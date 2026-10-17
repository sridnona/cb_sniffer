from collections import Counter, defaultdict
import gzip
import itertools
import pysam
import logging as logger
import sys
import argparse
#print(pysam.__version__)
#0.14.1
#python3
#3.6.5


class GenomicPosition:

	def __init__(self, line):
		"""
		:param line:
		 for each variant holds genomics positions
		 and other informations
		:raises ValueError: if the line has fewer than 7 tab-separated
		 columns or start/stop aren't integers
		"""
		if not line.strip():
			raise ValueError("blank line")
		region = line.rstrip('\r\n').split('\t')
		if len(region) < 7:
			raise ValueError("expected at least 7 tab-separated columns, found {}".format(len(region)))
		self.chrm = region[0]
		try:
			self.start = int(region[1])
			self.end = int(region[2])
		except ValueError:
			raise ValueError("start and stop must be integers")
		self.ref = region[3]
		self.alt = region[4]
		self.gene = region[5]  # gene name
		self.event = region[6]  #
		# print(self, self.chrm, self.start, self.ref, self.alt,self.event)
		self.normalize()

	def normalize(self):
		"""
		converts the variant to a minimal representation used for counting reads
		(the original start/ref/alt are kept for the output files)

		accepted inputs:
		  SNP / multi-base  T/C, ATC/GGC   start = first base
		  dash indels       C/-            start = first deleted base
		                    -/A            start = base before the insertion
		  VCF-style indels  ACT/A, A/AC    start = anchor base (VCF POS)

		sets:
		  self.kind     'snp', 'mnp', 'ins', 'del' or 'complex'
		  self.pos      1-based position to pile up on
		  self.var_ref  ref bases without padding ('' for an insertion)
		  self.var_alt  alt bases without padding ('' for a deletion)
		"""
		ref = '' if self.ref == '-' else self.ref.upper()
		alt = '' if self.alt == '-' else self.alt.upper()
		pos = self.start
		vcf_style = '-' not in (self.ref, self.alt)
		if vcf_style:
			# a VCF indel happens right after its anchor base, so drop the anchor
			# first; then trim shared trailing and leading bases. This keeps the
			# caller's position (A/AA stays after the anchor) while extra padding
			# across a repeat (ACCCC/ACCC) is trimmed back to the first base
			if len(ref) != len(alt) and ref[:1] == alt[:1]:
				ref, alt = ref[1:], alt[1:]
				pos += 1
			while ref and alt and ref[-1] == alt[-1]:
				ref, alt = ref[:-1], alt[:-1]
			while ref and alt and ref[0] == alt[0]:
				ref, alt = ref[1:], alt[1:]
				pos += 1
		if ref and not alt:
			kind = 'del'
		elif alt and not ref:
			kind = 'ins'
			if vcf_style:
				pos -= 1  # pile up on the base before the inserted bases
		elif ref and len(ref) == len(alt):
			kind = 'snp' if len(ref) == 1 else 'mnp'
		else:
			kind = 'complex'  # unequal-length substitution, or ref == alt
		self.kind, self.pos, self.var_ref, self.var_alt = kind, pos, ref, alt

	def classify(self):
		"""
		:return: 'snp', 'mnp', 'ins', 'del' or 'complex'
		"""
		return self.kind

	def read_allele(self, read, baseq):
		"""
		bases a read carries across a snp/mnp
		:return: base string, or None if the read doesn't cover every base
		without a gap or any base is below baseq
		"""
		n = len(self.var_alt)
		qpos = read.query_position
		if n > 1:
			aligned = dict((r, q) for q, r in read.alignment.get_aligned_pairs(matches_only=True))
			if any(aligned.get(self.pos - 1 + i) != qpos + i for i in range(n)):
				return None
		if min(read.alignment.query_qualities[qpos:qpos + n]) < int(baseq):
			return None
		return read.alignment.query_sequence[qpos:qpos + n]

	@classmethod
	def from_fields(cls, chrm, start, end, ref, alt, gene, event):
		return cls('\t'.join(str(f) for f in (chrm, start, end, ref, alt, gene, event)))

	@classmethod
	def good_barcodes(cls,f):
		"""

		:param f: good barcodes file, one per line, plain or gzipped
		:return: dict with good barcodes
		"""
		barc = {}
		with open_text(f) as barcodes:
			for line in barcodes:
				lines = line.strip()
				if lines:
					barc[lines] = 1
		return barc

	def count_barcodes(self, pile, bar, mapq=0, baseq=0, fasta=None):
		"""

		:param pile: open pysam.AlignmentFile
		:param bar: good barcodes dict
		:param mapq: mapping quality
		:param baseq: base quality (snp/mnp only)
		:param fasta: pysam.FastaFile used to left-align indels, or None
		:return: two dicts 1) CB:UB barcodes of all counted reads, by read name
		2) dict with ref and alt CB:UB barcodes
		"""
		barUcodes = defaultdict(list)  # UB
		bar_count = {'ref': [], 'alt': []}  # UB
		logger.info("processing variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(
			self.chrm,self.start, self.end, self.ref, self.alt, self.gene, self.event))
		mapq = int(mapq)

		def keep(aln):
			# cheap filters, checked before working out a read's allele
			if aln.mapping_quality < mapq:
				return False
			try:
				return aln.get_tag('CB') in bar and aln.has_tag('UB')  # only good barcodes
			except KeyError:  # no CB tag
				return False

		if self.kind in ('ins', 'del'):
			ref = None
			if fasta is not None:
				ref = RefWindow(fasta, self.chrm, self.pos - 1, self.pos + len(self.var_ref))
			reads = self.indel_reads(pile, ref, keep)
		else:
			reads = self.snp_reads(pile, baseq, keep)
		for aln, allele in reads:
			ustring = aln.get_tag('CB') + ':' + aln.get_tag('UB')
			barUcodes[aln.query_name].append(ustring)
			bar_count[allele].append(ustring)
		return barUcodes, bar_count

	def snp_reads(self, pile, baseq, keep):
		"""
		:param keep: function, False for reads to skip before checking their allele
		:return: generator of (alignment, 'ref' or 'alt') for a snp/mnp
		reads with a third allele (or N, or part of an mnp) are neither ref nor alt
		"""
		for pileupcolumn in pile.pileup(reference=self.chrm, start=self.pos - 1, end=self.pos,
										truncate=True, stepper="nofilter", max_depth=100000000,
										multiple_iterators=False,  # one pileup at a time; CRAM warns otherwise
										min_base_quality=0):
			for read in pileupcolumn.pileups:
				if read.is_del or read.is_refskip or not keep(read.alignment):
					continue
				allele = self.read_allele(read, baseq)
				if allele == self.var_ref:
					yield read.alignment, 'ref'
				elif allele == self.var_alt:
					yield read.alignment, 'alt'

	def indel_event(self):
		"""
		:return: the variant as (p, seq, ref_len), see left_align()
		"""
		if self.kind == 'del':
			return self.pos - 1, self.var_ref, len(self.var_ref)
		return self.pos, self.var_alt, 0

	def indel_reads(self, pile, ref, keep):
		"""
		:param ref: RefWindow around the variant, or None to match indels
		 only at the exact position given
		:param keep: function, False for reads to skip before checking their allele
		:return: generator of (alignment, 'ref' or 'alt') for an insertion/deletion

		The variant and each read's indels are left-aligned, so an indel the
		aligner placed elsewhere in a repeat still matches. A ref read must
		span the whole repeat (every position the indel could be placed at)
		with no indel or splice in it; other reads are neither ref nor alt.
		"""
		p, seq, ref_len = self.indel_event()
		if ref is not None and ref_len and ref.bases(p, ref_len) != seq:
			warn("deleted bases don't match the reference ({}): {}\t{}\t{}\t{}".format(
				ref.bases(p, ref_len), self.chrm, self.start, self.ref, self.alt))
		variant = left_align(p, seq, ref) + (ref_len,)
		right = right_align(p, seq, ref_len, ref)[0]
		lo, hi = variant[0] - 1, right + ref_len  # ref bases flanking the repeat
		for aln in pile.fetch(self.chrm, lo, hi + 1):
			if not keep(aln):
				continue
			allele = indel_allele(aln, lo, hi, variant, ref)
			if allele:
				yield aln, allele

	def consensus_calling(self,total_barcodes, ub_counts, wt, wu, wc):
		"""

		:param total_barcodes: list of all the barcodes at variant pos
		:param ub_counts: ref/alt barcode counts
		:param wt: file handle to write vaf information for each variant
		:param wu: file handle to write CB:UB information for each variant
		:param wc: file handle to write CB information for each variant
		:return: variant information

		"""
		uniqub_barcodes_raw = set()  # uni depth

		for tupru in total_barcodes.values():
			uniqub_barcodes_raw.update(tupru)
		# reads per CB:UB, counted once rather than per UMI
		ref_reads = Counter(ub_counts['ref'])
		alt_reads = Counter(ub_counts['alt'])
		# UB ====================================================================================================
		logger.info("Consensus calc variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(
			self.chrm, self.start, self.end, self.ref, self.alt, self.gene, self.event))
		d = {}
		UBuniq_barcodes_raw = []
		for utags in uniqub_barcodes_raw:
			unref = ref_reads[utags]
			unalt = alt_reads[utags]
			utotl = unref + unalt
			if unref and unalt:
				if unref / utotl < float(0.75):
					if unalt / utotl < float(0.75):
						continue
					else:
						un1alt = 1
						un1ref = 0
						u1totl = un1ref + un1alt
				else:
					if unalt / utotl < float(0.75):
						un1alt = 0
						un1ref = 1
						u1totl = un1ref + un1alt
					else:
						un1alt = 1  # this is sanity check in the final file
						un1ref = 1  # if both are one then the code block failed
						u1totl = un1ref + un1alt
			else:
				# uniq the ref UB barcode for vaf
				# if an UB is not seen in both ref and alt bin
				# meaning if a UB occurs more than once consider it only once
				# i.e if UB has 5 [alt] and 0[ref]
				# we would consider it 1[alt] and 0[ref]
				# (the UB file still gets the actual read counts, unref/unalt)
				un1ref = int(unref > 0)
				un1alt = int(unalt > 0)
				u1totl = un1ref + un1alt

			if un1alt:  # append only alt barcodes that have 1
				UBuniq_barcodes_raw.append(utags)

			ubtager = '{chrm}\t{st}\t{en}\t{ref}\t{alt}\t' \
					  '{type}\t{gene}\t{bar}\t{r}\t{v}\t{tot}\n'.format(
						chrm=self.chrm,
						st=self.start,
						en=self.end,
						ref=self.ref,
						alt=self.alt,
						type=self.event,
						gene=self.gene,
						bar=utags,
						r=unref,
						v=unalt, tot=utotl)
			# print(ubtager)
			wu.write(ubtager)

			# now work on CB barcodes
			# GGAAAGCCACCACGTG-2:GATGTCGCGC CB = GGAAAGCCACCACGTG-2 : UB = GATGTCGCGC
			if utags.split(':')[0] not in d:  # check if this already exists in dict(d)
				d[utags.split(':')[0]] = {'ref': un1ref, 'alt': un1alt}
			else:
				d[utags.split(':')[0]]['ref'] += un1ref
				d[utags.split(':')[0]]['alt'] += un1alt
		uni_alt = []
		tuni_alt = []
		CBuniq_barcodes_raw = []
		for i, v in d.items():
			# print i,v['alt']
			if v['alt'] >= 1:
				CBuniq_barcodes_raw.append(i)
			uni_alt.append(v['alt'])
			tot1 = v['alt'] + v['ref']
			tuni_alt.append(tot1)
			cbtager = '{chrm}\t{st}\t{en}\t{ref}\t{alt}\t' \
					  '{type}\t{gene}\t{bar}\t{ref1}\t{alt1}\t{tot}\n'.format(
						chrm=self.chrm,
						st=self.start,
						en=self.end,
						ref=self.ref,
						alt=self.alt,
						type=self.event,
						gene=self.gene,
						bar=i,
						alt1=v['alt'],
						ref1=v['ref'],
						tot=tot1)
			# print(cbtager)
			wc.write(cbtager)

		Uuni_alt_c = sum(uni_alt)
		tUuni_alt_c = sum(tuni_alt)  # total
		UBfin_umi = ','.join(UBuniq_barcodes_raw)
		CBfin_umi = ','.join(CBuniq_barcodes_raw)

		try:
			uvU = round(float(Uuni_alt_c) / float(tUuni_alt_c), 2)  # uni vaf
		except ZeroDivisionError:
			uvU = float(0.0)

		snp = '{chrm}\t{st}\t{en}\t{ref}\t{alt}\t{type}\t{gene}' \
			  '\t{Uunidepth}\t{Uunc}\t{Uuvaf}\t{umi}\t{Uumi}\n'.format(
				chrm=self.chrm,
				st=self.start,
				en=self.end,
				ref=self.ref,
				alt=self.alt,
				type=self.event,
				gene=self.gene,
				Uunidepth=tUuni_alt_c,
				Uumi=UBfin_umi,
				Uuvaf=uvU,
				umi=CBfin_umi,
				Uunc=Uuni_alt_c)
		# print(snp)
		wt.write(snp)
		logger.info("completed processing variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(
			self.chrm, self.start, self.end, self.ref, self.alt, self.gene, self.event))
		return snp


class RefWindow:
	"""reference bases around a variant, using 0-based coordinates"""

	def __init__(self, fasta, chrm, start, end, pad=1000):
		self.start = max(0, start - pad)
		self.seq = fasta.fetch(chrm, self.start, end + pad).upper()

	def base(self, i):
		i -= self.start
		return self.seq[i] if 0 <= i < len(self.seq) else None

	def bases(self, i, n):
		return self.seq[i - self.start:i - self.start + n]


def left_align(p, seq, ref):
	"""
	shift an indel left as far as it stays equivalent
	:param p: 0-based first reference base the indel affects: the first deleted
	 base for a deletion, the base after the inserted bases for an insertion
	:param seq: deleted or inserted bases
	:param ref: RefWindow, or None to leave the indel where it is
	:return: (p, seq) after shifting
	"""
	while ref is not None and seq and ref.base(p - 1) == seq[-1]:
		seq = seq[-1] + seq[:-1]
		p -= 1
	return p, seq


def right_align(p, seq, ref_len, ref):
	"""
	shift an indel right as far as it stays equivalent
	:param ref_len: len(seq) for a deletion, 0 for an insertion
	:return: (p, seq) after shifting, see left_align()
	"""
	while ref is not None and seq and ref.base(p + ref_len) == seq[0]:
		seq = seq[1:] + seq[0]
		p += 1
	return p, seq


def indel_allele(aln, lo, hi, variant, ref):
	"""
	:param aln: pysam.AlignedSegment
	:param lo, hi: 0-based reference bases flanking the repeat the indel is in
	:param variant: left-aligned (p, seq, ref_len) of the variant
	:param ref: RefWindow or None
	:return: 'alt' if the read's only indel between lo and hi is the variant,
	 'ref' if the read is aligned at lo and hi with no indel or splice
	 between them, otherwise None
	"""
	r, q = aln.reference_start, 0
	events = []
	covered = 0
	for op, n in aln.cigartuples:
		if op in (0, 7, 8):  # M = X
			covered += (r <= lo < r + n) + (r <= hi < r + n)
			r += n
			q += n
		elif op == 1:  # I
			if lo < r <= hi:
				if aln.query_sequence is None:
					return None
				events.append((r, aln.query_sequence[q:q + n], 0))
			q += n
		elif op in (2, 3):  # D N
			if r <= hi and r + n > lo:
				if op == 3:
					return None  # spliced across the repeat
				events.append((r, ref.bases(r, n) if ref is not None else None, n))
			r += n
		elif op == 4:  # S
			q += n
	if not events:
		return 'ref' if covered == 2 else None
	if len(events) > 1:
		return None
	r_p, r_seq, r_len = events[0]
	v_p, v_seq, v_len = variant
	if r_len != v_len or (r_len == 0 and len(r_seq) != len(v_seq)):
		return None  # different type or length
	if ref is None:
		# deleted bases aren't known without the reference
		same = r_p == v_p and (r_len > 0 or r_seq == v_seq)
	else:
		same = left_align(r_p, r_seq, ref) == (v_p, v_seq)
	return 'alt' if same else None


def warn(msg):
	logger.warning(msg)
	print('WARNING: ' + msg, file=sys.stderr)


def die(msg):
	logger.error(msg)
	print('ERROR: ' + msg, file=sys.stderr)
	sys.exit(1)


def is_cram(path):
	""":return: True if path is a CRAM file, judged by content"""
	with open(path, 'rb') as f:
		return f.read(4) == b'CRAM'


def is_gzipped(path):
	with open(path, 'rb') as f:
		return f.read(2) == b'\x1f\x8b'


def open_text(path):
	"""open a plain or gzipped (judged by content) text file for reading"""
	return gzip.open(path, 'rt') if is_gzipped(path) else open(path, 'r')


def is_vcf(path):
	"""
	:param path: variant file
	:return: True for VCF (plain, gzip or bgzip) or BCF, judged by content
	"""
	with (gzip.open(path, 'rb') if is_gzipped(path) else open(path, 'rb')) as f:
		head = f.read(16)
	return head.startswith(b'##fileformat=VCF') or head.startswith(b'BCF')


def is_tsv_header(line):
	"""
	:param line: first line of a tsv variant file
	:return: True if it's a header: it starts with '#', or its start column isn't an integer
	"""
	fields = line.rstrip('\r\n').split('\t')
	if line.startswith('#'):
		return True
	if len(fields) < 2:  # blank or not tab-separated: let the parser report it
		return False
	try:
		int(fields[1])
	except ValueError:
		return True
	return False


def read_variants(path, pass_only=False):
	"""
	:param path: tsv variant file (header line optional), or VCF/BCF
	:param pass_only: VCF only, skip records whose FILTER isn't PASS or '.'
	:return: generator of GenomicPosition, one per ALT allele for VCF records
	gene is taken from the VCF ID column and type is the variant kind
	"""
	if not is_vcf(path):
		with open(path, 'r') as regions:
			first = next(regions, None)
			if first is None:
				die("{}: variant file is empty".format(path))
			start = 2
			if not is_tsv_header(first):
				regions, start = itertools.chain([first], regions), 1
			for line_num, lines in enumerate(regions, start=start):
				try:
					x = GenomicPosition(lines)
				except ValueError as e:
					die("{}: line {}: {}: {!r}".format(path, line_num, e, lines.rstrip('\r\n')))
				yield x
		return
	verbosity = pysam.set_verbosity(0)  # no htslib error about a missing index
	try:
		vcf = pysam.VariantFile(path)
	finally:
		pysam.set_verbosity(verbosity)
	with vcf:
		for rec in vcf:
			if pass_only and list(rec.filter.keys()) not in ([], ['PASS']):
				logger.info("skipping filtered variant: {}\t{}\t{}\t{}".format(
					rec.chrom, rec.pos, rec.ref, ','.join(rec.alts or ())))
				continue
			for alt in rec.alts or ():
				if not alt or set(alt.upper()) - set('ACGTN'):
					warn("skipping symbolic/breakend allele: {}\t{}\t{}\t{}".format(
						rec.chrom, rec.pos, rec.ref, alt))
					continue
				x = GenomicPosition.from_fields(rec.chrom, rec.pos, rec.stop, rec.ref, alt,
												rec.id or '.', '.')
				x.event = x.kind
				yield x


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Parse CB barcodes from Single cell rna seq data')
	parser.add_argument('bam_file', help='BAM or CRAM file (CRAM requires -r)')
	parser.add_argument('variant_file', help='variants: tsv file (header line optional), or VCF/BCF (optionally gzipped)')
	parser.add_argument('barcodes', help='list of good barcodes file (may be gzipped)')
	parser.add_argument('upn', help='upn/sample name: will be used as prefix for out_file')
	parser.add_argument('-f', "--filter", type=int, default=0,
						help='number of reads required per barcode default: 0')
	parser.add_argument('-mq', "--mapq", type=int, default=0,
						help='Skip read with mapq smaller than default : 0')
	parser.add_argument('-bq', "--baseq", type=int, default=1,
						help='Skip bases with base quality less than default : 1')
	parser.add_argument('-r', '--reference',
						help='reference FASTA (indexed, may be bgzipped); used to left-align '
							 'indels so ones in repeats are found wherever the aligner put them. '
							 'Required for CRAM input, which is decoded with it')
	parser.add_argument("--pass-only", action='store_true',
						help='VCF input: skip records whose FILTER is not PASS or .')

	args = parser.parse_args()

	bam_file = args.bam_file
	variants = args.variant_file
	barcodes_good = args.barcodes
	outfile = args.upn
	counts = args.filter
	bq = args.baseq
	mq = args.mapq
	logger.basicConfig(filename=outfile + '.log', filemode='w+',
					   level=logger.DEBUG,
					   format='%(asctime)s %(levelname)s %(message)s')
	logger.info("Start process")
	cram = is_cram(bam_file)
	if cram and not args.reference:
		die("{}: CRAM input requires a reference FASTA (-r/--reference)".format(bam_file))
	variant_list = list(read_variants(variants, args.pass_only))
	logger.info("Number of variants:\t{}".format(len(variant_list)))
	fasta = None
	indel_chroms = {x.chrm for x in variant_list if x.kind in ('ins', 'del')}
	if args.reference:
		fasta = pysam.FastaFile(args.reference)
		missing = sorted(indel_chroms - set(fasta.references))
		if missing:
			die("{}: no sequence for chromosome(s) {}".format(args.reference, ', '.join(missing)))
	elif indel_chroms:
		warn("no --reference given: indels are only counted at the exact position given, "
			 "not elsewhere in a repeat")

	bars = GenomicPosition.good_barcodes(barcodes_good)
	logger.info("Number of good barcodes:\t{}".format(len(bars)))
	with open(outfile + '_AllCounts.tsv', 'w+') as var, \
			open(outfile + '_counts_CB.tsv', 'w+') as CB, \
			open(outfile + '_counts_UB.tsv', 'w+') as UB, \
			pysam.AlignmentFile(bam_file, 'rc' if cram else 'rb',  # opened once for all variants
								reference_filename=args.reference if cram else None) as bam:
		# write headers
		var_header = ['chr', 'start', 'end', 'ref', 'alt', 'type', 'gene',
					  'UB_DEPTH', 'UB_ALT', 'UB_VAF', 'CB_barcodes', 'CB:UB_barcodes']
		var_h = '\t'.join(var_header) + '\n'
		var.write(var_h)

		CB_header = ['chrm', 'start', 'end', 'ref', 'alt',
					 'type', 'gene', 'barcode', 'ref_count', 'alt_count', 'total_CB']
		CB_h = '\t'.join(CB_header) + '\n'
		CB.write(CB_h)

		UB_header = ['chrm', 'start', 'end', 'ref', 'alt', 'type',
					 'gene', 'barcode', 'ref_count', 'alt_count', 'total_UB']
		UB_h = '\t'.join(UB_header) + '\n'
		UB.write(UB_h)

		for x in variant_list:
			print(x.chrm, x.start, x.event, x.gene, x.classify())
			if x.kind == 'complex':
				warn("skipping unsupported variant (ref and alt differ in length "
					 "after trimming shared bases): {}\t{}\t{}\t{}".format(x.chrm, x.start, x.ref, x.alt))
				continue
			barcode_counts = x.count_barcodes(bam, bars, mq, bq, fasta)
			# print(check_classify[0])
			counters = x.consensus_calling(barcode_counts[0], barcode_counts[1], var, UB, CB)

logger.info("end process")

