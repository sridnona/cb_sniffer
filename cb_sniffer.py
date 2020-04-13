from collections import defaultdict
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
		"""
		region = line.strip().split('\t')
		if len(region) < 7:
			logger.debug("Insufficient columns fields please check variant file")
			sys.exit(1)
		self.chrm = region[0]
		self.start = int(region[1])
		self.end = int(region[2])
		self.ref = region[3]
		self.alt = region[4]
		self.gene = region[5]  # gene name
		self.event = region[6]  #
		# print(self, self.chrm, self.start, self.ref, self.alt,self.event)

	def classify(self):
		"""
		for each variant classifies if its a snp or variant
		:return:
		0 for snp
		int if its insertion or deletion
		"""
		v = ['A', 'T', 'G', 'C']
		if self.ref and self.alt in v:
			return 0

		else:
			if self.alt == '-':  # deletion
				return len(self.ref)
			elif self.ref == '-':  # insertion
				return len(self.alt)
			elif len(self.ref) == len(self.alt) > 1:
				return len(self.alt)

	@classmethod
	def good_barcodes(cls,f):
		"""

		:param f: good barcodes file
		:return: dict with good barcodes
		"""
		barc = {}
		with open(f, 'r') as barcodes:
			for line in barcodes:
				lines = line.strip()
				barc[lines] = 1
		return barc

	def count_barcodes(self, bam_fil, bar, indel, mapq=0, baseq=0):
		"""

		:param bam_fil: bam_file
		:param bar: good barcodes dict
		:param indel: SNP = 0 indel= int(window size) insertion or deletion
		:param mapq: mapping quality
		:param baseq: base quality
		:return: two dicts 1) list of all barcodes at given var pos
		2) dict with ref and alt barcodes
		"""
		barcodes = defaultdict(list)  # CB
		# ======================================
		barUcodes = defaultdict(list)  # UB
		bar_count = {'ref': [], 'alt': []}  # UB
		with pysam.AlignmentFile(bam_fil, 'rb') as pile:
			logger.info("processing variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(
				self.chrm,self.start, self.end, self.ref, self.alt, self.gene, self.event))

			for pileupcolumn in pile.pileup(reference=self.chrm, start=self.start - 1, end=self.start,
											truncate=True, stepper="nofilter", max_depth=100000000):
				for read in pileupcolumn.pileups:

					if read.alignment.has_tag('CB') and read.alignment.has_tag('UB'):
						if read.alignment.get_tag('CB') not in bar:  # filter reads with only good barcodes
							continue
							
						# print(read.query_position)
						if indel > 0:
							# print(int(read.alignment.query_qualities[read.query_position]))

							# print('{}'.format(read.query_position))
							# if read.query_position == None:
							# 	continue
							# if read.alignment.query_qualities[read.query_position + indel + 1] == None:
							# 	continue
							if read.is_refskip or \
							int(read.alignment.mapping_quality) < int(mapq):
								continue
						else:

							if read.is_del or read.is_refskip or \
							int(read.alignment.query_qualities[read.query_position]) < int(baseq) or \
							int(read.alignment.mapping_quality) < int(mapq):
								continue

						q_name = read.alignment.query_name
						q_tag = read.alignment.get_tag('CB')
						qU_tag = read.alignment.get_tag('UB')
						qstring = q_tag + ':' + qU_tag
						barcodes['{}'.format(q_name)].append(q_tag)
						barUcodes['{}'.format(q_name)].append(qstring)

						# barcode without mutations
						if indel > 0:
							
							# print('{}:{}:{}:{}:{}:{}'.format(read.indel, read.alignment.get_tag('CB'),read.alignment.cigarstring,
							# 	read.query_position,self.alt,pileupcolumn.pos))
							q_name = read.alignment.query_name
							q_tag = read.alignment.get_tag('CB')
							qU_tag = read.alignment.get_tag('UB')
							qstring = q_tag + ':' + qU_tag
							barcodes['{}'.format(q_name)].append(q_tag)
							barUcodes['{}'.format(q_name)].append(qstring)
							if self.alt == '-':
								if read.alignment.has_tag('CB') and read.query_position == None:
								# and read.indel == indel and \
								# 				read.alignment.query_alignment_sequence[read.query_position +
								# 						1:read.query_position + read.indel + 1] == self.alt:
									# print('indel length', read.indel)
									alt_tag = read.alignment.get_tag('CB')
									# UB
									altU_tag = read.alignment.get_tag('UB')
									ustringa = alt_tag + ':' + altU_tag
									bar_count['alt'].append(ustringa)
								else:  # ref reads
									# CB
									ref_tag = read.alignment.get_tag('CB')
									# UB
									refU_tag = read.alignment.get_tag('UB')
									ustring = ref_tag + ':' + refU_tag
									bar_count['ref'].append(ustring)
									# print(ustring)
									# print('{}\t{}\t{}'.format(ref_tag, refU_tag, 'ref'))
							else:
								if read.alignment.has_tag('CB') and read.indel == indel and \
									read.alignment.query_alignment_sequence[read.query_position +
										1:read.query_position + read.indel + 1] == self.alt:
									# print('indel length', read.indel)
									alt_tag = read.alignment.get_tag('CB')
									# UB
									altU_tag = read.alignment.get_tag('UB')
									ustringa = alt_tag + ':' + altU_tag
									bar_count['alt'].append(ustringa)
								else:  # ref reads
									# CB
									ref_tag = read.alignment.get_tag('CB')
									# UB
									refU_tag = read.alignment.get_tag('UB')
									ustring = ref_tag + ':' + refU_tag
									bar_count['ref'].append(ustring)
									# print(ustring)
									# print('{}\t{}\t{}'.format(ref_tag, refU_tag, 'ref'))
						else:
							#print(indel,read.query_position,read.alignment.query_sequence[read.query_position])
							if read.alignment.query_sequence[read.query_position] != self.alt:
								ref_tag = read.alignment.get_tag('CB')
								# UB
								refU_tag = read.alignment.get_tag('UB')
								ustring = ref_tag + ':' + refU_tag
								bar_count['ref'].append(ustring)
								# print('{}\t{}\t{}'.format(ref_tag, refU_tag, 'ref'))
							elif read.alignment.query_sequence[read.query_position] == self.alt:
								# CB
								alt_tag = read.alignment.get_tag('CB')
								# UB
								altU_tag = read.alignment.get_tag('UB')
								ustringa = alt_tag + ':' + altU_tag
								bar_count['alt'].append(ustringa)
		return barUcodes, bar_count

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
		ub_raw = []  # raw depth

		for tupru in total_barcodes.values():
			for sampleru in tupru:
				ub_raw.append(sampleru)
				uniqub_barcodes_raw.add(sampleru)
		# UB ====================================================================================================
		logger.info("Consensus calc variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(
			self.chrm, self.start, self.end, self.ref, self.alt, self.gene, self.event))
		d = {}
		UBuniq_barcodes_raw = []
		for utags in uniqub_barcodes_raw:
			if utags in ub_counts['alt'] and utags in ub_counts['ref']:
				# print 'both: ', utags)
				unref = ub_counts['ref'].count(utags)
				unalt = ub_counts['alt'].count(utags)
				utotl = unref + unalt
				un1ref = ub_counts['ref'].count(utags)
				un1alt = ub_counts['alt'].count(utags)
				ut1totl = un1ref + un1alt
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
				un1ref = list(set(ub_counts['ref'])).count(utags)
				# same stuff for alt UB barcodes
				un1alt = list(set(ub_counts['alt'])).count(utags)
				u1totl = un1ref + un1alt
				# although we Uniq the UB's above we need to print out the actual number
				# UB has 5 [alt] and 0[ref]
				# print out the same number
				unref = ub_counts['ref'].count(utags)
				unalt = ub_counts['alt'].count(utags)
				utotl = unref + unalt

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


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Parse CB barcodes from Single cell rna seq data')
	parser.add_argument('bam_file', help='BAM file')
	parser.add_argument('variant_file', help='variants file with header')
	parser.add_argument('barcodes', help='list of good barcodes file')
	parser.add_argument('upn', help='upn/sample name: will be used as prefix for out_file')
	parser.add_argument('-f', "--filter", type=int, default=0,
						help='number of reads required per barcode default: 0')
	parser.add_argument('-mq', "--mapq", type=int, default=0,
						help='Skip read with mapq smaller than default : 0')
	parser.add_argument('-bq', "--baseq", type=int, default=1,
						help='Skip bases with base quality less than default : 1')

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
	num_lines_variants = sum(1 for line in open(variants)) - 1  # minus header
	logger.info("Number of variants:\t{}".format(num_lines_variants))

	num_lines_barcode = sum(1 for lines in open(barcodes_good))
	logger.info("Number of good barcodes:\t{}".format(num_lines_barcode))
	bars = GenomicPosition.good_barcodes(barcodes_good)
	with open(outfile + '_AllCounts.tsv', 'w+') as var, \
			open(outfile + '_counts_CB.tsv', 'w+') as CB, \
			open(outfile + '_counts_UB.tsv', 'w+') as UB, \
			open(variants, 'r') as regions:
		# write headers
		var_header = ['chr', 'start', 'end', 'ref', 'alt', 'gene', 'type',
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

		next(regions)
		for lines in regions:
			x = GenomicPosition(lines)

			print(x.chrm, x.start, x.event, x.gene, x.classify())
			barcode_counts = x.count_barcodes(bam_file, bars, x.classify(), mq, bq)
			# print(check_classify[0])
			counters = x.consensus_calling(barcode_counts[0], barcode_counts[1], var, UB, CB)

logger.info("end process")

