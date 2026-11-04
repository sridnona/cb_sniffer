"""Tests for cb_sniffer.py using the small mouse 10x data set in this directory.

Run from the repository root with:  pytest test/
Tests marked xfail document known bugs listed in issues.md. They are strict,
so they start failing ("XPASS") once the bug is fixed, as a reminder to remove
the marker.
"""
import csv
import io
import os
import shutil
import subprocess
import sys

import pytest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(TEST_DIR)
SCRIPT = os.path.join(REPO_DIR, 'cb_sniffer.py')
BAM = os.path.join(TEST_DIR, 'test.bam')
VARIANTS = os.path.join(TEST_DIR, 'variants.tsv')
VARIANTS_VCF = os.path.join(TEST_DIR, 'variants_vcf.tsv')
VCF = os.path.join(TEST_DIR, 'variants.vcf')
REFERENCE = os.path.join(TEST_DIR, 'reference.fa.gz')  # N except +-3 kb around the sites
BARCODES = os.path.join(TEST_DIR, 'barcodes.tsv')

sys.path.insert(0, REPO_DIR)
from cb_sniffer import GenomicPosition, is_vcf, read_variants, left_align, right_align  # noqa: E402

DEL, INS, SNV, SNV_MQ, NOCOV, MNP = 12904213, 13562395, 13625741, 15268954, 10000001, 10066353
# start position of each site in variants_vcf.tsv -> the same site in variants.tsv
VCF_SITES = {12904212: DEL, 13562395: INS, 13625741: SNV, 10066352: MNP}
# (start, alt) of each usable record in variants.vcf -> the same site in variants.tsv
VCF_FILE_SITES = {(10000001, 'G'): NOCOV, (10066352, 'CCATA'): MNP, (12904212, 'ACCC'): DEL,
				  (13562395, 'AA'): INS, (13625741, 'C'): SNV, (15268954, 'C'): SNV_MQ}
COUNT_COLS = ('UB_DEPTH', 'UB_ALT', 'UB_VAF', 'CB_barcodes', 'CB:UB_barcodes')


def variant(ref, alt, start=100):
	return GenomicPosition('1\t{0}\t{0}\t{1}\t{2}\tgene\ttype\n'.format(start, ref, alt))


def read_tsv(path):
	with open(path) as f:
		return list(csv.DictReader(f, delimiter='\t'))


def run(tmp_path, *opts, variants=VARIANTS, prefix='out', bam=BAM, env=None):
	"""Run cb_sniffer.py and return its three output tables, keyed by output name."""
	cmd = [sys.executable, SCRIPT, *opts, bam, variants, BARCODES, prefix]
	subprocess.run(cmd, cwd=str(tmp_path), check=True, stdout=subprocess.DEVNULL, env=env)
	return {name: read_tsv(os.path.join(str(tmp_path), '{}_{}.tsv'.format(prefix, name)))
			for name in ('AllCounts', 'counts_CB', 'counts_UB')}


def by_site(rows):
	return {int(r['start']): r for r in rows}


def normalized(row):
	"""Summary row with barcode lists sorted; their order varies between runs (issues.md #9)."""
	row = dict(row)
	for col in ('CB_barcodes', 'CB:UB_barcodes'):
		row[col] = sorted(row[col].split(','))
	return row


def read_counts(ub_rows, site):
	"""Total (ref, alt) reads at a site, summed over all UMIs."""
	rows = [r for r in ub_rows if int(r['start']) == site]
	return (sum(int(r['ref_count']) for r in rows),
			sum(int(r['alt_count']) for r in rows))


@pytest.fixture(scope='module')
def default_run(tmp_path_factory):
	return run(tmp_path_factory.mktemp('default'))


# --- normalize() / classify() ---------------------------------------------------

@pytest.mark.parametrize('ref,alt,kind,pos,var_ref,var_alt', [
	('T', 'C', 'snp', 100, 'T', 'C'),
	('t', 'c', 'snp', 100, 'T', 'C'),           # lower case
	('ATC', 'GGC', 'mnp', 100, 'AT', 'GG'),     # shared trailing base trimmed
	('CTATG', 'CCATA', 'mnp', 101, 'TATG', 'CATA'),  # shared leading base trimmed
	('TG', 'CG', 'snp', 100, 'T', 'C'),         # padded SNV
	# dash notation
	('C', '-', 'del', 100, 'C', ''),            # start = first deleted base
	('CAG', '-', 'del', 100, 'CAG', ''),
	('-', 'A', 'ins', 100, '', 'A'),            # start = base before the insertion
	('-', 'AT', 'ins', 100, '', 'AT'),
	# VCF style: start = anchor base
	('ACT', 'A', 'del', 101, 'CT', ''),
	('A', 'AC', 'ins', 100, '', 'C'),
	('A', 'AA', 'ins', 100, '', 'A'),           # insertion stays after the anchor
	('ACCCC', 'ACCC', 'del', 101, 'C', ''),     # padding across a repeat is trimmed
	('GACT', 'GA', 'del', 102, 'CT', ''),       # more than one shared leading base
	('CA', 'A', 'del', 100, 'C', ''),           # no anchor base
	# not supported
	('ATC', 'GG', 'complex', 100, 'ATC', 'GG'),
	('A', 'A', 'complex', 100, '', ''),
])
def test_normalize(ref, alt, kind, pos, var_ref, var_alt):
	v = variant(ref, alt)
	assert (v.classify(), v.pos, v.var_ref, v.var_alt) == (kind, pos, var_ref, var_alt)
	assert (v.start, v.ref, v.alt) == (100, ref, alt)  # originals kept for output


# --- read_allele() ---------------------------------------------------------------

class FakePileupRead:
	def __init__(self, alignment, query_position):
		self.alignment, self.query_position = alignment, query_position


def make_read(cigar, quals=None, start=95):
	import pysam
	a = pysam.AlignedSegment()
	a.query_sequence = 'AAAAACATAGGGGG'
	a.reference_start = start - 1
	a.cigarstring = cigar
	a.query_qualities = pysam.qualitystring_to_array(quals or 'F' * 14)
	return FakePileupRead(a, 5)  # query position 5 is reference position 100


@pytest.mark.parametrize('cigar,quals,expected', [
	('14M', None, 'CATA'),
	('14M', 'FFFFFFF!FFFFFF', None),  # one low-quality base inside the MNP
	('7M50N7M', None, None),          # spliced inside the MNP
	('7M1D7M', None, None),           # deletion inside the MNP
	('7M1I6M', None, None),           # insertion inside the MNP
])
def test_read_allele_mnp(cigar, quals, expected):
	assert variant('TATG', 'CATA').read_allele(make_read(cigar, quals), 20) == expected


def test_read_allele_snp():
	assert variant('T', 'C').read_allele(make_read('14M'), 20) == 'C'
	assert variant('T', 'C').read_allele(make_read('14M', 'FFFFF!FFFFFFFF'), 20) is None


@pytest.mark.parametrize('line,message', [
	('\n', 'blank line'),
	('  \t \n', 'blank line'),
	('1\t100\t100\tA\tT\n', 'expected at least 7 tab-separated columns, found 5'),
	('1 100 100 A T g snp\n', 'expected at least 7 tab-separated columns, found 1'),
	('1\t1e2\t100\tA\tT\tg\tsnp\n', 'start and stop must be integers'),
])
def test_malformed_line_raises(line, message):
	with pytest.raises(ValueError, match=message):
		GenomicPosition(line)


def test_extra_columns_and_crlf_allowed():
	v = GenomicPosition('1\t100\t100\tA\tT\tg\tsnp\textra\r\n')
	assert (v.start, v.ref, v.alt, v.event) == (100, 'A', 'T', 'snp')
	assert GenomicPosition('1\t100\t100\tA\tT\tg\t\n').event == ''  # empty last column


def test_good_barcodes(tmp_path):
	f = tmp_path / 'bc.tsv'
	f.write_text('AAAC-1\nGGGT-1\n')
	assert GenomicPosition.good_barcodes(str(f)) == {'AAAC-1': 1, 'GGGT-1': 1}


# --- consensus_calling() --------------------------------------------------------

def call_consensus(ub_counts):
	"""Run consensus_calling on synthetic CB:UB read lists; return (AllCounts, CB, UB) text."""
	all_ub = set(ub_counts['ref']) | set(ub_counts['alt'])
	total = {'read{}'.format(i): [u] for i, u in enumerate(sorted(all_ub))}
	wt, wu, wc = io.StringIO(), io.StringIO(), io.StringIO()
	variant('T', 'C').consensus_calling(total, ub_counts, wt, wu, wc)
	return wt.getvalue().rstrip('\n').split('\t'), wc.getvalue(), wu.getvalue()


def test_consensus_umi_collapse():
	# cell1/umi1: 3 alt + 1 ref reads -> alt (>= 75%)
	# cell1/umi2: 1 alt + 1 ref       -> ambiguous, dropped
	# cell2/umi3: 4 ref reads         -> one ref UMI
	ub = {'alt': ['cell1:umi1'] * 3 + ['cell1:umi2'],
		  'ref': ['cell1:umi1', 'cell1:umi2'] + ['cell2:umi3'] * 4}
	summary, cb_text, ub_text = call_consensus(ub)
	depth, alt, vaf, cbs, umis = summary[7:12]
	assert (depth, alt, vaf) == ('2', '1', '0.5')
	assert cbs == 'cell1'
	assert umis == 'cell1:umi1'
	assert 'cell1:umi2' not in ub_text  # ambiguous UMIs are dropped from the UB table
	cb_rows = {line.split('\t')[7]: line.split('\t')[8:11] for line in cb_text.splitlines()}
	assert cb_rows == {'cell1': ['0', '1', '1'], 'cell2': ['1', '0', '1']}


def test_consensus_no_reads():
	summary, cb_text, ub_text = call_consensus({'ref': [], 'alt': []})
	assert summary[7:10] == ['0', '0', '0.0']
	assert cb_text == ub_text == ''


# --- end to end on test.bam ----------------------------------------------------

def test_output_headers(default_run):
	assert list(default_run['AllCounts'][0].keys()) == [
		'chr', 'start', 'end', 'ref', 'alt', 'type', 'gene',
		'UB_DEPTH', 'UB_ALT', 'UB_VAF', 'CB_barcodes', 'CB:UB_barcodes']
	for name, total in (('counts_CB', 'total_CB'), ('counts_UB', 'total_UB')):
		assert list(default_run[name][0].keys()) == [
			'chrm', 'start', 'end', 'ref', 'alt', 'type', 'gene',
			'barcode', 'ref_count', 'alt_count', total]


def test_type_and_gene_columns(default_run, vcf_run):
	row = by_site(default_run['AllCounts'])[DEL]
	assert (row['type'], row['gene']) == ('deletion', 'Slco5a1')
	row = [r for r in vcf_run if int(r['start']) == 12904212][0]
	assert (row['type'], row['gene']) == ('del', 'Slco5a1')
	for name in ('counts_CB', 'counts_UB'):
		row = [r for r in default_run[name] if int(r['start']) == DEL][0]
		assert (row['type'], row['gene']) == ('deletion', 'Slco5a1')


def test_one_summary_row_per_variant(default_run):
	assert [int(r['start']) for r in default_run['AllCounts']] == [DEL, INS, SNV, SNV_MQ, NOCOV, MNP]


@pytest.mark.parametrize('site,depth,alt,vaf', [
	(DEL, 3, 1, 0.33),
	(INS, 2, 1, 0.5),
	(SNV, 28, 14, 0.5),
	(SNV_MQ, 20, 13, 0.65),
	(NOCOV, 0, 0, 0.0),
	(MNP, 95, 73, 0.77),
])
def test_default_counts(default_run, site, depth, alt, vaf):
	row = by_site(default_run['AllCounts'])[site]
	assert int(row['UB_DEPTH']) == depth
	assert int(row['UB_ALT']) == alt
	assert float(row['UB_VAF']) == vaf


def test_alt_barcode_lists_match_alt_count(default_run):
	for row in default_run['AllCounts']:
		umis = [u for u in row['CB:UB_barcodes'].split(',') if u]
		assert len(umis) == int(row['UB_ALT'])
		cells = {u.split(':')[0] for u in umis}
		assert cells == {c for c in row['CB_barcodes'].split(',') if c}


def test_no_coverage_site_is_empty(default_run):
	row = by_site(default_run['AllCounts'])[NOCOV]
	assert row['CB_barcodes'] == row['CB:UB_barcodes'] == ''
	for name in ('counts_CB', 'counts_UB'):
		assert not [r for r in default_run[name] if int(r['start']) == NOCOV]


def test_only_whitelisted_barcodes_reported(default_run):
	good = set(open(BARCODES).read().split())
	cbs = {r['barcode'] for r in default_run['counts_CB']}
	cbs |= {r['barcode'].split(':')[0] for r in default_run['counts_UB']}
	assert cbs and cbs <= good
	# the BAM does contain reads from cells that are not in barcodes.tsv
	bam_cbs = set(subprocess.run(
		[sys.executable, '-c',
		 'import pysam,sys\nfor r in pysam.AlignmentFile(sys.argv[1]):\n'
		 ' print(r.get_tag("CB")) if r.has_tag("CB") else None', BAM],
		check=True, capture_output=True, text=True).stdout.split())
	assert bam_cbs - good


def test_cb_table_sums_match_summary(default_run):
	summary = by_site(default_run['AllCounts'])
	for site, row in summary.items():
		cb = [r for r in default_run['counts_CB'] if int(r['start']) == site]
		assert sum(int(r['alt_count']) for r in cb) == int(row['UB_ALT'])
		assert sum(int(r['total_CB']) for r in cb) == int(row['UB_DEPTH'])


def test_mapq_filter(tmp_path, default_run):
	mq = by_site(run(tmp_path, '-mq', '10')['AllCounts'])
	default = by_site(default_run['AllCounts'])
	# the multimapping reads (MAPQ 1 and 3) at SNV_MQ are removed
	assert (int(mq[SNV_MQ]['UB_DEPTH']), int(mq[SNV_MQ]['UB_ALT'])) == (16, 10)
	# all reads at the other sites have MAPQ 255, so nothing changes
	for site in (DEL, INS, SNV, NOCOV, MNP):
		assert normalized(mq[site]) == normalized(default[site])


def test_baseq_filter(tmp_path, default_run):
	# low-quality bases only remove reads, they don't change UMI totals here
	counts = {bq: read_counts(run(tmp_path, '-bq', bq, prefix='bq' + bq)['counts_UB'], SNV)
			  for bq in ('0', '30')}
	assert counts['0'] == (57, 52)    # BQ 0 bases are kept (needs min_base_quality=0, issues.md #2)
	assert read_counts(default_run['counts_UB'], SNV) == (52, 46)
	assert counts['30'] == (48, 42)


def test_deletion_ignores_baseq(tmp_path, default_run):
	# the indel branch has no base-quality filter
	high = run(tmp_path, '-bq', '60')['counts_UB']
	assert read_counts(high, DEL) == read_counts(default_run['counts_UB'], DEL) == (6, 9)


def test_insertion_ignores_baseq(tmp_path, default_run):
	high = run(tmp_path, '-bq', '60')['counts_UB']
	assert read_counts(high, INS) == read_counts(default_run['counts_UB'], INS)


def read_from_ref(refseq, cigar, start=80):
	"""read sequence matching refseq for a cigar, with A for inserted bases"""
	import pysam
	a = pysam.AlignedSegment()
	a.cigarstring = cigar
	seq, r = '', start
	for op, n in a.cigartuples:
		if op in (0, 7, 8):
			seq += refseq[r:r + n]
		elif op in (1, 4):
			seq += 'A' * n
		if op in (0, 2, 3, 7, 8):
			r += n
	return seq


def run_synthetic(tmp_path, reads, variant_line, refseq=None):
	"""
	Write reads to a BAM (each its own cell and UMI, all starting at 1-based
	position 81 on contig 1), run cb_sniffer on one variant and return its
	AllCounts row. reads: list of (cigar, sequence or None). A None sequence
	is all A's, or follows refseq when given; refseq is also passed as -r.
	"""
	import pysam
	bam = str(tmp_path / 'synthetic.bam')
	header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 1000}]}
	with pysam.AlignmentFile(bam, 'wb', header=header) as out:
		for i, (cigar, seq) in enumerate(reads):
			a = pysam.AlignedSegment()
			a.query_name, a.cigarstring = 'read{}'.format(i), cigar
			if seq is None:
				seq = read_from_ref(refseq, cigar) if refseq else 'A' * a.infer_query_length()
			a.query_sequence = seq
			a.reference_id, a.reference_start, a.mapping_quality = 0, 80, 255
			a.query_qualities = pysam.qualitystring_to_array('F' * len(a.query_sequence))
			a.set_tags([('CB', 'CELL{}-1'.format(i)), ('UB', 'UMI{}'.format(i))])
			out.write(a)
	pysam.index(bam)
	(tmp_path / 'bc.tsv').write_text(''.join('CELL{}-1\n'.format(i) for i in range(len(reads))))
	(tmp_path / 'v.tsv').write_text('chrm\tstart\tstop\tref\tvar\tgene\ttype\n' + variant_line + '\n')
	opts = []
	if refseq:
		(tmp_path / 'ref.fa').write_text('>1\n' + refseq + '\n')
		pysam.faidx(str(tmp_path / 'ref.fa'))
		opts = ['-r', 'ref.fa']
	subprocess.run([sys.executable, SCRIPT, *opts, bam, 'v.tsv', 'bc.tsv', 'out'],
				   cwd=str(tmp_path), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	return read_tsv(str(tmp_path / 'out_AllCounts.tsv'))[0]


def test_insertion_in_soft_clipped_read(tmp_path):
	"""issues.md #4: the inserted bases must be read from the full read sequence"""
	row = run_synthetic(tmp_path, [
		('20M2I20M', 'A' * 20 + 'GT' + 'A' * 20),               # alt
		('5S20M2I20M', 'CCCCC' + 'A' * 20 + 'GT' + 'A' * 20),   # alt, soft-clipped
		('5S40M', 'CCCCC' + 'A' * 40),                          # ref, soft-clipped
	], '1\t100\t100\t-\tGT\tg\tinsertion')
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('3', '2')
	assert sorted(row['CB_barcodes'].split(',')) == ['CELL0-1', 'CELL1-1']


@pytest.mark.parametrize('variant_line', [
	'1\t100\t100\tA\t-\tg\tdeletion',     # dash style: start = deleted base
	'1\t99\t100\tAA\tA\tg\tdeletion',     # VCF style: start = anchor base
])
def test_deletion_must_match(tmp_path, variant_line):
	"""issues.md #5: only a deletion with the variant's start and length is alt"""
	row = run_synthetic(tmp_path, [
		('19M1D20M', None),   # 1-bp deletion of base 100: alt
		('19M2D20M', None),   # 2-bp deletion starting at 100: skipped
		('18M2D20M', None),   # deletion starting at 99, covers 100: skipped
		('40M', None),        # no deletion: ref
		('20M1D20M', None),   # deletion at 101: doesn't span the site, skipped
	], variant_line)
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('2', '1')
	assert row['CB_barcodes'] == 'CELL0-1'


def base_at_100(bases):
	"""40-base read starting at 81 with the given bases from position 100"""
	return ('40M', 'A' * 19 + bases + 'A' * (21 - len(bases)))


def test_snp_other_alleles_skipped(tmp_path):
	"""issues.md #6: a read with neither the ref nor the alt base is not ref"""
	row = run_synthetic(tmp_path, [
		base_at_100('A'),   # ref
		base_at_100('C'),   # alt
		base_at_100('G'),   # third allele: skipped
		base_at_100('N'),   # N: skipped
	], '1\t100\t100\tA\tC\tg\tsnp')
	assert (row['UB_DEPTH'], row['UB_ALT'], row['UB_VAF']) == ('2', '1', '0.5')


def test_mnp_partial_alleles_skipped(tmp_path):
	"""issues.md #6: a read with only some of the mnp's alt bases is not ref"""
	row = run_synthetic(tmp_path, [
		base_at_100('AA'),  # ref
		base_at_100('CG'),  # alt
		base_at_100('CA'),  # first alt base only: skipped
		base_at_100('AG'),  # second alt base only: skipped
	], '1\t100\t101\tAA\tCG\tg\tmnp')
	assert (row['UB_DEPTH'], row['UB_ALT'], row['UB_VAF']) == ('2', '1', '0.5')


class FakeRef:
	"""stands in for RefWindow, over a plain string starting at 0"""
	def __init__(self, seq):
		self.seq = seq

	def base(self, i):
		return self.seq[i] if 0 <= i < len(self.seq) else None

	def bases(self, i, n):
		return self.seq[i:i + n]


@pytest.mark.parametrize('p,seq,ref_len,left,right', [
	# ref: TTAAAAATT; positions 2-6 are A
	(4, 'A', 1, (2, 'A'), (6, 'A')),         # 1-bp deletion in the run
	(6, 'A', 0, (2, 'A'), (7, 'A')),         # 1-bp insertion before base 6
	(8, 'T', 1, (7, 'T'), (8, 'T')),         # deletion after the run
	(3, 'AA', 2, (2, 'AA'), (5, 'AA')),      # 2-bp deletion
])
def test_left_right_align_homopolymer(p, seq, ref_len, left, right):
	ref = FakeRef('TTAAAAATT')
	assert left_align(p, seq, ref) == left
	assert right_align(p, seq, ref_len, ref) == right


def test_left_align_dinucleotide_repeat():
	ref = FakeRef('GGCACACAGG')  # CA repeat at 2-7
	assert left_align(6, 'CA', ref) == (2, 'CA')   # deletion of the last CA
	assert left_align(5, 'AC', ref) == (2, 'CA')   # same deletion, written shifted by one
	assert right_align(2, 'CA', 2, ref) == (6, 'CA')
	assert left_align(8, 'CA', ref) == (2, 'CA')   # insertion of CA after the repeat
	assert left_align(4, 'AC', ref) == (4, 'AC')   # no reference: unchanged
	assert left_align(4, 'AC', None) == (4, 'AC')


REPEAT_REF = 'T' * 95 + 'AAAAA' + 'T' * 900  # A run at 1-based 96-100


@pytest.mark.parametrize('use_ref,expected', [(True, ('3', '2')), (False, ('4', '1'))])
def test_deletion_in_repeat(tmp_path, use_ref, expected):
	"""issues.md #10: with a reference, a deletion anywhere in the run matches"""
	row = run_synthetic(tmp_path, [
		('15M1D24M', None),   # deletion of 96, the first A: alt
		('19M1D20M', None),   # deletion of 100, the last A: alt only with -r
		('40M', None),        # spans the run: ref
		('18M', None),        # ends at 98, inside the run: ref only without -r
		('15M2D23M', None),   # 2-bp deletion: neither
	], '1\t96\t96\tA\t-\tg\tdeletion', REPEAT_REF if use_ref else None)
	assert (row['UB_DEPTH'], row['UB_ALT']) == expected


@pytest.mark.parametrize('use_ref,expected', [(True, ('3', '2')), (False, ('3', '1'))])
def test_insertion_in_repeat(tmp_path, use_ref, expected):
	"""issues.md #10: with a reference, an insertion anywhere in the run matches"""
	reads = [
		('15M1I25M', None),   # A inserted after 95 (before the run): alt
		('20M1I20M', None),   # A inserted after 100 (end of the run): alt with -r, else ref
		('40M', None),        # ref
	]
	if not use_ref:
		reads = [(c, read_from_ref(REPEAT_REF, c)) for c, _ in reads]
	row = run_synthetic(tmp_path, reads, '1\t95\t95\t-\tA\tg\tinsertion',
						REPEAT_REF if use_ref else None)
	assert (row['UB_DEPTH'], row['UB_ALT']) == expected


def test_reference_insertion_in_repeat_real_data(tmp_path, default_run):
	"""the test insertion is in an AAA run; two soft-clipped reads put it at the start"""
	with_ref = run(tmp_path, '-r', REFERENCE)
	assert read_counts(default_run['counts_UB'], INS) == (16, 10)
	assert read_counts(with_ref['counts_UB'], INS) == (14, 12)
	# the deletion is already left-aligned in every read
	assert read_counts(with_ref['counts_UB'], DEL) == read_counts(default_run['counts_UB'], DEL)
	# snp/mnp counts don't use the reference
	for site in (SNV, SNV_MQ, NOCOV, MNP):
		assert normalized(by_site(with_ref['AllCounts'])[site]) == \
			normalized(by_site(default_run['AllCounts'])[site])


def test_reference_missing_chromosome_dies(tmp_path):
	(tmp_path / 'ref.fa').write_text('>2\nACGT\n')
	cmd = [sys.executable, SCRIPT, '-r', 'ref.fa', BAM, VARIANTS, BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True)
	assert result.returncode == 1
	assert 'no sequence for chromosome(s) 1' in result.stderr


def test_deletion_reference_mismatch_warns(tmp_path):
	variants = tmp_path / 'v.tsv'
	variants.write_text('chrm\tstart\tstop\tref\tvar\tgene\ttype\n1\t12904213\t12904213\tG\t-\tg\tdel\n')
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, BAM, str(variants), BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert "deleted bases don't match the reference (C)" in result.stderr


def test_no_reference_warns_for_indels(tmp_path):
	cmd = [sys.executable, SCRIPT, BAM, VARIANTS, BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert 'no --reference given' in result.stderr


def test_vcf_style_matches_dash(tmp_path, default_run):
	vcf = by_site(run(tmp_path, variants=VARIANTS_VCF)['AllCounts'])
	dash = by_site(default_run['AllCounts'])
	assert set(vcf) == set(VCF_SITES)
	cols = ('UB_DEPTH', 'UB_ALT', 'UB_VAF', 'CB_barcodes', 'CB:UB_barcodes')
	for vcf_site, dash_site in VCF_SITES.items():
		got, want = normalized(vcf[vcf_site]), normalized(dash[dash_site])
		assert [got[c] for c in cols] == [want[c] for c in cols], vcf_site


def test_complex_variant_is_skipped(tmp_path):
	variants = tmp_path / 'variants.tsv'
	lines = open(VARIANTS).read().splitlines()
	lines.insert(2, '1\t13625741\t13625743\tTGC\tCA\tLactb2\tcomplex')
	variants.write_text('\n'.join(lines) + '\n')
	cmd = [sys.executable, SCRIPT, BAM, str(variants), BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert 'skipping unsupported variant' in result.stderr
	rows = read_tsv(str(tmp_path / 'out_AllCounts.tsv'))
	assert [int(r['start']) for r in rows] == [DEL, INS, SNV, SNV_MQ, NOCOV, MNP]


def run_expecting_error(tmp_path, text):
	variants = tmp_path / 'variants.tsv'
	variants.write_text(text)
	cmd = [sys.executable, SCRIPT, BAM, str(variants), BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True)
	assert result.returncode == 1
	assert not list(tmp_path.glob('out_*.tsv'))  # fails before writing any output
	assert result.stderr.startswith('ERROR: ')
	assert result.stderr.strip()[len('ERROR: '):] in (tmp_path / 'out.log').read_text()
	return result.stderr


def test_trailing_blank_line_dies(tmp_path):
	"""issues.md #7"""
	err = run_expecting_error(tmp_path, open(VARIANTS).read() + '\n')
	assert 'line 8: blank line' in err


def test_malformed_line_dies(tmp_path):
	lines = open(VARIANTS).read().splitlines()
	lines[3] = lines[3].replace('\t', ' ')
	err = run_expecting_error(tmp_path, '\n'.join(lines) + '\n')
	assert 'line 4: expected at least 7 tab-separated columns, found 1' in err
	assert repr(lines[3]) in err


def test_empty_variant_file_dies(tmp_path):
	assert 'variant file is empty' in run_expecting_error(tmp_path, '')


@pytest.mark.parametrize('header', ['', '#chrm\tstart\tstop\tref\tvar\tgene\ttype\n'],
						 ids=['no_header', 'hash_header'])
def test_tsv_header_optional(tmp_path, default_run, header):
	lines = open(VARIANTS).read().splitlines(True)[1:]
	variants = tmp_path / 'variants.tsv'
	variants.write_text(header + ''.join(lines))
	result = run(tmp_path, variants=str(variants))
	assert [normalized(r) for r in result['AllCounts']] == [normalized(r) for r in default_run['AllCounts']]


def test_is_tsv_header():
	from cb_sniffer import is_tsv_header
	assert is_tsv_header(open(VARIANTS).readline())
	assert is_tsv_header('#1\t100\t100\tA\tC\tg\tt\n')
	assert not is_tsv_header('1\t100\t100\tA\tC\tg\tt\n')
	assert not is_tsv_header('\n')  # blank first line is still an error, not a header


def test_headerless_malformed_line_numbers(tmp_path):
	"""without a header, line numbers count from the first line"""
	lines = open(VARIANTS).read().splitlines()[1:]
	lines[2] = lines[2].replace('\t', ' ')
	err = run_expecting_error(tmp_path, '\n'.join(lines) + '\n')
	assert 'line 3: expected at least 7 tab-separated columns' in err


def test_blank_first_line_dies(tmp_path):
	assert 'line 1: blank line' in run_expecting_error(tmp_path, '\n' + open(VARIANTS).read())


def test_gzipped_barcodes(tmp_path, default_run):
	import gzip
	gz = str(tmp_path / 'barcodes.tsv.gz')
	with open(BARCODES, 'rb') as f, gzip.open(gz, 'wb') as out:
		out.write(f.read())
	cmd = [sys.executable, SCRIPT, BAM, VARIANTS, gz, 'out']
	subprocess.run(cmd, cwd=str(tmp_path), check=True, stdout=subprocess.DEVNULL)
	rows = read_tsv(str(tmp_path / 'out_AllCounts.tsv'))
	assert [normalized(r) for r in rows] == [normalized(r) for r in default_run['AllCounts']]


def test_output_order_is_reproducible(tmp_path):
	"""barcode lists and CB/UB row order don't depend on Python's per-process hash seed"""
	outputs = []
	for seed in ('1', '2'):
		d = tmp_path / seed
		d.mkdir()
		cmd = [sys.executable, SCRIPT, '-r', REFERENCE, BAM, VCF, BARCODES, 'out']
		subprocess.run(cmd, cwd=str(d), check=True, stdout=subprocess.DEVNULL,
					   stderr=subprocess.DEVNULL, env=dict(os.environ, PYTHONHASHSEED=seed))
		outputs.append([(d / 'out_{}.tsv'.format(n)).read_text()
						for n in ('AllCounts', 'counts_CB', 'counts_UB')])
	assert outputs[0] == outputs[1]
	ub = read_tsv(str(tmp_path / '1' / 'out_counts_UB.tsv'))
	site = [r['barcode'] for r in ub if r['start'] == str(SNV) and r['alt'] == 'C']
	assert site and site == sorted(site)


# --- CRAM input ----------------------------------------------------------------

# keep htslib from finding the reference anywhere but -r (no cache, no EBI download)
NO_REF_LOOKUP = dict(os.environ, REF_PATH='/nonexistent', REF_CACHE='/nonexistent')


@pytest.fixture(scope='module')
def cram(tmp_path_factory):
	"""test.bam as an indexed CRAM whose reference is not embedded, and whose
	header's UR: path no longer exists, so -r is the only way to decode it.
	The header is cut down to match the test reference (only contig 1, truncated);
	otherwise htslib embeds the reference in the CRAM."""
	import pysam
	d = tmp_path_factory.mktemp('cram')
	ref, path = str(d / 'gone.fa.gz'), str(d / 'test.cram')
	for ext in ('', '.fai', '.gzi'):
		shutil.copy(REFERENCE + ext, ref + ext)
	with pysam.FastaFile(REFERENCE) as fa:
		sq = [{'SN': '1', 'LN': fa.get_reference_length('1')}]
	with pysam.AlignmentFile(BAM) as src:
		header = src.header.to_dict()
		header['SQ'] = sq
		header = pysam.AlignmentHeader.from_dict(header)
		with pysam.AlignmentFile(path, 'wc', header=header, reference_filename=ref,
								 format_options=[b'embed_ref=0']) as out:
			for rec in src:
				fields = rec.to_string().split('\t')
				if fields[6] not in ('=', '*', '1'):  # mate on a contig not in the header
					fields[6:8] = ['*', '0']
				out.write(pysam.AlignedSegment.fromstring('\t'.join(fields), header))
	pysam.index(path)
	for ext in ('', '.fai', '.gzi'):
		os.remove(ref + ext)
	return path


def test_cram_requires_reference(tmp_path, cram):
	cmd = [sys.executable, SCRIPT, cram, VARIANTS, BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True, env=NO_REF_LOOKUP)
	assert result.returncode == 1
	assert 'CRAM input requires a reference FASTA' in result.stderr
	assert not list(tmp_path.glob('out_*.tsv'))


def test_cram_matches_bam(tmp_path, cram):
	from_bam = run(tmp_path, '-r', REFERENCE, prefix='bam')
	from_cram = run(tmp_path, '-r', REFERENCE, prefix='cram', bam=cram, env=NO_REF_LOOKUP)
	assert [normalized(r) for r in from_cram['AllCounts']] == [normalized(r) for r in from_bam['AllCounts']]
	for name in ('counts_CB', 'counts_UB'):
		key = lambda r: (r['start'], r['alt'], r['barcode'])
		assert sorted(from_cram[name], key=key) == sorted(from_bam[name], key=key)


def test_cram_no_multiple_iterators_warning(tmp_path, cram):
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, cram, VARIANTS, BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True, env=NO_REF_LOOKUP)
	assert result.returncode == 0
	assert 'multiple_iterators' not in result.stderr


# --- VCF input -----------------------------------------------------------------

@pytest.fixture(scope='module')
def vcf_copies(tmp_path_factory):
	"""variants.vcf as bgzipped VCF and as BCF (neither indexed)"""
	import pysam
	d = tmp_path_factory.mktemp('vcf')
	gz, bcf = str(d / 'v.vcf.gz'), str(d / 'v.bcf')
	pysam.tabix_compress(VCF, gz)
	with pysam.VariantFile(VCF) as src, pysam.VariantFile(bcf, 'wb', header=src.header) as out:
		for rec in src:
			out.write(rec)
	return {'vcf.gz': gz, 'bcf': bcf}


def test_is_vcf(tmp_path, vcf_copies):
	import gzip
	gz_tsv = str(tmp_path / 'variants.tsv.gz')
	with open(VARIANTS, 'rb') as f, gzip.open(gz_tsv, 'wb') as out:
		out.write(f.read())
	assert is_vcf(VCF) and is_vcf(vcf_copies['vcf.gz']) and is_vcf(vcf_copies['bcf'])
	assert not is_vcf(VARIANTS) and not is_vcf(gz_tsv)


def test_read_variants_tsv():
	assert [(v.start, v.ref, v.alt, v.gene, v.event) for v in read_variants(VARIANTS)][:2] == [
		(DEL, 'C', '-', 'Slco5a1', 'deletion'), (INS, '-', 'A', 'Tram1', 'insertion')]


def test_read_variants_vcf(capsys):
	got = [(v.start, v.end, v.ref, v.alt, v.gene, v.event) for v in read_variants(VCF)]
	assert got == [
		(10000001, 10000001, 'T', 'G', '.', 'snp'),           # no ID, FILTER '.'
		(10066352, 10066356, 'CTATG', 'CCATA', 'Cspp1', 'mnp'),
		(12904212, 12904216, 'ACCCC', 'ACCC', 'Slco5a1', 'del'),
		(13562395, 13562395, 'A', 'AA', 'Tram1', 'ins'),
		(13625741, 13625741, 'T', 'C', 'Lactb2', 'snp'),      # multi-allelic: one per ALT
		(13625741, 13625741, 'T', 'G', 'Lactb2', 'snp'),
		(15268954, 15268954, 'T', 'C', 'Gm10566', 'snp'),     # FILTER LowQual
	]
	assert 'skipping symbolic/breakend allele' in capsys.readouterr().err  # <DEL>


def test_read_variants_vcf_pass_only():
	assert 15268954 not in [v.start for v in read_variants(VCF, pass_only=True)]
	assert len(list(read_variants(VCF, pass_only=True))) == 6


@pytest.fixture(scope='module')
def vcf_run(tmp_path_factory):
	return run(tmp_path_factory.mktemp('vcf_run'), variants=VCF)['AllCounts']


def test_vcf_matches_tsv(vcf_run, default_run):
	vcf = {(int(r['start']), r['alt']): normalized(r) for r in vcf_run}
	tsv = {s: normalized(r) for s, r in by_site(default_run['AllCounts']).items()}
	assert set(vcf) == set(VCF_FILE_SITES) | {(13625741, 'G')}
	for key, site in VCF_FILE_SITES.items():
		assert [vcf[key][c] for c in COUNT_COLS] == [tsv[site][c] for c in COUNT_COLS], key
	# second ALT allele at the SNV: no reads carry G, and reads carrying C
	# are neither ref nor alt for it, so only the 14 ref UMIs are counted
	assert (vcf[(13625741, 'G')]['UB_DEPTH'], vcf[(13625741, 'G')]['UB_ALT']) == ('14', '0')


@pytest.mark.parametrize('fmt', ['vcf.gz', 'bcf'])
def test_compressed_vcf_matches_plain(tmp_path, vcf_copies, vcf_run, fmt):
	cmd = [sys.executable, SCRIPT, BAM, vcf_copies[fmt], BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert 'idx_find_and_load' not in result.stderr  # no htslib noise about the missing index
	rows = read_tsv(str(tmp_path / 'out_AllCounts.tsv'))
	assert [normalized(r) for r in rows] == [normalized(r) for r in vcf_run]


def test_pass_only_option(tmp_path):
	rows = run(tmp_path, '--pass-only', variants=VCF)['AllCounts']
	assert len(rows) == 6
	assert 15268954 not in [int(r['start']) for r in rows]
