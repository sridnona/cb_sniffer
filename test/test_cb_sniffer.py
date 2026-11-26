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

import random

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
from cb_sniffer import (GenomicPosition, is_vcf, read_variants, left_align, right_align,  # noqa: E402
						local_score, Haplotypes, RefWindow)

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
	"""Run cb_sniffer.py (with -r REFERENCE) and return its three output tables, keyed by output name."""
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, *opts, bam, variants, BARCODES, prefix]
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


def run_synthetic(tmp_path, reads, variant_line, refseq='A' * 1000):
	"""
	Write reads to a BAM (each its own cell and UMI, all starting at 1-based
	position 81 on contig 1), run cb_sniffer on one variant and return its
	AllCounts row. reads: list of (cigar, sequence or None). A None sequence
	follows refseq (all A by default), which is also passed as -r.
	"""
	import pysam
	bam = str(tmp_path / 'synthetic.bam')
	header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 1000}]}
	with pysam.AlignmentFile(bam, 'wb', header=header) as out:
		for i, (cigar, seq) in enumerate(reads):
			a = pysam.AlignedSegment()
			a.query_name, a.cigarstring = 'read{}'.format(i), cigar
			if seq is None:
				seq = read_from_ref(refseq, cigar)
			a.query_sequence = seq
			a.reference_id, a.reference_start, a.mapping_quality = 0, 80, 255
			a.query_qualities = pysam.qualitystring_to_array('F' * len(a.query_sequence))
			a.set_tags([('CB', 'CELL{}-1'.format(i)), ('UB', 'UMI{}'.format(i))])
			out.write(a)
	pysam.index(bam)
	(tmp_path / 'bc.tsv').write_text(''.join('CELL{}-1\n'.format(i) for i in range(len(reads))))
	(tmp_path / 'v.tsv').write_text('chrm\tstart\tstop\tref\tvar\tgene\ttype\n' + variant_line + '\n')
	(tmp_path / 'ref.fa').write_text('>1\n' + refseq + '\n')
	pysam.faidx(str(tmp_path / 'ref.fa'))
	subprocess.run([sys.executable, SCRIPT, '-r', 'ref.fa', bam, 'v.tsv', 'bc.tsv', 'out'],
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
	# 1-based 99-100 is AA, between C and G; deleting either A is the same deletion
	row = run_synthetic(tmp_path, [
		('19M1D20M', None),   # 1-bp deletion of base 100: alt
		('19M2D20M', None),   # 2-bp deletion starting at 100: skipped
		('18M2D20M', None),   # deletion starting at 99, covers 100: skipped
		('40M', None),        # no deletion: ref
		('20M1D20M', None),   # deletion of the G at 101: skipped
	], variant_line, MIXED_REF)
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('2', '1')
	assert row['CB_barcodes'] == 'CELL0-1'


_rng = random.Random(1)
MIXED_REF = ''.join(_rng.choice('ACGT') for _ in range(1000))
MIXED_REF = MIXED_REF[:97] + 'CAAG' + MIXED_REF[101:]  # 1-based 98-101


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
	assert left_align(4, 'AC', ref) == (4, 'AC')   # base 3 is A, not C: can't shift


REPEAT_REF = 'T' * 95 + 'AAAAA' + 'T' * 900  # A run at 1-based 96-100


def test_deletion_in_repeat(tmp_path):
	"""issues.md #10: a deletion anywhere in the run matches"""
	row = run_synthetic(tmp_path, [
		('15M1D24M', None),   # deletion of 96, the first A: alt
		('19M1D20M', None),   # deletion of 100, the last A: alt
		('40M', None),        # spans the run: ref
		('18M', None),        # ends at 98, inside the run: neither
		('15M2D23M', None),   # 2-bp deletion: neither
	], '1\t96\t96\tA\t-\tg\tdeletion', REPEAT_REF)
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('3', '2')


def test_insertion_in_repeat(tmp_path):
	"""issues.md #10: an insertion anywhere in the run matches"""
	row = run_synthetic(tmp_path, [
		('15M1I25M', None),   # A inserted after 95 (before the run): alt
		('20M1I20M', None),   # A inserted after 100 (end of the run): alt
		('40M', None),        # ref
	], '1\t95\t95\t-\tA\tg\tinsertion', REPEAT_REF)
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('3', '2')


def test_insertion_in_repeat_real_data(default_run):
	"""the test insertion is in an AAA run, and two soft-clipped reads put it
	at the start; they're alt too (only 16 / 10 reads without left-aligning).
	One more read is aligned to the end of the run and then soft-clipped, with
	the inserted A in the clip; realigning it makes it alt"""
	assert read_counts(default_run['counts_UB'], INS) == (14, 13)


# --- indels near soft clips, splices and mismatches --------------------------------
# A read's CIGAR is used only if it's aligned without a clip, splice or mismatch
# across the repeat the indel is in (and FLANK bases further where it's clipped or
# spliced). Other reads are realigned to the ref and alt haplotypes


def random_ref(seed, fixed):
	"""1,000 random bases, with fixed: {0-based index: base or bases}"""
	rng = random.Random(seed)
	ref = [rng.choice('ACGT') for _ in range(1000)]
	for i, bases in fixed.items():
		ref[i:i + len(bases)] = bases
	return ''.join(ref)


# 6-bp insertion after 1-based 100, not in a repeat
INS6 = 'TTGACA'
INS6_REF = random_ref(2, {99: 'C', 100: 'G'})
INS6_ALT = INS6_REF[:100] + INS6 + INS6_REF[100:]
INS6_LINE = '1\t100\t100\t-\t' + INS6 + '\tg\tinsertion'


@pytest.mark.parametrize('offset', range(9))
def test_insertion_soft_clipped_at_any_offset(tmp_path, offset):
	"""a read aligned up to the insertion, or partway into it, then soft-clipped
	is alt; a ref read clipped at the same place is ref"""
	cigar = '{}M{}S'.format(20 + offset, 26 - offset)
	row = run_synthetic(tmp_path, [(cigar, INS6_ALT[80:126]), (cigar, INS6_REF[80:126])],
						INS6_LINE, INS6_REF)
	assert (row['UB_DEPTH'], row['UB_ALT'], row['CB_barcodes']) == ('2', '1', 'CELL0-1')


def test_insertion_soft_clipped_on_left(tmp_path):
	"""the inserted bases are in a clip before the aligned part of the read"""
	ref = random_ref(2, {79: 'C', 80: 'G'})
	alt = ref[:80] + INS6 + ref[80:]
	row = run_synthetic(tmp_path, [
		('26S30M', alt[60:116]),   # clip holds 20 ref bases and the insertion: alt
		('20S30M', ref[60:110]),   # clip holds ref bases: ref
	], '1\t80\t80\t-\t' + INS6 + '\tg\tinsertion', ref)
	assert (row['UB_DEPTH'], row['UB_ALT'], row['CB_barcodes']) == ('2', '1', 'CELL0-1')


def test_tandem_duplication_6bp(tmp_path):
	"""a 6-bp tandem duplication: the aligner can put the inserted bases anywhere
	in the two copies, or clip them"""
	ref = random_ref(3, {})
	alt = ref[:100] + ref[94:100] + ref[100:]
	row = run_synthetic(tmp_path, [
		('20M6I20M', alt[80:126]),   # alt
		('14M6I26M', alt[80:126]),   # alt, inserted at the start of the first copy
		('20M26S', alt[80:126]),     # alt, clipped at the junction
		('26M20S', alt[80:126]),     # alt, second copy aligned to the reference
		('20M', ref[80:100]),        # ends before the junction: not counted
		('46M', ref[80:126]),        # ref
	], '1\t100\t100\t-\t' + ref[94:100] + '\tg\tinsertion', ref)
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('5', '4')
	assert row['CB_barcodes'] == 'CELL0-1,CELL1-1,CELL2-1,CELL3-1'


def test_tandem_duplication_40bp(tmp_path):
	"""a 40-bp tandem duplication of 1-based 101-140"""
	ref = random_ref(4, {99: 'A', 100: 'C', 139: 'G', 140: 'T'})
	alt = ref[:140] + ref[100:140] + ref[140:]
	row = run_synthetic(tmp_path, [
		('60M40I20M', alt[80:200]),   # alt
		('60M60S', alt[80:200]),      # alt, clipped at the junction
		('100M20S', alt[80:200]),     # alt, second copy aligned to the reference
		('50M', ref[80:130]),         # ends in the first copy: not counted
		('120M', ref[80:200]),        # ref
	], '1\t140\t140\t-\t' + ref[100:140] + '\tg\tinsertion', ref)
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('4', '3')
	assert row['CB_barcodes'] == 'CELL0-1,CELL1-1,CELL2-1'


def test_str_insertion(tmp_path):
	"""a CA inserted in (CA)6 at 1-based 91-102: reads that span the repeat are
	called, wherever the aligner puts the insertion or clips"""
	ref = random_ref(5, {89: 'G', 90: 'CA' * 6, 102: 'T'})
	alt = ref[:102] + 'CA' + ref[102:]
	row = run_synthetic(tmp_path, [
		('22M2I20M', alt[80:124]),   # alt
		('12M2I30M', alt[80:124]),   # alt, inserted partway through the repeat
		('22M22S', alt[80:124]),     # alt, clipped at the end of the repeat
		('16M28S', alt[80:124]),     # alt, clipped partway through the repeat
		('16M', ref[80:96]),         # ends in the repeat: not counted
		('44M', ref[80:124]),        # ref
	], '1\t102\t102\t-\tCA\tg\tinsertion', ref)
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('5', '4')
	assert row['CB_barcodes'] == 'CELL0-1,CELL1-1,CELL2-1,CELL3-1'


def test_spliced_reads_near_indel(tmp_path):
	row = run_synthetic(tmp_path, [
		('15M100N25M', None),                               # spliced across the site: skipped
		('20M100N20M', None),                               # spliced at the site: skipped
		('30M100N10M', None),                               # spliced 10 bases after it: ref
		('20M6I2M100N12M', INS6_ALT[80:108] + INS6_REF[202:214]),  # alt, spliced 2 bases after it
	], INS6_LINE, INS6_REF)
	assert (row['UB_DEPTH'], row['UB_ALT'], row['CB_barcodes']) == ('2', '1', 'CELL3-1')


def test_deletion_aligned_as_mismatches(tmp_path):
	"""a deletion near the end of a read can be aligned as mismatches instead;
	the mismatches make it be realigned, and it's alt"""
	ref = random_ref(6, {113: 'ACT'})
	row = run_synthetic(tmp_path, [
		('40M', ref[80:114] + ref[115:121]),   # base 115 deleted: alt
		('40M', None),                         # ref
	], '1\t115\t115\tC\t-\tg\tdeletion', ref)
	assert (row['UB_DEPTH'], row['UB_ALT'], row['CB_barcodes']) == ('2', '1', 'CELL0-1')


def test_local_score():
	assert local_score('ACGT', 'ACGT') == 4
	assert local_score('ACGT', 'TTACGTTT') == 4
	assert local_score('ANGT', 'ACGT') == 3              # N scores 0
	assert local_score('A' * 10 + 'C' * 10, 'A' * 10 + 'GG' + 'C' * 10) == 20 - 8  # 2-bp gap
	assert local_score('A' * 10 + 'G' + 'A' * 10, 'A' * 21) == 20 - 4  # mismatch
	assert local_score('', 'ACGT') == 0


def test_haplotypes():
	ref = RefWindow.__new__(RefWindow)
	ref.start, ref.seq = 0, 'ACGTACGTAC' * 20
	ref.end = len(ref.seq)
	h = Haplotypes(ref, (100, 'TT', 0), 99, 100)   # TT inserted before base 100
	assert h.alt == h.ref[:100 - h.start] + 'TT' + h.ref[100 - h.start:]
	assert h.alt_index(99) == 99 - h.start and h.alt_index(100) == 102 - h.start
	h = Haplotypes(ref, (100, 'AC', 2), 99, 102)   # 100-101 deleted
	assert h.alt == h.ref[:100 - h.start] + h.ref[102 - h.start:]
	assert h.alt_index(101) == 100 - h.start and h.alt_index(102) == 100 - h.start


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
	assert "ref allele doesn't match the reference (C), counting anyway" in result.stderr


def test_snv_reference_mismatch_warns(tmp_path):
	variants = tmp_path / 'v.tsv'
	variants.write_text('1\t13625741\t13625741\tG\tC\tg\tsnp\n')  # ref base is T
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, BAM, str(variants), BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert "ref allele doesn't match the reference (T)" in result.stderr
	assert read_tsv(str(tmp_path / 'out_AllCounts.tsv'))  # still counted


def test_reference_required(tmp_path):
	cmd = [sys.executable, SCRIPT, BAM, VARIANTS, BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True)
	assert result.returncode == 2
	assert 'required: -r/--reference' in result.stderr
	assert not list(tmp_path.glob('out_*.tsv'))


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
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, BAM, str(variants), BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert 'skipping unsupported variant' in result.stderr
	rows = read_tsv(str(tmp_path / 'out_AllCounts.tsv'))
	assert [int(r['start']) for r in rows] == [DEL, INS, SNV, SNV_MQ, NOCOV, MNP]


def run_expecting_error(tmp_path, text):
	variants = tmp_path / 'variants.tsv'
	variants.write_text(text)
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, BAM, str(variants), BARCODES, 'out']
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
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, BAM, VARIANTS, gz, 'out']
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


def test_cram_matches_bam(tmp_path, cram):
	from_bam = run(tmp_path, prefix='bam')
	from_cram = run(tmp_path, prefix='cram', bam=cram, env=NO_REF_LOOKUP)
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
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, BAM, vcf_copies[fmt], BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert 'idx_find_and_load' not in result.stderr  # no htslib noise about the missing index
	rows = read_tsv(str(tmp_path / 'out_AllCounts.tsv'))
	assert [normalized(r) for r in rows] == [normalized(r) for r in vcf_run]


def test_pass_only_option(tmp_path):
	rows = run(tmp_path, '--pass-only', variants=VCF)['AllCounts']
	assert len(rows) == 6
	assert 15268954 not in [int(r['start']) for r in rows]


# --- scATAC (--atac) -----------------------------------------------------------

def atac_read(name, cell, alt=False, flag=0, cigar='40M', ub=None):
	"""
	(name, cell, cigar, sequence, flag, ub) for write_atac_bam(): a 40-base read
	starting at 81 on contig 1, all A except base 100 (C when alt).
	For cigar 19M1D20M the sequence is all A (base 100 deleted).
	"""
	if 'D' in cigar:
		seq = 'A' * 39
	else:
		seq = 'A' * 19 + ('C' if alt else 'A') + 'A' * 20
	return name, cell, cigar, seq, flag, ub


def mates(name, cell, alt1=False, alt2=False, cigar='40M'):
	"""both reads of a pair, overlapping over their whole length"""
	return [atac_read(name, cell, alt1, 0x1 | 0x2 | 0x40, cigar),
			atac_read(name, cell, alt2, 0x1 | 0x2 | 0x80, cigar)]


def write_atac_bam(tmp_path, reads, cells=('CELL1-1', 'CELL2-1', 'CELL3-1')):
	"""write reads (see atac_read()) to an indexed BAM, plus a barcode file"""
	import pysam
	bam = str(tmp_path / 'atac.bam')
	header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 1000}]}
	with pysam.AlignmentFile(bam, 'wb', header=header) as out:
		for name, cell, cigar, seq, flag, ub in reads:
			a = pysam.AlignedSegment()
			a.query_name, a.cigarstring, a.query_sequence, a.flag = name, cigar, seq, flag
			a.reference_id, a.reference_start, a.mapping_quality = 0, 80, 60
			if flag & 0x1:
				a.next_reference_id, a.next_reference_start = 0, 80
			a.query_qualities = pysam.qualitystring_to_array('F' * len(seq))
			tags = [('CB', cell)] + ([('UB', ub)] if ub else [])
			a.set_tags(tags)
			out.write(a)
	pysam.index(bam)
	(tmp_path / 'bc.tsv').write_text(''.join(c + '\n' for c in cells))
	return bam


def run_atac(tmp_path, reads, variant_line='1\t100\t100\tA\tC\tg\tsnp', opts=('--atac',),
			 refseq='A' * 1000):
	"""run on a synthetic ATAC BAM; return (CompletedProcess, AllCounts rows or None)"""
	import pysam
	bam = write_atac_bam(tmp_path, reads)
	(tmp_path / 'v.tsv').write_text(variant_line + '\n')
	(tmp_path / 'ref.fa').write_text('>1\n' + refseq + '\n')
	pysam.faidx(str(tmp_path / 'ref.fa'))
	result = subprocess.run([sys.executable, SCRIPT, '-r', 'ref.fa', *opts, bam, 'v.tsv', 'bc.tsv', 'out'],
							cwd=str(tmp_path), capture_output=True, text=True)
	out = tmp_path / 'out_AllCounts.tsv'
	return result, read_tsv(str(out)) if out.exists() else None


ATAC_SNV_READS = (
	mates('fragA', 'CELL1-1', True, True)      # both mates alt: one alt fragment
	+ [atac_read('fragB', 'CELL1-1', True)]    # unpaired alt read
	+ mates('fragC', 'CELL2-1', False, False)  # both mates ref: one ref fragment
	+ mates('fragD', 'CELL2-1', False, True)   # mates disagree: dropped
	+ [atac_read('dup', 'CELL3-1', True, 0x400),     # duplicate: skipped
	   atac_read('fragB', 'CELL1-1', False, 0x800),  # supplementary of fragB: its primary wins
	   atac_read('supp', 'CELL3-1', True, 0x800),    # supplementary, primary elsewhere: counted
	   atac_read('sec', 'CELL3-1', True, 0x100),     # secondary: skipped
	   atac_read('qcfail', 'CELL3-1', True, 0x200)]  # QC fail: skipped
)


def test_no_umis_without_atac_dies(tmp_path):
	"""issue #17: an ATAC BAM used to give empty output with no error"""
	result, rows = run_atac(tmp_path, ATAC_SNV_READS, opts=())
	assert result.returncode == 1
	assert 'have a UB (UMI) tag' in result.stderr and '--atac' in result.stderr
	assert rows is None and not list(tmp_path.glob('out_*.tsv'))


def test_atac_counts_fragments(tmp_path):
	result, rows = run_atac(tmp_path, ATAC_SNV_READS)
	assert result.returncode == 0, result.stderr
	row = rows[0]
	assert (row['UB_DEPTH'], row['UB_ALT'], row['UB_VAF']) == ('4', '3', '0.75')
	assert row['CB_barcodes'] == 'CELL1-1,CELL3-1'
	assert row['CB:UB_barcodes'] == 'CELL1-1:fragA,CELL1-1:fragB,CELL3-1:supp'
	ub = read_tsv(str(tmp_path / 'out_counts_UB.tsv'))
	# read counts per fragment: both mates are counted, the fragment once;
	# fragD (mates disagree) fails consensus and is left out, like a discordant UMI.
	# fragB's ref supplementary isn't counted, so fragB isn't dropped as discordant
	assert {r['barcode']: (r['ref_count'], r['alt_count']) for r in ub} == {
		'CELL1-1:fragA': ('0', '2'), 'CELL1-1:fragB': ('0', '1'),
		'CELL2-1:fragC': ('2', '0'), 'CELL3-1:supp': ('0', '1')}
	cb = read_tsv(str(tmp_path / 'out_counts_CB.tsv'))
	assert {r['barcode']: (r['ref_count'], r['alt_count']) for r in cb} == {
		'CELL1-1': ('0', '2'), 'CELL2-1': ('1', '0'), 'CELL3-1': ('0', '1')}


def test_read_once():
	import pysam
	from cb_sniffer import read_once

	def aln(name, flag):
		a = pysam.AlignedSegment()
		a.query_name, a.flag = name, flag
		return a

	r1, r2 = 0x1 | 0x40, 0x1 | 0x80
	calls = [
		(aln('a', 0x800), 'ref'),        # supplementary first, then its primary: primary wins
		(aln('a', 0), 'alt'),
		(aln('b', r1), 'alt'),           # mates are different reads: both kept
		(aln('b', r2 | 0x800), 'ref'),
		(aln('c', 0x100), 'ref'),        # secondary loses to supplementary
		(aln('c', 0x800), 'alt'),
		(aln('d', 0x800), 'alt'),        # supplementary only: kept
	]
	assert [(a.query_name, a.flag, x) for a, x in read_once(calls)] == [
		('a', 0, 'alt'), ('b', r1, 'alt'), ('b', r2 | 0x800, 'ref'),
		('c', 0x800, 'alt'), ('d', 0x800, 'alt')]


# all A except a C at 1-based 90 and a G at 110, so a 40-base read at 81-120
# spans the A run with 9 bases to spare either side
ATAC_DEL_REF = 'A' * 89 + 'C' + 'A' * 19 + 'G' + 'A' * 890


def atac_del_read(name, cell, flag, deleted):
	"""a read of ATAC_DEL_REF at 81-120, with base 100 deleted if deleted"""
	seq = ATAC_DEL_REF[80:120]
	if deleted:
		return name, cell, '19M1D20M', seq[:19] + seq[20:], flag, None
	return name, cell, '40M', seq, flag, None


def test_atac_deletion_counts_fragments(tmp_path):
	pair = 0x1 | 0x2
	reads = [atac_del_read('fragA', 'CELL1-1', pair | 0x40, True),   # both mates carry the deletion
			 atac_del_read('fragA', 'CELL1-1', pair | 0x80, True),
			 atac_del_read('fragB', 'CELL2-1', pair | 0x40, False),  # ref
			 atac_del_read('fragB', 'CELL2-1', pair | 0x80, False),
			 atac_del_read('dup', 'CELL3-1', 0x400, True)]
	result, rows = run_atac(tmp_path, reads, '1\t100\t100\tA\t-\tg\tdeletion', refseq=ATAC_DEL_REF)
	assert result.returncode == 0, result.stderr
	assert (rows[0]['UB_DEPTH'], rows[0]['UB_ALT']) == ('2', '1')
	assert rows[0]['CB:UB_barcodes'] == 'CELL1-1:fragA'


def test_atac_with_umis_warns(tmp_path, default_run):
	"""--atac on a scRNA BAM warns; every read with a CB is counted as its own fragment"""
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, '--atac', BAM, VARIANTS, BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert '--atac given, but reads have UB (UMI) tags' in result.stderr


def test_rna_bam_passes_umi_check(tmp_path):
	cmd = [sys.executable, SCRIPT, '-r', REFERENCE, BAM, VARIANTS, BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), check=True, capture_output=True, text=True)
	assert 'UB (UMI) tag' not in result.stderr


# --- scATAC (--atac), real data ---------------------------------------------------
# atac.bam is a subset of the 10x PBMC 5k NextGEM scATAC BAM (Cell Ranger ATAC,
# GRCh38): reads overlapping the five sites in atac_variants.tsv. See README.md

ATAC_BAM = os.path.join(TEST_DIR, 'atac.bam')
ATAC_VARIANTS = os.path.join(TEST_DIR, 'atac_variants.tsv')
ATAC_BARCODES = os.path.join(TEST_DIR, 'atac_barcodes.tsv')
ATAC_REFERENCE = os.path.join(TEST_DIR, 'atac_reference.fa.gz')  # N except +-3 kb around each site
A_GAPDH, A_ACTB, A_HLAC, A_DEL, A_INS = 6534761, 5530212, 31271428, 6470901, 5532355
# fragment GCTTGC...-1:...34820 at A_INS: mate 2 has the insertion; mate 1 is aligned one
# base past it, then soft-clipped
A_SOFTCLIP_FRAG = 'GCTTGCTCAATGCACT-1:A00519:269:H7FM2DRXX:1:2123:32353:34820'


def run_real_atac(tmp_path, *opts, check=True):
	"""run on atac.bam; return (CompletedProcess, output tables keyed by name, or None)"""
	cmd = [sys.executable, SCRIPT, '-r', ATAC_REFERENCE, *opts, ATAC_BAM, ATAC_VARIANTS, ATAC_BARCODES, 'out']
	result = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True, check=check)
	if not (tmp_path / 'out_AllCounts.tsv').exists():
		return result, None
	return result, {name: read_tsv(str(tmp_path / 'out_{}.tsv'.format(name)))
					for name in ('AllCounts', 'counts_CB', 'counts_UB')}


@pytest.fixture(scope='module')
def atac_run(tmp_path_factory):
	return run_real_atac(tmp_path_factory.mktemp('atac'), '--atac')


def depth_alt(tables):
	return {int(r['start']): (int(r['UB_DEPTH']), int(r['UB_ALT']))
			for r in tables['AllCounts']}


def ub_rows(tables, start):
	return {r['barcode']: (int(r['ref_count']), int(r['alt_count']))
			for r in tables['counts_UB'] if int(r['start']) == start}


def test_real_atac_without_flag_dies(tmp_path):
	result, tables = run_real_atac(tmp_path, check=False)
	assert result.returncode == 1
	assert 'rerun with --atac' in result.stderr
	assert tables is None


def test_real_atac_no_umi_warning(atac_run):
	result, _ = atac_run
	assert 'UB (UMI) tag' not in result.stderr


# fragments (ref + alt) / alt fragments, checked against an independent pysam count
# of distinct non-duplicate read names per allele with a listed CB
# Every fragment with reads at the insertion is alt once soft-clipped reads are
# realigned. Two deletion fragments are lost compared with counting from the
# CIGAR alone: each has one read with a mismatch next to the deletion and only a
# base or two of evidence, so it isn't called
REAL_ATAC_EXPECTED = {
	(): {A_GAPDH: (240, 121), A_ACTB: (132, 68), A_HLAC: (320, 176),
		 A_DEL: (447, 182), A_INS: (10, 10)},
	# HLA-C has 46 alt fragments with MAPQ < 30
	('-mq', '30'): {A_GAPDH: (240, 121), A_ACTB: (132, 68), A_HLAC: (274, 130),
					A_DEL: (447, 182), A_INS: (10, 10)},
	('-bq', '30'): {A_GAPDH: (228, 119), A_ACTB: (126, 62), A_HLAC: (298, 156),
					A_DEL: (447, 182), A_INS: (10, 10)},
}


@pytest.mark.parametrize('opts', list(REAL_ATAC_EXPECTED), ids=lambda o: ' '.join(o) or 'default')
def test_real_atac_counts(tmp_path, atac_run, opts):
	tables = atac_run[1] if not opts else run_real_atac(tmp_path, '--atac', *opts)[1]
	assert depth_alt(tables) == REAL_ATAC_EXPECTED[opts]


def test_real_atac_one_count_per_fragment(atac_run):
	"""barcodes are CB:read name, and each fragment is one count in UB_DEPTH"""
	_, tables = atac_run
	listed = {l.strip() for l in open(ATAC_BARCODES)}
	for start, (depth, alt) in REAL_ATAC_EXPECTED[()].items():
		ub = ub_rows(tables, start)
		assert len(ub) == depth
		assert sum(1 for r, a in ub.values() if a) == alt
		for barcode in ub:
			cb, name = barcode.split(':', 1)
			assert cb in listed and name.startswith('A00519:')
	# where mates overlap the site both are counted in the UB file, but the
	# fragment is still one count
	gapdh = ub_rows(tables, A_GAPDH)
	assert sum(1 for r, a in gapdh.values() if r + a == 2) == 61
	assert max(r + a for r, a in gapdh.values()) == 2


def test_real_atac_skips_duplicates(atac_run):
	"""Cell Ranger ATAC marks (not removes) duplicates; none of them may be counted"""
	pysam = pytest.importorskip('pysam')
	with pysam.AlignmentFile(ATAC_BAM) as bam:
		dups = {r.get_tag('CB') + ':' + r.query_name for r in bam
				if r.is_duplicate and r.has_tag('CB')}
	assert len(dups) > 1000
	_, tables = atac_run
	counted = {r['barcode'] for r in tables['counts_UB']}
	assert counted and not counted & dups


def test_real_atac_softclip_mates(atac_run):
	"""mate 1 of A_SOFTCLIP_FRAG covers the base after the insertion, but its
	soft clip holds the inserted bases. It's realigned and called alt, so it
	agrees with mate 2 and the fragment counts (it used to be called ref, and the
	fragment was dropped)"""
	assert ub_rows(atac_run[1], A_INS)[A_SOFTCLIP_FRAG] == (0, 2)
