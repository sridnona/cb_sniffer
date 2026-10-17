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
BARCODES = os.path.join(TEST_DIR, 'barcodes.tsv')

sys.path.insert(0, REPO_DIR)
from cb_sniffer import GenomicPosition, is_vcf, read_variants  # noqa: E402

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


def run(tmp_path, *opts, variants=VARIANTS, prefix='out'):
	"""Run cb_sniffer.py and return its three output tables, keyed by output name."""
	cmd = [sys.executable, SCRIPT, *opts, BAM, variants, BARCODES, prefix]
	subprocess.run(cmd, cwd=str(tmp_path), check=True, stdout=subprocess.DEVNULL)
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


def test_too_few_columns_exits():
	with pytest.raises(SystemExit):
		GenomicPosition('1\t100\t100\tA\tT\n')


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


def run_synthetic(tmp_path, reads, variant_line):
	"""
	Write reads to a BAM (each its own cell and UMI, all starting at 1-based
	position 81 on contig 1), run cb_sniffer on one variant and return its
	AllCounts row. reads: list of (cigar, sequence or None for all A's).
	"""
	import pysam
	bam = str(tmp_path / 'synthetic.bam')
	header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 1000}]}
	with pysam.AlignmentFile(bam, 'wb', header=header) as out:
		for i, (cigar, seq) in enumerate(reads):
			a = pysam.AlignedSegment()
			a.query_name, a.cigarstring = 'read{}'.format(i), cigar
			a.query_sequence = seq or 'A' * a.infer_query_length()
			a.reference_id, a.reference_start, a.mapping_quality = 0, 80, 255
			a.query_qualities = pysam.qualitystring_to_array('F' * len(a.query_sequence))
			a.set_tags([('CB', 'CELL{}-1'.format(i)), ('UB', 'UMI{}'.format(i))])
			out.write(a)
	pysam.index(bam)
	(tmp_path / 'bc.tsv').write_text(''.join('CELL{}-1\n'.format(i) for i in range(len(reads))))
	(tmp_path / 'v.tsv').write_text('chrm\tstart\tstop\tref\tvar\tgene\ttype\n' + variant_line + '\n')
	subprocess.run([sys.executable, SCRIPT, bam, 'v.tsv', 'bc.tsv', 'out'],
				   cwd=str(tmp_path), check=True, stdout=subprocess.DEVNULL)
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
		('20M1D20M', None),   # deletion at 101, base 100 present: ref
	], variant_line)
	assert (row['UB_DEPTH'], row['UB_ALT']) == ('3', '1')
	assert row['CB_barcodes'] == 'CELL0-1'


@pytest.mark.parametrize('cigar,expected', [
	('19M1D20M', 1), ('19M2D20M', 2), ('18M2D20M', None), ('40M', None),
	('20M1D20M', None), ('5S19M1D20M', 1), ('10M5N4M1D20M', 1),
])
def test_deletion_length(cigar, expected):
	read = make_read(cigar, start=81)
	assert variant('A', '-').deletion_length(read) == expected


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


@pytest.mark.xfail(strict=True, reason='issues.md #7: a blank line in the variant file aborts the run')
def test_trailing_blank_line(tmp_path):
	variants = tmp_path / 'variants.tsv'
	shutil.copy(VARIANTS, str(variants))
	with open(str(variants), 'a') as f:
		f.write('\n')
	assert len(run(tmp_path, variants=str(variants))['AllCounts']) == 6


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
	# second ALT allele at the SNV: no reads carry G
	assert (vcf[(13625741, 'G')]['UB_DEPTH'], vcf[(13625741, 'G')]['UB_ALT']) == ('28', '0')


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
