# Test data

A tiny subset of a mouse (mm10) 10x Genomics scRNA-seq BAM (Cell Ranger 9.0.0, C57BL/6 bone marrow).
It contains only reads overlapping the six sites in `variants.tsv` (1,476 reads, ~100 KB).

| File | Contents |
|------|----------|
| `test.bam`, `test.bam.bai` | Reads overlapping the test sites on chr1. `@PG` header lines were removed. |
| `variants.tsv` | Six sites (see below), indels in dash notation |
| `variants_vcf.tsv` | Four of the same sites written VCF-style, with padding bases. Each must give the same counts as its `variants.tsv` counterpart |
| `variants.vcf` | The same sites as a VCF, plus a second ALT allele at 13625741, a `LowQual` record (15268954, for `--pass-only`), a record with no ID, and a symbolic `<DEL>` that should be skipped |
| `reference.fa.gz` (+ `.fai`, `.gzi`) | Contig `1` of mm10 up to 15.27 Mb, all N except ±3 kb around each site (~70 KB). Passed as `-r` (required) in every test |
| `barcodes.tsv` | 246 of the 327 cell barcodes in the BAM. About every 4th barcode is left out so the barcode filter has reads to drop. |

## Sites

| Site | Type | What it exercises |
|------|------|-------------------|
| `1:12904213 C>-` (Slco5a1) | 1-bp deletion in a C homopolymer | Deletion branch (reads with `is_del`) |
| `1:13562395 ->A` (Tram1) | 1-bp insertion in an A homopolymer | Insertion branch. The aligner places the insertion inconsistently in the `AAA` run, so left-aligning makes 2 more reads alt. One more read is aligned to the end of the run and then soft-clipped with the inserted A in the clip; realigning makes it alt (ref/alt reads 14/13; 16/10 from the CIGAR at the exact position alone) |
| `1:13625741 T>C` (Lactb2) | SNV, ~50% VAF | SNP branch. Many spliced reads span the site (`is_refskip`) and many bases have BQ 0 |
| `1:15268954 T>C` (Gm10566) | SNV | Includes multimapping reads with MAPQ 1 or 3, so `-mq` changes the result |
| `1:10000001 T>G` | No coverage | Zero-depth handling (VAF 0.0) |
| `1:10066353 TATG>CATA` (Cspp1) | Multi-base substitution (two SNVs 3 bp apart) | MNP branch. About half the reads are spliced partway through the MNP and are skipped |

The BAM also has reads with no `CB`/`UB` tags and reads whose barcode is not in `barcodes.tsv`, so both filters are exercised.

## scATAC test data (`--atac`)

A tiny subset of the 10x Genomics PBMC 5k NextGEM scATAC-seq BAM (`atac_pbmc_5k_nextgem_possorted_bam.bam`, Cell Ranger ATAC, GRCh38).
It contains only reads overlapping the five sites in `atac_variants.tsv` (4,486 reads, ~190 KB).

| File | Contents |
|------|----------|
| `atac.bam`, `atac.bam.bai` | Reads overlapping the test sites, made with `samtools view -M` over the sites and sorted. Every `@SQ` line is kept, and `@PG` lines were removed. Duplicates are flagged, not removed, as Cell Ranger ATAC writes them. 67 reads have no `CB` |
| `atac_variants.tsv` | Five sites (see below), with indels in dash notation |
| `atac_reference.fa.gz` (+ `.fai`, `.gzi`) | `chr6` up to 31.28 Mb, `chr7` up to 5.54 Mb and `chr12` up to 6.54 Mb of hg38, all N except ±3 kb around all five sites (from the UCSC REST API; ~170 KB). Passed as `-r` |
| `atac_barcodes.tsv` | 1,053 of the 1,403 cell barcodes in the BAM, with every 4th one left out |

| Site | Type | What it exercises |
|------|------|-------------------|
| `chr12:6534761 T>G` (GAPDH) | Het SNV | Overlapping mates: 61 fragments have both mates over the site |
| `chr7:5530212 C>A` (ACTB) | Het SNV | A second SNV, on another chromosome |
| `chr6:31271428 C>T` (HLA-C) | Het SNV | 46 alt fragments have MAPQ < 30, so `-mq 30` changes the result |
| `chr12:6470901 C>-` | 1-bp deletion, not in a repeat | Deletion branch at ~1,000 non-duplicate reads. Two fragments have only a read with a mismatch next to the deletion, with a base or two of evidence, so it isn't called (449 / 182 when the CIGAR alone was trusted) |
| `chr7:5532355 ->GCTACA` | 6-bp insertion, shifts 1 bp right | Many reads are soft-clipped at the insertion, with the inserted bases in the clip. They're realigned and called alt, so every fragment is alt (10 / 10). Counting from the CIGAR alone called them ref (10 / 2), and dropped one fragment whose mates disagreed |

Results with `--atac` (`UB_DEPTH / UB_ALT`, where each fragment counts once). The SNV counts match an independent pysam count of distinct non-duplicate read names per allele with a listed `CB`:

| Site | default | `-mq 30` | `-bq 30` | duplicates kept (pysam) |
|------|---------|----------|----------|-------------------------|
| GAPDH SNV | 240 / 121 | 240 / 121 | 228 / 119 | 504 / 255 |
| ACTB SNV | 132 / 68 | 132 / 68 | 126 / 62 | 291 / 151 |
| HLA-C SNV | 320 / 176 | 274 / 130 | 298 / 156 | 657 / 340 |
| deletion | 447 / 182 | 447 / 182 | 447 / 182 | 1024 / 395 |
| insertion | 10 / 10 | 10 / 10 | 10 / 10 | |

Each PCR duplicate has its own read name, so the last column shows that counting duplicates would roughly double every count.
Without `--atac`, cb_sniffer stops with an error because no read has a `UB` tag.

## Test suite

`test_cb_sniffer.py` checks these results with pytest (needs `pysam` and `pytest`):

```
pytest test/
# or, in the Docker image:
docker run --rm -v $PWD:/src -w /src cb_sniffer sh -c 'pip install -q pytest && pytest test/'
```

Tests marked `xfail` document known bugs from `issues.md`. They're strict, so once a bug is fixed they report XPASS, and you should remove the marker.

## Running by hand

```
python3 cb_sniffer.py -r test/reference.fa.gz test/test.bam test/variants.tsv test/barcodes.tsv default
python3 cb_sniffer.py -r test/reference.fa.gz -mq 10 test/test.bam test/variants.tsv test/barcodes.tsv mq10
python3 cb_sniffer.py -r test/reference.fa.gz -bq 30 test/test.bam test/variants.tsv test/barcodes.tsv bq30
python3 cb_sniffer.py -r test/reference.fa.gz test/test.bam test/variants_vcf.tsv test/barcodes.tsv vcf
python3 cb_sniffer.py -r test/reference.fa.gz --pass-only test/test.bam test/variants.vcf test/barcodes.tsv vcf_pass
python3 cb_sniffer.py -r test/atac_reference.fa.gz --atac test/atac.bam test/atac_variants.tsv test/atac_barcodes.tsv atac
```

Results with VCF-style/MNP support added (`UB_DEPTH / UB_ALT` from `*_AllCounts.tsv`):

| Site | default | `-mq 10` | `-bq 30` |
|------|---------|----------|----------|
| 12904213 del | 3 / 1 | 3 / 1 | 3 / 1 |
| 13562395 ins | 2 / 1 | 2 / 1 | 2 / 1 |
| 13625741 SNV | 28 / 14 | 28 / 14 | 28 / 14 |
| 15268954 SNV | 20 / 13 | 16 / 10 | 20 / 13 |
| 10000001 | 0 / 0 | 0 / 0 | 0 / 0 |
| 10066353 MNP | 95 / 73 | 95 / 73 | 91 / 72 |

At the SNVs, `-bq` doesn't change these UMI-level totals, because the low-quality reads are PCR duplicates of UMIs that also have good-quality reads. At the MNP, a read is dropped if any of its 4 bases is below the threshold, so a few UMIs are lost.
Its effect shows in the summed per-UMI read counts (`ref_count` + `alt_count` in `*_counts_UB.tsv`). At 13625741, ref/alt reads are 57/52 with `-bq 0`, 52/46 by default, and 48/42 with `-bq 30`.
Fixing the remaining bugs in `issues.md` will change some of these numbers.
