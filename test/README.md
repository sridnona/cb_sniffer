# Test data

A tiny subset of a mouse (mm10) 10x Genomics scRNA-seq BAM (Cell Ranger 9.0.0, C57BL/6 bone marrow).
It contains only reads overlapping the six sites in `variants.tsv` (1,476 reads, ~100 KB).

| File | Contents |
|------|----------|
| `test.bam`, `test.bam.bai` | Reads overlapping the test sites on chr1. `@PG` header lines were removed. |
| `variants.tsv` | Six sites (see below), indels in dash notation |
| `variants_vcf.tsv` | Four of the same sites written VCF-style, with padding bases. Each must give the same counts as its `variants.tsv` counterpart |
| `variants.vcf` | The same sites as a VCF, plus a second ALT allele at 13625741, a `LowQual` record (15268954, for `--pass-only`), a record with no ID, and a symbolic `<DEL>` that should be skipped |
| `barcodes.tsv` | 246 of the 327 cell barcodes in the BAM. About every 4th barcode is left out so the barcode filter has reads to drop. |

## Sites

| Site | Type | What it exercises |
|------|------|-------------------|
| `1:12904213 C>-` (Slco5a1) | 1-bp deletion in a C homopolymer | Deletion branch (reads with `is_del`) |
| `1:13562395 ->A` (Tram1) | 1-bp insertion in an A homopolymer | Insertion branch. The aligner places the insertion inconsistently in the repeat (issues.md #10) |
| `1:13625741 T>C` (Lactb2) | SNV, ~50% VAF | SNP branch. Many spliced reads span the site (`is_refskip`) and many bases have BQ 0 |
| `1:15268954 T>C` (Gm10566) | SNV | Includes multimapping reads with MAPQ 1 or 3, so `-mq` changes the result |
| `1:10000001 T>G` | No coverage | Zero-depth handling (VAF 0.0) |
| `1:10066353 TATG>CATA` (Cspp1) | Multi-base substitution (two SNVs 3 bp apart) | MNP branch. About half the reads are spliced partway through the MNP and are skipped |

The BAM also has reads with no `CB`/`UB` tags and reads whose barcode is not in `barcodes.tsv`, so both filters are exercised.

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
python3 cb_sniffer.py test/test.bam test/variants.tsv test/barcodes.tsv default
python3 cb_sniffer.py -mq 10 test/test.bam test/variants.tsv test/barcodes.tsv mq10
python3 cb_sniffer.py -bq 30 test/test.bam test/variants.tsv test/barcodes.tsv bq30
python3 cb_sniffer.py test/test.bam test/variants_vcf.tsv test/barcodes.tsv vcf
python3 cb_sniffer.py --pass-only test/test.bam test/variants.vcf test/barcodes.tsv vcf_pass
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
