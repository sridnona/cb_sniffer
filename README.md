cb_sniffer 
=============

Mutation barcode caller,  calls mutant and ref barcodes from 10x single cell data   
Bams aligned using 10x cellranger pipeline

Usage
----
```{shell}
python3 cb_sniffer.py -h
usage: cb_sniffer.py [-h] [-mq MAPQ] [-bq BASEQ] -r REFERENCE [--pass-only]
                     [--atac]
                     bam_file variant_file barcodes upn

Parse CB barcodes from Single cell rna seq data

positional arguments:
  bam_file              BAM or CRAM file
  variant_file          variants: tsv file (header line optional), or VCF/BCF
                        (optionally gzipped)
  barcodes              list of good barcodes file (may be gzipped)
  upn                   upn/sample name: will be used as prefix for out_file

options:
  -h, --help            show this help message and exit
  -mq, --mapq MAPQ      Skip read with mapq smaller than default : 0
  -bq, --baseq BASEQ    Skip bases with base quality less than default : 1
  -r, --reference REFERENCE
                        reference FASTA the reads were aligned to (indexed,
                        may be bgzipped). Used to check variant ref alleles,
                        to find indels wherever the aligner put them in a
                        repeat, to realign soft-clipped reads at indels, and
                        to decode CRAM input
  --pass-only           VCF input: skip records whose FILTER is not PASS or .
  --atac                scATAC BAM (no UB/UMI tags): count each fragment (read
                        pair) once in place of a UMI, and skip duplicate,
                        secondary and QC-fail reads

```
#### BAM or CRAM
An indexed BAM or CRAM from Cell Ranger, with `CB` and `UB` tags. For a CRAM, the `-r` FASTA must be the reference it was compressed against.

Each read is counted at most once per variant. Where a read is split into several alignments, the primary alignment's call is used, and a supplementary (0x800) or secondary (0x100) alignment counts only when the primary doesn't cover the site, as with split reads across a large deletion or a fusion junction. (CellRanger BAMs usually have none, but long-read single-cell BAMs from minimap2 often do)

The default is to assume the input bam is from scRNA and has CB and UB tags. If none of the reads with a listed barcode have a `UB` tag, the run stops and suggests `--atac`.

#### scATAC and Multiome (`--atac`)
scATAC BAMs (Cell Ranger ATAC or ARC) have `CB` tags but no `UB` UMIs. With `--atac`:
- Each fragment (read name) takes the place of a UMI. Both mates of a pair are counted, and the 75% consensus rule is applied per fragment, so overlapping mates count once and a fragment whose mates disagree is dropped.
- Reads flagged duplicate (0x400), secondary (0x100) or QC-fail (0x200) are skipped. Supplementary reads are handled as described above.
- The output columns are unchanged: in `_counts_UB.tsv` and the `CB:UB_barcodes` column, the barcode is `CB:read_name`.
- For Multiome (ARC) runs, Cell Ranger translates ATAC barcodes to their paired gene expression barcodes, so the same barcode list works for both BAMs.
- ATAC reads are aligned with bwa (MAPQ 0-60), so `-mq 30` is a reasonable filter.
- ATAC coverage at any one site is sparse, and outside open chromatin there may be no reads at all; expect far fewer cells per variant than with scRNA.

If `--atac` is given but the reads have `UB` tags, a warning is printed and reads are still counted per fragment, not per UMI.

#### Variant file
Either a VCF/BCF file (plain, gzipped or bgzipped, no index needed) or a tab-separated file. The format is detected from the file contents.

##### VCF
- Each ALT allele of a multi-allelic record is counted separately, as its own output row.
- Symbolic (`<DEL>`), breakend and `*` alleles are skipped with a warning.
- `--pass-only` skips records whose FILTER isn't `PASS` or `.`.
- In the output, `gene` is the record's `ID` (`.` if missing), and `type` is the variant kind: `snp`, `mnp`, `ins` or `del`.
- Sample columns are ignored.

##### TSV
```{shell}
chrm	start	stop	ref	var	gene_name	trv_type
19	4xxx	4xxx	G	A	gene_name	silent/misense/frame_shift  
```

The header line is optional. The first line is treated as a header if it starts with `#` or its `start` column isn't an integer; otherwise it's read as a variant.

Accepted `ref`/`var` notations (`start` is 1-based; `stop` isn't used):

| Event | Example ref / var | `start` |
|-------|-------------------|---------|
| SNV | `G` / `A` | the base |
| Multi-base substitution | `ATC` / `GGC` | first base |
| Deletion, dash style | `CT` / `-` | first deleted base |
| Insertion, dash style | `-` / `AC` | base before the insertion |
| Deletion, VCF style | `ACT` / `A` | anchor base (VCF POS) |
| Insertion, VCF style | `A` / `AAC` | anchor base (VCF POS) |

Shared padding bases are trimmed, so `TG` / `CG` is counted as the SNV `T` / `C`. Substitutions whose ref and alt still differ in length after trimming (e.g. `ATC` / `GG`) aren't supported and are skipped with a warning.
For multi-base substitutions, only reads covering every base without a gap are counted, and every base must pass `-bq`.

#### Reference (`-r`, required)
An indexed FASTA (`samtools faidx`; bgzipped is fine): the reference the reads were aligned to. `-r` is now required. Earlier versions made it optional, and without it matched indels only at the exact position given.

- The ref allele of every variant is checked against it, and a warning is printed if it doesn't match. A chromosome missing from the FASTA stops the run.
- Indels are left-aligned, both the variant and each read's indels, so an indel in a repeat matches wherever the aligner placed it in the repeat.
- A read is called from its CIGAR only if it's aligned with no soft clip, splice or mismatch across the whole repeat the indel is in, plus 5 bases further on any side where the alignment stops at a clip or splice. Such a read is ref if it has no indel there, and alt if its only indel there is the variant.
- Other reads near an indel, such as reads soft-clipped at an insertion or tandem duplication, or reads with a deletion near their end aligned as mismatches, are realigned. The part of the read near the indel, soft clips included, is aligned (Smith-Waterman, with BWA-MEM's scores) to the reference and to the reference with the indel applied. The read is called if it scores at least 3 better against one than the other. Reads that end inside the repeat score the same against both and aren't counted.
- A read spliced across the repeat isn't counted.

SNPs and multi-base substitutions use the reference only for the ref allele check.

#### barcodes
One barcode per line, plain or gzipped (e.g. Cell Ranger's `barcodes.tsv.gz`), detected from the file contents.
```{shell}
AACCTGAGAATGTTG-1
AAACCTGAGCTACCGC-1
AAACCTGAGCTGCCCA-1
AAACCTGAGGTCGGAT-1
AAACCTGAGTACGATA-1
AAACCTGAGTGCAAGC-1
AAACCTGAGTGCCATT-1
AAACCTGAGTGGAGTC-1
AAACCTGCAAAGTGCG-1
AAACCTGCATAAGACA-1
```

### Output files 
    1. variant list with barcodes supporting a variant :_AllCounts.tsv
    2. read counts for Chromium cellular barcode sequence(CB_tag) : _counts_CB.tsv
    3. read counts for Chromium molecular barcode(UB_tag): _counts_UB.tsv (with --atac: per fragment)


* Additional info on tags
[bam_tags](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam)



Dependencies
-------

python3  
pysam (version pinned in `requirements.txt`)  

Docker image
----------
[mgibio/cb_sniffer](https://hub.docker.com/r/mgibio/cb_sniffer) (built from the `Dockerfile` in this repo; linux/amd64). The image puts a `cb_sniffer` command on the `PATH`.

```{shell}
docker run --rm -v $PWD:/data -w /data mgibio/cb_sniffer \
    cb_sniffer -r ref.fa sample.bam variants.vcf barcodes.tsv.gz sample
```

On LSF:
```{shell}
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -Is -q research-hpc -a "docker(mgibio/cb_sniffer)" /bin/bash
```

To build it yourself:
```{shell}
docker build --platform linux/amd64 -t mgibio/cb_sniffer .
```
