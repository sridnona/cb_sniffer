cb_sniffer 
=============

Mutation barcode caller,  calls mutant and ref barcodes from 10x single cell data   
Bams aligned using 10x cellranger pipeline

Usage
----
```{shell}
python3 cb_sniffer.py -h
usage: cb_sniffer.py [-h] [-f FILTER] [-mq MAPQ] [-bq BASEQ] [-r REFERENCE]
                     [--pass-only]
                     bam_file variant_file barcodes upn

Parse CB barcodes from Single cell rna seq data

positional arguments:
  bam_file              BAM or CRAM file (CRAM requires -r)
  variant_file          variants: tsv file (header line optional), or VCF/BCF
                        (optionally gzipped)
  barcodes              list of good barcodes file (may be gzipped)
  upn                   upn/sample name: will be used as prefix for out_file

options:
  -h, --help            show this help message and exit
  -f FILTER, --filter FILTER
                        number of reads required per barcode default: 0
  -mq MAPQ, --mapq MAPQ
                        Skip read with mapq smaller than default : 0
  -bq BASEQ, --baseq BASEQ
                        Skip bases with base quality less than default : 1
  -r REFERENCE, --reference REFERENCE
                        reference FASTA (indexed, may be bgzipped); used to
                        left-align indels so ones in repeats are found
                        wherever the aligner put them. Required for CRAM
                        input, which is decoded with it
  --pass-only           VCF input: skip records whose FILTER is not PASS or .

```
#### BAM or CRAM
An indexed BAM or CRAM from Cell Ranger, with `CB` and `UB` tags. A CRAM needs `-r`, and the run stops with an error without it. The FASTA must be the reference the CRAM was compressed against, and it's also used for indel left-alignment.

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

#### Reference (`-r`)
An indexed FASTA (`samtools faidx`; bgzipped is fine). Required for CRAM input. It's used to left-align indels, both the variant and each read's indels, so an indel in a repeat matches wherever the aligner placed it in the repeat.
For a ref read, the read must span the whole repeat with no indel or splice in it; reads ending inside the repeat are neither ref nor alt.
Without `-r`, indels are matched only at the exact position given, and a warning is printed.
SNPs and multi-base substitutions don't use the reference.
With `-r`, the deleted bases of each deletion are also checked against the reference, and a warning is printed if they don't match.

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
    3. read counts for Chromium molecular barcode(UB_tag): _counts_UB.tsv  


* Additional info on tags
[bam_tags](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam)



Dependencies
-------

python3  
pysam  

Docker image
----------
* sridnona/python3:180925.v1

```{shell}
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -Is -q research-hpc -a "docker(sridnona/python3:180925.v1)" /bin/bash
```
