10x_singlecell 
=============

Mutation barcode caller,  calls mutant and ref barcodes from 10x single cell data   
Bams aligned using 10x cellranger pipeline

Usage
----
```{shell}
python3 cb_sniffer.py -h
usage: cb_sniffer.py [-h] [-f FILTER] [-mq MAPQ] [-bq BASEQ]
                     bam_file variant_file barcodes upn

Parse CB barcodes from Single cell rna seq data

positional arguments:
  bam_file              BAM file
  variant_file          variants file with header
  barcodes              list of good barcodes file
  upn                   upn/sample name: will be used as prefix for out_file

optional arguments:
  -h, --help            show this help message and exit
  -f FILTER, --filter FILTER
                        number of reads required per barcode default: 0
  -mq MAPQ, --mapq MAPQ
                        Skip read with mapq smaller than default : 0
  -bq BASEQ, --baseq BASEQ
                        Skip bases with base quality less than default : 1

```
#### Variant file 
```{shell}
chrm	start	stop	ref	var	gene_name	trv_type
19	4xxx	4xxx	G	A	gene_name	silent/misense/frame_shift  
```

#### barcodes
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


Dependencies
-------

python3  
pysam  

Docker image
----------
* sridnona/python3:180925.v1

```{shell}
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -Is -q research-hpc -a "docker(sridnona/python3:180924.v1)" /bin/bash
```
