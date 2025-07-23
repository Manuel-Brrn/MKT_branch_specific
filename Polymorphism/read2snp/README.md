# read2snp SNP Calling Pipeline

## Description
This pipeline uses reads2snp to identify SNPs and create consensus sequences from BAM alignments against a reference genome.

## Usage
```bash
sbatch read2snp_speltoides.sbatch <bam_list.txt> <reference.fasta> <output_directory>
```
#Example
#bash

```bash
sbatch read2snp_speltoides.sbatch \
    /path/to/bam_list.txt \
    /path/to/reference.fasta \
    /path/to/output_dir
```

#Required Arguments

    bam_list.txt: File containing paths to BAM files (one per line)

    reference.fasta: Reference genome in FASTA format

    output_directory: Where to store results

Default Parameters
Parameter       Value         Description
-th2            0.05          Minor allele frequency threshold
-fis            0.3           Inbreeding coefficient threshold
-nbth           4             Number of threads
-pre            0.000001      P-value threshold
-rlg            50            Read length
-bqt            20            Base quality threshold
-rqt            30            Read quality threshold
