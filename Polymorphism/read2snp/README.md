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

## Cleaning
Use the list of problematic cds detected by VESPA and the the list contig covered at 50% (a covered position is a position covered at a depth of 10 for at least 7 individuals), to keep only valuable contig sequences
```bash
grep -v -F -f problematic_cds.txt high_coverage_contigs_id_urartu.txt > filtered_high_coverage_ids.txt
```

    grep

    -v (invert match)

        Normally grep finds matching lines
        With -v, it does the opposite: keeps only lines that don't match the patterns
        This removes problematic IDs instead of keeping them

    -F (fixed strings)

        Treats the search patterns as literal strings (not regular expressions)
        Important because your gene IDs contain dots (.) which are special characters in regex
        Without -F, you'd need to escape the dots (e.g., TuG1812G0100000450\.01\.T03)

    -f problematic_cds.txt

        Reads search patterns from a file (problematic_cds.txt) instead of typing them manually
        Each line in the file is treated as a separate pattern to match against

## Extraction of valuable contigs

```bash
#!/bin/bash
#SBATCH --job-name=seqkit
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --partition=agap_normal


module load bioinfo-cirad
module load seqkit/2.8.1

seqkit grep -r -f filtered_high_coverage_ids.txt urartu_renamed.fas -o urartu_high_coverage_clean.fasta
```
