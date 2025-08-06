#!/bin/bash

module load bioinfo-cirad
module load samtools/1.14-bin

# Define output directory
################### ADAPT THE PATH ##########
OUTDIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/data_RNAseq/T_boeticum/1_sp_analysis/fastq_from_bams"

# Make sure output directory exists
mkdir -p "$OUTDIR"

# Loop over BAMs
while read BAM; do
    # Extract the sample name, 
    SAMPLE=$(basename "$BAM" .bam)

    echo "Processing $SAMPLE..."

    # Sort by query name (required for proper fastq conversion)
    samtools sort -n "$BAM" -o "${OUTDIR}/${SAMPLE}_namesorted.bam"

    # Convert to FASTQ (paired-end)
    samtools fastq -@ 4 \
        -1 "${OUTDIR}/${SAMPLE}_R1.fastq" \
        -2 "${OUTDIR}/${SAMPLE}_R2.fastq" \
        -0 /dev/null -s /dev/null -n \
        "${OUTDIR}/${SAMPLE}_namesorted.bam"

    # Optionally remove intermediate sorted BAM to save space
    rm "${OUTDIR}/${SAMPLE}_namesorted.bam"

done < bam_list.txt
