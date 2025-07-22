#!/bin/bash
#SBATCH --job-name=coverage_bams
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=9G
#SBATCH --time=1:00:00
#SBATCH --partition=agap_short

echo "=============================================="
echo "Starting Coverage Calculation Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "=============================================="

# Check if required parameters are provided
if [ $# -ne 4 ]; then
    echo "ERROR: Missing required parameters!"
    echo "Usage: $0 <BAM_DIR> <OUTPUT_PATH> <REF_FASTA> <READ_LENGTH>"
    echo "Example: $0 /path/to/bams /path/to/output.txt /path/to/reference.fasta 125"
    exit 1
fi

# Assign parameters
BAM_DIR="$1"
OUTPUT="$2"
REF_FASTA="$3"
READ_LENGTH="$4"

echo "=============================================="
echo "Parameters:"
echo "BAM Directory: $BAM_DIR"
echo "Output File: $OUTPUT"
echo "Reference FASTA: $REF_FASTA"
echo "Read Length: $READ_LENGTH"
echo "=============================================="

# Load required modules
echo "Loading required modules..."
module load bioinfo-cirad
module load seqkit/2.8.1
module load samtools/1.14-bin
echo "Modules loaded successfully."

# Get the total reference length
echo "Calculating reference length..."
REF_LENGTH=$(seqkit stats $REF_FASTA | awk 'NR==2 {print $5}' | tr -d ',')
echo "Reference length calculation complete."

# Check if reference length was retrieved successfully
if [[ -z "$REF_LENGTH" ]]; then
    echo "ERROR: Could not retrieve reference length from $REF_FASTA"
    exit 1
fi

echo "=============================================="
echo "Reference Statistics:"
echo "Total Reference Length: $REF_LENGTH"
echo "=============================================="

# Initialize output file
echo "Initializing output file..."
echo -e "Sample\tReads_mapped\tAverage_Coverage" > $OUTPUT
echo "Output file ready at $OUTPUT"

# Count number of BAM files to process
BAM_COUNT=$(ls $BAM_DIR/*.bam | wc -l)
echo "Found $BAM_COUNT BAM files to process"
echo "=============================================="

# Process each BAM file
COUNTER=1
for BAM in $BAM_DIR/*.bam; do
    SAMPLE=$(basename $BAM .bam)
    echo "Processing sample $COUNTER/$BAM_COUNT: $SAMPLE"

    # Count mapped reads
    echo " - Counting mapped reads..."
    READS=$(samtools view -c -F 4 $BAM)

    if [[ -z "$READS" ]]; then
        READS=0
        echo " - WARNING: No mapped reads found for $SAMPLE"
    fi

    # Calculate coverage
    echo " - Calculating coverage..."
    COV=$(awk -v reads=$READS -v rlen=$READ_LENGTH -v len=$REF_LENGTH 'BEGIN {printf "%.4f", (reads * rlen) / len}')

    # Print sample statistics
    echo " - Sample Statistics:"
    echo "   - Reads Mapped: $READS"
    echo "   - Read Length: $READ_LENGTH"
    echo "   - Coverage: $COV"

    # Write to output
    echo -e "$SAMPLE\t$READS\t$COV" >> $OUTPUT
    echo " - Results written to $OUTPUT"
    echo "----------------------------------------------"

    COUNTER=$((COUNTER+1))
done

echo "=============================================="
echo "Pipeline completed successfully!"
echo "Final results saved to: $OUTPUT"
echo "End Time: $(date)"
echo "=============================================="
