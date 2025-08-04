#!/bin/bash
#SBATCH --job-name=coverage_bams
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=9G
#SBATCH --time=1:00:00
#SBATCH --partition=agap_short

# -------------------------
# Initialization
# -------------------------
echo "=============================================="
echo " Starting Coverage Calculation Pipeline"
echo " Job ID       : $SLURM_JOB_ID"
echo " Start Time   : $(date)"
echo "=============================================="

# -------------------------
# Argument Parsing
# -------------------------
if [ $# -ne 3 ]; then
    echo "ERROR: Missing required arguments!"
    echo "Usage: $0 <BAM_DIR> <OUTPUT_FILE> <REF_FASTA>"
    echo "Example: $0 ./bams output.tsv genome.fasta"
    exit 1
fi

BAM_DIR="$1"
OUTPUT="$2"
REF_FASTA="$3"

echo "Parameters:"
echo " - BAM Directory   : $BAM_DIR"
echo " - Output File     : $OUTPUT"
echo " - Reference FASTA : $REF_FASTA"
echo "=============================================="

# -------------------------
# Module loading
# -------------------------
echo "Loading modules..."
module load bioinfo-cirad
module load samtools/1.14-bin
module load seqkit/2.8.1
echo "Modules loaded."
echo "=============================================="

# -------------------------
# Reference Stats
# -------------------------
echo "Calculating total reference length..."
REF_LENGTH=$(seqkit stats "$REF_FASTA" | awk 'NR==2 {gsub(",", "", $5); print $5}')

if [[ -z "$REF_LENGTH" || "$REF_LENGTH" -eq 0 ]]; then
    echo "ERROR: Invalid reference length from $REF_FASTA"
    exit 1
fi

echo "Total reference length: $REF_LENGTH bp"
echo "=============================================="

# -------------------------
# Prepare Output File
# -------------------------
echo -e "Sample\tReads_mapped\tMean_Read_Length\tCoverage" > "$OUTPUT"
echo "Initialized output file: $OUTPUT"
echo "=============================================="

# -------------------------
# BAM Processing
# -------------------------
BAM_FILES=("$BAM_DIR"/*.bam)
BAM_COUNT=${#BAM_FILES[@]}
COUNTER=1

for BAM in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)
    echo "[$COUNTER/$BAM_COUNT] Processing: $SAMPLE"

    # Detect paired-end data
    PAIRED_READS=$(samtools view -c -f 1 "$BAM" 2>/dev/null)

    if [[ "$PAIRED_READS" -gt 0 ]]; then
        echo " - Type: Paired-end"
        READS=$(samtools view -c -F 260 "$BAM")
        READS=$((READS * 2))  # Account for both mates
    else
        echo " - Type: Single-end"
        READS=$(samtools view -c -F 260 "$BAM")
    fi

    if [[ -z "$READS" || "$READS" -eq 0 ]]; then
        echo " - WARNING: No mapped reads"
        echo -e "$SAMPLE\t0\t0\t0.0000" >> "$OUTPUT"
        echo "----------------------------------------------"
        COUNTER=$((COUNTER + 1))
        continue
    fi

    # Get mean read length
    MEAN_LENGTH=$(samtools stats "$BAM" | awk '/^SN/ && /average length:/ {print int($4)}')

    # Calculate average coverage
    COVERAGE=$(awk -v r="$READS" -v l="$MEAN_LENGTH" -v g="$REF_LENGTH" \
        'BEGIN { printf "%.4f", (r * l) / g }')

    echo " - Mapped Reads     : $READS"
    echo " - Mean Read Length : $MEAN_LENGTH"
    echo " - Coverage         : $COVERAGE"

    # Write to output
    echo -e "$SAMPLE\t$READS\t$MEAN_LENGTH\t$COVERAGE" >> "$OUTPUT"
    echo "----------------------------------------------"

    COUNTER=$((COUNTER + 1))
done

# -------------------------
# Completion
# -------------------------
echo "=============================================="
echo "Coverage calculation completed!"
echo "Results saved to: $OUTPUT"
echo "End Time: $(date)"
echo "=============================================="

