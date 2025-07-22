#!/bin/bash
#SBATCH --job-name=depth_analysis
#SBATCH --output=depth_analysis_%j.log
#SBATCH --error=depth_analysis_%j.err
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=agap_short

# Check input arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/bam_files /path/to/output_directory"
    exit 1
fi

# Initialize variables
BAM_DIR="$1"
OUTPUT_DIR="$2"
FINAL_OUTPUT="$OUTPUT_DIR/depth_per_position_all_individuals_BAMs.depth"
TEMP_DIR="$OUTPUT_DIR/temp"
LOG_FILE="$OUTPUT_DIR/depth_analysis.log"

# Create directories
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR" || {
    echo "ERROR: Failed to create output directories" | tee -a "$LOG_FILE"
    exit 1
}

# Initialize log file
echo "=== Depth Analysis Pipeline ===" > "$LOG_FILE"
echo "Start time: $(date)" >> "$LOG_FILE"
echo "Input directory: $BAM_DIR" >> "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" >> "$LOG_FILE"

# Load required modules
echo "Loading modules..." | tee -a "$LOG_FILE"
module purge >> "$LOG_FILE" 2>&1
module load bioinfo-cirad samtools/1.14-bin >> "$LOG_FILE" 2>&1 || {
    echo "ERROR: Module loading failed" | tee -a "$LOG_FILE"
    exit 1
}

# Get list of BAM files
BAM_FILES=($(find "$BAM_DIR" -name "*.bam"))
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "ERROR: No BAM files found in $BAM_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM files to process" | tee -a "$LOG_FILE"

# Step 1: Generate depth files
echo -e "\n=== Step 1: Generating depth files ===" | tee -a "$LOG_FILE"

for bam_file in "${BAM_FILES[@]}"; do
    echo "Processing $bam_file" | tee -a "$LOG_FILE"
    sample=$(basename "$bam_file" .bam)
    samtools depth -aa "$bam_file" > "$TEMP_DIR/${sample}.depth" 2>> "$LOG_FILE"
done

# Verify depth files were created
DEPTH_FILES=($(find "$TEMP_DIR" -name "*.depth"))
if [ ${#DEPTH_FILES[@]} -ne ${#BAM_FILES[@]} ]; then
    echo "WARNING: Only ${#DEPTH_FILES[@]}/${#BAM_FILES[@]} depth files generated" | tee -a "$LOG_FILE"
fi

# Step 2: Merge depth files
echo -e "\n=== Step 2: Merging depth files ===" | tee -a "$LOG_FILE"

# Use the first file as reference
REF_FILE="${DEPTH_FILES[0]}"
if [ ! -s "$REF_FILE" ]; then
    echo "ERROR: Reference depth file is empty" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Using $REF_FILE for reference positions" | tee -a "$LOG_FILE"

# Build header
HEADER="contig\tpos"
for file in "${DEPTH_FILES[@]}"; do
    sample=$(basename "$file" .depth)
    HEADER+="\t$sample"
done
echo -e "$HEADER" > "$FINAL_OUTPUT"

# Paste and format all files
paste "${DEPTH_FILES[@]}" | awk '
BEGIN { OFS="\t" }
{
    contig = $1;
    pos = $2;
    printf "%s\t%s", contig, pos;
    for (i = 3; i <= NF; i += 3) {
        printf "\t%s", $(i);
    }
    printf "\n";
}' >> "$FINAL_OUTPUT"

# Final check
if [ ! -s "$FINAL_OUTPUT" ]; then
    echo "ERROR: Failed to create final output file" | tee -a "$LOG_FILE"
    exit 1
fi

# Cleanup
echo -e "\n=== Cleaning up ===" | tee -a "$LOG_FILE"
rm -rf "$TEMP_DIR"

# Done
echo -e "\n=== Analysis Complete ===" | tee -a "$LOG_FILE"
echo "Output file: $FINAL_OUTPUT" | tee -a "$LOG_FILE"
echo "Total positions: $(($(wc -l < "$FINAL_OUTPUT")-1))" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
