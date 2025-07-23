#!/bin/bash
#SBATCH --job-name=read2snp
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=70:00:00
#SBATCH --partition=agap_long

echo "=== Starting read2snp analysis ==="
echo "Job started at: $(date)"
echo "Running on host: $(hostname)"

# Load required modules
echo "Loading modules..."
module load bioinfo-cirad
module load reads2snp/2.0-bin

# Usage instructions
if [ $# -lt 3 ]; then
    echo "Usage: $0 <bam_list_path> <reference_fasta> <output_directory>"
    echo "Example:"
    echo "  sbatch $0 /path/to/bam_list.txt /path/to/reference.fasta /path/to/output"
    exit 1
fi

# Assign input parameters
BAM_LIST="$1"
REFERENCE="$2"
OUTPUT_DIR="$3"

# Analysis parameters (customize as needed)
TH2_VALUE=0.05
FIS_VALUE=0.9

echo "=== Analysis Parameters ==="
echo "BAM list: ${BAM_LIST}"
echo "Reference genome: ${REFERENCE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "TH2 threshold: ${TH2_VALUE}"
echo "FIS threshold: ${FIS_VALUE}"
echo ""

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Verify input files exist
echo "Verifying input files..."
if [ ! -f "${BAM_LIST}" ]; then
    echo "ERROR: BAM list file not found at ${BAM_LIST}" >&2
    exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference file not found at ${REFERENCE}" >&2
    exit 1
fi

# Run reads2snp
echo "Starting reads2snp analysis..."
reads2snp_2.0.64.bin \
    -bamlist "${BAM_LIST}" \
    -bamref "${REFERENCE}" \
    -th2 "${TH2_VALUE}" \
    -nbth 4 \
    -fis "${FIS_VALUE}" \
    -pre 0.000001 \
    -rlg 50 \
    -bqt 20 \
    -rqt 30 \
    -out "${OUTPUT_DIR}/urartu"

# Check exit status
if [ $? -eq 0 ]; then
    echo "SNP extraction completed successfully."
    echo "Results saved to: ${OUTPUT_DIR}"
else
    echo "ERROR: reads2snp failed with exit code $?" >&2
    exit 1
fi

echo "Job finished at: $(date)"
