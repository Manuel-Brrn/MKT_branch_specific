



**impMKT aligment:**

**Add the reference sequence of the focal species:**

**Extracts sequences from a FASTA file based on contig IDs using parallel processing.**
otbain : contigs_ID.txt ; list of genes to be extracted
```bash
grep ">" urartu_covered_cds.cds | cut -d "|" -f1 | sed 's/^>//' | uniq > contigs_ID.txt
```

## Usage
```bash
sbatch extract_sequences.sbatch <output_directory> <fasta_file_path>
```

extract_sequences_per_contig.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=sequences_extraction             # Job name for Slurm
#SBATCH --output=./log_%j_%x_out.txt                # Log file for standard output
#SBATCH --error=./log_%j_%x_err.txt                 # Log file for standard error
#SBATCH --nodes=1                                   # Number of nodes
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G                                   # Memory allocation
#SBATCH --time=10:00:00                             # Time limit
#SBATCH --partition=agap_normal                     # Job queue/partition

# Check if required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <output_directory> <fasta_file_path>"
    echo "Example: $0 /home/barrientosm/scratch/test/speltoides /path/to/modified_2.fasta"
    exit 1
fi

echo "Loading required modules..."
module load bioinfo-cirad
module load seqkit/2.8.1
module load parallel/20220922
echo "Modules loaded successfully."

# Define output directory and input FASTA file from command line arguments
export OUTPUT_DIR="$1"
export FASTA_FILE="$2"

echo "Output directory set to: $OUTPUT_DIR"
echo "FASTA file set to: $FASTA_FILE"

# Create the output directory if it doesn't exist
echo "Creating output directory if it doesn't exist..."
mkdir -p "$OUTPUT_DIR"
echo "Output directory ready."

# Extract unique contig IDs from FASTA headers and remove leading ">"
echo "Extracting unique contig IDs from FASTA headers..."
grep "^>" "$FASTA_FILE" | cut -d"|" -f1 | sed 's/^>//' | sort | uniq > contigs_ID.txt
echo "Contig IDs written to contigs_ID.txt"

# Define the function to process a single contig
process_contig() {
    contig="$1"
    echo "Processing contig: $contig"
    seqkit grep -r -p "$contig" "$FASTA_FILE" -o "${OUTPUT_DIR}/${contig}.fasta"
    if [ $? -eq 0 ]; then
        echo "Successfully processed contig: $contig"
    else
        echo "Error processing contig: $contig" >&2
    fi
}
export -f process_contig

echo "Starting parallel processing of contigs..."
# Run in parallel, processing max 20 contigs at the same time
unset PERL5LIB
cat contigs_ID.txt | parallel -j 20 --joblog "${OUTPUT_DIR}/parallel_joblog.txt" process_contig {}

# Check if parallel completed successfully
if [ $? -eq 0 ]; then
    echo "All contigs processed successfully."
else
    echo "Error: Some contigs failed to process. Check the logs for details."
    exit 1
fi

echo "Script completed."
```
**Add the reference CDS sequences to each file**
This SLURM script matches contig sequences from input FASTA files with their corresponding CDS sequences and combines them into output files.
It processes files in parallel for efficiency.

Input Data:
- **Input Directory**: Directory containing FASTA files (expected format: `*_1.fasta`)
- **CDS File**: A FASTA file containing all CDS sequences to match against
- **Output Directory**: Where processed files will be saved

Output:
- Processed FASTA files (original contig + matching CDS sequence) in output directory
- Job log file (`parallel_joblog.txt`) in output directory
- Console output with processing statistics

Processing Steps
1. For each input FASTA file:
   - Extracts matching CDS sequence from CDS file
   - Combines CDS sequence with original contig
   - Saves combined sequence to output directory
2. Runs in parallel using all allocated CPUs
3. Provides summary statistics upon completion

Usage
```bash
sbatch add_reference_sequence_all_genes.sbatch <input_fastas_dir> <cds_file> <output_dir>
```

add_reference_sequence_all_genes.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=add_reference
#SBATCH --output=add_reference_parallel_%j_out.txt
#SBATCH --error=add_reference_parallel_%j_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=20:00:00
#SBATCH --partition=agap_normal

module load bioinfo-cirad
module load seqkit/2.8.1
module load parallel/20220922

# Check if required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "ERROR: Missing arguments!"
    echo "Usage: $0 <input_fastas_dir> <cds_file> <output_dir>"
    echo "Example: $0 /path/to/fastas /path/to/species_CDS.fasta /path/to/output"
    exit 1
fi

# Input parameters from command line arguments
INPUT_DIR="$1"
CDS_FILE="$2"
OUTPUT_DIR="$3"

echo "Starting parallel contig matching"
echo "----------------------------------------"
echo "Input directory: $INPUT_DIR"
echo "CDS file: $CDS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo ""

# Verify input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Verify CDS file exists
if [ ! -f "$CDS_FILE" ]; then
    echo "ERROR: CDS file does not exist: $CDS_FILE"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create output directory: $OUTPUT_DIR"
    exit 1
fi

# Count number of fasta files to process
shopt -s nullglob
FASTA_FILES=("$INPUT_DIR"/*.fasta)
shopt -u nullglob  # Reset

if [ ${#FASTA_FILES[@]} -eq 0 ]; then
    echo "ERROR: No files matching '*.fasta' found in $INPUT_DIR"
    exit 1
fi

FASTA_COUNT=${#FASTA_FILES[@]}
echo "Found $FASTA_COUNT fasta files to process"
echo ""

# Create temporary directory for parallel processing
TMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TMP_DIR"
echo ""

# Function to process a single contig file
process_contig() {
    local CONTIG_FILE="$1"
    local CDS_FILE="$2"
    local OUTPUT_DIR="$3"
    local TMP_DIR="$4"

    # Extract just the filename without path
    local CONTIG_FILENAME=$(basename "$CONTIG_FILE")

    # Extract base filename without extension (now handles T01.fasta)
    local BASENAME=$(basename "$CONTIG_FILE" .fasta)
    local ID_PATTERN=$(echo "$BASENAME" | cut -d'|' -f1)  # Get the part before first |

    # Temporary files
    local TMP_CDS="${TMP_DIR}/${CONTIG_FILENAME}_cds.tmp"
    local TMP_OUT="${TMP_DIR}/${CONTIG_FILENAME}_out.tmp"

    # Extract matching CDS sequence using more flexible pattern matching
     seqkit grep -r -i -n -p "${ID_PATTERN}" "$CDS_FILE" > "$TMP_CDS" 2>/dev/null

    # Check if we found the CDS
    if [ ! -s "$TMP_CDS" ]; then
        echo "WARNING: No CDS found for $CONTIG_FILENAME (tried pattern: >${ID_PATTERN})"
        return 1
    fi

    # Combine CDS first, then original contig
    cat "$TMP_CDS" "$CONTIG_FILE" > "$TMP_OUT"
    mv "$TMP_OUT" "${OUTPUT_DIR}/${CONTIG_FILENAME}"
    rm -f "$TMP_CDS"

    echo "Processed: $CONTIG_FILENAME"
    return 0
}

export -f process_contig

# Process files in parallel
echo "Starting parallel processing with $SLURM_CPUS_PER_TASK jobs..."
echo ""

export PARALLEL_NO_TTY=1  # Prevents terminal-related errors in SLURM
START_TIME=$(date +%s)
unset PERL5LIB
printf "%s\n" "${FASTA_FILES[@]}" | parallel --will-cite -j $SLURM_CPUS_PER_TASK \
    --joblog "${OUTPUT_DIR}/parallel_joblog.txt" \
    --progress --eta \
    "process_contig {} '$CDS_FILE' '$OUTPUT_DIR' '$TMP_DIR'" 2>&1 | tee "${OUTPUT_DIR}/parallel.log"

PARALLEL_EXIT=$?
END_TIME=$(date +%s)
RUNTIME=$((END_TIME-START_TIME))

echo ""
echo "Parallel processing completed in $RUNTIME seconds"
echo ""

echo "First 5 files to process:"
printf "%s\n" "${FASTA_FILES[@]:0:5}"

# Count results
PROCESSED_COUNT=$(find "$OUTPUT_DIR" -type f -name "*_1.fasta" | wc -l)
SKIPPED_COUNT=$((FASTA_COUNT - PROCESSED_COUNT))
printf "%s\n" "${FASTA_FILES[@]:0:5}"
echo "----------------------------------------"
echo "Batch processing completed"
echo "Summary:"
echo "Total files found:    $FASTA_COUNT"
echo "Successfully processed: $PROCESSED_COUNT"
echo "Skipped files:        $SKIPPED_COUNT"
echo "Parallel exit status: $PARALLEL_EXIT"
echo ""

# Clean up temporary directory
rm -rf "$TMP_DIR"
echo "Cleaned up temporary directory: $TMP_DIR"
echo ""

# Verify at least some files were processed
if [ "$PROCESSED_COUNT" -eq 0 ]; then
    echo "ERROR: No files were successfully processed!"
    exit 1
fi

exit 0
```

**Change the reference sequence name to run MACSE:**
```bash
for file in *.fasta; do
    sed -i 's/^>T_urartu|\(.*\)/>\1|T_urartu/' "$file"
done
```
Adapt ">T_urartu" to the species name




**Clean alignments with Hmmcleaner**
```bash
module load bioinfo-cirad
module load hmmer/3.3.2-singularity
# Define the input directory
INPUT_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/MACSE/urartu_covered_impMKT"

# Process each aligned NT fasta file
for f in "$INPUT_DIR"/*_aligned_NT.fasta; do
    base=$(basename "$f" .fasta)  # includes _aligned_NT

    out_fasta="${INPUT_DIR}/${base}_hmm.fasta"
    out_log="${INPUT_DIR}/${base}_hmm.log"
    out_score="${INPUT_DIR}/${base}_hmm.score"

    if [[ -f "$out_fasta" && -f "$out_log" && -f "$out_score" ]]; then
        echo "=== Skipping $f (already cleaned) ==="
    else
        echo "=== Processing $f ==="
        env PERL5LIB=/home/barrientosm/my_perl_5.36.0/lib/perl5:/home/barrientosm/my_perl_5.36.0/lib/perl5/x86_64-linux \
            /home/barrientosm/my_perl_5.36.0/bin/perl -Ilib bin/HmmCleaner.pl "$f"
    fi
done > hmmcleaner_urartu_impMKT.log 2>&1
```
