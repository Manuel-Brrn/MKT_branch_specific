
**Alignment**

*dNdSpiNpiS input*
MACSE Alignment Pipeline

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


**Add the outgroup sequence to the multi fasta**
The script adds the outgroup sequence to each multi-fasta contig, and copy them to an input directory for aligning. Provides a log file with the contigs with no orthologs.

```bash
 ./add_outgroup.sh \
    /path/to/RBH_table.txt \
    /path/to/hvulgare_CDS.fasta \
    /path/to/output_dir/
```

```bash
#!/bin/bash

# ----------------------------
# Usage and argument parsing
# ----------------------------
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <RBH_table> <H_vulgare_fasta> <destination_directory>"
    echo "Example:"
    echo "  $0 RBH_table.txt hvulgare_CDS.fasta /path/to/output"
    exit 1
fi

RBH_TAB="$1"
H_VULGARE_FASTA="$2"
DEST_DIR="$3"

# ----------------------------
# Logging
# ----------------------------
MISSING_LOG="missing_orthologs.txt"
> "$MISSING_LOG"  # Clear log

module load bioinfo-cirad
module load seqkit/2.8.1
# ----------------------------
# Process FASTA files
# ----------------------------
for fasta in *.fasta; do
    gene_id=$(basename "$fasta" .fasta)
    echo "Processing $gene_id..."

    # Get EXACT H. vulgare ID from second column of RBH table (skip comments)
    h_ortholog=$(grep -v "^#" "$RBH_TAB" | awk -v id="$gene_id" '$1 == id { print $2 }')

    if [[ -z "$h_ortholog" ]]; then
        echo "No ortholog found in RBH for $gene_id" | tee -a "$MISSING_LOG"
        continue
    fi

    echo "  Hordeum ortholog ID: $h_ortholog"

    # Extract sequence using EXACT RBH-mapped ID
    if ! seqkit faidx "$H_VULGARE_FASTA" "$h_ortholog" > tmp_seq.fa 2>/dev/null; then
        echo "Ortholog '$h_ortholog' not found in FASTA for $gene_id" | tee -a "$MISSING_LOG"
        continue
    fi

    # Append ortholog to FASTA (only if tmp_seq.fa exists/non-empty)
    if [ -s tmp_seq.fa ]; then
        cat tmp_seq.fa >> "$fasta"
        echo "  Appended $h_ortholog to $fasta"
        cp "$fasta" "$DEST_DIR/"
    else
        echo "  Warning: tmp_seq.fa is empty for $h_ortholog" | tee -a "$MISSING_LOG"
    fi

    # Clean up
    rm -f tmp_seq.fa
done

echo "Done. Missing entries logged in: $MISSING_LOG"
```



Alignment with Macse: 

```bash
#!/bin/bash

# Load necessary modules
module load bioinfo-cirad
module load java/jre1.8.0_31

# Path to MACSE JAR
MACSE_JAR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/03_scripts/dn_ds_pipeline/MACSE/macse_v1.2.jar"

# Output directory
ALIGN_DIR="macse_alignments"
mkdir -p "$ALIGN_DIR"

# MACSE alignment parameters (you can change these if needed)
THREADS=4                         # Number of threads
GENETIC_CODE=1                   # NCBI genetic code table (1 = standard)
FRAMESHIFT_MODE="keep"          # How to handle frameshifts: keep | remove | correct
JAVA_MEM=2000m                  # Java max memory in MB

# Penalty costs
FSC=30                          # Frameshift cost
SC=100                          # Stop codon cost
GAPC=7                          # Gap creation cost
GAPE=1                          # Gap extension cost

# Less-reliable sequences penalty
FSC_LR=10
SC_LR=60

# Process each .fasta file
for file in *.fasta; do
    base=$(basename "$file" .fasta)

    echo "Aligning $file with MACSE..."

    java -Xmx"$JAVA_MEM" -jar "$MACSE_JAR" \
        -prog alignSequences \
        -seq "$file" \
        -out_NT "${ALIGN_DIR}/${base}_aligned_NT.fasta" \
        -out_AA "${ALIGN_DIR}/${base}_aligned_AA.fasta" \
        -thread "$THREADS" \
        -code "$GENETIC_CODE" \
        -shift "$FRAMESHIFT_MODE" \
        -fsc "$FSC" \
        -sc "$SC" \
        -gapc "$GAPC" \
        -gape "$GAPE" \
        -fsc_lr "$FSC_LR" \
        -sc_lr "$SC_LR" \
        -del yes

    echo "Alignment complete for: $file"
done

echo "All alignments stored in: $ALIGN_DIR"
```

**Core Configuration**
Parameter	Value	Explanation
-thread	4	Number of CPU threads for parallel processing.
-code	1	NCBI Genetic Code Table 1 (standard eukaryotic code). Other codes: 2 (vertebrate mitochondrial), 11 (bacterial), etc.
-shift	keep	How to handle frameshifts:
- keep: Preserve frameshifts as gaps (---).
- remove: Delete sequences with frameshifts.
- correct: Attempt to fix frameshifts.
-Xmx	2000m	Java heap memory (2GB). Increase (e.g., 4000m) for large datasets.

**Penalty Costs (Critical for Alignment Quality)**
A. Standard Penalties
Parameter	Value	Role in Alignment
-fsc	30	Frameshift cost: Penalty for introducing a frameshift (higher = fewer frameshifts).
-sc	100	Stop codon cost: Penalty for aligning a stop codon (*) in the middle of a sequence.
-gapc	7	Gap creation cost: Penalty for opening a gap (higher = fewer gaps).
-gape	1	Gap extension cost: Penalty for extending a gap (lower = longer gaps allowed).

B. Less-Reliable (LR) Sequence Penalties
Parameter	Value	Purpose
-fsc_lr	10	Lower frameshift penalty for unreliable sequences (e.g., low-quality data).
-sc_lr	60	Lower stop codon penalty for unreliable sequences.

    Standard penalties (-fsc, -sc) apply to high-confidence regions.
    LR penalties (-fsc_lr, -sc_lr) relax constraints for ambiguous regions, preventing over-penalization.

**Output Options**
Parameter	Output File	Content
-out_NT	*_aligned_NT.fasta	Nucleotide alignment (frameshifts as gaps).
-out_AA	*_aligned_AA.fasta	Amino acid alignment (stop codons as *).
-del	yes	Delete temporary files.

**Biological Interpretation of Penalties**
Frameshifts (-fsc)
    High cost (30): Favors alignments with fewer frameshifts.
    Example: A sequence with a frameshift mutation will be heavily penalized unless the mutation is biologically real.

Stop Codons (-sc)
    Very high cost (100): Strongly discourages internal stop codons (likely sequencing errors).
    Exception: Lower penalty (-sc_lr 60) for low-quality regions.

Gaps (-gapc, -gape)
    Gap creation (7) > extension (1): Favors fewer but longer gaps (common in evolutionary alignments).














*impMKT aligment:*

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
    echo "Example: $0 /path/to/fastas /path/to/speltoides_CDS.fasta /path/to/output"
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

**Merge the sequences to run MACSE:**
A Slurm batch script for concatenating multiple FASTA files into a single output file.
- Merges multiple FASTA files efficiently using `seqkit`

Usage:
```bash
sbatch merge_sequences.sbatch <input_directory> <output_filename>
```

merging_sequences.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=fasta_merger
#SBATCH --output=./%x_%j_out.txt
#SBATCH --error=./%x_%j_err.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G  # Reduced from 80G since we're using cat
#SBATCH --time=6:00:00  # Reduced time estimate
#SBATCH --partition=agap_normal

echo "========================================"
echo "Starting FASTA Merger (cat version)"
echo "Date: $(date)"
echo "========================================"

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "ERROR: Missing arguments!"
    echo "Usage: $0 <input_directory> <output_filename>"
    echo "Example: $0 /path/to/fastas merged.fasta"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_FILE="$2"

echo "Input directory: ${INPUT_DIR}"
echo "Output file: ${OUTPUT_FILE}"
echo ""

# Verify input directory
if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: Input directory does not exist!"
    exit 1
fi

# Count and verify input files
INPUT_COUNT=$(find "${INPUT_DIR}" -maxdepth 1 -name "*.fasta" | wc -l)
if [ "${INPUT_COUNT}" -eq 0 ]; then
    echo "ERROR: No .fasta files found in input directory!"
    exit 1
fi
echo "Found ${INPUT_COUNT} .fasta files to merge"
echo ""

# Merge using cat (memory efficient)
echo "Starting merge with cat..."
START_TIME=$(date +%s)

# Create empty output file first
> "${OUTPUT_FILE}"

# Process files one by one to avoid argument list too long errors
find "${INPUT_DIR}" -maxdepth 1 -name "*.fasta" -print0 | while IFS= read -r -d '' file; do
    echo "Adding ${file} to merge..."
    cat "${file}" >> "${OUTPUT_FILE}"
done

# Verify output
if [ $? -eq 0 ] && [ -s "${OUTPUT_FILE}" ]; then
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME-START_TIME))
    echo ""
    echo "SUCCESS: Merged ${INPUT_COUNT} files in ${RUNTIME} seconds!"
    echo "Output file: ${OUTPUT_FILE}"
    echo "Size: $(ls -lh "${OUTPUT_FILE}" | awk '{print $5}')"
    echo "Sequences: $(grep -c "^>" "${OUTPUT_FILE}")"
else
    echo ""
    echo "ERROR: Merge failed!"
    exit 1
fi

echo ""
echo "========================================"
echo "FASTA Merger completed"
echo "Date: $(date)"
echo "========================================"
```

