



**impMKT aligment:**

**Add the reference sequence of the focal species:**

**Extracts sequences from a FASTA file based on contig IDs**
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


**Alignment with Macse:**
Adapt the alignment output directory, the number of files to aligned and run the script on the input directory

```bash
#!/bin/bash
#SBATCH --job-name=alignment_impMKT_species
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --array=1-5630%17  # Adjust range based on number of files
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=50:00:00
#SBATCH --partition=agap_long

# ---------------------------
# Load necessary modules
# ---------------------------
module load bioinfo-cirad
module load java/jre1.8.0_31

# ---------------------------
# MACSE Configuration
# ---------------------------
MACSE_JAR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/03_scripts/dn_ds_pipeline/MACSE/macse_v1.2.jar"
ALIGN_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/MACSE/urartu_covered_impMKT"
mkdir -p "$ALIGN_DIR"

THREADS=4
GENETIC_CODE=1
FRAMESHIFT_MODE="keep"
JAVA_MEM=4000m

FSC=30
SC=100
GAPC=7
GAPE=1
FSC_LR=10
SC_LR=60

# ---------------------------
# Select the correct file based on SLURM_ARRAY_TASK_ID
# ---------------------------
FILES=(*.fasta)
TOTAL=${#FILES[@]}
INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if [ "$INDEX" -ge "$TOTAL" ]; then
    echo "SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds number of files ($TOTAL). Exiting."
    exit 1
fi

FILE="${FILES[$INDEX]}"
BASENAME=$(basename "$FILE" .fasta)

echo "SLURM Task ID: $SLURM_ARRAY_TASK_ID  Processing file: $FILE"

# ---------------------------
# Run MACSE alignment
# ---------------------------
java -Xmx"$JAVA_MEM" -jar "$MACSE_JAR" \
    -prog alignSequences \
    -seq "$FILE" \
    -out_NT "${ALIGN_DIR}/${BASENAME}_aligned_NT.fasta" \
    -out_AA "${ALIGN_DIR}/${BASENAME}_aligned_AA.fasta" \
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

echo "Alignment complete: $FILE $ALIGN_DIR/${BASENAME}_aligned_NT.fasta"
```
Core Configuration Parameter Value Explanation -thread 4 Number of CPU threads for parallel processing. -code 1 NCBI Genetic Code Table 1 (standard eukaryotic code). Other codes: 2 (vertebrate mitochondrial), 11 (bacterial), etc. -shift keep How to handle frameshifts:

    keep: Preserve frameshifts as gaps (---).
    remove: Delete sequences with frameshifts.
    correct: Attempt to fix frameshifts. -Xmx 2000m Java heap memory (2GB). Increase (e.g., 4000m) for large datasets.

Penalty Costs (Critical for Alignment Quality) A. Standard Penalties Parameter Value Role in Alignment -fsc 30 Frameshift cost: Penalty for introducing a frameshift (higher = fewer frameshifts). -sc 100 Stop codon cost: Penalty for aligning a stop codon (*) in the middle of a sequence. -gapc 7 Gap creation cost: Penalty for opening a gap (higher = fewer gaps). -gape 1 Gap extension cost: Penalty for extending a gap (lower = longer gaps allowed).

B. Less-Reliable (LR) Sequence Penalties Parameter Value Purpose -fsc_lr 10 Lower frameshift penalty for unreliable sequences (e.g., low-quality data). -sc_lr 60 Lower stop codon penalty for unreliable sequences.

Standard penalties (-fsc, -sc) apply to high-confidence regions.
LR penalties (-fsc_lr, -sc_lr) relax constraints for ambiguous regions, preventing over-penalization.

Output Options Parameter Output File Content -out_NT *_aligned_NT.fasta Nucleotide alignment (frameshifts as gaps). -out_AA *_aligned_AA.fasta Amino acid alignment (stop codons as *). -del yes Delete temporary files.

Biological Interpretation of Penalties Frameshifts (-fsc) High cost (30): Favors alignments with fewer frameshifts. Example: A sequence with a frameshift mutation will be heavily penalized unless the mutation is biologically real.

Stop Codons (-sc) Very high cost (100): Strongly discourages internal stop codons (likely sequencing errors). Exception: Lower penalty (-sc_lr 60) for low-quality regions.

Gaps (-gapc, -gape) Gap creation (7) > extension (1): Favors fewer but longer gaps (common in evolutionary alignments).

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

**Copy the cleaned alignment in the right directory**
```bash
cp *NT_hmm.fasta* /path/nucleotides_alignments_cleaned_hmm_cleaner/
```

**Ordonner les séquences avec la référence en première et l'outgroup en dernier:**
```bash
for f in *.fasta; do
    echo "Traitement de $f..."
    
    # 1. Extraire T_urartu et le placer en premier (si présent)
    seqkit grep -r -p "T_urartu" "$f" > tmp_first.fasta 2>/dev/null
    
    # 2. Extraire toutes les autres séquences (sauf H_vulgare et T_urartu)
    seqkit grep -v -r -p "H_vulgare" "$f" | seqkit grep -v -r -p "T_urartu" > tmp_rest.fasta 2>/dev/null
    
    # 3. Extraire H_vulgare pour le mettre à la fin
    seqkit grep -r -p "H_vulgare" "$f" > tmp_outgroup.fasta 2>/dev/null
    
    # 4. Fusionner dans l'ordre: T_urartu -> autres séquences -> H_vulgare
    cat tmp_first.fasta tmp_rest.fasta tmp_outgroup.fasta > "/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/MACSE/urartu_covered/nucleotides_alignments_cleaned_hmm_cleaner/reordered_alignments/${f%.fasta}_reordered.fasta"
    
    # 5. Nettoyage
    rm -f tmp_first.fasta tmp_rest.fasta tmp_outgroup.fasta
done
```


**Rename_sequences for dndspinpis**
```bash
#!/bin/bash

# =============================================
# CONFIGURABLE PARAMETERS (edit as needed)
# =============================================

# Species tags
OUTGROUP_OLD_TAG="H_vulgare"    # Original outgroup header pattern (keep original tag)
FOCAL_OLD_TAG="sp"              # Original focal species tag (to replace)
FOCAL_NEW_TAG="T_urartu"        # New focal species identifier

# RBH file settings
RBH_QUERY_COL=2                 # Column with outgroup IDs in RBH file
RBH_TARGET_COL=1                # Column with focal species IDs in RBH file

# =============================================
# FUNCTIONS (modular steps)
# =============================================

clean_headers() {
  local infile="$1"
  local outfile="$2"
  echo "Cleaning headers: removing text after first space..."
  awk '/^>/ { split($0, a, " "); print a[1]; next } { print }' "$infile" > "$outfile"
}

tag_outgroup() {
  local infile="$1"
  local outfile="$2"
  echo "Keeping original outgroup tag (${OUTGROUP_OLD_TAG}) and adding '|outgroup|' suffix..."
  awk -v oldtag="$OUTGROUP_OLD_TAG" '
    /^>/ {
      if ($0 ~ oldtag) {
        if ($0 !~ /\|outgroup\|$/) {
          print $0 "|outgroup|";
        } else {
          print $0;
        }
      } else {
        print $0;
      }
      next
    }
    { print }
  ' "$infile" > "$outfile"
}

map_rbh_contigs() {
  local infile="$1"
  local rbhfile="$2"
  local outfile="$3"
  echo "Mapping outgroup to focal contigs via RBH file..."
  awk -v rbhfile="$rbhfile" \
     -v query_col="$RBH_QUERY_COL" \
     -v target_col="$RBH_TARGET_COL" \
     -v outgroup_tag="$OUTGROUP_OLD_TAG" '
    BEGIN {
      FS="\t";
      while ((getline line < rbhfile) > 0) {
        if (line ~ /^#/) continue;
        split(line, fields, "\t");
        original_id = fields[query_col];
        b_to_a[original_id] = fields[target_col];
        modified_id = original_id;
        sub(/\.t/, "_t", modified_id);
        if (modified_id != original_id) {
          b_to_a[modified_id] = fields[target_col];
        }
      }
      close(rbhfile);
    }
    {
      if ($0 ~ /^>/ && $0 ~ outgroup_tag) {
        header = substr($0, 2);
        sub(/\|outgroup\|$/, "", header);
        split(header, arr, "|");
        hordeum_id = arr[1] "|" arr[2];
        if (hordeum_id in b_to_a) {
          print ">" b_to_a[hordeum_id] "|" header "|outgroup|";
        } else {
          print "WARNING: No RBH match for " hordeum_id > "/dev/stderr";
          print $0;
        }
      } else {
        print $0;
      }
    }' "$infile" > "$outfile"
}

replace_species_tag() {
  local infile="$1"
  local outfile="$2"
  echo "Replacing species tag (${FOCAL_OLD_TAG} to ${FOCAL_NEW_TAG})..."
  sed -E "s/^([^|]+\|)${FOCAL_OLD_TAG}(\|.*)/\1${FOCAL_NEW_TAG}\2/" "$infile" > "$outfile"
}

remove_trailing_pipes() {
  local infile="$1"
  local outfile="$2"
  echo "Removing trailing '|' characters..."
  sed -E '/^>/{ s/\|$// }' "$infile" > "$outfile"
}

# =============================================
# MAIN EXECUTION
# =============================================

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <alignment_fasta> <rbh_tab_file> <output_fasta>"
  exit 1
fi

for file in "$1" "$2"; do
  if [ ! -f "$file" ]; then
    echo "Error: File '$file' not found!" >&2
    exit 1
  fi
done

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

TMP1="$TMP_DIR/clean_headers.fasta"
TMP2="$TMP_DIR/outgroup_tagged.fasta"
TMP3="$TMP_DIR/rbh_mapped.fasta"
TMP4="$TMP_DIR/species_replaced.fasta"

clean_headers "$1" "$TMP1"
tag_outgroup "$TMP1" "$TMP2"
map_rbh_contigs "$TMP2" "$2" "$TMP3"
replace_species_tag "$TMP3" "$TMP4"
remove_trailing_pipes "$TMP4" "$3"

echo -e "\nAll steps completed. Output saved to: $3"
```

*Run the script*
```bash
OUTDIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/MACSE/urartu_covered/nucleotides_alignments_cleaned_hmm_cleaner/reordered_alignments/renamed_header"

mkdir -p "$OUTDIR"

for f in *.fasta; do
  ./reformat_header_alignment_file.sh "$f" /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/reciprocal_blast/urartu_hordeum/urartu_hordeum_vespa_cleaned/RBH_processed.tab "$OUTDIR/${f%.fasta}_renamed.fasta"
done
```
Must edit ; 
OUTGROUP_OLD_TAG=".."    # Original outgroup header pattern
OUTGROUP_NEW_TAG=.."         # New outgroup identifier
FOCAL_OLD_TAG="sp"                    # Original focal species tag (to replace)
FOCAL_NEW_TAG=".."              # New focal species identifier
OUTDIR path
RBH table path

**Delete empty alignments:**
```bash
find . -name "*.fasta" -type f -size 0 -delete
```

**Stastistic with AMAS:**
```bash
pip install amas
find ~/.local/lib/python3.9/site-packages -name "AMAS.py"
python3 /home/barrientosm/.local/lib/python3.9/site-packages/amas/AMAS.py summary -f fasta -d dna -i *.fasta
```
 
**Merge all alignments**
```bash
cat *.fasta > /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/dNdSpiNpiS/monococcum_covered/dNdSpiNpiS_input/monococcum_alignment.fasta
```
Must edit the path and name of the file

**Replace all x with N (unknown nucleotide)**
```bash
sed '/^>/! s/x/N/g' monococcum_alignment.fasta > monococcum_alignment_clean.fasta
```
Must edit the name's file


**Obtain daf and sfs from alignment**

```py

import sys
import time
import numpy as np
import pandas as pd
import pyfaidx as px
import pybedtools

def degenerancy(data,codonDict):

        #DEGENERANCY DICTIONARIES
        standardDict = {
                'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202',
                'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004',
                'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002',
                'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000',
                'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204',
                'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004',
                'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002',
                'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204',
                'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000',
                'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004',
                'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002',
                'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202',
                'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004',
                'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004',
                'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002',
                'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}

        if(codonDict == 'standard'):
                degenerateCodonTable = standardDict

        degenerancy = ''
        for i in range(0, len(data),3):
                codon = data[i:i+3]
                if('N' in codon or '-' in codon):
                        degenerancy += codon
                else:
                        degenerancy += degenerateCodonTable[codon]

        return(degenerancy)
def sequencesToMatrix(multiFasta,split=None):

        # Extract samples from fastas
        samples = list(multiFasta.keys())

        if(split is None):
                seqLen = len(multiFasta[samples[0]][:].seq)

                if((seqLen % 3) != 0):
                        print('cdsLength')
                        sys.exit('cdsLength')

                # Create empty array with ndimesions equal to multi-Fasta lines and length
                matrix = np.empty([len(samples),len(multiFasta[samples[0]][:].seq)],dtype='str')

                # List to append indexes if whole sequence at any population is len(seq) * 'N'
                deleteIndex = list()

                # Iter fasta to add sequence to matrix
                for i in range(1,len(samples),1):
                        # Extract each sample sequence
                        tmp = multiFasta[samples[i]][:].seq
                        if(len(tmp) != seqLen):
                                print('errorAlign')
                                sys.exit('errorAlign')
                        if('N' in tmp):
                                deleteIndex.append(i)
                        else:
                                matrix[i] = list(tmp)

                # Delete lines
                matrix = np.delete(matrix,deleteIndex,0)

                degenCode = degenerancy(multiFasta[samples[0]][:].seq,codonTable)
                # Put degenerancy in first ndarray element
                matrix[0] = list(degenCode)
                # NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer

                matrix = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')

                return(matrix)

        else:
                # Create empty array with ndimesions equal to multi-Fasta lines and length
                split = [split[0]+1,split[1]+1]
                matrix = np.empty([len(samples),(split[1]-split[0]+1)],dtype='str')

                # List to append indexes if whole sequence at any population is len(seq) * 'N'
                deleteIndex = list()

                # Iter fasta to add sequence to matrix
                for i in range(1,len(samples),1):
                        # Extract each sample sequence
                        tmp = multiFasta.get_spliced_seq(samples[i], [split]).seq
                        if(tmp == ('N' * len(tmp))):
                                deleteIndex.append(i)
                        else:
                                matrix[i] = list(tmp)

                # Delete lines
                matrix = np.delete(matrix,deleteIndex,0)

                degenCode = degenerancy( multiFasta.get_spliced_seq(samples[i], [split]).seq.upper(),codonTable)
                # Put degenerancy in first ndarray element
                matrix[0] = list(degenCode)
                # NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer

                matrix = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')

                return(matrix)

def uSfsFromFasta(sequenceMatrix):
        output = list()
        for x in np.nditer(sequenceMatrix, order='F',flags=['external_loop']):
                degen = x[0]
                AA = x[-1]

                # Undefined Ancestra Allele. Try to clean out of the loop
                if(AA == 'N' or AA == '-'):
                        next
                elif('N' in x[1:-1] or '-' in x[1:-1]):
                        next
                # Monomorphic sites. Try to clean out of the loop
                elif(np.unique(x[1:][np.where(x[1:]!='N')]).shape[0] == 1):
                        next
                else:
                        pol = x[1:-1]
                        if(degen == '4'):
                                functionalClass = '4fold'
                        else:
                                functionalClass = '0fold'

                        # Check if pol != AA and monomorphic
                        if((np.unique(pol).shape[0] == 1) and (np.unique(pol)[0] != AA)):
                                div = 1; AF = 0
                                tmp = [AF,div,functionalClass]
                                output.append(tmp)
                        else:
                                AN = x[1:-1].shape[0]
                                AC = pd.DataFrame(data=np.unique(x[1:-1], return_counts=True)[1],index=np.unique(x[1:-1], return_counts=True)[0])
                                div = 0
                                if(AA not in AC.index):
                                        next
                                else:
                                        AC = AC[AC.index!=AA]
                                        if(len(AC) == 0):
                                                next
                                        else:
                                                AF = AC.iloc[0]/AN
                                                AF = AF.iloc[0]
                                tmp = [AF,div,functionalClass]
                                output.append(tmp)
        return(output)
def formatSfs(sequenceMatrix,rawSfsOutput,dafFile,divFile,path,append=True):

        df = pd.DataFrame(rawSfsOutput)
        df['id'] = 'uploaded'
        df.columns = ['derivedAlleleFrequency','d','functionalClass','id']

        # Extract divergence data
        div = df[['id','functionalClass','d']]
        div = div[div['d']!=0]
        div = div.groupby(['id','functionalClass'])['d'].count().reset_index()
        div = div.pivot_table(index=['id'],columns=['functionalClass'],values='d').reset_index()
        try:
                div = div[['0fold','4fold']]
        except:
                if('4fold' in div.columns):
                        div = div[['4fold']]
                        div['0fold'] = 0
                elif('0fold' in div.columns):
                        div = div[['0fold']]
                        div['4fold'] = 0

        div['mi'] =  sequenceMatrix[0][np.where(sequenceMatrix[0]=='0')].shape[0]
        div['m0'] =  sequenceMatrix[0][np.where(sequenceMatrix[0]=='4')].shape[0]
        div.columns = ['Di','D0','mi','m0']
        # div = div.pivot_table(index=['functionalClass'],columns=['functionalClass'],values='div').reset_index()

        # Create SFS pd.DataFrame by functionClass and 20 frequency bin
        daf = df[df['d']!=1][['derivedAlleleFrequency','functionalClass','id']]
        bins = np.arange(0.025,1.05,0.05)
        labels = np.arange(0.025,1.0,0.05).tolist()
        daf['categories'] = pd.cut(daf['derivedAlleleFrequency'],bins=bins,labels=labels)
        daf = daf.groupby(['functionalClass','id','categories']).count().reset_index()
        sfs = pd.DataFrame({'daf':daf['categories'].unique(),'P0':daf[daf['functionalClass']=='4fold']['derivedAlleleFrequency'].reset_index(drop=True),'Pi':daf[daf['functionalClass']=='0fold']['derivedAlleleFrequency'].reset_index(drop=True)})

        sfs = sfs[['daf','P0','Pi']]
        sfs['P0'] = sfs['P0'].fillna(0)
        sfs['Pi'] = sfs['Pi'].fillna(0)
        sfs['daf'] = sfs['daf'].apply(lambda x: round(x,3))

        if(append is True):
                sfs.to_csv(path + dafFile,sep='\t',header=True,index=False,mode='a')
                div.to_csv(path + divFile,sep='\t',header=True,index=False,mode='a')
        else:
                sfs.to_csv(path + dafFile,sep='\t',header=True,index=False)
                div.to_csv(path + divFile,sep='\t',header=True,index=False)

start = time.time()
analysis=sys.argv[1]
multiFasta = sys.argv[2]
dafFile = sys.argv[3]
divFile = sys.argv[4]
codonTable = sys.argv[5]


# Open multi-Fasta
file = px.Fasta('/home/barrientosm/scratch/impMKT/' + multiFasta,duplicate_action='first',sequence_always_upper=True,read_long_names=True)

# Create ndarray with sequences
multiFastaMatrix = sequencesToMatrix(file)

if(multiFastaMatrix.shape[0] < 4):
        print('numberOfLines')
        sys.exit('numberOfLines')

# Estimating SFS
rawSfs = uSfsFromFasta(multiFastaMatrix)

formatSfs(multiFastaMatrix,rawSfs,dafFile,divFile,'/home/barrientosm/scratch/impMKT/')
```

**Run it for all genes**
```bash
#!/bin/bash
#SBATCH --job-name=sfs_from_fasta
#SBATCH --output=sfs_%A_%a.out
#SBATCH --error=sfs_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --array=1-N%10   # Remplace N par le nb de fasta
#SBATCH --time=48:00:00
#SBATCH --partition=agap_long

# === User-defined variables ===
indir="/path/to/input/fasta_dir"
outdir="/path/to/output_dir"
python_script="/path/to/sfsFromFasta.py"

# Load required modules
module load bioinfo-cirad
module load python/3.9.6

# Create output directory
mkdir -p "$outdir"

# Get list of .fasta files
fasta_files=($indir/*.fasta)
total_files=${#fasta_files[@]}

# Check if array ID is within bounds
if [ "$SLURM_ARRAY_TASK_ID" -le "$total_files" ]; then
  aln="${fasta_files[$SLURM_ARRAY_TASK_ID - 1]}"
  full=$(basename "$aln" .fasta)
  base=$(echo "$full" | sed 's/_NT_hmm_reordered_renamed//')

  echo "[$SLURM_ARRAY_TASK_ID/$total_files] Processing $aln"
  python "$python_script" "$base" "$aln" "${outdir}/${base}_daf.tsv" "${outdir}/${base}_div.tsv" standard
else
  echo "Array task $SLURM_ARRAY_TASK_ID exceeds number of files ($total_files). Exiting."
fi
```


**Merge results for all genes**
```py
import pandas as pd
from glob import glob

# --- SFS GLOBAL ---
sfs_files = sorted(glob("sfs_output*.tsv"))
sfs_global = None

for f in sfs_files:
    df = pd.read_csv(f, sep="\t")
    if sfs_global is None:
        sfs_global = df.copy()
    else:
        sfs_global["P0"] += df["P0"]
        sfs_global["Pi"] += df["Pi"]

sfs_global.to_csv("SFS_global.tsv", sep="\t", index=False)

# --- DIV GLOBAL ---
div_files = sorted(glob("div_output*.tsv"))
div_global = None

for f in div_files:
    df = pd.read_csv(f, sep="\t")
    if div_global is None:
        div_global = df.copy()
    else:
        div_global["Di"] += df["Di"]
        div_global["D0"] += df["D0"]
        div_global["mi"] += df["mi"]
        div_global["m0"] += df["m0"]

div_global.to_csv("DIV_global.tsv", sep="\t", index=False)


```



