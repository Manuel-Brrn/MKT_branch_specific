
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


**Clean alignments with Hmmcleaner**
```bash
module load bioinfo-cirad
module load hmmer/3.3.2-singularity

# Define the input directory
INPUT_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/MACSE/speltoides_covered"

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
done > hmmcleaner_speltoides.log 2>&1
```

**Copy the cleaned alignment in the right directory**
```bash
cp *NT_hmm.fasta* /path/nucleotides_alignments_cleaned_hmm_cleaner/
```

**Put the outgroup sequence at the end of the alignment**
```bash
for f in *.fasta; do
    seqkit grep -v -r -p "H_vulgare" "$f" > tmp_non_outgroup.fasta
    seqkit grep -r -p "H_vulgare" "$f" > tmp_outgroup.fasta
    cat tmp_non_outgroup.fasta tmp_outgroup.fasta > "/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/MACSE/urartu_covered/nucleotides_alignments_cleaned_hmm_cleaner/reordered_alignments/${f%.fasta}_reordered.fasta"
done
rm tmp_non_outgroup.fasta tmp_outgroup.fasta
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
 















