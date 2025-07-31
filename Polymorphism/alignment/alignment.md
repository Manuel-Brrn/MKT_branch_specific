# Alignment:

**Research of orthologs with focal species and outgroup**

# Reciprocal BLAST Pipeline
Performs reciprocal BLASTn comparisons between two FASTA files, using the input filenames to label outputs.

# Usage
```bash
sbatch reciprocal_blast_speltoides_hordeum_cds.sbatch \
    /path/to/query.fasta \
    /path/to/subject.fasta \
    /path/to/output_directory
```
reciprocal_blast.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=reciprocal_blast
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=9G
#SBATCH --time=1:00:00
#SBATCH --partition=agap_short

echo "Starting reciprocal BLAST analysis at $(date)"

### Check command-line arguments
if [ "$#" -ne 3 ]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 <query_fasta> <subject_fasta> <output_dir>"
    exit 1
fi

QUERY="$1"
SUBJECT="$2"
OUTPUT_DIR="$3"

# Extract nicknames from filenames
QUERY_NICKNAME=$(basename "$QUERY" | cut -d'.' -f1)
SUBJECT_NICKNAME=$(basename "$SUBJECT" | cut -d'.' -f1)

echo "Input parameters:"
echo "- Query: $QUERY ($QUERY_NICKNAME)"
echo "- Subject: $SUBJECT ($SUBJECT_NICKNAME)"
echo "- Output directory: $OUTPUT_DIR"

### Setup environment
echo "Setting up Conda environment..."
source /nfs/work/agap_id-bin/img/snakemake/7.32.4/bin/activate base

if ! conda env list | grep -q "blast_env"; then
    echo "Creating blast_env..."
    conda create -n blast_env -y blast -c bioconda
fi

echo "Activating blast_env..."
conda activate blast_env

### Run reciprocal BLAST
echo "Running BLAST: $QUERY_NICKNAME vs $SUBJECT_NICKNAME..."
blastn -query "$QUERY" \
       -subject "$SUBJECT" \
       -evalue 1e-50 \
       -out "${OUTPUT_DIR}/${QUERY_NICKNAME}_vs_${SUBJECT_NICKNAME}.tab" \
       -outfmt 6 \
       -num_threads 8

echo "Running BLAST: $SUBJECT_NICKNAME vs $QUERY_NICKNAME..."
blastn -query "$SUBJECT" \
       -subject "$QUERY" \
       -evalue 1e-50 \
       -out "${OUTPUT_DIR}/${SUBJECT_NICKNAME}_vs_${QUERY_NICKNAME}.tab" \
       -outfmt 6 \
       -num_threads 8

### Verify outputs
if [ -s "${OUTPUT_DIR}/${QUERY_NICKNAME}_vs_${SUBJECT_NICKNAME}.tab" ] && \
   [ -s "${OUTPUT_DIR}/${SUBJECT_NICKNAME}_vs_${QUERY_NICKNAME}.tab" ]; then
    echo "Successfully completed reciprocal BLAST at $(date)"
    echo "Output files:"
    echo "- ${QUERY_NICKNAME}_vs_${SUBJECT_NICKNAME}.tab"
    echo "- ${SUBJECT_NICKNAME}_vs_${QUERY_NICKNAME}.tab"
else
    echo "Error: BLAST outputs not generated properly!" >&2
    exit 1
fi
```

# Python script to obtain RBH
# Usage
```bash
python3 reciprocal_best_hits.py focal_vs_out.tab out_vs_focal.tab 1 2 11 low RBH.tab
```

 reciprocal_best_hits.py
```bash
#!/usr/bin/env python
"""Reciprocal Best Hit (RBH) using BLAST style tabular input

Takes seven command line options,
1. Tabular filename of A against B
2. Tabular filename of B against A
3. Query ID column number (assumed to be same for both input files), e.g. c1
4. Match ID column number (assumed to be same for both input files), e.g. c2
5. Score column number (assumed to be same for both input files), e.g. c12
6. Want higest or lowest score? (use string high or low)
7. Output filename

"""
import sys


def stop_err( msg ):
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)

#Parse Command Line
try:
    a_vs_b, b_vs_a, c_query, c_match, c_score, sort_order, out_file = sys.argv[1:]
except:
    stop_err("Expect 7 arguments: two input files, column settings, output file")


want_highest = want_lowest = False
if sort_order == "high":
    want_highest = True
elif sort_order == "low":
    want_lowest = True
else:
    stop_err("Sort order argument should be high or low")

if out_file in [a_vs_b, b_vs_a]:
    stop_err("Output file would overwrite an input file")

def get_col_index(col_str):
    if col_str[0]=="c":
        col_str = col_str[1:]
    return int(col_str)-1

c_query = get_col_index(c_query)
c_match = get_col_index(c_match)
c_score = get_col_index(c_score)
if len(set([c_query, c_match, c_score])) < 3:
    stop_err("Need three different column numbers!")

"""
def load_best(filename, id1col, id2col):
    best = dict()
    for line in open(a_vs_b):
        if line.startswith("#"): continue
        parts = line.rstrip("\n").split("\t")
        id1 = parts[id1col]
        id2 = parts[id2col]
        score = float(parts[c_score])
        if (a not in best_a_vs_b) \
        or (want_highest and score > best[id1][1]) \
        or (want_lowest and score < best[id1][1]):
            best[id1] = (id2, score)
    return best
best_a_vs_b = load_best(a_vs_b, c_query, c_match)
"""

best_a_vs_b = dict()
for line in open(a_vs_b):
    if line.startswith("#"): continue
    parts = line.rstrip("\n").split("\t")
    a = parts[c_query]
    b = parts[c_match]
    score = float(parts[c_score])
    if (a not in best_a_vs_b) \
    or (want_highest and score > best_a_vs_b[a][1]) \
    or (want_lowest and score < best_a_vs_b[a][1]):
        best_a_vs_b[a] = (b, score, parts[c_score])
b_short_list = set(b for (b,score, score_str) in best_a_vs_b.values())

best_b_vs_a = dict()
for line in open(b_vs_a):
    if line.startswith("#"): continue
    parts = line.rstrip("\n").split("\t")
    b = parts[c_query]
    a = parts[c_match]
    if a not in best_a_vs_b:
        continue
        #stop_err("The A-vs-B file does not have A-ID %r found in B-vs-A file" % a)
    if b not in b_short_list: continue
    score = float(parts[c_score])
    if (b not in best_b_vs_a) \
    or (want_highest and score > best_b_vs_a[b][1]) \
    or (want_lowest and score < best_b_vs_a[b][1]):
        best_b_vs_a[b] = (a, score, parts[c_score])
#TODO - Preserve order from A vs B?
a_short_list = sorted(set(a for (a,score,score_str) in best_b_vs_a.values()))

count = 0
outfile = open(out_file, 'w')
outfile.write("#A_id\tB_id\tA_vs_B\tB_vs_A\n")
for a in a_short_list:
    b = best_a_vs_b[a][0]
    if a == best_b_vs_a[b][0]:
        outfile.write("%s\t%s\t%s\t%s\n" % (a, b, best_a_vs_b[a][2], best_b_vs_a[b][2]))
        count += 1
outfile.close()
print ("Done, %i RBH found" % count)
```

**Processing of the RBH.tab (Tr.urartu - Ho.spontaneum)**
```bash
awk -F'\t' -v OFS='\t' 'NR>1 {gsub("T_urartu\\|","",$1); gsub("\\.","_",$1); print}' RBH.tab > RBH_processed.tab
```

**Alignment**

*dNdSpiNpiS alignment*
MACSE Alignment Pipeline

## Description
This pipeline performs pairwise alignment of ORF sequences using MACSE with reciprocal best hits (RBH) as input.

## Usage
```bash
 sbatch alignment.sbatch \
    <ORF_alignment_1.fasta> \
    <ORF_alignment_2.fasta> \
    <RBH_table.tab> \
    <output_directory>
```

#Arguments
Parameter           Description                 Example
ORF_alignment_1     First set of ORF sequences  speltoides_merged_bams.cds
ORF_alignment_2     Second set of ORF sequences hordeum_full_CDS_modified.fasta
RBH_table           Reciprocal best hits table  RBH_filtered_modified.tab
output_directory    Directory for results       ./alignment_results/


 Macse_alignment.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=alignment_MACSE
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=50:00:00
#SBATCH --partition=agap_long

echo "Starting alignment job at $(date)"
echo "Job ID: $SLURM_JOB_ID"

### Check command-line arguments
if [ "$#" -ne 4 ]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 <ORF_alignment_1> <ORF_alignment_2> <RBH_tab> <output_dir>"
    exit 1
fi

ORF_ALIGNMENT_1="$1"
ORF_ALIGNMENT_2="$2"
RBH_TAB="$3"
OUTPUT_DIR="$4"

echo "Input parameters:"
echo "- ORF alignment 1: $ORF_ALIGNMENT_1"
echo "- ORF alignment 2: $ORF_ALIGNMENT_2"
echo "- RBH table: $RBH_TAB"
echo "- Output directory: $OUTPUT_DIR"

### Setup working directory
echo "Setting up working directory..."
cd "$OUTPUT_DIR" || { echo "Error: Failed to change to output directory"; exit 1; }

### Load modules
echo "Loading required modules..."
module load bioinfo-cirad
module load snakemake/7.32.4-conda

### Conda environment setup
echo "Setting up Conda environment..."
source /nfs/work/agap_id-bin/img/snakemake/7.32.4/etc/profile.d/conda.sh
conda activate base

### Check/create environment
echo "Checking for bioperl_env..."
if ! conda env list | grep -q "bioperl_env"; then
    echo "Creating bioperl_env..."
    conda create -y -n bioperl_env perl-bioperl -c bioconda
fi

echo "Activating bioperl_env..."
conda activate bioperl_env

### Perl environment setup
echo "Setting up Perl environment..."
if ! command -v cpanm &> /dev/null; then
    echo "Installing cpanm..."
    curl -L https://cpanmin.us | perl - App::cpanminus
fi

if [[ ! "$PATH" =~ "$HOME/perl5/bin" ]]; then
    export PATH="$HOME/perl5/bin:$PATH"
fi

if ! perl -MSwitch -e 'print "Switch module is installed\n"' 2>/dev/null; then
    echo "Installing Switch module..."
    cpanm Switch
fi

### Verify BioPerl
echo "Verifying BioPerl installation..."
perl -MBio::SeqIO -e 'print "BioPerl works!\n"' || {
    echo "Error: BioPerl test failed!" >&2
    exit 1
}

### Set Perl library path
#export PERL5LIB="/storage/simple/users/barrientosm/.conda/envs/bioperl_env/lib/perl5/site_perl:$HOME/perl5/lib/perl5:$PERL5LIB"
export PERL5LIB="/home/barrientosm/my_perl_5.36.0/lib/site_perl/5.36.0:/home/barrientosm/my_perl_5.36.0/lib/site_perl/5.36.0/x86"
echo "PERL5LIB set to: $PERL5LIB"

### Run alignment
echo "Starting MACSE alignment..."
echo "Command: perl AlignTwoProfiles_wrapper_macse_v1_2_modified.pl -clean yes -p1 $ORF_ALIGNMENT_1 -p2 $ORF_ALIGNMENT_2 -rbh $RBH_TAB"

/storage/simple/users/barrientosm/.conda/envs/bioperl_env/bin/perl \
    /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/MACSE/AlignTwoProfiles_wrapper_macse_v1_2_modified.pl \
    -clean yes \
    -p1 "$ORF_ALIGNMENT_1" \
    -p2 "$ORF_ALIGNMENT_2" \
    -rbh "$RBH_TAB"

### Check exit status
if [ $? -eq 0 ]; then
    echo "Alignment completed successfully at $(date)"
    echo "Results saved in: $OUTPUT_DIR"
else
    echo "Error: Alignment failed!" >&2
    exit 1
fi
```

*impMKT aligment:*

Add the reference sequence of the focal species:
Extracts sequences from a FASTA file based on contig IDs using parallel processing.

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




