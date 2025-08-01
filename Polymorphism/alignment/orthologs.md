
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
