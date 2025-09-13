----------------------------------------------------------------------------------
|||||        Phase 3: Alignment Assessment & Phylogeny Reconstruction        |||||
----------------------------------------------------------------------------------
Details: This phase combines multiple third-party programs to automate the
assessment, selection, and phylogenetic reconstruction of protein MSAs.

Commands: metal_compare, prottest_setup, prottest_reader, mrbayes_setup
----------------------------------------------------------------------------------

*Alignment Assessment*

**Alignment**
```bash
#!/bin/bash
#SBATCH --job-name=mafft_align
#SBATCH --output=./macse_%j_out.txt
#SBATCH --error=./macse_%j_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --partition=agap_normal

INPUT_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/Single_Copy_Orthologue_Sequences"
OUTPUT_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/mafft"
THREADS=4  # Adjust

mkdir -p $OUTPUT_DIR

for FILE in $INPUT_DIR/N0.HOG*.fa; do
    BASE=$(basename $FILE .fa)
    mafft --auto --thread $THREADS $FILE > $OUTPUT_DIR/${BASE}_aln.fa 2>> mafft_errors.log
done

echo "Alignments completed. Results in $OUTPUT_DIR/"
```

**HmmCleaner**
```bash
cd /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/Bio-MUST-Apps-HmmCleaner-0.243280
module load hmmer/3.3.2-singularity
cpanm Bio::MUST::Apps::HmmCleaner
/home/barrientosm/my_perl_5.36.0/bin/cpanm -l ~/my_perl_5.36.0 \
    namespace::autoclean
/home/barrientosm/my_perl_5.36.0/bin/cpanm -l ~/my_perl_5.36.0 \
    PadWalker
/home/barrientosm/my_perl_5.36.0/bin/cpanm -l ~/my_perl_5.36.0 \
    Bio::MUST::Drivers

HmmCleaner.pl --version

# Set Perl environment
export PERL5LIB=/home/barrientosm/my_perl_5.36.0/lib/perl5:/home/barrientosm/my_perl_5.36.0/lib/perl5/x86_64-linux
export PATH=/home/barrientosm/my_perl_5.36.0/bin:$PATH

# Define the input directory
INPUT_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/MACSE/urartu_covered"

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
done > hmmcleaner_all.log 2>&1
```

**Renames orthologs**
```bash
for file in N0.*.fasta; do
    mv "$file" "${file#N0.}"
done
```

**Translate again in nucleotides**
```bash
#!/bin/bash
#SBATCH --job-name=translate
#SBATCH --output=log_%j_%x.out
#SBATCH --error=log_%j_%x.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --partition=agap_short

module load bioinfo-cirad
module load python/2.7.18

# === CONFIGURATION ===
INPUT_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/cleaned_alignments"
VESPA_SCRIPT="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/cleaned_alignments/vespa.py"
CDS_DATABASE="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/cds_sequences/database/cds_database.fasta"
OUTPUT_BASE="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/translated_alignments"

mkdir -p "$OUTPUT_BASE"
cd "$INPUT_DIR" || { echo "Erreur : impossible d'accéder à $INPUT_DIR"; exit 1; }

for fasta_file in *hmm.fasta; do
    [ ! -f "$fasta_file" ] && continue

    hog_id=$(echo "$fasta_file" | grep -o 'HOG[0-9]\+')
    hog_output_dir="$OUTPUT_BASE/$hog_id"
    mkdir -p "$hog_output_dir"

    # Nom temporaire sans extension
    tmp_name="$hog_id"
    cp "$fasta_file" "$hog_output_dir/$tmp_name"

    # Créer le répertoire attendu par vespa.py
    mkdir -p "$hog_output_dir/Map_Gaps_${hog_id}"

    echo "Traitement de $hog_id"

    (
        cd "$hog_output_dir" && \
        python2 "$VESPA_SCRIPT" map_alignments \
            -input="$tmp_name" \
            -database="$CDS_DATABASE"
    )

    if [ -f "$hog_output_dir/${hog_id}_aln_hmm.fasta" ]; then
        echo "OK pour $hog_id"
    else
        echo "ÉCHEC pour $hog_id"
    fi
done

echo "Traitement terminé. Résultats dans : $OUTPUT_BASE"
```


**Number of alignment lost**
```bash
grep -B 2 "Sequence" log_28749510_translate.out | grep "Traitement de" | awk '{print $3}' > extracted_hog_ids.txt
```

**Copy cleaned alignments translated in nucletoides**

```bash
#!/bin/bash

# Source and destination directories
SRC_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/translated_alignments"
DEST_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/translated_alignments_cleaned"
EXCLUDE_LIST="extracted_hog_ids.txt"

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Load the exclusion list into an associative array for fast lookup
declare -A EXCLUDE
while read -r line; do
    EXCLUDE["$line"]=1
done < "$EXCLUDE_LIST"

# Loop through all HOG directories in the source
for hog_dir in "$SRC_DIR"/*; do
    HOG_ID=$(basename "$hog_dir")

    # Skip if the HOG is in the exclusion list
    if [[ ${EXCLUDE["$HOG_ID"]+_} ]]; then
        echo "Skipping excluded HOG: $HOG_ID"
        continue
    fi

    # Define the path to the alignment file inside Map_Gaps_HOGXXXX
    FILE_TO_COPY="${hog_dir}/Map_Gaps_${HOG_ID}/${HOG_ID}"

    # Check if the file exists
    if [[ -f "$FILE_TO_COPY" ]]; then
        echo "Copying $FILE_TO_COPY to $DEST_DIR"
        cp "$FILE_TO_COPY" "$DEST_DIR/${HOG_ID}.fasta"
    else
        echo "File not found: $FILE_TO_COPY"
    fi
done
```

**Directory for codemld set_up**

**Gene tree inference function**
The -infer_genetree function is designed to automate the creation of the corresponding gene tree for a user-specified MSA. This is achieved by associating the taxa specified on a user-defined species tree with the headers created by ‘label_filename’ and infer_ensembl_species within the MSA. The function operates by first creating a copy of the species tree with the species names. The species tree is designated using the required species_tree option. The species names are then replaced with their associated MSA headers. If any species names remain after this phase, the taxa and their respective branches are removed from the tree to create the finished gene tree. It should be noted that the infer_genetree function incorporates the non-standard python library dendropy [Sukumaran et al., 2010].
The infer_genetree function in VESPA is mapping existing MSA sequences onto a predefined species tree through a process called tree reconciliation. 
Core Functionality
    Inputs:
        A multiple sequence alignment (MSA) with headers formatted by label_filename/infer_ensembl_species
        A user-provided species tree (divergence relationships of species/taxa)

    Output: A gene tree topology constrained by the species tree.

**Species tree**

Orthofinder repertory:
/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/cds_sequences/database/orthofinder_input/OrthoFinder/Results_Jul23_1/Species_Tree

 SpeciesTree_rooted.txt
(translated_hordeum_renamed:0.5,(translated_urartu_renamed:1,translated_monococcum_renamed:1)1:0.

```bash
sed -E 's/:[0-9.]+//g; s/[0-9]+//g' SpeciesTree_rooted.txt > SpeciesTree_topology_only.txt

#### remove problematic names
sed -i -e 's/Se_cereale_cds/S_cereale_cds/g' \
       -e 's/Tr_/T_/g' Simplified_SpeciesTree.txt
```



**Run infer_genetree**
```bash
#!/bin/bash

# Set species tree file
SPECIES_TREE="Simplified_SpeciesTree.txt"

# Loop through all fasta files in the directory
for fasta in *.fasta; do
    echo "Running vespa.py on $fasta..."
    python2 /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/03_scripts/dn_ds_pipeline/VESPA/VESPA-1.0.1/vespa.py infer_genetree -input="$fasta" -species_tree="$SPECIES_TREE"
done
```






