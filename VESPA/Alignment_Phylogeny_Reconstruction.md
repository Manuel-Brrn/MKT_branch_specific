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

export PERL5LIB=/home/barrientosm/my_perl_5.36.0/lib/site_perl/5.36.0:/home/barrientosm/my_perl_5.36.0/lib/site_perl/5.36.0/x86_64-linux:/home/barrientosm/my_perl_5.36.0/lib/5.36.0:/home/barrientosm/my_perl_5.36.0/lib/5.36.0/x86_64-linux
export PATH=/home/barrientosm/my_perl_5.36.0/bin:$PATH

 for f in /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/mafft/*.fa; do
     base=$(basename "$f" .fa);
cleaned="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/mafft/${base}_hmm.fasta";

if [ ! -f "$cleaned" ]; then
echo "=== Processing $f ===";
env  PERL5LIB=/home/barrientosm/my_perl_5.36.0/lib/perl5:/home/barrientosm/my_perl_5.36.0/lib/perl5/x86_64-linux \
     /home/barrientosm/my_perl_5.36.0/bin/perl -Ilib bin/HmmCleaner.pl "$f";
else
echo "=== Skipping $f (already cleaned) ===";
fi;
done > hmmcleaner_all.log 2>&1
```

**Renames orthologs**
```bash
for file in N0.*.fasta; do
    mv "$file" "${file#N0.}"
done
```

*Phylogeny Reconstruction for 4 or more species*


**Prottest**
```bash
cd /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest
tar -xzvf prottest-3.4.2-20160508.tar.gz
```

**Empirical model selection functions**
prottest_setup function: This function is designed to automate the process of identifying the best-fit model of amino acid replacement for a specified protein alignment using the third-party program ProtTest3 [Darriba et al., 2011]. The function is designed to test each amino acid substitution model in both the absence and presence of invariant sites, gamma categories, and a combination of the two.

Empirical_model_selection.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=prottest_setup
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --partition=agap_normal

# VESPA script path
VESPA_SCRIPT="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest/vespa.py"

# Process each fasta file
for fasta_file in *hmm.fasta; do
    [ ! -f "$fasta_file" ] && continue

    echo "Processing $fasta_file..."

    # Extract HOG ID (e.g. HOG0017898) from the filename
    hog_id=$(basename "$fasta_file" .fasta | sed 's/.*\(HOG[0-9]\+\).*/\1/')

    # Run vespa.py
    python2 "$VESPA_SCRIPT" prottest_setup -input="$fasta_file"

    # Rename ProtTest directory if exists
    prot_dir="ProtTest_Setup_${hog_id}_aln_hmm"
    [ -d "ProtTest_Setup_aln_hmm" ] && mv "ProtTest_Setup_aln_hmm" "$prot_dir"

    # Rename taskfarm file/dir if exists
    [ -e "setup_prottest_taskfarm" ] && mv "setup_prottest_taskfarm" "setup_prottest_taskfarm_${hog_id}"

done

echo "Processing complete."
```

**Preparing prottest**
Correction of paths and creation of directories to run prottest
prottest_preparation.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=prottest_prepare
#SBATCH --output=prepare_%j.log
#SBATCH --error=prepare_%j.err
#SBATCH --partition=agap_normal
#SBATCH --time=00:10:00
#SBATCH --mem=2G

# Configuration
PROTTEST_JAR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest/prottest-3.4.2/prottest-3.4.2.jar"
OUTPUT_BASE="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest/prottest_results"

# Process all taskfarm files
for tf_file in setup_prottest_taskfarm_*; do
    [ ! -f "$tf_file" ] && continue

    # Extract HOG ID
    hog_id=$(echo "$tf_file" | sed 's/.*\(HOG[0-9]\+\).*/\1/')
    
    # Define paths
    input_file="ProtTest_Setup_${hog_id}_aln_hmm"
    output_dir="${OUTPUT_BASE}/${hog_id}/models"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Write corrected command
    echo "java -jar $PROTTEST_JAR -i $input_file -o $output_dir -all-distributions" > "$tf_file"
    echo "Prepared $tf_file -> Output to $output_dir"
done

echo "Preparation completed. Ready to execute commands."
```

**Running prottest**
prottest_run.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=prottest_execute
#SBATCH --output=execute_%j.log
#SBATCH --error=execute_%j.err
#SBATCH --partition=agap_normal
#SBATCH --time=24:00:00
#SBATCH --mem=8G

module load bioinfo-cirad
module load java/jdk-24.0.2

# Execute all commands
for tf_file in setup_prottest_taskfarm_*; do
    [ ! -f "$tf_file" ] && continue
    
    echo "--------------------------------------------------"
    echo "Starting processing for $tf_file"
    echo "Command to execute:"
    cat "$tf_file"
    echo "--------------------------------------------------"
    
    # Execute the command and log output
    bash "$tf_file" >> "${tf_file}.log" 2>&1
    
    # Check exit status
    if [ $? -eq 0 ]; then
        echo "Successfully processed $tf_file"
    else
        echo "ERROR processing $tf_file - check ${tf_file}.log"
    fi
done

echo "All ProtTest executions completed"
```
**prottest_reader**
The prottest_reader function: 
This function automates the process of reading the output of ProtTest3. The function creates two output files: best_models.csv and best_supported_models.csv. The best models file reports the best-fit model of amino acid replacement (± rate-heterogeneity) reported by ProtTest3 whereas the best supported file reports the best-fit model of amino acid replacement (± rate-heterogeneity) supported by the third-party phylogenetic reconstruction program MrBayes [Ronquist and Huelsenbeck, 2003]. The two output files are given to enable the user to use different phylogenetic reconstruction software if desired.

*Phylogeny Reconstruction for 3 species*

map_alignments: The map_alignments function is designed to automate the conversion of protein MSAs to nucleotide. This process is mandatory as the codon substitution models of codeML require nucleotide alignments. Protein-MSA guided nucleotide MSAs are generated rather than directly generating nucleotide MSAs because: i) each column within the protein MSA represents aligned codons and therefore avoids aligning incomplete codons or frame-shift mutations, and ii) protein MSAs represent a comparison of the phenotype-producing elements of protein-coding sequences. The function begins by reading the protein MSA to map the non-gap position of each codon within the inferred nucleotide alignment. The sequence of the mapped codons is then inferred using the nucleotide dataset (preferably as a database) from earlier in the pipeline. If the mapping process results in no errors, the respective nucleotide MSA is created. All errors detected by the function will be returned within a separate log file. Please note that the map_alignments function requires the option -database to indicate the nucleotide dataset for correct sequence inference.

**Nucleotides database**
```bash
 sed -E 's/^>Hordeum_vulgare_([^ ]+) gene=.*/>H_vulgare|\1/' Cleaned_hordeum_full_CDS.fasta > Cleaned_hordeum_full_CDS_renamed.fasta

 sed -E 's/^>transcript:([^ ]+).*/>T_urartu|\1/' Cleaned_Triticum_urartu.IGDB.59.chr_cds.fasta > Cleaned_Triticum_urartu.IGDB.59.chr_cds_renamed.fasta

 sed -E 's/^>Tm\.(.+)/>T_monococcum|\1/' Cleaned_TA299-HC-cds_v1.1.fasta > Cleaned_TA299-HC-cds_v1.1_renamed.fasta
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
INPUT_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest/alignments_inputs"
VESPA_SCRIPT="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest/vespa.py"
CDS_DATABASE="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/cds_sequences/database/cds_database.fasta"
OUTPUT_BASE="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/gene_trees/Map_Gaps_cleaned"

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

**Number of alignment kept**
```bash
#!/bin/bash

BASE_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/gene_trees/Map_Gaps_cleaned"

total_hogs=0
with_files=0
with_3seq=0

echo "Analyzing HOG directories in $BASE_DIR..."

for hog_dir in "$BASE_DIR"/HOG*/; do
    ((total_hogs++))
    hog_id=$(basename "$hog_dir")

    #  Chemin correct vers le fichier généré par VESPA
    target_dir="${hog_dir}/Map_Gaps_${hog_id}"
    target_file="${target_dir}/${hog_id}"

    if [ -f "$target_file" ]; then
        ((with_files++))

        seq_count=$(grep -c '^>' "$target_file" 2>/dev/null || echo 0)
        if [ "$seq_count" -eq 3 ]; then
            ((with_3seq++))
        fi
    fi
done

echo "=== Analysis Results ==="
echo "Total HOG directories processed: $total_hogs"
echo "HOGs with alignment files: $with_files"
echo "HOGs with files containing exactly 3 sequences: $with_3seq"
echo "Detailed breakdown:"
echo "- HOGs with 0-2 sequences: $((with_files - with_3seq))"
echo "- HOGs with exactly 3 sequences: $with_3seq"
```

**Species Tree**
```bash
#!/bin/bash
#SBATCH --job-name=fasttree_3seq
#SBATCH --output=fasttree_%j.out
#SBATCH --error=fasttree_%j.err
#SBATCH --partition=agap_normal
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

module load bioinfo-cirad
module load fasttree/2.1.11

BASE_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/gene_trees/Map_Gaps_cleaned"
total_processed=0

echo "Processing 3-sequence alignments with FastTree..."

for hog_dir in "$BASE_DIR"/HOG*/; do
    hog_id=$(basename "$hog_dir")
    target_file="${hog_dir}/Map_Gaps_${hog_id}/${hog_id}"
    output_tree="${hog_dir}/Map_Gaps_${hog_id}/${hog_id}.tre"

    if [ -f "$target_file" ]; then
        seq_count=$(grep -c '^>' "$target_file" 2>/dev/null || echo 0)

        if [ "$seq_count" -eq 3 ]; then
            echo "Building tree for $hog_id..."
            
            # Run FastTree directly in the target directory
            (
                cd "${hog_dir}/Map_Gaps_${hog_id}" && \
                FastTree -wag "$hog_id" > "${hog_id}.tre" 2> "${hog_id}.fasttree.log"
            )
            
            # Verify output
            if [ -s "$output_tree" ]; then
                ((total_processed++))
                echo "Generated: ${hog_id}.tre"
            else
                echo "Failed: ${hog_id} (check Map_Gaps_${hog_id}/${hog_id}.fasttree.log)"
            fi
        fi
    fi
done

echo "=== FastTree Completion ==="
echo "Total 3-sequence alignments processed: $total_processed"
echo "Tree files (.tre) are in their respective Map_Gaps_HOG*/ directories"
```

########## infer gene tree
 pip2 install dendropy
# create the species tree
python2 vespa.py infer_genetree -input=nucleotides_MSA.fasta -species_tree=species_tree.nwk


