----------------------------------------------------------------------------------
|||||        Phase 3: Alignment Assessment & Phylogeny Reconstruction        |||||
----------------------------------------------------------------------------------
Details: This phase combines multiple third-party programs to automate the
assessment, selection, and phylogenetic reconstruction of protein MSAs.

Commands: metal_compare, prottest_setup, prottest_reader, mrbayes_setup
----------------------------------------------------------------------------------

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

**Prottest**
```bash
cd /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest
tar -xzvf prottest-3.4.2-20160508.tar.gz
```

**Empirical model selection functions**
prottest_setup function: This function is designed to automate the process of identifying the best-fit model of amino acid replacement for a specified protein alignment using the third-party program ProtTest3 [Darriba et al., 2011]. The function is designed to test each amino acid substitution model in both the absence and presence of invariant sites, gamma categories, and a combination of the two.

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
