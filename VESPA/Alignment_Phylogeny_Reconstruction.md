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

**Prottest**
```bash
cd /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest
tar -xzvf prottest-3.4.2-20160508.tar.gz
```

**Empirical model selection functions**
prottest_setup function: This function is designed to automate the process of identifying the best-fit model of amino acid replacement for a specified protein alignment using the third-party program ProtTest3 [Darriba et al., 2011]. The function is designed to test each amino acid substitution model in both the absence and presence of invariant sites, gamma categories, and a combination of the two.

```bash
#!/bin/bash
# Source directory with input fasta files
SOURCE_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/mafft"

# Target directory for ProtTest setup files
TARGET_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest/prottest_setup"

# VESPA script path
VESPA_SCRIPT="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/prottest/vespa.py"

# Create target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Process each hmm.fasta file
for fasta_file in "$SOURCE_DIR"/*hmm.fasta; do
    # Get just the filename without path
    filename=$(basename "$fasta_file")
    
    echo "Processing $filename..."
    
    # Run vespa.py command
    python2 "$VESPA_SCRIPT" prottest_setup -input="$fasta_file"
    
    # Copy the generated files to target directory
    cp -v "ProtTest_Setup_${filename}" "$TARGET_DIR/" 2>/dev/null
    cp -v "setup_prottest_taskfarm" "$TARGET_DIR/setup_prottest_taskfarm_${filename}" 2>/dev/null
    
    # Clean up working directory (optional)
    rm -f "ProtTest_Setup_${filename}" "setup_prottest_taskfarm"
done

echo "All files processed and copied to $TARGET_DIR"
```
