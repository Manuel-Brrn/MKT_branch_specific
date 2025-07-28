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
