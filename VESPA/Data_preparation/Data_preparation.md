
----------------------------------------------------------------------------------
|||||                        **Phase 1: Data Preparation**                       |||||
----------------------------------------------------------------------------------
This phase included for users new to bioinformatics. The phase prepares downloaded
genomes for homology searching.

Commands: clean, ensemble_clean, translate, create_database, gene_selection
----------------------------------------------------------------------------------

### **Clean cds**
Each sequence is confirmed as protein coding by using a conditional statement to verify that the nucleotide sequence contains only complete codons (i.e. the length of the sequence is exactly divisible by 3). And transcript with stop codon are eliminated.
```bash
python2 vespa.py clean –input=user_input
```

### **Traduction cds in protein**
The translate function translates nucleotide sequences that passed the QC filter of either clean function into amino acid sequences in the first reading frame forward only. The function operates by splitting the nucleotide sequence into codons and then translating them into their respective amino acids. The cleave_terminal option is enabled by default and is designed to cleave the terminal stop codon of each sequence 
```bash
python2 vespa.py translate -input=Cleaned_TA299-HC-cds_v1.1.renamed.fasta -cleave_terminal=True
```

### **Rename sequences names**
#### **Urartu**
```bash
 sed -E 's/^>transcript:([^ ]+).*/>T_urartu|\1/' Translated_Cleaned_Triticum_urartu.IGDB.59.chr_cds.fasta > translated_urartu_renamed.fasta
```
#### **Monoccocum**
```bash
 sed -E 's/^>Tm\.(.+)/>T_monococcum|\1/' Translated_Cleaned_TA299-HC-cds_v1.1.fasta > translated_monococcum_renamed.fasta
```
#### **Speltoides**
```bash
sed -E 's/^>([^ ]+).*/>Ae_speltoides|\1/' Cleaned_speltoides_CDS.fasta > Cleaned_speltoides_CDS_renamed.fasta
```
#### **Mutica**
```bash
sed -E 's/^>([^ ]+).*/>Ae_mutica|\1/' Cleaned_Ammut_EIv1.0.release_cds_clean.fasta > Cleaned_Ammut_EIv1.0.release_cds_clean_renamed.fasta
```
#### **Hordeum**
```bash
 sed -E 's/^>Hordeum_vulgare_([^ ]+) gene=.*/>H_vulgare|\1/' Translated_Cleaned_hordeum_full_CDS.fasta > translated_hordeum_renamed.fasta
```

**Cd-hit**
CD-HIT (Cluster Database at High Identity with Tolerance) is a widely-used and fast clustering program designed for comparing and clustering large sets of protein or nucleotide sequences. It is especially popular in bioinformatics pipelines for removing redundancy from datasets, such as after assembling transcriptomes or downloading large protein datasets from databases like UniProt.

Reduce Redundancy: Collapse highly similar sequences into representative clusters.
Speed Up Downstream Analyses: Reduces dataset size for alignment, annotation, etc.
Dereplication: Ensures only unique (or representative) sequences remain.
Preprocessing step in metagenomics, proteomics, transcriptomics, etc.

 Advantages of CD-HIT
    Speed:
        CD-HIT is extremely fast because of its short word filter and greedy incremental clustering method.
        Optimized for large-scale sequence datasets (e.g., millions of sequences).

    Memory Efficiency:
        Designed to handle large datasets even on machines with modest memory.

    Customizable Identity Thresholds:
        Users can specify sequence identity thresholds (e.g., -c 0.9 for 90% identity) to control clustering stringency.

    Supports Protein and Nucleotide Sequences:
        cd-hit for proteins.
        cd-hit-est for nucleotide sequences.

    Multiple Output Options:

        Clustered representative sequences.
        Cluster membership files for downstream analysis.

As mutica contains a lot of transcript I used cd-hit to remove any higly similar sequences

 cd_hit.sbatch
```bash
#!/bin/bash
#SBATCH --job-name=cd_hit
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --partition=agap_normal

module load bioinfo-cirad
module load cd-hit/4.8.1

cd-hit -i Cleaned_Ammut_EIv1.0.release_cds.fasta -o Cleaned_Ammut_EIv1.0.release_cds_clean.fasta -c 1.0
```

 ### **create_database for Blast**
The create_database function was designed for users to concatenate multiple genomes into the single database required for homology searching. The function operates by building the database a single sequence at a time.
Converts a FASTA file into a searchable BLAST database, allowing to:

    Run BLAST queries against your dataset (instead of NCBI’s databases).

    Accelerate sequence searches using indexed files.

```bash
cat *fasta* > translated_all_species.fasta
```

```bash
# Setup environment
source /nfs/work/agap_id-bin/img/snakemake/7.32.4/bin/activate base

# Create env if missing
if ! conda env list | grep -q "blast_env"; then
    conda create -n blast_env -y blast -c bioconda
fi

conda activate blast_env

makeblastdb \
  -in translated_all_species.fasta \
  -dbtype prot \
  -parse_seqids \
  -title "MultiSpeciesDB" \
  -out multispecies_prot_db \
  -blastdb_version 4
```

Parameter	Explanation:
- in translated_all_species.fasta	Input FASTA file containing protein sequences.
- dbtype prot	Specifies this is a protein database (use nucl for nucleotides).
- parse_seqids	Preserves original sequence IDs in the database (critical for tracking hits).
- title "MultiSpeciesDB"	Human-readable name for the database.
- out multispecies_prot_db	Base name for output files (will generate .phr, .pin, .psq, etc.).
- blastdb_version 4	Uses the newer BLAST database format (supports larger files).

 ### **Blast against the base**

 blast_PALM.sbatch
 ```bash
#!/bin/bash
#SBATCH --job-name=reciprocal_blast
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4  # ADDED
#SBATCH --mem=9G
#SBATCH --time=48:00:00
#SBATCH --partition=agap_normal
#SBATCH --array=1-3

echo "Starting reciprocal BLAST at $(date)"

### Setup environment
source /nfs/work/agap_id-bin/img/snakemake/7.32.4/bin/activate base

# Create env if missing
if ! conda env list | grep -q "blast_env"; then
    conda create -n blast_env -y blast -c bioconda
fi

conda activate blast_env

# Query directory
QUERY_DIR="/lustre/barrientosm/PALM/QUERIES"
QUERY_FILES=(
    "translated_hordeum.fasta"
    "translated_monoccocum.fasta"
    "translated_urartu.fasta"
)

# BLAST database PREFIX (without extensions)
BLAST_DB_PREFIX="/lustre/barrientosm/PALM/Blast_triticea/multispecies_prot_db"

# Get query file for this task
QUERY_FILE="${QUERY_DIR}/${QUERY_FILES[$SLURM_ARRAY_TASK_ID-1]}"

# Run BLAST with critical parameters
blastp -query "$QUERY_FILE" \
       -db "$BLAST_DB_PREFIX" \
       -out "${QUERY_FILE}.blastp.out" \
       -outfmt 6 \
       -evalue 1e-7 \
       -num_threads $SLURM_CPUS_PER_TASK \
       -max_target_seqs 5000 \
       -seg yes \
       -soft_masking true \
       2> "${QUERY_FILE}.blastp.err" \
       1> "${QUERY_FILE}.blastp.log"

echo "BLAST finished for $QUERY_FILE at $(date)"
```

### **Input for homology section 2**
```bash
cat *blastp.out* > BlastOutput_AllTriticeae.blastp.out
```
