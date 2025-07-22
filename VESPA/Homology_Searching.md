----------------------------------------------------------------------------------
|||||                       Phase 2: Homology Searching                      |||||
----------------------------------------------------------------------------------
Details: This phase is concerned with identifying groups of similar sequences from
either BLAST or HMMER homology searches.

Commands: similarity_groups, reciprocal_groups, best_reciprocal_groups
----------------------------------------------------------------------------------

**best reciprocal blast:** 
Best-reciprocal similarity” requires that the sequences pass two criteria:
(i) they are sequences from different species, and 
(ii) in the pair-wise connection each sequence finds no other sequence in the respective species with a lower E-value. 
These requirements limit identification to orthologs (non-orthologs may be identified due to identical E-values or the absence of a true ortholog).

```bash
mkdir Similarity_Groups
```

```bash
#!/bin/bash
#SBATCH --job-name=best_reciprocal_groups
#SBATCH --output=log_%j_%x.out
#SBATCH --error=log_%j_%x.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=agap_normal

echo "=== Starting VESPA best_reciprocal_groups  ==="
echo "Start time: $(date)"
echo "Running on node: $(hostname)"

# Load conda and activate environment
source /nfs/work/agap_id-bin/img/snakemake/7.32.4/bin/activate base

# Check and create env only if not existing
if ! conda env list | grep -q "blast_env"; then
    conda create -n blast_env -y blast -c bioconda
fi

conda activate blast_env

# Run VESPA best_reciprocal_groups
python2 vespa.py best_reciprocal_groups \
  -input=BlastOutput_AllTriticeae.blastp.out \
  -format=blast \
  -database=translated_all_species.fasta \
  -e_value=0.001 \
  -alignment_length=75 \
  -percent_identity=75

echo "=== Finished ==="
echo "End time: $(date)"
```

**List of orthologs:** 
```bash
for f in *.fasta; do
  if [ "$(grep -c '^>' "$f")" -eq 3 ]; then
    echo "$f"
  fi
done > files_with_3_seqs.txt
```

**Copy orthologs:**
```bash
while read f; do
  cp "$f" /home/barrientosm/scratch/PALM/Similarity_Groups_Triticeae/Best_Reciprocal_Groups/orthologs/
done < files_with_3_seqs.txt
```
