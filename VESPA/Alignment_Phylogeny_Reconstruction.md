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
#SBATCH --job-name=macse_align
#SBATCH --output=./macse_%j_out.txt
#SBATCH --error=./macse_%j_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=30G
#SBATCH --time=48:00:00
#SBATCH --partition=agap_normal

module load muscle/5.1  # ou la version disponible sur ton cluster

INPUT_DIR="/home/barrientosm/scratch/PALM/Similarity_Groups_Triticeae/Best_Reciprocal_Groups/orthologs"
OUTPUT_DIR="/home/barrientosm/scratch/PALM/Similarity_Groups_Triticeae/Best_Reciprocal_Groups/aligned"

mkdir -p "$OUTPUT_DIR"

for i in "$INPUT_DIR"/*.fasta
do
    base=$(basename "$i")
    muscle -in "$i" -out "$OUTPUT_DIR/${base%.fasta}.aln.fasta"
done
```
