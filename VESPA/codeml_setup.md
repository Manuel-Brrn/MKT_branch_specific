
**codeml setup**
Create a branch file, and renames the sequences to makes work the codeml set_up command, in order to run the branch model

```bash
#!/bin/bash

# Loop over all target directories
for hog_dir in Inferred_Genetree_HOG*/HOG*; do
    echo "Processing $hog_dir"
    
    # 1. Create branch_table.txt in HOG directory (not in codeml input dir)
    cat > "$hog_dir/branch_table.txt" <<EOF
Aespeltoides
Aemutica
Turartu
Tmonococcum
Aegilops_ancestor: Aespeltoides, Aemutica
Triticum_ancestor: Turartu, Tmonococcum
EOF

    # 2. Create codeml input directory with HOG* prefix
    codeml_dir="${hog_dir}/HOG_codeml_input"  # Changed to your requested name
    mkdir -p "$codeml_dir"

    # 3. Copy and clean fasta
    fasta_file=$(find "$hog_dir" -name "HOG*.fasta" | head -n 1)
    if [ -f "$fasta_file" ]; then
        cp "$fasta_file" "$codeml_dir/cleaned.fasta"
        sed -i 's/|.*//; s/_//g' "$codeml_dir/cleaned.fasta"
    else
        echo "No FASTA found in $hog_dir"
    fi

    # 4. Copy and clean tree
    tree_file=$(find "$hog_dir" -name "HOG*.tre" | head -n 1)
    if [ -f "$tree_file" ]; then
        cp "$tree_file" "$codeml_dir/cleaned.tre"
        sed -i 's/|[^),]*//g; s/_//g' "$codeml_dir/cleaned.tre"
    else
        echo "No TREE found in $hog_dir"
    fi

    # Verification
    echo "Created files:"
    echo "- Branch table: $hog_dir/branch_table.txt"
    echo "- CodeML inputs: $codeml_dir/{cleaned.fasta,cleaned.tre}"
done
```



**Script to run codeml_setup function in every directory**
```bash
#!/bin/bash

# Set PERL5LIB environment
export PERL5LIB=/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/03_scripts/dn_ds_pipeline/VESPA/VESPA-1.0.1:$PERL5LIB

# Base path to vespa.py
VESPA_SCRIPT="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/cleaned_alignments/vespa.py"

# Loop through each target HOG directory
for hog_dir in Inferred_Genetree_HOG*/HOG*; do
    echo "Entering directory: $hog_dir"

    cd "$hog_dir" || { echo "Failed to enter $hog_dir"; continue; }

    # Check inputs exist
    if [[ ! -f branch_table.txt ]]; then
        echo "Missing branch_table.txt in $hog_dir"
        cd - > /dev/null
        continue
    fi

    if [[ ! -d HOG_codeml_input ]]; then

        echo "Missing HOG_codeml_input directory in $hog_dir"
        cd - > /dev/null
        continue
    fi

    # Run VESPA codeml_setup
    echo "Running VESPA codeml_setup..."
    python2 "$VESPA_SCRIPT" codeml_setup -input=HOG_codeml_input -branch_file=branch_table.txt

    cd - > /dev/null
done
```

**Correct directories, codeml control file and create task_file to run codeml on every directories**
control files are based on : Beginner's Guide on the Use of PAML to Detect Positive Selection Sandra Alvarez-Carretero 2023
*Branch model*

codeml_set_up_branch.sh:
```bash
#!/bin/bash

ROOT="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/translated_alignments_cleaned"

TASKFARM="$ROOT/codeml_taskfarm_fullpath.sh"
> "$TASKFARM"

find "$ROOT" -type d -path "*/Inferred_Genetree_HOG*/HOG*/Codeml_Setup_HOG_codeml_input/cleaned_*" | while read codeml_dir; do
    echo "[INFO] Processing $codeml_dir"

    # si modelA existe → renommer
    if [ -d "$codeml_dir/modelA" ] && [ ! -d "$codeml_dir/model_branch" ]; then
        mv "$codeml_dir/modelA" "$codeml_dir/model_branch"
    fi

    if [ -d "$codeml_dir/model_branch" ]; then
        if [ -d "$codeml_dir/model_branch/Omega0" ]; then
            mv "$codeml_dir/model_branch/Omega0" "$codeml_dir/model_branch/Omega0_5"
        fi

        if [ -d "$codeml_dir/model_branch/Omega0_5" ]; then
            ctl="$codeml_dir/model_branch/Omega0_5/codeml.ctl"
            cat > "$ctl" <<EOF
seqfile = align.phy
treefile = tree
outfile = out

noisy = 3
verbose = 1
runmode = 0

seqtype = 1
ndata = 1
icode = 0
cleandata = 1

model = 2
NSsites = 0
CodonFreq = 7
estFreq = 0
clock = 0
fix_omega = 0
omega = 0.5
EOF

            echo "cd $codeml_dir/model_branch/Omega0_5; codeml" >> "$TASKFARM"
        fi
    fi

    #  Remove other directories except modelA, model_branch, Omega0, Omega0_5
    find "$codeml_dir" -mindepth 1 -maxdepth 1 -type d \
        ! -name "modelA" \
        ! -name "model_branch" \
        ! -name "Omega0" \
        ! -name "Omega0_5" \
        -exec rm -rf {} +
done
```

*Site model m0, m1 nearly neutral, m2 Selection, m7, m8*
codeml_set_up_site_model.sh:
```bash
codeml_set_up_site_model.sh
#!/bin/bash

ROOT="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/translated_alignments_cleaned"

TASKFARM="$ROOT/codeml_taskfarm_fullpath_site_model.sh"
> "$TASKFARM"

# ----------------------------
# Modèles à garder
# ----------------------------
MODELS=("m0" "m1Neutral" "m2Selection" "m7" "m8")

# Mapping modèle -> valeur codeml
declare -A MODEL_MAP
MODEL_MAP=( ["m0"]=0 ["m1Neutral"]=1 ["m2Selection"]=2 ["m7"]=7 ["m8"]=8 )

# ----------------------------
# Parcours des répertoires Codeml_Setup
# ----------------------------
find "$ROOT" -type d -name "cleaned" | while read codeml_dir; do
    echo "[INFO] Processing $codeml_dir"

    # Supprimer modèles non désirés
    for d in "$codeml_dir"/*/; do
        d=$(basename "$d")
        if [[ ! " ${MODELS[@]} " =~ " $d " ]]; then
            rm -rf "$codeml_dir/$d"
        fi
    done

    # Parcourir les modèles gardés
    for model in "${MODELS[@]}"; do
        if [ -d "$codeml_dir/$model" ]; then
            echo "  Keeping model: $model"

            # Supprimer Omega* sauf Omega0
            for omega in "$codeml_dir/$model"/Omega*; do
                if [ -d "$omega" ] && [ "$(basename "$omega")" != "Omega0" ]; then
                    rm -rf "$omega"
                fi
            done

            # Renommer Omega0 -> Omega0_5
            if [ -d "$codeml_dir/$model/Omega0" ]; then
                mv "$codeml_dir/$model/Omega0" "$codeml_dir/$model/Omega0_5"
            fi

            # Créer codeml.ctl
            if [ -d "$codeml_dir/$model/Omega0_5" ]; then
                ctl="$codeml_dir/$model/Omega0_5/codeml.ctl"
                cat > "$ctl" <<EOF
seqfile = align.phy
treefile = tree
outfile = out

noisy = 3
verbose = 1
runmode = 0

seqtype = 1
ndata = 1
icode = 0
cleandata = 1

model = ${MODEL_MAP[$model]}
NSsites = 0
CodonFreq = 7
estFreq = 0
clock = 0
fix_omega = 0
omega = 0.5
EOF

                # Ajouter au taskfarm
                echo "cd $codeml_dir/$model/Omega0_5; codeml" >> "$TASKFARM"
            fi
        fi
    done
done
```


Tasks files: 
- codeml_taskfarm_fullpath.sh (branch model)
- codeml_taskfarm_fullpath_site_model.sh (site model)

**Run Codeml for each alignments**
run_codeml.sbatch:
```bash
#!/bin/bash
#SBATCH --job-name=codeml
#SBATCH --output=codeml_%A_%a.out
#SBATCH --error=codeml_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --array=1-12342%10 #### adapt the number of jobs
#SBATCH --time=48:00:00
#SBATCH --partition=agap_long

module load bioinfo-cirad
module load paml/4.9.0

# Task file with one codeml command per line
TASKFILE="codeml_taskfarm_fullpath.sh"

# Get command for this array task
CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKFILE")

if [ -n "$CMD" ]; then
    echo "Running: $CMD"
    eval "$CMD"
else
    echo "No command found for task ID $SLURM_ARRAY_TASK_ID"
fi
```


