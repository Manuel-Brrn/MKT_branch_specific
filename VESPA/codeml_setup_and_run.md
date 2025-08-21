
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

**Parse Codeml output: dn_ds.py**
*Branch model and model m0*
```py
from Bio.Phylo.PAML import codeml
from Bio import Phylo
import io
import pandas as pd
import os
import re

# Obtenir le chemin actuel et extraire les informations
current_path = os.getcwd()

# Extraire le nom du HOG et de l'espèce depuis le chemin
hog_match = re.search(r'HOG(\d+)', current_path)
species_match = re.search(r'cleaned_(\w+)', current_path)

if hog_match and species_match:
    hog_number = hog_match.group(0)  # HOG0023328
    species_name = species_match.group(1)  # Turartu
    output_filename = f"{hog_number}_{species_name}_codeml.tab"
else:
    output_filename = "branch_results.csv"

# Charger le résultat codeml
results = codeml.read("out")
lnL = results["NSsites"][0]["lnL"]
branches = results["NSsites"][0]["parameters"]["branches"]
tree_newick = results["NSsites"][0]["tree"]

# Charger l'arbre Newick avec Biopython
tree = Phylo.read(io.StringIO(tree_newick), "newick")

# Associer chaque terminal (espèce) à sa longueur de branche
species_to_length = {term.name: term.branch_length for term in tree.get_terminals()}

# Fonction : trouver l'espèce correspondant à une branche
def match_species(branch, vals):
    # On compare la longueur de branche (t) dans "branches" et dans l'arbre
    for sp, length in species_to_length.items():
        if length is None:
            continue
        if abs(length - vals["t"]) < 1e-6:  # tolérance flottante
            return sp
    return f"branch_{branch}"

# Construire table
rows = []
for br, vals in branches.items():
    sp = match_species(br, vals)
    rows.append({
        "Species": sp,
        "Branch": br,
        "dN_branch": vals["dN"],
        "dS_branch": vals["dS"],
        "omega_branch": vals["omega"],
        "lnL": lnL
    })

df = pd.DataFrame(rows)
print(df.to_string(index=False))

# Sauvegarder avec le nom personnalisé
df.to_csv(output_filename, sep='\t', index=False)
print(f"\nResults saved to {output_filename}")
```

**Run: dn_ds.py in every directories branch and m0 model**
run_all_dn_ds.py:
```py
#!/usr/bin/env python3
import os
import glob
import subprocess

# ===== CONFIGURATION =====
BASE_DIR = "/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/translated_alignments_cleaned"
DN_DS_SCRIPT = "/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/translated_alignments_cleaned/dn_ds.py"
# =========================

# Check that dn_ds.py exists
if not os.path.exists(DN_DS_SCRIPT):
    print(f"ERROR: dn_ds.py not found at: {DN_DS_SCRIPT}")
    exit(1)

# Find all Omega0_5 directories - CORRECTION ICI
patterns = [
    os.path.join(BASE_DIR, "Inferred_Genetree_HOG*/HOG*/Codeml_Setup_HOG_codeml_input/cleaned_*/model_branch/Omega0_5"),
    os.path.join(BASE_DIR, "Inferred_Genetree_HOG*/HOG*/Codeml_Setup_HOG_codeml_input/cleaned/m0/Omega0_5")
]

directories = []
for pattern in patterns:
    directories.extend(glob.glob(pattern))

# Remove duplicates
directories = list(set(directories))

if not directories:
    print("No Omega0_5 directories found.")
    exit(0)

print(f"Found {len(directories)} Omega0_5 directories to process.")

# Run dn_ds.py inside each directory
for i, directory in enumerate(directories, 1):
    print(f"\n🔹 [{i}/{len(directories)}] Processing: {directory}")
    try:
        result = subprocess.run(
            ["python3", DN_DS_SCRIPT],
            cwd=directory,   # execute inside this Omega0_5 folder
            capture_output=True,
            text=True,
            timeout=600  # 10 min timeout per dir
        )

        if result.stdout:
            print("Output:\n", result.stdout.strip())
        if result.stderr:
            print("Errors:\n", result.stderr.strip())

    except subprocess.TimeoutExpired:
        print("Timeout: Script took too long")
    except Exception as e:
        print(f"Failed: {e}")

print("\nAll directories processed!")
```

**Create a summary table from codeml parsers**
Create two files one file with dn,ds, and omega for each branch and omega for m0 for every genes, and one file with every lnL of each branch and the lnL of the m0 model

summary_codeml.py
```py
summary_dn_ds_analysis.py:
#!/usr/bin/env python3
import os
import glob
import pandas as pd
import re

# ===== CONFIGURATION =====
BASE_DIR = "/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/translated_alignments_cleaned"
DNDS_OUTPUT = os.path.join(BASE_DIR, "summary_dn_ds_wide.tsv")
LNL_OUTPUT = os.path.join(BASE_DIR, "summary_lnL_wide.tsv")
# =========================

# Fichiers branch model (*.tab)
branch_files = glob.glob(
    os.path.join(BASE_DIR, "Inferred_Genetree_HOG*/HOG*/Codeml_Setup_HOG_codeml_input/cleaned_*/model_branch/Omega0_5/*_codeml.tab")
)

# Fichiers m0 model (branch_results.csv)
m0_files = glob.glob(
    os.path.join(BASE_DIR, "Inferred_Genetree_HOG*/HOG*/Codeml_Setup_HOG_codeml_input/cleaned/m0/Omega0_5/branch_results.csv")
)

print(f"Found {len(branch_files)} branch model files and {len(m0_files)} m0 model files.")

# Résultats
branch_results = {}  # {HOG: {species: {"dN":, "dS":, "omega":}}}
branch_lnl = {}      # {HOG: {species: lnL}}
m0_results = {}      # {HOG: {"omega_m0": val, "lnL_m0": val}}

# Charger les résultats du modèle branch
for f in branch_files:
    filename = os.path.basename(f)  # ex: HOG0023328_Turartu_codeml.tab
    match = re.match(r"(HOG\d+)_(\w+)_codeml\.tab", filename)
    if not match:
        continue
    hog, species = match.groups()

    df = pd.read_csv(f, sep="\t")
    # On prend la dernière ligne pour la branche testée
    row = df.iloc[-1]

    if hog not in branch_results:
        branch_results[hog] = {}
    if hog not in branch_lnl:
        branch_lnl[hog] = {}

    branch_results[hog][species] = {
        "dN": row["dN_branch"],
        "dS": row["dS_branch"],
        "omega": row["omega_branch"]
    }
    branch_lnl[hog][species] = row["lnL"]

# Charger les résultats du modèle m0
for f in m0_files:
    match = re.search(r"(HOG\d+)", f)
    if not match:
        continue
    hog = match.group(1)

    df = pd.read_csv(f, sep="\t")
    row = df.iloc[0]  # m0 a une seule valeur
    m0_results[hog] = {"omega_m0": row["omega_branch"], "lnL_m0": row["lnL"]}

# Construire les tableaux finaux
all_hogs = sorted(set(branch_results.keys()) | set(m0_results.keys()))
dn_ds_rows = []
lnl_rows = []

# Liste de toutes les espèces pour ordonner les colonnes
all_species = sorted(set(
    species for hog_data in branch_results.values()
    for species in hog_data.keys()
))

for hog in all_hogs:
    dn_ds_row = {"HOG": hog}
    lnl_row = {"HOG": hog}

    # Ajouter les données pour chaque espèce (dans l'ordre)
    for species in all_species:
        if hog in branch_results and species in branch_results[hog]:
            dn_ds_row[f"{species}_dN"] = branch_results[hog][species]["dN"]
            dn_ds_row[f"{species}_dS"] = branch_results[hog][species]["dS"]
            dn_ds_row[f"{species}_omega"] = branch_results[hog][species]["omega"]

        if hog in branch_lnl and species in branch_lnl[hog]:
            lnl_row[f"{species}_lnL"] = branch_lnl[hog][species]

    # Ajouter le modèle m0
    if hog in m0_results:
        dn_ds_row["omega_m0"] = m0_results[hog]["omega_m0"]
        lnl_row["lnL_m0"] = m0_results[hog]["lnL_m0"]

    dn_ds_rows.append(dn_ds_row)
    lnl_rows.append(lnl_row)

# Créer les DataFrames
dn_ds_df = pd.DataFrame(dn_ds_rows)
lnl_df = pd.DataFrame(lnl_rows)

# Ordonner les colonnes pour dn_ds
dn_ds_cols = ["HOG"]
for species in all_species:
    dn_ds_cols.extend([f"{species}_dN", f"{species}_dS", f"{species}_omega"])
dn_ds_cols.append("omega_m0")

# Ordonner les colonnes pour lnl
lnl_cols = ["HOG"]
for species in all_species:
    lnl_cols.append(f"{species}_lnL")
lnl_cols.append("lnL_m0")

# Garder seulement les colonnes qui existent
dn_ds_df = dn_ds_df[[col for col in dn_ds_cols if col in dn_ds_df.columns]]
lnl_df = lnl_df[[col for col in lnl_cols if col in lnl_df.columns]]

# Sauvegarder
dn_ds_df.to_csv(DNDS_OUTPUT, sep="\t", index=False)
lnl_df.to_csv(LNL_OUTPUT, sep="\t", index=False)

print(f"Saved DN/DS file: {DNDS_OUTPUT}")
print(f"Saved lnL file: {LNL_OUTPUT}")
print(f"Species order: {all_species}")
```
