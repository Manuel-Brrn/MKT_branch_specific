
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
        echo "❌ No FASTA found in $hog_dir"
    fi

    # 4. Copy and clean tree
    tree_file=$(find "$hog_dir" -name "HOG*.tre" | head -n 1)
    if [ -f "$tree_file" ]; then
        cp "$tree_file" "$codeml_dir/cleaned.tre"
        sed -i 's/|[^),]*//g; s/_//g' "$codeml_dir/cleaned.tre"
    else
        echo "❌ No TREE found in $hog_dir"
    fi

    # Verification
    echo "Created files:"
    echo "- Branch table: $hog_dir/branch_table.txt"
    echo "- CodeML inputs: $codeml_dir/{cleaned.fasta,cleaned.tre}"
done
```



Script to run codeml_setup function in every directory
```bash
#!/bin/bash

# Set PERL5LIB environment
export PERL5LIB=/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/03_scripts/dn_ds_pipeline/VESPA/VESPA-1.0.1:$PERL5LIB

# Base path to vespa.py
VESPA_SCRIPT="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/alignment/cleaned_alignments/vespa.py"

# Loop through each target HOG directory
for hog_dir in Inferred_Genetree_HOG*/HOG*; do
    echo "Entering directory: $hog_dir"

    cd "$hog_dir" || { echo "❌ Failed to enter $hog_dir"; continue; }

    # Check inputs exist
    if [[ ! -f branch_table.txt ]]; then
        echo "❌ Missing branch_table.txt in $hog_dir"
        cd - > /dev/null
        continue
    fi

    if [[ ! -d HOG_codeml_input ]]; then

        echo "❌ Missing HOG_codeml_input directory in $hog_dir"
        cd - > /dev/null
        continue
    fi

    # Run VESPA codeml_setup
    echo "Running VESPA codeml_setup..."
    python2 "$VESPA_SCRIPT" codeml_setup -input=HOG_codeml_input -branch_file=branch_table.txt

    cd - > /dev/null
done
```

**Modify the codeml task file**
```bash
for dir in Inferred_Genetree_HOG*/HOG*/; do
  if [ -f "$dir/codeml_taskfarm.txt" ]; then
    # Crée un nouveau fichier avec extension .sh
    outfile="$dir/codeml_taskfarm_fullpath.sh"
    
    sed "s|cd Codeml_Setup|cd $(realpath "$dir")/Codeml_Setup|g" "$dir/codeml_taskfarm.txt" \
      > "$outfile"
    
    # Rendre le fichier exécutable
    chmod +x "$outfile"
    
    echo "Created and made executable: $outfile"
  fi
done
```
