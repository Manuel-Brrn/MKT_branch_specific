#!/bin/bash

# Configuration - Set these absolute paths
INDIVIDUALS_INFO="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/01_raw_data/Accessions_table/individuals_info.txt"
FASTQ_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/01_raw_data/Trans_raw_fastq"

# Validate paths
if [ ! -f "$INDIVIDUALS_INFO" ]; then
    echo "Error: Individuals info file not found at $INDIVIDUALS_INFO" >&2
    exit 1
fi

if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: FASTQ directory not found at $FASTQ_DIR" >&2
    exit 1
fi

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 \"<species_name>\" <target_directory>" >&2
    echo "Note: Species name must match exactly as shown below" >&2
    echo -e "\nAvailable species names:"
    awk 'NR>1 {print $2}' "$INDIVIDUALS_INFO" | sort -u | while read -r species; do
        echo "  \"$species\""
    done
    exit 1
fi

species_name="$1"
target_dir="$2"

# Create target directory
mkdir -p "$target_dir" || { echo "Failed to create target directory"; exit 1; }

# Extract arcad_codes for exact species name match
echo -e "\nSearching for species: '$species_name'"
arcad_codes=$(awk -F'\t' -v species="$species_name" '
NR>1 && $2 == species {print $1}
' "$INDIVIDUALS_INFO" | sort -u)

if [ -z "$arcad_codes" ]; then
    echo "Error: No exact match found for species: '$species_name'" >&2
    echo -e "\nAvailable species names (must match exactly):"
    awk 'NR>1 {print $2}' "$INDIVIDUALS_INFO" | sort -u | while read -r species; do
        echo "  \"$species\""
    done
    exit 1
fi

echo -e "Found accessions:\n$arcad_codes"

# Copy FASTQ files
total_copied=0
for code in $arcad_codes; do
    echo -e "\nProcessing $code..."
    while IFS= read -r -d $'\0' file; do
        echo "Copying $(basename "$file")"
        cp -v "$file" "$target_dir/" && ((total_copied++))
    done < <(find "$FASTQ_DIR" -maxdepth 1 -type f -iname "${code}_*.fastq*" -print0)
done

echo -e "\nDone! Copied $total_copied files to $target_dir"
ls -lh "$target_dir"
