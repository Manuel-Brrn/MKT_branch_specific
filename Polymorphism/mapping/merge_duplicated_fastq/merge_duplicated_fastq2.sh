#!/bin/bash

# Arguments check
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
mkdir -p "$OUTPUT_DIR"

echo "Searching for split FASTQ files in: $INPUT_DIR"

# --- Merge R1 files ---
echo -e "\nProcessing R1 files..."
find "$INPUT_DIR" -type f -name "*_R1_*.fastq.bz2" \
    | sed -E 's|.*/([^/]+)_R1_[0-9]+\.fastq\.bz2|\1|' \
    | sort -u \
    | while read -r base; do
        files=($(find "$INPUT_DIR" -type f -name "${base}_R1_*.fastq.bz2" | sort -V))
        echo "Merging R1 for $base"
        printf '  %s\n' "${files[@]##*/}"
        bzcat "${files[@]}" | bzip2 > "$OUTPUT_DIR/${base}_R1.fastq.bz2"
done

# --- Merge R2 files ---
echo -e "\nProcessing R2 files..."
find "$INPUT_DIR" -type f -name "*_R2_*.fastq.bz2" \
    | sed -E 's|.*/([^/]+)_R2_[0-9]+\.fastq\.bz2|\1|' \
    | sort -u \
    | while read -r base; do
        files=($(find "$INPUT_DIR" -type f -name "${base}_R2_*.fastq.bz2" | sort -V))
        echo "Merging R2 for $base"
        printf '  %s\n' "${files[@]##*/}"
        bzcat "${files[@]}" | bzip2 > "$OUTPUT_DIR/${base}_R2.fastq.bz2"
done

# --- Copy non-split files ---
echo -e "\nCopying non-split files..."
find "$INPUT_DIR" -type f -name "*.fastq.bz2" | grep -vE '_R[12]_[0-9]+\.fastq\.bz2' | while read -r file; do
    echo "Copying: $(basename "$file")"
    cp "$file" "$OUTPUT_DIR/"
done

echo -e "\nAll done. Merged files are in: $OUTPUT_DIR"
