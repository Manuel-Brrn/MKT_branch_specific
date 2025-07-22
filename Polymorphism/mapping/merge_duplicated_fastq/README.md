# FASTQ Merge Script

## Description
This script merges split FASTQ files (with `_00`/`_01` suffixes) while preserving the original naming convention. It handles both R1 and R2 files separately and copies non-split files unchanged.

## Usage
```bash
./merge_duplicated_fastq.sh <input_directory> <output_directory>
```

## Example
```bash
# Merge files in Urartu_FatsQ directory:
./merge_duplicated_fastq.sh \
  /home/user/Urartu_FatsQ \
  /home/user/merged_fastqs
```

## File Naming Requirements
- Split files must follow pattern: `TrXXX_XXXXX_XXX_00.R1.fastq.bz2` and `TrXXX_XXXXX_XXX_01.R1.fastq.bz2`
- Non-split files: `TrXXX_XXXXX_XXX.R1.fastq.bz2`

## What It Does
1. Finds all split file pairs (`_00` + `_01`) for both R1 and R2
2. Merges them using `bzcat` and `bzip2`
3. Copies non-split files as-is
4. Preserves original naming format (`TrXXX_XXXXX_XXX.R1.fastq.bz2`)

## Output Structure
```
merged_output/
Tr309_TTAGGC_L005.R1.fastq.bz2   # Merged from _00 + _01
Tr309_TTAGGC_L005.R2.fastq.bz2
Tr206_GATCAG_L002.R1.fastq.bz2   # Copied as-is (non-split)
```

## Requirements
- Linux/Unix environment
- `bzip2` installed
- Files must be named consistently as shown above

## Troubleshooting
If no files are merged:
1. Verify your filenames match the expected patterns
2. Check file permissions in both input/output directories
3. Ensure files have both `_00` and `_01` pairs

> **Note**: The script will skip files that don't have matching pairs.
