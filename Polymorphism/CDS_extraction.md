*CDS extraction*

**UTR Length Extractor**

A script to extract 5' and 3' UTR lengths from GFF3 files, with transcript IDs and strand information.

Features:
- Extracts UTR lengths for all mRNA features
- Handles both 5' and 3' UTRs
- Preserves strand information
- Identifies transcripts missing UTR annotations
- Provides validation statistics

**Usage**
```bash
./get_utr_lengths.sh <input.gff3> <output.tab>
```

utr_lenght.sh
```bash
#!/bin/bash

# Usage: ./get_utr_lengths.sh <input_gff> <output_file>
# Example: ./get_utr_lengths.sh /path/to/Triticum_urartu.IGDB.59.chr.gff3 /path/to/utr_lengths.tab

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <input_gff> <output_file>"
    exit 1
fi

INPUT_GFF="$1"
OUTPUT_FILE="$2"

echo "========================================"
echo "Extracting UTR lengths from: $INPUT_GFF"
echo "Output will be saved to: $OUTPUT_FILE"
echo "========================================"

# Step 1: Extract UTR lengths
echo "Processing GFF file to calculate UTR lengths..."
awk '
BEGIN {
    FS = "\t";
    OFS = "\t";
    print "contigID", "five_prime_UTR_len", "three_prime_UTR_len", "strand";
}
$3 == "mRNA" {
    split($9, attr, ";");
    for (i in attr) {
        if (attr[i] ~ /^ID=transcript:/) {
            split(attr[i], id, ":");
            transcript = id[2];
            strand = $7;
            utr5[transcript] = 0;
            utr3[transcript] = 0;
            transcript_strand[transcript] = strand;
        }
    }
}
$3 == "five_prime_UTR" {
    split($9, attr, ";");
    for (i in attr) {
        if (attr[i] ~ /^Parent=transcript:/) {
            split(attr[i], parent, ":");
            transcript = parent[2];
            len = $5 - $4 + 1;
            utr5[transcript] += len;
        }
    }
}
$3 == "three_prime_UTR" {
    split($9, attr, ";");
    for (i in attr) {
        if (attr[i] ~ /^Parent=transcript:/) {
            split(attr[i], parent, ":");
            transcript = parent[2];
            len = $5 - $4 + 1;
            utr3[transcript] += len;
        }
    }
}
END {
    for (t in utr5) {
        print t, utr5[t] + 0, utr3[t] + 0, transcript_strand[t];
    }
}' "$INPUT_GFF" > "$OUTPUT_FILE"

echo "UTR extraction complete. Validating output..."

# Step 2: Verify counts
TOTAL_MRNA=$(grep -P "\tmRNA\t" "$INPUT_GFF" | wc -l)
OUTPUT_MRNA=$(tail -n +2 "$OUTPUT_FILE" | wc -l)

echo "Total mRNAs in GFF: $TOTAL_MRNA"
echo "mRNAs in output: $OUTPUT_MRNA"

if [ "$TOTAL_MRNA" -ne "$OUTPUT_MRNA" ]; then
    echo "Warning: Output is missing $((TOTAL_MRNA - OUTPUT_MRNA)) mRNAs."

    # Identify missing contigs
    echo "Identifying missing contigs..."
    grep -P "\tmRNA\t" "$INPUT_GFF" | \
    awk -F'\t' '{split($9, a, /[:;]/); print a[2]}' | \
    sort -u > all_mrna_ids.txt

    tail -n +2 "$OUTPUT_FILE" | cut -f1 | sort -u > output_contigs.txt

    comm -23 all_mrna_ids.txt output_contigs.txt > missing_contigs.txt

    MISSING_COUNT=$(wc -l < missing_contigs.txt)
    echo "Missing contigs written to: missing_contigs.txt ($MISSING_COUNT entries)"
else
    echo "Success: All mRNAs accounted for in output."
fi

echo "========================================"
echo "Script completed successfully"
echo "========================================"
```

**Fastas stats**
Pipeline for analyzing FASTA/FASTA.GZ files with SeqKit statistics, N-count analysis, and sequence composition metrics.

Features:
- **Multi-file processing**: Handles multiple input files in one run
- **Advanced statistics**: Uses SeqKit for comprehensive sequence analysis
- **N/X counting**: Detailed ambiguous base counting

Usage:
```bash
sbatch transcriptome_analysis.sh <output_directory> <fasta_file1> [<fasta_file2> ...]
```
Input Requirements
    File formats:
        Uncompressed FASTA (.fasta, .fa)
        Gzipped FASTA (.fasta.gz, .fa.gz)

Output Files
File Pattern                            Description
transcriptome_analysis_table.txt        Main summary table
*.seqkit_stats.txt                      Comprehensive SeqKit statistics
*.per_sequence_stats.txt                Per-contig/chromosome stats
*.n_counts.txt                          N/X counts per sequence
all_samples.*.txt                       Consolidated reports

```bash
#!/bin/bash
#SBATCH --job-name=sequences_stats
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G  # Increased for large files
#SBATCH --time=1:00:00  # Extended runtime
#SBATCH --partition=agap_short

echo "=== Starting Transcriptome Analysis ==="
echo "Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "Current directory: $(pwd)"

# Load required modules
module purge
module load bioinfo-cirad
module load seqkit/2.8.1

# Verify SeqKit is loaded
if ! command -v seqkit &> /dev/null; then
    echo "ERROR: SeqKit failed to load"
    module avail seqkit 2>&1
    exit 1
fi

# Usage instructions
if [ $# -lt 2 ]; then
    echo "Usage: $0 <output_directory> <fasta_file1> [<fasta_file2> ...]"
    echo "Example:"
    echo "  sbatch $0 ./results /path/to/file1.fasta /path/to/file2.fasta.gz"
    exit 1
fi

# Initialize directories
OUTPUT_DIR="$1"
mkdir -p "$OUTPUT_DIR" || { echo "Failed to create output directory"; exit 1; }

echo "=== Analysis Parameters ==="
echo "Output directory: ${OUTPUT_DIR}"
echo "Input files to process:"
printf '  - %s\n' "${@:2}"
echo ""

# Main output file
output_file="${OUTPUT_DIR}/transcriptome_analysis_table.txt"
echo -e "Sample\tContigs\tTotal_Length\tN_Count\tGC_Percent" > "$output_file"

# Process each FASTA file
for fasta_file in "${@:2}"; do
    base_name=$(basename "$fasta_file" .fasta | sed 's/.gz$//')
    echo "=== Processing: ${base_name} ==="

    # Verify file exists and is readable
    if [ ! -f "$fasta_file" ]; then
        echo "  ERROR: File not found - skipping"
        echo -e "${base_name}\tFile_not_found\tNA\tNA\tNA" >> "$output_file"
        continue
    fi

    if [ ! -r "$fasta_file" ]; then
        echo "  ERROR: File not readable - skipping"
        echo -e "${base_name}\tFile_not_readable\tNA\tNA\tNA" >> "$output_file"
        continue
    fi

    # Generate file-specific outputs
    stats_file="${OUTPUT_DIR}/${base_name}.seqkit_stats.txt"
    per_seq_file="${OUTPUT_DIR}/${base_name}.per_sequence_stats.txt"
    ncount_file="${OUTPUT_DIR}/${base_name}.n_counts.txt"

    # Run SeqKit stats (without path in output)
    echo "  Running comprehensive analysis..."
    if [[ "$fasta_file" == *.gz ]]; then
        seqkit stats -a -T <(zcat "$fasta_file") | awk -v bn="$base_name" 'NR==1 || NR==2 {$1=bn; print}' > "$stats_file" 2>> "${OUTPUT_DIR}/seqkit_errors.log"
    else
        seqkit stats -a -T "$fasta_file" | awk -v bn="$base_name" 'NR==1 || NR==2 {$1=bn; print}' > "$stats_file" 2>> "${OUTPUT_DIR}/seqkit_errors.log"
    fi

    # Generate per-sequence stats
    echo "  Generating per-sequence statistics..."
    if [[ "$fasta_file" == *.gz ]]; then
        seqkit fx2tab -n -l -g -i <(zcat "$fasta_file") > "$per_seq_file" 2>> "${OUTPUT_DIR}/seqkit_errors.log"
    else
        seqkit fx2tab -n -l -g -i "$fasta_file" > "$per_seq_file" 2>> "${OUTPUT_DIR}/seqkit_errors.log"
    fi

    # Count N/X bases
    echo "  Counting ambiguous bases..."
    if [[ "$fasta_file" == *.gz ]]; then
        zcat "$fasta_file" | seqkit fx2tab | awk -F'\t' 'BEGIN{OFS="\t"; total=0}
        {count = gsub(/[NnXx]/, "&", $2); print $1, count; total += count}
        END {print "Total", total > "'"${ncount_file}.total"'"}'> "$ncount_file"
    else
        seqkit fx2tab "$fasta_file" | awk -F'\t' 'BEGIN{OFS="\t"; total=0}
        {count = gsub(/[NnXx]/, "&", $2); print $1, count; total += count}
        END {print "Total", total > "'"${ncount_file}.total"'"}'> "$ncount_file"
    fi

    # Extract summary statistics
    echo "  Extracting summary metrics..."
    contig_count=$(grep -c "^>" <(if [[ "$fasta_file" == *.gz ]]; then zcat "$fasta_file"; else cat "$fasta_file"; fi) 2>/dev/null || echo "NA")
    total_length=$(awk 'NR==2 {print $5}' "$stats_file" 2>/dev/null || echo "NA")
    gc_content=$(awk 'NR==2 {print $NF}' "$stats_file" 2>/dev/null || echo "NA")
    n_count=$(awk '{sum += $2} END {print sum}' "$ncount_file" 2>/dev/null || echo "NA")

    # Append results
    echo -e "${base_name}\t${contig_count}\t${total_length}\t${n_count}\t${gc_content}" >> "$output_file"

    echo "  Completed: ${base_name}"
    echo "  Contigs: ${contig_count} | Length: ${total_length} | N-count: ${n_count} | GC%: ${gc_content}"
    echo ""
done

echo "=== Job Summary ==="
echo "Processed $(( $# - 1 )) files"
echo "Main results: ${output_file}"
echo "Detailed reports in: ${OUTPUT_DIR}/"
echo "Error logs: ${OUTPUT_DIR}/seqkit_errors.log"
echo "Job finished at: $(date)"
```

**CDS Position Calculator**

Description:
This Bash script calculates Coding Sequence (CDS) positions for transcripts based on:
- Transcript length statistics
- Annotated UTR (Untranslated Region) lengths
- Strand orientation information

Key Features
- **Strand-aware calculations**: Correctly handles + and - strand transcripts
- **Position validation**: Checks for invalid CDS ranges (start ≥ end)
- **Error reporting**: Identifies transcripts missing from UTR annotations
- **Flexible input**: Handles transcript IDs with pipe (`|`) delimiters
- **Header preservation**: Maintains original full transcript IDs

Usage
```bash
./CDS_positions.sh <transcript_stats> <utr_lengths> <output_file>
```

```bash
#!/bin/bash

# Usage: ./CDS_positions.sh <transcript_stats> <utr_lengths> <output_file>
# Example:
# ./CDS_positions.sh cleaned_sequence_stats.txt utr_lengths.tab cds_positions.tab

# Input files
STATS_FILE="$1"
UTR_FILE="$2"
OUTPUT_FILE="$3"

# Check args
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <transcript_stats> <utr_lengths> <output_file>"
    exit 1
fi

echo "Generating CDS positions table..."
echo "Transcript stats: $STATS_FILE"
echo "UTR lengths:      $UTR_FILE"
echo "Output file:      $OUTPUT_FILE"

awk -F'\t' '
BEGIN {
    OFS = "\t";
    print "full_transcript_id", "transcript_id", "cds_start", "cds_end", "strand", "five_prime_utr", "three_prime_utr", "original_length";
}
NR==FNR {
    # Read UTR file (skip header)
    if (FNR > 1) {
        utr[$1] = 1;
        five_utr[$1] = $2;
        three_utr[$1] = $3;
        strand[$1] = $4;
    }
    next;
}
{
    # Read transcript stats file
    full_id = $1;
    transcript_id = $1;  # Now matches UTR file IDs exactly
    original_length = $2;

    # Skip if transcript not in UTR file
    if (!(transcript_id in utr)) {
        print "Warning: Transcript " transcript_id " not found in UTR file" > "/dev/stderr";
        next;
    }

    # Calculate CDS positions
        cds_start = 200 + five_utr[transcript_id] + 1;
        cds_end = original_length - 200 - three_utr[transcript_id];

    # Verify positions are valid
    if (cds_start >= cds_end) {
        print "Warning: Invalid CDS positions for " transcript_id " (start >= end)" > "/dev/stderr";
        next;
    }

    # Output result
    print full_id, transcript_id, cds_start, cds_end, strand[transcript_id],
          five_utr[transcript_id], three_utr[transcript_id], original_length;
}' "$UTR_FILE" "$STATS_FILE" > "$OUTPUT_FILE"

echo "CDS positions saved to: $OUTPUT_FILE"

```

**CDS table for ORF_extractor**

```bash
awk 'NR>1 {gsub(/\./, "_", $1); print $1, $3, $4}' table_cds.txt > new_file.txt
sed 's/ /\t/g' new_file.txt > ORF_tab-delimited.tab
```

**Replace dots by underscore on the multi fasta for ORF_extractor**
```bash
awk '/^>/ {gsub(/\./, "_"); print; next} {print}' monococcum_high_coverage.fasta > monococcum_high_coverage_underscores.fasta
```

**ORF_extractor**
```bash
#!/bin/bash
#SBATCH --job-name=ORF_extraction
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --time=6:00:00
#SBATCH --partition=agap_normal

# ====== Use your personal Perl environment ======
#export PATH=$HOME/my_perl_5.36.0/bin:$PATH
#export PERL5LIB=$HOME/my_perl_5.36.0/lib/site_perl/5.36.0:$HOME/my_perl_5.36.0/lib/5.36.0

export PATH=$HOME/my_perl_5.36.0/bin:$PATH
export PERL5LIB=$HOME/my_perl_5.36.0/lib/perl5:$HOME/my_perl_5.36.0/lib/site_perl/5.36.0:$HOME/my_perl_5.36.0/lib/5.36.0

# ====== Arguments ======
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <ORF_table_path> <multi_fasta_path> <output_prefix>"
    exit 1
fi

ORF_TABLE="$1"
MULTI_FASTA="$2"
OUTPUT_PREFIX="$3"

echo "Starting ORF extraction job at $(date)"
echo "Using Perl: $(which perl)"
perl -MBio::SeqIO -e 'print "BioPerl is working in job environment\n"' || {
    echo "Error: BioPerl not found in job environment"
    exit 1
}

# ====== Run ORF extractor ======
$HOME/my_perl_5.36.0/bin/perl \
    /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/orf_extractor/ORF_extractor.pl \
    -orf "$ORF_TABLE" \
    -otype auto \
    -seq "$MULTI_FASTA" \
    -stype auto \
    -out "$OUTPUT_PREFIX"

if [ $? -eq 0 ]; then
    echo "ORF extraction completed successfully at $(date)"
else
    echo "Error: ORF extraction failed!"
    exit 1
fi
```

