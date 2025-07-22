#!/bin/bash

# Check for required arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/input_file.depth /path/to/output_directory"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_DIR="$2"

# Define output filenames
FILTERED_OUTPUT="${OUTPUT_DIR}/depth_per_position_covered_all_individuals.depth"
COUNT_OUTPUT="${OUTPUT_DIR}/contig_percentage_covered_positions.txt"
SUMMARY_OUTPUT="${OUTPUT_DIR}/number_covered_position.txt"

# Parameters (customize if needed)
MIN_DEPTH=10
MIN_INDIVIDUALS=7

# Validate input
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR" || {
    echo "Error: Failed to create output directory: $OUTPUT_DIR" >&2
    exit 1
}

echo "=== Starting Depth Coverage Analysis ==="
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Output files:"
echo "1. Filtered positions: $(basename "$FILTERED_OUTPUT")"
echo "2. Coverage statistics: $(basename "$COUNT_OUTPUT")"
echo "3. Summary table: $(basename "$SUMMARY_OUTPUT")"
echo "Parameters: Minimum depth = $MIN_DEPTH, Minimum individuals = $MIN_INDIVIDUALS"
echo ""

# Initialize output files
echo "Initializing output files..."
echo -e "contig\tpos\t$(head -1 "$INPUT_FILE" | cut -f3-)" > "$FILTERED_OUTPUT" || {
    echo "Error: Failed to initialize filtered output file" >&2
    exit 1
}

echo -e "Contig\tPositions_covered\tTotal_positions\tPercentage_covered" > "$COUNT_OUTPUT" || {
    echo "Error: Failed to initialize count output file" >&2
    exit 1
}

echo -e "total_positions\tcovered_positions\tpercentage_covered" > "$SUMMARY_OUTPUT" || {
    echo "Error: Failed to initialize summary output file" >&2
    exit 1
}

# Process with awk
echo "Processing depth data..."
awk -v filtered="$FILTERED_OUTPUT" \
    -v counts="$COUNT_OUTPUT" \
    -v min_depth="$MIN_DEPTH" \
    -v min_ind="$MIN_INDIVIDUALS" \
'
BEGIN {
    OFS = "\t";
    current_contig = "";
    contig_count = 0;
    total_positions_per_contig = 0;
    total_filtered_positions = 0;
    total_all_positions = 0;
}
NR == 1 { next } # Skip header
{
    total_all_positions++;
    contig = $1;
    pos = $2;
    count = 0;

    # Count individuals with sufficient depth
    for (i = 3; i <= NF; i++) {
        if ($i >= min_depth) count++;
    }

    # Handle contig changes
    if (contig != current_contig) {
        if (current_contig != "") {
            percentage = (contig_count / total_positions_per_contig) * 100;
            print current_contig, contig_count, total_positions_per_contig, percentage >> counts;
        }
        current_contig = contig;
        contig_count = 0;
        total_positions_per_contig = 0;
    }

    # Count positions
    total_positions_per_contig++;
    if (count >= min_ind) {
        contig_count++;
        total_filtered_positions++;
        print $0 >> filtered;
    }
}
END {
    # Write data for last contig
    if (current_contig != "") {
        percentage = (contig_count / total_positions_per_contig) * 100;
        print current_contig, contig_count, total_positions_per_contig, percentage >> counts;
    }

    # Calculate overall coverage percentage
    percentage_covered = (total_filtered_positions / total_all_positions) * 100;

    # Print summary statistics to a temporary file
    print "TOTAL_POSITIONS", total_all_positions > "summary.tmp";
    print "COVERED_POSITIONS", total_filtered_positions >> "summary.tmp";
    print "PERCENTAGE_COVERED", percentage_covered >> "summary.tmp";

    print "Processed " total_all_positions " total positions" > "/dev/stderr";
    print "Found " total_filtered_positions " covered positions (" percentage_covered "%)" > "/dev/stderr";
}
' "$INPUT_FILE" || {
    echo "Error: Failed to process input file" >&2
    exit 1
}

# Generate the summary table
echo "Generating summary table..."
total_positions=$(awk '/TOTAL_POSITIONS/{print $2}' summary.tmp)
covered_positions=$(awk '/COVERED_POSITIONS/{print $2}' summary.tmp)
percentage_covered=$(awk '/PERCENTAGE_COVERED/{printf "%.2f", $2}' summary.tmp)

echo -e "${total_positions}\t${covered_positions}\t${percentage_covered}" >> "$SUMMARY_OUTPUT"
rm -f summary.tmp

echo ""
echo "=== Analysis Complete ==="
echo "Results saved to:"
echo "1. Filtered positions: $FILTERED_OUTPUT ($(wc -l < "$FILTERED_OUTPUT") lines)"
echo "2. Coverage statistics: $COUNT_OUTPUT ($(wc -l < "$COUNT_OUTPUT") lines)"
echo "3. Summary table: $SUMMARY_OUTPUT"
echo ""
echo "Final coverage:"
echo "- Total positions: $total_positions"
echo "- Covered positions: $covered_positions"
echo "- Percentage covered: ${percentage_covered}%"
