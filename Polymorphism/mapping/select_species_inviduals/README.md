# FASTQ File Selector Script

**Description**
This script copies FASTQ files for specified species from a source directory to a target folder,
using data from `individuals_info.txt`.

**Usage**
./select_individuals.sh "<exact_species_name>" "/path/to/target_directory"

**Examples**
For species without spaces:
./select_individuals.sh speltoides ./speltoides_samples

For species with spaces (MUST use quotes):
./select_individuals.sh "Aegilops bicornis" ./bicornis_samples
./select_individuals.sh "Triticum monococcum ssp. boeoticum" ./boeoticum_samples
