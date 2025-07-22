#!/bin/bash
#SBATCH --job-name=bz2_to_gz
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --partition=agap_short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=15G
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-44%10
echo  [$(date +'%d-%m-%Y')]-[$(date +'%H:%M:%S')]-[$SLURM_JOB_ID]-[Start] $0 $@ >> '.logs'
if [ -d ./your-repository/.git ]; then git add .logs; git commit -m 'Update .logs'; fi
# --------------------------------------------------------------------------------- #
# Set Up variables
folder_path=$1
finalpath=$2
echo "Folder : $folder_path ; finalpath : $finalpath"

# Generate an array
mapfile -t files_in_folder < <(ls "$folder_path")

# Find target
filename="$folder_path/${files_in_folder[$SLURM_ARRAY_TASK_ID - 1]}"

# Do the target is a file ?
if [ ! -f "$filename" ]; then
    echo "$filename is not a regular file"
    exit 1
fi

# Do the target have the right format
if [[ ! "$filename" == *.bz2 ]]; then
   if [[ "$filename" == *.gz ]]; then
       cp $filename $finalpath
       echo ".gz found, cp only"
       exit 0
   fi

   echo "$filename should end with .bz2"
   exit 2
fi

# Set up more variabes
dir_path=$(dirname "$(pwd "$filename")")
base_name=$(basename "$filename" .bz2)
intermediary="$finalpath/$base_name.temp"
newbz2="$intermediary.bz2"
tempfile="$intermediary.verif"
finalpath="$finalpath/$base_name.gz"

# Do the name used by the result is used
if [ -f "$finalpath" ]; then
    echo "$finalpath already exists"
    exit 3
fi

echo "INI: $filename's weight: $(stat -c %s "$filename")"

# Do an intermediary allready exist ?
if [ -f "$intermediary" ]; then
    echo "$intermediary already exists." >&2
    exit 5
else
    bunzip2 -kcd "$filename" > "$intermediary"
fi

# Build the intermdiary
echo "TEMP: $intermediary 's weight: $(stat -c %s "$intermediary")"
gzip -c "$intermediary" > "$finalpath"

# Build the gz file
echo "TEMP: $finalpath content inspection"

echo "DONE: $finalpath's weight: $(stat -c %s "$finalpath")"

# Inspect the content of the gz fule
zcat "$finalpath" > "$tempfile"

echo "TEMP: $tempfile 's weight: $(stat -c %s "$tempfile")"

if diff "$intermediary" "$tempfile" >/dev/null; then
   echo "TEMP: FastQ: Good"
else
   echo "TEMP: FastQ: WARNING UNINDENTICAL FILE !!
   echo "TEMP: FastQ: WARNING UNINDENTICAL FILE !! >&2
   exit 10
fi

rm "$tempfile"

# Verify that the intermediary has correctly been extracted
bzip2 -c "$intermediary" > "$newbz2"

echo "TEMP: $newbz2 's weight: $(stat -c %s "$newbz2")"

if diff "$newbz2" "$filename" >/dev/null; then
   echo "END: bz2: Good"
else
   echo "END: bz2:  WARNING UNINDENTICAL FILE !!"
   echo "END: bz2:  WARNING UNINDENTICAL FILE !!" >&2
   exit 11
fi

rm $newbz2
rm $intermediary

echo  [$(date +'%d-%m-%Y')]-[$(date +'%H:%M:%S')]-[$SLURM_JOB_ID]-[Start] $0 $@
