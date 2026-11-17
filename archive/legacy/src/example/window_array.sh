#!/bin/bash

#SBATCH --job-name=extract-window
#SBATCH --array=0-XXX
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00-1:0:00
#SBATCH --mem=10000M
#SBATCH --partition=low_p

FILE_IN="pietzner_2021_PMID34648354/files.txt"
FILE_OUT="pietzner_2021_PMID34648354/window/"
WINDOW=1000000
POS_COLUMN=3

mkdir -p "$FILE_OUT"

# Get the current line from files.txt (skip header with NR==SLURM_ARRAY_TASK_ID+1)
line=$(awk -v n=$SLURM_ARRAY_TASK_ID 'NR==n+1 {print}' "$FILE_IN")

# Parse columns from the line
file_path=$(echo "$line" | cut -f1 | tr -d '"')
ID=$(echo "$line"        | cut -f2)
CHR=$(echo "$line"       | cut -f3 | tr -d '"')
START=$(echo "$line"     | cut -f4)
END=$(echo "$line"       | cut -f5)

# Skip invalid entries
if [[ -z "$file_path" || -z "$CHR" || -z "$START" || -z "$END" || ! "$CHR" =~ ^[0-9XY]+$ ]]; then
    echo "Invalid or header line for task $SLURM_ARRAY_TASK_ID. Skipping."
    exit 0
fi

if ! [[ "$START" =~ ^[0-9]+$ ]] || ! [[ "$END" =~ ^[0-9]+$ ]]; then
    echo "Invalid START/END values on line $SLURM_ARRAY_TASK_ID. Skipping."
    exit 0
fi

if [[ ! -f "$file_path" ]]; then
    echo "File not found: $file_path"
    exit 0
fi

# Expand window
START=$((START - WINDOW))
END=$((END + WINDOW))

# Construct output filename
base_filename=$(basename "$file_path" .gz)
base_filename="${base_filename%%.*}"
output_file="$FILE_OUT/${base_filename}_UNIPROT-${ID}.txt"

# Extract header
header=$(zcat "$file_path" | head -n 1)
echo "$header" > "$output_file"

# Filter
zcat "$file_path" | awk -v CHR="$CHR" -v start="$START" -v end="$END" -v pos_col="$POS_COLUMN" \
    'BEGIN {OFS="\t"} NR > 1 {if ($1 == CHR && $pos_col >= start && $pos_col <= end) print $0}' >> "$output_file"

# Report
if [[ -s "$output_file" ]]; then
    line_count=$(($(wc -l < "$output_file") - 1))
    echo "Task $SLURM_ARRAY_TASK_ID: Saved $line_count rows for $base_filename to $output_file"
else
    echo "Task $SLURM_ARRAY_TASK_ID: No data found. Removing $output_file"
    rm -f "$output_file"
fi
