#!/bin/bash

#SBATCH --job-name=extract-window
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00-1:0:00
#SBATCH --mem=10000M
#SBATCH --partition=low_p

FILE_IN="pietzner_2021_PMID34648354/files.txt"
FILE_OUT="pietzner_2021_PMID34648354/window/"
WINDOW=1000000
POS_COLUMN=3

# Ensure output directory exists
mkdir -p "$FILE_OUT"

# Check if the input file exists and is readable
if [[ ! -f "$FILE_IN" ]]; then
    echo "Input file $FILE_IN not found."
    exit 1
fi

# Read each line from the input file, skipping the header
while IFS=$'\t' read -r file_path ID CHR START END; do
    # Remove quotes from the 'CHR' column if present
    CHR=$(echo "$CHR" | tr -d '"')

    # Debugging statement for each line read
    echo "# START"
    echo "## FILE: $file_path"
    echo "## UNIPROT: $ID"
    echo "## CHR: $CHR"
    echo "## START: $START"
    echo "## END: $END"
    echo "## WINDOW: $WINDOW"

    # Skip the header or any invalid entries
    if [[ "$file_path" == "file" || -z "$CHR" || -z "$START" || -z "$END" || ! "$CHR" =~ ^[0-9XY]+$ ]]; then
        echo "** Skipping header or invalid entry **"
        continue
    fi
    
    # Ensure that START and END are integers
    if ! [[ "$START" =~ ^[0-9]+$ ]] || ! [[ "$END" =~ ^[0-9]+$ ]]; then
    echo "** Invalid START or END values. Skipping line. **"
    continue
    fi

    # Remove any quotes from the file path if present
    file_path=$(echo "$file_path" | tr -d '"')

    # Check if the file exists
    if [[ -f "$file_path" ]]; then
        # Extend the range by 1,000,000 on both sides
        START=$((START - WINDOW))
        END=$((END + WINDOW))

        # Extract the base filename (without path and extension)
        base_filename=$(basename "$file_path" .gz)
        # Remove everything after the first dot (.)
        base_filename="${base_filename%%.*}"
        # Construct the output file path
        output_file="$FILE_OUT/${base_filename}_UNIPROT-${ID}.txt"

        # Display progress
        echo "## PROCESSING"

        # Capture the header and filter the file based on CHR, START, and END
        # The first line (header) will be saved separately
        header=$(zcat "$file_path" | head -n 1)

        # Write the header to the output file first
        echo "$header" > "$output_file"

        # Filter the file and append the results to the output file
        zcat "$file_path" | awk -v CHR="$CHR" -v start="$START" -v end="$END" -v pos_col="$POS_COLUMN" \
            'BEGIN {OFS="\t"} NR > 1 {if ($1 == CHR && $pos_col >= start && $pos_col <= end) print $0}' >> "$output_file"

        # Check if the output file contains data
        echo "## SAVING"
        if [[ -s "$output_file" ]]; then
            # Count the number of lines in the new file (excluding the header)
            line_count=$(($(wc -l < "$output_file") - 1))
            echo "** Filtered rows for $base_filename saved to $output_file **"
            echo "** N rows: $line_count **"
        else
            echo "** No data found for $file_path within the specified range. Output file not created. **"
            rm -f "$output_file"  # Remove empty file
        fi
    else
        echo "File not found: $file_path"
    fi
done < "$FILE_IN"

echo "## DONE"
