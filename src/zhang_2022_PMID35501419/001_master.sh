#!/bin/bash
#SBATCH --job-name=finemap-zhang
#SBATCH --array=0-4657
#SBATCH --time=03-00:0:00
#SBATCH --mem=15000M
#SBATCH --cpus-per-task=1
#SBATCH --partition=low_p
#SBATCH --output=log/%A_%a.out
#SBATCH --error=log/%A_%a.err

# Path variables
VAR_study="zhang_2022_PMID35501419"
FILE_LIST="/data/IARC_Biostat/work/finemapping/src/${VAR_study}/files.txt"

# Get the file corresponding to this array task
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")

# Skip if no file
if [ -z "$FILE" ]; then
  echo "No file for task ID ${SLURM_ARRAY_TASK_ID}, exiting."
  exit 0
fi

echo "Processing file: $FILE"

# Run your processing script
/opt/R/4.4.1/bin/Rscript "/data/IARC_Biostat/work/finemapping/src/${VAR_study}/001_master.R" "$FILE"
