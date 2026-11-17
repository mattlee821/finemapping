#!/bin/bash
#SBATCH --job-name=gwasfm_reference
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
# Install reference data as a job (a whole-population Pan-UKB index takes hours). Submit from the repository root:
#   sbatch scripts/install_reference_data.sbatch [options of install_reference_data.sh]
# 1 population per job runs the 6 in parallel:
#   for p in AFR AMR CSA EAS EUR MID; do sbatch --job-name=gwasfm_$p scripts/install_reference_data.sbatch --only panukb --populations $p; done
# The final check covers only what the job was asked to install (its --only and --populations).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
bash scripts/install_reference_data.sh install "$@"
check_args=()
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
  case "${args[i]}" in --only|--populations) check_args+=("${args[i]}" "${args[i + 1]:-}") ;; esac
done
bash scripts/install_reference_data.sh check ${check_args[@]+"${check_args[@]}"}
