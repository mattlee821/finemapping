#!/usr/bin/env bash
# Install the reference data a run needs (Pan-UKB variant indexes, GRCh37 FASTA, liftover chain, gene annotation).
#
# Usage:
#   bash scripts/install_reference_data.sh                      # everything, all 6 Pan-UKB populations
#   bash scripts/install_reference_data.sh --only genome,annotation
#   bash scripts/install_reference_data.sh --only panukb --populations EUR,AFR
#   bash scripts/install_reference_data.sh check                # report what is missing or changed
#   bash scripts/install_reference_data.sh paths --nextflow     # the installed files as nextflow run options
# `install` is the default command. The reference root is GWASFM_REFERENCE_ROOT from .env (default
# ~/.cache/gwasfinemapping/reference). Needs the tools environment (bash scripts/source.sh --tools-only first) and
# network access to UCSC, Ensembl, GitHub and the public Pan-UKB S3 bucket. See docs/reference-data.md.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case "${1:-}" in -h|--help) sed -n '2,12p' "${BASH_SOURCE[0]}"; exit 0 ;; esac
# shellcheck disable=SC1091
source "${ROOT_DIR}/scripts/source.sh"   # loads .env and activates the tools environment
command="install"
case "${1:-}" in install|check|paths) command="$1"; shift ;; esac
exec python3 "${ROOT_DIR}/scripts/reference_data.py" "$command" "$@"
