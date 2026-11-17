#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 <population> <output.parquet> [GRCh37-region]" >&2
  echo "Population must be one of: AFR AMR CSA EAS EUR MID" >&2
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage
  exit 2
fi

population=$1
output=$2
region="${3:-}"

case "$population" in
  AFR|AMR|CSA|EAS|EUR|MID) ;;
  *)
    echo "Unsupported Pan-UKB population: $population" >&2
    usage
    exit 2
    ;;
esac

if [[ "$output" != *.parquet ]]; then
  echo "Output path must end in .parquet: $output" >&2
  exit 2
fi
if [[ -n "$region" && ! "$region" =~ ^[^:]+:[0-9]+-[0-9]+$ ]]; then
  echo "Invalid GRCh37 region: $region" >&2
  usage
  exit 2
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required to build the Pan-UKB variant index." >&2
  exit 1
fi

for path in "$output" "${output}.tmp" "${output}.export.tsv.bgz"; do
  if [[ -e "$path" ]]; then
    echo "Refusing to overwrite existing output or temporary file: $path" >&2
    exit 1
  fi
done

mkdir -p "$(dirname "$output")"
builder="$(dirname "${BASH_SOURCE[0]}")/build_panukb_variant_index.py"
if [[ -n "$region" ]]; then
  python3 "$builder" "$population" "$output" "$region"
else
  python3 "$builder" "$population" "$output"
fi
