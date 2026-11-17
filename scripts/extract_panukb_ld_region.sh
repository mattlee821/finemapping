#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 <population> <variant-index.parquet> <GRCh37-region> <output-prefix>" >&2
  echo "Population must be one of: AFR AMR CSA EAS EUR MID" >&2
  echo "Region format: contig:start-end (for example 22:25000000-25500000)" >&2
}

if [[ $# -ne 4 ]]; then
  usage
  exit 2
fi

population=$1
variant_index=$2
region=$3
output_prefix=$4

case "$population" in
  AFR|AMR|CSA|EAS|EUR|MID) ;;
  *)
    echo "Unsupported Pan-UKB population: $population" >&2
    usage
    exit 2
    ;;
esac

if [[ ! -f "$variant_index" ]]; then
  echo "Variant index does not exist: $variant_index" >&2
  exit 1
fi
if [[ ! "$region" =~ ^[^:]+:[0-9]+-[0-9]+$ ]]; then
  echo "Invalid GRCh37 region: $region" >&2
  usage
  exit 2
fi
if ! command -v ldcov >/dev/null 2>&1; then
  echo "ldcov is required; install it with its S3 extra (ldcov[s3])." >&2
  exit 1
fi

for suffix in .npz .variants.tsv .bcor .bcor.idx; do
  if [[ -e "${output_prefix}${suffix}" ]]; then
    echo "Refusing to overwrite existing output: ${output_prefix}${suffix}" >&2
    exit 1
  fi
done

mkdir -p "$(dirname "$output_prefix")"
block_matrix="s3://pan-ukb-us-east-1/ld_release/UKBB.${population}.ldadj.bm"
ldcov --ld-bm \
  --bm "$block_matrix" \
  --variant-index "$variant_index" \
  --region "$region" \
  --out "$output_prefix" \
  --output-format npz
