#!/usr/bin/env bash
#SBATCH --job-name=gwasfm
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=slurm-gwasfm-%j.out
#
# Unified local / SLURM entrypoint for GWASfinemapping (mirrors GWASprocessing/src/standardise.sh).
# Can be invoked directly via `bash scripts/finemap.sh ...` or submitted via `sbatch scripts/finemap.sh ...`.
#
# Examples:
#   bash scripts/finemap.sh --vcf trait.vcf.gz --study-id ukb --gwas-id bmi --ld-pop EUR --outdir results/
#   sbatch scripts/finemap.sh --manifest manifest.tsv --ld-pop EUR --outdir results/ --slurm

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/scripts/lib/env.sh" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}/scripts"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
# shellcheck source=lib/env.sh
source "${SCRIPT_DIR}/lib/env.sh"
gwasfm_activate_runtime

usage() {
  cat <<'EOF'
Usage:
  bash scripts/finemap.sh [options]
  sbatch scripts/finemap.sh [options] --slurm

Required (one of --vcf or --manifest, plus --outdir):
  --vcf PATH                 Single GRCh38 GWAS-VCF (.vcf.gz)
  --manifest PATH            TSV manifest with column 'vcf' (and optional study_id, gwas_id, metadata, qc, gene_id)
  --outdir DIR               Output directory

Single-GWAS metadata options:
  --study-id ID              Study identifier (required with --vcf unless --metadata supplies it)
  --gwas-id ID               GWAS identifier (required with --vcf unless --metadata supplies it)
  --metadata PATH            Optional GWASprocessing metadata.json
  --qc PATH                  Optional GWASprocessing qc.json
  --trait NAME               Optional phenotype label
  --ancestry CODE            Optional GWAS ancestry check label
  --gene-id ID               Optional UniProt ID or HGNC gene symbol for cis_gene window

Fine-mapping & reference options:
  --ld-pop POP               Pan-UKB LD population: AFR|AMR|CSA|EAS|EUR|MID (default: EUR)
  --reference-root DIR       Reference data root (default: $GWASFM_REFERENCE_ROOT)
  --exclude-mhc              Exclude extended MHC (chr6:25.5-34.0 Mb GRCh38)
  --significance-lp FLOAT    Regional significance gate in -log10(P) (default: 7.30103 = P <= 5e-8)
  --cis-significance-lp FLT  Optional cis_gene significance gate override
  --region-flank-bp INT      Greedy lead flank in bp (default: 1500000 = 3 Mb window)
  --max-loci INT             Cap analyzed regions per GWAS (0 = all; default: 0)
  --keep-intermediates       Retain scratch per-locus LD/VCF/TSV files for debugging

Execution options:
  --slurm                    Submit tasks via Nextflow SLURM profile
  --container ENGINE         Run tasks in container: apptainer|docker
  --max-forks INT            Max concurrent GWAS jobs (default: $GWASFM_MAX_FORKS or 2)
  --resume                   Pass -resume to Nextflow
EOF
}

VCF=""
MANIFEST=""
OUTDIR=""
STUDY_ID=""
GWAS_ID=""
METADATA=""
QC=""
TRAIT=""
ANCESTRY=""
GENE_ID=""
LD_POP="EUR"
REF_ROOT="${GWASFM_REFERENCE_ROOT}"
EXCLUDE_MHC="false"
SIG_LP=""
CIS_SIG_LP=""
FLANK_BP=""
MAX_LOCI=""
KEEP_INTERMEDIATES="false"
LOCAL_LD_PANEL=""
USE_SLURM=0
CONTAINER_ENGINE=""
MAX_FORKS=""
RESUME=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2 ;;
    --manifest) MANIFEST="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --study-id|--study_id) STUDY_ID="$2"; shift 2 ;;
    --gwas-id|--gwas_id) GWAS_ID="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --qc) QC="$2"; shift 2 ;;
    --trait) TRAIT="$2"; shift 2 ;;
    --ancestry) ANCESTRY="$2"; shift 2 ;;
    --gene-id|--gene_id) GENE_ID="$2"; shift 2 ;;
    --ld-pop|--population) LD_POP="$2"; shift 2 ;;
    --reference-root) REF_ROOT="$2"; shift 2 ;;
    --exclude-mhc|--exclude_mhc) EXCLUDE_MHC="true"; shift ;;
    --significance-lp) SIG_LP="$2"; shift 2 ;;
    --cis-significance-lp) CIS_SIG_LP="$2"; shift 2 ;;
    --region-flank-bp) FLANK_BP="$2"; shift 2 ;;
    --max-loci|--max_loci) MAX_LOCI="$2"; shift 2 ;;
    --keep-intermediates) KEEP_INTERMEDIATES="true"; shift ;;
    --local-ld-panel) LOCAL_LD_PANEL="$2"; shift 2 ;;
    --slurm) USE_SLURM=1; shift ;;
    --container) CONTAINER_ENGINE="$2"; shift 2 ;;
    --max-forks|--max_forks) MAX_FORKS="$2"; shift 2 ;;
    --resume|-resume) RESUME=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$OUTDIR" || ( -z "$VCF" && -z "$MANIFEST" ) ]]; then
  echo "Error: supply --outdir and either --vcf or --manifest." >&2
  usage >&2
  exit 2
fi

CHAIN="${REF_ROOT}/genome/hg38ToHg19.over.chain.gz"
FASTA="${REF_ROOT}/genome/GRCh37.fa"
FAI="${REF_ROOT}/genome/GRCh37.fa.fai"
GENE_ANNOT="${REF_ROOT}/annotation/gene_regions_GRCh38.tsv"
PANUKB_DIR="${REF_ROOT}/panukb"

for req in "$CHAIN" "$FASTA" "$FAI"; do
  if [[ ! -f "$req" ]]; then
    echo "Missing reference file: $req (run 'bash scripts/install_reference_data.sh install' or pass --reference-root)" >&2
    exit 1
  fi
done

PROFILES=("local")
if [[ "$USE_SLURM" -eq 1 ]]; then
  PROFILES=("slurm")
fi
if [[ -n "$CONTAINER_ENGINE" ]]; then
  PROFILES+=("$CONTAINER_ENGINE")
elif [[ -n "${GWASFM_CONTAINER:-}" ]]; then
  PROFILES+=("${GWASFM_CONTAINER}")
fi
PROFILE_CSV="$(IFS=,; echo "${PROFILES[*]}")"

NF_ARGS=(
  run "${GWASFM_ROOT}/finemapping.nf"
  -profile "$PROFILE_CSV"
  --outdir "$OUTDIR"
  --chain "$CHAIN"
  --grch37_fasta "$FASTA"
  --grch37_fai "$FAI"
  --population "$LD_POP"
  --exclude_mhc "$EXCLUDE_MHC"
  --keep_intermediates "$KEEP_INTERMEDIATES"
)

[[ -f "$GENE_ANNOT" ]] && NF_ARGS+=(--gene_annotation "$GENE_ANNOT")
[[ -d "$PANUKB_DIR" ]] && NF_ARGS+=(--panukb_dir "$PANUKB_DIR")
[[ -n "$VCF" ]] && NF_ARGS+=(--vcf "$VCF")
[[ -n "$MANIFEST" ]] && NF_ARGS+=(--manifest "$MANIFEST")
[[ -n "$STUDY_ID" ]] && NF_ARGS+=(--study_id "$STUDY_ID")
[[ -n "$GWAS_ID" ]] && NF_ARGS+=(--gwas_id "$GWAS_ID")
[[ -n "$METADATA" ]] && NF_ARGS+=(--metadata "$METADATA")
[[ -n "$QC" ]] && NF_ARGS+=(--qc "$QC")
[[ -n "$TRAIT" ]] && NF_ARGS+=(--trait "$TRAIT")
[[ -n "$ANCESTRY" ]] && NF_ARGS+=(--ancestry "$ANCESTRY")
[[ -n "$GENE_ID" ]] && NF_ARGS+=(--gene_id "$GENE_ID")
[[ -n "$SIG_LP" ]] && NF_ARGS+=(--significance_lp "$SIG_LP")
[[ -n "$CIS_SIG_LP" ]] && NF_ARGS+=(--cis_significance_lp "$CIS_SIG_LP")
[[ -n "$FLANK_BP" ]] && NF_ARGS+=(--region_flank_bp "$FLANK_BP")
[[ -n "$MAX_LOCI" ]] && NF_ARGS+=(--max_loci "$MAX_LOCI")
[[ -n "$MAX_FORKS" ]] && NF_ARGS+=(--max_forks "$MAX_FORKS")
[[ -n "$LOCAL_LD_PANEL" ]] && NF_ARGS+=(--local_ld_panel "$LOCAL_LD_PANEL")
[[ "$RESUME" -eq 1 ]] && NF_ARGS+=(-resume)

exec nextflow "${NF_ARGS[@]}"
