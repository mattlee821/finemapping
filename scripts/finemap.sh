#!/usr/bin/env bash
#SBATCH --job-name=GWASfinemapping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=slurm-gwasfm-%j.out
#
# Unified local / SLURM entrypoint for GWASfinemapping (mirrors GWASprocessing/src/standardise.sh).
# Usage:
#   bash scripts/finemap.sh [options]
#   sbatch scripts/finemap.sh [options] --slurm

set -euo pipefail

find_root() {
  if [[ -n "${GWASFM_ROOT_OVERRIDE:-}" ]]; then
    (cd "${GWASFM_ROOT_OVERRIDE}" && pwd); return
  fi
  local script_root candidate
  script_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  if [[ -f "${script_root}/finemapping.nf" ]]; then
    printf '%s\n' "${script_root}"; return
  fi
  for candidate in "${SLURM_SUBMIT_DIR:-}" "${PWD}"; do
    if [[ -n "${candidate}" && -f "${candidate}/finemapping.nf" ]]; then
      (cd "${candidate}" && pwd); return
    fi
  done
  echo "Cannot locate GWASfinemapping root (finemapping.nf)." >&2
  exit 2
}

ROOT_DIR="$(find_root)"
# shellcheck source=lib/env.sh
source "${ROOT_DIR}/scripts/lib/env.sh"
gwasfm_activate_runtime

usage() {
  cat <<'EOF'
Usage:
  bash scripts/finemap.sh [options]
  sbatch scripts/finemap.sh [options]

Inputs (supply one of --input-dir, --manifest, or --vcf):
  --input-dir DIR            GWASprocessing root, study standardise/ folder, or single GWAS folder
  --manifest PATH            TSV manifest with 'vcf' or 'gwas_dir' column
  --vcf PATH                 Single GRCh38 GWAS-VCF (.vcf.gz)
  --outdir DIR               Output directory (default: $GWASFM_WORK_ROOT -> <repo>/finemapping)
  --temp-dir DIR             Scratch/temp root for Nextflow workdirs, run logs, locus matrix.npy, and TMPDIR (default: $GWASFM_TEMP_ROOT)

Single-GWAS / filter options:
  --only-study STUDY_ID      Filter --input-dir or --manifest to one study_id
  --study-id ID              Study identifier override for --vcf
  --gwas-id ID               GWAS identifier override for --vcf
  --metadata PATH            Optional metadata.json override for --vcf
  --qc PATH                  Optional qc.json override for --vcf
  --trait NAME               Optional phenotype label
  --ancestry CODE            Optional GWAS ancestry code
  --gene-id ID               Optional UniProt ID or HGNC gene symbol for cis_gene

Fine-mapping & reference options:
  --methods LIST             Methods to run, space/comma-separated (default: "p-ld susie-rss susie-rrss carma polyfun")
  --output-format FMT        Method file format: rda|rds|both (default: rda)
  --compress COMP            Compression for saved method files: xz|gzip (default: xz)
  --ld-pop POP               Pan-UKB LD population: AFR|AMR|CSA|EAS|EUR|MID (default: EUR)
  --reference-root DIR       Reference root (default: $GWASFM_REFERENCE_ROOT)
  --bm-path DIR              Optional local BlockMatrix directory override (defaults to live S3)
  --exclude-mhc              Exclude extended MHC (chr6:25.5-34.0 Mb GRCh38)
  --significance-lp FLOAT    Significance gate in -log10(P) (default: 7.30103 = P <= 5e-8)
  --cis-significance-lp FLT  Optional cis_gene significance gate override
  --region-flank-bp INT      Greedy lead flank in bp (default: 1500000 = 3 Mb window)
  --max-loci INT             Cap analyzed regions per GWAS (0 = all; default: 0)

Batch control & execution (matching GWASprocessing):
  --row-limit INT            Eligible rows to run (0 = all; default: 0)
  --force                    Re-run even when complete outputs exist for current settings
  --dry-run                  Validate and summarize inputs without running fine-mapping
  --report true|false        Render per-GWAS report.html (default: true)
  --keep-intermediates       Retain per-locus scratch LD/VCF/TSV files
  --allow-ancestry-mismatch  Allow mismatched single-ancestry LD
  --profile NAME             Nextflow profile (default: slurm under sbatch, else local)
  --slurm                    Alias for --profile slurm
  --slurm-exclude NODES      Comma-separated SLURM nodes to exclude (default: $GWASFM_SLURM_EXCLUDE)
  --container ENGINE         Run tasks inside container: apptainer|docker
  --max-forks INT            Max concurrent GWAS jobs (default: $GWASFM_MAX_FORKS or 2)
  --resume                   Pass -resume to Nextflow
EOF
}

INPUT_DIR=""
MANIFEST=""
VCF=""
OUTDIR="${GWASFM_WORK_ROOT}"
TEMP_DIR_OPT=""
ONLY_STUDY=""
STUDY_ID=""
GWAS_ID=""
METADATA=""
QC=""
TRAIT=""
ANCESTRY=""
GENE_ID=""
METHODS_SPEC="p_ld,susie_rss,susie_rrss,carma,polyfun"
OUTPUT_FORMAT="rda"
COMPRESS_MODE="xz"
LD_POP="EUR"
REF_ROOT="${GWASFM_REFERENCE_ROOT}"
BM_PATH=""
EXCLUDE_MHC="false"
SIG_LP=""
CIS_SIG_LP=""
FLANK_BP=""
MAX_LOCI=""
ROW_LIMIT="0"
FORCE="false"
DRY_RUN="false"
REPORT="true"
KEEP_INTERMEDIATES="false"
ALLOW_ANC_MISMATCH="false"
LOCAL_LD_PANEL=""
PROFILE=""
SLURM_EXCLUDE="${GWASFM_SLURM_EXCLUDE:-}"
CONTAINER_ENGINE=""
MAX_FORKS=""
RESUME=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir|--input_dir) INPUT_DIR="$2"; shift 2 ;;
    --manifest) MANIFEST="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --temp-dir|--temp_dir) TEMP_DIR_OPT="$2"; shift 2 ;;
    --only-study|--only_study) ONLY_STUDY="$2"; shift 2 ;;
    --study-id|--study_id) STUDY_ID="$2"; shift 2 ;;
    --gwas-id|--gwas_id) GWAS_ID="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --qc) QC="$2"; shift 2 ;;
    --trait) TRAIT="$2"; shift 2 ;;
    --ancestry) ANCESTRY="$2"; shift 2 ;;
    --gene-id|--gene_id) GENE_ID="$2"; shift 2 ;;
    --methods) METHODS_SPEC="$2"; shift 2 ;;
    --output-format|--output_format) OUTPUT_FORMAT="$2"; shift 2 ;;
    --compress) COMPRESS_MODE="$2"; shift 2 ;;
    --ld-pop|--population) LD_POP="$2"; shift 2 ;;
    --reference-root) REF_ROOT="$2"; shift 2 ;;
    --bm-path|--bm_path) BM_PATH="$2"; shift 2 ;;
    --exclude-mhc|--exclude_mhc) EXCLUDE_MHC="true"; shift ;;
    --significance-lp) SIG_LP="$2"; shift 2 ;;
    --cis-significance-lp) CIS_SIG_LP="$2"; shift 2 ;;
    --region-flank-bp) FLANK_BP="$2"; shift 2 ;;
    --max-loci|--max_loci) MAX_LOCI="$2"; shift 2 ;;
    --row-limit|--row_limit) ROW_LIMIT="$2"; shift 2 ;;
    --force) FORCE="true"; shift ;;
    --dry-run|--dry_run) DRY_RUN="true"; shift ;;
    --report) REPORT="$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]')"; shift 2 ;;
    --keep-intermediates) KEEP_INTERMEDIATES="true"; shift ;;
    --allow-ancestry-mismatch) ALLOW_ANC_MISMATCH="true"; shift ;;
    --local-ld-panel) LOCAL_LD_PANEL="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --slurm) PROFILE="slurm"; shift ;;
    --slurm-exclude|--slurm_exclude) SLURM_EXCLUDE="$2"; shift 2 ;;
    --container) CONTAINER_ENGINE="$2"; shift 2 ;;
    --max-forks|--max_forks) MAX_FORKS="$2"; shift 2 ;;
    --resume|-resume) RESUME=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$VCF" && -z "$MANIFEST" && -z "$INPUT_DIR" ]]; then
  echo "Error: supply one of --input-dir, --manifest, or --vcf." >&2
  usage >&2
  exit 2
fi
if [[ "$REPORT" != "true" && "$REPORT" != "false" ]]; then
  echo "Error: --report must be true or false (got $REPORT)." >&2
  exit 2
fi

resolve_abs() {
  [[ "$1" = /* ]] && printf '%s\n' "$1" || printf '%s\n' "${PWD}/$1"
}

OUTDIR="$(resolve_abs "$OUTDIR")"
[[ -n "$TEMP_DIR_OPT" ]] && TEMP_DIR_OPT="$(resolve_abs "$TEMP_DIR_OPT")"
REF_ROOT="$(resolve_abs "$REF_ROOT")"
[[ -n "$INPUT_DIR" ]] && INPUT_DIR="$(resolve_abs "$INPUT_DIR")"
[[ -n "$MANIFEST" ]] && MANIFEST="$(resolve_abs "$MANIFEST")"
[[ -n "$VCF" ]] && VCF="$(resolve_abs "$VCF")"
[[ -n "$METADATA" ]] && METADATA="$(resolve_abs "$METADATA")"
[[ -n "$QC" ]] && QC="$(resolve_abs "$QC")"
[[ -n "$LOCAL_LD_PANEL" ]] && LOCAL_LD_PANEL="$(resolve_abs "$LOCAL_LD_PANEL")"
[[ -n "$BM_PATH" ]] && BM_PATH="$(resolve_abs "$BM_PATH")"

CHAIN="${REF_ROOT}/genome/hg38ToHg19.over.chain.gz"
FASTA="${REF_ROOT}/genome/GRCh37.fa"
FAI="${REF_ROOT}/genome/GRCh37.fa.fai"
GENE_ANNOT="${REF_ROOT}/annotation/gene_regions_GRCh38.tsv"
PANUKB_DIR="${REF_ROOT}/panukb"
POLYFUN_DIR=""
for pf_cand in "${REF_ROOT}/polyfun" "${GWASFM_CACHE_DIR}/reference/polyfun" "${HOME}/.cache/gwasfinemapping/reference/polyfun"; do
  if [[ -d "${pf_cand}/by_chrom" ]]; then
    POLYFUN_DIR="${pf_cand}"
    break
  fi
done

for req in "$CHAIN" "$FASTA" "$FAI"; do
  if [[ ! -f "$req" ]]; then
    echo "Missing reference file: $req (run 'bash scripts/install_reference_data.sh install' or pass --reference-root)" >&2
    exit 1
  fi
done

if [[ -z "$PROFILE" ]]; then
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    PROFILE="slurm"
  else
    PROFILE="${GWASFM_PROFILE:-local}"
  fi
fi
if [[ -n "$CONTAINER_ENGINE" ]]; then
  PROFILE="${PROFILE},${CONTAINER_ENGINE}"
elif [[ -n "${GWASFM_CONTAINER:-}" ]]; then
  PROFILE="${PROFILE},${GWASFM_CONTAINER}"
fi

RUN_ID="${SLURM_JOB_ID:-$(date +%Y%m%d-%H%M%S)-$$}"
TEMP_ROOT="${TEMP_DIR_OPT:-${GWASFM_TEMP_ROOT:-}}"
if [[ -n "$TEMP_ROOT" ]]; then
  TEMP_ROOT="$(resolve_abs "$TEMP_ROOT")"
  SCRATCH_ROOT="${TEMP_ROOT}/scratch"
  LOG_ROOT="${TEMP_ROOT}/logs/runs/${RUN_ID}"
  RUN_TMP_DIR="${TEMP_ROOT}/tmp/${RUN_ID}"
  BM_CACHE_DIR="${GWASFM_BM_CACHE_DIR:-${TEMP_ROOT}/bm_cache}"
  if [[ -n "$TEMP_DIR_OPT" ]]; then
    export GWASFM_NXF_WORK="${TEMP_ROOT}/nextflow-work"
  fi
  export TMPDIR="${RUN_TMP_DIR}"
  export NXF_TEMP="${RUN_TMP_DIR}"
  export NXF_CACHE_DIR="${GWASFM_NXF_WORK}/nxf-cache"
  mkdir -p "$OUTDIR" "$LOG_ROOT" "${GWASFM_NXF_WORK}" "${NXF_CACHE_DIR}" "${SCRATCH_ROOT}" "${RUN_TMP_DIR}" "${BM_CACHE_DIR}"
else
  SCRATCH_ROOT="${GWASFM_SCRATCH_DIR:-${GWASFM_CACHE_DIR}/scratch}"
  LOG_ROOT="${GWASFM_CACHE_DIR}/logs/runs/${RUN_ID}"
  RUN_TMP_DIR="${GWASFM_TMP_DIR:-}"
  BM_CACHE_DIR="${GWASFM_BM_CACHE_DIR:-}"
  export NXF_CACHE_DIR="${GWASFM_NXF_WORK}/nxf-cache"
  mkdir -p "$OUTDIR" "$LOG_ROOT" "${GWASFM_NXF_WORK}" "${NXF_CACHE_DIR}"
fi
NXF_LOG="${LOG_ROOT}/nextflow.${RUN_ID}.log"
export NXF_OPTS="${NXF_OPTS:--Xms2g -Xmx10g}"

# Fast preflight path on host: skip Nextflow JVM startup for --dry-run or when all rows are already complete
SKIP_NXF=0
if [[ -z "$CONTAINER_ENGINE" && -z "${GWASFM_CONTAINER:-}" ]]; then
  PREP_CMD=(
    python3 "${ROOT_DIR}/scripts/prepare_manifest.py" prepare
    --outdir "$OUTDIR"
    --state-file "${OUTDIR}/.finemap_state.tsv"
    --output "${LOG_ROOT}/run_manifest.tsv"
    --summary "${LOG_ROOT}/run_manifest_summary.json"
    --rows-dir "${LOG_ROOT}/rows"
    --chain "$CHAIN"
    --grch37-fasta "$FASTA"
    --cache-dir "${GWASFM_CACHE_DIR}"
    --population "$LD_POP"
    --methods "$METHODS_SPEC"
    --output-format "$OUTPUT_FORMAT"
    --significance-lp "${SIG_LP:-7.301029995664}"
    --region-flank-bp "${FLANK_BP:-1500000}"
    --max-region-width-bp 9500000
    --gene-window-bp 1000000
    --exclude-mhc "$EXCLUDE_MHC"
    --max-loci "${MAX_LOCI:-0}"
    --max-causal 10
    --max-iter 50
    --row-limit "$ROW_LIMIT"
    --force "$FORCE"
    --dry-run "$DRY_RUN"
    --report "$REPORT"
    --keep-intermediates "$KEEP_INTERMEDIATES"
    --allow-ancestry-mismatch "$ALLOW_ANC_MISMATCH"
  )
  [[ -f "$GENE_ANNOT" ]] && PREP_CMD+=(--gene-annotation "$GENE_ANNOT")
  [[ -d "$PANUKB_DIR" ]] && PREP_CMD+=(--panukb-dir "$PANUKB_DIR")
  [[ -n "$POLYFUN_DIR" ]] && PREP_CMD+=(--polyfun-dir "$POLYFUN_DIR")
  [[ -n "$INPUT_DIR" ]] && PREP_CMD+=(--input-dir "$INPUT_DIR")
  [[ -n "$MANIFEST" ]] && PREP_CMD+=(--manifest "$MANIFEST")
  [[ -n "$VCF" ]] && PREP_CMD+=(--vcf "$VCF")
  [[ -n "$ONLY_STUDY" ]] && PREP_CMD+=(--only-study "$ONLY_STUDY")
  [[ -n "$STUDY_ID" ]] && PREP_CMD+=(--study-id "$STUDY_ID")
  [[ -n "$GWAS_ID" ]] && PREP_CMD+=(--gwas-id "$GWAS_ID")
  [[ -n "$METADATA" ]] && PREP_CMD+=(--metadata "$METADATA")
  [[ -n "$QC" ]] && PREP_CMD+=(--qc "$QC")
  [[ -n "$TRAIT" ]] && PREP_CMD+=(--trait "$TRAIT")
  [[ -n "$ANCESTRY" ]] && PREP_CMD+=(--ancestry "$ANCESTRY")
  [[ -n "$GENE_ID" ]] && PREP_CMD+=(--gene-id "$GENE_ID")
  [[ -n "$CIS_SIG_LP" ]] && PREP_CMD+=(--cis-significance-lp "$CIS_SIG_LP")
  [[ -n "$LOCAL_LD_PANEL" ]] && PREP_CMD+=(--local-ld-panel "$LOCAL_LD_PANEL")
  [[ -n "$BM_PATH" ]] && PREP_CMD+=(--bm-path "$BM_PATH")

  if "${PREP_CMD[@]}" >/dev/null 2>&1; then
    if [[ "$DRY_RUN" == "true" ]]; then
      SKIP_NXF=1
    elif python3 -c "import json,sys; d=json.load(open(sys.argv[1])); raise SystemExit(0 if d.get('rows_to_run',1)==0 and not d.get('invalid') else 1)" "${LOG_ROOT}/run_manifest_summary.json"; then
      python3 "${ROOT_DIR}/scripts/prepare_manifest.py" update-state \
        --state-file "${OUTDIR}/.finemap_state.tsv" \
        --outdir "$OUTDIR" \
        --output "${LOG_ROOT}/state_update_summary.json" >/dev/null 2>&1 || true
      SKIP_NXF=1
    fi
  fi
fi

[[ -n "$MAX_FORKS" ]] && export GWASFM_MAX_FORKS="$MAX_FORKS"
NF_CMD=(
  nextflow
  -log "$NXF_LOG"
  run "${ROOT_DIR}/finemapping.nf"
  -profile "$PROFILE"
  -work-dir "${GWASFM_NXF_WORK}"
  -with-trace "${LOG_ROOT}/trace.tsv"
  -with-report "${LOG_ROOT}/nextflow-report.html"
  -with-timeline "${LOG_ROOT}/nextflow-timeline.html"
  --outdir "$OUTDIR"
  --log_root "$LOG_ROOT"
  --scratch_dir "$SCRATCH_ROOT"
  --chain "$CHAIN"
  --grch37_fasta "$FASTA"
  --grch37_fai "$FAI"
  --cache_dir "${GWASFM_CACHE_DIR}"
  --population "$LD_POP"
  --methods "$METHODS_SPEC"
  --output_format "$OUTPUT_FORMAT"
  --compress "$COMPRESS_MODE"
  --exclude_mhc "$EXCLUDE_MHC"
  --row_limit "$ROW_LIMIT"
  --force "$FORCE"
  --dry_run "$DRY_RUN"
  --report "$REPORT"
  --keep_intermediates "$KEEP_INTERMEDIATES"
  --allow_ancestry_mismatch "$ALLOW_ANC_MISMATCH"
)

[[ -f "$GENE_ANNOT" ]] && NF_CMD+=(--gene_annotation "$GENE_ANNOT")
[[ -d "$PANUKB_DIR" ]] && NF_CMD+=(--panukb_dir "$PANUKB_DIR")
[[ -n "$POLYFUN_DIR" ]] && NF_CMD+=(--polyfun_dir "$POLYFUN_DIR")
[[ -n "$INPUT_DIR" ]] && NF_CMD+=(--input_dir "$INPUT_DIR")
[[ -n "$MANIFEST" ]] && NF_CMD+=(--manifest "$MANIFEST")
[[ -n "$VCF" ]] && NF_CMD+=(--vcf "$VCF")
[[ -n "$ONLY_STUDY" ]] && NF_CMD+=(--only_study "$ONLY_STUDY")
[[ -n "$STUDY_ID" ]] && NF_CMD+=(--study_id "$STUDY_ID")
[[ -n "$GWAS_ID" ]] && NF_CMD+=(--gwas_id "$GWAS_ID")
[[ -n "$METADATA" ]] && NF_CMD+=(--metadata "$METADATA")
[[ -n "$QC" ]] && NF_CMD+=(--qc "$QC")
[[ -n "$TRAIT" ]] && NF_CMD+=(--trait "$TRAIT")
[[ -n "$ANCESTRY" ]] && NF_CMD+=(--ancestry "$ANCESTRY")
[[ -n "$GENE_ID" ]] && NF_CMD+=(--gene_id "$GENE_ID")
[[ -n "$SIG_LP" ]] && NF_CMD+=(--significance_lp "$SIG_LP")
[[ -n "$CIS_SIG_LP" ]] && NF_CMD+=(--cis_significance_lp "$CIS_SIG_LP")
[[ -n "$FLANK_BP" ]] && NF_CMD+=(--region_flank_bp "$FLANK_BP")
[[ -n "$MAX_LOCI" ]] && NF_CMD+=(--max_loci "$MAX_LOCI")
[[ -n "$MAX_FORKS" ]] && NF_CMD+=(--max_forks "$MAX_FORKS")
[[ -n "$SLURM_EXCLUDE" ]] && NF_CMD+=(--slurm_exclude "$SLURM_EXCLUDE")
[[ -n "$LOCAL_LD_PANEL" ]] && NF_CMD+=(--local_ld_panel "$LOCAL_LD_PANEL")
[[ -n "$BM_PATH" ]] && NF_CMD+=(--bm_path "$BM_PATH")
[[ -n "$RUN_TMP_DIR" ]] && NF_CMD+=(--temp_dir "$RUN_TMP_DIR")
[[ -n "$BM_CACHE_DIR" ]] && NF_CMD+=(--bm_cache_dir "$BM_CACHE_DIR")
[[ "$RESUME" -eq 1 ]] && NF_CMD+=(-resume)

nf_status=0
if [[ "$SKIP_NXF" -eq 0 ]]; then
  "${NF_CMD[@]}" || nf_status=$?
fi

summary_status=0
python3 - "$LOG_ROOT" "$DRY_RUN" <<'PY' || summary_status=$?
import json, sys
from pathlib import Path

log_root = Path(sys.argv[1])
dry_run = sys.argv[2] == "true"
prep_file = log_root / "run_manifest_summary.json"
if not prep_file.is_file():
    raise SystemExit(0)

prep = json.loads(prep_file.read_text(encoding="utf-8"))
state_file = log_root / "state_update_summary.json"
upd = json.loads(state_file.read_text(encoding="utf-8")) if state_file.is_file() else {}
invalid = prep.get("invalid") or []
failed_rows = upd.get("failed_rows") or []

print("\nGWASfinemapping summary")
print(f"  Source:                 {prep.get('source', '')}")
print(f"  Output directory:       {prep.get('outdir', '')}")
print(f"  Rows to run:            {prep.get('rows_to_run', 0)}")
print(f"  Rows skipped (complete):{prep.get('rows_skipped_complete', 0)}")
if prep.get("rows_deferred_due_to_limit", 0):
    print(f"  Deferred by row limit:  {prep['rows_deferred_due_to_limit']}")
if invalid:
    print(f"  Failed preflight:       {len(invalid)}")
    for item in invalid:
        print(f"    - {item.get('study_id')}/{item.get('gwas_id')}: {'; '.join(item.get('errors', []))}")
if not dry_run:
    rows_expected = int(prep.get("rows_to_run", 0) or 0)
    updated_n = int(upd.get("updated", 0) or 0)
    unaccounted_n = max(0, rows_expected - updated_n - len(failed_rows))
    print(f"  Completed this run:     {updated_n}")
    if failed_rows:
        print(f"  Failed during run:      {len(failed_rows)}")
        for fr in failed_rows:
            print(f"    - {fr.get('row_key')}: {fr.get('error')}")
    if unaccounted_n > 0:
        print(f"  Failed upstream (screen/task): {unaccounted_n}")
    if rows_expected == 0 and prep.get("rows_skipped_complete", 0) > 0 and not invalid:
        print("  All selected rows are already complete for current settings (use --force to re-run).")

if (invalid or failed_rows or (not dry_run and unaccounted_n > 0)) and not dry_run:
    raise SystemExit(1)
PY

if [[ "$DRY_RUN" != "true" && -f "${LOG_ROOT}/run_manifest.tsv" ]]; then
  slurm_out=""
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_out="${SLURM_SUBMIT_DIR:-${PWD}}/slurm-gwasfm-${SLURM_JOB_ID}.out"
  fi
  python3 "${ROOT_DIR}/scripts/prepare_manifest.py" copy-logs \
    --outdir "$OUTDIR" \
    --run-manifest "${LOG_ROOT}/run_manifest.tsv" \
    --nextflow-log "$NXF_LOG" \
    --trace-file "${LOG_ROOT}/trace.tsv" \
    ${slurm_out:+--slurm-out "$slurm_out"} || true
fi

if [[ "$RESUME" -eq 0 && "$KEEP_INTERMEDIATES" != "true" ]]; then
  # Note: ${SCRATCH_ROOT}/by_gwas/<study>__<gwas>__<hash> is cleaned per-GWAS on completion
  # by finalise_gwas.py so unfinished GWAS keep already-completed loci across runs.
  rm -rf "${SCRATCH_ROOT}/${RUN_ID}" 2>/dev/null || true
  if [[ -n "$TEMP_ROOT" && -n "$RUN_TMP_DIR" ]]; then
    rm -rf "${RUN_TMP_DIR}" 2>/dev/null || true
  fi
  if [[ "$nf_status" -eq 0 && "$summary_status" -eq 0 ]]; then
    rm -rf "${GWASFM_NXF_WORK}"/[0-9a-f][0-9a-f] 2>/dev/null || true
  fi
fi

[[ "$nf_status" -ne 0 ]] && exit "$nf_status"
exit "$summary_status"
