#!/usr/bin/env bash
# Validate an installation by fine-mapping the two reproducible benchmark GWAS
# (Locke 2015 BMI + Sun 2023 UKB-PPP NPPB P16860_OID20049) against public Pan-UKB EUR LD
# and comparing variant-by-variant results against benchmark/expected/ (mirrors GWASprocessing/src/benchmark.sh).
#
# Usage:
#   bash scripts/benchmark.sh                              # real-data benchmark (skips already-completed rows unless --force)
#   bash scripts/benchmark.sh --force                      # force re-run from scratch
#   bash scripts/benchmark.sh --gwas-root /path/to/GWAS    # explicit GWASprocessing output root
#   bash scripts/benchmark.sh --update-expected            # write/update benchmark/expected/*.json
#   bash scripts/benchmark.sh --offline                    # fast zero-network synthetic benchmark (for CI)
#   bash scripts/benchmark.sh --server-panukb <vcf.gz>     # 3-point PyArrow-cut vs Hail parity + peak RSS check

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/env.sh
source "${SCRIPT_DIR}/lib/env.sh"
gwasfm_activate_runtime
cd "${GWASFM_ROOT}"

FORCE_FLAG=()
UPDATE_EXPECTED=0
OFFLINE=0
SERVER_VCF=""
GWAS_ROOT="${GWASFM_BENCHMARK_GWAS_ROOT:-}"
BENCH_OUTDIR="${GWASFM_BENCHMARK_OUTDIR:-${GWASFM_CACHE_DIR}/benchmark-output}"
MAX_LOCI=2
LD_POP="EUR"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE_FLAG=(--force); shift ;;
    --update-expected) UPDATE_EXPECTED=1; shift ;;
    --offline) OFFLINE=1; shift ;;
    --server-panukb) SERVER_VCF="$2"; shift 2 ;;
    --gwas-root) GWAS_ROOT="$2"; shift 2 ;;
    --outdir) BENCH_OUTDIR="$2"; shift 2 ;;
    --max-loci) MAX_LOCI="$2"; shift 2 ;;
    --ld-pop) LD_POP="$2"; shift 2 ;;
    -h|--help) sed -n '2,13p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

measure_rss_and_run() {
  local rss_file="$1"; shift
  if [[ "$(uname -s)" == "Darwin" ]] && [[ -x /usr/bin/time ]]; then
    /usr/bin/time -l "$@" 2> >(tee "${rss_file}" >&2)
  elif [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v "$@" 2> >(tee "${rss_file}" >&2)
  else
    "$@"
  fi
}

extract_peak_mb() {
  local rss_file="$1"
  [[ -f "$rss_file" ]] || { echo "0"; return; }
  python3 - "$rss_file" <<'PY'
import re, sys
txt = open(sys.argv[1], errors="ignore").read()
m_linux = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", txt)
if m_linux:
    print(f"{int(m_linux.group(1)) / 1024.0:.1f}")
    raise SystemExit(0)
m_mac = re.search(r"(\d+)\s+maximum resident set size", txt)
if m_mac:
    val = int(m_mac.group(1))
    mb = val / (1024.0 * 1024.0) if val > 10_000_000 else val / 1024.0
    print(f"{mb:.1f}")
    raise SystemExit(0)
print("0")
PY
}

# Mode 1: Server 3-point PyArrow-cut vs Hail parity check (update.md section 3)
if [[ -n "$SERVER_VCF" ]]; then
  echo "=== Live Server Pan-UKB Parity & Memory Benchmark (${SERVER_VCF}) ==="
  bash scripts/install_reference_data.sh check
  TMP_SRV="$(mktemp -d "${TMPDIR:-/tmp}/gwasfm-srvbench-XXXXXX")"
  trap 'rm -rf "${TMP_SRV}"' EXIT

  REF_ROOT="${GWASFM_REFERENCE_ROOT}"
  CHAIN="${REF_ROOT}/genome/hg38ToHg19.over.chain.gz"
  FASTA="${REF_ROOT}/genome/GRCh37.fa"
  GENE_ANNOT="${REF_ROOT}/annotation/gene_regions_GRCh38.tsv"
  PANUKB_DIR="${REF_ROOT}/panukb"

  mkdir -p "${TMP_SRV}/cache_cut" "${TMP_SRV}/out_cut" "${TMP_SRV}/empty_panukb" "${TMP_SRV}/cache_hail" "${TMP_SRV}/out_hail"
  measure_rss_and_run "${TMP_SRV}/cut.time" \
    python3 scripts/run_finemapping.py \
      --vcf "${SERVER_VCF}" \
      --study-id srv_bench \
      --gwas-id cut_mode \
      --chain "${CHAIN}" \
      --grch37-fasta "${FASTA}" \
      --gene-annotation "${GENE_ANNOT}" \
      --panukb-dir "${PANUKB_DIR}" \
      --cache-dir "${TMP_SRV}/cache_cut" \
      --output-root "${TMP_SRV}/out_cut" \
      --population "${LD_POP}" \
      --max-loci 2 > "${TMP_SRV}/cut.stdout"

  if grep -qiE "Initializing Hail|SparkUI|Running on Apache Spark" "${TMP_SRV}/cut.time" "${TMP_SRV}/cut.stdout"; then
    echo "FAIL: Hail/Spark was unexpectedly invoked during cut_regional_index mode!" >&2
    exit 1
  fi

  python3 scripts/run_finemapping.py \
    --vcf "${SERVER_VCF}" \
    --study-id srv_bench \
    --gwas-id hail_mode \
    --chain "${CHAIN}" \
    --grch37-fasta "${FASTA}" \
    --gene-annotation "${GENE_ANNOT}" \
    --panukb-dir "${TMP_SRV}/empty_panukb" \
    --cache-dir "${TMP_SRV}/cache_hail" \
    --output-root "${TMP_SRV}/out_hail" \
    --population "${LD_POP}" \
    --max-loci 2 > "${TMP_SRV}/hail.stdout"

  PEAK_CUT_MB="$(extract_peak_mb "${TMP_SRV}/cut.time")"
  python3 - "${TMP_SRV}" "${LD_POP}" "${PEAK_CUT_MB}" <<'PY'
import json, sys
from pathlib import Path

root = Path(sys.argv[1]); pop = sys.argv[2]; peak_mb = sys.argv[3]
qc_cut = json.loads((root / "out_cut/srv_bench/cut_mode/qc.json").read_text())
qc_hail = json.loads((root / "out_hail/srv_bench/hail_mode/qc.json").read_text())
cut_pqs = sorted((root / "cache_cut/panukb").glob(f"UKBB.{pop}.chr*.variant_index.parquet"))
assert len(cut_pqs) >= 1, "No regional cut parquet files found"
for pq_f in cut_pqs:
    assert Path(f"{pq_f}.json").is_file(), f"Missing JSON provenance sidecar for {pq_f}"
for loc_c, loc_h in zip(qc_cut["loci"], qc_hail["loci"]):
    assert loc_c["region_id"] == loc_h["region_id"]
    assert loc_c["variants_ld_matched"] == loc_h["variants_ld_matched"], (
        f"MISMATCH in match_variants for {loc_c['region_id']}: cut={loc_c['variants_ld_matched']} vs hail={loc_h['variants_ld_matched']}"
    )
    assert loc_c["index_cut_seconds"] < 30.0
    print(f"PASS {loc_c['region_id']}: matched={loc_c['variants_ld_matched']} | cut_index={loc_c['index_cut_seconds']}s")
print(f"Peak RSS (cut mode): {peak_mb} MB")
print("=== ALL 3 SERVER PARITY CHECKS PASSED ===")
PY
  exit 0
fi

# Mode 2: Offline synthetic benchmark (used in container CI where drive1/S3 are absent)
if [[ "$OFFLINE" -eq 1 ]]; then
  echo "=== Deterministic Offline Synthetic Benchmark ==="
  TMP_BENCH="$(mktemp -d "${TMPDIR:-/tmp}/gwasfm-bench-XXXXXX")"
  trap 'rm -rf "${TMP_BENCH}"' EXIT
  python3 scripts/build_synthetic_fixture.py --out "${TMP_BENCH}/fixture" >/dev/null
  GWASFM_NXF_WORK="${TMP_BENCH}/nxf-work" \
  GWASFM_CACHE_DIR="${TMP_BENCH}/cache" \
  bash scripts/finemap.sh \
    --manifest "${TMP_BENCH}/fixture/manifest.tsv" \
    --reference-root "${TMP_BENCH}/fixture/reference" \
    --local-ld-panel "${TMP_BENCH}/fixture/reference/local_eur_panel.npz" \
    --ld-pop EUR \
    --outdir "${TMP_BENCH}/out"
  MODE_ARG="check"
  [[ "$UPDATE_EXPECTED" -eq 1 ]] && MODE_ARG="write"
  Rscript scripts/check_benchmark.R --mode "$MODE_ARG" \
    --gwas-dir "${TMP_BENCH}/out/sun_synthetic/P00001_olink" \
    --expected "${GWASFM_ROOT}/benchmark/expected/offline_P00001_olink.json"
  Rscript scripts/check_benchmark.R --mode "$MODE_ARG" \
    --gwas-dir "${TMP_BENCH}/out/cancer_synthetic/colon_rectal_sim" \
    --expected "${GWASFM_ROOT}/benchmark/expected/offline_colon_rectal_sim.json"
  echo "  Result: PASS (offline synthetic benchmark)."
  exit 0
fi

# Mode 3 (Default): Real-data benchmark on Locke 2015 BMI + Sun 2023 UKB-PPP P16860_OID20049 (NPPB)
if [[ -z "$GWAS_ROOT" ]]; then
  for cand in \
    "/Volumes/drive1/research/GWASprocessing/data/GWAS" \
    "${HOME}/GWASprocessing/GWAS" \
    "${HOME}/drive1/research/GWASprocessing/data/GWAS"; do
    if [[ -d "$cand/locke_2015_25673413/standardise" && -d "$cand/sun_2023_37794186/standardise" ]]; then
      GWAS_ROOT="$cand"
      break
    fi
  done
fi

if [[ -z "$GWAS_ROOT" || ! -d "$GWAS_ROOT" ]]; then
  echo "Error: could not find standardised Locke 2015 and Sun 2023 GWAS directories." >&2
  echo "Pass --gwas-root /path/to/GWASprocessing/data/GWAS (or --offline for synthetic CI check)." >&2
  exit 1
fi

LOCKE_DIR="${GWAS_ROOT}/locke_2015_25673413/standardise/GWAScatalog_GCST002783"
SUN_DIR="${GWAS_ROOT}/sun_2023_37794186/standardise/P16860_OID20049"
for d in "$LOCKE_DIR" "$SUN_DIR"; do
  if [[ ! -d "$d" ]]; then
    echo "Missing required benchmark GWAS folder: $d" >&2
    exit 1
  fi
done

TMP_MAN="$(mktemp -d "${TMPDIR:-/tmp}/gwasfm-realbench-XXXXXX")"
trap 'rm -rf "${TMP_MAN}"' EXIT
BENCH_MANIFEST="${TMP_MAN}/benchmark_manifest.tsv"
cat > "${BENCH_MANIFEST}" <<EOF
study_id	gwas_id	gwas_dir
locke_2015_25673413	GWAScatalog_GCST002783	${LOCKE_DIR}
sun_2023_37794186	P16860_OID20049	${SUN_DIR}
EOF

echo "=== GWASfinemapping Real-Data Benchmark (Locke 2015 BMI + Sun 2023 UKB-PPP NPPB, --max-loci ${MAX_LOCI}) ==="
T0="$(python3 -c 'import time; print(time.time())')"
measure_rss_and_run "${TMP_MAN}/bench.time" \
  bash scripts/finemap.sh \
    --manifest "${BENCH_MANIFEST}" \
    --ld-pop "${LD_POP}" \
    --max-loci "${MAX_LOCI}" \
    --max-forks 1 \
    --outdir "${BENCH_OUTDIR}" \
    ${FORCE_FLAG[@]+"${FORCE_FLAG[@]}"}
T1="$(python3 -c 'import time; print(time.time())')"
PEAK_MB="$(extract_peak_mb "${TMP_MAN}/bench.time")"

MODE_ARG="check"
[[ "$UPDATE_EXPECTED" -eq 1 ]] && MODE_ARG="write"

failed=0
echo ""
echo "GWASfinemapping benchmark checks:"
Rscript scripts/check_benchmark.R \
  --mode "$MODE_ARG" \
  --gwas-dir "${BENCH_OUTDIR}/locke_2015_25673413/GWAScatalog_GCST002783" \
  --expected "${GWASFM_ROOT}/benchmark/expected/locke_2015_GWAScatalog_GCST002783.json" || failed=1

Rscript scripts/check_benchmark.R \
  --mode "$MODE_ARG" \
  --gwas-dir "${BENCH_OUTDIR}/sun_2023_37794186/P16860_OID20049" \
  --expected "${GWASFM_ROOT}/benchmark/expected/sun_2023_P16860_OID20049.json" || failed=1

python3 - "${T0}" "${T1}" "${PEAK_MB}" "${BENCH_OUTDIR}" "${failed}" <<'PY'
import sys
wall = float(sys.argv[2]) - float(sys.argv[1])
peak = float(sys.argv[3])
outdir = sys.argv[4]
failed = int(sys.argv[5])
print(f"  Wall time: {wall:.1f}s | Peak RSS: {peak:.1f} MB")
print(f"  Manifest report: {outdir}/manifest.html")
if failed:
    print("  Result: FAIL.")
    raise SystemExit(1)
print("  Result: PASS. This install reproduces the stored Locke 2015 BMI + Sun 2023 UKB-PPP benchmark.")
PY
