#!/usr/bin/env bash
# Reproducible benchmark & server Pan-UKB parity validator for GWASfinemapping
# (mirrors GWASprocessing/src/benchmark.sh + implements update.md section 3 server checks).
#
# Usage:
#   bash scripts/benchmark.sh                              # run deterministic offline benchmark & compare to benchmark/expected/
#   bash scripts/benchmark.sh --update-expected            # regenerate benchmark/expected/expected_metrics.json
#   bash scripts/benchmark.sh --server-panukb <vcf.gz>     # run live server PyArrow-cut vs Hail parity + MHC memory check

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/env.sh
source "${SCRIPT_DIR}/lib/env.sh"
gwasfm_activate_runtime

cd "${GWASFM_ROOT}"

UPDATE_EXPECTED=0
SERVER_VCF=""
LD_POP="EUR"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --update-expected) UPDATE_EXPECTED=1; shift ;;
    --server-panukb) SERVER_VCF="$2"; shift 2 ;;
    --ld-pop) LD_POP="$2"; shift 2 ;;
    -h|--help) sed -n '2,8p' "${BASH_SOURCE[0]}"; exit 0 ;;
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
    # macOS /usr/bin/time -l reports bytes on modern macOS
    mb = val / (1024.0 * 1024.0) if val > 10_000_000 else val / 1024.0
    print(f"{mb:.1f}")
    raise SystemExit(0)
print("0")
PY
}

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

  echo "--- [Check 1] PyArrow cut_regional_index run (--max-loci 2) ---"
  mkdir -p "${TMP_SRV}/cache_cut" "${TMP_SRV}/out_cut" "${TMP_SRV}/empty_panukb"
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
      --max-loci 2 \
      > "${TMP_SRV}/cut.stdout"

  if grep -qiE "Initializing Hail|SparkUI|Running on Apache Spark" "${TMP_SRV}/cut.time" "${TMP_SRV}/cut.stdout"; then
    echo "FAIL: Hail/Spark was unexpectedly invoked during cut_regional_index mode!" >&2
    exit 1
  fi

  echo "--- [Check 2] Fallback regional Hail S3 build (--panukb-dir empty, --max-loci 2) ---"
  mkdir -p "${TMP_SRV}/cache_hail" "${TMP_SRV}/out_hail"
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
    --max-loci 2 \
    > "${TMP_SRV}/hail.stdout"

  PEAK_CUT_MB="$(extract_peak_mb "${TMP_SRV}/cut.time")"
  python3 - "${TMP_SRV}" "${LD_POP}" "${PEAK_CUT_MB}" <<'PY'
import json, sys
from pathlib import Path

root = Path(sys.argv[1])
pop = sys.argv[2]
peak_mb = sys.argv[3]
qc_cut = json.loads((root / "out_cut/srv_bench/cut_mode/qc.json").read_text())
qc_hail = json.loads((root / "out_hail/srv_bench/hail_mode/qc.json").read_text())

cut_pqs = sorted((root / "cache_cut/panukb").glob(f"UKBB.{pop}.chr*.variant_index.parquet"))
assert len(cut_pqs) >= 1, "No regional cut parquet files found in cache_cut/panukb/"
for pq_f in cut_pqs:
    assert Path(f"{pq_f}.json").is_file(), f"Missing JSON provenance sidecar for {pq_f}"

for loc_c, loc_h in zip(qc_cut["loci"], qc_hail["loci"]):
    assert loc_c["region_id"] == loc_h["region_id"], (loc_c["region_id"], loc_h["region_id"])
    assert loc_c["variants_ld_matched"] == loc_h["variants_ld_matched"], (
        f"MISMATCH in match_variants for {loc_c['region_id']}: cut={loc_c['variants_ld_matched']} vs hail={loc_h['variants_ld_matched']}"
    )
    assert loc_c["index_cut_seconds"] < 30.0, f"Index cut took {loc_c['index_cut_seconds']}s (>30s)"
    print(
        f"PASS locus {loc_c['region_id']}: matched={loc_c['variants_ld_matched']} "
        f"(exact parity with Hail) | cut_index={loc_c['index_cut_seconds']}s | locus_wall={loc_c['wall_seconds']}s"
    )
print(f"Peak RSS (cut mode): {peak_mb} MB")
print("=== ALL 3 SERVER PARITY CHECKS PASSED ===")
PY
  exit 0
fi

echo "=== Deterministic Offline Benchmark ==="
TMP_BENCH="$(mktemp -d "${TMPDIR:-/tmp}/gwasfm-bench-XXXXXX")"
trap 'rm -rf "${TMP_BENCH}"' EXIT

python3 scripts/build_synthetic_fixture.py --out "${TMP_BENCH}/fixture" >/dev/null

T0="$(python3 -c 'import time; print(time.time())')"
GWASFM_NXF_WORK="${TMP_BENCH}/nxf-work" \
GWASFM_CACHE_DIR="${TMP_BENCH}/cache" \
measure_rss_and_run "${TMP_BENCH}/bench.time" \
  bash scripts/finemap.sh \
    --manifest "${TMP_BENCH}/fixture/manifest.tsv" \
    --reference-root "${TMP_BENCH}/fixture/reference" \
    --local-ld-panel "${TMP_BENCH}/fixture/reference/local_eur_panel.npz" \
    --ld-pop EUR \
    --outdir "${TMP_BENCH}/out"
T1="$(python3 -c 'import time; print(time.time())')"
PEAK_MB="$(extract_peak_mb "${TMP_BENCH}/bench.time")"

EXPECTED_JSON="${GWASFM_ROOT}/benchmark/expected/expected_metrics.json"
mkdir -p "$(dirname "${EXPECTED_JSON}")"

Rscript - "${TMP_BENCH}/out" "${TMP_BENCH}/metrics.json" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1L]]
metrics_out <- args[[2L]]

collect_gwas <- function(rel) {
  gdir <- file.path(out_dir, rel)
  qc <- jsonlite::fromJSON(file.path(gdir, "qc.json"), simplifyVector = FALSE)
  e_rss <- new.env(parent = emptyenv())
  e_rrss <- new.env(parent = emptyenv())
  load(file.path(gdir, "methods/susie-rss.rda"), envir = e_rss)
  load(file.path(gdir, "methods/susie-rrss.rda"), envir = e_rrss)

  regs <- list()
  for (rid in sort(names(e_rss$susie_rss$regions))) {
    v_rss <- e_rss$susie_rss$regions[[rid]]$result$variants
    v_rrss <- e_rrss$susie_rrss$regions[[rid]]$result$variants
    top_i <- which.max(v_rss$PIP)
    regs[[rid]] <- list(
      variants = nrow(v_rss),
      cs_rss = nrow(e_rss$susie_rss$regions[[rid]]$result$credible_sets),
      cs_rrss = nrow(e_rrss$susie_rrss$regions[[rid]]$result$credible_sets),
      top_variant = v_rss$GWAS_VARIANT_ID[[top_i]],
      top_pip_rss = round(as.numeric(v_rss$PIP[[top_i]]), 5),
      top_pip_rrss = round(as.numeric(v_rrss$PIP[[top_i]]), 5)
    )
  }
  total_bytes <- sum(file.info(list.files(gdir, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE)
  list(
    overall_status = qc$overall_status,
    analyzed_regions = qc$analyzed_regions_count,
    total_cs_rss = qc$total_credible_sets_susie_rss,
    total_cs_rrss = qc$total_credible_sets_susie_rrss,
    published_bytes = total_bytes,
    regions = regs
  )
}

res <- list(
  P00001_olink = collect_gwas("sun_synthetic/P00001_olink"),
  colon_rectal_sim = collect_gwas("cancer_synthetic/colon_rectal_sim")
)
jsonlite::write_json(res, metrics_out, auto_unbox = TRUE, pretty = TRUE)
RS

python3 - "${TMP_BENCH}/metrics.json" "${EXPECTED_JSON}" "${UPDATE_EXPECTED}" "${T0}" "${T1}" "${PEAK_MB}" <<'PY'
import json, math, sys
from pathlib import Path

actual = json.loads(Path(sys.argv[1]).read_text())
expected_path = Path(sys.argv[2])
update = int(sys.argv[3]) == 1
wall_s = float(sys.argv[5]) - float(sys.argv[4])
peak_mb = float(sys.argv[6])

# Strip runtime-variable byte counts from deterministic statistical snapshot, but enforce <= 5 MB budget per synthetic GWAS
deterministic = {}
for k, v in actual.items():
    assert v["published_bytes"] < 5 * 1024 * 1024, f"Published folder {k} exceeded 5 MB budget: {v['published_bytes']}"
    deterministic[k] = { kk: vv for kk, vv in v.items() if kk != "published_bytes" }

if update or not expected_path.is_file():
    expected_path.write_text(json.dumps(deterministic, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote expected benchmark snapshot to {expected_path}")
else:
    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    assert set(deterministic) == set(expected), (set(deterministic), set(expected))
    for gwas_k in expected:
        exp_g = expected[gwas_k]
        act_g = deterministic[gwas_k]
        for field in ("overall_status", "analyzed_regions", "total_cs_rss", "total_cs_rrss"):
            assert act_g[field] == exp_g[field], f"{gwas_k}.{field}: {act_g[field]} != {exp_g[field]}"
        assert set(act_g["regions"]) == set(exp_g["regions"]), f"{gwas_k} regions mismatch"
        for rid, exp_r in exp_g["regions"].items():
            act_r = act_g["regions"][rid]
            for f in ("variants", "cs_rss", "cs_rrss", "top_variant"):
                assert act_r[f] == exp_r[f], f"{gwas_k}.{rid}.{f}: {act_r[f]} != {exp_r[f]}"
            for f in ("top_pip_rss", "top_pip_rrss"):
                assert math.isclose(act_r[f], exp_r[f], abs_tol=1e-4), f"{gwas_k}.{rid}.{f}: {act_r[f]} != {exp_r[f]}"

for k, v in actual.items():
    print(
        f"  {k:20s} | status={v['overall_status']} | regions={v['analyzed_regions']} "
        f"| CS={v['total_cs_rss']}/{v['total_cs_rrss']} | disk={v['published_bytes']/1024:.1f} KB"
    )
print(f"Wall time: {wall_s:.2f}s | Peak RSS: {peak_mb:.1f} MB")
print("=== BENCHMARK CONCORDANCE PASSED ===")
PY
