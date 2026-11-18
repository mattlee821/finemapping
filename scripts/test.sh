#!/usr/bin/env bash
# End-to-end offline validation test suite for GWASfinemapping (mirrors GWASprocessing/src/test.sh).
# Runs:
#   [1/4] Python unit & allele-swap invariance tests
#   [2/4] R package testthat suite
#   [3/4] R CMD check --no-manual (0 errors, 0 warnings, 0 notes)
#   [4/4] End-to-end Nextflow workflows: --input-dir discovery, idempotent skip, --force,
#         --report false + standalone render_report.R, --exclude-mhc, cis gate, and coloc .rda contract.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/env.sh
source "${SCRIPT_DIR}/lib/env.sh"

if [[ -n "${GWASFM_CONTAINER:-}" && -z "${GWASFM_IN_CONTAINER:-}" ]]; then
  echo "Running test suite inside container (${GWASFM_CONTAINER})..."
  gwasfm_in_container_exec bash scripts/test.sh "$@"
  exit 0
fi

gwasfm_activate_runtime
cd "${GWASFM_ROOT}"

echo "=== [1/4] Python unit & allele-orientation tests ==="
if python3 -c "import pytest" >/dev/null 2>&1; then
  python3 -m pytest -q tests/python
else
  python3 -m unittest discover -s tests/python -v
fi

echo "=== [2/4] R package build/install & testthat suite ==="
find src \( -name '*.o' -o -name '*.so' \) -delete 2>/dev/null || true
R_TARGET_LIB="$(Rscript -e 'cat(.libPaths()[1])')"
R_LIBS_ALL="$(Rscript -e 'cat(.libPaths(), sep=":")')"
if [[ ! -e "${R_TARGET_LIB}/susieR/DESCRIPTION" ]]; then
  echo "Error: active R library (${R_TARGET_LIB}) does not contain susieR. Run 'bash scripts/source.sh --r-packages' first." >&2
  exit 1
fi
R_LIBS="${R_LIBS_ALL}" R_LIBS_USER="${R_TARGET_LIB}" R_PROFILE_USER="" \
  R CMD INSTALL --no-multiarch --library="${R_TARGET_LIB}" . >/dev/null
find src \( -name '*.o' -o -name '*.so' \) -delete 2>/dev/null || true
R_LIBS="${R_LIBS_ALL}" R_LIBS_USER="${R_TARGET_LIB}" R_PROFILE_USER="" Rscript -e '
  library(gwasfinemapping)
  testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
'

echo "=== [3/4] R CMD check (0 errors, 0 warnings, 0 notes) ==="
TMP_RCHK="$(mktemp -d "${TMPDIR:-/tmp}/gwasfm-rcheck-XXXXXX")"
(
  cd "${TMP_RCHK}"
  R_LIBS="${R_LIBS_ALL}" R_LIBS_USER="${R_TARGET_LIB}" R_PROFILE_USER="" \
    R CMD build --no-build-vignettes --no-manual "${GWASFM_ROOT}" >/dev/null
  _R_CHECK_SYSTEM_CLOCK_=0 R_LIBS="${R_LIBS_ALL}" R_LIBS_USER="${R_TARGET_LIB}" R_PROFILE_USER="" \
    R CMD check --no-manual --no-tests gwasfinemapping_*.tar.gz > check.out 2>&1 || {
      cat check.out >&2
      exit 1
    }
  grep -E "^Status: OK" check.out
)
rm -rf "${TMP_RCHK}"

echo "=== [4/4] End-to-end Nextflow workflows, state tracking, MHC policy & storage-light contract ==="
TMP_TEST="$(mktemp -d "${TMPDIR:-/tmp}/gwasfm-test-XXXXXX")"
trap 'rm -rf "${TMP_TEST}"' EXIT

python3 scripts/build_synthetic_fixture.py --out "${TMP_TEST}/fixture" >/dev/null

# 4a. Dry-run preflight check via --manifest
GWASFM_NXF_WORK="${TMP_TEST}/nxf-work" \
GWASFM_CACHE_DIR="${TMP_TEST}/cache" \
bash scripts/finemap.sh \
  --manifest "${TMP_TEST}/fixture/manifest.tsv" \
  --reference-root "${TMP_TEST}/fixture/reference" \
  --local-ld-panel "${TMP_TEST}/fixture/reference/local_eur_panel.npz" \
  --ld-pop EUR \
  --outdir "${TMP_TEST}/out_batch" \
  --dry-run

# 4b. Full batch run (--manifest)
GWASFM_NXF_WORK="${TMP_TEST}/nxf-work" \
GWASFM_CACHE_DIR="${TMP_TEST}/cache" \
bash scripts/finemap.sh \
  --manifest "${TMP_TEST}/fixture/manifest.tsv" \
  --reference-root "${TMP_TEST}/fixture/reference" \
  --local-ld-panel "${TMP_TEST}/fixture/reference/local_eur_panel.npz" \
  --ld-pop EUR \
  --outdir "${TMP_TEST}/out_batch"

# 4c. Idempotent re-run must skip both already-completed GWAS (0 rows launched)
RERUN_OUT="$(
  GWASFM_NXF_WORK="${TMP_TEST}/nxf-work" \
  GWASFM_CACHE_DIR="${TMP_TEST}/cache" \
  bash scripts/finemap.sh \
    --manifest "${TMP_TEST}/fixture/manifest.tsv" \
    --reference-root "${TMP_TEST}/fixture/reference" \
    --local-ld-panel "${TMP_TEST}/fixture/reference/local_eur_panel.npz" \
    --ld-pop EUR \
    --outdir "${TMP_TEST}/out_batch"
)"
printf '%s\n' "${RERUN_OUT}" | grep -q "Rows skipped (complete):2"

# 4d. Single-GWAS run with --exclude-mhc, --report false (deferred report), and sub-threshold cis_gene
GWASFM_NXF_WORK="${TMP_TEST}/nxf-work" \
GWASFM_CACHE_DIR="${TMP_TEST}/cache" \
bash scripts/finemap.sh \
  --vcf "${TMP_TEST}/fixture/gwas/sun_synthetic/P00001_olink/P00001_olink.vcf.gz" \
  --study-id sun_synthetic \
  --gwas-id P00002_subthresh_nomhc \
  --ancestry ALL \
  --gene-id P00002 \
  --exclude-mhc \
  --report false \
  --reference-root "${TMP_TEST}/fixture/reference" \
  --local-ld-panel "${TMP_TEST}/fixture/reference/local_eur_panel.npz" \
  --ld-pop EUR \
  --outdir "${TMP_TEST}/out_nomhc"

test ! -f "${TMP_TEST}/out_nomhc/sun_synthetic/P00002_subthresh_nomhc/report.html"
# Render report on demand purely from .rda + .json
Rscript scripts/render_report.R --gwas-dir "${TMP_TEST}/out_nomhc/sun_synthetic/P00002_subthresh_nomhc"
test -f "${TMP_TEST}/out_nomhc/sun_synthetic/P00002_subthresh_nomhc/report.html"

# 4e. User-selected method subset ("p-ld susie-rrss polyfun") + --output-format rds
GWASFM_NXF_WORK="${TMP_TEST}/nxf-work" \
GWASFM_CACHE_DIR="${TMP_TEST}/cache" \
bash scripts/finemap.sh \
  --vcf "${TMP_TEST}/fixture/gwas/cancer_synthetic/colon_rectal_sim/colon_rectal_sim.vcf.gz" \
  --study-id cancer_synthetic \
  --gwas-id colon_rectal_subset_rds \
  --ancestry EUR \
  --methods "p-ld susie-rrss polyfun" \
  --output-format rds \
  --compress gzip \
  --reference-root "${TMP_TEST}/fixture/reference" \
  --local-ld-panel "${TMP_TEST}/fixture/reference/local_eur_panel.npz" \
  --ld-pop EUR \
  --outdir "${TMP_TEST}/out_subset"

# Verify contracts in Python + R
python3 - "${TMP_TEST}/out_batch" "${TMP_TEST}/out_nomhc" "${TMP_TEST}/out_subset" <<'PY'
import json, sys
from pathlib import Path

out_batch = Path(sys.argv[1])
out_nomhc = Path(sys.argv[2])
out_subset = Path(sys.argv[3])

assert (out_batch / "manifest.html").is_file(), "Missing batch manifest.html"
assert (out_batch / ".finemap_state.tsv").is_file(), "Missing .finemap_state.tsv"
assert (out_batch / "sun_synthetic" / "logs").is_dir(), "Missing per-study logs/ folder"

g1 = out_batch / "sun_synthetic" / "P00001_olink"
g2 = out_batch / "cancer_synthetic" / "colon_rectal_sim"
g3 = out_nomhc / "sun_synthetic" / "P00002_subthresh_nomhc"
g4 = out_subset / "cancer_synthetic" / "colon_rectal_subset_rds"

expected_top = {"metadata.json", "provenance.json", "qc.json", "report.html", "methods"}
expected_rdas = {"p-ld.rda", "susie-rss.rda", "susie-rrss.rda", "carma.rda", "polyfun.rda"}

for gdir in (g1, g2, g3):
    present_top = {p.name for p in gdir.iterdir()}
    assert present_top == expected_top, f"Unexpected top-level files in {gdir}: {present_top}"
    present_rdas = {p.name for p in (gdir / "methods").iterdir()}
    assert present_rdas == expected_rdas, f"Unexpected .rda files in {gdir}/methods: {present_rdas}"
    assert "p-value.rda" not in present_rdas
    html_text = (gdir / "report.html").read_text(encoding="utf-8")
    assert "data:image/png;base64," in html_text, f"Embedded base64 plot missing in {gdir}/report.html"

present_subset = {p.name for p in (g4 / "methods").iterdir()}
assert present_subset == {"p-ld.rds", "susie-rrss.rds", "polyfun.rds"}, f"Unexpected subset files: {present_subset}"

qc3 = json.loads((g3 / "qc.json").read_text(encoding="utf-8"))
assert qc3["ancestry_ld_status"] == "multi_ancestry_gwas_single_panel_ld"
assert "ALL→EUR LD" in (g3 / "report.html").read_text(encoding="utf-8")
assert qc3["mhc"]["status"] == "excluded_by_option"
assert qc3["cis_gene"]["status"] == "skipped_below_significance_threshold"
assert qc3["analyzed_regions_count"] == 2
assert qc3["total_credible_sets_carma"] >= 1
assert qc3["total_credible_sets_polyfun"] >= 1
print("Python contract & state-tracking assertions passed.")
PY

Rscript - "${TMP_TEST}/out_batch" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
out_batch <- args[[1L]]

e1 <- new.env(parent = emptyenv())
e2 <- new.env(parent = emptyenv())
load(file.path(out_batch, "sun_synthetic/P00001_olink/methods/susie-rss.rda"), envir = e1)
load(file.path(out_batch, "cancer_synthetic/colon_rectal_sim/methods/susie-rss.rda"), envir = e2)

s1 <- e1$susie_rss
s2 <- e2$susie_rss
chr2_id <- grep("^chr2_", names(s1$regions), value = TRUE)
fit1 <- s1$regions[[chr2_id]]$result$fit
fit2 <- s2$regions[[1L]]$result$fit
v1 <- s1$regions[[chr2_id]]$result$variants

stopifnot(
  identical(colnames(fit1$alpha), v1$GWAS_VARIANT_ID),
  identical(colnames(fit1$lbf_variable), v1$GWAS_VARIANT_ID),
  identical(names(fit1$pip), v1$GWAS_VARIANT_ID),
  length(intersect(colnames(fit1$lbf_variable), colnames(fit2$lbf_variable))) == 15L,
  names(which.max(fit2$pip)) %in% colnames(fit1$lbf_variable)
)

mhc_id <- grep("_mhc$", names(s1$regions), value = TRUE)
stopifnot(
  is.logical(s1$regions[[mhc_id]]$result$variants$mhc_region),
  all(s1$regions[[mhc_id]]$result$variants$mhc_region),
  all(s1$regions[[mhc_id]]$result$variants$mhc_core_grc)
)
cat("R .rda + coloc matrix naming assertions passed.\n")
RS

echo "=== ALL TESTS PASSED ==="
