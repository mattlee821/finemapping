#!/usr/bin/env bash
# End-to-end validation test suite for GWASfinemapping (mirrors GWASprocessing/src/test.sh).
# Runs: (1) Python unit tests, (2) R testthat suite, (3) Nextflow batch + --exclude-mhc + cis-gate smoke tests.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/env.sh
source "${SCRIPT_DIR}/lib/env.sh"
gwasfm_activate_runtime

cd "${GWASFM_ROOT}"

echo "=== [1/3] Python unit tests ==="
if python3 -c "import pytest" >/dev/null 2>&1; then
  python3 -m pytest -q tests/python
else
  python3 -m unittest discover -s tests/python -v
fi

echo "=== [2/3] R package testthat suite ==="
Rscript -e '
  for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
  res <- testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
'

echo "=== [3/3] End-to-end Nextflow manifest, MHC policy, cis gate & storage-light contract ==="
TMP_TEST="$(mktemp -d "${TMPDIR:-/tmp}/gwasfm-test-XXXXXX")"
trap 'rm -rf "${TMP_TEST}"' EXIT

python3 scripts/build_synthetic_fixture.py --out "${TMP_TEST}/fixture" >/dev/null

# 3a. Batch manifest run (protein GWAS with cis + merged trans + MHC, plus non-protein colon/rectal sim GWAS)
GWASFM_NXF_WORK="${TMP_TEST}/nxf-work" \
GWASFM_CACHE_DIR="${TMP_TEST}/cache" \
bash scripts/finemap.sh \
  --manifest "${TMP_TEST}/fixture/manifest.tsv" \
  --reference-root "${TMP_TEST}/fixture/reference" \
  --local-ld-panel "${TMP_TEST}/fixture/reference/local_eur_panel.npz" \
  --ld-pop EUR \
  --outdir "${TMP_TEST}/out_batch"

# 3b. Single-GWAS run with --exclude-mhc AND sub-threshold cis_gene (P00002 on chr3 has max LP ~ 2.7 < 7.30103)
GWASFM_NXF_WORK="${TMP_TEST}/nxf-work" \
GWASFM_CACHE_DIR="${TMP_TEST}/cache" \
bash scripts/finemap.sh \
  --vcf "${TMP_TEST}/fixture/gwas/sun_synthetic/P00001_olink/P00001_olink.vcf.gz" \
  --study-id sun_synthetic \
  --gwas-id P00002_subthresh_nomhc \
  --gene-id P00002 \
  --exclude-mhc \
  --reference-root "${TMP_TEST}/fixture/reference" \
  --local-ld-panel "${TMP_TEST}/fixture/reference/local_eur_panel.npz" \
  --ld-pop EUR \
  --outdir "${TMP_TEST}/out_nomhc"

# Verify storage-light contract, JSON provenance/QC, MHC flags, cis gate, report.html, manifest.html, and coloc-ready .rda
python3 - "${TMP_TEST}/out_batch" "${TMP_TEST}/out_nomhc" <<'PY'
import json
import sys
from pathlib import Path

out_batch = Path(sys.argv[1])
out_nomhc = Path(sys.argv[2])

assert (out_batch / "manifest.html").is_file(), "Missing batch manifest.html"
g1 = out_batch / "sun_synthetic" / "P00001_olink"
g2 = out_batch / "cancer_synthetic" / "colon_rectal_sim"
g3 = out_nomhc / "sun_synthetic" / "P00002_subthresh_nomhc"

expected_top = {"metadata.json", "provenance.json", "qc.json", "report.html", "methods"}
expected_rdas = {"p-value.rda", "p-ld.rda", "susie-rss.rda", "susie-rrss.rda"}

for gdir in (g1, g2, g3):
    present_top = {p.name for p in gdir.iterdir()}
    assert present_top == expected_top, f"Unexpected top-level files in {gdir}: {present_top}"
    present_rdas = {p.name for p in (gdir / "methods").iterdir()}
    assert present_rdas == expected_rdas, f"Unexpected .rda files in {gdir}/methods: {present_rdas}"
    html_text = (gdir / "report.html").read_text(encoding="utf-8")
    assert "data:image/png;base64," in html_text, f"Embedded base64 plot missing in {gdir}/report.html"

qc1 = json.loads((g1 / "qc.json").read_text(encoding="utf-8"))
assert qc1["analyzed_regions_count"] == 3, f"Expected 3 regions in P00001_olink, got {qc1['analyzed_regions_count']}"
assert qc1["cis_gene"]["status"] == "activated"
assert qc1["mhc"]["status"] == "included_single_region"
assert qc1["total_credible_sets_susie_rss"] >= 4, f"Expected >=4 CSs across 3 loci, got {qc1['total_credible_sets_susie_rss']}"

qc2 = json.loads((g2 / "qc.json").read_text(encoding="utf-8"))
assert qc2["cis_gene"]["status"] == "not_applicable", qc2["cis_gene"]
assert qc2["analyzed_regions_count"] == 1

qc3 = json.loads((g3 / "qc.json").read_text(encoding="utf-8"))
assert qc3["mhc"]["status"] == "excluded_by_option", qc3["mhc"]
assert qc3["cis_gene"]["status"] == "skipped_below_significance_threshold", qc3["cis_gene"]
assert qc3["analyzed_regions_count"] == 2  # chr1 and chr2 genome-wide loci only (MHC excluded, chr3 cis skipped)
print("Python contract assertions passed.")
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
stopifnot(length(s1$regions) == 3L, length(s2$regions) == 1L)

chr2_id <- grep("^chr2_", names(s1$regions), value = TRUE)
stopifnot(length(chr2_id) == 1L)
fit1 <- s1$regions[[chr2_id]]$result$fit
fit2 <- s2$regions[[1L]]$result$fit
v1 <- s1$regions[[chr2_id]]$result$variants

stopifnot(
  !is.null(colnames(fit1$alpha)),
  !is.null(colnames(fit1$lbf_variable)),
  identical(colnames(fit1$alpha), v1$GWAS_VARIANT_ID),
  identical(names(fit1$pip), v1$GWAS_VARIANT_ID),
  all(c("AF", "AFS", "HM", "LD_R2_TO_LEAD", "EXPECTED_Z", "mhc_region", "mhc_core_grc") %in% names(v1)),
  all(!is.na(v1$AF)),
  identical(unique(v1$HM), "direct"),
  length(intersect(colnames(fit1$lbf_variable), colnames(fit2$lbf_variable))) == 15L,
  names(which.max(fit2$pip)) %in% colnames(fit1$lbf_variable)
)

mhc_id <- grep("_mhc$", names(s1$regions), value = TRUE)
stopifnot(
  length(mhc_id) == 1L,
  all(s1$regions[[mhc_id]]$result$variants$mhc_region),
  all(s1$regions[[mhc_id]]$result$variants$mhc_core_grc)
)
cat("R .rda + coloc matrix naming assertions passed.\n")
RS

echo "=== ALL TESTS PASSED ==="
