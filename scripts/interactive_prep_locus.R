#!/usr/bin/env Rscript
# ==============================================================================
# Interactive Walkthrough: Stage 1 (Screen GWAS) & Stage 2 (Prepare Locus)
# ==============================================================================
# Step through this script line-by-line (Cmd+Enter) in RStudio / Positron / R
# to run and inspect the P <= 5e-8 screen and shared locus preparation on NPPB
# before opening any standalone method script in scripts/methods/*.R.
# ==============================================================================

# --- 1. Paths & configuration -------------------------------------------------
repo_root   <- normalizePath(getwd(), mustWork = TRUE)
home_dir    <- Sys.getenv("HOME")
nppb_dir    <- file.path(home_dir, "research-local/GWASprocessing/data/GWAS/sun_2023_37794186/standardise/P16860_OID20049")
work_dir    <- file.path(home_dir, "research-local/GWASfinemapping/test_nppb")
local_bm    <- file.path(home_dir, "research-local/GWASfinemapping/nppb_test_bm/UKBB.EUR.ldadj.bm")
python_bin  <- file.path(home_dir, ".cache/gwasfinemapping/env/bin/python")
if (!file.exists(python_bin)) python_bin <- Sys.which("python3")

# --- 2. Stage 1: Genome-wide P <= 5e-8 screen & region construction -----------
screen_json <- file.path(work_dir, "screen_summary.json")
status_screen <- system2(
  python_bin,
  args = c(
    file.path(repo_root, "scripts/screen_gwas.py"),
    "--vcf",      file.path(nppb_dir, "P16860_OID20049.gwas.vcf.gz"),
    "--metadata", file.path(nppb_dir, "P16860_OID20049.metadata.json"),
    "--qc",       file.path(nppb_dir, "P16860_OID20049.qc.json"),
    "--study-id", "sun_2023_37794186",
    "--gwas-id",  "P16860_OID20049",
    "--max-loci", "2",
    "--out",      screen_json
  )
)
stopifnot(status_screen == 0L)

screen <- jsonlite::fromJSON(screen_json, simplifyVector = TRUE)
cat("Candidate regions found:", screen$candidate_regions_total, "\n")
print(screen$selected_regions[, c("region_id", "region_type", "CHROM", "start", "end", "significant_variant_count")])

# --- 3. Stage 2: Prepare shared locus inputs (liftover + Pan-UKB LD) ----------
# Uses local NPPB test BlockMatrix if present; omit --bm-path to fetch live from AWS S3.
bm_args <- if (dir.exists(local_bm)) c("--bm-path", local_bm) else character()

for (idx in c(0L, 1L)) {
  st <- system2(
    python_bin,
    args = c(
      file.path(repo_root, "scripts/prep_locus.py"),
      "--screen-json", screen_json,
      "--locus-index", as.character(idx),
      "--outdir", file.path(work_dir, "loci", screen$selected_regions$region_id[[idx + 1L]]),
      bm_args
    )
  )
  stopifnot(st == 0L)
}

# --- 4. Inspect prepared cis_NPPB inputs directly in R ------------------------
if (requireNamespace("gwasfinemapping", quietly = TRUE)) library(gwasfinemapping)
cis_dir <- file.path(work_dir, "loci", screen$selected_regions$region_id[[1L]])
prep_qc <- jsonlite::fromJSON(file.path(cis_dir, "locus_prep_qc.json"))
variants <- utils::read.delim(file.path(cis_dir, "matched_variants.tsv"), stringsAsFactors = FALSE)
R <- read_ld_matrix(file.path(cis_dir, "matrix.npy"))

cat(sprintf(
  "Prepared %s: %d variants, R dim = %d x %d, lead = %s (LP = %.2f)\n",
  prep_qc$region_id, nrow(variants), nrow(R), ncol(R),
  prep_qc$lead_variant_id, prep_qc$lead_lp
))
