#!/usr/bin/env Rscript
# ==============================================================================
# Standalone Method Script: SuSiE-rRSS (Robust Finite-Sample / Mismatch SuSiE)
# ==============================================================================
# Designed for BOTH:
#   1. Interactive line-by-line execution in RStudio / Positron / R console
#      (uses the default NPPB test locus paths in Section 1 when no CLI flags
#      are passed), AND
#   2. Automated CLI / Nextflow execution:
#      Rscript scripts/methods/run_susie_rrss.R \
#        --variants /path/to/matched_variants.tsv \
#        --matrix   /path/to/matrix.npy \
#        --region-id cis_NPPB_chr1_10857464_12858945 \
#        --region-type cis_gene \
#        --reference-sample-size 420542 \
#        --out-rds  /path/to/susie_rrss/region_result.rds \
#        --qc-json  /path/to/susie_rrss/method_qc.json
# ==============================================================================

t_script_start <- proc.time()
invisible(gc(reset = TRUE))

# --- 1. Configuration (editable interactive defaults + CLI overrides) ---------
default_locus_dir <- file.path(
  Sys.getenv("HOME"),
  "research-local/GWASfinemapping/test_nppb/loci/cis_NPPB_chr1_10857464_12858945"
)

variants_path         <- file.path(default_locus_dir, "matched_variants.tsv")
matrix_path           <- file.path(default_locus_dir, "matrix.npy")
matrix_ids            <- NULL
region_id             <- "cis_NPPB_chr1_10857464_12858945"
region_type           <- "cis_gene"
reference_sample_size <- 420542
effective_n           <- NULL
max_causal            <- 10L
max_iter              <- 100L
susie_refine          <- FALSE
coverage              <- 0.95
min_abs_corr          <- 0.5
r_mismatch            <- "eb_mix"
out_rds               <- file.path(default_locus_dir, "methods/susie_rrss/region_result.rds")
qc_json               <- file.path(default_locus_dir, "methods/susie_rrss/method_qc.json")
rds_compress          <- FALSE
allow_skip            <- FALSE

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) > 0L) {
  if (length(cli_args) %% 2L != 0L) {
    stop("CLI options must be --key value pairs.", call. = FALSE)
  }
  for (i in seq.int(1L, length(cli_args), by = 2L)) {
    flag <- cli_args[[i]]
    val  <- cli_args[[i + 1L]]
    if (flag == "--variants")        variants_path <- val
    else if (flag == "--matrix")     matrix_path   <- val
    else if (flag == "--matrix-ids") matrix_ids    <- val
    else if (flag == "--region-id")  region_id     <- val
    else if (flag == "--region-type") region_type  <- val
    else if (flag == "--reference-sample-size") reference_sample_size <- as.numeric(val)
    else if (flag == "--effective-n") effective_n  <- as.numeric(val)
    else if (flag == "--max-causal") max_causal    <- as.integer(val)
    else if (flag == "--max-iter")   max_iter      <- as.integer(val)
    else if (flag == "--susie-refine") susie_refine <- tolower(val) %in% c("true", "1", "yes")
    else if (flag == "--coverage")   coverage      <- as.numeric(val)
    else if (flag == "--min-abs-corr") min_abs_corr <- as.numeric(val)
    else if (flag == "--r-mismatch") r_mismatch    <- val
    else if (flag == "--out-rds")    out_rds       <- val
    else if (flag == "--qc-json")    qc_json       <- val
    else if (flag == "--allow-skipped-locus") {
      allow_skip <- tolower(val) %in% c("true", "1", "yes")
    }
    else if (flag == "--compress") {
      rds_compress <- if (tolower(val) %in% c("false", "0", "none")) FALSE else val
    } else {
      stop("Unknown CLI option: ", flag, call. = FALSE)
    }
  }
}

if (allow_skip && (!file.exists(variants_path) || !file.exists(matrix_path))) {
  quit(save = "no", status = 0L)
}

# --- 2. Load gwasfinemapping functions ----------------------------------------
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
repo_root <- if (length(script_arg) > 0L) {
  dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))))
} else {
  getwd()
}
if (requireNamespace("gwasfinemapping", quietly = TRUE)) {
  library(gwasfinemapping)
}
if (dir.exists(file.path(repo_root, "R"))) {
  for (f in sort(list.files(file.path(repo_root, "R"), pattern = "[.]R$", full.names = TRUE))) {
    source(f, local = FALSE)
  }
}

# --- 3. Read & normalise locus variant table ----------------------------------
t_io_start <- proc.time()
if (!file.exists(variants_path)) {
  stop("Matched variants table not found: ", variants_path, call. = FALSE)
}
variants <- utils::read.delim(
  variants_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
for (bcol in intersect(c("mhc_region", "mhc_core_grc", "PANEL_ALLELE_SWAPPED"), names(variants))) {
  if (!is.logical(variants[[bcol]])) {
    variants[[bcol]] <- tolower(as.character(variants[[bcol]])) == "true"
  }
}

# --- 4. Read Fortran-order binary LD matrix (.npy) scaled in-place to XtX -----
# Pre-scaling R -> (n - 1) * R on load avoids allocating extra M x M copies for
# XtX and sum(R^2), cutting peak RAM by >65% with exact bit-for-bit PIPs.
ss_non_na <- if ("SS" %in% names(variants)) as.numeric(variants$SS[!is.na(variants$SS)]) else numeric()
n_eff_load <- if (!is.null(effective_n)) {
  as.numeric(effective_n)
} else if (length(ss_non_na) > 0L) {
  min(ss_non_na)
} else {
  NULL
}
R <- read_ld_matrix(matrix_path, ids_path = matrix_ids, scale = if (is.null(n_eff_load)) 1 else (n_eff_load - 1))
io_seconds <- (proc.time() - t_io_start)[["elapsed"]]

# --- 5. Fit SuSiE-rRSS (R_finite + R_mismatch) & extract 95% credible sets ----
t_fit_start <- proc.time()
result <- run_susie_rrss(
  variants = variants,
  R = R,
  reference_sample_size = reference_sample_size,
  effective_n = effective_n,
  L = max_causal,
  coverage = coverage,
  min_abs_corr = min_abs_corr,
  max_iter = max_iter,
  refine = susie_refine,
  R_mismatch = r_mismatch,
  allow_missing_n = TRUE
)
fit_seconds <- (proc.time() - t_fit_start)[["elapsed"]]

# --- 6. Save standalone region_result.rds & method_qc.json --------------------
t_save_start <- proc.time()
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    region_id = region_id,
    region_type = region_type,
    result = result
  ),
  file = out_rds,
  compress = rds_compress
)
save_seconds <- (proc.time() - t_save_start)[["elapsed"]]

proc_delta <- proc.time() - t_script_start
gc_mat <- gc(FALSE)
peak_ram_mb <- round(sum(gc_mat[, 6L] * c(56, 8)) / (1024^2), 1)

diag_compact <- result$diagnostics
diag_compact$R_finite_diagnostics <- NULL

qc_payload <- c(
  list(
    method_id = "susie_rrss",
    region_id = region_id,
    region_type = region_type,
    variants = nrow(result$variants),
    credible_sets = nrow(result$credible_sets),
    max_pip = max(result$variants$PIP, na.rm = TRUE),
    min_cs_purity = if (nrow(result$credible_sets) > 0L) {
      min(result$credible_sets$min_abs_corr, na.rm = TRUE)
    } else {
      NA_real_
    }
  ),
  diag_compact,
  list(
    timing = list(
      io_seconds = round(io_seconds, 3),
      fit_seconds = round(fit_seconds, 3),
      save_seconds = round(save_seconds, 3),
      wall_seconds = round(proc_delta[["elapsed"]], 3),
      user_seconds = round(proc_delta[["user.self"]] + proc_delta[["user.child"]], 3),
      sys_seconds = round(proc_delta[["sys.self"]] + proc_delta[["sys.child"]], 3),
      peak_r_heap_mb = peak_ram_mb
    )
  )
)

if (!is.null(qc_json) && nzchar(qc_json)) {
  dir.create(dirname(qc_json), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(qc_payload, qc_json, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
}

if (interactive()) {
  message(sprintf(
    "[susie-rrss] %s (%d variants): %d CS, converged=%s, R_rel_flag=%s | fit=%.2fs, total=%.2fs",
    region_id, nrow(result$variants), nrow(result$credible_sets),
    as.character(result$diagnostics$converged),
    as.character(result$diagnostics$R_reliability_flag),
    fit_seconds, proc_delta[["elapsed"]]
  ))
}
