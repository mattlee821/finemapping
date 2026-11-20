#!/usr/bin/env Rscript
# ==============================================================================
# Standalone Method Script: PolyFun + SuSiE (Functionally Informed Fine-Mapping)
# ==============================================================================
# Designed for BOTH:
#   1. Interactive line-by-line execution in RStudio / Positron / R console
#      (uses the default NPPB test locus paths in Section 1 when no CLI flags
#      are passed), AND
#   2. Automated CLI / Nextflow execution:
#      Rscript scripts/methods/run_polyfun.R \
#        --variants /path/to/matched_variants.tsv \
#        --matrix   /path/to/matrix.npy \
#        --region-id cis_NPPB_chr1_10857464_12858945 \
#        --region-type cis_gene \
#        --out-rds  /path/to/polyfun/region_result.rds \
#        --qc-json  /path/to/polyfun/method_qc.json
# ==============================================================================

t_script_start <- proc.time()
invisible(gc(reset = TRUE))

# --- 1. Configuration (editable interactive defaults + CLI overrides) ---------
default_locus_dir <- file.path(
  Sys.getenv("HOME"),
  "research-local/GWASfinemapping/test_nppb/loci/cis_NPPB_chr1_10857464_12858945"
)

variants_path        <- file.path(default_locus_dir, "matched_variants.tsv")
matrix_path          <- file.path(default_locus_dir, "matrix.npy")
matrix_ids           <- NULL
region_id            <- "cis_NPPB_chr1_10857464_12858945"
region_type          <- "cis_gene"
polyfun_priors       <- NULL
missing_prior_policy <- "baseline_floor"
effective_n          <- NULL
max_causal           <- 10L
max_iter             <- 100L
susie_refine         <- FALSE
coverage             <- 0.95
min_abs_corr         <- 0.5
out_rds              <- file.path(default_locus_dir, "methods/polyfun/region_result.rds")
qc_json              <- file.path(default_locus_dir, "methods/polyfun/method_qc.json")
rds_compress         <- FALSE
allow_skip           <- FALSE

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) > 0L) {
  if (length(cli_args) %% 2L != 0L) {
    stop("CLI options must be --key value pairs.", call. = FALSE)
  }
  for (i in seq.int(1L, length(cli_args), by = 2L)) {
    flag <- cli_args[[i]]
    val  <- cli_args[[i + 1L]]
    if (flag == "--variants")             variants_path <- val
    else if (flag == "--matrix")          matrix_path   <- val
    else if (flag == "--matrix-ids")      matrix_ids    <- val
    else if (flag == "--region-id")       region_id     <- val
    else if (flag == "--region-type")     region_type   <- val
    else if (flag == "--polyfun-priors")  polyfun_priors <- val
    else if (flag == "--missing-prior-policy") missing_prior_policy <- val
    else if (flag == "--effective-n")     effective_n   <- as.numeric(val)
    else if (flag == "--max-causal")      max_causal    <- as.integer(val)
    else if (flag == "--max-iter")        max_iter      <- as.integer(val)
    else if (flag == "--susie-refine")    susie_refine  <- tolower(val) %in% c("true", "1", "yes")
    else if (flag == "--coverage")        coverage      <- as.numeric(val)
    else if (flag == "--min-abs-corr")    min_abs_corr  <- as.numeric(val)
    else if (flag == "--reference-sample-size") {
      # Accepted and ignored so Nextflow can pass a uniform Bayesian interface
    }
    else if (flag == "--out-rds")         out_rds       <- val
    else if (flag == "--qc-json")         qc_json       <- val
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

# --- 4. Resolve PolyFun functional priors BEFORE loading large LD matrix ------
# Resolving Parquet priors while R's heap is ~40 MB (rather than after loading a
# multi-GB LD matrix) avoids fork() copy-on-write page pressure on 8 GB hosts.
resolve_pf_fn <- if (exists("resolve_polyfun_priors", mode = "function")) {
  resolve_polyfun_priors
} else {
  getFromNamespace("resolve_polyfun_priors", "gwasfinemapping")
}
priors_pre <- resolve_pf_fn(
  variants = variants,
  polyfun_priors = polyfun_priors,
  missing_prior_policy = missing_prior_policy
)

# --- 5. Read Fortran-order binary LD matrix (.npy) scaled in-place to XtX -----
ss_non_na <- if ("SS" %in% names(variants)) as.numeric(variants$SS[!is.na(variants$SS)]) else numeric()
n_eff_load <- if (!is.null(effective_n)) {
  as.numeric(effective_n)
} else if (length(ss_non_na) > 0L) {
  min(ss_non_na)
} else {
  NULL
}
R <- if (identical(priors_pre$status, "skipped_unsupported_chromosome")) {
  matrix(1, nrow = 1L, ncol = 1L, dimnames = list("skip", "skip"))
} else {
  read_ld_matrix(matrix_path, ids_path = matrix_ids, scale = if (is.null(n_eff_load)) 1 else (n_eff_load - 1))
}
io_seconds <- (proc.time() - t_io_start)[["elapsed"]]

# --- 6. Fit functionally informed PolyFun + SuSiE -----------------------------
t_fit_start <- proc.time()
if (!identical(priors_pre$status, "skipped_unsupported_chromosome") &&
    all(priors_pre$keep_mask)) {
  variants$POLYFUN_SNPVAR <- priors_pre$snpvar
  result <- run_polyfun_susie(
    variants = variants,
    R = R,
    prior_weights = priors_pre$prior_weights,
    effective_n = effective_n,
    L = max_causal,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter,
    refine = susie_refine,
    allow_missing_n = TRUE
  )
  result$variants$SNPVAR <- priors_pre$snpvar
  result$variants$PRIOR_MATCHED <- priors_pre$matched
  result$diagnostics$prior_source <- priors_pre$source
  result$diagnostics$prior_matched_variants <- as.integer(sum(priors_pre$matched))
  result$diagnostics$prior_baseline_assigned_variants <- as.integer(sum(!priors_pre$matched))
  result$diagnostics$prior_match_rate <- round(mean(priors_pre$matched), 4)
  result$provenance$parameters$prior_source <- priors_pre$source
  result$provenance$parameters$missing_prior_policy <- missing_prior_policy
  result$provenance$parameters$prior_floor_value <- priors_pre$floor_value
} else {
  result <- run_polyfun_susie(
    variants = variants,
    R = R,
    polyfun_priors = polyfun_priors,
    missing_prior_policy = missing_prior_policy,
    effective_n = effective_n,
    L = max_causal,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter,
    refine = susie_refine,
    allow_missing_n = TRUE
  )
}
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

max_pip_val <- if (all(is.na(result$variants$PIP))) {
  NA_real_
} else {
  max(result$variants$PIP, na.rm = TRUE)
}

qc_payload <- c(
  list(
    method_id = "polyfun",
    region_id = region_id,
    region_type = region_type,
    variants = nrow(result$variants),
    credible_sets = nrow(result$credible_sets),
    max_pip = max_pip_val,
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
    "[polyfun] %s (%d variants, %.1f%% prior matched): %d CS, converged=%s, max_pip=%.4f | fit=%.2fs, total=%.2fs",
    region_id, nrow(result$variants), 100 * (result$diagnostics$prior_match_rate %||% 0),
    nrow(result$credible_sets), as.character(result$diagnostics$converged),
    max_pip_val, fit_seconds, proc_delta[["elapsed"]]
  ))
}
