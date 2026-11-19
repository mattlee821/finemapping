#!/usr/bin/env Rscript
# ==============================================================================
# Single-Load Multi-Method Runner (`p-ld`, `carma`, `susie-rss`, `susie-rrss`, `polyfun`)
# ==============================================================================
# Loads `matched_variants.tsv` and `matrix.npy` ONCE into R memory per locus,
# optionally unlinks the multi-GB `matrix.npy` immediately after read to free
# scratch/tmp disk space, runs unscaled-R methods (`p_ld`, `carma`), scales R
# in-place in C (`0` bytes copied) to `XtX = (n_eff - 1) * R`, and runs scaled
# SuSiE methods (`susie_rss`, `susie_rrss`, `polyfun`) back-to-back.
# ==============================================================================

t_total_start <- proc.time()

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
repo_root <- if (length(script_arg) > 0L) {
  dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))))
} else {
  getwd()
}

variants_path         <- NULL
matrix_path           <- NULL
matrix_ids            <- NULL
region_id             <- NULL
region_type           <- "lead_locus"
methods_csv           <- "p_ld,susie_rss,susie_rrss,carma,polyfun"
out_methods_dir       <- NULL
log_dir               <- NULL
safe_region           <- NULL
effective_n           <- NULL
reference_sample_size <- 420542L
max_causal            <- 10L
max_iter              <- 50L
coverage              <- 0.95
min_abs_corr          <- 0.5
r_mismatch            <- "eb_mix"
polyfun_priors        <- NULL
missing_prior_policy  <- "baseline"
rds_compress          <- "gzip"
unlink_matrix         <- FALSE
allow_skip            <- FALSE

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) %% 2L != 0L) {
  stop("CLI options must be --key value pairs.", call. = FALSE)
}
for (i in seq.int(1L, length(cli_args), by = 2L)) {
  flag <- cli_args[[i]]
  val  <- cli_args[[i + 1L]]
  if (flag == "--variants")             variants_path <- val
  else if (flag == "--matrix")          matrix_path <- val
  else if (flag == "--matrix-ids")      matrix_ids <- val
  else if (flag == "--region-id")       region_id <- val
  else if (flag == "--region-type")     region_type <- val
  else if (flag == "--methods")         methods_csv <- val
  else if (flag == "--out-methods-dir") out_methods_dir <- val
  else if (flag == "--log-dir")         log_dir <- val
  else if (flag == "--safe-region")     safe_region <- val
  else if (flag == "--effective-n")     effective_n <- as.numeric(val)
  else if (flag == "--reference-sample-size") reference_sample_size <- as.numeric(val)
  else if (flag == "--max-causal")      max_causal <- as.integer(val)
  else if (flag == "--max-iter")        max_iter <- as.integer(val)
  else if (flag == "--coverage")        coverage <- as.numeric(val)
  else if (flag == "--min-abs-corr")    min_abs_corr <- as.numeric(val)
  else if (flag == "--r-mismatch")      r_mismatch <- val
  else if (flag == "--polyfun-priors")  polyfun_priors <- val
  else if (flag == "--missing-prior-policy") missing_prior_policy <- val
  else if (flag == "--compress") {
    rds_compress <- if (tolower(val) %in% c("false", "0", "none")) FALSE else val
  } else if (flag == "--unlink-matrix") {
    unlink_matrix <- tolower(val) %in% c("true", "1", "yes")
  } else if (flag == "--allow-skipped-locus") {
    allow_skip <- tolower(val) %in% c("true", "1", "yes")
  } else {
    stop("Unknown CLI option: ", flag, call. = FALSE)
  }
}

if (allow_skip && (!file.exists(variants_path) || !file.exists(matrix_path))) {
  quit(save = "no", status = 0L)
}

methods_req <- unique(Filter(nzchar, trimws(strsplit(methods_csv, ",")[[1L]])))
if (length(methods_req) == 0L) {
  quit(save = "no", status = 0L)
}
if (is.null(safe_region) || !nzchar(safe_region)) {
  safe_region <- gsub("[^A-Za-z0-9._-]+", "_", region_id)
}

# Load compiled gwasfinemapping C routines + latest R helpers
if (requireNamespace("gwasfinemapping", quietly = TRUE)) {
  library(gwasfinemapping)
}
# Source repo R/ files if running against an uninstalled or pre-update build
if (!exists("scale_ld_matrix_inplace", mode = "function") &&
    !(requireNamespace("gwasfinemapping", quietly = TRUE) &&
      exists("scale_ld_matrix_inplace", where = asNamespace("gwasfinemapping"), mode = "function"))) {
  for (f in list.files(file.path(repo_root, "R"), pattern = "[.]R$", full.names = TRUE)) {
    source(f)
  }
}

write_step_log <- function(method_id, msg) {
  if (!is.null(log_dir) && nzchar(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_p <- file.path(log_dir, sprintf("03_method_%s_%s.log", safe_region, method_id))
    tryCatch(writeLines(msg, con = log_p), error = function(e) NULL)
  }
}

peak_heap_mb_now <- function() {
  gc_mat <- gc(FALSE)
  round(sum(gc_mat[, 6L] * c(56, 8)) / (1024^2), 1)
}

# --- 1. Read variants & pre-resolve PolyFun Parquet priors while heap is small ---
t_io_start <- proc.time()
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

priors_pre <- NULL
if ("polyfun" %in% methods_req) {
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
}

# --- 2. Load matrix.npy ONCE (unscaled if p_ld/carma needed, else pre-scaled XtX) ---
n_eff_load <- if (!is.null(effective_n)) {
  as.numeric(effective_n)
} else {
  min(variants$SS[!is.na(variants$SS)])
}
need_unscaled <- any(c("p_ld", "carma") %in% methods_req)
need_scaled <- any(c("susie_rss", "susie_rrss", "polyfun") %in% methods_req)
only_polyfun_skip <- identical(methods_req, "polyfun") &&
  !is.null(priors_pre) &&
  identical(priors_pre$status, "skipped_unsupported_chromosome")

initial_scale <- if (need_unscaled) 1 else (n_eff_load - 1)
R <- if (only_polyfun_skip) {
  matrix(1, nrow = 1L, ncol = 1L, dimnames = list("skip", "skip"))
} else {
  read_ld_matrix(matrix_path, ids_path = matrix_ids, scale = initial_scale)
}

if (!"LD_R2_TO_LEAD" %in% names(variants) && nrow(R) == nrow(variants)) {
  lead_idx <- which.max(ifelse(is.na(variants$LP), -Inf, as.numeric(variants$LP)))
  if (length(lead_idx) == 1L) {
    variants$LD_R2_TO_LEAD <- as.numeric((R[, lead_idx] / initial_scale)^2)
  }
}

if (unlink_matrix && file.exists(matrix_path)) {
  unlink(matrix_path)
}
shared_io_seconds <- (proc.time() - t_io_start)[["elapsed"]]
io_recorded <- FALSE

take_io_sec <- function() {
  if (!io_recorded) {
    io_recorded <<- TRUE
    shared_io_seconds
  } else {
    0
  }
}

save_bayesian_outputs <- function(method_id, result, io_sec, fit_sec, proc_start) {
  out_dir <- file.path(out_methods_dir, method_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_rds <- file.path(out_dir, "region_result.rds")
  qc_json <- file.path(out_dir, "method_qc.json")

  t_sv <- proc.time()
  saveRDS(
    list(region_id = region_id, region_type = region_type, result = result),
    file = out_rds,
    compress = rds_compress
  )
  save_sec <- (proc.time() - t_sv)[["elapsed"]]
  p_delta <- proc.time() - proc_start

  diag_compact <- result$diagnostics
  diag_compact$R_finite_diagnostics <- NULL
  max_pip_val <- if (all(is.na(result$variants$PIP))) {
    NA_real_
  } else {
    max(result$variants$PIP, na.rm = TRUE)
  }

  qc_payload <- c(
    list(
      method_id = method_id,
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
        io_seconds = round(io_sec, 3),
        fit_seconds = round(fit_sec, 3),
        save_seconds = round(save_sec, 3),
        wall_seconds = round(p_delta[["elapsed"]] + io_sec, 3),
        user_seconds = round(p_delta[["user.self"]] + p_delta[["user.child"]], 3),
        sys_seconds = round(p_delta[["sys.self"]] + p_delta[["sys.child"]], 3),
        peak_r_heap_mb = peak_heap_mb_now()
      )
    )
  )
  jsonlite::write_json(qc_payload, qc_json, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  write_step_log(
    method_id,
    sprintf("[%s] %s: %d variants, %d CS, fit=%.2fs, io=%.2fs, save=%.2fs",
            method_id, region_id, nrow(result$variants), nrow(result$credible_sets),
            fit_sec, io_sec, save_sec)
  )
}

# --- 3. Run unscaled-R methods (`p_ld`, `carma`) ------------------------------
if ("p_ld" %in% methods_req) {
  t_m0 <- proc.time()
  io_s <- take_io_sec()
  r2_thresholds <- c(0.1, 0.01, 0.001)
  t_f0 <- proc.time()
  pld_df <- prioritize_pld(variants, R, r2_thresholds)
  fit_s <- (proc.time() - t_f0)[["elapsed"]]

  selected_by_r2 <- stats::setNames(
    lapply(r2_thresholds, function(t) as.integer(sum(pld_df$selected[pld_df$r2_threshold == t]))),
    as.character(r2_thresholds)
  )
  out_dir <- file.path(out_methods_dir, "p_ld")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  t_sv <- proc.time()
  saveRDS(
    list(region_id = region_id, region_type = region_type, result = pld_df),
    file = file.path(out_dir, "region_result.rds"),
    compress = rds_compress
  )
  save_s <- (proc.time() - t_sv)[["elapsed"]]
  p_delta <- proc.time() - t_m0

  qc_pld <- list(
    method_id = "p_ld",
    region_id = region_id,
    region_type = region_type,
    variants = nrow(variants),
    selected_by_r2 = selected_by_r2,
    timing = list(
      io_seconds = round(io_s, 3),
      fit_seconds = round(fit_s, 3),
      save_seconds = round(save_s, 3),
      wall_seconds = round(p_delta[["elapsed"]] + io_s, 3),
      user_seconds = round(p_delta[["user.self"]] + p_delta[["user.child"]], 3),
      sys_seconds = round(p_delta[["sys.self"]] + p_delta[["sys.child"]], 3),
      peak_r_heap_mb = peak_heap_mb_now()
    )
  )
  jsonlite::write_json(qc_pld, file.path(out_dir, "method_qc.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  write_step_log(
    "p_ld",
    sprintf("[p_ld] %s: %d variants, fit=%.2fs, io=%.2fs", region_id, nrow(variants), fit_s, io_s)
  )
  rm(pld_df, qc_pld)
  invisible(gc(FALSE))
}

if ("carma" %in% methods_req) {
  t_m0 <- proc.time()
  io_s <- take_io_sec()
  t_f0 <- proc.time()
  res_carma <- run_carma(
    variants = variants,
    R = R,
    lambda = 1.0,
    effect_size_prior = "Spike-slab",
    tau = 0.04,
    L = max_causal,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    outlier_switch = TRUE,
    outlier_bf_index = 1 / 3.2,
    seed = 1L
  )
  fit_s <- (proc.time() - t_f0)[["elapsed"]]
  save_bayesian_outputs("carma", res_carma, io_s, fit_s, t_m0)
  rm(res_carma)
  invisible(gc(FALSE))
}

# --- 4. Scale R -> XtX in-place in C (0-copy) & run SuSiE methods -------------
if (need_unscaled && need_scaled && !only_polyfun_skip) {
  scale_fn <- if (exists("scale_ld_matrix_inplace", mode = "function")) {
    scale_ld_matrix_inplace
  } else {
    getFromNamespace("scale_ld_matrix_inplace", "gwasfinemapping")
  }
  R <- scale_fn(R, n_eff_load - 1)
}

if ("susie_rss" %in% methods_req) {
  t_m0 <- proc.time()
  io_s <- take_io_sec()
  t_f0 <- proc.time()
  res_rss <- run_susie_rss(
    variants = variants,
    R = R,
    effective_n = effective_n,
    L = max_causal,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter
  )
  fit_s <- (proc.time() - t_f0)[["elapsed"]]
  save_bayesian_outputs("susie_rss", res_rss, io_s, fit_s, t_m0)
  rm(res_rss)
  invisible(gc(FALSE))
}

if ("susie_rrss" %in% methods_req) {
  t_m0 <- proc.time()
  io_s <- take_io_sec()
  t_f0 <- proc.time()
  res_rrss <- run_susie_rrss(
    variants = variants,
    R = R,
    reference_sample_size = reference_sample_size,
    effective_n = effective_n,
    L = max_causal,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter,
    R_mismatch = r_mismatch
  )
  fit_s <- (proc.time() - t_f0)[["elapsed"]]
  save_bayesian_outputs("susie_rrss", res_rrss, io_s, fit_s, t_m0)
  rm(res_rrss)
  invisible(gc(FALSE))
}

if ("polyfun" %in% methods_req) {
  t_m0 <- proc.time()
  io_s <- take_io_sec()
  t_f0 <- proc.time()
  if (!identical(priors_pre$status, "skipped_unsupported_chromosome") &&
      all(priors_pre$keep_mask)) {
    variants_pf <- variants
    variants_pf$POLYFUN_SNPVAR <- priors_pre$snpvar
    res_pf <- run_polyfun_susie(
      variants = variants_pf,
      R = R,
      prior_weights = priors_pre$prior_weights,
      effective_n = effective_n,
      L = max_causal,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      max_iter = max_iter
    )
    res_pf$variants$SNPVAR <- priors_pre$snpvar
    res_pf$variants$PRIOR_MATCHED <- priors_pre$matched
    res_pf$diagnostics$prior_source <- priors_pre$source
    res_pf$diagnostics$prior_matched_variants <- as.integer(sum(priors_pre$matched))
    res_pf$diagnostics$prior_baseline_assigned_variants <- as.integer(sum(!priors_pre$matched))
    res_pf$diagnostics$prior_match_rate <- round(mean(priors_pre$matched), 4)
    res_pf$provenance$parameters$prior_source <- priors_pre$source
    res_pf$provenance$parameters$missing_prior_policy <- missing_prior_policy
    res_pf$provenance$parameters$prior_floor_value <- priors_pre$floor_value
    rm(variants_pf)
  } else {
    res_pf <- run_polyfun_susie(
      variants = variants,
      R = R,
      polyfun_priors = polyfun_priors,
      missing_prior_policy = missing_prior_policy,
      effective_n = effective_n,
      L = max_causal,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      max_iter = max_iter
    )
  }
  fit_s <- (proc.time() - t_f0)[["elapsed"]]
  save_bayesian_outputs("polyfun", res_pf, io_s, fit_s, t_m0)
  rm(res_pf)
  invisible(gc(FALSE))
}
