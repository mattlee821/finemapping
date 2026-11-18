#!/usr/bin/env Rscript
# ==============================================================================
# Standalone Method Script: p-LD Prioritisation (r^2 = 0.1, 0.01, 0.001)
# ==============================================================================
# Designed for BOTH:
#   1. Interactive line-by-line execution in RStudio / Positron / R console
#      (uses the default NPPB test locus paths in Section 1 when no CLI flags
#      are passed), AND
#   2. Automated CLI / Nextflow execution:
#      Rscript scripts/methods/run_p_ld.R \
#        --variants /path/to/matched_variants.tsv \
#        --matrix   /path/to/matrix.npy \
#        --region-id cis_NPPB_chr1_10857464_12858945 \
#        --region-type cis_gene \
#        --out-rds  /path/to/p_ld/region_result.rds \
#        --qc-json  /path/to/p_ld/method_qc.json
# ==============================================================================

t_script_start <- proc.time()
invisible(gc(reset = TRUE))

# --- 1. Configuration (editable interactive defaults + CLI overrides) ---------
default_locus_dir <- file.path(
  Sys.getenv("HOME"),
  "research-local/GWASfinemapping/test_nppb/loci/cis_NPPB_chr1_10857464_12858945"
)

variants_path <- file.path(default_locus_dir, "matched_variants.tsv")
matrix_path   <- file.path(default_locus_dir, "matrix.npy")
matrix_ids    <- NULL
region_id     <- "cis_NPPB_chr1_10857464_12858945"
region_type   <- "cis_gene"
r2_thresholds <- c(0.1, 0.01, 0.001)
out_rds       <- file.path(default_locus_dir, "methods/p_ld/region_result.rds")
qc_json       <- file.path(default_locus_dir, "methods/p_ld/method_qc.json")
rds_compress  <- FALSE
allow_skip    <- FALSE

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
    else if (flag == "--r2-thresholds") {
      r2_thresholds <- as.numeric(strsplit(val, ",", fixed = TRUE)[[1L]])
    }
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
if (requireNamespace("gwasfinemapping", quietly = TRUE)) {
  library(gwasfinemapping)
} else {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  repo_root <- if (length(script_arg) > 0L) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))))
  } else {
    getwd()
  }
  for (f in list.files(file.path(repo_root, "R"), pattern = "[.]R$", full.names = TRUE)) {
    source(f)
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

# --- 4. Read Fortran-order binary LD matrix (.npy) ----------------------------
R <- read_ld_matrix(matrix_path, ids_path = matrix_ids)
if (!"LD_R2_TO_LEAD" %in% names(variants) && nrow(R) == nrow(variants)) {
  lead_idx <- which.max(ifelse(is.na(variants$LP), -Inf, as.numeric(variants$LP)))
  if (length(lead_idx) == 1L) {
    variants$LD_R2_TO_LEAD <- as.numeric(R[, lead_idx]^2)
  }
}
io_seconds <- (proc.time() - t_io_start)[["elapsed"]]

# --- 5. Run greedy p-LD selection across r^2 thresholds -----------------------
t_fit_start <- proc.time()
pld_df <- prioritize_pld(variants, R, r2_thresholds)
fit_seconds <- (proc.time() - t_fit_start)[["elapsed"]]

selected_by_r2 <- stats::setNames(
  lapply(r2_thresholds, function(t) {
    as.integer(sum(pld_df$selected[pld_df$r2_threshold == t]))
  }),
  as.character(r2_thresholds)
)

# --- 6. Save standalone region_result.rds & method_qc.json --------------------
t_save_start <- proc.time()
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
region_payload <- list(
  region_id = region_id,
  region_type = region_type,
  result = pld_df
)
saveRDS(region_payload, file = out_rds, compress = rds_compress)
save_seconds <- (proc.time() - t_save_start)[["elapsed"]]

proc_delta <- proc.time() - t_script_start
gc_mat <- gc(FALSE)
peak_ram_mb <- round(sum(gc_mat[, 6L] * c(56, 8)) / (1024^2), 1)

qc_payload <- list(
  method_id = "p_ld",
  region_id = region_id,
  region_type = region_type,
  variants = nrow(variants),
  selected_by_r2 = selected_by_r2,
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

if (!is.null(qc_json) && nzchar(qc_json)) {
  dir.create(dirname(qc_json), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(qc_payload, qc_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

if (interactive()) {
  message(sprintf(
    "[p-ld] %s (%d variants): selected r2<0.1=%d, r2<0.01=%d, r2<0.001=%d | fit=%.2fs, total=%.2fs",
    region_id, nrow(variants),
    selected_by_r2[["0.1"]], selected_by_r2[["0.01"]], selected_by_r2[["0.001"]],
    fit_seconds, proc_delta[["elapsed"]]
  ))
}
