#!/usr/bin/env Rscript

parse_options <- function(args) {
  if (length(args) %% 2L != 0L) {
    stop("Options must be supplied as --name value pairs.", call. = FALSE)
  }
  options <- list()
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- gsub("-", "_", sub("^--", "", args[[i]]), fixed = TRUE)
    if (!startsWith(args[[i]], "--") || key %in% names(options)) {
      stop("Invalid or repeated option: ", args[[i]], call. = FALSE)
    }
    options[[key]] <- args[[i + 1L]]
  }
  options
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

parse_bool <- function(val, default = TRUE) {
  if (is.null(val)) return(default)
  tolower(as.character(val)) %in% c("1", "true", "yes", "y")
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
required <- c(
  "variants", "outdir", "methods", "region_id", "region_type"
)
missing <- setdiff(required, names(options))
if (length(missing) > 0L) {
  stop("Missing required option(s): ", paste(missing, collapse = ", "),
       call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
project_root <- dirname(dirname(script_path))
for (file in list.files(file.path(project_root, "R"), pattern = "\\.R$",
                        full.names = TRUE)) {
  source(file)
}

variants <- utils::read.delim(
  options$variants,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
for (bcol in intersect(
  c("mhc_region", "mhc_core_grc", "PANEL_ALLELE_SWAPPED"),
  names(variants)
)) {
  if (!is.logical(variants[[bcol]])) {
    variants[[bcol]] <- tolower(as.character(variants[[bcol]])) == "true"
  }
}
if (nrow(variants) == 0L) {
  stop("The region has no LD-matched variants.", call. = FALSE)
}

methods <- strsplit(options$methods, ",", fixed = TRUE)[[1L]]
allowed_methods <- c("p_value", "p_ld", "susie_rss", "susie_rrss")
if (any(!methods %in% allowed_methods) || anyDuplicated(methods)) {
  stop("Unsupported or repeated method in --methods.", call. = FALSE)
}
if (!nzchar(options$region_id)) {
  stop("--region-id must not be empty.", call. = FALSE)
}
allowed_region_types <- c(
  "cis_gene", "lead_locus", "mhc", "gene_region", "genome_wide_locus"
)
if (!options$region_type %in% allowed_region_types) {
  stop("Unsupported --region-type: ", options$region_type, call. = FALSE)
}

matrix_path <- options$matrix
if (any(methods %in% c("p_ld", "susie_rss", "susie_rrss")) &&
    (is.null(matrix_path) || !file.exists(matrix_path))) {
  stop("The selected methods require --matrix.", call. = FALSE)
}

R <- NULL
if (!is.null(matrix_path)) {
  R <- read_ld_matrix(matrix_path, ids_path = options$matrix_ids)
  if (!"LD_R2_TO_LEAD" %in% names(variants) && nrow(R) == nrow(variants)) {
    lead_idx <- which.max(
      ifelse(is.na(variants$LP), -Inf, as.numeric(variants$LP))
    )
    if (length(lead_idx) == 1L) {
      variants$LD_R2_TO_LEAD <- as.numeric(R[, lead_idx]^2)
    }
  }
}

outdir <- options$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
emit_tsv <- parse_bool(options$write_tsv, default = TRUE)

thresholds <- if (is.null(options$r2_thresholds)) {
  c(0.1, 0.01, 0.001)
} else {
  as.numeric(strsplit(options$r2_thresholds, ",", fixed = TRUE)[[1L]])
}
reference_sample_size <- if (is.null(options$reference_sample_size)) {
  NULL
} else {
  as.numeric(options$reference_sample_size)
}
effective_n <- if (is.null(options$effective_n)) {
  NULL
} else {
  as.numeric(options$effective_n)
}
max_causal <- if (is.null(options$max_causal)) {
  10L
} else {
  as.integer(options$max_causal)
}
max_iter <- if (is.null(options$max_iter)) {
  200L
} else {
  as.integer(options$max_iter)
}
min_abs_corr <- if (is.null(options$min_abs_corr)) {
  0.5
} else {
  as.numeric(options$min_abs_corr)
}
coverage <- if (is.null(options$coverage)) {
  0.95
} else {
  as.numeric(options$coverage)
}

summary_rows <- vector("list", length(methods))
method_qc <- list()

for (i in seq_along(methods)) {
  method <- methods[[i]]
  method_dir <- file.path(outdir, method)
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

  if (method == "p_value") {
    result <- prioritize_p(variants)
    if (emit_tsv) {
      write_tsv(result, file.path(method_dir, "variants.tsv"))
    }
    summary_rows[[i]] <- data.frame(
      method_id = method,
      variants = nrow(result),
      selected = sum(result$selected),
      credible_sets = NA_integer_,
      stringsAsFactors = FALSE
    )
    method_qc[[method]] <- list(
      variants = nrow(result),
      selected = sum(result$selected)
    )
  } else if (method == "p_ld") {
    result <- prioritize_pld(variants, R, thresholds)
    if (emit_tsv) {
      write_tsv(result, file.path(method_dir, "variants.tsv"))
      jsonlite::write_json(
        list(r2_thresholds = thresholds, definition = "independent greedy runs"),
        file.path(method_dir, "diagnostics.json"),
        auto_unbox = TRUE,
        pretty = TRUE
      )
    }
    by_thresh <- stats::setNames(
      lapply(thresholds, function(t) {
        sum(result$selected[result$r2_threshold == t])
      }),
      as.character(thresholds)
    )
    summary_rows[[i]] <- data.frame(
      method_id = method,
      variants = nrow(result),
      selected = sum(result$selected),
      credible_sets = NA_integer_,
      stringsAsFactors = FALSE
    )
    method_qc[[method]] <- list(
      variants = nrow(variants),
      selected_by_r2 = by_thresh
    )
  } else {
    args <- list(
      variants = variants,
      R = R,
      effective_n = effective_n,
      L = max_causal,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      max_iter = max_iter
    )
    if (method == "susie_rss") {
      result <- do.call(run_susie_rss, args)
    } else {
      if (is.null(reference_sample_size)) {
        stop("SuSiE-rRSS requires --reference-sample-size.", call. = FALSE)
      }
      result <- do.call(
        run_susie_rrss,
        c(args, list(reference_sample_size = reference_sample_size))
      )
    }

    if (emit_tsv) {
      write_tsv(result$variants, file.path(method_dir, "variants.tsv"))
      write_tsv(result$credible_sets, file.path(method_dir, "credible_sets.tsv"))
      write_tsv(
        result$credible_set_membership,
        file.path(method_dir, "credible_set_membership.tsv")
      )
      jsonlite::write_json(
        result$diagnostics,
        file.path(method_dir, "diagnostics.json"),
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null"
      )
      jsonlite::write_json(
        result$provenance,
        file.path(method_dir, "provenance.json"),
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null"
      )
    }
    summary_rows[[i]] <- data.frame(
      method_id = method,
      variants = nrow(result$variants),
      selected = NA_integer_,
      credible_sets = nrow(result$credible_sets),
      stringsAsFactors = FALSE
    )
    diag_compact <- result$diagnostics
    diag_compact$R_finite_diagnostics <- NULL
    method_qc[[method]] <- c(
      list(
        variants = nrow(result$variants),
        credible_sets = nrow(result$credible_sets),
        max_pip = max(result$variants$PIP, na.rm = TRUE),
        min_cs_purity = if (nrow(result$credible_sets) > 0L) {
          min(result$credible_sets$min_abs_corr, na.rm = TRUE)
        } else {
          NA_real_
        }
      ),
      diag_compact
    )
  }
  saveRDS(
    list(
      region_id = options$region_id,
      region_type = options$region_type,
      result = result
    ),
    file.path(method_dir, "region_result.rds")
  )
  rm(result)
  invisible(gc(FALSE))
}

if (emit_tsv) {
  write_tsv(do.call(rbind, summary_rows), file.path(outdir, "method_summary.tsv"))
}
jsonlite::write_json(
  method_qc,
  file.path(outdir, "method_qc.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  na = "null"
)
