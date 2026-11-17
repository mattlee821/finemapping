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
if (!options$region_type %in% c("gene_region", "genome_wide_locus")) {
  stop("Unsupported --region-type.", call. = FALSE)
}

matrix_path <- options$matrix
if (any(methods %in% c("p_ld", "susie_rss", "susie_rrss")) &&
    (is.null(matrix_path) || !file.exists(matrix_path))) {
  stop("The selected methods require --matrix.", call. = FALSE)
}

R <- NULL
if (!is.null(matrix_path)) {
  matrix_frame <- utils::read.delim(
    matrix_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ids <- as.character(matrix_frame[[1L]])
  R <- as.matrix(matrix_frame[-1L])
  storage.mode(R) <- "double"
  rownames(R) <- ids
  colnames(R) <- names(matrix_frame)[-1L]
}

outdir <- options$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
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
summary_rows <- vector("list", length(methods))

for (i in seq_along(methods)) {
  method <- methods[[i]]
  method_dir <- file.path(outdir, method)
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

  if (method == "p_value") {
    result <- prioritize_p(variants)
    write_tsv(result, file.path(method_dir, "variants.tsv"))
    summary_rows[[i]] <- data.frame(
      method_id = method,
      variants = nrow(result),
      selected = sum(result$selected),
      credible_sets = NA_integer_
    )
  } else if (method == "p_ld") {
    result <- prioritize_pld(variants, R, thresholds)
    write_tsv(result, file.path(method_dir, "variants.tsv"))
    summary_rows[[i]] <- data.frame(
      method_id = method,
      variants = nrow(result),
      selected = sum(result$selected),
      credible_sets = NA_integer_
    )
    jsonlite::write_json(
      list(r2_thresholds = thresholds, definition = "independent greedy runs"),
      file.path(method_dir, "diagnostics.json"),
      auto_unbox = TRUE,
      pretty = TRUE
    )
  } else {
    args <- list(variants = variants, R = R)
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
    summary_rows[[i]] <- data.frame(
      method_id = method,
      variants = nrow(result$variants),
      selected = NA_integer_,
      credible_sets = nrow(result$credible_sets)
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
}

write_tsv(do.call(rbind, summary_rows), file.path(outdir, "method_summary.tsv"))
