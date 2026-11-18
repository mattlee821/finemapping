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

parse_bool <- function(val, default = FALSE) {
  if (is.null(val)) return(default)
  tolower(as.character(val)) %in% c("1", "true", "yes", "y")
}

normalize_method_id <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("-", "_", x, fixed = TRUE)
  if (x == "polyfun_susie") x <- "polyfun"
  x
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
if (!"output_root" %in% names(options) || !dir.exists(options$output_root)) {
  stop("A valid --output-root is required.", call. = FALSE)
}
allow_empty <- parse_bool(options$allow_empty, default = FALSE)
keep_rds <- parse_bool(options$keep_rds, default = FALSE)
out_format <- tolower(options$format %||% "rda")
`%||%` <- function(x, y) if (is.null(x) || !nzchar(x)) y else x
out_format <- tolower(if (is.null(options$format) || !nzchar(options$format)) "rda" else options$format)
if (!out_format %in% c("rda", "rds", "both")) {
  stop("--format must be 'rda', 'rds', or 'both'.", call. = FALSE)
}

search_root <- if (!is.null(options$input_root) && nzchar(options$input_root)) {
  options$input_root
} else {
  options$output_root
}

allowed_methods <- c("p_ld", "susie_rss", "susie_rrss", "carma", "polyfun")
default_active <- c("p_ld", "susie_rss", "susie_rrss", "carma", "polyfun")
requested_methods <- if (!is.null(options$methods) && nzchar(options$methods)) {
  raw_m <- unlist(strsplit(options$methods, "[,[:space:]]+"))
  raw_m <- raw_m[nzchar(raw_m)]
  unique(vapply(raw_m, normalize_method_id, character(1)))
} else {
  NULL
}
if (!is.null(requested_methods) && any(!requested_methods %in% allowed_methods)) {
  bad <- setdiff(requested_methods, allowed_methods)
  stop("Unsupported method(s): ", paste(bad, collapse = ", "), call. = FALSE)
}

result_files <- if (dir.exists(search_root)) {
  list.files(
    search_root,
    pattern = "^region_result\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
} else {
  character()
}

if (!is.null(requested_methods) && length(result_files) > 0L) {
  m_dirs <- vapply(result_files, function(f) normalize_method_id(basename(dirname(f))), character(1))
  result_files <- result_files[m_dirs %in% requested_methods]
}

if (length(result_files) == 0L && !allow_empty) {
  stop("No per-region method results were found.", call. = FALSE)
}

allowed_region_types <- c(
  "cis_gene", "lead_locus", "mhc", "gene_region", "genome_wide_locus"
)

by_method <- list()
if (allow_empty) {
  seed_methods <- if (!is.null(requested_methods)) requested_methods else default_active
  for (m in seed_methods) {
    by_method[[m]] <- stats::setNames(list(), character())
  }
}

for (result_file in result_files) {
  method_id <- normalize_method_id(basename(dirname(result_file)))
  if (!method_id %in% allowed_methods) {
    stop("Unexpected method result path: ", result_file, call. = FALSE)
  }
  region <- readRDS(result_file)
  valid_region <- is.list(region) &&
    all(c("region_id", "region_type", "result") %in% names(region)) &&
    is.character(region$region_id) && length(region$region_id) == 1L &&
    nzchar(region$region_id) &&
    is.character(region$region_type) && length(region$region_type) == 1L &&
    region$region_type %in% allowed_region_types &&
    !is.null(region$result)
  if (!valid_region) {
    stop("Malformed per-region result: ", result_file, call. = FALSE)
  }
  if (method_id %in% c("p_value", "p_ld")) {
    if (!is.data.frame(region$result) ||
        !"method_id" %in% names(region$result) ||
        !identical(unique(region$result$method_id), method_id)) {
      stop("Method ID mismatch in per-region result: ", result_file,
           call. = FALSE)
    }
  } else {
    if (!is.list(region$result) ||
        !identical(region$result$provenance$method_id, method_id)) {
      stop("Method ID mismatch in per-region result: ", result_file,
           call. = FALSE)
    }
  }
  if (is.null(by_method[[method_id]])) {
    by_method[[method_id]] <- list()
  }
  if (region$region_id %in% names(by_method[[method_id]])) {
    stop(
      "Duplicate region ID for ", method_id, ": ", region$region_id,
      call. = FALSE
    )
  }
  by_method[[method_id]][[region$region_id]] <- region
}

output_dir <- file.path(options$output_root, "methods")
if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) &&
    !dir.exists(output_dir)) {
  stop("Could not create aggregate methods directory.", call. = FALSE)
}

method_stems <- c(
  p_value = "p-value",
  p_ld = "p-ld",
  susie_rss = "susie-rss",
  susie_rrss = "susie-rrss",
  carma = "carma",
  polyfun = "polyfun"
)
compression <- if (!is.null(options$compress) && nzchar(options$compress)) {
  options$compress
} else {
  "xz"
}

save_one_method <- function(method_id) {
  regions <- by_method[[method_id]]
  method_object <- list(
    schema_version = 1L,
    method_id = method_id,
    regions = regions
  )
  stem <- method_stems[[method_id]]
  if (out_format %in% c("rda", "both")) {
    environment <- new.env(parent = emptyenv())
    assign(method_id, method_object, envir = environment)
    save(
      list = method_id,
      envir = environment,
      file = file.path(output_dir, paste0(stem, ".rda")),
      compress = compression
    )
  }
  if (out_format %in% c("rds", "both")) {
    saveRDS(
      method_object,
      file = file.path(output_dir, paste0(stem, ".rds")),
      compress = compression
    )
  }
  TRUE
}

m_names <- names(by_method)
if (.Platform$OS.type == "unix" && length(m_names) > 1L) {
  res_par <- parallel::mclapply(
    m_names,
    save_one_method,
    mc.cores = min(length(m_names), 4L)
  )
  bad_par <- Filter(function(x) inherits(x, "try-error") || !isTRUE(x), res_par)
  if (length(bad_par) > 0L) {
    stop("Parallel method save failed: ", as.character(bad_par[[1L]]), call. = FALSE)
  }
} else {
  lapply(m_names, save_one_method)
}

if (!keep_rds && length(result_files) > 0L) {
  unlink_failures <- unlink(result_files)
  if (unlink_failures > 0L || any(file.exists(result_files))) {
    stop("Could not remove temporary per-region RDS files.", call. = FALSE)
  }
}
