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

options <- parse_options(commandArgs(trailingOnly = TRUE))
if (!"output_root" %in% names(options) || !dir.exists(options$output_root)) {
  stop("A valid --output-root is required.", call. = FALSE)
}

result_files <- list.files(
  options$output_root,
  pattern = "^region_result\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(result_files) == 0L) {
  stop("No per-region method results were found.", call. = FALSE)
}

by_method <- list()
for (result_file in result_files) {
  method_id <- basename(dirname(result_file))
  if (!method_id %in% c("p_value", "p_ld", "susie_rss", "susie_rrss")) {
    stop("Unexpected method result path: ", result_file, call. = FALSE)
  }
  region <- readRDS(result_file)
  valid_region <- is.list(region) &&
    all(c("region_id", "region_type", "result") %in% names(region)) &&
    is.character(region$region_id) && length(region$region_id) == 1L &&
    nzchar(region$region_id) &&
    is.character(region$region_type) && length(region$region_type) == 1L &&
    region$region_type %in% c("gene_region", "genome_wide_locus") &&
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

method_files <- c(
  p_value = "p-value.rda",
  p_ld = "p-ld.rda",
  susie_rss = "susie-rss.rda",
  susie_rrss = "susie-rrss.rda"
)
for (method_id in names(by_method)) {
  regions <- by_method[[method_id]]
  method_object <- list(
    schema_version = 1L,
    method_id = method_id,
    regions = regions
  )
  environment <- new.env(parent = emptyenv())
  assign(method_id, method_object, envir = environment)
  save(
    list = method_id,
    envir = environment,
    file = file.path(output_dir, method_files[[method_id]])
  )
}

unlink_failures <- unlink(result_files)
if (unlink_failures > 0L || any(file.exists(result_files))) {
  stop("Could not remove temporary per-region RDS files.", call. = FALSE)
}
