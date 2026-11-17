extract_info_field <- function(info, field) {
  tokens <- strsplit(as.character(info), ";", fixed = TRUE)
  vapply(tokens, function(parts) {
    key <- paste0(field, "=")
    match <- which(startsWith(parts, key))
    if (length(match) > 1L) {
      stop("Repeated INFO/", field, " value.", call. = FALSE)
    }
    if (length(match) == 1L) {
      return(substring(parts[match], nchar(key) + 1L))
    }
    if (field %in% parts) {
      return("TRUE")
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}

parse_format_fields <- function(records, format_column, sample_column) {
  if (nrow(records) == 0L) {
    return(records)
  }

  format_values <- as.character(records[[format_column]])
  sample_values <- as.character(records[[sample_column]])
  format_groups <- unique(format_values)
  reserved_columns <- unique(c(names(records), "SNPID", "HM"))
  format_output_sources <- character()

  for (format_value in format_groups) {
    rows <- which(format_values == format_value)
    if (is.na(format_value) || format_value %in% c("", ".")) {
      next
    }
    tags <- strsplit(format_value, ":", fixed = TRUE)[[1]]
    if (anyDuplicated(tags)) {
      stop("Duplicate FORMAT tag in VCF record(s).", call. = FALSE)
    }
    values <- strsplit(sample_values[rows], ":", fixed = TRUE)
    if (any(lengths(values) > length(tags))) {
      stop("Sample field contains more values than its FORMAT declaration.",
           call. = FALSE)
    }

    for (column in tags) {
      output_column <- if (column %in% reserved_columns) {
        paste0("FORMAT_", column)
      } else {
        column
      }
      if (output_column %in% reserved_columns) {
        stop("FORMAT tag cannot be represented without colliding with VCF fields: ",
             column, call. = FALSE)
      }
      if (output_column %in% names(format_output_sources) &&
          format_output_sources[[output_column]] != column) {
        stop("Multiple FORMAT tags map to output field: ", output_column,
             call. = FALSE)
      }
      format_output_sources[output_column] <- column
      if (!output_column %in% names(records)) {
        records[[output_column]] <- NA_character_
      }
      tag_index <- match(column, tags)
      records[[output_column]][rows] <- vapply(values, function(record_values) {
        if (length(record_values) < tag_index ||
            record_values[[tag_index]] %in% c("", ".")) {
          return(NA_character_)
        }
        record_values[[tag_index]]
      }, character(1), USE.NAMES = FALSE)
    }
  }

  records
}

read_gwas_vcf <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    stop("`path` must be one non-empty file path.", call. = FALSE)
  }
  if (!file.exists(path) || dir.exists(path)) {
    stop("GWAS-VCF file does not exist: ", path, call. = FALSE)
  }

  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  metadata_lines <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      stop("VCF header is missing the #CHROM column line.", call. = FALSE)
    }
    if (startsWith(line, "##")) {
      metadata_lines <- c(metadata_lines, line)
      next
    }
    if (startsWith(line, "#CHROM\t")) {
      columns <- strsplit(substring(line, 2L), "\t", fixed = TRUE)[[1]]
      break
    }
    stop("Unexpected line before the VCF #CHROM header.", call. = FALSE)
  }

  fixed_columns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                     "INFO", "FORMAT")
  if (!all(fixed_columns %in% columns)) {
    stop("VCF header must contain the standard fixed fields and FORMAT.",
         call. = FALSE)
  }
  if (anyDuplicated(columns)) {
    stop("VCF header contains duplicate column names.", call. = FALSE)
  }
  sample_columns <- columns[seq_along(columns) > match("FORMAT", columns)]
  if (length(sample_columns) != 1L) {
    stop("Exactly one GWAS-VCF sample column is supported per input file.",
         call. = FALSE)
  }

  records <- utils::read.delim(
    con,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    colClasses = "character",
    na.strings = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(records) != length(columns)) {
    stop("VCF data column count does not match the #CHROM header.",
         call. = FALSE)
  }
  names(records) <- columns

  if (nrow(records) > 0L) {
    records <- parse_format_fields(records, "FORMAT", sample_columns[[1]])
    if ("SNPID" %in% names(records) || "HM" %in% names(records)) {
      stop("INFO-derived SNPID/HM would collide with a FORMAT tag.",
           call. = FALSE)
    }
    records$SNPID <- extract_info_field(records$INFO, "SNPID")
    records$HM <- extract_info_field(records$INFO, "HM")
  }

  numeric_fields <- intersect(c("POS", "ES", "SE", "LP", "AF", "SS", "NC",
                                "EZ", "SI"), names(records))
  for (field in numeric_fields) {
    records[[field]] <- parse_numeric_field(records[[field]], field)
  }
  if (nrow(records) > 0L) {
    records <- validate_gwas_stats(records)
  }

  meta_value <- function(key) {
    prefix <- paste0("##", key, "=")
    matching <- metadata_lines[startsWith(metadata_lines, prefix)]
    if (length(matching) > 1L) {
      stop("Repeated VCF metadata key: ", key, call. = FALSE)
    }
    if (length(matching) == 0L) {
      return(NA_character_)
    }
    substring(matching, nchar(prefix) + 1L)
  }

  structure(
    list(
      variants = records,
      metadata = list(
        header_lines = metadata_lines,
        fileformat = meta_value("fileformat"),
        reference = meta_value("reference"),
        source = meta_value("source")
      ),
      sample_id = sample_columns[[1]],
      path = normalizePath(path, mustWork = TRUE)
    ),
    class = "gwas_vcf"
  )
}
