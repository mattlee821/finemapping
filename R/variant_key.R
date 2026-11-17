variant_key <- function(variants) {
  required <- c("CHROM", "POS", "REF", "ALT")
  missing <- setdiff(required, names(variants))
  if (length(missing) > 0) {
    stop("Missing variant identity field(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  if ("LD_VARIANT_ID" %in% names(variants)) {
    ld_ids <- as.character(variants$LD_VARIANT_ID)
    if (anyNA(ld_ids) || any(!nzchar(ld_ids))) {
      stop("LD_VARIANT_ID must be present for every variant when supplied.",
           call. = FALSE)
    }
    ids <- ld_ids
  } else {
    ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT,
                 sep = ":")
  }
  ids
}

parse_numeric_field <- function(x, field) {
  x <- as.character(x)
  missing <- is.na(x) | x %in% c("", ".")
  numeric_pattern <- paste0(
    "^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
  )
  valid <- missing | grepl(numeric_pattern, x, perl = TRUE)
  if (any(!valid)) {
    bad_values <- unique(x[!valid])
    stop("Invalid numeric value(s) in ", field, ": ",
         paste(utils::head(bad_values, 5L), collapse = ", "), call. = FALSE)
  }

  values <- rep(NA_real_, length(x))
  values[!missing] <- as.numeric(x[!missing])
  if (any(!is.finite(values[!missing]))) {
    stop("Non-finite numeric value in ", field, call. = FALSE)
  }
  values
}
