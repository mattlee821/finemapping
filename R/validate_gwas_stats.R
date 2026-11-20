validate_gwas_stats <- function(variants, require_p = FALSE, require_z = FALSE) {
  if (!is.data.frame(variants)) {
    stop("`variants` must be a data frame.", call. = FALSE)
  }

  required <- c("CHROM", "POS", "REF", "ALT")
  missing <- setdiff(required, names(variants))
  if (length(missing) > 0) {
    stop("Missing required variant field(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  out <- variants
  if (nrow(out) == 0L) {
    return(out)
  }

  for (field in c("CHROM", "REF", "ALT")) {
    value <- as.character(out[[field]])
    if (anyNA(value) || any(!nzchar(value)) || any(value == ".")) {
      stop("Missing or invalid variant identity in ", field, ".", call. = FALSE)
    }
    out[[field]] <- value
  }

  out$POS <- parse_numeric_field(out$POS, "POS")
  if (anyNA(out$POS) || any(out$POS < 1) || any(out$POS != floor(out$POS)) ||
      any(out$POS > .Machine$integer.max)) {
    stop("POS must contain positive integer coordinates.", call. = FALSE)
  }
  out$POS <- as.integer(out$POS)

  if (any(grepl(",", out$ALT, fixed = TRUE))) {
    stop("Multiallelic ALT records are not supported; split them upstream.",
         call. = FALSE)
  }

  ids <- variant_key(out)
  duplicate_ids <- unique(ids[duplicated(ids)])
  if (length(duplicate_ids) > 0) {
    stop("Duplicate allele-aware variant key(s): ",
         paste(utils::head(duplicate_ids, 5L), collapse = ", "), call. = FALSE)
  }

  for (field in intersect(c("P", "LP", "EZ", "Z", "ES", "SE", "SS",
                            "AF", "NC", "SI"), names(out))) {
    if (is.logical(out[[field]]) && all(is.na(out[[field]]))) {
      out[[field]] <- as.numeric(out[[field]])
    }
    if (!is.numeric(out[[field]])) {
      stop(field, " must be numeric.", call. = FALSE)
    }
    finite_or_missing <- is.finite(out[[field]]) | is.na(out[[field]])
    if (any(!finite_or_missing)) {
      stop(field, " contains non-finite values.", call. = FALSE)
    }
  }

  if ("P" %in% names(out) &&
      any(!is.na(out$P) & (out$P < 0 | out$P > 1))) {
    stop("P must be between 0 and 1.", call. = FALSE)
  }
  if ("LP" %in% names(out) &&
      any(!is.na(out$LP) & out$LP < 0)) {
    stop("LP must be non-negative.", call. = FALSE)
  }
  if ("SE" %in% names(out) &&
      any(!is.na(out$SE) & out$SE <= 0)) {
    stop("SE must be positive when supplied.", call. = FALSE)
  }
  if ("SS" %in% names(out) &&
      any(!is.na(out$SS) & out$SS <= 0)) {
    stop("SS must be positive when supplied.", call. = FALSE)
  }
  if ("AF" %in% names(out) &&
      any(!is.na(out$AF) & (out$AF < 0 | out$AF > 1))) {
    stop("AF must be between 0 and 1.", call. = FALSE)
  }

  z <- rep(NA_real_, nrow(out))
  if ("EZ" %in% names(out)) {
    z <- out$EZ
  } else if ("Z" %in% names(out)) {
    z <- out$Z
  }
  if (anyNA(z) && all(c("ES", "SE") %in% names(out))) {
    derive <- is.na(z) & !is.na(out$ES) & !is.na(out$SE)
    z[derive] <- out$ES[derive] / out$SE[derive]
  }
  if (any(!is.na(z) & !is.finite(z))) {
    stop("Derived or supplied Z-scores must be finite.", call. = FALSE)
  }
  if (!"EZ" %in% names(out) && any(!is.na(z))) {
    out$EZ <- z
  } else if ("EZ" %in% names(out)) {
    out$EZ[is.na(out$EZ)] <- z[is.na(out$EZ)]
  }

  lp <- if ("LP" %in% names(out)) out$LP else rep(NA_real_, nrow(out))
  if ("P" %in% names(out)) {
    derive <- is.na(lp) & !is.na(out$P) & out$P > 0
    lp[derive] <- -log10(out$P[derive])
  }
  derive_from_z <- is.na(lp) & !is.na(z)
  if (any(derive_from_z)) {
    lp[derive_from_z] <- -stats::pchisq(
      z[derive_from_z]^2,
      df = 1,
      lower.tail = FALSE,
      log.p = TRUE
    ) / log(10)
  }
  if (any(!is.na(lp) & (!is.finite(lp) | lp < 0))) {
    stop("LP must be finite and non-negative.", call. = FALSE)
  }
  if (!"LP" %in% names(out)) {
    out$LP <- lp
  } else {
    out$LP[is.na(out$LP)] <- lp[is.na(out$LP)]
  }

  if (require_p && anyNA(out$LP)) {
    stop("A finite P-value, LP, or Z-score is required for every variant.",
         call. = FALSE)
  }
  if (require_z && anyNA(z)) {
    stop("A signed Z-score is required for every variant; supply EZ or valid ES/SE.",
         call. = FALSE)
  }

  out
}
