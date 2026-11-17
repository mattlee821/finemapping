validate_ld_matrix <- function(R, variant_ids = NULL, check_psd = FALSE,
                               tolerance = 1e-8, psd_max_variants = 5000L) {
  if (!is.matrix(R) || !is.numeric(R)) {
    stop("LD must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(R) != ncol(R)) {
    stop("LD matrix must be square.", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("`tolerance` must be one positive finite number.", call. = FALSE)
  }
  if (length(psd_max_variants) != 1L || !is.numeric(psd_max_variants) ||
      is.na(psd_max_variants) || psd_max_variants < 1L) {
    stop("`psd_max_variants` must be one positive integer.", call. = FALSE)
  }
  if (nrow(R) == 0L) {
    return(invisible("empty"))
  }

  row_ids <- rownames(R)
  col_ids <- colnames(R)
  if (is.null(row_ids) || is.null(col_ids) || anyNA(row_ids) ||
      anyNA(col_ids) || any(!nzchar(row_ids)) || any(!nzchar(col_ids))) {
    stop("LD matrix must have non-empty row and column variant identifiers.",
         call. = FALSE)
  }
  if (anyDuplicated(row_ids) || anyDuplicated(col_ids)) {
    stop("LD matrix variant identifiers must be unique.", call. = FALSE)
  }
  if (!identical(row_ids, col_ids)) {
    stop("LD row and column variant order must be identical.", call. = FALSE)
  }
  if (!is.null(variant_ids) && !identical(as.character(variant_ids), row_ids)) {
    stop("LD variant order does not match the supplied variant order.",
         call. = FALSE)
  }
  if (any(!is.finite(R))) {
    stop("LD contains missing or non-finite values; unavailable LD is not zero.",
         call. = FALSE)
  }
  if (max(abs(R - t(R))) > tolerance) {
    stop("LD matrix is not symmetric within tolerance.", call. = FALSE)
  }
  if (any(abs(diag(R) - 1) > tolerance)) {
    stop("LD diagonal values must be approximately 1.", call. = FALSE)
  }
  if (any(abs(R) > 1 + tolerance)) {
    stop("LD correlations must be within [-1, 1].", call. = FALSE)
  }

  psd_status <- "not_requested"
  if (check_psd) {
    if (nrow(R) <= as.integer(psd_max_variants)) {
      min_eigenvalue <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
      if (min_eigenvalue < -tolerance) {
        stop("LD matrix is not positive semidefinite; no repair was applied.",
             call. = FALSE)
      }
      psd_status <- "passed"
    } else {
      psd_status <- "skipped_large_m"
    }
  }
  res <- TRUE
  attr(res, "psd_status") <- psd_status
  invisible(res)
}
