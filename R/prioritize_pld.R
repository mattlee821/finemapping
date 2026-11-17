prioritize_pld <- function(variants, R, r2_thresholds) {
  if (!is.numeric(r2_thresholds) || length(r2_thresholds) == 0L ||
      anyNA(r2_thresholds) || any(!is.finite(r2_thresholds)) ||
      any(r2_thresholds <= 0 | r2_thresholds > 1) ||
      anyDuplicated(r2_thresholds)) {
    stop("`r2_thresholds` must be unique finite values in (0, 1].",
         call. = FALSE)
  }

  variants <- validate_gwas_stats(variants)
  if (nrow(variants) == 0L) {
    result <- variants[rep(NA_integer_, 0L), , drop = FALSE]
    result$method_id <- character()
    result$r2_threshold <- numeric()
    result$rank <- integer()
    result$selected <- logical()
    result$excluded_by_variant_id <- character()
    return(result)
  }
  if (anyNA(variants$LP)) {
    stop("A finite P-value, LP, or Z-score is required for every variant.",
         call. = FALSE)
  }

  ids <- variant_key(variants)
  validate_ld_matrix(R, variant_ids = ids, check_psd = FALSE)
  ordering <- order(-variants$LP, ids, method = "radix")
  rank <- integer(nrow(variants))
  rank[ordering] <- seq_along(ordering)
  ranked <- variants[ordering, , drop = FALSE]
  ranked_ids <- ids[ordering]

  results <- lapply(r2_thresholds, function(threshold) {
    n_var <- nrow(ranked)
    selected <- integer(n_var)
    n_sel <- 0L
    excluded_by <- rep(NA_character_, n_var)
    keep <- logical(n_var)

    for (i in seq_len(n_var)) {
      if (n_sel == 0L) {
        keep[[i]] <- TRUE
        n_sel <- 1L
        selected[[1L]] <- i
        next
      }

      sel_pos <- selected[seq_len(n_sel)]
      # R is validated symmetric; column access is contiguous in column-major R
      r_col <- R[ordering[sel_pos], ordering[[i]], drop = TRUE]
      r2 <- r_col * r_col
      high_ld <- which(r2 >= threshold)
      if (length(high_ld) == 0L) {
        keep[[i]] <- TRUE
        n_sel <- n_sel + 1L
        selected[[n_sel]] <- i
      } else {
        max_r2 <- max(r2[high_ld])
        first_max <- high_ld[which(r2[high_ld] == max_r2)[[1L]]]
        excluded_by[[i]] <- ranked_ids[sel_pos[[first_max]]]
      }
    }

    result <- ranked
    result$method_id <- "p_ld"
    result$r2_threshold <- threshold
    result$rank <- seq_len(nrow(result))
    result$selected <- keep
    result$excluded_by_variant_id <- excluded_by
    result
  })

  do.call(rbind, results)
}
