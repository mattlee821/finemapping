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
    best_r2 <- rep(-1.0, n_var)
    best_sel_rank <- integer(n_var)
    excluded_by <- rep(NA_character_, n_var)
    keep <- logical(n_var)

    # Selected-push sieve: only selected variants read their contiguous LD column
    for (i in seq_len(n_var)) {
      if (best_r2[[i]] >= threshold) {
        excluded_by[[i]] <- ranked_ids[[best_sel_rank[[i]]]]
        next
      }
      keep[[i]] <- TRUE
      if (i == n_var) break

      r_col <- R[, ordering[[i]], drop = TRUE]
      tail_idx <- (i + 1L):n_var
      r_tail <- r_col[ordering[tail_idx]]
      r2_tail <- r_tail * r_tail
      upd <- which(r2_tail >= threshold & r2_tail > best_r2[tail_idx])
      if (length(upd) > 0L) {
        pos_upd <- tail_idx[upd]
        best_r2[pos_upd] <- r2_tail[upd]
        best_sel_rank[pos_upd] <- i
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
