prioritize_p <- function(variants, top_n = NULL) {
  if (!is.null(top_n) &&
      (length(top_n) != 1L || is.na(top_n) || !is.numeric(top_n) ||
       top_n < 1 || top_n != floor(top_n))) {
    stop("`top_n` must be NULL or one positive integer.", call. = FALSE)
  }

  variants <- validate_gwas_stats(variants)
  if (nrow(variants) == 0L) {
    variants$rank <- integer()
    variants$method_id <- character()
    variants$selected <- logical()
    variants$selection_reason <- character()
    return(variants)
  }
  if (anyNA(variants$LP)) {
    stop("A finite P-value, LP, or Z-score is required for every variant.",
         call. = FALSE)
  }

  ids <- variant_key(variants)
  ordering <- order(-variants$LP, ids, method = "radix")
  result <- variants[ordering, , drop = FALSE]
  result$rank <- seq_len(nrow(result))
  result$method_id <- "p_value"
  result$selected <- if (is.null(top_n)) {
    rep(TRUE, nrow(result))
  } else {
    result$rank <= top_n
  }
  result$selection_reason <- if (is.null(top_n)) {
    rep("all_ranked_variants", nrow(result))
  } else {
    ifelse(result$selected, "top_n", "outside_top_n")
  }
  result
}
