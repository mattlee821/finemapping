resolve_effective_n <- function(variants, effective_n) {
  sample_sizes <- if ("SS" %in% names(variants)) {
    variants$SS[!is.na(variants$SS)]
  } else {
    numeric()
  }
  summary <- if (length(sample_sizes) == 0L) {
    NULL
  } else {
    c(
      n_variants_with_ss = length(sample_sizes),
      min_ss = min(sample_sizes),
      median_ss = stats::median(sample_sizes),
      max_ss = max(sample_sizes)
    )
  }

  if (!is.null(effective_n)) {
    if (length(effective_n) != 1L || !is.numeric(effective_n) ||
        is.na(effective_n) || !is.finite(effective_n) || effective_n <= 0) {
      stop("`effective_n` must be NULL or one positive finite number.",
           call. = FALSE)
    }
    return(list(
      value = effective_n,
      source = "explicit_override",
      summary = summary
    ))
  }

  if (!"SS" %in% names(variants)) {
    stop("SS is unavailable; supply an explicit `effective_n` override.",
         call. = FALSE)
  }
  if (length(sample_sizes) == 0L) {
    stop("No valid per-variant SS values; supply an explicit `effective_n` override.",
         call. = FALSE)
  }

  list(
    value = min(sample_sizes),
    source = "minimum_variant_ss",
    summary = summary
  )
}

resolve_gwas_variant_ids <- function(variants, fallback_ids) {
  if ("GWAS_VARIANT_ID" %in% names(variants) &&
      all(nzchar(as.character(variants$GWAS_VARIANT_ID))) &&
      !anyNA(variants$GWAS_VARIANT_ID)) {
    return(as.character(variants$GWAS_VARIANT_ID))
  }
  if (all(c("GWAS_CHROM", "GWAS_POS", "GWAS_REF", "GWAS_ALT") %in%
          names(variants))) {
    chrom <- sub("^chr", "", as.character(variants$GWAS_CHROM),
                 ignore.case = TRUE)
    return(paste(
      chrom,
      as.integer(variants$GWAS_POS),
      toupper(as.character(variants$GWAS_REF)),
      toupper(as.character(variants$GWAS_ALT)),
      sep = ":"
    ))
  }
  fallback_ids
}

name_susie_fit_variants <- function(fit, variant_ids) {
  p <- length(variant_ids)
  if (!is.null(fit$pip) && length(fit$pip) == p) {
    names(fit$pip) <- variant_ids
  }
  if (!is.null(fit$XtXr) && length(fit$XtXr) == p) {
    names(fit$XtXr) <- variant_ids
  }
  for (slot in c("alpha", "lbf_variable", "mu", "mu2")) {
    mat <- fit[[slot]]
    if (is.matrix(mat) && ncol(mat) == p) {
      colnames(mat) <- variant_ids
      fit[[slot]] <- mat
    }
  }
  fit$snp <- variant_ids
  if (!is.null(fit$sets) && !is.null(fit$sets$cs) && length(fit$sets$cs) > 0L) {
    for (i in seq_along(fit$sets$cs)) {
      idx <- as.integer(fit$sets$cs[[i]])
      names(idx) <- variant_ids[idx]
      fit$sets$cs[[i]] <- idx
    }
  }
  fit
}

run_susie_model <- function(variants, R, method_id, effective_n = NULL,
                            reference_sample_size = NULL, L = 10L,
                            coverage = 0.95, min_abs_corr = 0.5,
                            max_iter = 200L, refine = FALSE,
                            R_mismatch = "eb_mix",
                            psd_max_variants = 5000L) {
  if (!method_id %in% c("susie_rss", "susie_rrss")) {
    stop("Unsupported SuSiE method identifier.", call. = FALSE)
  }
  if (length(L) != 1L || !is.numeric(L) || is.na(L) || L < 1 ||
      L != floor(L)) {
    stop("`L` must be one positive integer.", call. = FALSE)
  }
  if (length(coverage) != 1L || !is.numeric(coverage) ||
      is.na(coverage) || !is.finite(coverage) || coverage <= 0 ||
      coverage > 1) {
    stop("`coverage` must be in (0, 1].", call. = FALSE)
  }
  if (!is.null(min_abs_corr) &&
      (length(min_abs_corr) != 1L || !is.numeric(min_abs_corr) ||
       is.na(min_abs_corr) || !is.finite(min_abs_corr) ||
       min_abs_corr < 0 || min_abs_corr > 1)) {
    stop("`min_abs_corr` must be NULL or a number in [0, 1].",
         call. = FALSE)
  }
  if (length(max_iter) != 1L || !is.numeric(max_iter) || is.na(max_iter) ||
      max_iter < 1 || max_iter != floor(max_iter)) {
    stop("`max_iter` must be one positive integer.", call. = FALSE)
  }
  if (!is.logical(refine) || length(refine) != 1L || is.na(refine)) {
    stop("`refine` must be one non-missing logical value.", call. = FALSE)
  }

  variants <- validate_gwas_stats(variants, require_p = TRUE, require_z = TRUE)
  if (nrow(variants) == 0L) {
    stop("SuSiE cannot fit a region with no variants.", call. = FALSE)
  }
  ids <- variant_key(variants)
  gwas_ids <- resolve_gwas_variant_ids(variants, ids)
  ld_check <- validate_ld_matrix(
    R,
    variant_ids = ids,
    check_psd = TRUE,
    psd_max_variants = psd_max_variants
  )
  psd_status <- attr(ld_check, "psd_status") %||% "passed"
  n <- resolve_effective_n(variants, effective_n)

  if (method_id == "susie_rrss") {
    if (length(reference_sample_size) != 1L ||
        !is.numeric(reference_sample_size) ||
        is.na(reference_sample_size) || !is.finite(reference_sample_size) ||
        reference_sample_size <= 0) {
      stop("SuSiE-rRSS requires a positive LD reference sample size B.",
           call. = FALSE)
    }
    if (reference_sample_size <= 2) {
      stop("LD reference sample size B must be greater than 2.",
           call. = FALSE)
    }
    if (length(R_mismatch) != 1L || !is.character(R_mismatch) ||
        is.na(R_mismatch) ||
        !R_mismatch %in% c("none", "eb", "eb_mix")) {
      stop("`R_mismatch` must be one of 'none', 'eb', or 'eb_mix'.",
           call. = FALSE)
    }
    R_finite <- reference_sample_size
  } else {
    R_finite <- NULL
    R_mismatch <- "none"
  }

  fit <- susieR::susie_rss(
    z = variants$EZ,
    R = R,
    n = n$value,
    L = as.integer(L),
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = as.integer(max_iter),
    refine = refine,
    R_finite = R_finite,
    R_mismatch = R_mismatch
  )
  cs <- susieR::susie_get_cs(
    fit,
    Xcorr = R,
    coverage = coverage,
    min_abs_corr = min_abs_corr
  )
  fit$sets <- cs
  fit <- name_susie_fit_variants(fit, gwas_ids)

  credible_sets <- if (length(cs$cs) == 0L) {
    data.frame(
      credible_set_id = character(),
      effect_index = integer(),
      size = integer(),
      coverage = numeric(),
      min_abs_corr = numeric(),
      mean_abs_corr = numeric(),
      median_abs_corr = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      credible_set_id = names(cs$cs),
      effect_index = as.integer(cs$cs_index),
      size = as.integer(lengths(cs$cs)),
      coverage = as.numeric(cs$coverage),
      min_abs_corr = cs$purity$min.abs.corr,
      mean_abs_corr = cs$purity$mean.abs.corr,
      median_abs_corr = cs$purity$median.abs.corr,
      stringsAsFactors = FALSE
    )
  }
  membership <- if (length(cs$cs) == 0L) {
    data.frame(
      credible_set_id = character(),
      effect_index = integer(),
      variant_id = character(),
      gwas_variant_id = character(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, lapply(seq_along(cs$cs), function(i) {
      members <- as.integer(cs$cs[[i]])
      data.frame(
        credible_set_id = names(cs$cs)[[i]],
        effect_index = as.integer(cs$cs_index[[i]]),
        variant_id = ids[members],
        gwas_variant_id = gwas_ids[members],
        stringsAsFactors = FALSE
      )
    }))
  }

  cs_lookup <- rep(NA_character_, nrow(variants))
  if (length(cs$cs) > 0L) {
    for (i in seq_along(cs$cs)) {
      members <- as.integer(cs$cs[[i]])
      cs_lookup[members] <- names(cs$cs)[[i]]
    }
  }

  diagnostics <- list(
    converged = isTRUE(fit$converged),
    convergence_status = if (is.null(fit$converged)) NA else fit$converged,
    convergence_reason = fit$convergence_reason %||% NA_character_,
    iterations = fit$niter %||% NA_integer_,
    final_elbo = if (length(fit$elbo) > 0L) utils::tail(fit$elbo, 1L) else NA_real_,
    psd_check = psd_status,
    R_reliability_flag = fit$R_finite_diagnostics$R_reliability_flag %||% NA,
    SER_fallback = fit$R_finite_diagnostics$SER_fallback %||% NA,
    artifact_flag = fit$R_finite_diagnostics$artifact_flag %||% NA,
    R_sensitivity_flag = fit$R_sensitivity_flag %||% NA,
    R_finite_diagnostics = fit$R_finite_diagnostics
  )

  lead_idx <- which.max(
    ifelse(is.na(variants$LP), -Inf, as.numeric(variants$LP))
  )
  ld_r2_to_lead <- if (length(lead_idx) == 1L && nrow(R) == nrow(variants)) {
    as.numeric(R[, lead_idx]^2)
  } else {
    rep(NA_real_, nrow(variants))
  }
  expected_z <- if (!is.null(fit$XtXr) && length(fit$XtXr) == nrow(variants) &&
                    n$value > 1) {
    as.numeric(fit$XtXr / sqrt(n$value - 1))
  } else {
    rep(NA_real_, nrow(variants))
  }

  variant_results <- variants
  if (!"GWAS_VARIANT_ID" %in% names(variant_results)) {
    variant_results$GWAS_VARIANT_ID <- gwas_ids
  }
  if (!"LD_R2_TO_LEAD" %in% names(variant_results)) {
    variant_results$LD_R2_TO_LEAD <- ld_r2_to_lead
  }
  variant_results$method_id <- method_id
  variant_results$PIP <- as.numeric(fit$pip)
  variant_results$EXPECTED_Z <- expected_z
  variant_results$credible_set_id <- cs_lookup

  list(
    variants = variant_results,
    credible_sets = credible_sets,
    credible_set_membership = membership,
    diagnostics = diagnostics,
    provenance = list(
      method_id = method_id,
      method_version = as.character(utils::packageVersion("susieR")),
      parameters = list(
        L = as.integer(L),
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        max_iter = as.integer(max_iter),
        refine = refine,
        effective_n = n$value,
        effective_n_source = n$source,
        effective_n_summary = n$summary,
        R_finite = R_finite,
        R_mismatch = R_mismatch,
        psd_max_variants = as.integer(psd_max_variants)
      )
    ),
    fit = fit
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

run_susie_rss <- function(variants, R, effective_n = NULL, L = 10L,
                          coverage = 0.95, min_abs_corr = 0.5,
                          max_iter = 200L, refine = FALSE,
                          psd_max_variants = 5000L) {
  run_susie_model(
    variants = variants,
    R = R,
    method_id = "susie_rss",
    effective_n = effective_n,
    L = L,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter,
    refine = refine,
    psd_max_variants = psd_max_variants
  )
}

run_susie_rrss <- function(variants, R, reference_sample_size,
                           effective_n = NULL, L = 10L, coverage = 0.95,
                           min_abs_corr = 0.5, max_iter = 200L,
                           refine = FALSE, R_mismatch = "eb_mix",
                           psd_max_variants = 5000L) {
  run_susie_model(
    variants = variants,
    R = R,
    method_id = "susie_rrss",
    effective_n = effective_n,
    reference_sample_size = reference_sample_size,
    L = L,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter,
    refine = refine,
    R_mismatch = R_mismatch,
    psd_max_variants = psd_max_variants
  )
}
