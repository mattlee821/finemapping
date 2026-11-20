resolve_effective_n <- function(variants, effective_n, allow_missing = FALSE) {
  sample_sizes <- if ("SS" %in% names(variants) && is.numeric(variants$SS)) {
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

  if (length(sample_sizes) == 0L) {
    if (isTRUE(allow_missing)) {
      return(list(
        value = NULL,
        source = "z_only_unreported_n",
        summary = NULL
      ))
    }
    if (!"SS" %in% names(variants)) {
      stop("SS is unavailable; supply an explicit `effective_n` override.",
           call. = FALSE)
    }
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

compute_credible_sets_fast <- function(fit, R, coverage = 0.95,
                                       min_abs_corr = 0.5,
                                       max_exact_purity = 500L,
                                       ld_scale = 1) {
  raw_cs <- susieR::susie_get_cs(
    fit,
    Xcorr = NULL,
    coverage = coverage,
    min_abs_corr = NULL,
    check_symmetric = FALSE
  )
  if (is.null(raw_cs$cs) || length(raw_cs$cs) == 0L) {
    return(list(
      cs = NULL,
      coverage = NULL,
      requested_coverage = coverage
    ))
  }
  cs_list <- raw_cs$cs
  effect_indices <- as.integer(sub("^L", "", names(cs_list)))
  claimed_cov <- as.numeric(raw_cs$coverage)

  purity_mat <- matrix(NA_real_, nrow = length(cs_list), ncol = 3L)
  for (i in seq_along(cs_list)) {
    pos <- as.integer(cs_list[[i]])
    if (length(pos) <= 1L) {
      purity_mat[i, ] <- c(1, 1, 1)
    } else {
      pos_eval <- if (length(pos) > max_exact_purity) {
        pos[unique(round(seq.int(1L, length(pos), length.out = 200L)))]
      } else {
        pos
      }
      sub_r <- abs(R[pos_eval, pos_eval, drop = FALSE]) / ld_scale
      vals <- sub_r[upper.tri(sub_r, diag = FALSE)]
      purity_mat[i, ] <- c(
        min(vals),
        sum(vals) / length(vals),
        stats::median(vals)
      )
    }
  }
  purity_df <- as.data.frame(purity_mat)
  colnames(purity_df) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")

  keep <- if (!is.null(min_abs_corr)) {
    which(purity_df[, 1L] >= min_abs_corr)
  } else {
    seq_len(nrow(purity_df))
  }
  if (length(keep) == 0L) {
    return(list(
      cs = NULL,
      coverage = NULL,
      requested_coverage = coverage
    ))
  }

  cs_list <- cs_list[keep]
  purity_df <- purity_df[keep, , drop = FALSE]
  claimed_cov <- claimed_cov[keep]
  effect_indices <- effect_indices[keep]

  row_names <- paste0("L", effect_indices)
  names(cs_list) <- row_names
  rownames(purity_df) <- row_names

  ord <- order(purity_df[, 1L], decreasing = TRUE)
  list(
    cs = cs_list[ord],
    purity = purity_df[ord, , drop = FALSE],
    cs_index = effect_indices[ord],
    coverage = claimed_cov[ord],
    requested_coverage = coverage
  )
}

fit_susie_from_prescaled_xtx <- function(z, XtX, n, L, coverage, max_iter,
                                         refine, prior_weights = NULL,
                                         R_finite = NULL, R_mismatch = "none",
                                         frob_sq = NULL, ld_scale = 1) {
  ns <- asNamespace("susieR")
  apply_pve <- get("apply_pve_adjustment", envir = ns)
  working_fn <- get("summary_stats_working_quantities", envir = ns)
  norm_pw <- get("normalize_prior_weights", envir = ns)
  val_params <- get("validate_and_override_params", envir = ns)
  workhorse <- get("susie_workhorse", envir = ns)
  mix_p_to_z <- get("mixture_reference_p_to_z", envir = ns)

  z_clean <- as.numeric(z)
  z_clean[is.na(z_clean)] <- 0
  p <- length(z_clean)
  if (!is.null(colnames(XtX))) {
    names(z_clean) <- colnames(XtX)
  }
  pve <- apply_pve(z_clean, n, "wald")
  wk <- working_fn(
    z = pve$z,
    n = n,
    pve_adjustment = pve$adjustment,
    scaled_prior_variance = 0.2
  )
  if (!is.null(colnames(XtX))) {
    names(wk$Xty) <- colnames(XtX)
  }

  if (is.null(attr(XtX, "d"))) {
    d_vec <- diag(XtX)
    attr(XtX, "d") <- d_vec
    attr(XtX, "scaled:scale") <- rep(1, p)
  }

  pw <- norm_pw(prior_weights, p)
  conv_method <- if (!is.null(R_finite) || R_mismatch != "none") "pip" else "elbo"

  params_obj <- list(
    L = min(as.integer(L), p),
    scaled_prior_variance = wk$scaled_prior_variance,
    residual_variance = NULL,
    prior_weights = pw,
    null_weight = 0,
    estimate_residual_variance = FALSE,
    estimate_residual_method = "MoM",
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    prior_variance_grid = NULL,
    mixture_weights = NULL,
    unmappable_effects = "none",
    check_null_threshold = 0,
    prior_tol = 1e-9,
    residual_variance_upperbound = Inf,
    model_init = NULL,
    coverage = coverage,
    min_abs_corr = NULL,
    median_abs_corr = NULL,
    compute_univariate_zscore = FALSE,
    max_iter = as.integer(max_iter),
    tol = NULL,
    convergence_method = conv_method,
    verbose = FALSE,
    track_fit = FALSE,
    residual_variance_lowerbound = 0,
    refine = refine,
    n_purity = "auto",
    cs_extension_corr = NULL,
    alpha0 = NULL,
    beta0 = NULL,
    n = wk$n,
    use_NIG = FALSE,
    intercept = FALSE,
    standardize = FALSE,
    check_prior = TRUE,
    slot_prior = NULL,
    L_greedy = NULL,
    greedy_lbf_cutoff = 0.1
  )
  params_obj <- val_params(params_obj)
  params_obj$coverage <- NULL
  params_obj$min_abs_corr <- NULL

  data_obj <- structure(
    list(
      XtX = XtX,
      X = NULL,
      Xty = wk$Xty,
      yty = wk$yty,
      n = wk$n,
      p = p,
      X_colmeans = rep(NA, p),
      y_mean = NA
    ),
    class = "ss"
  )
  ss <- list(data = data_obj, params = params_obj)

  if (!is.null(R_finite) || R_mismatch != "none") {
    B <- if (is.null(R_finite)) Inf else R_finite
    if (is.null(frob_sq) || !is.finite(frob_sq)) {
      if (is.loaded("gwasfm_frob_and_scale_inplace", PACKAGE = "gwasfinemapping")) {
        frob_sq <- .Call("gwasfm_frob_and_scale_inplace", XtX, 1.0,
                         PACKAGE = "gwasfinemapping") / (ld_scale * ld_scale)
      } else {
        frob_sq <- sum((XtX / ld_scale)^2)
      }
    }
    r_db <- if (is.infinite(B)) frob_sq else (B * frob_sq - p^2) / (B + 1)
    eff_r <- p^2 / max(r_db, 1)
    ss$data$R_finite_B <- B
    ss$data$R_finite_diagnostics <- list(
      B = B,
      p = p,
      R_frob_sq_debiased = r_db,
      effective_rank = eff_r,
      r_over_B = eff_r / B,
      Rhat_diag_deviation = abs(diag(XtX) / ld_scale - 1)
    )
    ss$data$R_mismatch <- R_mismatch
    if (R_mismatch != "none") {
      eig <- eigen(XtX, symmetric = TRUE)
      eig$values <- eig$values / ld_scale
      ss$data$eigen_R <- eig
    }
    ss$params$R_mismatch <- R_mismatch
    ss$params$R_mismatch_method <- "mle"
    ss$params$eb_mix_ref <- 5e-8
    ss$data$eb_mix_ref <- 5e-8
    ss$params$eb_mix_z_ref <- mix_p_to_z(5e-8)
    ss$data$eb_mix_z_ref <- ss$params$eb_mix_z_ref
    ss$data$z_marginal <- pve$z
    ss$params$eig_delta_rel <- 1e-3
    ss$params$eig_delta_abs <- 0
    ss$params$artifact_threshold <- 0.1
    ss$params$R_sensitivity_threshold <- 30
  }

  model <- workhorse(ss$data, ss$params)
  if (!is.null(model$R_finite_diagnostics)) {
    model$R_finite_diagnostics$SER_fallback <-
      isTRUE(model$R_finite_diagnostics$R_reliability_flag)
  }
  model
}

run_susie_model <- function(variants, R, method_id, effective_n = NULL,
                            reference_sample_size = NULL, L = 10L,
                            coverage = 0.95, min_abs_corr = 0.5,
                            max_iter = 50L, refine = FALSE,
                            R_mismatch = "eb_mix",
                            psd_max_variants = 5000L,
                            eigen_max_variants = 5000L,
                            prior_weights = NULL,
                            allow_missing_n = FALSE) {
  if (!method_id %in% c("susie_rss", "susie_rrss", "polyfun")) {
    stop("Unsupported SuSiE method identifier.", call. = FALSE)
  }
  if (length(eigen_max_variants) != 1L || !is.numeric(eigen_max_variants) ||
      is.na(eigen_max_variants) || eigen_max_variants < 1L) {
    stop("`eigen_max_variants` must be one positive integer.", call. = FALSE)
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
  n <- resolve_effective_n(variants, effective_n, allow_missing = allow_missing_n)

  r_mismatch_requested <- if (method_id == "susie_rrss") R_mismatch else "none"
  q_art_eigen_status <- "not_applicable"
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
    if (R_mismatch != "none") {
      if (nrow(R) <= as.integer(eigen_max_variants)) {
        q_art_eigen_status <- "computed"
      } else {
        R_mismatch <- "none"
        q_art_eigen_status <- "skipped_large_m"
      }
    }
  } else {
    R_finite <- NULL
    R_mismatch <- "none"
  }

  ld_scale <- attr(R, "ld_scale")
  if (is.null(ld_scale) || !is.finite(ld_scale) || ld_scale <= 0) {
    ld_scale <- 1
  }
  frob_sq_attr <- attr(R, "frob_sq")

  # validate_ld_matrix already verified unit diagonal and symmetry:
  # standardize = FALSE avoids 4 redundant M x M copies, and min_abs_corr = NULL
  # avoids susieR 2.0.0's un-subsampled get_purity() O(M^2) allocation on diffuse slots.
  # When R was pre-scaled in-place to XtX = (n - 1) * R at load time, bypass
  # susie_rss's second M x M allocation cleanly via sufficient_stats_constructor.
  expected_scale <- if (is.null(n$value)) 1 else (n$value - 1)
  if (abs(ld_scale - expected_scale) <= 1e-6 * max(1, expected_scale)) {
    fit <- fit_susie_from_prescaled_xtx(
      z = variants$EZ,
      XtX = R,
      n = n$value,
      L = as.integer(L),
      coverage = coverage,
      max_iter = as.integer(max_iter),
      refine = refine,
      prior_weights = prior_weights,
      R_finite = R_finite,
      R_mismatch = R_mismatch,
      frob_sq = frob_sq_attr,
      ld_scale = ld_scale
    )
  } else {
    R_unscaled <- if (ld_scale == 1) R else R / ld_scale
    fit <- susieR::susie_rss(
      z = variants$EZ,
      R = R_unscaled,
      n = n$value,
      L = as.integer(L),
      prior_weights = prior_weights,
      standardize = FALSE,
      coverage = coverage,
      min_abs_corr = NULL,
      max_iter = as.integer(max_iter),
      refine = refine,
      R_finite = R_finite,
      R_mismatch = R_mismatch
    )
  }
  cs <- compute_credible_sets_fast(
    fit,
    R = R,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    ld_scale = ld_scale
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
    q_art_eigen = q_art_eigen_status,
    r_mismatch_effective = R_mismatch,
    R_reliability_flag = fit$R_finite_diagnostics$R_reliability_flag %||% NA,
    SER_fallback = fit$R_finite_diagnostics$SER_fallback %||% NA,
    artifact_flag = fit$R_finite_diagnostics$artifact_flag %||% NA,
    R_sensitivity_flag = fit$R_sensitivity_flag %||% NA,
    R_finite_diagnostics = fit$R_finite_diagnostics
  )

  lead_idx <- which.max(
    ifelse(is.na(variants$LP), -Inf, as.numeric(variants$LP))
  )
  ld_r2_to_lead <- if ("LD_R2_TO_LEAD" %in% names(variants) &&
                       !anyNA(variants$LD_R2_TO_LEAD)) {
    as.numeric(variants$LD_R2_TO_LEAD)
  } else if (length(lead_idx) == 1L && nrow(R) == nrow(variants)) {
    (as.numeric(R[, lead_idx]) / ld_scale)^2
  } else {
    rep(NA_real_, nrow(variants))
  }
  expected_z <- if (!is.null(fit$XtXr) && length(fit$XtXr) == nrow(variants) &&
                    !is.null(n$value) && n$value > 1) {
    as.numeric(fit$XtXr / sqrt(n$value - 1))
  } else if (!is.null(fit$XtXr) && length(fit$XtXr) == nrow(variants) && is.null(n$value)) {
    as.numeric(fit$XtXr)
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
        R_mismatch_requested = r_mismatch_requested,
        psd_max_variants = as.integer(psd_max_variants),
        eigen_max_variants = as.integer(eigen_max_variants)
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
                          max_iter = 50L, refine = FALSE,
                          psd_max_variants = 5000L,
                          allow_missing_n = FALSE) {
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
    psd_max_variants = psd_max_variants,
    allow_missing_n = allow_missing_n
  )
}

run_susie_rrss <- function(variants, R, reference_sample_size,
                           effective_n = NULL, L = 10L, coverage = 0.95,
                           min_abs_corr = 0.5, max_iter = 50L,
                           refine = FALSE, R_mismatch = "eb_mix",
                           psd_max_variants = 5000L,
                           eigen_max_variants = 5000L,
                           allow_missing_n = FALSE) {
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
    psd_max_variants = psd_max_variants,
    eigen_max_variants = eigen_max_variants,
    allow_missing_n = allow_missing_n
  )
}
