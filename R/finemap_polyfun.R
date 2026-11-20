resolve_polyfun_priors <- function(variants, polyfun_priors = NULL,
                                   prior_weights = NULL,
                                   missing_prior_policy = c("baseline_floor", "drop", "error"),
                                   truncation_q = 100) {
  missing_prior_policy <- match.arg(missing_prior_policy)
  m <- nrow(variants)

  if (!is.null(prior_weights)) {
    if (!is.numeric(prior_weights) || length(prior_weights) != m ||
        anyNA(prior_weights) || any(!is.finite(prior_weights)) ||
        any(prior_weights < 0) || sum(prior_weights) <= 0) {
      stop("`prior_weights` must be a non-negative finite numeric vector of length nrow(variants) with positive sum.",
           call. = FALSE)
    }
    pw <- as.numeric(prior_weights) / sum(prior_weights)
    return(list(
      status = "ok",
      source = "explicit_vector",
      snpvar = as.numeric(prior_weights),
      prior_weights = pw,
      matched = rep(TRUE, m),
      keep_mask = rep(TRUE, m),
      floor_value = min(prior_weights[prior_weights > 0])
    ))
  }

  if ("POLYFUN_SNPVAR" %in% names(variants) && is.null(polyfun_priors)) {
    sv <- suppressWarnings(as.numeric(variants$POLYFUN_SNPVAR))
    matched <- !is.na(sv) & is.finite(sv) & sv > 0
    if (!any(matched)) {
      stop("`POLYFUN_SNPVAR` column contains no positive finite values.", call. = FALSE)
    }
    floor_val <- max(sv[matched]) / truncation_q
    sv[matched & sv < floor_val] <- floor_val
    if (any(!matched)) {
      if (missing_prior_policy == "error") {
        stop(sum(!matched), " variants lack positive POLYFUN_SNPVAR priors.", call. = FALSE)
      } else if (missing_prior_policy == "baseline_floor") {
        sv[!matched] <- floor_val
      }
    }
    keep <- if (missing_prior_policy == "drop") matched else rep(TRUE, m)
    sv_keep <- sv[keep]
    return(list(
      status = "ok",
      source = "variants_column",
      snpvar = sv_keep,
      prior_weights = sv_keep / sum(sv_keep),
      matched = matched[keep],
      keep_mask = keep,
      floor_value = floor_val
    ))
  }

  ld_chrom <- if ("LD_CHROM" %in% names(variants)) {
    sub("^chr", "", as.character(variants$LD_CHROM), ignore.case = TRUE)
  } else {
    sub("^chr", "", as.character(variants$CHROM), ignore.case = TRUE)
  }
  ld_pos <- as.integer(if ("LD_POS" %in% names(variants)) variants$LD_POS else variants$POS)
  ld_ref <- toupper(as.character(if ("LD_REF" %in% names(variants)) variants$LD_REF else variants$REF))
  ld_alt <- toupper(as.character(if ("LD_ALT" %in% names(variants)) variants$LD_ALT else variants$ALT))

  chrom_u <- unique(ld_chrom)
  if (length(chrom_u) != 1L) {
    stop("PolyFun locus variants must belong to a single chromosome.", call. = FALSE)
  }
  chrom <- chrom_u[[1L]]

  ref_root <- Sys.getenv(
    "GWASFM_REFERENCE_ROOT",
    path.expand("~/.cache/gwasfinemapping/reference")
  )
  priors_target <- if (!is.null(polyfun_priors) && nzchar(polyfun_priors)) {
    path.expand(polyfun_priors)
  } else if (dir.exists(file.path(ref_root, "polyfun"))) {
    file.path(ref_root, "polyfun")
  } else {
    path.expand("~/.cache/gwasfinemapping/reference/polyfun")
  }

  # Directory containing installed per-chromosome Parquet slices
  if (dir.exists(priors_target)) {
    if (!chrom %in% as.character(1:22)) {
      return(list(
        status = "skipped_unsupported_chromosome",
        source = "ukb_meta_15traits",
        reason = paste0("PolyFun precomputed priors cover autosomes 1-22 only (got chromosome ", chrom, ").")
      ))
    }
    cpath <- file.path(priors_target, "by_chrom", paste0("chr", chrom, ".parquet"))
    meta_path <- file.path(priors_target, "polyfun_meta.json")
    if (!file.exists(cpath)) {
      stop("PolyFun reference Parquet not found: ", cpath,
           " (run `bash scripts/install_reference_data.sh --only polyfun`).",
           call. = FALSE)
    }
    meta <- if (file.exists(meta_path) && requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::fromJSON(meta_path, simplifyVector = TRUE)
    } else {
      list(min_snpvar_scaled = NA_real_)
    }

    tmp_in <- tempfile(pattern = "gwasfm_pf_in_", fileext = ".tsv")
    tmp_out <- tempfile(pattern = "gwasfm_pf_out_", fileext = ".tsv")
    on.exit(unlink(c(tmp_in, tmp_out), force = TRUE), add = TRUE)
    utils::write.table(
      data.frame(BP = ld_pos, REF = ld_ref, ALT = ld_alt, stringsAsFactors = FALSE),
      file = tmp_in, sep = "\t", quote = FALSE, row.names = FALSE
    )

    py_bin <- Sys.which("python3")
    if (!nzchar(py_bin)) py_bin <- "python3"
    py_code <- paste(
      "import sys, pandas as pd, pyarrow.parquet as pq",
      "cpath, qpath, opath = sys.argv[1:4]",
      "q = pd.read_csv(qpath, sep='\\t', dtype={'BP': int, 'REF': str, 'ALT': str})",
      "pmin, pmax = int(q['BP'].min()), int(q['BP'].max())",
      "ref = pq.read_table(cpath, filters=[('BP', '>=', pmin), ('BP', '<=', pmax)]).to_pandas()",
      "k1 = ref['BP'].astype(str) + ':' + ref['A1'].str.upper() + ':' + ref['A2'].str.upper()",
      "k2 = ref['BP'].astype(str) + ':' + ref['A2'].str.upper() + ':' + ref['A1'].str.upper()",
      "lkp = dict(zip(k1, ref['SNPVAR']))",
      "for k, v in zip(k2, ref['SNPVAR']): lkp.setdefault(k, v)",
      "qk = q['BP'].astype(str) + ':' + q['REF'] + ':' + q['ALT']",
      "sv = qk.map(lkp)",
      "pd.DataFrame({'SNPVAR': sv}).to_csv(opath, sep='\\t', index=False, na_rep='NA')",
      sep = "\n"
    )
    status_py <- system2(py_bin, c("-c", shQuote(py_code), shQuote(cpath), shQuote(tmp_in), shQuote(tmp_out)))
    if (status_py != 0L || !file.exists(tmp_out)) {
      stop("Failed to query PolyFun Parquet priors via python3/pyarrow.", call. = FALSE)
    }
    sv <- suppressWarnings(as.numeric(utils::read.delim(tmp_out, stringsAsFactors = FALSE)$SNPVAR))
    matched <- !is.na(sv) & is.finite(sv) & sv > 0
    floor_val <- as.numeric(meta$min_snpvar_scaled)
    if (is.na(floor_val) || !is.finite(floor_val) || floor_val <= 0) {
      floor_val <- if (any(matched)) max(sv[matched]) / truncation_q else 1.0
    }
    if (!any(matched) && missing_prior_policy != "baseline_floor") {
      stop("Zero variants in locus matched PolyFun reference priors.", call. = FALSE)
    }
    if (any(!matched)) {
      if (missing_prior_policy == "error") {
        stop(sum(!matched), " locus variants were not found in PolyFun snpvar_meta.", call. = FALSE)
      } else if (missing_prior_policy == "baseline_floor") {
        sv[!matched] <- floor_val
      }
    }
    keep <- if (missing_prior_policy == "drop") matched else rep(TRUE, m)
    if (!any(keep)) {
      stop("All locus variants were dropped by missing_prior_policy = 'drop'.", call. = FALSE)
    }
    sv_keep <- sv[keep]
    return(list(
      status = "ok",
      source = "ukb_meta_15traits",
      snpvar = sv_keep,
      prior_weights = sv_keep / sum(sv_keep),
      matched = matched[keep],
      keep_mask = keep,
      floor_value = floor_val
    ))
  }

  # Tab/whitespace-delimited custom PolyFun file (e.g. *.snpvar_ridge_constrained.gz or .tsv)
  if (!file.exists(priors_target)) {
    stop("PolyFun priors path does not exist: ", priors_target, call. = FALSE)
  }
  tbl <- utils::read.table(priors_target, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  req_cols <- c("CHR", "BP", "A1", "A2")
  if (!all(req_cols %in% names(tbl))) {
    stop("Custom PolyFun prior table must include columns CHR, BP, A1, A2, and SNPVAR (or snpvar_bin).",
         call. = FALSE)
  }
  sv_col <- if ("SNPVAR" %in% names(tbl)) "SNPVAR" else if ("snpvar_bin" %in% names(tbl)) "snpvar_bin" else NULL
  if (is.null(sv_col)) {
    stop("Custom PolyFun prior table lacks SNPVAR or snpvar_bin column.", call. = FALSE)
  }
  tbl_chr <- sub("^chr", "", as.character(tbl$CHR), ignore.case = TRUE)
  tbl <- tbl[tbl_chr == chrom, , drop = FALSE]
  if (nrow(tbl) == 0L) {
    return(list(
      status = "skipped_unsupported_chromosome",
      source = "custom_polyfun_table",
      reason = paste0("Custom PolyFun table contains no variants on chromosome ", chrom, ".")
    ))
  }
  raw_sv <- as.numeric(tbl[[sv_col]])
  floor_val <- max(raw_sv, na.rm = TRUE) / truncation_q
  raw_sv[is.na(raw_sv) | raw_sv < floor_val] <- floor_val

  k_exact <- paste(as.integer(tbl$BP), toupper(tbl$A1), toupper(tbl$A2), sep = ":")
  k_swap <- paste(as.integer(tbl$BP), toupper(tbl$A2), toupper(tbl$A1), sep = ":")
  q_keys <- paste(ld_pos, ld_ref, ld_alt, sep = ":")

  idx1 <- match(q_keys, k_exact)
  idx2 <- match(q_keys, k_swap)
  idx <- ifelse(!is.na(idx1), idx1, idx2)
  matched <- !is.na(idx)
  sv <- ifelse(matched, raw_sv[idx], NA_real_)
  if (any(!matched)) {
    if (missing_prior_policy == "error") {
      stop(sum(!matched), " locus variants were not found in custom PolyFun table.", call. = FALSE)
    } else if (missing_prior_policy == "baseline_floor") {
      sv[!matched] <- floor_val
    }
  }
  keep <- if (missing_prior_policy == "drop") matched else rep(TRUE, m)
  if (!any(keep)) {
    stop("All locus variants were dropped by missing_prior_policy = 'drop'.", call. = FALSE)
  }
  sv_keep <- sv[keep]
  list(
    status = "ok",
    source = "custom_polyfun_table",
    snpvar = sv_keep,
    prior_weights = sv_keep / sum(sv_keep),
    matched = matched[keep],
    keep_mask = keep,
    floor_value = floor_val
  )
}

run_polyfun_susie <- function(variants, R, prior_weights = NULL,
                              polyfun_priors = NULL,
                              missing_prior_policy = c("baseline_floor", "drop", "error"),
                              effective_n = NULL, L = 10L,
                              coverage = 0.95, min_abs_corr = 0.5,
                              max_iter = 100L, refine = FALSE,
                              psd_max_variants = 5000L,
                              allow_missing_n = FALSE) {
  missing_prior_policy <- match.arg(missing_prior_policy)
  variants <- validate_gwas_stats(variants, require_p = TRUE, require_z = TRUE)
  if (nrow(variants) == 0L) {
    stop("PolyFun + SuSiE cannot fit a region with no variants.", call. = FALSE)
  }

  priors_info <- resolve_polyfun_priors(
    variants = variants,
    polyfun_priors = polyfun_priors,
    prior_weights = prior_weights,
    missing_prior_policy = missing_prior_policy
  )

  if (identical(priors_info$status, "skipped_unsupported_chromosome")) {
    ids <- variant_key(variants)
    gwas_ids <- resolve_gwas_variant_ids(variants, ids)
    variant_results <- variants
    if (!"GWAS_VARIANT_ID" %in% names(variant_results)) {
      variant_results$GWAS_VARIANT_ID <- gwas_ids
    }
    variant_results$method_id <- "polyfun"
    variant_results$SNPVAR <- NA_real_
    variant_results$PRIOR_WEIGHT <- NA_real_
    variant_results$PRIOR_MATCHED <- FALSE
    variant_results$PIP <- NA_real_
    variant_results$EXPECTED_Z <- NA_real_
    variant_results$credible_set_id <- NA_character_
    return(list(
      variants = variant_results,
      credible_sets = data.frame(
        credible_set_id = character(), effect_index = integer(),
        size = integer(), coverage = numeric(), min_abs_corr = numeric(),
        mean_abs_corr = numeric(), median_abs_corr = numeric(),
        prior_variance_V = numeric(), lbf = numeric(),
        ld_mismatch_warning = logical(),
        stringsAsFactors = FALSE
      ),
      credible_set_membership = data.frame(
        credible_set_id = character(), effect_index = integer(),
        variant_id = character(), gwas_variant_id = character(),
        stringsAsFactors = FALSE
      ),
      diagnostics = list(
        status = "skipped_unsupported_chromosome",
        skip_reason = priors_info$reason,
        converged = NA,
        prior_source = priors_info$source,
        prior_matched_variants = 0L,
        prior_baseline_assigned_variants = 0L,
        prior_match_rate = 0
      ),
      provenance = list(
        method_id = "polyfun",
        method_version = paste0("polyfun_priors+susieR_", utils::packageVersion("susieR")),
        parameters = list(
          L = as.integer(L),
          coverage = coverage,
          min_abs_corr = min_abs_corr,
          max_iter = as.integer(max_iter),
          prior_source = priors_info$source,
          missing_prior_policy = missing_prior_policy
        )
      ),
      fit = NULL
    ))
  }

  keep_mask <- priors_info$keep_mask
  if (!all(keep_mask)) {
    ld_scale_orig <- attr(R, "ld_scale")
    preval_orig <- attr(R, "gwasfm_prevalidated")
    variants <- variants[keep_mask, , drop = FALSE]
    R <- R[keep_mask, keep_mask, drop = FALSE]
    if (!is.null(ld_scale_orig)) {
      attr(R, "ld_scale") <- ld_scale_orig
      d_vec <- rep(ld_scale_orig, nrow(R))
      names(d_vec) <- colnames(R)
      attr(R, "d") <- d_vec
      attr(R, "scaled:scale") <- rep(1, nrow(R))
    }
    if (isTRUE(preval_orig)) {
      attr(R, "gwasfm_prevalidated") <- TRUE
    }
  }

  res <- run_susie_model(
    variants = variants,
    R = R,
    method_id = "polyfun",
    effective_n = effective_n,
    L = L,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter,
    refine = refine,
    psd_max_variants = psd_max_variants,
    prior_weights = priors_info$prior_weights,
    allow_missing_n = allow_missing_n
  )

  n_matched <- sum(priors_info$matched)
  n_total <- length(priors_info$matched)
  n_floor <- n_total - n_matched
  pw <- priors_info$prior_weights
  max_min_ratio <- max(pw) / min(pw[pw > 0])

  res$variants$SNPVAR <- priors_info$snpvar
  res$variants$PRIOR_WEIGHT <- pw
  res$variants$PRIOR_MATCHED <- priors_info$matched

  res$diagnostics$prior_source <- priors_info$source
  res$diagnostics$prior_matched_variants <- as.integer(n_matched)
  res$diagnostics$prior_baseline_assigned_variants <- as.integer(n_floor)
  res$diagnostics$prior_match_rate <- round(n_matched / n_total, 4)
  res$diagnostics$prior_max_min_ratio <- round(max_min_ratio, 4)
  res$provenance$method_version <- paste0(
    "polyfun_priors+susieR_", as.character(utils::packageVersion("susieR"))
  )
  res$provenance$parameters$prior_source <- priors_info$source
  res$provenance$parameters$missing_prior_policy <- missing_prior_policy
  res$provenance$parameters$prior_floor_value <- priors_info$floor_value
  res
}
