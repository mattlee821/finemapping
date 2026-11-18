carma_credible_sets_fast <- function(pip, R, rho = 0.95, min_abs_corr = 0.5,
                                     ld_scale = 1) {
  start_r <- max(0.5, as.numeric(min_abs_corr %||% 0.5))
  grid_coarse <- seq(from = 0.5, to = 0.95, by = 0.05)
  grid_fine <- seq(from = 0.96, to = 0.99, by = 0.01)
  candidate_r <- c(grid_coarse, grid_fine)
  candidate_r <- candidate_r[candidate_r >= start_r & candidate_r <= 0.99]
  if (length(candidate_r) == 0L) {
    candidate_r <- start_r
  }

  p <- length(pip)
  pip_clean <- ifelse(is.na(pip), 0, as.numeric(pip))
  if (sum(pip_clean) <= rho) {
    return(list())
  }

  snp_list <- vector("list", length(candidate_r))
  base_order <- order(pip_clean, decreasing = TRUE)

  for (r_idx in seq_along(candidate_r)) {
    cor_thr <- candidate_r[[r_idx]] * ld_scale
    avail <- rep(TRUE, p)
    pip_order <- base_order
    rem_sum <- sum(pip_clean)
    cs_r <- list()
    s <- 1L
    while (s <= length(pip_order) && rem_sum > rho) {
      lead_i <- pip_order[[s]]
      cor_group <- which(avail & abs(R[, lead_i]) > cor_thr)
      if (sum(pip_clean[cor_group]) > rho) {
        g_pip <- pip_clean[cor_group]
        ord_g <- order(g_pip, decreasing = TRUE)
        sorted_pip <- g_pip[ord_g]
        k_take <- which(cumsum(sorted_pip) > rho)[[1L]]
        snp_index <- cor_group[ord_g[seq_len(k_take)]]
        cs_r[[length(cs_r) + 1L]] <- as.integer(snp_index)
        avail[snp_index] <- FALSE
        rem_Before <- sum(pip_clean[pip_order[s:length(pip_order)]])
        pip_order <- pip_order[!pip_order %in% snp_index]
        rem_sum <- if (s <= length(pip_order)) {
          sum(pip_clean[pip_order[s:length(pip_order)]])
        } else {
          0
        }
      } else {
        rem_sum <- rem_sum - pip_clean[[lead_i]]
        s <- s + 1L
      }
    }
    snp_list[[r_idx]] <- cs_r
  }

  counts <- lengths(snp_list)
  if (max(counts) == 0L) {
    return(list())
  }
  best_r_idx <- max(which(counts == max(counts)))
  snp_list[[best_r_idx]]
}

carma_model_keys <- function(model_mat) {
  if (nrow(model_mat) == 0L) return(character())
  k <- ncol(model_mat)
  if (k == 1L) return(as.character(model_mat[, 1L]))
  cols <- lapply(seq_len(k), function(j) model_mat[, j])
  do.call(paste, c(cols, list(sep = ",")))
}

carma_eval_spike_slab_batch <- function(model_mat, R, z, tau, ld_scale = 1) {
  n_mod <- nrow(model_mat)
  if (n_mod == 0L) return(numeric())
  k <- ncol(model_mat)

  if (is.loaded("gwasfm_carma_spike_slab_batch", PACKAGE = "gwasfinemapping")) {
    storage.mode(model_mat) <- "integer"
    return(.Call(
      "gwasfm_carma_spike_slab_batch",
      model_mat, R, as.numeric(z), as.numeric(tau), as.numeric(ld_scale),
      PACKAGE = "gwasfinemapping"
    ))
  }

  if (k == 1L) {
    idx <- as.integer(model_mat[, 1L])
    zi <- z[idx]
    return(0.5 * log(tau) - 0.5 * log(1 + tau) + 0.5 * (zi * zi) / (1 + tau))
  }

  ns_carma <- asNamespace("CARMA")
  ind_fn <- get("ind_Normal_fixed_sigma_marginal", envir = ns_carma)
  z_mat <- matrix(as.numeric(z), ncol = 1L)
  out <- numeric(n_mod)
  seq_k <- seq_len(k)
  for (m in seq_len(n_mod)) {
    idx <- as.integer(model_mat[m, ])
    sub_r <- R[idx, idx, drop = FALSE] / ld_scale
    sub_z <- z_mat[idx, , drop = FALSE]
    out[[m]] <- ind_fn(seq_k, sub_r, sub_z, tau, k, 1)
  }
  out
}

carma_eval_cauchy_batch <- function(model_mat, R, z, tau_sample, y_var = 1,
                                    ld_scale = 1) {
  n_mod <- nrow(model_mat)
  if (n_mod == 0L) return(numeric())
  k <- ncol(model_mat)
  ns_carma <- asNamespace("CARMA")
  cauchy_fn <- get("Cauchy_fixed_sigma_marginal", envir = ns_carma)
  z_mat <- matrix(as.numeric(z), ncol = 1L)
  out <- numeric(n_mod)
  seq_k <- seq_len(k)
  for (m in seq_len(n_mod)) {
    idx <- as.integer(model_mat[m, ])
    sub_r <- R[idx, idx, drop = FALSE] / ld_scale
    sub_z <- z_mat[idx, , drop = FALSE]
    out[[m]] <- cauchy_fn(seq_k, sub_r, sub_z, tau_sample, k, y_var)
  }
  out
}

carma_build_neighborhood <- function(p, input_S, conditional_S = integer()) {
  working_S <- if (length(conditional_S) > 0L) {
    setdiff(input_S, conditional_S)
  } else {
    input_S
  }
  excluded <- unique(c(working_S, conditional_S))
  S_sub <- if (length(excluded) > 0L) seq_len(p)[-excluded] else seq_len(p)
  k <- length(working_S)
  n_sub <- length(S_sub)

  if (k == 0L) {
    return(list(
      working_S = integer(),
      gamma_minus = matrix(integer(), nrow = 0L, ncol = 0L),
      gamma_plus = matrix(as.integer(S_sub), ncol = 1L),
      gamma_equal = matrix(integer(), nrow = 0L, ncol = 0L)
    ))
  }

  if (k == 1L) {
    s1 <- working_S[[1L]]
    # Vectorized sorted 2-column matrix for (s1, S_sub)
    g_plus <- cbind(pmin.int(s1, S_sub), pmax.int(s1, S_sub))
    storage.mode(g_plus) <- "integer"
    return(list(
      working_S = as.integer(working_S),
      gamma_minus = matrix(integer(), nrow = 1L, ncol = 0L),
      gamma_plus = g_plus,
      gamma_equal = matrix(as.integer(S_sub), ncol = 1L)
    ))
  }

  # k >= 2: vectorized row-sorted construction without apply() closures
  ws_sorted <- sort(as.integer(working_S))
  # gamma_minus: k rows, each dropping one element of ws_sorted (already sorted!)
  g_minus <- matrix(0L, nrow = k, ncol = k - 1L)
  for (r in seq_len(k)) {
    g_minus[r, ] <- ws_sorted[-r]
  }

  # Insert sorted vector `v` into sorted vector `base` of length m in O(m * length(v))
  insert_sorted_batch <- function(base_vec, cand_vec) {
    m_len <- length(base_vec)
    n_c <- length(cand_vec)
    if (n_c == 0L) return(matrix(integer(), nrow = 0L, ncol = m_len + 1L))
    pos <- findInterval(cand_vec, base_vec) # 0..m_len
    mat <- matrix(0L, nrow = n_c, ncol = m_len + 1L)
    for (col_j in seq_len(m_len + 1L)) {
      before <- pos >= col_j
      equal_pos <- pos == (col_j - 1L)
      after <- pos < (col_j - 1L)
      if (col_j <= m_len && any(before)) {
        mat[before, col_j] <- base_vec[[col_j]]
      }
      if (any(equal_pos)) {
        mat[equal_pos, col_j] <- cand_vec[equal_pos]
      }
      if (col_j >= 2L && any(after)) {
        mat[after, col_j] <- base_vec[[col_j - 1L]]
      }
    }
    mat
  }

  g_plus <- insert_sorted_batch(ws_sorted, S_sub)
  eq_chunks <- lapply(seq_len(k), function(r) {
    insert_sorted_batch(g_minus[r, ], S_sub)
  })
  g_equal <- do.call(rbind, eq_chunks)

  list(
    working_S = ws_sorted,
    gamma_minus = g_minus,
    gamma_plus = g_plus,
    gamma_equal = g_equal
  )
}

carma_shotgun_module <- function(z, R, ld_scale = 1, lambda = 1,
                                 effect_size_prior = c("Spike-slab", "Cauchy"),
                                 tau = 0.04, y_var = 1,
                                 num_causal = 10L, max_model_dim = 200000L,
                                 inner_all_iter = 10L, epsilon = 1e-3,
                                 outlier_switch = TRUE,
                                 outlier_bf_index = 1 / 3.2,
                                 visited_init = NULL,
                                 conditional_S_init = integer(),
                                 seed = 1L) {
  effect_size_prior <- match.arg(effect_size_prior)
  p <- length(z)
  max_k <- max(1L, as.integer(num_causal))
  ns_carma <- asNamespace("CARMA")
  outlier_fn <- if (effect_size_prior == "Spike-slab") {
    get("outlier_ind_Normal_marginal", envir = ns_carma)
  } else {
    get("outlier_Cauchy_fixed_sigma_marginal", envir = ns_carma)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))
  tau_sample <- if (effect_size_prior == "Cauchy") {
    stats::rgamma(1e5L, 0.5, rate = 0.5)
  } else {
    tau
  }
  outlier_tau <- tau_sample

  poisson_log_prior <- function(dim_model) {
    dim_model * log(lambda) + lgamma(p - dim_model + 1) - lgamma(p + 1)
  }
  null_margin <- poisson_log_prior(0)
  running_best_lp <- null_margin

  # Per-dimension visited model store (avoids millions of R list/string allocations)
  visited_by_k <- if (!is.null(visited_init)) {
    visited_init
  } else {
    lapply(seq_len(max_k), function(k) {
      list(mat = matrix(integer(), nrow = 0L, ncol = k), lp = numeric())
    })
  }

  # Exact 1-SNP log-posterior table (computed once in O(p))
  k1_mat <- matrix(seq_len(p), ncol = 1L)
  k1_ml <- if (effect_size_prior == "Spike-slab") {
    carma_eval_spike_slab_batch(k1_mat, R, z, tau = tau, ld_scale = ld_scale)
  } else {
    carma_eval_cauchy_batch(k1_mat, R, z, tau_sample = tau_sample,
                            y_var = y_var, ld_scale = ld_scale)
  }
  k1_lp <- k1_ml + poisson_log_prior(1L)
  running_best_lp <- max(running_best_lp, max(k1_lp, na.rm = TRUE))

  record_visited <- function(mat, log_post) {
    n_r <- nrow(mat)
    k <- ncol(mat)
    if (n_r == 0L || k == 0L || k > length(visited_by_k)) return(invisible(NULL))
    max_cur <- max(log_post, na.rm = TRUE)
    if (is.finite(max_cur) && max_cur > running_best_lp) {
      running_best_lp <<- max_cur
    }
    # Retain all models within 40 log-units of best (exp(-40) = 4.2e-18 < machine eps)
    # plus up to top 256 models from this batch for q25 convergence stability.
    keep <- which(is.finite(log_post) & log_post >= (running_best_lp - 40))
    if (length(keep) < min(256L, n_r)) {
      ord_b <- order(log_post, decreasing = TRUE)
      keep <- unique(c(keep, ord_b[seq_len(min(256L, n_r))]))
    }
    if (length(keep) == 0L) return(invisible(NULL))
    sub_m <- mat[keep, , drop = FALSE]
    sub_lp <- log_post[keep]
    cur <- visited_by_k[[k]]
    comb_m <- rbind(cur$mat, sub_m)
    comb_lp <- c(cur$lp, sub_lp)
    dup <- duplicated(comb_m)
    if (any(dup)) {
      comb_m <- comb_m[!dup, , drop = FALSE]
      comb_lp <- comb_lp[!dup]
    }
    if (length(comb_lp) > 25000L) {
      ord_c <- order(comb_lp, decreasing = TRUE)[seq_len(25000L)]
      comb_m <- comb_m[ord_c, , drop = FALSE]
      comb_lp <- comb_lp[ord_c]
    }
    visited_by_k[[k]] <<- list(mat = comb_m, lp = comb_lp)
    invisible(NULL)
  }

  eval_candidates <- function(mat, record = TRUE) {
    if (nrow(mat) == 0L) return(numeric())
    k <- ncol(mat)
    if (k == 0L) return(rep(null_margin, nrow(mat)))
    log_post <- if (k == 1L) {
      k1_lp[mat[, 1L]]
    } else {
      ml <- if (effect_size_prior == "Spike-slab") {
        carma_eval_spike_slab_batch(mat, R, z, tau = tau, ld_scale = ld_scale)
      } else {
        carma_eval_cauchy_batch(mat, R, z, tau_sample = tau_sample,
                                y_var = y_var, ld_scale = ld_scale)
      }
      ml + poisson_log_prior(k)
    }
    if (record) {
      record_visited(mat, log_post)
    }
    log_post
  }

  sample_softmax <- function(scores) {
    aa <- scores - max(scores, na.rm = TRUE)
    aa[is.nan(aa) | is.na(aa)] <- min(aa[is.finite(aa)], 0)
    probs <- exp(aa)
    if (length(probs) == 1L) return(1L)
    sample.int(length(probs), size = 1L, prob = probs)
  }

  conditional_S <- unique(as.integer(conditional_S_init))
  S <- conditional_S
  stored_bf <- 0
  z_col <- matrix(as.numeric(z), ncol = 1L)

  collect_top_models <- function() {
    lps <- null_margin
    models <- list(integer())
    for (k in seq_along(visited_by_k)) {
      vk <- visited_by_k[[k]]
      if (length(vk$lp) == 0L) next
      mat_k <- vk$mat
      lp_k <- vk$lp
      if (length(conditional_S) > 0L) {
        bad <- matrix(mat_k %in% conditional_S, nrow = nrow(mat_k), ncol = k)
        keep <- rowSums(bad) == 0L
        mat_k <- mat_k[keep, , drop = FALSE]
        lp_k <- lp_k[keep]
      }
      if (length(lp_k) == 0L) next
      lps <- c(lps, lp_k)
      row_list <- split(mat_k, row(mat_k))
      names(row_list) <- NULL
      models <- c(models, row_list)
    }
    ord <- order(lps, decreasing = TRUE)
    if (length(ord) > max_model_dim) {
      ord <- ord[seq_len(max_model_dim)]
    }
    list(log_post = lps[ord], models = models[ord])
  }

  for (l in seq_len(inner_all_iter)) {
    for (h in seq_len(10L)) {
      nb <- carma_build_neighborhood(p, S, conditional_S)
      working_S <- nb$working_S
      k_ws <- length(working_S)

      current_log_margin <- if (k_ws == 0L) {
        null_margin
      } else {
        eval_candidates(matrix(working_S, nrow = 1L), record = FALSE)[[1L]]
      }

      if (k_ws == 0L) {
        g_plus <- nb$gamma_plus
        m_plus <- eval_candidates(g_plus)
        aa <- m_plus - current_log_margin
        aa <- aa - max(aa, na.rm = TRUE)
        aa[is.nan(aa) | is.na(aa)] <- min(aa[is.finite(aa)], 0)
        top_n <- min(length(aa), max(1L, floor(p / 2)))
        ord_top <- order(exp(aa), decreasing = TRUE)[seq_len(top_n)]
        pick <- if (length(ord_top) == 1L) {
          ord_top[[1L]]
        } else {
          ord_top[[sample.int(length(ord_top), 1L, prob = exp(aa[ord_top]))]]
        }
        S <- unique(c(as.integer(g_plus[pick, ]), conditional_S))
      } else {
        sets <- list(nb$gamma_minus, nb$gamma_plus, nb$gamma_equal)
        margins <- lapply(sets, eval_candidates)

        star_idx <- rep(NA_integer_, 3L)
        star_margin <- rep(-Inf, 3L)

        if (length(margins[[1L]]) > 0L) {
          idx1 <- sample_softmax(margins[[1L]] - current_log_margin)
          star_idx[[1L]] <- idx1
          star_margin[[1L]] <- margins[[1L]][[idx1]]
        }

        active_sets <- if (k_ws >= num_causal) c(3L) else c(2L, 3L)
        for (i_set in active_sets) {
          cand_mat <- sets[[i_set]]
          cand_mar <- margins[[i_set]]
          while (nrow(cand_mat) > 0L) {
            pick_i <- sample_softmax(cand_mar - current_log_margin)
            test_S <- as.integer(cand_mat[pick_i, ])
            if (outlier_switch && length(test_S) > 1L) {
              k_t <- length(test_S)
              seq_kt <- seq_len(k_t)
              modi_ld <- R[test_S, test_S, drop = FALSE] / ld_scale
              z_sub <- z_col[test_S, , drop = FALSE]
              ridge_fn <- function(xi) {
                tmp_ld <- xi * modi_ld + (1 - xi) * diag(k_t)
                outlier_fn(seq_kt, tmp_ld, z_sub, outlier_tau, k_t, 1)
              }
              opizer <- stats::optimize(ridge_fn, interval = c(0, 1), maximum = TRUE)
              opt_ld <- opizer$maximum * modi_ld + (1 - opizer$maximum) * diag(k_t)
              raw_bf <- outlier_fn(seq_kt, modi_ld, z_sub, outlier_tau, k_t, 1) -
                        outlier_fn(seq_kt, opt_ld, z_sub, outlier_tau, k_t, 1)
              test_log_bf <- -abs(raw_bf)
              if (is.finite(test_log_bf) && exp(test_log_bf) < outlier_bf_index) {
                cand_mat <- cand_mat[-pick_i, , drop = FALSE]
                cand_mar <- cand_mar[-pick_i]
                new_out <- test_S[is.na(match(test_S, working_S))]
                conditional_S <- unique(c(conditional_S, new_out))
                next
              }
            }
            star_idx[[i_set]] <- pick_i
            star_margin[[i_set]] <- cand_mar[[pick_i]]
            sets[[i_set]] <- cand_mat
            margins[[i_set]] <- cand_mar
            break
          }
        }

        valid_choices <- which(is.finite(star_margin))
        if (k_ws >= num_causal) {
          valid_choices <- setdiff(valid_choices, 2L)
        }
        if (length(valid_choices) > 0L) {
          chosen_i <- if (length(valid_choices) == 1L) {
            valid_choices[[1L]]
          } else {
            valid_choices[[sample_softmax(star_margin[valid_choices] - current_log_margin)]]
          }
          chosen_row <- sets[[chosen_i]][star_idx[[chosen_i]], ]
          S <- unique(c(as.integer(chosen_row), conditional_S))
        }
      }
    }

    top_res <- collect_top_models()
    q25 <- max(1L, round(0.25 * length(top_res$log_post)))
    cur_bf <- mean(top_res$log_post[seq_len(q25)])
    diff_bf <- abs(cur_bf - stored_bf)
    if (diff_bf < epsilon) {
      break
    }
    stored_bf <- cur_bf
  }

  top_res <- collect_top_models()
  lps <- top_res$log_post
  models <- top_res$models

  if (num_causal == 1L) {
    aa <- k1_lp - max(k1_lp, na.rm = TRUE)
    pip <- exp(aa) / sum(exp(aa))
  } else if (length(lps) == 0L) {
    pip <- rep(0, p)
  } else {
    aa <- lps - max(lps, na.rm = TRUE)
    w <- exp(aa)
    w <- w / sum(w)
    pip <- numeric(p)
    lens <- lengths(models)
    non_empty <- which(lens > 0L)
    if (length(non_empty) > 0L) {
      all_snps <- unlist(models[non_empty], use.names = FALSE)
      all_w <- rep.int(w[non_empty], lens[non_empty])
      agg <- rowsum(all_w, all_snps, reorder = FALSE)
      pip[as.integer(rownames(agg))] <- as.numeric(agg[, 1L])
    }
  }

  list(
    pip = pip,
    log_post = lps,
    models = models,
    outliers = conditional_S,
    visited_by_k = visited_by_k,
    inner_iterations = l
  )
}

run_carma <- function(variants, R, lambda = 1,
                      effect_size_prior = c("Spike-slab", "Cauchy"),
                      tau = 0.04, L = 10L, coverage = 0.95,
                      min_abs_corr = 0.5, outlier_switch = TRUE,
                      outlier_bf_index = 1 / 3.2,
                      max_model_dim = 200000L, all_iter = 3L,
                      all_inner_iter = 10L, epsilon_threshold = 1e-5,
                      psd_max_variants = 5000L, seed = 1L) {
  if (!requireNamespace("CARMA", quietly = TRUE)) {
    stop("Package 'CARMA' (ZikunY/CARMA) is required for run_carma().",
         call. = FALSE)
  }
  effect_size_prior <- match.arg(effect_size_prior)
  variants <- validate_gwas_stats(variants, require_p = TRUE, require_z = TRUE)
  if (nrow(variants) == 0L) {
    stop("CARMA cannot fit a region with no variants.", call. = FALSE)
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
  ld_scale <- attr(R, "ld_scale")
  if (is.null(ld_scale) || !is.finite(ld_scale) || ld_scale <= 0) {
    ld_scale <- 1
  }

  p <- nrow(variants)
  z <- as.numeric(variants$EZ)
  z[is.na(z)] <- 0
  eps_locus <- epsilon_threshold * p

  # Burn-in Shotgun step (matches CARMA.R lines 922-934)
  mod_res <- carma_shotgun_module(
    z = z,
    R = R,
    ld_scale = ld_scale,
    lambda = lambda,
    effect_size_prior = effect_size_prior,
    tau = tau,
    num_causal = as.integer(L),
    max_model_dim = as.integer(max_model_dim),
    inner_all_iter = as.integer(all_inner_iter),
    epsilon = eps_locus,
    outlier_switch = outlier_switch,
    outlier_bf_index = outlier_bf_index,
    visited_init = NULL,
    conditional_S_init = integer(),
    seed = seed
  )

  converged <- FALSE
  outer_used <- 0L
  for (g in seq_len(as.integer(all_iter))) {
    outer_used <- g
    q25_prev <- max(1L, round(0.25 * length(mod_res$log_post)))
    prev_bf <- mean(mod_res$log_post[seq_len(q25_prev)])

    mod_res <- carma_shotgun_module(
      z = z,
      R = R,
      ld_scale = ld_scale,
      lambda = lambda,
      effect_size_prior = effect_size_prior,
      tau = tau,
      num_causal = as.integer(L),
      max_model_dim = as.integer(max_model_dim),
      inner_all_iter = as.integer(all_inner_iter),
      epsilon = eps_locus,
      outlier_switch = outlier_switch,
      outlier_bf_index = outlier_bf_index,
      visited_init = mod_res$visited_by_k,
      conditional_S_init = mod_res$outliers,
      seed = if (is.null(seed)) NULL else as.integer(seed) + g
    )

    q25_cur <- max(1L, round(0.25 * length(mod_res$log_post)))
    cur_bf <- mean(mod_res$log_post[seq_len(q25_cur)])
    if (abs(prev_bf - cur_bf) < eps_locus) {
      converged <- TRUE
      break
    }
  }
  if (all_iter <= 1L) converged <- TRUE

  pip <- as.numeric(mod_res$pip)
  cs_raw <- carma_credible_sets_fast(
    pip = pip,
    R = R,
    rho = coverage,
    min_abs_corr = min_abs_corr,
    ld_scale = ld_scale
  )

  if (length(cs_raw) == 0L) {
    credible_sets <- data.frame(
      credible_set_id = character(),
      effect_index = integer(),
      size = integer(),
      coverage = numeric(),
      min_abs_corr = numeric(),
      mean_abs_corr = numeric(),
      median_abs_corr = numeric(),
      stringsAsFactors = FALSE
    )
    membership <- data.frame(
      credible_set_id = character(),
      effect_index = integer(),
      variant_id = character(),
      gwas_variant_id = character(),
      stringsAsFactors = FALSE
    )
    cs_lookup <- rep(NA_character_, p)
  } else {
    cs_ids <- paste0("L", seq_along(cs_raw))
    purity_mat <- matrix(NA_real_, nrow = length(cs_raw), ncol = 3L)
    cov_vec <- numeric(length(cs_raw))
    cs_lookup <- rep(NA_character_, p)
    for (i in seq_along(cs_raw)) {
      pos <- as.integer(cs_raw[[i]])
      cov_vec[[i]] <- sum(pip[pos])
      cs_lookup[pos] <- cs_ids[[i]]
      if (length(pos) <= 1L) {
        purity_mat[i, ] <- c(1, 1, 1)
      } else {
        pos_eval <- if (length(pos) > 500L) {
          pos[unique(round(seq.int(1L, length(pos), length.out = 200L)))]
        } else {
          pos
        }
        sub_r <- abs(R[pos_eval, pos_eval, drop = FALSE]) / ld_scale
        vals <- sub_r[upper.tri(sub_r, diag = FALSE)]
        purity_mat[i, ] <- c(min(vals), sum(vals) / length(vals), stats::median(vals))
      }
    }
    credible_sets <- data.frame(
      credible_set_id = cs_ids,
      effect_index = seq_along(cs_raw),
      size = as.integer(lengths(cs_raw)),
      coverage = cov_vec,
      min_abs_corr = purity_mat[, 1L],
      mean_abs_corr = purity_mat[, 2L],
      median_abs_corr = purity_mat[, 3L],
      stringsAsFactors = FALSE
    )
    membership <- do.call(rbind, lapply(seq_along(cs_raw), function(i) {
      members <- as.integer(cs_raw[[i]])
      data.frame(
        credible_set_id = cs_ids[[i]],
        effect_index = as.integer(i),
        variant_id = ids[members],
        gwas_variant_id = gwas_ids[members],
        stringsAsFactors = FALSE
      )
    }))
  }

  lead_idx <- which.max(ifelse(is.na(variants$LP), -Inf, as.numeric(variants$LP)))
  ld_r2_to_lead <- if ("LD_R2_TO_LEAD" %in% names(variants) &&
                       !anyNA(variants$LD_R2_TO_LEAD)) {
    as.numeric(variants$LD_R2_TO_LEAD)
  } else if (length(lead_idx) == 1L && nrow(R) == p) {
    (as.numeric(R[, lead_idx]) / ld_scale)^2
  } else {
    rep(NA_real_, p)
  }

  outlier_flag <- rep(FALSE, p)
  if (length(mod_res$outliers) > 0L) {
    outlier_flag[mod_res$outliers] <- TRUE
  }

  variant_results <- variants
  if (!"GWAS_VARIANT_ID" %in% names(variant_results)) {
    variant_results$GWAS_VARIANT_ID <- gwas_ids
  }
  if (!"LD_R2_TO_LEAD" %in% names(variant_results)) {
    variant_results$LD_R2_TO_LEAD <- ld_r2_to_lead
  }
  variant_results$method_id <- "carma"
  variant_results$PIP <- pip
  variant_results$CARMA_OUTLIER <- outlier_flag
  variant_results$credible_set_id <- cs_lookup

  diagnostics <- list(
    converged = converged,
    convergence_status = converged,
    outer_iterations = outer_used,
    visited_models = length(mod_res$log_post),
    psd_check = psd_status,
    outlier_switch = outlier_switch,
    outlier_bf_index = outlier_bf_index,
    outlier_count = as.integer(length(mod_res$outliers)),
    outlier_gwas_variant_ids = if (length(mod_res$outliers) > 0L) {
      gwas_ids[mod_res$outliers]
    } else {
      character()
    }
  )

  list(
    variants = variant_results,
    credible_sets = credible_sets,
    credible_set_membership = membership,
    diagnostics = diagnostics,
    provenance = list(
      method_id = "carma",
      method_version = as.character(utils::packageVersion("CARMA")),
      parameters = list(
        lambda = lambda,
        effect_size_prior = effect_size_prior,
        tau = tau,
        L = as.integer(L),
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        outlier_switch = outlier_switch,
        outlier_bf_index = outlier_bf_index,
        max_model_dim = as.integer(max_model_dim),
        all_iter = as.integer(all_iter),
        all_inner_iter = as.integer(all_inner_iter),
        epsilon_threshold = epsilon_threshold,
        seed = seed
      )
    )
  )
}
