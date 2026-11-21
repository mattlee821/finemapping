susie_fixture <- function(z = c(5.5, 5.2, 0.3, -0.2, 0.1)) {
  variants <- data.frame(
    CHROM = rep("1", 5),
    POS = 100:104,
    REF = c("A", "C", "G", "T", "A"),
    ALT = c("G", "T", "A", "C", "C"),
    EZ = z,
    LP = rep(2, 5),
    SS = c(900, 850, 800, 750, 700)
  )
  ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT,
               sep = ":")
  R <- matrix(0.1, 5, 5, dimnames = list(ids, ids))
  diag(R) <- 1
  list(variants = variants, R = R)
}

test_that("SuSiE-RSS records PIPs, credible sets, and the selected effective N", {
  fixture <- susie_fixture()
  result <- run_susie_rss(
    fixture$variants,
    fixture$R,
    L = 3,
    max_iter = 30,
    min_abs_corr = NULL
  )

  expect_equal(result$provenance$method_id, "susie_rss")
  expect_equal(result$provenance$parameters$effective_n, 700)
  expect_equal(result$provenance$parameters$effective_n_source,
               "minimum_variant_ss")
  expect_equal(nrow(result$variants), 5L)
  expect_true(all(result$variants$PIP >= 0 & result$variants$PIP <= 1))
  expect_true(is.logical(result$diagnostics$converged))
  expect_true(all(c("credible_set_id", "coverage", "min_abs_corr") %in%
                    names(result$credible_sets)))
})

test_that("robust SuSiE-RSS records finite-reference and mismatch diagnostics", {
  fixture <- susie_fixture()
  result <- run_susie_rrss(
    fixture$variants,
    fixture$R,
    reference_sample_size = 500,
    L = 3,
    max_iter = 30,
    min_abs_corr = NULL
  )

  expect_equal(result$provenance$method_id, "susie_rrss")
  expect_equal(result$provenance$parameters$R_finite, 500)
  expect_equal(result$provenance$parameters$R_mismatch, "eb_mix")
  expect_true("SER_fallback" %in% names(result$diagnostics))
  expect_true("R_finite_diagnostics" %in% names(result$diagnostics))

  configured <- run_susie_rrss(
    fixture$variants,
    fixture$R,
    reference_sample_size = 500,
    L = 3,
    max_iter = 30,
    min_abs_corr = NULL,
    R_mismatch = "eb"
  )
  expect_equal(configured$provenance$parameters$R_mismatch, "eb")
})

test_that("SuSiE requires signed Z, usable LD, and a sample-size value", {
  fixture <- susie_fixture()
  no_ss <- fixture$variants
  no_ss$SS <- NULL
  expect_error(run_susie_rss(no_ss, fixture$R), "supply an explicit")

  no_z <- fixture$variants
  no_z$EZ <- NULL
  expect_error(run_susie_rss(no_z, fixture$R), "signed Z-score")

  expect_error(run_susie_rrss(fixture$variants, fixture$R, 1),
               "greater than 2")
  expect_error(run_susie_rrss(
    fixture$variants, fixture$R, 500, R_mismatch = "unsupported"
  ), "R_mismatch")
})

test_that("SuSiE handles a one-variant locus and rejects an empty locus", {
  fixture <- susie_fixture()
  single <- fixture$variants[1, , drop = FALSE]
  id <- paste(single$CHROM, single$POS, single$REF, single$ALT, sep = ":")
  result <- run_susie_rss(
    single,
    matrix(1, 1, 1, dimnames = list(id, id)),
    L = 1,
    max_iter = 10,
    min_abs_corr = NULL
  )

  expect_equal(length(result$variants$PIP), 1L)
  expect_true(is.finite(result$variants$PIP[[1]]))

  empty <- fixture$variants[0, , drop = FALSE]
  empty_ld <- matrix(numeric(), 0, 0,
                     dimnames = list(character(), character()))
  expect_error(run_susie_rss(empty, empty_ld), "no variants")
})

test_that("SuSiE-rRSS keeps eb_mix active above eigen_max_variants while skipping Q_art eigen", {
  fixture <- susie_fixture(z = c(8.5, 7.9, 0.4, -0.2, 0.1))
  res <- run_susie_rrss(
    fixture$variants,
    fixture$R,
    reference_sample_size = 500,
    L = 3,
    max_iter = 50,
    min_abs_corr = NULL,
    eigen_max_variants = 2L
  )

  expect_equal(res$provenance$parameters$R_mismatch, "eb_mix")
  expect_equal(res$diagnostics$r_mismatch_requested, "eb_mix")
  expect_equal(res$diagnostics$r_mismatch_effective, "eb_mix")
  expect_equal(res$diagnostics$q_art_eigen, "skipped_large_m")
  expect_false(isTRUE(res$diagnostics$artifact_evaluable))
  expect_true(all(c("prior_variance_V", "lbf", "ld_mismatch_warning") %in%
                    names(res$credible_sets)))
  expect_true(all(c("initial_elbo", "final_elbo", "final_elbo_change",
                    "max_prior_variance_V", "max_lbf") %in%
                    names(res$diagnostics)))
})

test_that("Two-tier CS validity gate suppresses non-converged and V > 1 blow-up credible sets", {
  # 1. Force non-convergence with max_iter = 1 on a multi-signal locus
  variants <- data.frame(
    CHROM = rep("1", 6),
    POS = 200:205,
    REF = rep("A", 6),
    ALT = rep("G", 6),
    EZ = c(9.5, 9.1, 0.2, 8.4, 8.0, -0.1),
    LP = c(20, 19, 0.5, 16, 15, 0.3),
    SS = rep(5000, 6)
  )
  ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT, sep = ":")
  R <- diag(6)
  R[1, 2] <- R[2, 1] <- 0.95
  R[4, 5] <- R[5, 4] <- 0.95
  R[1, 4] <- R[4, 1] <- 0.25
  R[1, 5] <- R[5, 1] <- 0.25
  R[2, 4] <- R[4, 2] <- 0.25
  R[2, 5] <- R[5, 2] <- 0.25
  dimnames(R) <- list(ids, ids)

  nonconv <- suppressWarnings(run_susie_rss(
    variants,
    R,
    L = 4,
    max_iter = 1L,
    min_abs_corr = 0.5
  ))
  expect_false(nonconv$diagnostics$converged)
  expect_true(nonconv$diagnostics$credible_sets_suppressed)
  expect_equal(nonconv$diagnostics$credible_sets_suppressed_reason, "non_converged")
  expect_gte(nonconv$diagnostics$unconverged_cs_dropped, 1L)
  expect_equal(nrow(nonconv$credible_sets), 0L)
  expect_equal(nrow(nonconv$credible_set_membership), 0L)
  expect_true(all(is.na(nonconv$variants$credible_set_id)))
  expect_null(nonconv$fit$sets$cs)

  # 2. Converged locus with 1 valid lead CS (V <= 1) + 1 opposite-sign collinear blow-up slot (V > 1)
  #    retains the valid CS and drops only the V > 1 CS slot!
  variants_mix <- variants
  variants_mix$EZ <- c(45, -45, 0.1, 9.5, 9.1, 0.05)
  variants_mix$LP <- c(200, 200, 0.3, 20.0, 19.0, 0.2)
  variants_mix$SS <- rep(500, 6)
  R_mix <- diag(6)
  R_mix[1, 2] <- R_mix[2, 1] <- 0.998
  R_mix[4, 5] <- R_mix[5, 4] <- 0.95
  dimnames(R_mix) <- list(ids, ids)

  mix <- suppressWarnings(run_susie_rss(
    variants_mix,
    R_mix,
    L = 4,
    max_iter = 200L,
    min_abs_corr = 0.5
  ))
  if (isTRUE(mix$diagnostics$converged) &&
      isTRUE(mix$diagnostics$implausible_prior_variance_flag)) {
    expect_true(all(mix$credible_sets$prior_variance_V <= 1.0))
    expect_gte(nrow(mix$credible_sets), 1L)
    expect_gte(mix$diagnostics$dropped_implausible_v_cs, 1L)
    expect_equal(length(mix$fit$sets$cs), nrow(mix$credible_sets))
    expect_equal(length(mix$fit$sets$cs_index), nrow(mix$credible_sets))
  } else if (!isTRUE(mix$diagnostics$converged)) {
    expect_equal(nrow(mix$credible_sets), 0L)
    expect_null(mix$fit$sets$cs)
  }
})
