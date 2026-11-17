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
