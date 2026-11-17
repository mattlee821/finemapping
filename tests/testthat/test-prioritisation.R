test_that("P-value prioritization is stable, configurable, and non-Bayesian", {
  variants <- data.frame(
    CHROM = c("1", "1", "1"),
    POS = c(300, 100, 200),
    REF = c("G", "A", "C"),
    ALT = c("A", "G", "T"),
    ID = c("rs3", "rs1", "rs2"),
    LP = c(5, 5, 3)
  )

  result <- prioritize_p(variants, top_n = 2)
  expect_equal(result$ID, c("rs1", "rs3", "rs2"))
  expect_equal(result$rank, 1:3)
  expect_equal(result$selected, c(TRUE, TRUE, FALSE))
  expect_false(any(c("PIP", "credible_set_id") %in% names(result)))
})

test_that("p-LD creates separate greedy selections for each r-squared threshold", {
  variants <- data.frame(
    CHROM = rep("1", 3),
    POS = 100:102,
    REF = c("A", "C", "G"),
    ALT = c("G", "T", "A"),
    LP = c(8, 7, 6)
  )
  ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT,
               sep = ":")
  R <- matrix(c(
    1, 0.2, 0.9,
    0.2, 1, 0.5,
    0.9, 0.5, 1
  ), nrow = 3, byrow = TRUE, dimnames = list(ids, ids))

  result <- prioritize_pld(variants, R, c(0.1, 0.01, 0.001))
  selected_counts <- vapply(c(0.1, 0.01, 0.001), function(threshold) {
    sum(result$selected[result$r2_threshold == threshold])
  }, integer(1))

  expect_equal(selected_counts, c(2L, 1L, 1L))
  expect_equal(result$excluded_by_variant_id[
    result$r2_threshold == 0.1 & result$rank == 3
  ], ids[[1]])
  expect_false(any(c("PIP", "credible_set_id") %in% names(result)))
})

test_that("p-LD rejects unavailable or misordered LD instead of treating it as zero", {
  variants <- data.frame(
    CHROM = c("1", "1"),
    POS = c(100, 101),
    REF = c("A", "C"),
    ALT = c("G", "T"),
    LP = c(8, 7)
  )
  ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT,
               sep = ":")
  R <- diag(2)
  dimnames(R) <- list(ids, ids)
  R[1, 2] <- R[2, 1] <- NA_real_
  expect_error(prioritize_pld(variants, R, 0.1), "unavailable LD")

  R[1, 2] <- R[2, 1] <- 0.1
  dimnames(R) <- list(rev(ids), rev(ids))
  expect_error(prioritize_pld(variants, R, 0.1), "variant order")
})

test_that("LD validation checks symmetry, diagonal, range, and PSD when requested", {
  ids <- c("v1", "v2")
  R <- matrix(c(1, 0.2, 0.2, 1), 2, dimnames = list(ids, ids))
  expect_true(isTRUE(validate_ld_matrix(R, ids, check_psd = TRUE)))

  asymmetric <- R
  asymmetric[1, 2] <- 0.3
  expect_error(validate_ld_matrix(asymmetric), "symmetric")

  bad_diagonal <- R
  bad_diagonal[1, 1] <- 0.9
  expect_error(validate_ld_matrix(bad_diagonal), "diagonal")

  non_psd <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1),
                    nrow = 3, byrow = TRUE,
                    dimnames = list(paste0("v", 1:3), paste0("v", 1:3)))
  expect_error(validate_ld_matrix(non_psd, check_psd = TRUE), "positive semidefinite")
})
