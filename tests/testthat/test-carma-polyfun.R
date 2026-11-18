test_that("read_ld_matrix in-place XtX scaling is bit-for-bit identical to unscaled SuSiE", {
  set.seed(101)
  n <- 10000
  p <- 60
  X <- scale(matrix(rnorm(200 * p), 200, p))
  R <- stats::cor(X)
  ids <- paste0("1:", seq_len(p) * 1000L, ":A:G")
  dimnames(R) <- list(ids, ids)

  z <- as.numeric(R[, 12] * 7.5 + R[, 40] * 5.5 + rnorm(p, sd = 0.3))
  variants <- data.frame(
    CHROM = "1",
    POS = seq_len(p) * 1000L,
    ID = paste0("rs", seq_len(p)),
    REF = "A",
    ALT = "G",
    GWAS_VARIANT_ID = paste0("1:", seq_len(p) * 1000L + 500L, ":A:G"),
    ES = z * 0.01,
    SE = 0.01,
    LP = -stats::pchisq(z^2, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10),
    SS = n,
    EZ = z,
    stringsAsFactors = FALSE
  )

  tmp_tsv <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_tsv, force = TRUE), add = TRUE)
  out_df <- cbind(data.frame(variant_id = ids, stringsAsFactors = FALSE), as.data.frame(R))
  utils::write.table(out_df, tmp_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  R_unscaled <- read_ld_matrix(tmp_tsv, scale = 1)
  R_scaled <- read_ld_matrix(tmp_tsv, scale = n - 1)

  res_u <- run_susie_rss(variants, R_unscaled, L = 5L)
  res_s <- run_susie_rss(variants, R_scaled, L = 5L)
  expect_equal(res_s$variants$PIP, res_u$variants$PIP, tolerance = 1e-14)
  expect_equal(res_s$credible_sets, res_u$credible_sets, tolerance = 1e-14)

  rr_u <- run_susie_rrss(variants, R_unscaled, reference_sample_size = 50000, L = 5L)
  rr_s <- run_susie_rrss(variants, R_scaled, reference_sample_size = 50000, L = 5L)
  expect_equal(rr_s$variants$PIP, rr_u$variants$PIP, tolerance = 1e-13)
})

test_that("run_carma recovers causal variants and flags discrepancy outliers", {
  skip_if_not_installed("CARMA")
  set.seed(202)
  n <- 1500
  p <- 80
  X <- scale(matrix(rnorm(n * p), n, p))
  X[, 16] <- scale(0.96 * X[, 15] + sqrt(1 - 0.96^2) * rnorm(n))
  R <- stats::cor(X)
  beta <- numeric(p)
  beta[15] <- 0.28
  beta[55] <- 0.22
  y <- as.numeric(X %*% beta + rnorm(n))
  z <- as.numeric(crossprod(X, y) / sqrt(n * sum(y^2) / n))
  z[16] <- -7.2

  variants <- data.frame(
    CHROM = "2",
    POS = seq_len(p) * 2000L,
    ID = paste0("rs_c_", seq_len(p)),
    REF = "C",
    ALT = "T",
    GWAS_VARIANT_ID = paste0("2:", seq_len(p) * 2000L, ":C:T"),
    ES = z * 0.02,
    SE = 0.02,
    LP = -stats::pchisq(z^2, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10),
    SS = n,
    EZ = z,
    stringsAsFactors = FALSE
  )
  ids <- paste0("2:", variants$POS, ":C:T")
  dimnames(R) <- list(ids, ids)

  res <- run_carma(variants, R, coverage = 0.95, outlier_switch = TRUE, seed = 1L)
  expect_identical(res$provenance$method_id, "carma")
  expect_true(res$variants$CARMA_OUTLIER[[16L]])
  expect_gte(res$diagnostics$outlier_count, 1L)
  top2 <- head(order(-res$variants$PIP), 2L)
  expect_true(all(c(15L, 55L) %in% top2))
  expect_gte(nrow(res$credible_sets), 2L)
})

test_that("run_polyfun_susie prioritizes functionally enriched SNPs and skips chrX cleanly", {
  set.seed(303)
  n <- 20000
  p <- 40
  X <- scale(matrix(rnorm(400 * p), 400, p))
  # Make variants 10 and 11 nearly perfectly correlated (r > 0.995)
  X[, 11] <- scale(0.998 * X[, 10] + sqrt(1 - 0.998^2) * rnorm(400))
  R <- stats::cor(X)

  z <- as.numeric(R[, 11] * 8.5 + rnorm(p, sd = 0.1))
  z[10] <- 8.55
  z[11] <- 8.54

  variants <- data.frame(
    CHROM = "10",
    POS = seq_len(p) * 1000L,
    ID = paste0("rs_pf_", seq_len(p)),
    REF = "G",
    ALT = "A",
    GWAS_VARIANT_ID = paste0("10:", seq_len(p) * 1000L + 100L, ":G:A"),
    LD_CHROM = "10",
    LD_POS = seq_len(p) * 1000L,
    LD_REF = "G",
    LD_ALT = "A",
    ES = z * 0.01,
    SE = 0.01,
    LP = -stats::pchisq(z^2, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10),
    SS = n,
    EZ = z,
    stringsAsFactors = FALSE
  )
  ids <- paste0("10:", variants$POS, ":G:A")
  dimnames(R) <- list(ids, ids)

  # Custom PolyFun prior table with swapped allele on variant 11 (100x enriched)
  tmp_priors <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_priors, force = TRUE), add = TRUE)
  prior_tbl <- data.frame(
    CHR = rep("10", p - 2L),
    BP = variants$LD_POS[seq_len(p - 2L)],
    A1 = c(rep("G", 10L), "A", rep("G", p - 13L)),
    A2 = c(rep("A", 10L), "G", rep("A", p - 13L)),
    SNPVAR = ifelse(seq_len(p - 2L) == 11L, 4e-7, 4e-9),
    stringsAsFactors = FALSE
  )
  utils::write.table(prior_tbl, tmp_priors, sep = "\t", quote = FALSE, row.names = FALSE)

  unif_res <- run_susie_rss(variants, R, L = 5L)
  pf_res <- run_polyfun_susie(
    variants, R,
    polyfun_priors = tmp_priors,
    missing_prior_policy = "baseline_floor",
    L = 5L
  )

  expect_identical(pf_res$provenance$method_id, "polyfun")
  expect_equal(pf_res$diagnostics$prior_matched_variants, p - 2L)
  expect_equal(pf_res$diagnostics$prior_baseline_assigned_variants, 2L)
  expect_gt(pf_res$variants$PIP[[11L]], 0.90)
  expect_gt(pf_res$variants$PIP[[11L]], unif_res$variants$PIP[[11L]] + 0.30)
  expect_identical(colnames(pf_res$fit$lbf_variable), variants$GWAS_VARIANT_ID)

  # Chromosome X is unsupported by precomputed autosomal PolyFun priors
  variants_x <- variants
  variants_x$CHROM <- "X"
  variants_x$LD_CHROM <- "X"
  ids_x <- paste0("X:", variants_x$POS, ":G:A")
  R_x <- R
  dimnames(R_x) <- list(ids_x, ids_x)
  tmp_dir <- tempfile("pf_dir_")
  dir.create(file.path(tmp_dir, "by_chrom"), recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  res_x <- run_polyfun_susie(variants_x, R_x, polyfun_priors = tmp_dir)
  expect_identical(res_x$diagnostics$status, "skipped_unsupported_chromosome")
  expect_true(all(is.na(res_x$variants$PIP)))
})
