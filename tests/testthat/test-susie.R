susie_fixture <- function(z = c(5.5, 5.2, 0.3, -0.2, 0.1)) {
  variants <- data.frame(
    CHROM = rep("1", 5),
    POS = 100:104,
    GWAS_CHROM = rep("1", 5),
    GWAS_POS = 1100:1104,
    GWAS_REF = c("A", "C", "G", "T", "A"),
    GWAS_ALT = c("G", "T", "A", "C", "C"),
    REF = c("A", "C", "G", "T", "A"),
    ALT = c("G", "T", "A", "C", "C"),
    EZ = z,
    LP = c(8.0, 7.5, 0.5, 0.3, 0.2),
    SS = c(900, 850, 800, 750, 700),
    stringsAsFactors = FALSE
  )
  ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT,
               sep = ":")
  R <- matrix(0.1, 5, 5, dimnames = list(ids, ids))
  diag(R) <- 1
  list(variants = variants, R = R)
}

test_that("SuSiE-RSS records PIPs, coloc-ready named matrices, and effective N", {
  fixture <- susie_fixture()
  result <- run_susie_rss(
    fixture$variants,
    fixture$R,
    L = 3,
    max_iter = 50,
    min_abs_corr = NULL
  )

  expected_gwas_ids <- paste(
    "1", 1100:1104, c("A", "C", "G", "T", "A"), c("G", "T", "A", "C", "C"),
    sep = ":"
  )
  expect_equal(result$provenance$method_id, "susie_rss")
  expect_equal(result$provenance$parameters$effective_n, 700)
  expect_equal(result$provenance$parameters$effective_n_source,
               "minimum_variant_ss")
  expect_equal(nrow(result$variants), 5L)
  expect_equal(result$variants$GWAS_VARIANT_ID, expected_gwas_ids)
  expect_equal(colnames(result$fit$alpha), expected_gwas_ids)
  expect_equal(colnames(result$fit$lbf_variable), expected_gwas_ids)
  expect_equal(names(result$fit$pip), expected_gwas_ids)
  expect_equal(result$fit$snp, expected_gwas_ids)
  expect_true(all(is.finite(result$variants$LD_R2_TO_LEAD)))
  expect_true(all(is.finite(result$variants$EXPECTED_Z)))
  expect_true(all(result$variants$PIP >= 0 & result$variants$PIP <= 1))
  expect_true(is.logical(result$diagnostics$converged))
  expect_equal(result$diagnostics$psd_check, "passed")
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
    max_iter = 50,
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
    max_iter = 50,
    min_abs_corr = NULL,
    R_mismatch = "eb"
  )
  expect_equal(configured$provenance$parameters$R_mismatch, "eb")
})

test_that("validate_ld_matrix skips O(M^3) eigen PSD check above psd_max_variants", {
  fixture <- susie_fixture()
  chk <- validate_ld_matrix(
    fixture$R,
    check_psd = TRUE,
    psd_max_variants = 3L
  )
  expect_true(isTRUE(chk))
  expect_equal(attr(chk, "psd_status"), "skipped_large_m")
})

test_that("read_ld_matrix loads binary float64 .npy + .ids.txt accurately", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  npy_path <- file.path(tmp, "mat.npy")
  ids_path <- file.path(tmp, "mat.ids.txt")
  writeLines(c("1:100:A:G", "1:200:C:T"), ids_path)

  # Construct minimal NumPy v1.0 header for 2x2 '<f8' C-order matrix
  dict_str <- "{'descr': '<f8', 'fortran_order': False, 'shape': (2, 2), }"
  pad_len <- 64L - ((10L + nchar(dict_str) + 1L) %% 64L)
  header_txt <- paste0(dict_str, strrep(" ", pad_len), "\n")
  con <- file(npy_path, "wb")
  writeBin(as.raw(c(0x93, charToRaw("NUMPY"), 0x01, 0x00)), con)
  writeBin(as.integer(nchar(header_txt)), con, size = 2L, endian = "little")
  writeBin(charToRaw(header_txt), con)
  writeBin(c(1.0, 0.35, 0.35, 1.0), con, size = 8L, endian = "little")
  close(con)

  loaded <- read_ld_matrix(npy_path)
  expect_equal(dim(loaded), c(2L, 2L))
  expect_equal(rownames(loaded), c("1:100:A:G", "1:200:C:T"))
  expect_equal(as.numeric(loaded[1, 2]), 0.35, tolerance = 1e-12)
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
