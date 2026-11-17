test_that("GWAS-VCF fields and header metadata are preserved", {
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path))
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##reference=GRCh38",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait1",
    "1\t100\trs1\tA\tG\t.\tPASS\tSNPID=source1;HM=1\tES:SE:LP:AF:SS:EZ\t0.2:0.1:2:0.2:1000:2"
  ), path)

  result <- read_gwas_vcf(path)
  expect_equal(result$metadata$reference, "GRCh38")
  expect_equal(result$sample_id, "trait1")
  expect_equal(result$variants$SNPID, "source1")
  expect_equal(result$variants$HM, "1")
  expect_equal(result$variants$ES, 0.2)
  expect_equal(result$variants$SS, 1000)
  expect_equal(result$variants$LP, 2)
})

test_that("FORMAT fields colliding with fixed and INFO fields are preserved", {
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path))
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait1",
    "1\t100\trs_fixed\tA\tG\t.\tPASS\tSNPID=source1;HM=forward\tID:SNPID:HM:ES:SE:LP\tid_format:snp_format:hm_format:0.2:0.1:2"
  ), path)

  result <- read_gwas_vcf(path)
  expect_equal(result$variants$ID, "rs_fixed")
  expect_equal(result$variants$FORMAT_ID, "id_format")
  expect_equal(result$variants$SNPID, "source1")
  expect_equal(result$variants$FORMAT_SNPID, "snp_format")
  expect_equal(result$variants$HM, "forward")
  expect_equal(result$variants$FORMAT_HM, "hm_format")
})

test_that("GWAS-VCF reader rejects multiple sample columns and malformed values", {
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path))
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait1\ttrait2"
  ), path)
  expect_error(read_gwas_vcf(path), "Exactly one")

  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait1",
    "1\t100\trs1\tA\tG\t.\tPASS\t.\tES:SE:LP\t0.2:bad:2"
  ), path)
  expect_error(read_gwas_vcf(path), "Invalid numeric value")
})

test_that("statistics are derived only when valid and remain numerically stable", {
  variants <- data.frame(
    CHROM = c("1", "1"),
    POS = c(100, 200),
    REF = c("A", "C"),
    ALT = c("G", "T"),
    ES = c(0.2, -0.3),
    SE = c(0.1, 0.1),
    P = c(0.01, 0.001)
  )

  result <- validate_gwas_stats(variants, require_p = TRUE, require_z = TRUE)
  expect_equal(result$LP, c(2, 3))
  expect_equal(result$EZ, c(2, -3))

  extreme <- variants[1, c("CHROM", "POS", "REF", "ALT")]
  extreme$LP <- 400
  expect_equal(prioritize_p(extreme)$rank, 1L)
})

test_that("invalid and duplicate variants fail explicitly", {
  variants <- data.frame(
    CHROM = c("1", "1"),
    POS = c(100, 100),
    REF = c("A", "A"),
    ALT = c("G", "G"),
    LP = c(2, 1)
  )
  expect_error(validate_gwas_stats(variants), "Duplicate")

  variants$POS <- c(0, 1)
  expect_error(validate_gwas_stats(variants), "positive integer")

  variants$POS <- c(100, 101)
  variants$ALT <- c("G,T", "C")
  expect_error(validate_gwas_stats(variants), "Multiallelic")
})

test_that("empty regions return empty prioritization results", {
  empty <- data.frame(
    CHROM = character(),
    POS = integer(),
    REF = character(),
    ALT = character(),
    LP = numeric()
  )

  expect_equal(nrow(prioritize_p(empty)), 0L)
  expect_equal(nrow(prioritize_pld(empty, matrix(numeric(), 0, 0),
                                   r2_thresholds = 0.1)), 0L)
})
