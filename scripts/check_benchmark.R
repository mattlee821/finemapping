#!/usr/bin/env Rscript
# Extract or verify variant-by-variant reproducible benchmark results from .rda + .json outputs
# (mirrors GWASprocessing/src/benchmark/check_install.py).

parse_options <- function(args) {
  opts <- list()
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- gsub("-", "_", sub("^--", "", args[[i]]), fixed = TRUE)
    opts[[key]] <- args[[i + 1L]]
  }
  opts
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

summarize_gwas_dir <- function(gwas_dir) {
  qc_path <- file.path(gwas_dir, "qc.json")
  rss_path <- file.path(gwas_dir, "methods", "susie-rss.rda")
  rrss_path <- file.path(gwas_dir, "methods", "susie-rrss.rda")
  pld_path <- file.path(gwas_dir, "methods", "p-ld.rda")
  carma_path <- file.path(gwas_dir, "methods", "carma.rda")
  polyfun_path <- file.path(gwas_dir, "methods", "polyfun.rda")

  for (req in c(qc_path, rss_path, rrss_path, pld_path, carma_path, polyfun_path)) {
    if (!file.exists(req)) {
      stop("Missing expected output file: ", req, call. = FALSE)
    }
  }
  if (file.exists(file.path(gwas_dir, "methods", "p-value.rda"))) {
    stop("Deprecated output file p-value.rda must not be generated.", call. = FALSE)
  }

  qc <- jsonlite::fromJSON(qc_path, simplifyVector = FALSE)
  e_rss <- new.env(parent = emptyenv()); load(rss_path, envir = e_rss)
  e_rrss <- new.env(parent = emptyenv()); load(rrss_path, envir = e_rrss)
  e_pld <- new.env(parent = emptyenv()); load(pld_path, envir = e_pld)
  e_carma <- new.env(parent = emptyenv()); load(carma_path, envir = e_carma)
  e_pf <- new.env(parent = emptyenv()); load(polyfun_path, envir = e_pf)

  s_rss <- e_rss$susie_rss
  s_rrss <- e_rrss$susie_rrss
  s_pld <- e_pld$p_ld
  s_carma <- e_carma$carma
  s_pf <- e_pf$polyfun

  regions <- list()
  for (rid in sort(names(s_rss$regions))) {
    r_rss <- s_rss$regions[[rid]]$result
    r_rrss <- s_rrss$regions[[rid]]$result
    r_carma <- s_carma$regions[[rid]]$result
    r_pf <- s_pf$regions[[rid]]$result
    df_pld <- s_pld$regions[[rid]]$result
    v_rss <- r_rss$variants
    v_rrss <- r_rrss$variants
    v_carma <- r_carma$variants
    v_pf <- r_pf$variants

    # Capture top 25 variants by SuSiE-RSS PIP (deterministic tie-break by GWAS_VARIANT_ID)
    ord <- order(-v_rss$PIP, -v_rss$LP, v_rss$GWAS_VARIANT_ID)
    keep <- utils::head(ord, 25L)
    sentinel_df <- data.frame(
      gwas_variant_id = as.character(v_rss$GWAS_VARIANT_ID[keep]),
      rsid = as.character(v_rss$ID[keep]),
      ez = round(as.numeric(v_rss$EZ[keep]), 6),
      lp = round(as.numeric(v_rss$LP[keep]), 6),
      pip_rss = round(as.numeric(v_rss$PIP[keep]), 6),
      pip_rrss = round(as.numeric(v_rrss$PIP[keep]), 6),
      pip_carma = round(as.numeric(v_carma$PIP[keep]), 6),
      pip_polyfun = round(as.numeric(v_pf$PIP[keep]), 6),
      cs_rss = ifelse(is.na(v_rss$credible_set_id[keep]), "", as.character(v_rss$credible_set_id[keep])),
      cs_rrss = ifelse(is.na(v_rrss$credible_set_id[keep]), "", as.character(v_rrss$credible_set_id[keep])),
      cs_carma = ifelse(is.na(v_carma$credible_set_id[keep]), "", as.character(v_carma$credible_set_id[keep])),
      cs_polyfun = ifelse(is.na(v_pf$credible_set_id[keep]), "", as.character(v_pf$credible_set_id[keep])),
      stringsAsFactors = FALSE
    )

    pld_sel <- stats::setNames(
      as.list(vapply(c(0.1, 0.01, 0.001), function(t) {
        as.integer(sum(df_pld$selected[df_pld$r2_threshold == t]))
      }, integer(1))),
      c("0.1", "0.01", "0.001")
    )

    regions[[rid]] <- list(
      region_type = s_rss$regions[[rid]]$region_type,
      variants_matched = nrow(v_rss),
      pld_selected = pld_sel,
      cs_count_rss = nrow(r_rss$credible_sets),
      cs_count_rrss = nrow(r_rrss$credible_sets),
      cs_count_carma = nrow(r_carma$credible_sets),
      cs_count_polyfun = nrow(r_pf$credible_sets),
      converged_rss = isTRUE(r_rss$diagnostics$converged),
      converged_rrss = isTRUE(r_rrss$diagnostics$converged),
      converged_carma = isTRUE(r_carma$diagnostics$converged),
      converged_polyfun = isTRUE(r_pf$diagnostics$converged),
      top_variants = sentinel_df
    )
  }

  list(
    study_id = qc$study_id,
    gwas_id = qc$gwas_id,
    overall_status = qc$overall_status,
    ancestry_ld_status = qc$ancestry_ld_status,
    genomewide_significant_variants = qc$genomewide_significant_variants,
    analyzed_regions_count = qc$analyzed_regions_count,
    cis_status = qc$cis_gene$status %||% "not_applicable",
    mhc_status = qc$mhc$status %||% "no_significant_variants",
    regions = regions
  )
}

compare_summaries <- function(actual, expected, pip_tol = 1e-3, z_tol = 1e-4) {
  diffs <- character()
  scalar_fields <- c(
    "study_id", "gwas_id", "overall_status", "ancestry_ld_status",
    "genomewide_significant_variants", "analyzed_regions_count",
    "cis_status", "mhc_status"
  )
  for (f in scalar_fields) {
    if (!identical(as.character(actual[[f]]), as.character(expected[[f]]))) {
      diffs <- c(diffs, sprintf("%s: actual=%s vs expected=%s", f, actual[[f]], expected[[f]]))
    }
  }
  act_regs <- names(actual$regions)
  exp_regs <- names(expected$regions)
  if (!identical(act_regs, exp_regs)) {
    diffs <- c(diffs, sprintf(
      "region IDs differ: actual=[%s] vs expected=[%s]",
      paste(act_regs, collapse = ", "), paste(exp_regs, collapse = ", ")
    ))
    return(diffs)
  }

  for (rid in exp_regs) {
    ar <- actual$regions[[rid]]
    er <- expected$regions[[rid]]
    for (rf in c("region_type", "variants_matched", "cs_count_rss", "cs_count_rrss", "cs_count_carma", "cs_count_polyfun")) {
      if (!identical(as.integer_or_char <- as.character(ar[[rf]]), as.character(er[[rf]]))) {
        diffs <- c(diffs, sprintf("%s.%s: actual=%s vs expected=%s", rid, rf, ar[[rf]], er[[rf]]))
      }
    }
    at <- ar$top_variants
    et <- jsonlite::fromJSON(jsonlite::toJSON(er$top_variants))
    if (nrow(at) != nrow(et) || !identical(as.character(at$gwas_variant_id), as.character(et$gwas_variant_id))) {
      diffs <- c(diffs, sprintf("%s.top_variants: variant ordering changed", rid))
      next
    }
    max_dz <- max(abs(as.numeric(at$ez) - as.numeric(et$ez)))
    max_dpip_rss <- max(abs(as.numeric(at$pip_rss) - as.numeric(et$pip_rss)))
    max_dpip_rrss <- max(abs(as.numeric(at$pip_rrss) - as.numeric(et$pip_rrss)))
    max_dpip_carma <- max(abs(as.numeric(at$pip_carma) - as.numeric(et$pip_carma)))
    max_dpip_pf <- max(abs(as.numeric(at$pip_polyfun) - as.numeric(et$pip_polyfun)))
    if (!is.finite(max_dz) || max_dz > z_tol) {
      diffs <- c(diffs, sprintf("%s: max |dZ| = %.6f > %.6f", rid, max_dz, z_tol))
    }
    if (!is.finite(max_dpip_rss) || max_dpip_rss > pip_tol) {
      diffs <- c(diffs, sprintf("%s: max |dPIP_rss| = %.6f > %.6f", rid, max_dpip_rss, pip_tol))
    }
    if (!is.finite(max_dpip_rrss) || max_dpip_rrss > pip_tol) {
      diffs <- c(diffs, sprintf("%s: max |dPIP_rrss| = %.6f > %.6f", rid, max_dpip_rrss, pip_tol))
    }
    if (!is.finite(max_dpip_carma) || max_dpip_carma > pip_tol) {
      diffs <- c(diffs, sprintf("%s: max |dPIP_carma| = %.6f > %.6f", rid, max_dpip_carma, pip_tol))
    }
    if (!is.finite(max_dpip_pf) || max_dpip_pf > pip_tol) {
      diffs <- c(diffs, sprintf("%s: max |dPIP_polyfun| = %.6f > %.6f", rid, max_dpip_pf, pip_tol))
    }
  }
  diffs
}

opts <- parse_options(commandArgs(trailingOnly = TRUE))
mode <- opts$mode %||% "check"
gwas_dir <- opts$gwas_dir
expected_file <- opts$expected

if (is.null(gwas_dir) || is.null(expected_file)) {
  stop("Supply --gwas-dir <dir> and --expected <json>.", call. = FALSE)
}

actual <- summarize_gwas_dir(gwas_dir)
if (mode == "write") {
  dir.create(dirname(expected_file), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(actual, expected_file, auto_unbox = TRUE, pretty = TRUE)
  total_v <- sum(vapply(actual$regions, function(r) r$variants_matched, integer(1)))
  cat(sprintf(
    "WROTE %s: %d region(s), %d matched variants\n",
    expected_file, actual$analyzed_regions_count, total_v
  ))
  quit(status = 0L)
}

if (!file.exists(expected_file)) {
  stop("Expected benchmark file does not exist: ", expected_file, call. = FALSE)
}
expected <- jsonlite::fromJSON(expected_file, simplifyVector = FALSE)
diffs <- compare_summaries(actual, expected)
total_v <- sum(vapply(actual$regions, function(r) r$variants_matched, integer(1)))
if (length(diffs) > 0L) {
  cat(sprintf(
    "FAIL %s/%s (%d region(s), %d matched variants):\n%s\n",
    actual$study_id, actual$gwas_id, actual$analyzed_regions_count, total_v,
    paste0("  - ", diffs, collapse = "\n")
  ))
  quit(status = 1L)
}
cat(sprintf(
  "PASS %s/%s: %d expected region(s), %d LD-matched variants, exact PIP/CS concordance\n",
  actual$study_id, actual$gwas_id, actual$analyzed_regions_count, total_v
))
