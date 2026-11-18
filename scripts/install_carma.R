#!/usr/bin/env Rscript
# Install pinned ZikunY/CARMA from source with portable gsl-config + R BLAS/LAPACK flags.
# Removes upstream's precompiled arm64 .o/.so artifacts and hardcoded /opt/homebrew paths.

carma_commit <- Sys.getenv(
  "GWASFM_CARMA_COMMIT",
  "856f0fb48cc42f1ad895600f0a4d80ad8451c2ff"
)
force_reinstall <- tolower(Sys.getenv("GWASFM_FORCE_CARMA", "false")) %in% c("1", "true", "yes")

target_lib <- commandArgs(trailingOnly = TRUE)
target_lib <- if (length(target_lib) >= 1L && nzchar(target_lib[[1L]])) {
  target_lib[[1L]]
} else {
  .libPaths()[[1L]]
}

if (!force_reinstall && requireNamespace("CARMA", lib.loc = target_lib, quietly = TRUE)) {
  cat("CARMA already installed in", target_lib, "\n")
  quit(save = "no", status = 0L)
}

cran_deps <- c("Rcpp", "RcppArmadillo", "RcppGSL", "glmnet", "dplyr", "MASS", "Matrix")
missing_deps <- cran_deps[!vapply(cran_deps, function(pkg) {
  requireNamespace(pkg, lib.loc = c(target_lib, .libPaths()), quietly = TRUE)
}, logical(1))]

if (length(missing_deps) > 0L) {
  cat("Installing CRAN dependencies for CARMA:", paste(missing_deps, collapse = ", "), "\n")
  utils::install.packages(
    missing_deps,
    lib = target_lib,
    repos = "https://cloud.r-project.org",
    quiet = TRUE
  )
}

url <- sprintf("https://github.com/ZikunY/CARMA/archive/%s.tar.gz", carma_commit)
work_dir <- tempfile(pattern = "gwasfm_carma_src_")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

tar_path <- file.path(work_dir, "carma.tar.gz")
utils::download.file(url, destfile = tar_path, mode = "wb", quiet = TRUE)
utils::untar(tar_path, exdir = work_dir)

src_dirs <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
pkg_dir <- src_dirs[grepl("^CARMA-", basename(src_dirs))][[1L]]

# Strip upstream's checked-in .o / .so binaries and .DS_Store / Rd2pdf clutter
stale_bins <- list.files(
  file.path(pkg_dir, "src"),
  pattern = "\\.(o|so|dll)$",
  full.names = TRUE
)
if (length(stale_bins) > 0L) unlink(stale_bins, force = TRUE)

gsl_cflags <- tryCatch(
  trimws(system2("gsl-config", "--cflags", stdout = TRUE)),
  error = function(e) ""
)
gsl_libs <- tryCatch(
  trimws(system2("gsl-config", "--libs", stdout = TRUE)),
  error = function(e) "-lgsl -lgslcblas -lm"
)

makevars_lines <- c(
  "CXX_STD = CXX14",
  paste("PKG_CPPFLAGS =", gsl_cflags),
  paste("PKG_LIBS =", gsl_libs, "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)")
)
writeLines(makevars_lines, con = file.path(pkg_dir, "src", "Makevars"))

env_vars <- c(
  paste0("R_LIBS=", paste(c(target_lib, .libPaths()), collapse = .Platform$path.sep)),
  "R_PROFILE_USER="
)
status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", target_lib), pkg_dir),
  env = env_vars
)
if (status != 0L) {
  stop("Failed to compile and install CARMA (", carma_commit, ").", call. = FALSE)
}
cat("Successfully installed pinned CARMA (", substr(carma_commit, 1L, 8L), ") into ", target_lib, "\n", sep = "")
