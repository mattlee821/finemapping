if (!nzchar(Sys.getenv("GWASFM_RPROFILE_LOADED", unset = ""))) {
  Sys.setenv(GWASFM_RPROFILE_LOADED = "1")
  fast_lib <- Sys.getenv("GWASFM_DIRECT_RLIB", unset = "")
  if (!nzchar(fast_lib)) {
    lib_root <- Sys.getenv(
      "RENV_PATHS_LIBRARY",
      unset = file.path(Sys.getenv("HOME"), ".cache", "gwasfinemapping", "renv-library")
    )
    cands <- Sys.glob(file.path(lib_root, "*", "*", "*"))
    cands <- cands[dir.exists(file.path(cands, "susieR"))]
    if (length(cands) > 0L) fast_lib <- cands[[1L]]
  }
  if (nzchar(fast_lib) && dir.exists(fast_lib) &&
      !identical(Sys.getenv("GWASFM_FORCE_RENV", unset = ""), "1")) {
    .libPaths(unique(c(normalizePath(fast_lib), .libPaths())))
  } else {
    profile <- Sys.getenv("R_PROFILE_USER", unset = "")
    if (nzchar(profile) && file.exists(profile)) {
      project <- dirname(normalizePath(profile))
    } else {
      project <- Sys.getenv("RENV_PROJECT", unset = getwd())
    }
    if (file.exists(file.path(project, "renv", "activate.R"))) {
      project <- normalizePath(project)
      Sys.setenv(RENV_PROJECT = project)
      source(file.path(project, "renv", "activate.R"))
    }
  }
}
