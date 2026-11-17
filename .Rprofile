if (!nzchar(Sys.getenv("GWASFM_RPROFILE_LOADED", unset = ""))) {
  Sys.setenv(GWASFM_RPROFILE_LOADED = "1")
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
