profile <- Sys.getenv("R_PROFILE_USER", unset = "")
if (nzchar(profile)) {
  project <- dirname(normalizePath(profile))
} else {
  project <- Sys.getenv("RENV_PROJECT", unset = getwd())
}
project <- normalizePath(project)
Sys.setenv(RENV_PROJECT = project)
source(file.path(project, "renv", "activate.R"))
