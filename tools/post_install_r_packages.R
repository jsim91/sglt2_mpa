args <- commandArgs(trailingOnly = TRUE)
env_name <- if (length(args) >= 1) args[[1]] else Sys.getenv("CONDA_DEFAULT_ENV", unset = "")

message("Installing supplemental R packages")
if (nzchar(env_name)) {
  message("Conda environment: ", env_name)
}

cran_repo <- Sys.getenv("R_CRAN_MIRROR", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo))

ensure_cran_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

ensure_github_package <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(repo, upgrade = "never", dependencies = TRUE)
  }
}

ensure_cran_package("remotes")
ensure_cran_package("trelliscopejs")
ensure_github_package("sccomp", "MangiolaLaboratory/sccomp")

message("Supplemental R package install complete")