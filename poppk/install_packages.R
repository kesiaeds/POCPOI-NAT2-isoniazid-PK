# Install nlmixr2 + dependencies into a user-private library on Sherlock.
# This compilation is heavy (~30-60 min). Run once on an sdev node, not the
# login node, otherwise you'll be killed for CPU usage.
#
# Recommended sdev allocation:
#   sdev -t 03:00:00 -c 8 -m 32G
#
# Then:
#   ml R/4.3 gcc/12 cmake openblas
#   Rscript scripts/install_packages.R

LIB <- Sys.getenv("R_LIBS_USER",
                  file.path("~", "R", "x86_64-pc-linux-gnu-library", "4.3"))
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, .libPaths()))

options(
  repos      = c(CRAN = "https://cloud.r-project.org"),
  Ncpus      = max(1L, parallel::detectCores() - 1L),
  install.packages.check.source = "no"
)

cran <- c(
  # core tidyverse-ish + plotting
  "dplyr", "tidyr", "readr", "ggplot2", "scales", "stringr",
  # PK modelling stack
  "rxode2", "nlmixr2", "nlmixr2plot", "nlmixr2extra",
  "tidyvpc",
  # tables/labels
  "table1", "gtsummary"
)

cat("Installing into:", LIB, "\n")
cat("Using", getOption("Ncpus"), "cores\n\n")

for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("--- installing", p, "---\n")
    try(install.packages(p, dependencies = TRUE))
  } else {
    cat(p, "already installed\n")
  }
}

cat("\n--- sessionInfo ---\n"); print(sessionInfo())
cat("\nDONE.\n")
