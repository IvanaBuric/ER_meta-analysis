# =========================================================
# 00_setup.R
# Purpose: Reproduce the project's R package environment.
# Run this ONCE before sourcing scripts 01–08.
# =========================================================

# ---- 1. Use a fixed CRAN mirror ----
# Prevents the interactive "choose a CRAN mirror" menu on a
# fresh machine.
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- 2. Prefer pre-compiled binaries ----
# Avoids building packages from source (and the compiler/
# Xcode errors that can cause) on macOS and Windows.
if (.Platform$OS.type %in% c("unix", "windows")) {
  options(pkgType = "binary")
}

# ---- 3. Bootstrap renv ----
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# ---- 4. Restore the recorded package environment ----
# Installs the exact package versions stored in renv.lock.
# This is the command others should run to reproduce the
# analysis. It does NOT modify the lockfile.
renv::restore(prompt = FALSE)

# ---- 5. Load core packages ----
library(metafor)
library(clubSandwich)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(here)
library(janitor)

# ---- Analysis mode: data-integrity exclusions ----
# TRUE  = PRIMARY analysis: excludes Atta (Study_ID 6) and Hosseinian (Study_ID 30)
#         for documented data-integrity reasons -> 59 studies.
# FALSE = SENSITIVITY analysis: includes all 61 studies.
# Outputs go to folders suffixed "_primary59" or "_all61", so the two runs
# do not overwrite each other. You can also flip this in the console before
# sourcing scripts 03-08 (e.g., EXCLUDE_INTEGRITY <- FALSE).
EXCLUDE_INTEGRITY <- TRUE
studies_excluded  <- c(6, 30)

# ---- Note for maintainers only ----
# To RECORD new/updated package versions into the lockfile
# (e.g., after adding a package), run renv::snapshot() manually.
# Do NOT call snapshot() automatically here, so that running
# this script never overwrites the lockfile.
