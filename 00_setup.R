# =========================================================
# 00_setup.R
# Purpose: Set up a reproducible R environment for the
# meta-analysis project and install/load core packages
# =========================================================

# ---- 1. Install renv if needed ----
# renv creates a project-specific package library so that
# your analysis can be reproduced later on another computer
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# ---- 2. Initialise renv for this project ----
# This sets up the reproducible project environment
renv::init()

# ---- 3. Install core packages needed for this project ----
# metafor = meta-analysis
# readr   = importing data
# dplyr   = data manipulation
# tidyr   = tidying data
# stringr = working with text/labels
# here    = stable file paths within the project
# janitor = cleaning variable names and quick tabulations
install.packages(c(
  "metafor",
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "here",
  "janitor"
))

# ---- 4. Load the packages ----
library(metafor)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(janitor)

# ---- 5. Save the package versions used in this project ----
# This creates/updates the lockfile for reproducibility
renv::snapshot()

