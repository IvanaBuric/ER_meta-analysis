# Emotion Regulation Meta-Analysis

This repository contains the R scripts, cleaned dataset, and RStudio project file used to reproduce the analyses for this meta-analysis.

The repository is organised so that a user can download it, open it in RStudio, and run the scripts in order to reproduce the main results.

## What is included in this repository

This repository contains:

- 8 R scripts for data preparation, analysis, and reporting
- the cleaned dataset used for the analyses
- an `.Rproj` file for opening the project in RStudio

The repository does **not** include derived output files. These are created automatically when the scripts are run.

## Before you start

To reproduce the analyses, you will need:

- [R](https://cran.r-project.org/)
- [RStudio Desktop](https://posit.co/download/rstudio-desktop/)

You do not need extensive R experience. The main steps are simply to open the project and run the scripts in order.

## Important note about the folder structure

You do **not** need to create any folders manually.

After downloading this repository, keep all files exactly as they are. Do not rename files, move scripts into different folders, or change the location of the dataset, as this may prevent the scripts from running correctly.

## How to reproduce the analyses

### Step 1. Download the repository

Download this repository from GitHub to your computer and unzip it if needed.

### Step 2. Open the project in RStudio

Locate the `.Rproj` file in the downloaded folder and double-click it.

This will open the project in RStudio with the correct working directory.

### Step 3. Run the scripts in order

Run the scripts in the following order:

1. `00_setup.R`
2. `01_import_and_check_data.R`
3. `02_compute_effect_sizes.R`
4. `03_meta_analysis_models.R`
5. `04_sensitivity_analysis_ER.R`
6. `05_moderator_analysis_ER.R`
7. `06_er_md_association.R`
8. `07_risk_of_bias_and_study_quality_ER.R`
9. `08_final_reporting_tables_and_figures_ER.R`

It is important to run them in this order because later scripts depend on files created by earlier scripts.

### Step 4. Wait for package installation if needed

The setup script and other scripts will install required R packages if they are not already available on your computer.

### Step 5. Check the generated outputs

When the scripts have finished running, output files created by the analysis will be saved automatically in the appropriate project folders.

## Recommended way to run each script

Some scripts may prompt you to confirm actions (e.g., by typing `y` or `n` in the console).

For this reason, it is recommended to:

1. Open the first script
2. Click **Source** to run the script
3. Monitor the Console for any prompts
4. If prompted, type the requested input (e.g., `y` or `n`) and press Enter
5. Wait until the script finishes
6. Open the next script and repeat

Do not leave the scripts running unattended, as some steps may require user input.

## Overview of the scripts

### `00_setup.R`
Sets up the R environment and installs or loads required packages.

### `01_import_and_check_data.R`
Imports the cleaned dataset and performs initial data checks.

### `02_compute_effect_sizes.R`
Computes the effect sizes used in the meta-analysis.

### `03_meta_analysis_models.R`
Runs the main meta-analytic models.

### `04_sensitivity_analysis_ER.R`
Runs sensitivity analyses for emotion regulation outcomes.

### `05_moderator_analysis_ER.R`
Runs moderator analyses.

### `06_er_md_association.R`
Examines the association between emotion regulation and mental distress effects.

### `07_risk_of_bias_and_study_quality_ER.R`
Runs analyses related to risk of bias and study quality.

### `08_final_reporting_tables_and_figures_ER.R`
Creates the final reporting tables and figures.

## Troubleshooting

### The scripts do not run
First check that you opened the project by clicking the `.Rproj` file, rather than opening individual scripts directly.

### A package is missing
Run `00_setup.R` again and allow time for package installation to finish.

### A file cannot be found
Make sure you have kept the repository structure unchanged after downloading it.

### An error appears in a later script
Make sure all earlier scripts were run successfully and in the correct order.

## Reproducibility note

This repository is intended to allow full reproduction of the analyses from the cleaned dataset onward. All file paths in the scripts are relative, so the project should run on a new computer as long as the repository structure is kept unchanged.

## Citation

If you use materials from this repository, please cite the associated paper/project.

## Contact

For questions about the repository, please contact the repository owner.
