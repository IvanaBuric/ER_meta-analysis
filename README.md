# Emotion Regulation Meta-Analysis

This repository contains the R scripts, cleaned dataset, and RStudio project file used to reproduce the analyses for this meta-analysis of the effects of mindfulness-based interventions on emotion regulation.

The repository is organised so that a user can download it, open it in RStudio, and run the scripts in order to reproduce the results.

## What is included in this repository

- 9 R scripts (`00`–`08`) for data preparation, analysis, and reporting
- the cleaned dataset used for the analyses (`Merged ER META Data Extraction_CLEAN.xlsx`)
- `study_labels.csv`, used by script `08` to label studies in the forest plot
- an `.Rproj` file for opening the project in RStudio

The repository does **not** include derived output files. These are created automatically when the scripts are run.

## Primary (59-study) and sensitivity (61-study) analyses

The dataset contains 61 studies. Two of these (Atta et al., 2024 and Hosseinian et al., 2016) were excluded from the primary analysis for documented data-integrity reasons, leaving 59 studies.

This is controlled by a single switch at the top of `00_setup.R`:

```r
EXCLUDE_INTEGRITY <- TRUE    # TRUE  = PRIMARY analysis (59 studies)
                             # FALSE = SENSITIVITY analysis (all 61 studies)
studies_excluded  <- c(6, 30)
```

- `TRUE` (default) reproduces the **primary** analysis (59 studies). Outputs are written to folders suffixed `_primary59`.
- `FALSE` reproduces the **sensitivity** analysis (all 61 studies). Outputs are written to folders suffixed `_all61`.

Because the two runs write to different folders, they do not overwrite each other. To produce both sets of results, run the pipeline once with `EXCLUDE_INTEGRITY <- TRUE` and once with `FALSE`.

## Before you start

To reproduce the analyses you will need:

- R
- RStudio Desktop

You do not need extensive R experience. The main steps are simply to open the project and run the scripts in order.

## Important note about the folder structure

You do not need to create any folders manually. After downloading this repository, keep all files exactly as they are. Do not rename files, move scripts into different folders, or change the location of the dataset, as this may prevent the scripts from running correctly. In particular, the dataset must remain at `Merged ER META Data Extraction_CLEAN.xlsx`, and `study_labels.csv` must remain in the project root.

## How to reproduce the analyses

**Step 1. Download the repository** from GitHub to your computer and unzip it if needed.

**Step 2. Open the project in RStudio** by double-clicking the `.Rproj` file. This opens the project with the correct working directory.

**Step 3. Run the scripts in order:**

```
00_setup.R
01_import_and_check_data.R
02_compute_effect_sizes.R
03_meta_analysis_models.R
04_sensitivity_analysis_ER.R
05_moderator_analysis_ER.R
06_er_md_association.R
07_risk_of_bias_and_study_quality_ER.R
08_final_reporting_tables_and_figures_ER.R
```

It is important to run them in this order because later scripts depend on files created by earlier scripts.

**Step 4. Wait for package installation if needed.** The setup script will install required R packages if they are not already available.

**Step 5. Check the generated outputs.** Output files are saved automatically in the appropriate project folders (suffixed `_primary59` or `_all61`, depending on the switch in `00_setup.R`).

## Recommended way to run each script

Some scripts may prompt you to confirm actions (e.g., by typing `y` or `n` in the console). For this reason it is recommended to:

1. Open the first script
2. Click **Source** to run it
3. Monitor the Console for any prompts
4. If prompted, type the requested input and press Enter
5. Wait until the script finishes, then open the next script and repeat

Do not leave the scripts running unattended, as some steps may require user input.

## Overview of the scripts

- `00_setup.R` — Sets up the R environment, installs/loads required packages, and defines the primary/sensitivity analysis switch.
- `01_import_and_check_data.R` — Imports the cleaned dataset and performs initial data checks.
- `02_compute_effect_sizes.R` — Computes the effect sizes used in the meta-analysis.
- `03_meta_analysis_models.R` — Runs the main meta-analytic models.
- `04_sensitivity_analysis_ER.R` — Runs sensitivity analyses for emotion regulation outcomes.
- `05_moderator_analysis_ER.R` — Runs moderator analyses.
- `06_er_md_association.R` — Examines the association between emotion regulation and mental distress effects.
- `07_risk_of_bias_and_study_quality_ER.R` — Runs analyses related to risk of bias and study quality.
- `08_final_reporting_tables_and_figures_ER.R` — Creates the final reporting tables and figures (uses `study_labels.csv`).

## Troubleshooting

- **The scripts do not run.** Check that you opened the project by clicking the `.Rproj` file, rather than opening individual scripts directly.
- **A package is missing.** Run `00_setup.R` again and allow time for package installation to finish.
- **A file cannot be found.** Make sure you have kept the repository structure unchanged, including `Merged ER META Data Extraction_CLEAN.xlsx` and `study_labels.csv`.
- **An error appears in a later script.** Make sure all earlier scripts were run successfully and in the correct order.

## Reproducibility note

This repository is intended to allow full reproduction of the analyses from the cleaned dataset onward. All file paths in the scripts are relative, so the project should run on a new computer as long as the repository structure is kept unchanged.

## Citation

If you use materials from this repository, please cite the associated paper/project.

## Contact

For questions about the repository, please contact the repository owner.
