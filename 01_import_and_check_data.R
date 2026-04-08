############################################################
# PROJECT: ER meta-analysis_final
# SCRIPT: 01_import_and_check_data.R
# PURPOSE:
#   1. Import the cleaned Excel dataset
#   2. Check that the data were imported correctly
#   3. Run transparent and reproducible data checks
#   4. Save checked files for later analysis scripts
#
# AUTHOR NOTES:
#   - This script does NOT compute effect sizes yet.
#   - This script is intended to be shareable and reproducible.
#   - Run this script before any analysis or meta-analysis scripts.
############################################################


############################
# 0. Install and Load required packages
############################

# These packages are used for:
# - readxl: importing Excel files
# - dplyr: data wrangling
# - stringr: cleaning text variables
# - janitor: quick tabulations and cleaner outputs
# - readr: exporting CSV files
# - tibble: neat data display
required_packages <- c(
  "readxl",
  "dplyr",
  "stringr",
  "janitor",
  "readr",
  "tibble"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed)) {
    install.packages(pkg)
  }
}

install.packages("readxl")
############################################################
# Load packages
############################################################

library(readxl)
library(dplyr)
library(stringr)
library(janitor)
library(readr)
library(tibble)

#####################################
# 1. Define file paths for the project
#####################################

# Input file: cleaned Excel dataset
input_file <- "data/Merged ER META Data Extraction_CLEAN.xlsx"

# Output folder for files created by this script
output_folder <- "data/derived"

# Create output folder if it does not already exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


########################
# 2. Import the dataset
########################

# Read the first sheet from the Excel workbook
# If your data are located on another sheet, replace sheet = 1
# with the correct sheet name or sheet number
dat_raw <- read_excel(input_file, sheet = 1)

# Create a working copy
dat <- dat_raw


#################################################
# 3. Basic inspection of the imported dataset
#################################################

cat("\n--- DATA DIMENSIONS ---\n")
print(dim(dat))

cat("\n--- VARIABLE NAMES ---\n")
print(names(dat))

cat("\n--- DATA STRUCTURE ---\n")
str(dat)

cat("\n--- FIRST 6 ROWS ---\n")
print(head(dat))


#############################################################
# 4. Remove fully empty rows and clean text-type variables
#############################################################

# Remove rows that are completely empty
dat <- dat %>%
  filter(!if_all(everything(), ~ is.na(.)))

# Remove extra spaces from text variables
dat <- dat %>%
  mutate(across(where(is.character), ~ str_squish(.)))


##########################################################
# 5. Convert variables to the intended type where needed
##########################################################

# List of variables that should be numeric
numeric_vars <- c(
  "Effect_ID",
  "Study_ID",
  "Higher_Score_More_Problems",
  "Exp_N_Pre",
  "Exp_N_Post",
  "Exp_Pre_Mean",
  "Exp_Pre_SD",
  "Exp_Post_Mean",
  "Exp_Post_SD",
  "Ctrl_N_Pre",
  "Ctrl_N_Post",
  "Ctrl_Pre_Mean",
  "Ctrl_Pre_SD",
  "Ctrl_Post_Mean",
  "Ctrl_Post_SD",
  "Preregistered",
  "Population",
  "Recruitment_Source",
  "Age_Mean",
  "%_Female",
  "Ethnicity",
  "MBI_type",
  "MBI_Duration_Weeks",
  "MBI_Frequency_Per_Week",
  "MBI_Session_Length_Minutes",
  "MBI_Delivery_Mode",
  "Control_Type",
  "ROB_D1",
  "ROB_D2",
  "ROB_D3",
  "ROB_D4",
  "ROB_D5"
)


# First, turn empty strings and text "NA" into real missing values
dat <- dat %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  mutate(across(where(is.character), ~ na_if(., "NA")))

# Then convert selected variables to numeric
dat <- dat %>%
  mutate(across(any_of(numeric_vars), as.numeric))

library(dplyr)

dat <- dat %>%
  mutate(
    # Recode Population into broader groups
    Population_collapsed = case_when(
      Population %in% c(0, 1)        ~ "Nonclinical",
      Population %in% c(2, 3)        ~ "Clinical",
      Population == 4                ~ "AtRisk",
      TRUE                           ~  NA_character_  # treat unknown or missing as NA 
    ),
    
    # Recode Control_Type into Passive vs Active
    Control_collapsed = case_when(
      Control_Type %in% c(0, 1)      ~ "Passive",
      Control_Type %in% c(2, 3, 4)    ~ "Active",
      TRUE                           ~ NA_character_  # treat unknown or missing as NA 
    )
  )

# Make domain, population and control type categorical

dat <- dat %>%
  mutate(
    Outcome_Domain = factor(
      Outcome_Domain,
      levels = c("ER_Adaptive", "ER_Maladaptive", "ER_Total", "MD")
    ),
    Population_collapsed = factor(
      Population_collapsed,
      levels = c("Nonclinical", "Clinical", "AtRisk")
    ),
    Control_collapsed = factor(
      Control_collapsed,
      levels = c("Passive", "Active")
    )
  )
# Check the distribution of collapsed categories
cat("\nPopulation_collapsed distribution:\n")
print(dat %>% count(Population_collapsed) %>% mutate(pct = n / sum(n) * 100))

cat("\nControl_collapsed distribution:\n")
print(dat %>% count(Control_collapsed) %>% mutate(pct = n / sum(n) * 100))

# Keep grouping and text variables as character
dat <- dat %>%
  mutate(
    Comparison_ID = as.character(Comparison_ID),
    Outcome_Domain = as.character(Outcome_Domain),
    Outcome_Measure = as.character(Outcome_Measure),
    Group = as.character(Group),
  )


#############################################
# 6. Check whether required columns are present
#############################################

required_vars <- c(
  "Effect_ID",
  "Study_ID",
  "Comparison_ID",
  "Outcome_Domain",
  "Higher_Score_More_Problems",
  "Outcome_Measure",
  "Group",
  "Exp_N_Pre",
  "Exp_N_Post",
  "Exp_Pre_Mean",
  "Exp_Pre_SD",
  "Exp_Post_Mean",
  "Exp_Post_SD",
  "Ctrl_N_Pre",
  "Ctrl_N_Post",
  "Ctrl_Pre_Mean",
  "Ctrl_Pre_SD",
  "Ctrl_Post_Mean",
  "Ctrl_Post_SD"
)

missing_required_vars <- setdiff(required_vars, names(dat))

cat("\n--- REQUIRED COLUMN CHECK ---\n")
if (length(missing_required_vars) == 0) {
  cat("All required columns are present.\n")
} else {
  cat("Missing required columns:\n")
  print(missing_required_vars)
}


###########################################
# 7. Check key identifiers and row counts
###########################################

cat("\n--- IDENTIFIER CHECKS ---\n")
cat("Number of rows:", nrow(dat), "\n")
cat("Number of unique studies:", n_distinct(dat$Study_ID), "\n")
cat("Number of unique comparisons:", n_distinct(dat$Comparison_ID), "\n")
cat("Number of unique effect IDs:", n_distinct(dat$Effect_ID), "\n")

# Check whether Effect_ID is unique
effect_id_duplicates <- dat %>%
  count(Effect_ID) %>%
  filter(n > 1)

cat("\nDuplicate Effect_ID values:\n")
print(effect_id_duplicates)

# Check whether Comparison_ID is unique
comparison_id_duplicates <- dat %>%
  count(Comparison_ID) %>%
  filter(n > 1)

cat("\nDuplicate Comparison_ID values:\n")
print(comparison_id_duplicates)


##########################################
# 8. Check coding of key analysis variables
##########################################

cat("\n--- CODING CHECKS ---\n")

cat("\nOutcome_Domain frequencies:\n")
print(tabyl(dat$Outcome_Domain))

cat("\nHigher_Score_More_Problems frequencies:\n")
print(tabyl(dat$Higher_Score_More_Problems))

cat("\nGroup frequencies:\n")
print(tabyl(dat$Group))


##################################################
# 9. Summarise missing values across all variables
##################################################

missing_summary <- tibble(
  variable = names(dat),
  n_missing = sapply(dat, function(x) sum(is.na(x))),
  pct_missing = round(100 * n_missing / nrow(dat), 2)
)

cat("\n--- MISSING VALUE SUMMARY ---\n")
print(missing_summary)

# Save missing-value summary
write_csv(
  missing_summary,
  file.path(output_folder, "missing_value_summary.csv")
)


#############################################################
# 10. Flag impossible or suspicious values in key variables
#############################################################

# These checks flag rows for review.
# They do NOT automatically delete or modify data.

dat_checks <- dat %>%
  mutate(
    flag_missing_ids =
      is.na(Effect_ID) | is.na(Study_ID) | is.na(Comparison_ID),
    
    flag_invalid_direction =
      !Higher_Score_More_Problems %in% c(0, 1),
    
    flag_nonpositive_sample_size =
      (!is.na(Exp_N_Pre)  & Exp_N_Pre  <= 0) |
      (!is.na(Exp_N_Post) & Exp_N_Post <= 0) |
      (!is.na(Ctrl_N_Pre) & Ctrl_N_Pre <= 0) |
      (!is.na(Ctrl_N_Post)& Ctrl_N_Post <= 0),
    
    flag_noninteger_sample_size =
      (!is.na(Exp_N_Pre)  & Exp_N_Pre  != round(Exp_N_Pre)) |
      (!is.na(Exp_N_Post) & Exp_N_Post != round(Exp_N_Post)) |
      (!is.na(Ctrl_N_Pre) & Ctrl_N_Pre != round(Ctrl_N_Pre)) |
      (!is.na(Ctrl_N_Post)& Ctrl_N_Post!= round(Ctrl_N_Post)),
    
    flag_negative_sd =
      (!is.na(Exp_Pre_SD)   & Exp_Pre_SD   < 0) |
      (!is.na(Exp_Post_SD)  & Exp_Post_SD  < 0) |
      (!is.na(Ctrl_Pre_SD)  & Ctrl_Pre_SD  < 0) |
      (!is.na(Ctrl_Post_SD) & Ctrl_Post_SD < 0),
    
    flag_zero_sd =
      (!is.na(Exp_Pre_SD)   & Exp_Pre_SD   == 0) |
      (!is.na(Exp_Post_SD)  & Exp_Post_SD  == 0) |
      (!is.na(Ctrl_Pre_SD)  & Ctrl_Pre_SD  == 0) |
      (!is.na(Ctrl_Post_SD) & Ctrl_Post_SD == 0),
    
    flag_missing_core_stats =
      is.na(Exp_Pre_Mean)  | is.na(Exp_Post_Mean) |
      is.na(Ctrl_Pre_Mean) | is.na(Ctrl_Post_Mean) |
      is.na(Exp_Pre_SD)    | is.na(Exp_Post_SD)   |
      is.na(Ctrl_Pre_SD)   | is.na(Ctrl_Post_SD),
    
    flag_post_n_greater_than_pre =
      (!is.na(Exp_N_Pre)  & !is.na(Exp_N_Post)  & Exp_N_Post  > Exp_N_Pre) |
      (!is.na(Ctrl_N_Pre) & !is.na(Ctrl_N_Post) & Ctrl_N_Post > Ctrl_N_Pre),
    
    # Large SD relative to mean can be real, but it is worth reviewing
    flag_suspicious_sd_ratio =
      (!is.na(Exp_Pre_Mean)  & !is.na(Exp_Pre_SD)  & abs(Exp_Pre_Mean)  > 0 & Exp_Pre_SD  / abs(Exp_Pre_Mean)  > 3) |
      (!is.na(Exp_Post_Mean) & !is.na(Exp_Post_SD) & abs(Exp_Post_Mean) > 0 & Exp_Post_SD / abs(Exp_Post_Mean) > 3) |
      (!is.na(Ctrl_Pre_Mean) & !is.na(Ctrl_Pre_SD) & abs(Ctrl_Pre_Mean) > 0 & Ctrl_Pre_SD / abs(Ctrl_Pre_Mean) > 3) |
      (!is.na(Ctrl_Post_Mean)& !is.na(Ctrl_Post_SD)& abs(Ctrl_Post_Mean)> 0 & Ctrl_Post_SD/ abs(Ctrl_Post_Mean)> 3)
  ) %>%
  mutate(
    any_data_flag =
      flag_missing_ids |
      flag_invalid_direction |
      flag_nonpositive_sample_size |
      flag_noninteger_sample_size |
      flag_negative_sd |
      flag_zero_sd |
      flag_missing_core_stats |
      flag_post_n_greater_than_pre |
      flag_suspicious_sd_ratio
  )


#######################################################
# 11. Review rows that were flagged by automatic checks
#######################################################

# Count how many rows triggered any flag
cat("\n--- NUMBER OF FLAGGED ROWS ---\n")
cat(sum(dat_checks$any_data_flag), "rows flagged for review\n")

# Extract rows with any flag
flagged_rows <- dat_checks %>%
  filter(any_data_flag) %>%
  select(
    Effect_ID,
    Study_ID,
    Comparison_ID,
    Outcome_Domain,
    Outcome_Measure,
    Group,
    starts_with("flag_"),
  )

cat("\n--- FLAGGED ROWS ---\n")
print(flagged_rows, n = nrow(flagged_rows))

# Save flagged rows for manual inspection
write_csv(
  flagged_rows,
  file.path(output_folder, "flagged_rows.csv")
)

########################################################
# 12. Check for possible duplicate substantive entries
########################################################

# This checks whether the same study/comparison/outcome/group
# appears more than once
possible_duplicate_rows <- dat %>%
  count(Study_ID, Comparison_ID, Outcome_Domain, Outcome_Measure, Group) %>%
  filter(n > 1)

cat("\n--- POSSIBLE DUPLICATE COMPARISONS ---\n")
print(possible_duplicate_rows)

write_csv(
  possible_duplicate_rows,
  file.path(output_folder, "possible_duplicate_comparisons.csv")
)


############################################################
# 13. Check whether studies contain multiple comparisons
############################################################

study_comparison_summary <- dat %>%
  group_by(Study_ID) %>%
  summarise(
    n_rows = n(),
    n_comparisons = n_distinct(Comparison_ID),
    n_outcomes = n_distinct(Outcome_Measure),
    .groups = "drop"
  ) %>%
  arrange(Study_ID)

cat("\n--- STUDY-LEVEL SUMMARY ---\n")
print(study_comparison_summary, n = nrow(study_comparison_summary))

write_csv(
  study_comparison_summary,
  file.path(output_folder, "study_comparison_summary.csv")
)


############################################################
# 14. Save checked datasets for later analysis scripts
############################################################

# Save a clean checked dataset
write_csv(
  dat,
  file.path(output_folder, "er_meta_data_checked.csv")
)

# Save dataset with automatic flags included
write_csv(
  dat_checks,
  file.path(output_folder, "er_meta_data_checked_with_flags.csv")
)


#########################################
# 15. Final console message
#########################################

cat("\n====================================================\n")
cat("Data import and validation script finished.\n")
cat("Files saved in: data/derived/\n")
cat("- er_meta_data_checked.csv\n")
cat("- er_meta_data_checked_with_flags.csv\n")
cat("- missing_value_summary.csv\n")
cat("- flagged_rows.csv\n")
cat("- possible_duplicate_comparisons.csv\n")
cat("- study_comparison_summary.csv\n")
cat("====================================================\n")

