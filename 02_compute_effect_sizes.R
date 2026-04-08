############################################################
# Script: 02_compute_effect_sizes.R
# Project: ER meta-analysis_final
# Purpose:
#   Compute effect sizes for a meta-analysis of mindfulness-
#   based interventions (MBIs) on emotion regulation.
#
# Effect size approach:
#   - Standardised mean difference based on pre-post change
#     between intervention and control groups
#   - Hedges' g small-sample correction applied
#   - Positive values always indicate improvement in favour
#     of the MBI group
#
# Main assumptions:
#   - Post-intervention sample sizes are used when attrition
#     occurs
#   - Assumed pre-post correlation r = .50 for primary
#     analyses
#   - Sensitivity analyses also computed for r = .30 and .70
#
# Input:
#   data/derived/er_meta_data_checked.csv
#
# Output:
#   data/derived/effect_sizes/
#
# Author: Ivana Buric
############################################################


############################
# 1. Load required packages
############################

required_packages <- c(
  "readr",
  "dplyr",
  "tibble",
  "stringr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(readr)
library(dplyr)
library(tibble)
library(stringr)


############################
# 2. Define file locations
############################

input_file <- "data/derived/er_meta_data_checked.csv"
output_folder <- "data/derived/effect_sizes"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


############################
# 3. Read checked dataset
############################

dat <- read_csv(input_file, show_col_types = FALSE)

cat("\n----------------------------------------\n")
cat("Effect size script started\n")
cat("----------------------------------------\n")
cat("Rows read:", nrow(dat), "\n")
cat("Columns read:", ncol(dat), "\n")


#########################################
# 4. Check required columns are available
#########################################

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

if (length(missing_required_vars) > 0) {
  stop(
    "The following required columns are missing:\n",
    paste(missing_required_vars, collapse = ", ")
  )
}

cat("All required columns are present.\n")


##################################################
# 5. Additional checks before computing effect sizes
##################################################

# Check uniqueness of Effect_ID
if (anyDuplicated(dat$Effect_ID) > 0) {
  stop("Effect_ID is not unique. Please fix duplicate Effect_ID values before continuing.")
}

# Check allowed coding of Higher_Score_More_Problems
invalid_direction <- dat %>%
  filter(!Higher_Score_More_Problems %in% c(0, 1))

if (nrow(invalid_direction) > 0) {
  stop("Higher_Score_More_Problems contains values other than 0 and 1.")
}

cat("Effect_ID is unique.\n")
cat("Higher_Score_More_Problems coding is valid.\n")


###########################################################
# 6. Define function to compute change-score Hedges' g
###########################################################

# Notes on the computation:
# - The effect size compares pre-post change in the MBI group
#   against pre-post change in the control group.
# - Positive values are defined as improvement in favour of MBI.
# - The standardiser is the pooled pretest SD across groups.
# - The sampling variance uses the assumed pre-post correlation.
#
# This is a transparent and commonly used approach for
# pretest-posttest-control group designs in meta-analysis.
# Sensitivity analyses are included for alternative values of r.

compute_change_smd <- function(data, r_value) {
  
  data %>%
    mutate(
      # Use post-intervention sample sizes for effect-size variance
      n_e = Exp_N_Post,
      n_c = Ctrl_N_Post,
      
      # Direction multiplier:
      # If higher score = more problems, then improvement is a decrease.
      # If higher score = fewer problems / more adaptive ER, then
      # improvement is an increase.
      direction_multiplier = if_else(
        Higher_Score_More_Problems == 1, -1, 1
      ),
      
      # Direction-corrected change scores
      change_e = direction_multiplier * (Exp_Post_Mean - Exp_Pre_Mean),
      change_c = direction_multiplier * (Ctrl_Post_Mean - Ctrl_Pre_Mean),
      
      # Pooled pretest SD
      sd_pooled_pre = sqrt(
        (
          ((Exp_N_Pre - 1) * Exp_Pre_SD^2) +
            ((Ctrl_N_Pre - 1) * Ctrl_Pre_SD^2)
        ) /
          (Exp_N_Pre + Ctrl_N_Pre - 2)
      ),
      
      # Uncorrected standardised mean difference in change
      d_unbiased = (change_e - change_c) / sd_pooled_pre,
      
      # Degrees of freedom based on post-intervention Ns
      df_post = n_e + n_c - 2,
      
      # Small-sample correction
      J = 1 - (3 / (4 * df_post - 1)),
      
      # Hedges' g
      yi = J * d_unbiased,
      
      # Approximate sampling variance for standardised change difference
      vi = (2 * (1 - r_value) * ((n_e + n_c) / (n_e * n_c))) + 
        ((yi^2) / (2 * df_post)),
      
      sei = sqrt(vi),
      ci_lb = yi - 1.96 * sei,
      ci_ub = yi + 1.96 * sei,
      assumed_r = r_value
    )
}


##############################################
# 7. Compute primary and sensitivity datasets
##############################################

dat_es_r50 <- compute_change_smd(dat, r_value = 0.50)
dat_es_r30 <- compute_change_smd(dat, r_value = 0.30)
dat_es_r70 <- compute_change_smd(dat, r_value = 0.70)

cat("Effect sizes computed for r = .50, .30, and .70.\n")


#############################################################
# 8. Basic checks of computed effect sizes and variances
#############################################################

es_check_summary <- tibble(
  dataset = c("r_50", "r_30", "r_70"),
  n_rows = c(nrow(dat_es_r50), nrow(dat_es_r30), nrow(dat_es_r70)),
  n_missing_yi = c(sum(is.na(dat_es_r50$yi)),
                   sum(is.na(dat_es_r30$yi)),
                   sum(is.na(dat_es_r70$yi))),
  n_missing_vi = c(sum(is.na(dat_es_r50$vi)),
                   sum(is.na(dat_es_r30$vi)),
                   sum(is.na(dat_es_r70$vi))),
  n_nonpositive_vi = c(sum(dat_es_r50$vi <= 0, na.rm = TRUE),
                       sum(dat_es_r30$vi <= 0, na.rm = TRUE),
                       sum(dat_es_r70$vi <= 0, na.rm = TRUE))
)

cat("\n--- EFFECT SIZE CHECK SUMMARY ---\n")
print(es_check_summary)

# Stop if the primary dataset has missing or invalid variances
if (sum(is.na(dat_es_r50$yi)) > 0) {
  stop("Primary effect size dataset contains missing yi values.")
}

if (sum(is.na(dat_es_r50$vi)) > 0) {
  stop("Primary effect size dataset contains missing vi values.")
}

if (sum(dat_es_r50$vi <= 0, na.rm = TRUE) > 0) {
  stop("Primary effect size dataset contains non-positive variances.")
}


###########################################################
# 9. Create useful analysis variables for later scripts
###########################################################

dat_es_r50 <- dat_es_r50 %>%
  mutate(
    effect_size_type = "Change-score Hedges g",
    effect_direction = "Positive values indicate improvement in favour of MBI",
    analysis_set = case_when(
      Outcome_Domain %in% c("ER_Total", "ER_Adaptive", "ER_Maladaptive") ~ "ER",
      Outcome_Domain == "MD" ~ "MD",
      TRUE ~ "Other"
    )
  )

dat_es_r30 <- dat_es_r30 %>%
  mutate(
    effect_size_type = "Change-score Hedges g",
    effect_direction = "Positive values indicate improvement in favour of MBI",
    analysis_set = case_when(
      Outcome_Domain %in% c("ER_Total", "ER_Adaptive", "ER_Maladaptive") ~ "ER",
      Outcome_Domain == "MD" ~ "MD",
      TRUE ~ "Other"
    )
  )

dat_es_r70 <- dat_es_r70 %>%
  mutate(
    effect_size_type = "Change-score Hedges g",
    effect_direction = "Positive values indicate improvement in favour of MBI",
    analysis_set = case_when(
      Outcome_Domain %in% c("ER_Total", "ER_Adaptive", "ER_Maladaptive") ~ "ER",
      Outcome_Domain == "MD" ~ "MD",
      TRUE ~ "Other"
    )
  )


###########################################################
# 10. Inspect the largest and smallest effect sizes
###########################################################

extreme_effects <- dat_es_r50 %>%
  arrange(desc(abs(yi))) %>%
  select(
    Effect_ID,
    Study_ID,
    Comparison_ID,
    Outcome_Domain,
    Outcome_Measure,
    Group,
    yi,
    vi,
    ci_lb,
    ci_ub
  ) %>%
  slice(1:20)

cat("\n--- 20 MOST EXTREME EFFECT SIZES (r = .50) ---\n")
print(extreme_effects, n = 20)


###########################################################
# 11. Split into ER and MD datasets for later analyses
###########################################################

dat_es_r50_er <- dat_es_r50 %>%
  filter(analysis_set == "ER")

dat_es_r50_md <- dat_es_r50 %>%
  filter(analysis_set == "MD")

cat("\nPrimary dataset counts:\n")
cat("ER rows:", nrow(dat_es_r50_er), "\n")
cat("MD rows:", nrow(dat_es_r50_md), "\n")
cat("ER studies:", n_distinct(dat_es_r50_er$Study_ID), "\n")
cat("MD studies:", n_distinct(dat_es_r50_md$Study_ID), "\n")


###########################################################
# 12. Save outputs
###########################################################
# ER subset:
# Used for the primary multilevel meta-analysis of mindfulness
# interventions on emotion regulation outcomes.

# MD subset:
# Not part of the primary ER meta-analysis.
# Retained for later exploratory analyses examining the
# association between ER and mental distress effects.

# Primary analysis dataset
write_csv(
  dat_es_r50,
  file.path(output_folder, "er_meta_effect_sizes_r50.csv")
)

# Sensitivity datasets
write_csv(
  dat_es_r30,
  file.path(output_folder, "er_meta_effect_sizes_r30.csv")
)

write_csv(
  dat_es_r70,
  file.path(output_folder, "er_meta_effect_sizes_r70.csv")
)

# ER-only and MD-only files from primary dataset
write_csv(
  dat_es_r50_er,
  file.path(output_folder, "er_meta_effect_sizes_r50_ER_only.csv")
)

write_csv(
  dat_es_r50_md,
  file.path(output_folder, "er_meta_effect_sizes_r50_MD_only.csv")
)

# QC outputs
write_csv(
  es_check_summary,
  file.path(output_folder, "effect_size_check_summary.csv")
)

write_csv(
  extreme_effects,
  file.path(output_folder, "extreme_effect_sizes_r50.csv")
)


###########################################################
# 13. Final console output
###########################################################

cat("\n====================================================\n")
cat("Effect size computation finished.\n")
cat("Files saved in: data/derived/effect_sizes/\n")
cat("- er_meta_effect_sizes_r50.csv\n")
cat("- er_meta_effect_sizes_r30.csv\n")
cat("- er_meta_effect_sizes_r70.csv\n")
cat("- er_meta_effect_sizes_r50_ER_only.csv\n")
cat("- er_meta_effect_sizes_r50_MD_only.csv\n")
cat("- effect_size_check_summary.csv\n")
cat("- extreme_effect_sizes_r50.csv\n")
cat("====================================================\n")

############################################################
#14. Check extreme effecs
# Purpose:
#   Inspect the most extreme effect sizes and compare them
#   with the original extracted data.
############################################################

library(readr)
library(dplyr)

# File paths
extreme_file <- "data/derived/effect_sizes/extreme_effect_sizes_r50.csv"
full_file <- "data/derived/effect_sizes/er_meta_effect_sizes_r50.csv"

# Load datasets
extreme_es <- read_csv(extreme_file, show_col_types = FALSE)
full_es <- read_csv(full_file, show_col_types = FALSE)

cat("\n========================================\n")
cat("Extreme effect sizes inspection\n")
cat("========================================\n")

# Print all extreme rows
print(extreme_es, n = nrow(extreme_es))


############################################################
# Function to inspect a single effect size in detail
############################################################

inspect_effect <- function(effect_id) {
  
  row <- full_es %>%
    filter(Effect_ID == effect_id)
  
  cat("\n----------------------------------------\n")
  cat("Inspecting Effect_ID:", effect_id, "\n")
  cat("----------------------------------------\n")
  
  print(row, width = Inf)
  
  cat("\nKey values:\n")
  
  print(row %>%
          select(
            Study_ID,
            Comparison_ID,
            Outcome_Domain,
            Outcome_Measure,
            Exp_N_Pre,
            Exp_N_Post,
            Exp_Pre_Mean,
            Exp_Post_Mean,
            Ctrl_N_Pre,
            Ctrl_N_Post,
            Ctrl_Pre_Mean,
            Ctrl_Post_Mean,
            yi,
            vi,
            ci_lb,
            ci_ub
          ))
}


############################################################
# Loop through extreme effects automatically
############################################################

for (id in extreme_es$Effect_ID) {
  
  inspect_effect(id)
  
  readline(prompt = "Press [enter] to inspect next effect size...")
}

###########################################################
# 15. Formal outlier and influence screening
###########################################################

# Rationale:
# This block does not delete any studies automatically.
# It creates transparent flags for unusually large effects,
# unusually small standardisers, and potentially influential
# rows that should be checked before fitting the final models.

outlier_influence_screen <- dat_es_r50 %>%
  mutate(
    abs_yi = abs(yi),
    
    # Flag unusually large standardised effects
    flag_abs_yi_gt_2 = abs_yi > 2,
    flag_abs_yi_gt_3 = abs_yi > 3,
    
    # Flag very small pooled baseline SDs
    flag_sd_pooled_pre_lt_0_10 = !is.na(sd_pooled_pre) & sd_pooled_pre < 0.10,
    flag_sd_pooled_pre_lt_0_20 = !is.na(sd_pooled_pre) & sd_pooled_pre < 0.20,
    flag_sd_pooled_pre_lt_0_50 = !is.na(sd_pooled_pre) & sd_pooled_pre < 0.50,
    
    # Flag very small post-intervention sample sizes
    flag_small_post_n = !is.na(n_e) & !is.na(n_c) & (n_e < 20 | n_c < 20),
    
    # Flag strong attrition in either arm
    flag_attrition_exp_gt_20pct =
      !is.na(Exp_N_Pre) & !is.na(Exp_N_Post) &
      Exp_N_Pre > 0 &
      ((Exp_N_Pre - Exp_N_Post) / Exp_N_Pre) > 0.20,
    
    flag_attrition_ctrl_gt_20pct =
      !is.na(Ctrl_N_Pre) & !is.na(Ctrl_N_Post) &
      Ctrl_N_Pre > 0 &
      ((Ctrl_N_Pre - Ctrl_N_Post) / Ctrl_N_Pre) > 0.20,
    
    # Flag rows where the approximate sampling variance is very small
    flag_vi_below_10th_percentile =
      vi < quantile(vi, probs = 0.10, na.rm = TRUE),
    
    flag_vi_below_5th_percentile =
      vi < quantile(vi, probs = 0.05, na.rm = TRUE),
    
    # Overall flag for manual review
    any_outlier_influence_flag =
      flag_abs_yi_gt_2 |
      flag_abs_yi_gt_3 |
      flag_sd_pooled_pre_lt_0_10 |
      flag_sd_pooled_pre_lt_0_20 |
      flag_sd_pooled_pre_lt_0_50 |
      flag_small_post_n |
      flag_attrition_exp_gt_20pct |
      flag_attrition_ctrl_gt_20pct |
      flag_vi_below_10th_percentile |
      flag_vi_below_5th_percentile
  )

cat("\n--- OUTLIER / INFLUENCE SCREENING SUMMARY ---\n")

outlier_flag_summary <- tibble(
  flag = c(
    "abs_yi_gt_2",
    "abs_yi_gt_3",
    "sd_pooled_pre_lt_0_10",
    "sd_pooled_pre_lt_0_20",
    "sd_pooled_pre_lt_0_50",
    "small_post_n",
    "attrition_exp_gt_20pct",
    "attrition_ctrl_gt_20pct",
    "vi_below_10th_percentile",
    "vi_below_5th_percentile",
    "any_outlier_influence_flag"
  ),
  n_flagged = c(
    sum(outlier_influence_screen$flag_abs_yi_gt_2, na.rm = TRUE),
    sum(outlier_influence_screen$flag_abs_yi_gt_3, na.rm = TRUE),
    sum(outlier_influence_screen$flag_sd_pooled_pre_lt_0_10, na.rm = TRUE),
    sum(outlier_influence_screen$flag_sd_pooled_pre_lt_0_20, na.rm = TRUE),
    sum(outlier_influence_screen$flag_sd_pooled_pre_lt_0_50, na.rm = TRUE),
    sum(outlier_influence_screen$flag_small_post_n, na.rm = TRUE),
    sum(outlier_influence_screen$flag_attrition_exp_gt_20pct, na.rm = TRUE),
    sum(outlier_influence_screen$flag_attrition_ctrl_gt_20pct, na.rm = TRUE),
    sum(outlier_influence_screen$flag_vi_below_10th_percentile, na.rm = TRUE),
    sum(outlier_influence_screen$flag_vi_below_5th_percentile, na.rm = TRUE),
    sum(outlier_influence_screen$any_outlier_influence_flag, na.rm = TRUE)
  )
)

print(outlier_flag_summary, n = nrow(outlier_flag_summary))

flagged_outlier_rows <- outlier_influence_screen %>%
  filter(any_outlier_influence_flag) %>%
  select(
    Effect_ID,
    Study_ID,
    Comparison_ID,
    Outcome_Domain,
    Outcome_Measure,
    Group,
    yi,
    vi,
    ci_lb,
    ci_ub,
    sd_pooled_pre,
    n_e,
    n_c,
    Exp_N_Pre,
    Exp_N_Post,
    Ctrl_N_Pre,
    Ctrl_N_Post,
    starts_with("flag_")
  ) %>%
  arrange(desc(abs(yi)))

cat("\n--- FLAGGED OUTLIER / INFLUENCE ROWS ---\n")
print(flagged_outlier_rows, n = nrow(flagged_outlier_rows), width = Inf)

effects_abs_gt_3 <- outlier_influence_screen %>%
  filter(flag_abs_yi_gt_3) %>%
  select(
    Effect_ID,
    Study_ID,
    Comparison_ID,
    Outcome_Domain,
    Outcome_Measure,
    Group,
    yi,
    vi,
    ci_lb,
    ci_ub,
    sd_pooled_pre,
    Exp_N_Pre,
    Exp_N_Post,
    Ctrl_N_Pre,
    Ctrl_N_Post
  ) %>%
  arrange(desc(abs(yi)))

cat("\n--- EFFECTS WITH |g| > 3 ---\n")
print(effects_abs_gt_3, n = nrow(effects_abs_gt_3), width = Inf)

write_csv(
  outlier_flag_summary,
  file.path(output_folder, "outlier_influence_flag_summary.csv")
)

write_csv(
  flagged_outlier_rows,
  file.path(output_folder, "flagged_outlier_influence_rows.csv")
)

write_csv(
  effects_abs_gt_3,
  file.path(output_folder, "effects_abs_gt_3.csv")
)

cat("- outlier_influence_flag_summary.csv\n")
cat("- flagged_outlier_influence_rows.csv\n")
cat("- effects_abs_gt_3.csv\n")


