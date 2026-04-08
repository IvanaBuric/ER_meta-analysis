############################################################
# Script: 08_final_reporting_tables_and_figures_ER.R
# Project: ER meta-analysis_final
# Purpose:
#   Create final reporting tables and figures for the
#   multilevel meta-analysis of mindfulness-based
#   interventions (MBIs) on emotion regulation (ER).
#
# This script collates outputs from previous scripts:
#   - Script 03: Main ER meta-analytic models
#   - Script 04: ER sensitivity analyses
#   - Script 05: ER moderator analyses
#   - Script 06: ER-MD association analyses
#   - Script 07: Risk-of-bias analyses
#
# Main reporting outputs:
#   - Final primary-results summary table
#   - Final sensitivity summary table
#   - Final moderator summary table
#   - Final risk-of-bias summary table
#   - Final ER-MD association summary table
#   - Study-level forest plot for aggregated ER effects
#   - Combined reporting .RData object
#
# Input:
#   data/derived/meta_models/
#   data/derived/sensitivity_ER/
#   data/derived/moderators_ER/
#   data/derived/er_md_association/
#   data/derived/risk_of_bias_ER/
#
# Output:
#   data/derived/final_reporting_ER/
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
  "stringr",
  "ggplot2"
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
library(ggplot2)


############################
# 2. Define file locations
############################

input_meta_folder <- "data/derived/meta_models"
input_sensitivity_folder <- "data/derived/sensitivity_ER"
input_moderators_folder <- "data/derived/moderators_ER"
input_association_folder <- "data/derived/er_md_association"
input_rob_folder <- "data/derived/risk_of_bias_ER"

output_folder <- "data/derived/final_reporting_ER"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


############################
# 3. Helper functions
############################

save_object_summary <- function(object, path, title = NULL) {
  sink(path)
  if (!is.null(title)) {
    cat(title, "\n")
    cat(strrep("=", nchar(title)), "\n\n")
  }
  print(object)
  sink()
}

assert_required_columns <- function(data, cols, object_name) {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", object_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}


############################
# 4. Read required input files
############################

overall_er_model_summary <- read_csv(
  file.path(input_meta_folder, "overall_ER_model_summary.csv"),
  show_col_types = FALSE
)

domain_specific_model_summary <- read_csv(
  file.path(input_meta_folder, "domain_specific_model_summary.csv"),
  show_col_types = FALSE
)

study_level_er <- read_csv(
  file.path(input_meta_folder, "study_level_aggregated_ER_effects.csv"),
  show_col_types = FALSE
)

all_sensitivity_results <- read_csv(
  file.path(input_sensitivity_folder, "all_sensitivity_results_ER.csv"),
  show_col_types = FALSE
)

sensitivity_n_table <- read_csv(
  file.path(input_sensitivity_folder, "sensitivity_sample_sizes_ER.csv"),
  show_col_types = FALSE
)

moderator_overview <- read_csv(
  file.path(input_moderators_folder, "moderator_overview_ER.csv"),
  show_col_types = FALSE
)

domain_interaction_overview <- read_csv(
  file.path(input_moderators_folder, "domain_interaction_overview_ER.csv"),
  show_col_types = FALSE
)

er_md_meta_regression_summary <- read_csv(
  file.path(input_association_folder, "er_md_meta_regression_summary.csv"),
  show_col_types = FALSE
)

descriptive_correlations_er_md <- read_csv(
  file.path(input_association_folder, "descriptive_correlations_ER_MD.csv"),
  show_col_types = FALSE
)

sensitivity_er_md_results <- read_csv(
  file.path(input_association_folder, "sensitivity_ER_MD_results.csv"),
  show_col_types = FALSE
)

sensitivity_er_md_excluding_extreme <- read_csv(
  file.path(input_association_folder, "sensitivity_ER_MD_excluding_extreme_ER_studies.csv"),
  show_col_types = FALSE
)

risk_of_bias_overview <- read_csv(
  file.path(input_rob_folder, "risk_of_bias_overview_ER.csv"),
  show_col_types = FALSE
)

risk_of_bias_sensitivity <- read_csv(
  file.path(input_rob_folder, "sensitivity_excluding_high_risk_studies_ER.csv"),
  show_col_types = FALSE
)

cat("\n----------------------------------------\n")
cat("Final reporting script started\n")
cat("----------------------------------------\n")
cat("All required files were read successfully.\n")


############################
# 5. Check critical columns
############################

assert_required_columns(
  overall_er_model_summary,
  c("estimate", "se_robust", "ci_lb_robust", "ci_ub_robust",
    "p_value_robust", "k_effect_sizes", "n_studies", "n_comparisons",
    "tau2_level_3", "tau2_level_2"),
  "overall_er_model_summary"
)

assert_required_columns(
  domain_specific_model_summary,
  c("Outcome_Domain", "estimate", "se_robust", "ci_lb_robust",
    "ci_ub_robust", "p_value_robust", "k_effect_sizes", "n_studies",
    "n_comparisons", "tau2_level_3", "tau2_level_2"),
  "domain_specific_model_summary"
)

assert_required_columns(
  study_level_er,
  c("Study_ID", "yi_agg", "vi_agg", "sei_agg", "ci_lb", "ci_ub", "k_effects"),
  "study_level_er"
)

assert_required_columns(
  all_sensitivity_results,
  c("analysis", "estimate", "se_robust", "p_value_robust",
    "ci_lb_robust", "ci_ub_robust", "k_effect_sizes", "n_studies",
    "n_comparisons", "tau2_level_3", "tau2_level_2", "QE", "QE_df", "QE_p"),
  "all_sensitivity_results"
)

assert_required_columns(
  moderator_overview,
  c("moderator", "label", "type", "status", "k_effect_sizes",
    "n_studies", "n_comparisons", "omnibus_stat", "omnibus_df_num",
    "omnibus_df_den", "omnibus_p"),
  "moderator_overview"
)

assert_required_columns(
  domain_interaction_overview,
  c("moderator", "label", "type", "status", "k_effect_sizes",
    "n_studies", "n_comparisons", "interaction_stat",
    "interaction_df_num", "interaction_df_den", "interaction_p"),
  "domain_interaction_overview"
)

assert_required_columns(
  risk_of_bias_overview,
  c("moderator", "label", "status", "k_effect_sizes",
    "n_studies", "n_comparisons", "omnibus_stat",
    "omnibus_df_num", "omnibus_df_den", "omnibus_p"),
  "risk_of_bias_overview"
)

assert_required_columns(
  risk_of_bias_sensitivity,
  c("moderator", "label", "analysis", "status", "estimate",
    "se_robust", "ci_lb_robust", "ci_ub_robust", "p_value_robust",
    "k_effect_sizes", "n_studies", "n_comparisons", "n_studies_total",
    "n_studies_high", "tau2_level_3", "tau2_level_2", "QE", "QE_df", "QE_p"),
  "risk_of_bias_sensitivity"
)

assert_required_columns(
  er_md_meta_regression_summary,
  c("intercept", "slope_ER", "intercept_se", "slope_se",
    "intercept_ci_lb", "intercept_ci_ub", "slope_ci_lb", "slope_ci_ub",
    "intercept_p", "slope_p", "tau2", "QE", "QE_df", "QE_p",
    "QM", "QM_df", "QM_p", "k_studies"),
  "er_md_meta_regression_summary"
)

assert_required_columns(
  descriptive_correlations_er_md,
  c("statistic", "estimate", "p_value", "ci_lb", "ci_ub", "n_studies"),
  "descriptive_correlations_er_md"
)


###########################################################
# 6. Final primary-results summary tables
###########################################################

final_primary_summary <- overall_er_model_summary %>%
  mutate(
    analysis = "Primary overall ER multilevel model"
  ) %>%
  select(
    analysis,
    estimate,
    se_robust,
    ci_lb_robust,
    ci_ub_robust,
    p_value_robust,
    k_effect_sizes,
    n_studies,
    n_comparisons,
    tau2_level_3,
    tau2_level_2,
    QE,
    QE_df,
    QE_p
  )

final_domain_summary <- domain_specific_model_summary %>%
  mutate(
    analysis = case_when(
      Outcome_Domain == "ER_Total" ~ "Emotion dysregulation",
      Outcome_Domain == "ER_Adaptive" ~ "Adaptive emotion regulation strategies",
      Outcome_Domain == "ER_Maladaptive" ~ "Maladaptive emotion regulation strategies",
      TRUE ~ Outcome_Domain
    )
  ) %>%
  select(
    analysis,
    Outcome_Domain,
    estimate,
    se_robust,
    ci_lb_robust,
    ci_ub_robust,
    p_value_robust,
    k_effect_sizes,
    n_studies,
    n_comparisons,
    tau2_level_3,
    tau2_level_2
  )

cat("\n--- FINAL PRIMARY SUMMARY ---\n")
print(final_primary_summary)

cat("\n--- FINAL DOMAIN SUMMARY ---\n")
print(final_domain_summary, n = nrow(final_domain_summary))


###########################################################
# 7. Final sensitivity summary table
###########################################################

final_sensitivity_summary <- all_sensitivity_results %>%
  select(
    analysis,
    estimate,
    se_robust,
    ci_lb_robust,
    ci_ub_robust,
    p_value_robust,
    k_effect_sizes,
    n_studies,
    n_comparisons,
    tau2_level_3,
    tau2_level_2,
    QE,
    QE_df,
    QE_p
  )

cat("\n--- FINAL SENSITIVITY SUMMARY ---\n")
print(final_sensitivity_summary, n = nrow(final_sensitivity_summary), width = Inf)


###########################################################
# 8. Final moderator summary tables
###########################################################

final_moderator_summary <- moderator_overview %>%
  mutate(
    interpretation = case_when(
      is.na(omnibus_p) ~ "Not estimated or not applicable",
      omnibus_p < .05 ~ "Evidence that moderator explains between-study variation",
      TRUE ~ "No clear evidence that moderator explains between-study variation"
    )
  ) %>%
  select(
    moderator,
    label,
    type,
    status,
    k_effect_sizes,
    n_studies,
    n_comparisons,
    omnibus_stat,
    omnibus_df_num,
    omnibus_df_den,
    omnibus_p,
    interpretation
  )

final_domain_interaction_summary <- domain_interaction_overview %>%
  mutate(
    interpretation = case_when(
      is.na(interaction_p) ~ "Not estimated or not applicable",
      interaction_p < .05 ~ "Evidence that moderator differs across ER domains",
      TRUE ~ "No clear evidence that moderator differs across ER domains"
    )
  ) %>%
  select(
    moderator,
    label,
    type,
    status,
    k_effect_sizes,
    n_studies,
    n_comparisons,
    interaction_stat,
    interaction_df_num,
    interaction_df_den,
    interaction_p,
    interpretation
  )

cat("\n--- FINAL MODERATOR SUMMARY ---\n")
print(final_moderator_summary, n = nrow(final_moderator_summary), width = Inf)

cat("\n--- FINAL DOMAIN-INTERACTION SUMMARY ---\n")
print(final_domain_interaction_summary, n = nrow(final_domain_interaction_summary), width = Inf)


###########################################################
# 9. Final risk-of-bias summary tables
###########################################################

final_rob_summary <- risk_of_bias_overview %>%
  mutate(
    interpretation = case_when(
      is.na(omnibus_p) ~ "Not estimated or not applicable",
      omnibus_p < .05 ~ "Evidence that this ROB2 domain moderates effect sizes",
      TRUE ~ "No clear evidence that this ROB2 domain moderates effect sizes"
    )
  ) %>%
  select(
    moderator,
    label,
    status,
    k_effect_sizes,
    n_studies,
    n_comparisons,
    omnibus_stat,
    omnibus_df_num,
    omnibus_df_den,
    omnibus_p,
    interpretation
  )

final_rob_sensitivity_summary <- risk_of_bias_sensitivity %>%
  select(
    moderator,
    label,
    analysis,
    status,
    estimate,
    se_robust,
    ci_lb_robust,
    ci_ub_robust,
    p_value_robust,
    k_effect_sizes,
    n_studies,
    n_comparisons,
    n_studies_total,
    n_studies_high,
    tau2_level_3,
    tau2_level_2,
    QE,
    QE_df,
    QE_p
  )

cat("\n--- FINAL RISK-OF-BIAS SUMMARY ---\n")
print(final_rob_summary, n = nrow(final_rob_summary), width = Inf)

cat("\n--- FINAL RISK-OF-BIAS SENSITIVITY SUMMARY ---\n")
print(final_rob_sensitivity_summary, n = nrow(final_rob_sensitivity_summary), width = Inf)


###########################################################
# 10. Final ER-MD association summary tables
###########################################################

final_er_md_summary <- er_md_meta_regression_summary %>%
  mutate(
    analysis = "Primary ER-MD meta-regression"
  ) %>%
  select(
    analysis,
    intercept,
    slope_ER,
    intercept_se,
    slope_se,
    intercept_ci_lb,
    intercept_ci_ub,
    slope_ci_lb,
    slope_ci_ub,
    intercept_p,
    slope_p,
    tau2,
    QE,
    QE_df,
    QE_p,
    QM,
    QM_df,
    QM_p,
    k_studies
  )

final_er_md_correlations <- descriptive_correlations_er_md %>%
  mutate(
    analysis = "Primary overlap dataset"
  ) %>%
  select(
    analysis,
    statistic,
    estimate,
    p_value,
    ci_lb,
    ci_ub,
    n_studies
  )

final_er_md_sensitivity <- bind_rows(
  sensitivity_er_md_results %>%
    mutate(sensitivity_type = "Restricted to studies with >1 ER and >1 MD effect"),
  sensitivity_er_md_excluding_extreme %>%
    mutate(sensitivity_type = "Exclude pre-identified extreme overlap studies")
)

cat("\n--- FINAL ER-MD META-REGRESSION SUMMARY ---\n")
print(final_er_md_summary, width = Inf)

cat("\n--- FINAL ER-MD CORRELATION SUMMARY ---\n")
print(final_er_md_correlations, width = Inf)

cat("\n--- FINAL ER-MD SENSITIVITY SUMMARY ---\n")
print(final_er_md_sensitivity, n = nrow(final_er_md_sensitivity), width = Inf)


###########################################################
# 11. Publication-style forest plot using metafor
###########################################################


library(metafor)
library(dplyr)

# ----------------------------------
# 1. Prepare data for plotting
# ----------------------------------

forest_dat <- study_level_er %>%
  filter(yi_agg <= 3) %>%   # exclude extreme displayed studies only
  mutate(
    study_label = paste0("Study ", Study_ID),
    weight_raw = 1 / vi_agg
  ) %>%
  arrange(desc(yi_agg)) %>%
  mutate(
    weight_pct = 100 * weight_raw / sum(weight_raw)
  )

# Pooled robust result from your main summary file
pooled_est   <- overall_er_model_summary$estimate[1]
pooled_ci_lb <- overall_er_model_summary$ci_lb_robust[1]
pooled_ci_ub <- overall_er_model_summary$ci_ub_robust[1]

# Scale point sizes by weight
psize_vals <- 0.8 + 2.2 * sqrt(forest_dat$weight_pct / max(forest_dat$weight_pct))

# Row positions
k_plot <- nrow(forest_dat)
row_spacing <- 1.5

rows_studies <- seq(
  from = k_plot * row_spacing + 1,
  to = 2,
  by = -row_spacing
)

row_pooled <- 1

# ----------------------------------
# 2. Open graphics device
# ----------------------------------


png(
  filename = file.path(output_folder, "forest_plot_study_level_ER_standard_weighted.png"),
  width = 3400,
  height = 4000,
  res = 300
)

par(
  mar = c(5.2, 9, 1.2, 1.8),   # bottom, left, top, right
  cex = 0.85
)


# ----------------------------------
# 3. Draw standard forest plot
# ----------------------------------

forest(
  x = forest_dat$yi_agg,
  vi = forest_dat$vi_agg,
  slab = forest_dat$study_label,
  rows = rows_studies,
  xlim = c(-1.95, 4.00),
  alim = c(-1.5, 3.8),
  at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3.0, 3.5),
  psize = psize_vals,
  efac = 1,
  refline = 0,
  lty = "dashed",
  lwd = 0.9,
  cex = 0.78,
  xlab = "Hedges' g",
  mlab = "",
  header = "Study",
  annotate = FALSE
)

# ----------------------------------
# 4. Add pooled effect diamond
# ----------------------------------

addpoly(
  x = pooled_est,
  ci.lb = pooled_ci_lb,
  ci.ub = pooled_ci_ub,
  row = 1,
  mlab = "Pooled effect",
  col = "black",
  border = "black",
  cex = 0.8
)

# ---------------------------------------------------------
# 5. Optional footnote
# ---------------------------------------------------------

dev.off()

cat("\n--- CLEAN WEIGHTED METAFOR FOREST PLOT CREATED ---\n")
cat("Saved as: forest_plot_study_level_ER_standard_weighted.png\n")

###########################################################
# 12. Manuscript-style key findings table
###########################################################

key_findings_table <- tibble(
  section = c(
    "Primary ER effect",
    "ER domains",
    "Sensitivity analyses",
    "A priori moderators",
    "Moderator-by-domain interactions",
    "Risk of bias domains",
    "ER-MD association"
  ),
  main_result = c(
    "Overall pooled ER effect from the primary multilevel model",
    "Separate pooled effects for ER_Total, ER_Adaptive, and ER_Maladaptive",
    "Primary ER result compared across prespecified sensitivity analyses",
    "Separate moderator models for study-level moderators",
    "Tests of whether moderator effects differ across ER domains",
    "Separate ROB2 domain moderator models and exclusion-of-high-risk sensitivity analyses",
    "Study-level association between aggregated ER and MD effects"
  ),
  primary_output_file = c(
    "overall_ER_model_summary.csv",
    "domain_specific_model_summary.csv",
    "all_sensitivity_results_ER.csv",
    "moderator_overview_ER.csv",
    "domain_interaction_overview_ER.csv",
    "risk_of_bias_overview_ER.csv",
    "er_md_meta_regression_summary.csv"
  )
)

cat("\n--- KEY FINDINGS TABLE ---\n")
print(key_findings_table, width = Inf)


###########################################################
# 13. Save outputs
###########################################################

write_csv(
  final_primary_summary,
  file.path(output_folder, "final_primary_summary_ER.csv")
)

write_csv(
  final_domain_summary,
  file.path(output_folder, "final_domain_summary_ER.csv")
)

write_csv(
  final_sensitivity_summary,
  file.path(output_folder, "final_sensitivity_summary_ER.csv")
)

write_csv(
  final_moderator_summary,
  file.path(output_folder, "final_moderator_summary_ER.csv")
)

write_csv(
  final_domain_interaction_summary,
  file.path(output_folder, "final_domain_interaction_summary_ER.csv")
)

write_csv(
  final_rob_summary,
  file.path(output_folder, "final_risk_of_bias_summary_ER.csv")
)

write_csv(
  final_rob_sensitivity_summary,
  file.path(output_folder, "final_risk_of_bias_sensitivity_summary_ER.csv")
)

write_csv(
  final_er_md_summary,
  file.path(output_folder, "final_er_md_meta_regression_summary.csv")
)

write_csv(
  final_er_md_correlations,
  file.path(output_folder, "final_er_md_correlations_summary.csv")
)

write_csv(
  final_er_md_sensitivity,
  file.path(output_folder, "final_er_md_sensitivity_summary.csv")
)

write_csv(
  key_findings_table,
  file.path(output_folder, "key_findings_table_ER.csv")
)


###########################################################
# 14. Save R objects
###########################################################

save(
  overall_er_model_summary,
  domain_specific_model_summary,
  study_level_er,
  all_sensitivity_results,
  sensitivity_n_table,
  moderator_overview,
  domain_interaction_overview,
  er_md_meta_regression_summary,
  descriptive_correlations_er_md,
  sensitivity_er_md_results,
  sensitivity_er_md_excluding_extreme,
  risk_of_bias_overview,
  risk_of_bias_sensitivity,
  final_primary_summary,
  final_domain_summary,
  final_sensitivity_summary,
  final_moderator_summary,
  final_domain_interaction_summary,
  final_rob_summary,
  final_rob_sensitivity_summary,
  final_er_md_summary,
  final_er_md_correlations,
  final_er_md_sensitivity,
  forest_plot_dat,
  forest_plot_pub,
  key_findings_table,
  file = file.path(output_folder, "final_reporting_ER.RData")
)


###########################################################
# 15. Save text summaries
###########################################################

save_object_summary(
  final_primary_summary,
  file.path(output_folder, "final_primary_summary_ER.txt"),
  "Final primary ER summary"
)

save_object_summary(
  final_sensitivity_summary,
  file.path(output_folder, "final_sensitivity_summary_ER.txt"),
  "Final sensitivity ER summary"
)

save_object_summary(
  final_moderator_summary,
  file.path(output_folder, "final_moderator_summary_ER.txt"),
  "Final moderator ER summary"
)

save_object_summary(
  final_rob_summary,
  file.path(output_folder, "final_risk_of_bias_summary_ER.txt"),
  "Final risk-of-bias ER summary"
)

save_object_summary(
  final_er_md_summary,
  file.path(output_folder, "final_er_md_meta_regression_summary.txt"),
  "Final ER-MD association summary"
)


###########################################################
# 16. Final console output
###########################################################

cat("\n====================================================\n")
cat("Final reporting script finished.\n")
cat("Files saved in: data/derived/final_reporting_ER/\n")
cat("- final_primary_summary_ER.csv\n")
cat("- final_domain_summary_ER.csv\n")
cat("- final_sensitivity_summary_ER.csv\n")
cat("- final_moderator_summary_ER.csv\n")
cat("- final_domain_interaction_summary_ER.csv\n")
cat("- final_risk_of_bias_summary_ER.csv\n")
cat("- final_risk_of_bias_sensitivity_summary_ER.csv\n")
cat("- final_er_md_meta_regression_summary.csv\n")
cat("- final_er_md_correlations_summary.csv\n")
cat("- final_er_md_sensitivity_summary.csv\n")
cat("- key_findings_table_ER.csv\n")
cat("- forest_plot_study_level_ER_publication.png\n")
cat("- final_reporting_ER.RData\n")
cat("====================================================\n")