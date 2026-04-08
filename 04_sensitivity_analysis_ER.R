############################################################
# Script: 04_sensitivity_analysis_ER.R
# Project: ER meta-analysis_final
# Purpose:
#   Run sensitivity analyses for the primary emotion
#   regulation (ER) meta-analysis.
#
# Sensitivity analyses included:
#   1. Alternative assumed pre-post correlations:
#      r = .30, .50, .70
#   2. Exclusion of very extreme effect sizes:
#      |Hedges' g| > 3
#   3. Exclusion of flagged small-standardiser cases:
#      pooled pretest SD < 0.10 / 0.20 / 0.50
#
# Main modelling approach:
#   - Multilevel random-effects meta-analysis
#   - Effect sizes nested within studies, with a lower-level
#     random effect for unique effect/comparison rows within study
#   - REML estimation
#   - Cluster-robust inference at the study level (CR2)
#   - Robust confidence intervals based on Satterthwaite-
#     adjusted degrees of freedom
#
# Input:
#   data/derived/effect_sizes/er_meta_effect_sizes_r30.csv
#   data/derived/effect_sizes/er_meta_effect_sizes_r50.csv
#   data/derived/effect_sizes/er_meta_effect_sizes_r70.csv
#
# Output:
#   data/derived/sensitivity_ER/
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
  "metafor",
  "clubSandwich",
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
library(metafor)
library(clubSandwich)
library(ggplot2)


############################
# 2. Define file locations
############################

input_r30 <- "data/derived/effect_sizes/er_meta_effect_sizes_r30.csv"
input_r50 <- "data/derived/effect_sizes/er_meta_effect_sizes_r50.csv"
input_r70 <- "data/derived/effect_sizes/er_meta_effect_sizes_r70.csv"

output_folder <- "data/derived/sensitivity_ER"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


############################
# 3. Read datasets
############################

dat_r30 <- read_csv(input_r30, show_col_types = FALSE)
dat_r50 <- read_csv(input_r50, show_col_types = FALSE)
dat_r70 <- read_csv(input_r70, show_col_types = FALSE)

cat("\n----------------------------------------\n")
cat("ER sensitivity analysis script started\n")
cat("----------------------------------------\n")
cat("Rows read (r = .30):", nrow(dat_r30), "\n")
cat("Rows read (r = .50):", nrow(dat_r50), "\n")
cat("Rows read (r = .70):", nrow(dat_r70), "\n")


###########################################################
# 4. Helper functions
###########################################################

prepare_er_data <- function(data) {
  data %>%
    filter(analysis_set == "ER") %>%
    mutate(
      Study_ID = as.factor(Study_ID),
      Comparison_ID = as.factor(Comparison_ID),
      Outcome_Domain = factor(
        Outcome_Domain,
        levels = c("ER_Total", "ER_Adaptive", "ER_Maladaptive")
      )
    )
}

fit_multilevel_model <- function(data) {
  rma.mv(
    yi = yi,
    V = vi,
    random = ~ 1 | Study_ID/Comparison_ID,
    data = data,
    method = "REML",
    test = "t"
  )
}

get_robust_tests <- function(model, data) {
  coef_test(
    model,
    vcov = "CR2",
    cluster = data$Study_ID
  )
}

extract_model_summary <- function(model, data, analysis_label) {
  robust <- get_robust_tests(model, data)
  robust_row <- as.data.frame(robust)[1, ]
  t_crit <- qt(0.975, df = robust_row$df_Satt)
  
  tibble(
    analysis = analysis_label,
    estimate = as.numeric(coef(model)[1]),
    se_robust = robust_row$SE,
    t_robust = robust_row$tstat,
    df_robust = robust_row$df_Satt,
    p_value_robust = robust_row$p_Satt,
    ci_lb_robust = robust_row$beta - t_crit * robust_row$SE,
    ci_ub_robust = robust_row$beta + t_crit * robust_row$SE,
    k_effect_sizes = model$k,
    n_studies = n_distinct(data$Study_ID),
    n_comparisons = n_distinct(data$Comparison_ID),
    tau2_level_3 = model$sigma2[1],
    tau2_level_2 = model$sigma2[2],
    QE = model$QE,
    QE_df = model$k - model$p,
    QE_p = model$QEp
  )
}

save_model_output <- function(object, file_path, title_text = NULL) {
  sink(file_path)
  if (!is.null(title_text)) {
    cat(title_text, "\n")
    cat(paste(rep("=", nchar(title_text)), collapse = ""), "\n\n")
  }
  print(object)
  sink()
}


###########################################################
# 5. Prepare ER datasets
###########################################################

dat_r30_er <- prepare_er_data(dat_r30)
dat_r50_er <- prepare_er_data(dat_r50)
dat_r70_er <- prepare_er_data(dat_r70)

cat("\n--- ER DATASET COUNTS ---\n")
cat("ER rows (r = .30):", nrow(dat_r30_er), "\n")
cat("ER rows (r = .50):", nrow(dat_r50_er), "\n")
cat("ER rows (r = .70):", nrow(dat_r70_er), "\n")


###########################################################
# 6. Sensitivity analysis 1:
#    Alternative pre-post correlations
###########################################################

model_r30 <- fit_multilevel_model(dat_r30_er)
model_r50 <- fit_multilevel_model(dat_r50_er)
model_r70 <- fit_multilevel_model(dat_r70_er)

robust_r30 <- get_robust_tests(model_r30, dat_r30_er)
robust_r50 <- get_robust_tests(model_r50, dat_r50_er)
robust_r70 <- get_robust_tests(model_r70, dat_r70_er)

results_r <- bind_rows(
  extract_model_summary(model_r30, dat_r30_er, "Primary model with assumed r = .30"),
  extract_model_summary(model_r50, dat_r50_er, "Primary model with assumed r = .50"),
  extract_model_summary(model_r70, dat_r70_er, "Primary model with assumed r = .70")
)

cat("\n--- SENSITIVITY RESULTS: DIFFERENT PRE-POST CORRELATIONS ---\n")
print(results_r)


###########################################################
# 7. Sensitivity analysis 2:
#    Excluding very extreme effect sizes |g| > 3
###########################################################

dat_r50_no_abs_gt_3 <- dat_r50_er %>%
  filter(abs(yi) <= 3)

model_no_abs_gt_3 <- fit_multilevel_model(dat_r50_no_abs_gt_3)
robust_no_abs_gt_3 <- get_robust_tests(model_no_abs_gt_3, dat_r50_no_abs_gt_3)

results_abs_gt_3 <- extract_model_summary(
  model_no_abs_gt_3,
  dat_r50_no_abs_gt_3,
  "Primary model excluding effect sizes with |g| > 3"
)

cat("\n--- SENSITIVITY RESULTS: EXCLUDING |g| > 3 ---\n")
print(results_abs_gt_3)


###########################################################
# 8. Sensitivity analysis 3:
#    Excluding rows with very small pooled baseline SD
###########################################################

dat_r50_no_sd_lt_010 <- dat_r50_er %>%
  filter(is.na(sd_pooled_pre) | sd_pooled_pre >= 0.10)

dat_r50_no_sd_lt_020 <- dat_r50_er %>%
  filter(is.na(sd_pooled_pre) | sd_pooled_pre >= 0.20)

dat_r50_no_sd_lt_050 <- dat_r50_er %>%
  filter(is.na(sd_pooled_pre) | sd_pooled_pre >= 0.50)

model_no_sd_lt_010 <- fit_multilevel_model(dat_r50_no_sd_lt_010)
model_no_sd_lt_020 <- fit_multilevel_model(dat_r50_no_sd_lt_020)
model_no_sd_lt_050 <- fit_multilevel_model(dat_r50_no_sd_lt_050)

robust_no_sd_lt_010 <- get_robust_tests(model_no_sd_lt_010, dat_r50_no_sd_lt_010)
robust_no_sd_lt_020 <- get_robust_tests(model_no_sd_lt_020, dat_r50_no_sd_lt_020)
robust_no_sd_lt_050 <- get_robust_tests(model_no_sd_lt_050, dat_r50_no_sd_lt_050)

results_small_sd <- bind_rows(
  extract_model_summary(
    model_no_sd_lt_010,
    dat_r50_no_sd_lt_010,
    "Primary model excluding pooled pretest SD < 0.10"
  ),
  extract_model_summary(
    model_no_sd_lt_020,
    dat_r50_no_sd_lt_020,
    "Primary model excluding pooled pretest SD < 0.20"
  ),
  extract_model_summary(
    model_no_sd_lt_050,
    dat_r50_no_sd_lt_050,
    "Primary model excluding pooled pretest SD < 0.50"
  )
)

cat("\n--- SENSITIVITY RESULTS: EXCLUDING SMALL POOLED PRETEST SD CASES ---\n")
print(results_small_sd)


###########################################################
# 9. Combine all sensitivity results
###########################################################

main_result <- extract_model_summary(
  model_r50,
  dat_r50_er,
  "Primary model with assumed r = .50 (reference)"
)

all_sensitivity_results <- bind_rows(
  main_result,
  results_r %>% filter(analysis != "Primary model with assumed r = .50"),
  results_abs_gt_3,
  results_small_sd
)

cat("\n--- ALL SENSITIVITY RESULTS ---\n")
print(all_sensitivity_results)


###########################################################
# 10. Create manuscript-ready table and comparison plot
###########################################################

sensitivity_results_pub <- all_sensitivity_results %>%
  mutate(
    estimate = round(estimate, 2),
    se_robust = round(se_robust, 2),
    ci_lb_robust = round(ci_lb_robust, 2),
    ci_ub_robust = round(ci_ub_robust, 2),
    p_value_robust = ifelse(
      p_value_robust < .001,
      "< .001",
      sprintf("%.3f", p_value_robust)
    ),
    `95% CI` = paste0("[", ci_lb_robust, ", ", ci_ub_robust, "]")
  ) %>%
  select(
    analysis,
    estimate,
    `95% CI`,
    se_robust,
    p_value_robust,
    k_effect_sizes,
    n_studies,
    tau2_level_3,
    tau2_level_2
  )

cat("\n--- MANUSCRIPT-READY SENSITIVITY TABLE ---\n")
print(sensitivity_results_pub)

sensitivity_plot <- ggplot(
  all_sensitivity_results,
  aes(
    x = reorder(analysis, estimate),
    y = estimate
  )
) +
  geom_point() +
  geom_errorbar(
    aes(ymin = ci_lb_robust, ymax = ci_ub_robust),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    title = "Sensitivity analyses for the primary ER model",
    x = NULL,
    y = "Pooled Hedges' g (robust 95% CI)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(output_folder, "sensitivity_results_ER_plot.png"),
  plot = sensitivity_plot,
  width = 10,
  height = 6,
  dpi = 300
)


###########################################################
# 11. Save outputs
###########################################################

write_csv(
  results_r,
  file.path(output_folder, "sensitivity_results_prepost_correlation.csv")
)

write_csv(
  results_abs_gt_3,
  file.path(output_folder, "sensitivity_results_excluding_abs_gt_3.csv")
)

write_csv(
  results_small_sd,
  file.path(output_folder, "sensitivity_results_excluding_small_sd.csv")
)

write_csv(
  all_sensitivity_results,
  file.path(output_folder, "all_sensitivity_results_ER.csv")
)

write_csv(
  sensitivity_results_pub,
  file.path(output_folder, "sensitivity_results_manuscript_ready.csv")
)

write_csv(
  dat_r50_no_abs_gt_3,
  file.path(output_folder, "dataset_ER_excluding_abs_gt_3.csv")
)

write_csv(
  dat_r50_no_sd_lt_010,
  file.path(output_folder, "dataset_ER_excluding_sd_pooled_pre_lt_010.csv")
)

write_csv(
  dat_r50_no_sd_lt_020,
  file.path(output_folder, "dataset_ER_excluding_sd_pooled_pre_lt_020.csv")
)

write_csv(
  dat_r50_no_sd_lt_050,
  file.path(output_folder, "dataset_ER_excluding_sd_pooled_pre_lt_050.csv")
)

# Model-based summaries kept for transparency
save_model_output(
  summary(model_r30),
  file.path(output_folder, "model_r30_summary.txt"),
  "Primary ER model with assumed pre-post correlation r = .30"
)

save_model_output(
  summary(model_r50),
  file.path(output_folder, "model_r50_summary.txt"),
  "Primary ER model with assumed pre-post correlation r = .50"
)

save_model_output(
  summary(model_r70),
  file.path(output_folder, "model_r70_summary.txt"),
  "Primary ER model with assumed pre-post correlation r = .70"
)

save_model_output(
  summary(model_no_abs_gt_3),
  file.path(output_folder, "model_excluding_abs_gt_3_summary.txt"),
  "Primary ER model excluding effect sizes with |g| > 3"
)

save_model_output(
  summary(model_no_sd_lt_010),
  file.path(output_folder, "model_excluding_sd_lt_010_summary.txt"),
  "Primary ER model excluding pooled pretest SD < 0.10"
)

save_model_output(
  summary(model_no_sd_lt_020),
  file.path(output_folder, "model_excluding_sd_lt_020_summary.txt"),
  "Primary ER model excluding pooled pretest SD < 0.20"
)

save_model_output(
  summary(model_no_sd_lt_050),
  file.path(output_folder, "model_excluding_sd_lt_050_summary.txt"),
  "Primary ER model excluding pooled pretest SD < 0.50"
)

# Robust outputs for manuscript reporting
save_model_output(
  robust_r30,
  file.path(output_folder, "model_r30_robust_tests.txt"),
  "Cluster-robust tests for primary ER model with assumed pre-post correlation r = .30"
)

save_model_output(
  robust_r50,
  file.path(output_folder, "model_r50_robust_tests.txt"),
  "Cluster-robust tests for primary ER model with assumed pre-post correlation r = .50"
)

save_model_output(
  robust_r70,
  file.path(output_folder, "model_r70_robust_tests.txt"),
  "Cluster-robust tests for primary ER model with assumed pre-post correlation r = .70"
)

save_model_output(
  robust_no_abs_gt_3,
  file.path(output_folder, "model_excluding_abs_gt_3_robust_tests.txt"),
  "Cluster-robust tests for primary ER model excluding effect sizes with |g| > 3"
)

save_model_output(
  robust_no_sd_lt_010,
  file.path(output_folder, "model_excluding_sd_lt_010_robust_tests.txt"),
  "Cluster-robust tests for primary ER model excluding pooled pretest SD < 0.10"
)

save_model_output(
  robust_no_sd_lt_020,
  file.path(output_folder, "model_excluding_sd_lt_020_robust_tests.txt"),
  "Cluster-robust tests for primary ER model excluding pooled pretest SD < 0.20"
)

save_model_output(
  robust_no_sd_lt_050,
  file.path(output_folder, "model_excluding_sd_lt_050_robust_tests.txt"),
  "Cluster-robust tests for primary ER model excluding pooled pretest SD < 0.50"
)

write_csv(
  as_tibble(robust_r30, rownames = "term"),
  file.path(output_folder, "model_r30_robust_tests.csv")
)

write_csv(
  as_tibble(robust_r50, rownames = "term"),
  file.path(output_folder, "model_r50_robust_tests.csv")
)

write_csv(
  as_tibble(robust_r70, rownames = "term"),
  file.path(output_folder, "model_r70_robust_tests.csv")
)

write_csv(
  as_tibble(robust_no_abs_gt_3, rownames = "term"),
  file.path(output_folder, "model_excluding_abs_gt_3_robust_tests.csv")
)

write_csv(
  as_tibble(robust_no_sd_lt_010, rownames = "term"),
  file.path(output_folder, "model_excluding_sd_lt_010_robust_tests.csv")
)

write_csv(
  as_tibble(robust_no_sd_lt_020, rownames = "term"),
  file.path(output_folder, "model_excluding_sd_lt_020_robust_tests.csv")
)

write_csv(
  as_tibble(robust_no_sd_lt_050, rownames = "term"),
  file.path(output_folder, "model_excluding_sd_lt_050_robust_tests.csv")
)

save(
  dat_r30_er,
  dat_r50_er,
  dat_r70_er,
  dat_r50_no_abs_gt_3,
  dat_r50_no_sd_lt_010,
  dat_r50_no_sd_lt_020,
  dat_r50_no_sd_lt_050,
  model_r30,
  model_r50,
  model_r70,
  model_no_abs_gt_3,
  model_no_sd_lt_010,
  model_no_sd_lt_020,
  model_no_sd_lt_050,
  robust_r30,
  robust_r50,
  robust_r70,
  robust_no_abs_gt_3,
  robust_no_sd_lt_010,
  robust_no_sd_lt_020,
  robust_no_sd_lt_050,
  all_sensitivity_results,
  sensitivity_results_pub,
  file = file.path(output_folder, "sensitivity_analysis_ER.RData")
)


###########################################################
# 12. Create and save sensitivity sample-size table
###########################################################

sensitivity_n_table <- tibble(
  analysis = c(
    "Primary ER dataset",
    "Exclude |g| > 3",
    "Exclude sd_pooled_pre < 0.10",
    "Exclude sd_pooled_pre < 0.20",
    "Exclude sd_pooled_pre < 0.50"
  ),
  k_effect_sizes = c(
    nrow(dat_r50_er),
    nrow(dat_r50_no_abs_gt_3),
    nrow(dat_r50_no_sd_lt_010),
    nrow(dat_r50_no_sd_lt_020),
    nrow(dat_r50_no_sd_lt_050)
  ),
  n_studies = c(
    dplyr::n_distinct(dat_r50_er$Study_ID),
    dplyr::n_distinct(dat_r50_no_abs_gt_3$Study_ID),
    dplyr::n_distinct(dat_r50_no_sd_lt_010$Study_ID),
    dplyr::n_distinct(dat_r50_no_sd_lt_020$Study_ID),
    dplyr::n_distinct(dat_r50_no_sd_lt_050$Study_ID)
  )
)

print(sensitivity_n_table)

write_csv(
  sensitivity_n_table,
  file.path(output_folder, "sensitivity_sample_sizes_ER.csv")
)


###########################################################
# 13. Final console output
###########################################################

cat("\n====================================================\n")
cat("ER sensitivity analysis script finished.\n")
cat("Files saved in: data/derived/sensitivity_ER/\n")
cat("- sensitivity_results_prepost_correlation.csv\n")
cat("- sensitivity_results_excluding_abs_gt_3.csv\n")
cat("- sensitivity_results_excluding_small_sd.csv\n")
cat("- all_sensitivity_results_ER.csv\n")
cat("- sensitivity_results_manuscript_ready.csv\n")
cat("- dataset_ER_excluding_abs_gt_3.csv\n")
cat("- dataset_ER_excluding_sd_pooled_pre_lt_010.csv\n")
cat("- dataset_ER_excluding_sd_pooled_pre_lt_020.csv\n")
cat("- dataset_ER_excluding_sd_pooled_pre_lt_050.csv\n")
cat("- model_r30_summary.txt\n")
cat("- model_r50_summary.txt\n")
cat("- model_r70_summary.txt\n")
cat("- model_excluding_abs_gt_3_summary.txt\n")
cat("- model_excluding_sd_lt_010_summary.txt\n")
cat("- model_excluding_sd_lt_020_summary.txt\n")
cat("- model_excluding_sd_lt_050_summary.txt\n")
cat("- model_r30_robust_tests.txt\n")
cat("- model_r50_robust_tests.txt\n")
cat("- model_r70_robust_tests.txt\n")
cat("- model_excluding_abs_gt_3_robust_tests.txt\n")
cat("- model_excluding_sd_lt_010_robust_tests.txt\n")
cat("- model_excluding_sd_lt_020_robust_tests.txt\n")
cat("- model_excluding_sd_lt_050_robust_tests.txt\n")
cat("- model_r30_robust_tests.csv\n")
cat("- model_r50_robust_tests.csv\n")
cat("- model_r70_robust_tests.csv\n")
cat("- model_excluding_abs_gt_3_robust_tests.csv\n")
cat("- model_excluding_sd_lt_010_robust_tests.csv\n")
cat("- model_excluding_sd_lt_020_robust_tests.csv\n")
cat("- model_excluding_sd_lt_050_robust_tests.csv\n")
cat("- sensitivity_results_ER_plot.png\n")
cat("- sensitivity_sample_sizes_ER.csv\n")
cat("- sensitivity_analysis_ER.RData\n")
cat("====================================================\n")