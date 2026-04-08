############################################################
# Script: 03_meta_analysis_models.R
# Project: ER meta-analysis_final
# Purpose:
#   Fit the main multilevel meta-analytic models for the
#   meta-analysis of mindfulness-based interventions (MBIs)
#   on emotion regulation.
#
# Main analysis dataset:
#   - Primary effect sizes computed with assumed pre-post
#     correlation r = .50
#   - Emotion regulation outcomes only
#
# Main modelling approach:
#   - Multilevel random-effects meta-analysis
#   - Effect sizes nested within studies, with a lower-level
#     random effect for unique comparisons/effects within study
#   - Models estimated with REML
#   - Cluster-robust inference at the study level
#
# Additional analyses:
#   - Overall ER model
#   - Moderator model by ER domain
#   - Domain-specific pooled models
#   - Study-level aggregation for descriptive summaries
#   - Funnel plot and Egger's regression test on aggregated
#     study-level effects
#   - Leave-one-study-out sensitivity analysis
#
# Input:
#   data/derived/effect_sizes/er_meta_effect_sizes_r50.csv
#
# Output:
#   data/derived/meta_models/
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

input_file <- "data/derived/effect_sizes/er_meta_effect_sizes_r50.csv"
output_folder <- "data/derived/meta_models"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


############################
# 3. Read effect size dataset
############################

dat <- read_csv(input_file, show_col_types = FALSE)

cat("\n----------------------------------------\n")
cat("Meta-analysis model script started\n")
cat("----------------------------------------\n")
cat("Rows read:", nrow(dat), "\n")
cat("Columns read:", ncol(dat), "\n")


###########################################################
# 4. Check required variables are present
###########################################################

required_vars <- c(
  "Effect_ID",
  "Study_ID",
  "Comparison_ID",
  "Outcome_Domain",
  "Outcome_Measure",
  "yi",
  "vi",
  "analysis_set"
)

missing_required_vars <- setdiff(required_vars, names(dat))

if (length(missing_required_vars) > 0) {
  stop(
    "The following required columns are missing:\n",
    paste(missing_required_vars, collapse = ", ")
  )
}

cat("All required columns are present.\n")


###########################################################
# 5. Prepare analysis dataset
###########################################################

# Main inferential models focus on emotion regulation outcomes.
# Mental distress outcomes will be used later for the separate
# association analyses.

dat_er <- dat %>%
  filter(analysis_set == "ER") %>%
  mutate(
    Study_ID = as.factor(Study_ID),
    Comparison_ID = as.factor(Comparison_ID),
    Outcome_Domain = factor(
      Outcome_Domain,
      levels = c("ER_Total", "ER_Adaptive", "ER_Maladaptive")
    )
  )

cat("\n--- MAIN ANALYSIS DATASET ---\n")
cat("ER rows:", nrow(dat_er), "\n")
cat("Unique studies:", n_distinct(dat_er$Study_ID), "\n")
cat("Unique comparisons:", n_distinct(dat_er$Comparison_ID), "\n")
cat("Comparison_ID unique within ER dataset:",
    n_distinct(dat_er$Comparison_ID) == nrow(dat_er), "\n")
cat("Outcome domains:\n")
print(table(dat_er$Outcome_Domain, useNA = "ifany"))

if (nrow(dat_er) == 0) {
  stop("No ER rows found in the primary effect size dataset.")
}

if (sum(is.na(dat_er$yi)) > 0) {
  stop("ER dataset contains missing yi values.")
}

if (sum(is.na(dat_er$vi)) > 0) {
  stop("ER dataset contains missing vi values.")
}

if (sum(dat_er$vi <= 0, na.rm = TRUE) > 0) {
  stop("ER dataset contains non-positive vi values.")
}


###########################################################
# 6. Helper function to save model summaries to text files
###########################################################

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
# 7. Fit main multilevel random-effects model
###########################################################

# Primary model:
# Effect sizes are nested within studies, with a lower-level
# random effect for unique comparisons/effects within studies.
# In this dataset, Comparison_ID is effectively unique per row,
# so the lower level functions as an effect/comparison-within-
# study random effect rather than a shared comparison cluster.

model_er_overall <- rma.mv(
  yi = yi,
  V = vi,
  random = ~ 1 | Study_ID/Comparison_ID,
  data = dat_er,
  method = "REML",
  test = "t"
)

cat("\n--- OVERALL MULTILEVEL MODEL (ER) ---\n")
print(summary(model_er_overall))


###########################################################
# 8. Cluster-robust inference for the main model
###########################################################

# Robust variance estimation is used to obtain study-clustered
# standard errors and tests that are less sensitive to residual
# dependence structures within studies.

model_er_overall_robust <- coef_test(
  model_er_overall,
  vcov = "CR2",
  cluster = dat_er$Study_ID
)

cat("\n--- ROBUST TESTS: OVERALL MULTILEVEL MODEL (ER) ---\n")
print(model_er_overall_robust)


###########################################################
# 9. Variance components and model information
###########################################################

variance_components <- tibble(
  component = c("study_level", "effect_or_comparison_level"),
  sigma2 = model_er_overall$sigma2
)

model_fit_summary <- tibble(
  k_effect_sizes = model_er_overall$k,
  n_studies = n_distinct(dat_er$Study_ID),
  n_comparisons = n_distinct(dat_er$Comparison_ID),
  logLik = as.numeric(logLik(model_er_overall)),
  AIC = AIC(model_er_overall),
  BIC = BIC(model_er_overall),
  QE = model_er_overall$QE,
  QE_df = model_er_overall$k - model_er_overall$p,
  QE_p = model_er_overall$QEp
)

cat("\n--- VARIANCE COMPONENTS ---\n")
print(variance_components)

cat("\n--- MODEL FIT SUMMARY ---\n")
print(model_fit_summary)


###########################################################
# 10. Moderator model by emotion regulation domain
###########################################################

# This model examines whether pooled effects differ across
# ER_Total, ER_Adaptive, and ER_Maladaptive.

model_er_domain <- rma.mv(
  yi = yi,
  V = vi,
  mods = ~ Outcome_Domain,
  random = ~ 1 | Study_ID/Comparison_ID,
  data = dat_er,
  method = "REML",
  test = "t"
)

cat("\n--- MODERATOR MODEL: OUTCOME DOMAIN ---\n")
print(summary(model_er_domain))

model_er_domain_robust <- coef_test(
  model_er_domain,
  vcov = "CR2",
  cluster = dat_er$Study_ID
)

cat("\n--- ROBUST TESTS: MODERATOR MODEL (OUTCOME DOMAIN) ---\n")
print(model_er_domain_robust)

domain_coef_idx <- grep("^Outcome_Domain", names(coef(model_er_domain)))

domain_omnibus_robust <- Wald_test(
  model_er_domain,
  constraints = constrain_zero(domain_coef_idx),
  vcov = "CR2",
  cluster = dat_er$Study_ID
)

cat("\n--- ROBUST OMNIBUS TEST OF DOMAIN MODERATOR ---\n")
print(domain_omnibus_robust)

###########################################################
# 11. Domain-specific pooled models
###########################################################

fit_domain_model <- function(domain_name) {
  
  dat_sub <- dat_er %>%
    filter(Outcome_Domain == domain_name)
  
  model_sub <- rma.mv(
    yi = yi,
    V = vi,
    random = ~ 1 | Study_ID/Comparison_ID,
    data = dat_sub,
    method = "REML",
    test = "t"
  )
  
  robust_sub <- coef_test(
    model_sub,
    vcov = "CR2",
    cluster = dat_sub$Study_ID
  )
  
  list(
    data = dat_sub,
    model = model_sub,
    robust = robust_sub
  )
}

domain_total <- fit_domain_model("ER_Total")
domain_adaptive <- fit_domain_model("ER_Adaptive")
domain_maladaptive <- fit_domain_model("ER_Maladaptive")

cat("\n--- DOMAIN-SPECIFIC MODEL: ER_TOTAL ---\n")
print(summary(domain_total$model))
print(domain_total$robust)

cat("\n--- DOMAIN-SPECIFIC MODEL: ER_ADAPTIVE ---\n")
print(summary(domain_adaptive$model))
print(domain_adaptive$robust)

cat("\n--- DOMAIN-SPECIFIC MODEL: ER_MALADAPTIVE ---\n")
print(summary(domain_maladaptive$model))
print(domain_maladaptive$robust)


###########################################################
# 12. Create study-level aggregated dataset
###########################################################
# These files are useful for descriptive reporting and for
# exploratory small-study-effect analyses on aggregated
# study-level effects rather than all dependent effects.

aggregate_study_effects <- function(data, domain_name = NULL) {
  
  dat_use <- data
  
  if (!is.null(domain_name)) {
    dat_use <- dat_use %>%
      filter(Outcome_Domain == domain_name)
  }
  
  dat_use %>%
    group_by(Study_ID) %>%
    summarise(
      k_effects = n(),
      yi_agg = weighted.mean(yi, w = 1 / vi, na.rm = TRUE),
      vi_agg = 1 / sum(1 / vi, na.rm = TRUE),
      sei_agg = sqrt(vi_agg),
      Outcome_Domains = paste(sort(unique(Outcome_Domain)), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lb = yi_agg - 1.96 * sei_agg,
      ci_ub = yi_agg + 1.96 * sei_agg
    ) %>%
    arrange(as.numeric(as.character(Study_ID)))
}

study_level_er <- aggregate_study_effects(dat_er)
study_level_er_total <- aggregate_study_effects(dat_er, "ER_Total")
study_level_er_adaptive <- aggregate_study_effects(dat_er, "ER_Adaptive")
study_level_er_maladaptive <- aggregate_study_effects(dat_er, "ER_Maladaptive")

cat("\n--- STUDY-LEVEL AGGREGATED ER EFFECTS ---\n")
cat("Overall pooled ER studies:", nrow(study_level_er), "\n")
cat("ER_Total studies:", nrow(study_level_er_total), "\n")
cat("ER_Adaptive studies:", nrow(study_level_er_adaptive), "\n")
cat("ER_Maladaptive studies:", nrow(study_level_er_maladaptive), "\n")

###########################################################
# 13. Small-study effect analyses
###########################################################

# Exploratory small-study-effect analyses are run on aggregated
# study-level effects. These results should be interpreted
# cautiously because standard funnel-plot and Egger approaches
# are primarily designed for independent effect sizes and do not
# fully accommodate multilevel dependence.

run_egger_analysis <- function(data, analysis_label) {
  
  model <- rma(
    yi = yi_agg,
    vi = vi_agg,
    data = data,
    method = "REML",
    test = "t"
  )
  
  test <- regtest(
    model,
    model = "rma",
    predictor = "sei"
  )
  
  cat("\n--- EXPLORATORY EGGER MODEL:", analysis_label, "---\n")
  print(summary(model))
  
  cat("\n--- EXPLORATORY EGGER TEST:", analysis_label, "---\n")
  print(test)
  
  list(
    model = model,
    test = test
  )
}

egger_overall <- run_egger_analysis(
  study_level_er,
  "Overall ER (all domains pooled)"
)

egger_total <- run_egger_analysis(
  study_level_er_total,
  "ER_Total"
)

egger_adaptive <- run_egger_analysis(
  study_level_er_adaptive,
  "ER_Adaptive"
)

egger_maladaptive <- run_egger_analysis(
  study_level_er_maladaptive,
  "ER_Maladaptive"
)
############################################################
# 14. Leave-one-study-out sensitivity analysis
###########################################################

# This is a transparent influence analysis at the study level.
# It refits the main model while excluding one study at a time.
# Inference is based on cluster-robust variance estimation
# (CR2) with Satterthwaite-adjusted degrees of freedom.

loo_results <- vector("list", length = n_distinct(dat_er$Study_ID))
study_ids <- sort(unique(dat_er$Study_ID))

for (i in seq_along(study_ids)) {
  
  current_study <- study_ids[i]
  
  dat_minus_one <- dat_er %>%
    filter(Study_ID != current_study)
  
  model_minus_one <- rma.mv(
    yi = yi,
    V = vi,
    random = ~ 1 | Study_ID/Comparison_ID,
    data = dat_minus_one,
    method = "REML",
    test = "t"
  )
  
  robust_minus_one <- coef_test(
    model_minus_one,
    vcov = "CR2",
    cluster = dat_minus_one$Study_ID
  )
  
  robust_row <- as.data.frame(robust_minus_one)[1, ]
  t_crit <- qt(0.975, df = robust_row$df_Satt)
  
  loo_results[[i]] <- tibble(
    omitted_study = as.character(current_study),
    estimate = as.numeric(coef(model_minus_one)[1]),
    se_robust = robust_row$SE,
    df_robust = robust_row$df_Satt,
    p_value_robust = robust_row$p_Satt,
    ci_lb_robust = robust_row$beta - t_crit * robust_row$SE,
    ci_ub_robust = robust_row$beta + t_crit * robust_row$SE,
    tau2_level_3 = model_minus_one$sigma2[1],
    tau2_level_2 = model_minus_one$sigma2[2],
    k_remaining = model_minus_one$k
  )
}

full_estimate <- as.numeric(coef(model_er_overall)[1])

loo_results <- bind_rows(loo_results) %>%
  mutate(
    delta_estimate = estimate - full_estimate
  ) %>%
  arrange(desc(abs(delta_estimate)))

cat("\n--- LEAVE-ONE-STUDY-OUT RESULTS ---\n")
print(loo_results, n = nrow(loo_results))


###########################################################
# Helper function to extract robust inference for manuscript
# summary tables
###########################################################

extract_robust_row <- function(robust_table, term_name) {
  robust_df <- as.data.frame(robust_table)
  robust_df$term <- rownames(robust_df)
  row <- robust_df %>% filter(term == term_name)
  
  if (nrow(row) != 1) {
    stop("Could not uniquely identify robust row for term: ", term_name)
  }
  
  t_crit <- qt(0.975, df = row$df_Satt)
  
  tibble(
    term = term_name,
    estimate_robust = row$beta,
    se_robust = row$SE,
    df_robust = row$df_Satt,
    p_value_robust = row$p_Satt,
    ci_lb_robust = row$beta - t_crit * row$SE,
    ci_ub_robust = row$beta + t_crit * row$SE
  )
}
###########################################################
# 15. Create publication-ready summary tables
###########################################################

overall_robust_row <- extract_robust_row(model_er_overall_robust, "intrcpt")

overall_summary_table <- tibble(
  model = "Overall ER multilevel model",
  estimate = as.numeric(coef(model_er_overall)[1]),
  se_robust = overall_robust_row$se_robust,
  ci_lb_robust = overall_robust_row$ci_lb_robust,
  ci_ub_robust = overall_robust_row$ci_ub_robust,
  df_robust = overall_robust_row$df_robust,
  p_value_robust = overall_robust_row$p_value_robust,
  k_effect_sizes = model_er_overall$k,
  n_studies = n_distinct(dat_er$Study_ID),
  n_comparisons = n_distinct(dat_er$Comparison_ID),
  tau2_level_3 = model_er_overall$sigma2[1],
  tau2_level_2 = model_er_overall$sigma2[2],
  QE = model_er_overall$QE,
  QE_df = model_er_overall$k - model_er_overall$p,
  QE_p = model_er_overall$QEp
)

domain_total_robust_row <- extract_robust_row(domain_total$robust, "intrcpt")
domain_adaptive_robust_row <- extract_robust_row(domain_adaptive$robust, "intrcpt")
domain_maladaptive_robust_row <- extract_robust_row(domain_maladaptive$robust, "intrcpt")

domain_summary_table <- bind_rows(
  tibble(
    Outcome_Domain = "ER_Total",
    estimate = as.numeric(coef(domain_total$model)[1]),
    se_robust = domain_total_robust_row$se_robust,
    ci_lb_robust = domain_total_robust_row$ci_lb_robust,
    ci_ub_robust = domain_total_robust_row$ci_ub_robust,
    df_robust = domain_total_robust_row$df_robust,
    p_value_robust = domain_total_robust_row$p_value_robust,
    k_effect_sizes = domain_total$model$k,
    n_studies = n_distinct(domain_total$data$Study_ID),
    n_comparisons = n_distinct(domain_total$data$Comparison_ID),
    tau2_level_3 = domain_total$model$sigma2[1],
    tau2_level_2 = domain_total$model$sigma2[2]
  ),
  tibble(
    Outcome_Domain = "ER_Adaptive",
    estimate = as.numeric(coef(domain_adaptive$model)[1]),
    se_robust = domain_adaptive_robust_row$se_robust,
    ci_lb_robust = domain_adaptive_robust_row$ci_lb_robust,
    ci_ub_robust = domain_adaptive_robust_row$ci_ub_robust,
    df_robust = domain_adaptive_robust_row$df_robust,
    p_value_robust = domain_adaptive_robust_row$p_value_robust,
    k_effect_sizes = domain_adaptive$model$k,
    n_studies = n_distinct(domain_adaptive$data$Study_ID),
    n_comparisons = n_distinct(domain_adaptive$data$Comparison_ID),
    tau2_level_3 = domain_adaptive$model$sigma2[1],
    tau2_level_2 = domain_adaptive$model$sigma2[2]
  ),
  tibble(
    Outcome_Domain = "ER_Maladaptive",
    estimate = as.numeric(coef(domain_maladaptive$model)[1]),
    se_robust = domain_maladaptive_robust_row$se_robust,
    ci_lb_robust = domain_maladaptive_robust_row$ci_lb_robust,
    ci_ub_robust = domain_maladaptive_robust_row$ci_ub_robust,
    df_robust = domain_maladaptive_robust_row$df_robust,
    p_value_robust = domain_maladaptive_robust_row$p_value_robust,
    k_effect_sizes = domain_maladaptive$model$k,
    n_studies = n_distinct(domain_maladaptive$data$Study_ID),
    n_comparisons = n_distinct(domain_maladaptive$data$Comparison_ID),
    tau2_level_3 = domain_maladaptive$model$sigma2[1],
    tau2_level_2 = domain_maladaptive$model$sigma2[2]
  )
)

cat("\n--- OVERALL SUMMARY TABLE ---\n")
print(overall_summary_table)

cat("\n--- DOMAIN SUMMARY TABLE ---\n")
print(domain_summary_table)


###########################################################
# 16. Save plots
###########################################################

# Funnel plot helper
save_funnel_plot <- function(model_object, file_name, title_text) {
  
  png(
    filename = file.path(output_folder, file_name),
    width = 1800,
    height = 1600,
    res = 220
  )
  
  pooled_effect <- as.numeric(coef(model_object)[1])
  
  funnel(
    model_object,
    refline = pooled_effect,
    lty = 2,
    main = title_text,
    xlab = "Aggregated Hedges' g",
    ylab = "Standard Error"
  )
  
  dev.off()
}

# Overall pooled ER funnel plot
save_funnel_plot(
  egger_overall$model,
  "funnel_plot_study_level_ER_overall.png",
  "Funnel Plot: Study-Level Aggregated ER Effects (All Domains)"
)

# Domain-specific funnel plots
save_funnel_plot(
  egger_total$model,
  "funnel_plot_study_level_ER_total.png",
  "Funnel Plot: Study-Level Aggregated ER_Total Effects"
)

save_funnel_plot(
  egger_adaptive$model,
  "funnel_plot_study_level_ER_adaptive.png",
  "Funnel Plot: Study-Level Aggregated ER_Adaptive Effects"
)

save_funnel_plot(
  egger_maladaptive$model,
  "funnel_plot_study_level_ER_maladaptive.png",
  "Funnel Plot: Study-Level Aggregated ER_Maladaptive Effects"
)

# Leave-one-study-out plot

# Leave-one-study-out plot

loo_plot <- ggplot(
  loo_results,
  aes(x = reorder(omitted_study, estimate), y = estimate)
) +
  geom_point() +
  geom_errorbar(
    aes(ymin = ci_lb_robust, ymax = ci_ub_robust),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    title = "Leave-One-Study-Out Sensitivity Analysis",
    x = "Omitted Study ID",
    y = "Pooled effect estimate"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(output_folder, "leave_one_study_out_plot.png"),
  plot = loo_plot,
  width = 9,
  height = 10,
  dpi = 300
)

###########################################################
# 17. Save model objects and output tables
###########################################################

# Save core tables
write_csv(
  overall_summary_table,
  file.path(output_folder, "overall_ER_model_summary.csv")
)

write_csv(
  domain_summary_table,
  file.path(output_folder, "domain_specific_model_summary.csv")
)

write_csv(
  variance_components,
  file.path(output_folder, "overall_ER_variance_components.csv")
)

write_csv(
  model_fit_summary,
  file.path(output_folder, "overall_ER_model_fit_summary.csv")
)

write_csv(
  study_level_er,
  file.path(output_folder, "study_level_aggregated_ER_effects.csv")
)

write_csv(
  loo_results,
  file.path(output_folder, "leave_one_study_out_results.csv")
)

# Save robust coefficient test tables
write_csv(
  as_tibble(model_er_overall_robust, rownames = "term"),
  file.path(output_folder, "overall_ER_model_robust_tests.csv")
)

write_csv(
  as_tibble(model_er_domain_robust, rownames = "term"),
  file.path(output_folder, "domain_moderator_model_robust_tests.csv")
)

write_csv(
  as_tibble(domain_total$robust, rownames = "term"),
  file.path(output_folder, "ER_Total_robust_tests.csv")
)

write_csv(
  as_tibble(domain_adaptive$robust, rownames = "term"),
  file.path(output_folder, "ER_Adaptive_robust_tests.csv")
)

write_csv(
  as_tibble(domain_maladaptive$robust, rownames = "term"),
  file.path(output_folder, "ER_Maladaptive_robust_tests.csv")
)

write_csv(
  study_level_er_total,
  file.path(output_folder, "study_level_aggregated_ER_total_effects.csv")
)

write_csv(
  study_level_er_adaptive,
  file.path(output_folder, "study_level_aggregated_ER_adaptive_effects.csv")
)

write_csv(
  study_level_er_maladaptive,
  file.path(output_folder, "study_level_aggregated_ER_maladaptive_effects.csv")
)

# Save text summaries for transparency and reproducibility
save_model_output(
  summary(model_er_overall),
  file.path(output_folder, "overall_ER_model_summary.txt"),
  "Overall multilevel random-effects model for emotion regulation"
)

save_model_output(
  model_er_overall_robust,
  file.path(output_folder, "overall_ER_model_robust_tests.txt"),
  "Cluster-robust tests for overall ER model"
)

save_model_output(
  summary(model_er_domain),
  file.path(output_folder, "domain_moderator_model_summary.txt"),
  "Moderator model by emotion regulation domain"
)

save_model_output(
  model_er_domain_robust,
  file.path(output_folder, "domain_moderator_model_robust_tests.txt"),
  "Cluster-robust tests for moderator model by domain"
)

save_model_output(
  domain_omnibus_robust,
  file.path(output_folder, "domain_moderator_omnibus_test.txt"),
  "Robust omnibus test of the emotion regulation domain moderator"
)

save_model_output(
  summary(domain_total$model),
  file.path(output_folder, "ER_Total_model_summary.txt"),
  "Domain-specific model: ER_Total"
)

save_model_output(
  summary(domain_adaptive$model),
  file.path(output_folder, "ER_Adaptive_model_summary.txt"),
  "Domain-specific model: ER_Adaptive"
)

save_model_output(
  summary(domain_maladaptive$model),
  file.path(output_folder, "ER_Maladaptive_model_summary.txt"),
  "Domain-specific model: ER_Maladaptive"
)

save_model_output(
  summary(egger_model),
  file.path(output_folder, "egger_model_summary.txt"),
  "Egger model on aggregated study-level ER effects"
)

save_model_output(
  egger_test,
  file.path(output_folder, "egger_test_results.txt"),
  "Egger test on aggregated study-level ER effects"
)

save_model_output(
  summary(egger_overall$model),
  file.path(output_folder, "egger_model_summary_overall_ER.txt"),
  "Egger model on aggregated study-level ER effects (all domains pooled)"
)

save_model_output(
  egger_overall$test,
  file.path(output_folder, "egger_test_results_overall_ER.txt"),
  "Egger test on aggregated study-level ER effects (all domains pooled)"
)

save_model_output(
  summary(egger_total$model),
  file.path(output_folder, "egger_model_summary_ER_total.txt"),
  "Egger model on aggregated study-level ER_Total effects"
)

save_model_output(
  egger_total$test,
  file.path(output_folder, "egger_test_results_ER_total.txt"),
  "Egger test on aggregated study-level ER_Total effects"
)

save_model_output(
  summary(egger_adaptive$model),
  file.path(output_folder, "egger_model_summary_ER_adaptive.txt"),
  "Egger model on aggregated study-level ER_Adaptive effects"
)

save_model_output(
  egger_adaptive$test,
  file.path(output_folder, "egger_test_results_ER_adaptive.txt"),
  "Egger test on aggregated study-level ER_Adaptive effects"
)

save_model_output(
  summary(egger_maladaptive$model),
  file.path(output_folder, "egger_model_summary_ER_maladaptive.txt"),
  "Egger model on aggregated study-level ER_Maladaptive effects"
)

save_model_output(
  egger_maladaptive$test,
  file.path(output_folder, "egger_test_results_ER_maladaptive.txt"),
  "Egger test on aggregated study-level ER_Maladaptive effects"
)

###########################################################
# 18. Save R objects for later scripts
###########################################################

save(
  dat_er,
  study_level_er,
  study_level_er_total,
  study_level_er_adaptive,
  study_level_er_maladaptive,
  model_er_overall,
  model_er_overall_robust,
  model_er_domain,
  model_er_domain_robust,
  domain_omnibus_robust,
  domain_total,
  domain_adaptive,
  domain_maladaptive,
  egger_overall,
  egger_total,
  egger_adaptive,
  egger_maladaptive,
  loo_results,
  file = file.path(output_folder, "meta_analysis_models_ER.RData")
)

###########################################################
# 19. Final console output
###########################################################

cat("\n====================================================\n")
cat("Meta-analysis model script finished.\n")
cat("Files saved in: data/derived/meta_models/\n")
cat("- overall_ER_model_summary.csv\n")
cat("- domain_specific_model_summary.csv\n")
cat("- overall_ER_variance_components.csv\n")
cat("- overall_ER_model_fit_summary.csv\n")
cat("- study_level_aggregated_ER_effects.csv\n")
cat("- leave_one_study_out_results.csv\n")
cat("- overall_ER_model_robust_tests.csv\n")
cat("- domain_moderator_model_robust_tests.csv\n")
cat("- ER_Total_robust_tests.csv\n")
cat("- ER_Adaptive_robust_tests.csv\n")
cat("- ER_Maladaptive_robust_tests.csv\n")
cat("- overall_ER_model_summary.txt\n")
cat("- overall_ER_model_robust_tests.txt\n")
cat("- domain_moderator_model_summary.txt\n")
cat("- domain_moderator_model_robust_tests.txt\n")
cat("- domain_moderator_omnibus_test.txt\n")
cat("- ER_Total_model_summary.txt\n")
cat("- ER_Adaptive_model_summary.txt\n")
cat("- ER_Maladaptive_model_summary.txt\n")
cat("- egger_model_summary.txt\n")
cat("- egger_test_results.txt\n")
cat("- funnel_plot_study_level_ER.png\n")
cat("- leave_one_study_out_plot.png\n")
cat("- meta_analysis_models_ER.RData\n")
cat("- study_level_aggregated_ER_total_effects.csv\n")
cat("- study_level_aggregated_ER_adaptive_effects.csv\n")
cat("- study_level_aggregated_ER_maladaptive_effects.csv\n")
cat("- egger_model_summary_overall_ER.txt\n")
cat("- egger_test_results_overall_ER.txt\n")
cat("- egger_model_summary_ER_total.txt\n")
cat("- egger_test_results_ER_total.txt\n")
cat("- egger_model_summary_ER_adaptive.txt\n")
cat("- egger_test_results_ER_adaptive.txt\n")
cat("- egger_model_summary_ER_maladaptive.txt\n")
cat("- egger_test_results_ER_maladaptive.txt\n")
cat("- funnel_plot_study_level_ER_overall.png\n")
cat("- funnel_plot_study_level_ER_total.png\n")
cat("- funnel_plot_study_level_ER_adaptive.png\n")
cat("- funnel_plot_study_level_ER_maladaptive.png\n")
cat("====================================================\n")

