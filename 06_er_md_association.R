############################################################
# Script: 06_er_md_association.R
# Project: ER meta-analysis_final
# Purpose:
#   Examine whether study-level effects on emotion regulation
#   (ER) are associated with study-level effects on mental
#   distress (MD) across mindfulness-based intervention (MBI)
#   studies.
#
# Conceptual note:
#   - This script tests an ecological association at the
#     study level.
#   - It does NOT test mediation or within-person mechanism.
#   - It asks whether studies showing larger improvements in
#     ER also tend to show larger improvements in MD.
#
# Main analysis dataset:
#   - Primary effect sizes computed with assumed pre-post
#     correlation r = .50
#   - Study-level aggregated ER and MD effects
#
# Main modelling approach:
#   - Aggregate ER effects within study
#   - Aggregate MD effects within study
#   - Merge studies contributing both ER and MD outcomes
#   - Fit random-effects meta-regression predicting MD effect
#     size from ER effect size
#   - Use inverse-variance weighting based on aggregated MD
#     sampling variance
#   - Provide robust descriptive summaries, scatterplot, and
#     leave-one-study-out sensitivity analysis
#
# Input:
#   data/derived/effect_sizes/er_meta_effect_sizes_r50.csv
#
# Output:
#   data/derived/er_md_association/
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
library(ggplot2)


############################
# 2. Define file locations
############################

input_file <- "data/derived/effect_sizes/er_meta_effect_sizes_r50.csv"
output_folder <- "data/derived/er_md_association"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


############################
# 3. Read effect size dataset
############################

dat <- read_csv(input_file, show_col_types = FALSE)

cat("\n----------------------------------------\n")
cat("ER-MD association script started\n")
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

dat <- dat %>%
  mutate(
    Study_ID = as.factor(Study_ID),
    Comparison_ID = as.factor(Comparison_ID)
  )

if (sum(is.na(dat$yi)) > 0) {
  stop("Dataset contains missing yi values.")
}

if (sum(is.na(dat$vi)) > 0) {
  stop("Dataset contains missing vi values.")
}

if (sum(dat$vi <= 0, na.rm = TRUE) > 0) {
  stop("Dataset contains non-positive vi values.")
}


###########################################################
# 6. Helper function to save text output
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
# 7. Aggregate ER and MD effects to the study level
###########################################################

# Rationale:
# Studies contribute multiple dependent effect sizes.
# For this cross-domain association analysis, the unit of
# analysis is the study.
#
# Aggregation approach:
# - Within each study, aggregate effect sizes using
#   inverse-variance weighted means.
# - This is used here for descriptive ecological association
#   analysis, not as a replacement for the multilevel models
#   in Scripts 03-05.

aggregate_study_effects <- function(data, analysis_set_value, value_name) {
  
  data %>%
    filter(analysis_set == analysis_set_value) %>%
    group_by(Study_ID) %>%
    summarise(
      k_effects = n(),
      yi_agg = weighted.mean(yi, w = 1 / vi, na.rm = TRUE),
      vi_agg = 1 / sum(1 / vi, na.rm = TRUE),
      sei_agg = sqrt(vi_agg),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lb = yi_agg - 1.96 * sei_agg,
      ci_ub = yi_agg + 1.96 * sei_agg
    ) %>%
    rename(
      !!paste0(value_name, "_k_effects") := k_effects,
      !!paste0(value_name, "_yi") := yi_agg,
      !!paste0(value_name, "_vi") := vi_agg,
      !!paste0(value_name, "_sei") := sei_agg,
      !!paste0(value_name, "_ci_lb") := ci_lb,
      !!paste0(value_name, "_ci_ub") := ci_ub
    ) %>%
    arrange(as.numeric(as.character(Study_ID)))
}

study_level_er <- aggregate_study_effects(
  data = dat,
  analysis_set_value = "ER",
  value_name = "ER"
)

study_level_md <- aggregate_study_effects(
  data = dat,
  analysis_set_value = "MD",
  value_name = "MD"
)

cat("\n--- STUDY-LEVEL AGGREGATED DATASETS ---\n")
cat("Studies with ER data:", nrow(study_level_er), "\n")
cat("Studies with MD data:", nrow(study_level_md), "\n")


###########################################################
# 8. Merge studies contributing both ER and MD outcomes
###########################################################

association_dat <- study_level_er %>%
  inner_join(study_level_md, by = "Study_ID") %>%
  arrange(as.numeric(as.character(Study_ID)))

cat("\n--- MERGED ER-MD STUDY-LEVEL DATASET ---\n")
cat("Studies with both ER and MD data:", nrow(association_dat), "\n")

if (nrow(association_dat) < 10) {
  stop("Fewer than 10 studies contribute both ER and MD outcomes. Association analysis not stable enough to proceed.")
}

print(association_dat, n = nrow(association_dat))


###########################################################
# 9. Descriptive correlations
###########################################################

# These correlations are descriptive and unweighted.
# The meta-regression below is the main inferential analysis.

pearson_cor <- cor.test(
  association_dat$ER_yi,
  association_dat$MD_yi,
  method = "pearson"
)

spearman_cor <- cor.test(
  association_dat$ER_yi,
  association_dat$MD_yi,
  method = "spearman",
  exact = FALSE
)

correlation_summary <- tibble(
  statistic = c("Pearson r", "Spearman rho"),
  estimate = c(
    unname(pearson_cor$estimate),
    unname(spearman_cor$estimate)
  ),
  p_value = c(
    pearson_cor$p.value,
    spearman_cor$p.value
  ),
  ci_lb = c(
    pearson_cor$conf.int[1],
    NA_real_
  ),
  ci_ub = c(
    pearson_cor$conf.int[2],
    NA_real_
  ),
  n_studies = nrow(association_dat)
)

cat("\n--- DESCRIPTIVE CORRELATIONS ---\n")
print(correlation_summary)


###########################################################
# 10. Main meta-regression model
###########################################################

# Main inferential model:
# Predict aggregated MD effect size from aggregated ER effect
# size across studies.
#
# Interpretation:
# - Positive slope: studies with larger ER improvements also
#   tend to show larger MD improvements.
# - This is a study-level association, not a causal mediation
#   test.

model_er_md <- rma(
  yi = MD_yi,
  vi = MD_vi,
  mods = ~ ER_yi,
  data = association_dat,
  method = "REML",
  test = "t"
)

cat("\n--- MAIN ER-MD META-REGRESSION MODEL ---\n")
print(summary(model_er_md))

er_md_model_summary <- tibble(
  model = "Random-effects meta-regression: MD ~ ER",
  intercept = as.numeric(coef(model_er_md)[1]),
  slope_ER = as.numeric(coef(model_er_md)[2]),
  intercept_se = model_er_md$se[1],
  slope_se = model_er_md$se[2],
  intercept_ci_lb = model_er_md$ci.lb[1],
  intercept_ci_ub = model_er_md$ci.ub[1],
  slope_ci_lb = model_er_md$ci.lb[2],
  slope_ci_ub = model_er_md$ci.ub[2],
  intercept_p = model_er_md$pval[1],
  slope_p = model_er_md$pval[2],
  tau2 = model_er_md$tau2,
  QE = model_er_md$QE,
  QE_df = model_er_md$k - model_er_md$p,
  QE_p = model_er_md$QEp,
  QM = model_er_md$QM,
  QM_df = model_er_md$m,
  QM_p = model_er_md$QMp,
  k_studies = model_er_md$k
)

cat("\n--- ER-MD META-REGRESSION SUMMARY TABLE ---\n")
print(er_md_model_summary)


###########################################################
# 11. Influence and leave-one-study-out analysis
###########################################################

loo_results <- vector("list", length = nrow(association_dat))
study_ids <- association_dat$Study_ID

for (i in seq_along(study_ids)) {
  
  current_study <- study_ids[i]
  
  dat_minus_one <- association_dat %>%
    filter(Study_ID != current_study)
  
  model_minus_one <- rma(
    yi = MD_yi,
    vi = MD_vi,
    mods = ~ ER_yi,
    data = dat_minus_one,
    method = "REML",
    test = "t"
  )
  
  loo_results[[i]] <- tibble(
    omitted_study = as.character(current_study),
    intercept = as.numeric(coef(model_minus_one)[1]),
    slope_ER = as.numeric(coef(model_minus_one)[2]),
    slope_se = model_minus_one$se[2],
    slope_ci_lb = model_minus_one$ci.lb[2],
    slope_ci_ub = model_minus_one$ci.ub[2],
    slope_p = model_minus_one$pval[2],
    tau2 = model_minus_one$tau2,
    k_remaining = model_minus_one$k
  )
}

loo_results <- bind_rows(loo_results) %>%
  arrange(desc(abs(slope_ER)))

cat("\n--- LEAVE-ONE-STUDY-OUT RESULTS ---\n")
print(loo_results, n = nrow(loo_results))


###########################################################
# 12. Optional sensitivity analysis:
#     exclude studies with only 1 ER or 1 MD effect
###########################################################

# This checks whether the association is robust when retaining
# only studies with more than one contributing effect in both
# domains before aggregation.

association_dat_multi <- association_dat %>%
  filter(ER_k_effects > 1, MD_k_effects > 1)

if (nrow(association_dat_multi) >= 10) {
  
  model_er_md_multi <- rma(
    yi = MD_yi,
    vi = MD_vi,
    mods = ~ ER_yi,
    data = association_dat_multi,
    method = "REML",
    test = "t"
  )
  
  sensitivity_summary <- tibble(
    analysis = c(
      "Primary ER-MD association model",
      "Sensitivity: studies with >1 ER and >1 MD effect"
    ),
    slope_ER = c(
      as.numeric(coef(model_er_md)[2]),
      as.numeric(coef(model_er_md_multi)[2])
    ),
    slope_se = c(
      model_er_md$se[2],
      model_er_md_multi$se[2]
    ),
    slope_ci_lb = c(
      model_er_md$ci.lb[2],
      model_er_md_multi$ci.lb[2]
    ),
    slope_ci_ub = c(
      model_er_md$ci.ub[2],
      model_er_md_multi$ci.ub[2]
    ),
    slope_p = c(
      model_er_md$pval[2],
      model_er_md_multi$pval[2]
    ),
    tau2 = c(
      model_er_md$tau2,
      model_er_md_multi$tau2
    ),
    k_studies = c(
      model_er_md$k,
      model_er_md_multi$k
    )
  )
  
  cat("\n--- SENSITIVITY ANALYSIS SUMMARY ---\n")
  print(sensitivity_summary)
  
} else {
  
  model_er_md_multi <- NULL
  
  sensitivity_summary <- tibble(
    analysis = "Sensitivity: studies with >1 ER and >1 MD effect",
    slope_ER = NA_real_,
    slope_se = NA_real_,
    slope_ci_lb = NA_real_,
    slope_ci_ub = NA_real_,
    slope_p = NA_real_,
    tau2 = NA_real_,
    k_studies = nrow(association_dat_multi),
    note = "Not estimated: fewer than 10 studies remained."
  )
  
  cat("\n--- SENSITIVITY ANALYSIS SUMMARY ---\n")
  print(sensitivity_summary)
}


###########################################################
# 13. Create plots
###########################################################

scatter_plot <- ggplot(
  association_dat,
  aes(x = ER_yi, y = MD_yi)
) +
  geom_point(aes(size = 1 / MD_vi), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Association between study-level ER and MD effects",
    x = "Aggregated ER effect (Hedges' g)",
    y = "Aggregated MD effect (Hedges' g)",
    size = "Precision\n(1 / MD variance)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(output_folder, "scatterplot_ER_MD_association.png"),
  plot = scatter_plot,
  width = 8,
  height = 6,
  dpi = 300
)

loo_plot <- ggplot(
  loo_results,
  aes(x = reorder(omitted_study, slope_ER), y = slope_ER)
) +
  geom_point() +
  geom_errorbar(
    aes(ymin = slope_ci_lb, ymax = slope_ci_ub),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    title = "Leave-one-study-out analysis for ER-MD slope",
    x = "Omitted Study ID",
    y = "Slope of ER predicting MD"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(output_folder, "leave_one_study_out_ER_MD_slope.png"),
  plot = loo_plot,
  width = 9,
  height = 10,
  dpi = 300
)


###########################################################
# 14. Save output tables
###########################################################

write_csv(
  study_level_er,
  file.path(output_folder, "study_level_aggregated_ER_effects.csv")
)

write_csv(
  study_level_md,
  file.path(output_folder, "study_level_aggregated_MD_effects.csv")
)

write_csv(
  association_dat,
  file.path(output_folder, "study_level_ER_MD_overlap_dataset.csv")
)

write_csv(
  correlation_summary,
  file.path(output_folder, "descriptive_correlations_ER_MD.csv")
)

write_csv(
  er_md_model_summary,
  file.path(output_folder, "er_md_meta_regression_summary.csv")
)

write_csv(
  loo_results,
  file.path(output_folder, "leave_one_study_out_ER_MD_results.csv")
)

write_csv(
  sensitivity_summary,
  file.path(output_folder, "sensitivity_ER_MD_results.csv")
)


###########################################################
# 15. Save text summaries
###########################################################

save_model_output(
  summary(model_er_md),
  file.path(output_folder, "er_md_meta_regression_summary.txt"),
  "Random-effects meta-regression of study-level MD effects on study-level ER effects"
)

save_model_output(
  pearson_cor,
  file.path(output_folder, "pearson_correlation_ER_MD.txt"),
  "Pearson correlation between aggregated ER and MD effects"
)

save_model_output(
  spearman_cor,
  file.path(output_folder, "spearman_correlation_ER_MD.txt"),
  "Spearman correlation between aggregated ER and MD effects"
)

if (!is.null(model_er_md_multi)) {
  save_model_output(
    summary(model_er_md_multi),
    file.path(output_folder, "er_md_meta_regression_sensitivity_summary.txt"),
    "Sensitivity meta-regression for studies with >1 ER and >1 MD effect"
  )
}


###########################################################
# 16. Save R objects
###########################################################

save(
  dat,
  study_level_er,
  study_level_md,
  association_dat,
  pearson_cor,
  spearman_cor,
  model_er_md,
  model_er_md_multi,
  loo_results,
  sensitivity_summary,
  file = file.path(output_folder, "er_md_association_analysis.RData")
)


###########################################################
# 17. Final console output
###########################################################

cat("\n====================================================\n")
cat("ER-MD association script finished.\n")
cat("Files saved in: data/derived/er_md_association/\n")
cat("- study_level_aggregated_ER_effects.csv\n")
cat("- study_level_aggregated_MD_effects.csv\n")
cat("- study_level_ER_MD_overlap_dataset.csv\n")
cat("- descriptive_correlations_ER_MD.csv\n")
cat("- er_md_meta_regression_summary.csv\n")
cat("- leave_one_study_out_ER_MD_results.csv\n")
cat("- sensitivity_ER_MD_results.csv\n")
cat("- er_md_meta_regression_summary.txt\n")
cat("- pearson_correlation_ER_MD.txt\n")
cat("- spearman_correlation_ER_MD.txt\n")
cat("- scatterplot_ER_MD_association.png\n")
cat("- leave_one_study_out_ER_MD_slope.png\n")
cat("- er_md_association_analysis.RData\n")
cat("====================================================\n")

###########################################################
# Additional sensitivity: exclude pre-identified extreme ER studies
###########################################################

# Rationale:
# Study 43 clearly contributes an extreme ER effect.
# Study 30 also contributes a very large ER effect.
# These analyses test whether the ER-MD association remains
# when those overlap studies are excluded.

extreme_er_studies <- c("43", "30")

association_dat_no_extreme <- association_dat %>%
  filter(!as.character(Study_ID) %in% extreme_er_studies)

cat("\n--- SENSITIVITY DATASET: EXCLUDING PRE-IDENTIFIED EXTREME ER STUDIES ---\n")
cat("Studies remaining:", nrow(association_dat_no_extreme), "\n")
cat("Excluded Study_IDs:", paste(extreme_er_studies, collapse = ", "), "\n")

if (nrow(association_dat_no_extreme) >= 10) {
  
  pearson_cor_no_extreme <- cor.test(
    association_dat_no_extreme$ER_yi,
    association_dat_no_extreme$MD_yi,
    method = "pearson"
  )
  
  spearman_cor_no_extreme <- cor.test(
    association_dat_no_extreme$ER_yi,
    association_dat_no_extreme$MD_yi,
    method = "spearman",
    exact = FALSE
  )
  
  model_er_md_no_extreme <- rma(
    yi = MD_yi,
    vi = MD_vi,
    mods = ~ ER_yi,
    data = association_dat_no_extreme,
    method = "REML",
    test = "t"
  )
  
  sensitivity_extreme_summary <- tibble(
    analysis = c(
      "Primary ER-MD association model",
      "Sensitivity: exclude pre-identified extreme ER studies (43, 30)"
    ),
    slope_ER = c(
      as.numeric(coef(model_er_md)[2]),
      as.numeric(coef(model_er_md_no_extreme)[2])
    ),
    slope_se = c(
      model_er_md$se[2],
      model_er_md_no_extreme$se[2]
    ),
    slope_ci_lb = c(
      model_er_md$ci.lb[2],
      model_er_md_no_extreme$ci.lb[2]
    ),
    slope_ci_ub = c(
      model_er_md$ci.ub[2],
      model_er_md_no_extreme$ci.ub[2]
    ),
    slope_p = c(
      model_er_md$pval[2],
      model_er_md_no_extreme$pval[2]
    ),
    tau2 = c(
      model_er_md$tau2,
      model_er_md_no_extreme$tau2
    ),
    k_studies = c(
      model_er_md$k,
      model_er_md_no_extreme$k
    )
  )
  
  sensitivity_extreme_correlations <- tibble(
    analysis = c(
      "Primary overlap dataset",
      "Exclude pre-identified extreme ER studies (43, 30)"
    ),
    pearson_r = c(
      unname(pearson_cor$estimate),
      unname(pearson_cor_no_extreme$estimate)
    ),
    pearson_p = c(
      pearson_cor$p.value,
      pearson_cor_no_extreme$p.value
    ),
    spearman_rho = c(
      unname(spearman_cor$estimate),
      unname(spearman_cor_no_extreme$estimate)
    ),
    spearman_p = c(
      spearman_cor$p.value,
      spearman_cor_no_extreme$p.value
    ),
    n_studies = c(
      nrow(association_dat),
      nrow(association_dat_no_extreme)
    )
  )
  
  cat("\n--- SENSITIVITY SUMMARY: EXCLUDING PRE-IDENTIFIED EXTREME ER STUDIES ---\n")
  print(sensitivity_extreme_summary)
  
  cat("\n--- CORRELATION SUMMARY: EXCLUDING PRE-IDENTIFIED EXTREME ER STUDIES ---\n")
  print(sensitivity_extreme_correlations)
  
} else {
  
  pearson_cor_no_extreme <- NULL
  spearman_cor_no_extreme <- NULL
  model_er_md_no_extreme <- NULL
  
  sensitivity_extreme_summary <- tibble(
    analysis = "Sensitivity: exclude pre-identified extreme ER studies (43, 30)",
    slope_ER = NA_real_,
    slope_se = NA_real_,
    slope_ci_lb = NA_real_,
    slope_ci_ub = NA_real_,
    slope_p = NA_real_,
    tau2 = NA_real_,
    k_studies = nrow(association_dat_no_extreme),
    note = "Not estimated: fewer than 10 studies remained."
  )
  
  sensitivity_extreme_correlations <- tibble(
    analysis = "Exclude pre-identified extreme ER studies (43, 30)",
    pearson_r = NA_real_,
    pearson_p = NA_real_,
    spearman_rho = NA_real_,
    spearman_p = NA_real_,
    n_studies = nrow(association_dat_no_extreme),
    note = "Not estimated: fewer than 10 studies remained."
  )
  
  cat("\n--- SENSITIVITY SUMMARY: EXCLUDING PRE-IDENTIFIED EXTREME ER STUDIES ---\n")
  print(sensitivity_extreme_summary)
}
write_csv(
  association_dat_no_extreme,
  file.path(output_folder, "study_level_ER_MD_overlap_dataset_no_extreme_ER_studies.csv")
)

write_csv(
  sensitivity_extreme_summary,
  file.path(output_folder, "sensitivity_ER_MD_excluding_extreme_ER_studies.csv")
)

write_csv(
  sensitivity_extreme_correlations,
  file.path(output_folder, "descriptive_correlations_ER_MD_excluding_extreme_ER_studies.csv")
)

if (!is.null(model_er_md_no_extreme)) {
  save_model_output(
    summary(model_er_md_no_extreme),
    file.path(output_folder, "er_md_meta_regression_excluding_extreme_ER_studies.txt"),
    "Sensitivity meta-regression excluding pre-identified extreme ER studies (43, 30)"
  )
}

if (!is.null(pearson_cor_no_extreme)) {
  save_model_output(
    pearson_cor_no_extreme,
    file.path(output_folder, "pearson_correlation_ER_MD_excluding_extreme_ER_studies.txt"),
    "Pearson correlation excluding pre-identified extreme ER studies (43, 30)"
  )
}

if (!is.null(spearman_cor_no_extreme)) {
  save_model_output(
    spearman_cor_no_extreme,
    file.path(output_folder, "spearman_correlation_ER_MD_excluding_extreme_ER_studies.txt"),
    "Spearman correlation excluding pre-identified extreme ER studies (43, 30)"
  )
}

cat("- study_level_ER_MD_overlap_dataset_no_extreme_ER_studies.csv\n")
cat("- sensitivity_ER_MD_excluding_extreme_ER_studies.csv\n")
cat("- descriptive_correlations_ER_MD_excluding_extreme_ER_studies.csv\n")
cat("- er_md_meta_regression_excluding_extreme_ER_studies.txt\n")
cat("- pearson_correlation_ER_MD_excluding_extreme_ER_studies.txt\n")
cat("- spearman_correlation_ER_MD_excluding_extreme_ER_studies.txt\n")

