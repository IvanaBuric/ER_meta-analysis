############################################################
# Script: 07_risk_of_bias_and_study_quality_ER.R
# Project: ER meta-analysis_final
# Purpose:
#   Examine risk of bias (RoB 2) domains as study-level
#   moderators of the pooled emotion regulation (ER) effect
#   of mindfulness-based interventions (MBIs), and run
#   sensitivity analyses excluding high-risk studies.
#
# Main analysis dataset:
#   - Primary effect sizes computed with assumed pre-post
#     correlation r = .50
#   - Emotion regulation outcomes only
#
# Main modelling approach:
#   - Multilevel random-effects meta-analysis
#   - Effect sizes nested within comparisons nested within
#     studies
#   - Models estimated with REML
#   - Cluster-robust inference at the study level (CR2)
#   - Robust confidence intervals based on Satterthwaite-
#     adjusted degrees of freedom
#   - Each ROB2 domain tested in a separate model
#   - Outcome_Domain retained as an adjustment covariate
#
# Risk of bias strategy:
#   - ROB2_D1 to ROB2_D5 are treated as categorical,
#     study-level moderators
#   - Only domains with sufficient usable variability are
#     estimated
#   - Additional sensitivity analyses exclude studies rated
#     as high risk (ROB2 = 2) on each domain separately
#
# ROB2 coding:
#   0 = Low risk
#   1 = Some concerns
#   2 = High risk
#
# Input:
#   data/derived/effect_sizes/er_meta_effect_sizes_r50.csv
#
# Output:
#   data/derived/risk_of_bias_ER/
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
  "purrr",
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
library(purrr)
library(metafor)
library(clubSandwich)
library(ggplot2)


############################
# 2. Define file locations
############################

input_file <- "data/derived/effect_sizes/er_meta_effect_sizes_r50.csv"
output_folder <- "data/derived/risk_of_bias_ER"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


############################
# 3. Read effect size dataset
############################

dat <- read_csv(input_file, show_col_types = FALSE)

cat("\n----------------------------------------\n")
cat("Risk of bias analysis script started\n")
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
  "analysis_set",
  "ROB2_D1",
  "ROB2_D2",
  "ROB2_D3",
  "ROB2_D4",
  "ROB2_D5"
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
# 5. Prepare ER analysis dataset
###########################################################

# Assumption:
# Missing / unclear risk-of-bias codes have already been
# converted to NA in earlier data-cleaning steps.

dat_er <- dat %>%
  filter(analysis_set == "ER") %>%
  mutate(
    Study_ID = as.factor(Study_ID),
    Comparison_ID = as.factor(Comparison_ID),
    Outcome_Domain = factor(
      Outcome_Domain,
      levels = c("ER_Total", "ER_Adaptive", "ER_Maladaptive")
    ),
    ROB2_D1 = factor(
      ROB2_D1,
      levels = c(0, 1, 2),
      labels = c("Low risk", "Some concerns", "High risk")
    ),
    ROB2_D2 = factor(
      ROB2_D2,
      levels = c(0, 1, 2),
      labels = c("Low risk", "Some concerns", "High risk")
    ),
    ROB2_D3 = factor(
      ROB2_D3,
      levels = c(0, 1, 2),
      labels = c("Low risk", "Some concerns", "High risk")
    ),
    ROB2_D4 = factor(
      ROB2_D4,
      levels = c(0, 1, 2),
      labels = c("Low risk", "Some concerns", "High risk")
    ),
    ROB2_D5 = factor(
      ROB2_D5,
      levels = c(0, 1, 2),
      labels = c("Low risk", "Some concerns", "High risk")
    )
  )

cat("\n--- MAIN RISK OF BIAS DATASET ---\n")
cat("ER rows:", nrow(dat_er), "\n")
cat("Unique studies:", n_distinct(dat_er$Study_ID), "\n")
cat("Unique comparisons:", n_distinct(dat_er$Comparison_ID), "\n")

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
# 6. Helper functions
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

get_robust_tests <- function(model, data) {
  coef_test(
    model,
    vcov = "CR2",
    cluster = data$Study_ID
  )
}

extract_robust_rows <- function(robust_object) {
  robust_tbl <- as_tibble(robust_object, rownames = "term")
  t_crit <- qt(0.975, df = robust_tbl$df_Satt)
  
  robust_tbl %>%
    mutate(
      ci_lb_robust = beta - t_crit * SE,
      ci_ub_robust = beta + t_crit * SE
    ) %>%
    rename(
      estimate_robust = beta,
      se_robust = SE,
      t_stat_robust = tstat,
      df_robust = df_Satt,
      p_value_robust = p_Satt
    )
}

extract_intercept_robust_summary <- function(model, data, analysis_label, moderator, label) {
  robust <- get_robust_tests(model, data)
  robust_row <- as.data.frame(robust)[1, ]
  t_crit <- qt(0.975, df = robust_row$df_Satt)
  
  tibble(
    moderator = moderator,
    label = label,
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


###########################################################
# 7. Define ROB2 specifications
###########################################################

rob_specs <- tribble(
  ~moderator, ~label,
  "ROB2_D1", "ROB2 Domain 1: Randomisation process",
  "ROB2_D2", "ROB2 Domain 2: Deviations from intended interventions",
  "ROB2_D3", "ROB2 Domain 3: Missing outcome data",
  "ROB2_D4", "ROB2 Domain 4: Measurement of the outcome",
  "ROB2_D5", "ROB2 Domain 5: Selection of the reported result"
)

# Pragmatic thresholds for estimation
min_total_studies_for_model <- 10
min_studies_per_factor_level <- 4


###########################################################
# 8. Prepare moderator-specific ROB2 datasets
###########################################################

prepare_rob_dataset <- function(data, moderator) {
  
  dat_mod <- data %>%
    filter(!is.na(.data[[moderator]])) %>%
    mutate(
      moderator_value = droplevels(as.factor(.data[[moderator]]))
    ) %>%
    filter(!is.na(moderator_value))
  
  level_counts <- dat_mod %>%
    group_by(moderator_value) %>%
    summarise(
      k_effects = n(),
      n_studies = n_distinct(Study_ID),
      .groups = "drop"
    ) %>%
    arrange(moderator_value)
  
  valid_levels <- level_counts %>%
    filter(n_studies >= min_studies_per_factor_level) %>%
    pull(moderator_value) %>%
    as.character()
  
  dat_mod <- dat_mod %>%
    filter(as.character(moderator_value) %in% valid_levels) %>%
    mutate(
      moderator_value = droplevels(moderator_value)
    )
  
  level_counts_final <- dat_mod %>%
    group_by(moderator_value) %>%
    summarise(
      k_effects = n(),
      n_studies = n_distinct(Study_ID),
      .groups = "drop"
    ) %>%
    arrange(moderator_value)
  
  list(
    data = dat_mod,
    level_counts = level_counts_final
  )
}


###########################################################
# 9. Fit single ROB2 moderator models
###########################################################

fit_single_rob_moderator <- function(data, moderator, label) {
  
  prepared <- prepare_rob_dataset(data, moderator)
  dat_mod <- prepared$data
  level_counts <- prepared$level_counts
  
  n_studies <- n_distinct(dat_mod$Study_ID)
  n_effects <- nrow(dat_mod)
  n_comparisons <- n_distinct(dat_mod$Comparison_ID)
  
  if (n_effects == 0) {
    return(list(
      moderator = moderator,
      label = label,
      status = "Not estimated: no usable rows after filtering missing values / sparse levels.",
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      omnibus = NULL
    ))
  }
  
  if (n_studies < min_total_studies_for_model) {
    return(list(
      moderator = moderator,
      label = label,
      status = paste0(
        "Not estimated: fewer than ", min_total_studies_for_model,
        " studies with usable data."
      ),
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      omnibus = NULL
    ))
  }
  
  if (n_distinct(dat_mod$moderator_value) < 2) {
    return(list(
      moderator = moderator,
      label = label,
      status = "Not estimated: fewer than 2 factor levels after filtering sparse levels.",
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      omnibus = NULL
    ))
  }
  
  model <- rma.mv(
    yi = yi,
    V = vi,
    mods = ~ Outcome_Domain + moderator_value,
    random = ~ 1 | Study_ID/Comparison_ID,
    data = dat_mod,
    method = "REML",
    test = "t"
  )
  
  robust <- get_robust_tests(model, dat_mod)
  
  mod_coef_idx <- grep("^moderator_value", names(coef(model)))
  
  omnibus <- Wald_test(
    model,
    constraints = constrain_zero(mod_coef_idx),
    vcov = "CR2",
    cluster = dat_mod$Study_ID
  )
  
  list(
    moderator = moderator,
    label = label,
    status = "Estimated successfully.",
    data = dat_mod,
    level_counts = level_counts,
    model = model,
    robust = robust,
    omnibus = omnibus
  )
}


###########################################################
# 10. Run all ROB2 moderator analyses
###########################################################

rob_results <- pmap(
  rob_specs,
  function(moderator, label) {
    fit_single_rob_moderator(
      data = dat_er,
      moderator = moderator,
      label = label
    )
  }
)

names(rob_results) <- rob_specs$moderator


###########################################################
# 11. Create ROB2 overview table
###########################################################

rob_overview <- map_dfr(rob_results, function(x) {
  
  if (is.null(x$model)) {
    return(tibble(
      moderator = x$moderator,
      label = x$label,
      status = x$status,
      k_effect_sizes = nrow(x$data),
      n_studies = n_distinct(x$data$Study_ID),
      n_comparisons = n_distinct(x$data$Comparison_ID),
      omnibus_stat = NA_real_,
      omnibus_df_num = NA_real_,
      omnibus_df_den = NA_real_,
      omnibus_p = NA_real_
    ))
  }
  
  omnibus_object <- x$omnibus
  
  omnibus_p <- tryCatch(
    if ("p_val" %in% names(omnibus_object)) omnibus_object$p_val else omnibus_object$QMp,
    error = function(e) NA_real_
  )
  
  omnibus_stat <- tryCatch(
    if ("Fstat" %in% names(omnibus_object)) omnibus_object$Fstat else omnibus_object$QM,
    error = function(e) NA_real_
  )
  
  omnibus_df_num <- tryCatch(
    if ("df_num" %in% names(omnibus_object)) omnibus_object$df_num else omnibus_object$m,
    error = function(e) NA_real_
  )
  
  omnibus_df_den <- tryCatch(
    if ("df_denom" %in% names(omnibus_object)) omnibus_object$df_denom else NA_real_,
    error = function(e) NA_real_
  )
  
  tibble(
    moderator = x$moderator,
    label = x$label,
    status = x$status,
    k_effect_sizes = x$model$k,
    n_studies = n_distinct(x$data$Study_ID),
    n_comparisons = n_distinct(x$data$Comparison_ID),
    omnibus_stat = omnibus_stat,
    omnibus_df_num = omnibus_df_num,
    omnibus_df_den = omnibus_df_den,
    omnibus_p = omnibus_p
  )
})

cat("\n--- RISK OF BIAS MODERATOR OVERVIEW ---\n")
print(rob_overview, n = nrow(rob_overview), width = Inf)


###########################################################
# 12. Create model-based coefficient table (transparency only)
###########################################################

rob_coefficients <- map_dfr(rob_results, function(x) {
  
  if (is.null(x$model)) {
    return(tibble(
      moderator = x$moderator,
      label = x$label,
      term = NA_character_,
      estimate = NA_real_,
      se_model = NA_real_,
      ci_lb_model = NA_real_,
      ci_ub_model = NA_real_,
      p_value_model = NA_real_
    ))
  }
  
  tibble(
    moderator = x$moderator,
    label = x$label,
    term = names(coef(x$model)),
    estimate = as.numeric(coef(x$model)),
    se_model = x$model$se,
    ci_lb_model = x$model$ci.lb,
    ci_ub_model = x$model$ci.ub,
    p_value_model = x$model$pval
  )
})

cat("\n--- RISK OF BIAS MODERATOR COEFFICIENT TABLE ---\n")
print(rob_coefficients, n = nrow(rob_coefficients), width = Inf)


###########################################################
# 13. Create robust coefficient table (manuscript-facing)
###########################################################

rob_robust_tests <- map_dfr(rob_results, function(x) {
  
  if (is.null(x$robust)) {
    return(tibble(
      moderator = x$moderator,
      label = x$label,
      term = NA_character_,
      estimate_robust = NA_real_,
      se_robust = NA_real_,
      ci_lb_robust = NA_real_,
      ci_ub_robust = NA_real_,
      t_stat_robust = NA_real_,
      df_robust = NA_real_,
      p_value_robust = NA_real_
    ))
  }
  
  extract_robust_rows(x$robust) %>%
    mutate(
      moderator = x$moderator,
      label = x$label,
      .before = 1
    ) %>%
    select(
      moderator,
      label,
      term,
      estimate_robust,
      se_robust,
      ci_lb_robust,
      ci_ub_robust,
      t_stat_robust,
      df_robust,
      p_value_robust
    )
})

cat("\n--- RISK OF BIAS ROBUST TEST TABLE ---\n")
print(rob_robust_tests, n = nrow(rob_robust_tests), width = Inf)


###########################################################
# 14. Sensitivity analyses excluding high-risk studies
###########################################################

fit_overall_er_model <- function(data) {
  rma.mv(
    yi = yi,
    V = vi,
    random = ~ 1 | Study_ID/Comparison_ID,
    data = data,
    method = "REML",
    test = "t"
  )
}

rob_sensitivity_results <- vector("list", length = nrow(rob_specs))

for (i in seq_len(nrow(rob_specs))) {
  
  current_moderator <- rob_specs$moderator[i]
  current_label <- rob_specs$label[i]
  
  dat_current <- dat_er %>%
    filter(!is.na(.data[[current_moderator]]))
  
  n_studies_total <- n_distinct(dat_current$Study_ID)
  
  n_studies_high <- dat_current %>%
    filter(.data[[current_moderator]] == "High risk") %>%
    summarise(n = n_distinct(Study_ID)) %>%
    pull(n)
  
  dat_no_high <- dat_current %>%
    filter(.data[[current_moderator]] != "High risk")
  
  n_studies_no_high <- n_distinct(dat_no_high$Study_ID)
  
  if (nrow(dat_current) == 0) {
    rob_sensitivity_results[[i]] <- tibble(
      moderator = current_moderator,
      label = current_label,
      analysis = "Exclude high-risk studies",
      status = "Not estimated: no usable rows for this ROB2 domain.",
      estimate = NA_real_,
      se_robust = NA_real_,
      t_robust = NA_real_,
      df_robust = NA_real_,
      p_value_robust = NA_real_,
      ci_lb_robust = NA_real_,
      ci_ub_robust = NA_real_,
      k_effect_sizes = NA_integer_,
      n_studies = n_studies_no_high,
      n_comparisons = NA_integer_,
      tau2_level_3 = NA_real_,
      tau2_level_2 = NA_real_,
      QE = NA_real_,
      QE_df = NA_real_,
      QE_p = NA_real_,
      n_studies_total = n_studies_total,
      n_studies_high = n_studies_high
    )
    next
  }
  
  if (n_studies_high == 0) {
    rob_sensitivity_results[[i]] <- tibble(
      moderator = current_moderator,
      label = current_label,
      analysis = "Exclude high-risk studies",
      status = "Not estimated: no studies rated High risk on this domain.",
      estimate = NA_real_,
      se_robust = NA_real_,
      t_robust = NA_real_,
      df_robust = NA_real_,
      p_value_robust = NA_real_,
      ci_lb_robust = NA_real_,
      ci_ub_robust = NA_real_,
      k_effect_sizes = NA_integer_,
      n_studies = n_studies_no_high,
      n_comparisons = NA_integer_,
      tau2_level_3 = NA_real_,
      tau2_level_2 = NA_real_,
      QE = NA_real_,
      QE_df = NA_real_,
      QE_p = NA_real_,
      n_studies_total = n_studies_total,
      n_studies_high = n_studies_high
    )
    next
  }
  
  if (n_studies_no_high < min_total_studies_for_model) {
    rob_sensitivity_results[[i]] <- tibble(
      moderator = current_moderator,
      label = current_label,
      analysis = "Exclude high-risk studies",
      status = paste0(
        "Not estimated: fewer than ", min_total_studies_for_model,
        " studies remain after excluding High risk."
      ),
      estimate = NA_real_,
      se_robust = NA_real_,
      t_robust = NA_real_,
      df_robust = NA_real_,
      p_value_robust = NA_real_,
      ci_lb_robust = NA_real_,
      ci_ub_robust = NA_real_,
      k_effect_sizes = NA_integer_,
      n_studies = n_studies_no_high,
      n_comparisons = NA_integer_,
      tau2_level_3 = NA_real_,
      tau2_level_2 = NA_real_,
      QE = NA_real_,
      QE_df = NA_real_,
      QE_p = NA_real_,
      n_studies_total = n_studies_total,
      n_studies_high = n_studies_high
    )
    next
  }
  
  model_no_high <- fit_overall_er_model(dat_no_high)
  
  rob_sensitivity_results[[i]] <- extract_intercept_robust_summary(
    model = model_no_high,
    data = dat_no_high,
    analysis_label = "Exclude high-risk studies",
    moderator = current_moderator,
    label = current_label
  ) %>%
    mutate(
      status = "Estimated successfully.",
      n_studies_total = n_studies_total,
      n_studies_high = n_studies_high
    )
}

rob_sensitivity_results <- bind_rows(rob_sensitivity_results)

cat("\n--- SENSITIVITY ANALYSES EXCLUDING HIGH-RISK STUDIES ---\n")
print(rob_sensitivity_results, n = nrow(rob_sensitivity_results), width = Inf)


###########################################################
# 15. Create comparison table against primary overall ER model
###########################################################

primary_overall_model <- fit_overall_er_model(dat_er)

primary_overall_summary <- extract_intercept_robust_summary(
  model = primary_overall_model,
  data = dat_er,
  analysis_label = "Primary overall ER model",
  moderator = "Primary",
  label = "Primary overall ER model"
)

rob_sensitivity_comparison <- rob_sensitivity_results %>%
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
    tau2_level_3,
    tau2_level_2,
    n_studies_total,
    n_studies_high
  )

cat("\n--- PRIMARY OVERALL ER MODEL SUMMARY ---\n")
print(primary_overall_summary, width = Inf)


###########################################################
# 16. Plot high-risk exclusion sensitivity results
###########################################################

plot_data <- rob_sensitivity_results %>%
  filter(status == "Estimated successfully.")

if (nrow(plot_data) > 0) {
  
  sensitivity_plot <- ggplot(
    plot_data,
    aes(
      x = reorder(label, estimate),
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
      title = "Sensitivity analyses excluding high-risk studies",
      x = NULL,
      y = "Pooled Hedges' g (robust 95% CI)"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(output_folder, "sensitivity_excluding_high_risk_studies_plot.png"),
    plot = sensitivity_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
}


###########################################################
# 17. Save level counts and per-domain text outputs
###########################################################

for (x in rob_results) {
  
  if (!is.null(x$level_counts)) {
    write_csv(
      x$level_counts,
      file.path(
        output_folder,
        paste0("rob_level_counts_", x$moderator, ".csv")
      )
    )
  }
  
  coverage_tbl <- tibble(
    moderator = x$moderator,
    label = x$label,
    status = x$status,
    k_effect_sizes = nrow(x$data),
    n_studies = n_distinct(x$data$Study_ID),
    n_comparisons = n_distinct(x$data$Comparison_ID)
  )
  
  write_csv(
    coverage_tbl,
    file.path(output_folder, paste0("coverage_", x$moderator, ".csv"))
  )
  
  if (!is.null(x$model)) {
    save_model_output(
      summary(x$model),
      file.path(output_folder, paste0("model_summary_", x$moderator, ".txt")),
      paste0("Risk-of-bias moderator model: ", x$label)
    )
    
    save_model_output(
      x$robust,
      file.path(output_folder, paste0("robust_tests_", x$moderator, ".txt")),
      paste0("Cluster-robust tests: ", x$label)
    )
    
    save_model_output(
      x$omnibus,
      file.path(output_folder, paste0("omnibus_test_", x$moderator, ".txt")),
      paste0("Robust omnibus moderator test: ", x$label)
    )
  } else {
    save_model_output(
      coverage_tbl,
      file.path(output_folder, paste0("model_summary_", x$moderator, ".txt")),
      paste0("Risk-of-bias moderator model not estimated: ", x$label)
    )
  }
}


###########################################################
# 18. Save main outputs
###########################################################

write_csv(
  rob_overview,
  file.path(output_folder, "risk_of_bias_overview_ER.csv")
)

write_csv(
  rob_coefficients,
  file.path(output_folder, "risk_of_bias_coefficients_ER.csv")
)

write_csv(
  rob_robust_tests,
  file.path(output_folder, "risk_of_bias_robust_tests_ER.csv")
)

write_csv(
  primary_overall_summary,
  file.path(output_folder, "primary_overall_ER_model_for_rob_reference.csv")
)

write_csv(
  rob_sensitivity_results,
  file.path(output_folder, "sensitivity_excluding_high_risk_studies_ER.csv")
)

write_csv(
  rob_sensitivity_comparison,
  file.path(output_folder, "sensitivity_excluding_high_risk_studies_ER_summary.csv")
)


###########################################################
# 19. Save R objects for later scripts
###########################################################

save(
  dat_er,
  rob_specs,
  rob_results,
  rob_overview,
  rob_coefficients,
  rob_robust_tests,
  primary_overall_model,
  primary_overall_summary,
  rob_sensitivity_results,
  rob_sensitivity_comparison,
  file = file.path(output_folder, "risk_of_bias_analysis_ER.RData")
)


###########################################################
# 20. Final console output
###########################################################

cat("\n====================================================\n")
cat("Risk of bias analysis script finished.\n")
cat("Files saved in: data/derived/risk_of_bias_ER/\n")
cat("- risk_of_bias_overview_ER.csv\n")
cat("- risk_of_bias_coefficients_ER.csv\n")
cat("- risk_of_bias_robust_tests_ER.csv\n")
cat("- primary_overall_ER_model_for_rob_reference.csv\n")
cat("- sensitivity_excluding_high_risk_studies_ER.csv\n")
cat("- sensitivity_excluding_high_risk_studies_ER_summary.csv\n")
cat("- coverage_*.csv\n")
cat("- model_summary_*.txt\n")
cat("- robust_tests_*.txt\n")
cat("- omnibus_test_*.txt\n")
cat("- rob_level_counts_*.csv\n")
cat("- sensitivity_excluding_high_risk_studies_plot.png (if estimable)\n")
cat("- risk_of_bias_analysis_ER.RData\n")
cat("====================================================\n")