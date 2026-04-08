############################################################
# Script: 05_moderator_analysis_ER.R
# Project: ER meta-analysis_final
# Purpose:
#   Test a priori study-level moderators of the pooled
#   emotion regulation (ER) effect of mindfulness-based
#   interventions (MBIs).
#
# Main analysis dataset:
#   - Primary effect sizes computed with assumed pre-post
#     correlation r = .50
#   - Emotion regulation outcomes only
#
# Main modelling approach:
#   - Multilevel random-effects meta-analysis
#   - Effect sizes nested within studies, with a lower‑level
#     random effect for unique comparisons/effect rows
#   - Models estimated with REML
#   - Cluster-robust inference at the study level
#   - Each moderator tested in a separate model
#   - Outcome_Domain retained as an adjustment covariate
#
# A priori moderators tested:
#   - Population type
#   - Age_Mean
#   - MBI_Duration_Weeks
#   - Control_Type
#
# Input:
#   data/derived/effect_sizes/er_meta_effect_sizes_r50.csv
#
# Output:
#   data/derived/moderators_ER/
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
  "clubSandwich"
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


############################
# 2. Define file locations
############################

input_file <- "data/derived/effect_sizes/er_meta_effect_sizes_r50.csv"
output_folder <- "data/derived/moderators_ER"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


############################
# 3. Read effect size dataset
############################

dat <- read_csv(input_file, show_col_types = FALSE)

cat("\n----------------------------------------\n")
cat("Moderator analysis script started\n")
cat("----------------------------------------\n")
cat("Rows read:", nrow(dat), "\n")
cat("Columns read:", ncol(dat), "\n")


###########################################################
# 4. Check required variables are present
###########################################################

required_vars <- c(
  "Effect_ID", "Study_ID", "Comparison_ID", "Outcome_Domain",
  "yi", "vi", "analysis_set",
  "Population_collapsed",   # use collapsed population
  "Age_Mean",
  "MBI_Duration_Weeks",
  "Control_collapsed"       # use collapsed control
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
# 5. Prepare moderator analysis dataset
###########################################################

# Assumption:
# Missing or unclear study-level codes (e.g., 99) have already
# been converted to NA in earlier data-cleaning steps.

dat_er <- dat %>%
  filter(analysis_set == "ER") %>%
  mutate(
    Study_ID = as.factor(Study_ID),
    Comparison_ID = as.factor(Comparison_ID),
    Outcome_Domain = factor(
      Outcome_Domain,
      levels = c("ER_Total", "ER_Adaptive", "ER_Maladaptive")
    ),
    Population_collapsed = factor(Population_collapsed),
    Control_collapsed = factor(Control_collapsed)
  )
    
   
cat("\n--- MAIN MODERATOR DATASET ---\n")
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
# 7. Define moderator specifications
###########################################################

moderator_specs <- tribble(
  ~moderator,              ~type,        ~label,
  "Population_collapsed",  "categorical","Population type (collapsed)",
  "Age_Mean",              "continuous", "Mean age",
  "MBI_Duration_Weeks",    "continuous", "Intervention duration (weeks)",
  "Control_collapsed",     "categorical","Control type (collapsed)"
)

# Practical thresholds for estimation
min_total_studies_for_model <- 10
min_studies_per_factor_level <- 4


###########################################################
# 8. Prepare moderator-specific datasets
###########################################################

prepare_moderator_dataset <- function(data, moderator, type) {
  
  dat_mod <- data %>%
    filter(!is.na(.data[[moderator]]))
  
  if (type == "categorical") {
    
    dat_mod <- dat_mod %>%
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
    
    return(list(
      data = dat_mod,
      level_counts = level_counts_final
    ))
  }
  
  if (type == "continuous") {
    
    dat_mod <- dat_mod %>%
      mutate(
        moderator_value = as.numeric(.data[[moderator]])
      ) %>%
      filter(!is.na(moderator_value))
    
    return(list(
      data = dat_mod,
      level_counts = NULL
    ))
  }
  
  stop("Unknown moderator type.")
}


###########################################################
# 9. Fit single-moderator models
###########################################################

fit_single_moderator <- function(data, moderator, type, label) {
  
  prepared <- prepare_moderator_dataset(data, moderator, type)
  dat_mod <- prepared$data
  level_counts <- prepared$level_counts
  
  n_studies <- n_distinct(dat_mod$Study_ID)
  n_effects <- nrow(dat_mod)
  n_comparisons <- n_distinct(dat_mod$Comparison_ID)
  
  if (n_effects == 0) {
    return(list(
      moderator = moderator,
      label = label,
      type = type,
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
      type = type,
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
  
  if (type == "categorical" && n_distinct(dat_mod$moderator_value) < 2) {
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Not estimated: fewer than 2 factor levels after filtering sparse levels.",
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      omnibus = NULL
    ))
  }
  
  if (type == "continuous" && length(unique(dat_mod$moderator_value)) < 2) {
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Not estimated: moderator has no usable variation.",
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      omnibus = NULL
    ))
  }
  
  if (type == "continuous") {
    
    dat_mod <- dat_mod %>%
      mutate(
        moderator_value_c = moderator_value - mean(moderator_value, na.rm = TRUE)
      )
    
    model <- rma.mv(
      yi = yi,
      V = vi,
      mods = ~ Outcome_Domain + moderator_value_c,
      random = ~ 1 | Study_ID/Comparison_ID,
      data = dat_mod,
      method = "REML",
      test = "t"
    )
    
    robust <- coef_test(
      model,
      vcov = "CR2",
      cluster = dat_mod$Study_ID
    )
    
    omnibus <- Wald_test(
      model,
      constraints = constrain_zero(length(coef(model))),
      vcov = "CR2",
      cluster = dat_mod$Study_ID
    )
    
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Estimated successfully.",
      data = dat_mod,
      level_counts = level_counts,
      model = model,
      robust = robust,
      omnibus = omnibus
    ))
  }
  
  if (type == "categorical") {
    
    model <- rma.mv(
      yi = yi,
      V = vi,
      mods = ~ Outcome_Domain + moderator_value,
      random = ~ 1 | Study_ID/Comparison_ID,
      data = dat_mod,
      method = "REML",
      test = "t"
    )
    
    robust <- coef_test(
      model,
      vcov = "CR2",
      cluster = dat_mod$Study_ID
    )
    
    mod_coef_idx <- grep("^moderator_value", names(coef(model)))
    omnibus <- Wald_test(
      model,
      constraints = constrain_zero(mod_coef_idx),
      vcov = "CR2",
      cluster = dat_mod$Study_ID
    )
    
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Estimated successfully.",
      data = dat_mod,
      level_counts = level_counts,
      model = model,
      robust = robust,
      omnibus = omnibus
    ))
  }
}


###########################################################
# 10. Run all moderator analyses
###########################################################

moderator_results <- pmap(
  moderator_specs,
  function(moderator, type, label) {
    fit_single_moderator(
      data = dat_er,
      moderator = moderator,
      type = type,
      label = label
    )
  }
)

names(moderator_results) <- moderator_specs$moderator


###########################################################
# 11. Create overview table
###########################################################

moderator_overview <- map_dfr(moderator_results, function(x) {
  
  if (is.null(x$model)) {
    return(tibble(
      moderator = x$moderator,
      label = x$label,
      type = x$type,
      status = x$status,
      k_effect_sizes = nrow(x$data),
      n_studies = n_distinct(x$data$Study_ID),
      n_comparisons = n_distinct(x$data$Comparison_ID),
      estimate = NA_real_,
      se_model = NA_real_,
      ci_lb_model = NA_real_,
      ci_ub_model = NA_real_,
      p_value_model = NA_real_,
      omnibus_stat = NA_real_,
      omnibus_df_num = NA_real_,
      omnibus_df_den = NA_real_,
      omnibus_p = NA_real_,
      ci_lb_robust = NA_real_,
      ci_ub_robust = NA_real_
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
  
  coef_index <- if (x$type == "continuous") length(coef(x$model)) else NA_integer_
  
  estimate <- if (x$type == "continuous") as.numeric(coef(x$model)[coef_index]) else NA_real_
  se_model <- if (x$type == "continuous") x$model$se[coef_index] else NA_real_
  ci_lb_model <- if (x$type == "continuous") x$model$ci.lb[coef_index] else NA_real_
  ci_ub_model <- if (x$type == "continuous") x$model$ci.ub[coef_index] else NA_real_
  p_value_model <- if (x$type == "continuous") x$model$pval[coef_index] else NA_real_
  
  ci_lb_robust <- NA_real_
  ci_ub_robust <- NA_real_
  if (x$type == "continuous" && !is.null(x$robust)) {
    # Convert the robust test object into a data frame
    robust_df <- as.data.frame(x$robust)
    # Identify the last coefficient term (the moderator slope)
    slope_term  <- names(coef(x$model))[length(coef(x$model))]
    slope_row   <- robust_df[rownames(robust_df) == slope_term, ]
    # Compute critical t-value using the Satterthwaite df
    t_crit  <- qt(0.975, df = slope_row$df_Satt)
    # Calculate the robust CI
    ci_lb_robust <- slope_row$beta - t_crit * slope_row$SE
    ci_ub_robust <- slope_row$beta + t_crit * slope_row$SE
  }
  
  tibble(
    moderator = x$moderator,
    label = x$label,
    type = x$type,
    status = x$status,
    k_effect_sizes = x$model$k,
    n_studies = n_distinct(x$data$Study_ID),
    n_comparisons = n_distinct(x$data$Comparison_ID),
    estimate = estimate,
    se_model = se_model,
    ci_lb_model = ci_lb_model,
    ci_ub_model = ci_ub_model,
    p_value_model = p_value_model,
    omnibus_stat = omnibus_stat,
    omnibus_df_num = omnibus_df_num,
    omnibus_df_den = omnibus_df_den,
    omnibus_p = omnibus_p,
    ci_lb_robust = ci_lb_robust,
    ci_ub_robust = ci_ub_robust
  )
  })

cat("\n--- MODERATOR OVERVIEW ---\n")
print(moderator_overview, n = nrow(moderator_overview))


###########################################################
# 12. Create coefficient table
###########################################################

moderator_coefficients <- map_dfr(moderator_results, function(x) {
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
cat("\n--- MODERATOR COEFFICIENT TABLE ---\n")
print(moderator_coefficients, n = nrow(moderator_coefficients))

###########################################################
# 13. Create robust coefficient table
###########################################################

moderator_robust_tests <- map_dfr(moderator_results, function(x) {
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
  robust_tbl <- as_tibble(x$robust, rownames = "term")
  # compute critical t and CIs for each term
  t_crit <- qt(0.975, df = robust_tbl$df_Satt)
  ci_lb <- robust_tbl$beta - t_crit * robust_tbl$SE
  ci_ub <- robust_tbl$beta + t_crit * robust_tbl$SE
  tibble(
    moderator = x$moderator,
    label = x$label,
    term = robust_tbl$term,
    estimate_robust = robust_tbl$beta,
    se_robust = robust_tbl$SE,
    ci_lb_robust = ci_lb,
    ci_ub_robust = ci_ub,
    t_stat_robust = robust_tbl$tstat,
    df_robust = robust_tbl$df_Satt,
    p_value_robust = robust_tbl$p_Satt
  )
})

cat("\n--- ROBUST MODERATOR TEST TABLE ---\n")
print(moderator_robust_tests, n = nrow(moderator_robust_tests))

###########################################################
# 14. Subgroup pooled effects (overall + domain-specific)
###########################################################

# Purpose:
# To compute pooled effect sizes (Hedges' g) within each level
# of categorical moderators, both:
#   (1) overall across all ER domains
#   (2) separately for each ER domain
#
# This enables reporting of:
# - overall subgroup effects
# - domain-specific subgroup effects
# - alongside omnibus moderator tests
#
# Model structures:
# (1) Overall effect:
#     yi ~ 1
#
# (2) Domain-specific effects:
#     yi ~ Outcome_Domain - 1
#
# Key features:
# - Multilevel random-effects model
# - REML estimation
# - Cluster-robust inference (CR2)
#
# Output:
# - One row per moderator level × domain
# - Domain = "Overall", "ER_Total", "ER_Adaptive", "ER_Maladaptive"
###########################################################


############################
# 14.1 Function to extract subgroup effects
############################

extract_subgroup_effects_full <- function(data, moderator_name, label) {
  
  data <- data %>% filter(!is.na(.data[[moderator_name]]))
  
  data <- data %>%
    mutate(subgroup = as.factor(.data[[moderator_name]]))
  
  subgroups <- levels(data$subgroup)
  
  results <- map_dfr(subgroups, function(g) {
    
    dat_sub <- data %>% filter(subgroup == g)
    
    # -----------------------------
    # (1) Overall pooled effect
    # -----------------------------
    
    model_overall <- rma.mv(
      yi = yi,
      V = vi,
      mods = ~ 1,
      random = ~ 1 | Study_ID/Comparison_ID,
      data = dat_sub,
      method = "REML",
      test = "t"
    )
    
    robust_overall <- coef_test(
      model_overall,
      vcov = "CR2",
      cluster = dat_sub$Study_ID
    )
    
    row_overall <- as.data.frame(robust_overall)[1, ]
    t_crit_overall <- qt(0.975, df = row_overall$df_Satt)
    
    overall_row <- tibble(
      moderator = moderator_name,
      label = label,
      level = g,
      domain = "Overall",
      k_effects = nrow(dat_sub),
      n_studies = n_distinct(dat_sub$Study_ID),
      estimate = row_overall$beta,
      se = row_overall$SE,
      ci_lb = row_overall$beta - t_crit_overall * row_overall$SE,
      ci_ub = row_overall$beta + t_crit_overall * row_overall$SE,
      p_value = row_overall$p_Satt
    )
    
    # -----------------------------
    # (2) Domain-specific effects
    # -----------------------------
    
    model_domain <- rma.mv(
      yi = yi,
      V = vi,
      mods = ~ Outcome_Domain - 1,
      random = ~ 1 | Study_ID/Comparison_ID,
      data = dat_sub,
      method = "REML",
      test = "t"
    )
    
    robust_domain <- coef_test(
      model_domain,
      vcov = "CR2",
      cluster = dat_sub$Study_ID
    )
    
    robust_df <- as.data.frame(robust_domain)
    
    domain_rows <- map_dfr(rownames(robust_df), function(term) {
      
      row <- robust_df[term, ]
      t_crit <- qt(0.975, df = row$df_Satt)
      
      domain_name <- gsub("^Outcome_Domain", "", term)
      
      tibble(
        moderator = moderator_name,
        label = label,
        level = g,
        domain = domain_name,
        k_effects = nrow(dat_sub %>% filter(Outcome_Domain == domain_name)),
        n_studies = dat_sub %>%
          filter(Outcome_Domain == domain_name) %>%
          summarise(n = n_distinct(Study_ID)) %>%
          pull(n),
        estimate = row$beta,
        se = row$SE,
        ci_lb = row$beta - t_crit * row$SE,
        ci_ub = row$beta + t_crit * row$SE,
        p_value = row$p_Satt
      )
    })
    
    bind_rows(overall_row, domain_rows)
  })
  
  return(results)
}

############################
# 14.2 Run subgroup analyses
############################

categorical_mods <- moderator_specs %>%
  filter(type == "categorical")

subgroup_results_full <- map2_dfr(
  categorical_mods$moderator,
  categorical_mods$label,
  ~ extract_subgroup_effects_full(dat_er, .x, .y)
)

cat("\n--- SUBGROUP EFFECTS (OVERALL + DOMAIN) ---\n")
print(subgroup_results_full, n = nrow(subgroup_results_full))


############################
# 14.3 Combine with omnibus tests
############################

omnibus_results <- moderator_overview %>%
  filter(type == "categorical") %>%
  select(
    moderator,
    omnibus_stat,
    omnibus_df_num,
    omnibus_df_den,
    omnibus_p
  ) %>%
  distinct()

final_moderator_table <- subgroup_results_full %>%
  left_join(omnibus_results, by = "moderator") %>%
  arrange(moderator, level, domain)

cat("\n--- FINAL MODERATOR TABLE ---\n")
print(final_moderator_table, n = nrow(final_moderator_table))


############################
# 14.4 Save
############################

write_csv(
  final_moderator_table,
  file.path(output_folder, "moderator_subgroup_effects_FULL_ER.csv")
)


###########################################################
# 15. Domain-specific moderator tests via interactions
###########################################################

# Purpose:
# Test whether each moderator operates differently across
# ER_Total, ER_Adaptive, and ER_Maladaptive.
#
# Model structure:
# yi ~ Outcome_Domain * moderator
#
# Interpretation:
# - Continuous moderator:
#   tests whether the slope differs by outcome domain
# - Categorical moderator:
#   tests whether category differences differ by outcome domain

fit_domain_interaction_moderator <- function(data, moderator, type, label) {
  
  prepared <- prepare_moderator_dataset(data, moderator, type)
  dat_mod <- prepared$data
  level_counts <- prepared$level_counts
  
  n_studies <- n_distinct(dat_mod$Study_ID)
  n_effects <- nrow(dat_mod)
  n_comparisons <- n_distinct(dat_mod$Comparison_ID)
  
  if (n_effects == 0) {
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Not estimated: no usable rows after filtering missing values / sparse levels.",
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      interaction_test = NULL
    ))
  }
  
  if (n_studies < min_total_studies_for_model) {
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = paste0(
        "Not estimated: fewer than ", min_total_studies_for_model,
        " studies with usable data."
      ),
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      interaction_test = NULL
    ))
  }
  
  if (type == "categorical" && n_distinct(dat_mod$moderator_value) < 2) {
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Not estimated: fewer than 2 factor levels after filtering sparse levels.",
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      interaction_test = NULL
    ))
  }
  
  if (type == "continuous" && length(unique(dat_mod$moderator_value)) < 2) {
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Not estimated: moderator has no usable variation.",
      data = dat_mod,
      level_counts = level_counts,
      model = NULL,
      robust = NULL,
      interaction_test = NULL
    ))
  }
  
  if (type == "continuous") {
    
    dat_mod <- dat_mod %>%
      mutate(
        moderator_value_c = moderator_value - mean(moderator_value, na.rm = TRUE)
      )
    
    model <- rma.mv(
      yi = yi,
      V = vi,
      mods = ~ Outcome_Domain * moderator_value_c,
      random = ~ 1 | Study_ID/Comparison_ID,
      data = dat_mod,
      method = "REML",
      test = "t"
    )
    
    robust <- coef_test(
      model,
      vcov = "CR2",
      cluster = dat_mod$Study_ID
    )
    
    interaction_terms <- grep(
      "Outcome_Domain.*:moderator_value_c|moderator_value_c:Outcome_Domain",
      names(coef(model))
    )
    
    interaction_test <- if (length(interaction_terms) > 0) {
      Wald_test(
        model,
        constraints = constrain_zero(interaction_terms),
        vcov = "CR2",
        cluster = dat_mod$Study_ID
      )
    } else {
      NULL
    }
    
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Estimated successfully.",
      data = dat_mod,
      level_counts = level_counts,
      model = model,
      robust = robust,
      interaction_test = interaction_test
    ))
  }
  
  if (type == "categorical") {
    
    model <- rma.mv(
      yi = yi,
      V = vi,
      mods = ~ Outcome_Domain * moderator_value,
      random = ~ 1 | Study_ID/Comparison_ID,
      data = dat_mod,
      method = "REML",
      test = "t"
    )
    
    robust <- coef_test(
      model,
      vcov = "CR2",
      cluster = dat_mod$Study_ID
    )
    
    interaction_terms <- grep("Outcome_Domain.*:moderator_value|moderator_value.*:Outcome_Domain",
                              names(coef(model)))
    
    interaction_test <- if (length(interaction_terms) > 0) {
      Wald_test(
        model,
        constraints = constrain_zero(interaction_terms),
        vcov = "CR2",
        cluster = dat_mod$Study_ID
      )
    } else {
      NULL
    }
    
    return(list(
      moderator = moderator,
      label = label,
      type = type,
      status = "Estimated successfully.",
      data = dat_mod,
      level_counts = level_counts,
      model = model,
      robust = robust,
      interaction_test = interaction_test
    ))
  }
}

domain_interaction_results <- purrr::pmap(
  moderator_specs,
  function(moderator, type, label) {
    fit_domain_interaction_moderator(
      data = dat_er,
      moderator = moderator,
      type = type,
      label = label
    )
  }
)

names(domain_interaction_results) <- moderator_specs$moderator

cat("\n--- DOMAIN-SPECIFIC MODERATOR INTERACTION MODELS FITTED ---\n")
cat("Number of interaction models:", length(domain_interaction_results), "\n")

###########################################################
# Summarise domain-specific moderator interaction tests
###########################################################

domain_interaction_overview <- purrr::map_dfr(domain_interaction_results, function(x) {
  if (is.null(x$model)) {
    return(tibble(
      moderator = x$moderator,
      label = x$label,
      type = x$type,
      status = x$status,
      k_effect_sizes = nrow(x$data),
      n_studies = dplyr::n_distinct(x$data$Study_ID),
      n_comparisons = dplyr::n_distinct(x$data$Comparison_ID),
      interaction_stat = NA_real_,
      interaction_df_num = NA_real_,
      interaction_df_den = NA_real_,
      interaction_p = NA_real_
    ))
  }
  
  interaction_stat <- NA_real_
  interaction_df_num <- NA_real_
  interaction_df_den <- NA_real_
  interaction_p <- NA_real_
  
  if (!is.null(x$interaction_test)) {
    if ("Fstat" %in% names(x$interaction_test)) {
      interaction_stat <- x$interaction_test$Fstat
      interaction_df_num <- x$interaction_test$df_num
      interaction_df_den <- x$interaction_test$df_denom
      interaction_p <- x$interaction_test$p_val
    } else {
      interaction_stat <- x$interaction_test$QM
      interaction_df_num <- x$interaction_test$m
      interaction_df_den <- NA_real_
      interaction_p <- x$interaction_test$QMp
    }
  }
  
  tibble(
    moderator = x$moderator,
    label = x$label,
    type = x$type,
    status = x$status,
    k_effect_sizes = x$model$k,
    n_studies = dplyr::n_distinct(x$data$Study_ID),
    n_comparisons = dplyr::n_distinct(x$data$Comparison_ID),
    interaction_stat = interaction_stat,
    interaction_df_num = interaction_df_num,
    interaction_df_den = interaction_df_den,
    interaction_p = interaction_p
  )
})
cat("\n--- DOMAIN-SPECIFIC MODERATOR INTERACTION OVERVIEW ---\n")
print(domain_interaction_overview, n = nrow(domain_interaction_overview), width = Inf)
###########################################################
#  Extract coefficients from interaction models
###########################################################

domain_interaction_coefficients <- purrr::map_dfr(domain_interaction_results, function(x) {
  
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

cat("\n--- DOMAIN-SPECIFIC MODERATOR INTERACTION COEFFICIENTS ---\n")
print(domain_interaction_coefficients, n = nrow(domain_interaction_coefficients), width = Inf)


###########################################################
# Save interaction-model outputs
###########################################################

write_csv(
  domain_interaction_overview,
  file.path(output_folder, "domain_interaction_overview_ER.csv")
)

write_csv(
  domain_interaction_coefficients,
  file.path(output_folder, "domain_interaction_coefficients_ER.csv")
)

for (x in domain_interaction_results) {
  
  if (!is.null(x$model)) {
    
    save_model_output(
      summary(x$model),
      file.path(output_folder, paste0("interaction_model_summary_", x$moderator, ".txt")),
      paste0("Domain-specific moderator interaction model: ", x$label)
    )
    
    save_model_output(
      x$robust,
      file.path(output_folder, paste0("interaction_robust_tests_", x$moderator, ".txt")),
      paste0("Cluster-robust tests for domain-specific moderator interaction: ", x$label)
    )
    
    save_model_output(
      x$interaction_test,
      file.path(output_folder, paste0("interaction_omnibus_test_", x$moderator, ".txt")),
      paste0("Omnibus interaction test with Outcome_Domain: ", x$label)
    )
  }
}

###########################################################
# 15. Save level counts for categorical moderators
###########################################################

for (x in moderator_results) {
  if (!is.null(x$level_counts)) {
    write_csv(
      x$level_counts,
      file.path(
        output_folder,
        paste0("moderator_level_counts_", x$moderator, ".csv")
      )
    )
  }
}


###########################################################
# 16. Save per-moderator outputs
###########################################################

for (x in moderator_results) {
  
  coverage_tbl <- tibble(
    moderator = x$moderator,
    label = x$label,
    type = x$type,
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
      paste0("Moderator model: ", x$label)
    )
    
    save_model_output(
      x$robust,
      file.path(output_folder, paste0("robust_tests_", x$moderator, ".txt")),
      paste0("Cluster-robust tests: ", x$label)
    )
    
    save_model_output(
      x$omnibus,
      file.path(output_folder, paste0("omnibus_test_", x$moderator, ".txt")),
      paste0("Omnibus moderator test: ", x$label)
    )
    
  } else {
    
    save_model_output(
      coverage_tbl,
      file.path(output_folder, paste0("model_summary_", x$moderator, ".txt")),
      paste0("Moderator model not estimated: ", x$label)
    )
  }
}


###########################################################
# 17. Save summary tables and R objects
###########################################################

write_csv(
  moderator_overview,
  file.path(output_folder, "moderator_overview_ER.csv")
)

write_csv(
  moderator_coefficients,
  file.path(output_folder, "moderator_coefficients_ER.csv")
)

write_csv(
  moderator_robust_tests,
  file.path(output_folder, "moderator_robust_tests_ER.csv")
)

save(
  dat_er,
  moderator_specs,
  moderator_results,
  moderator_overview,
  moderator_coefficients,
  moderator_robust_tests,
  domain_interaction_results,
  domain_interaction_overview,
  domain_interaction_coefficients,
  file = file.path(output_folder, "moderator_analysis_ER.RData")
)
###########################################################
# 18. Final console output
###########################################################

cat("\n====================================================\n")
cat("Moderator analysis script finished.\n")
cat("Files saved in: data/derived/moderators_ER/\n")
cat("- moderator_overview_ER.csv\n")
cat("- moderator_coefficients_ER.csv\n")
cat("- moderator_robust_tests_ER.csv\n")
cat("- coverage_*.csv\n")
cat("- model_summary_*.txt\n")
cat("- robust_tests_*.txt\n")
cat("- omnibus_test_*.txt\n")
cat("- moderator_level_counts_*.csv (for categorical moderators)\n")
cat("- moderator_analysis_ER.RData\n")
cat("- domain_interaction_overview_ER.csv\n")
cat("- domain_interaction_coefficients_ER.csv\n")
cat("- interaction_model_summary_*.txt\n")
cat("- interaction_robust_tests_*.txt\n")
cat("- interaction_omnibus_test_*.txt\n")
cat("====================================================\n")

library(dplyr)
library(stringr)

final_table_pub <- final_moderator_table_all %>%
  
  # Rename columns
  rename(
    Moderator = label,
    Category = level,
    `Outcome domain` = domain,
    `k (effects)` = k_effects,
    `k (studies)` = n_studies,
    Effect = estimate,
    SE = se,
    `CI lower` = ci_lb,
    `CI upper` = ci_ub,
    p = p_value
  ) %>%
  
  # Label continuous moderators clearly
  mutate(
    Category = ifelse(is.na(Category), "Linear effect", Category),
    
    # Clean domain labels
    `Outcome domain` = recode(
      `Outcome domain`,
      "Overall" = "Overall",
      "ER_Total" = "ER total",
      "ER_Adaptive" = "Adaptive strategies",
      "ER_Maladaptive" = "Maladaptive strategies"
    ),
    
    # Combine CI
    `95% CI` = paste0(
      "[", round(`CI lower`, 2), ", ", round(`CI upper`, 2), "]"
    ),
    
    # Round values
    Effect = round(Effect, 2),
    SE = round(SE, 2),
    p = ifelse(p < .001, "< .001", sprintf("%.3f", p))
  ) %>%
  
  # Keep only needed columns
  select(
    Moderator,
    Category,
    `Outcome domain`,
    `k (effects)`,
    `k (studies)`,
    Effect,
    `95% CI`,
    p
  ) %>%
  
  # Order rows nicely
  mutate(
    `Outcome domain` = factor(
      `Outcome domain`,
      levels = c("Overall", "ER total", "Adaptive strategies", "Maladaptive strategies")
    )
  ) %>%
  arrange(Moderator, Category, `Outcome domain`)

final_table_pub <- final_table_pub %>%
  mutate(
    Category = ifelse(
      duplicated(Moderator),
      paste0("  ", Category),
      Category
    )
  )

write.csv(final_table_pub, "moderator_table_CPR.csv", row.names = FALSE)