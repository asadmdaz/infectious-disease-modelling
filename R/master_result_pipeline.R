#################################################################################################################
# ============================================================
# MASTER RESULTS PIPELINE
# SEIR + Dynamic Network Paper
# Produces tables and figures for MLE and Bayesian MCMC results
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(knitr)
  library(kableExtra)
  library(scales)
})

set.seed(20260603)
dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# USER SETTINGS
# ============================================================
complete_mle_dir <- "output"
posterior_file   <- "output/fit_results_all_parameter_bayesian_posterior_draws.csv"

# Truth values used in simulation
truth_map <- tibble::tribble(
  ~parameter, ~truth,
  "beta", 0.30, "xi", 0.05, "kappa", 0.40, "gamma", 0.25,
  "p_E", 0.60, "p_I", 0.85, "s", 0.90, "c", 0.95,
  "eta_SS", 0.08, "omega_SS", 0.03,
  "eta_SE", 0.08, "omega_SE", 0.03,
  "eta_SI", 0.05, "omega_SI", 0.06,
  "eta_SR", 0.08, "omega_SR", 0.03,
  "eta_EE", 0.08, "omega_EE", 0.03,
  "eta_EI", 0.08, "omega_EI", 0.03,
  "eta_ER", 0.08, "omega_ER", 0.03,
  "eta_II", 0.04, "omega_II", 0.07,
  "eta_IR", 0.08, "omega_IR", 0.03,
  "eta_RR", 0.08, "omega_RR", 0.03
)

# ============================================================
# HELPERS
# ============================================================
parse_complete_mle <- function(file) {
  x <- read_csv(file, show_col_types = FALSE)
  
  if (!any(startsWith(names(x), "fit_summary."))) return(NULL)
  
  fit <- x %>%
    select(starts_with("fit_summary.")) %>%
    rename_with(~ str_remove(.x, "^fit_summary\\.")) %>%
    distinct()
  
  fit$source_file <- basename(file)
  fit$replicate <- str_extract(basename(file), "\\d+")
  fit
}

safe_reorder <- function(df, value_col) {
  ord <- df %>%
    filter(!is.na(.data[[value_col]]), !is.na(parameter)) %>%
    arrange(.data[[value_col]]) %>%
    pull(parameter)
  df %>% mutate(parameter = factor(parameter, levels = unique(ord)))
}

latex_table <- function(df, caption, label, align = NULL) {
  if (is.null(align)) align <- paste(rep("l", ncol(df)), collapse = "")
  kbl(df, format = "latex", booktabs = TRUE, caption = caption, label = label, align = align) %>%
    kable_styling(latex_options = c("hold_position", "striped"), font_size = 9)
}

# ============================================================
# 1) COMPLETE MLE RESULTS
# ============================================================
mle_files <- list.files(complete_mle_dir, pattern = "\\.csv$", full.names = TRUE)

mle_all <- map_dfr(mle_files, parse_complete_mle)

if (nrow(mle_all) > 0) {
  mle_all <- mle_all %>%
    left_join(truth_map, by = "parameter") %>%
    mutate(
      estimate = as.numeric(estimate),
      lower = as.numeric(lower),
      upper = as.numeric(upper)
    )
  
  mle_summary <- mle_all %>%
    group_by(parameter) %>%
    summarise(
      truth = first(truth),
      estimate = mean(estimate, na.rm = TRUE),
      se = mean(se, na.rm = TRUE),
      lower = mean(lower, na.rm = TRUE),
      upper = mean(upper, na.rm = TRUE),
      bias = mean(estimate - truth, na.rm = TRUE),
      rmse = sqrt(mean((estimate - truth)^2, na.rm = TRUE)),
      coverage = mean(lower <= truth & truth <= upper, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(parameter)
  
  write_csv(mle_all, "output/tables/mle_replicate_results.csv")
  write_csv(mle_summary, "output/tables/mle_summary.csv")
  
  tab_mle <- mle_summary %>%
    mutate(
      across(where(is.numeric), ~ round(.x, 4)),
      coverage = round(coverage, 3)
    )
  save_kable(
    latex_table(tab_mle, "Complete-data MLE summaries.", "tab:mle_summary"),
    "output/tables/tab_mle_summary.tex"
  )
  
  p_mle_truth <- mle_all %>%
    filter(!is.na(truth)) %>%
    ggplot(aes(x = truth, y = estimate)) +
    geom_point(alpha = 0.55, size = 1.7, color = "#2c7fb8") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray35") +
    facet_wrap(~ parameter, scales = "free") +
    theme_classic(base_size = 12) +
    labs(x = "True value", y = "MLE") +
    theme(strip.background = element_blank())
  
  ggsave("output/figures/mle_truth_vs_estimate.png", p_mle_truth, width = 12, height = 8, dpi = 300)
  
  p_bias <- mle_summary %>%
    filter(!is.na(bias)) %>%
    mutate(parameter = factor(parameter, levels = parameter[order(bias)])) %>%
    ggplot(aes(x = parameter, y = bias)) +
    geom_col(fill = "#2c7fb8", width = 0.7) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = "Bias")
  
  ggsave("output/figures/mle_bias.png", p_bias, width = 10, height = 8, dpi = 300)
  
  p_rmse <- mle_summary %>%
    filter(!is.na(rmse)) %>%
    mutate(parameter = factor(parameter, levels = parameter[order(rmse)])) %>%
    ggplot(aes(x = parameter, y = rmse)) +
    geom_col(fill = "#41ab5d", width = 0.7) +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = "RMSE")
  
  ggsave("output/figures/mle_rmse.png", p_rmse, width = 10, height = 8, dpi = 300)
  
  p_cov <- mle_summary %>%
    filter(!is.na(coverage)) %>%
    mutate(parameter = factor(parameter, levels = parameter[order(coverage)])) %>%
    ggplot(aes(x = parameter, y = coverage)) +
    geom_col(fill = "#f03b20", width = 0.7) +
    coord_flip() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = "Coverage")
  
  ggsave("output/figures/mle_coverage.png", p_cov, width = 10, height = 8, dpi = 300)
}

# ============================================================
# 2) BAYESIAN POSTERIOR SUMMARIES
# ============================================================
if (file.exists(posterior_file)) {
  post <- read_csv(posterior_file, show_col_types = FALSE)
  names(post) <- str_remove(names(post), "^posterior_draws\\.")
  post <- post %>% select(where(is.numeric))
  
  bayes_summary <- post %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "draw") %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(draw, na.rm = TRUE),
      median = median(draw, na.rm = TRUE),
      sd = sd(draw, na.rm = TRUE),
      mcse = sd(draw, na.rm = TRUE) / sqrt(sum(!is.na(draw))),
      lower = quantile(draw, 0.025, na.rm = TRUE),
      upper = quantile(draw, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(truth_map, by = "parameter") %>%
    mutate(
      bias = mean - truth,
      rmse = abs(mean - truth),
      cover = lower <= truth & truth <= upper
    ) %>%
    arrange(parameter)
  
  write_csv(bayes_summary, "output/tables/bayes_summary.csv")
  
  tab_bayes <- bayes_summary %>%
    mutate(
      across(where(is.numeric), ~ round(.x, 4)),
      cover = ifelse(cover, "Yes", "No")
    )
  
  save_kable(
    latex_table(tab_bayes, "Posterior summaries from Bayesian MCMC.", "tab:bayes_summary"),
    "output/tables/tab_bayes_summary.tex"
  )
  
  p_bayes_truth <- bayes_summary %>%
    filter(!is.na(truth)) %>%
    ggplot(aes(x = truth, y = mean)) +
    geom_point(size = 2, alpha = 0.7, color = "#7b3294") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01, alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray35") +
    facet_wrap(~ parameter, scales = "free") +
    theme_classic(base_size = 12) +
    labs(x = "True value", y = "Posterior mean") +
    theme(strip.background = element_blank())
  
  ggsave("output/figures/bayes_truth_vs_mean.png", p_bayes_truth, width = 12, height = 8, dpi = 300)
  
  post_long <- post %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "draw")
  
  p_density <- post_long %>%
    ggplot(aes(x = draw, fill = parameter)) +
    geom_density(alpha = 0.25, color = NA) +
    facet_wrap(~ parameter, scales = "free", ncol = 4) +
    theme_classic(base_size = 11) +
    labs(x = "Posterior draw", y = "Density") +
    theme(legend.position = "none", strip.background = element_blank())
  
  ggsave("output/figures/bayes_density.png", p_density, width = 14, height = 10, dpi = 300)
}

# ============================================================
# 3) COMBINED MLE vs BAYES COMPARISON
# ============================================================
if (exists("mle_summary") && exists("bayes_summary")) {
  comp <- mle_summary %>%
    select(parameter, mle_estimate = estimate, mle_lower = lower, mle_upper = upper, mle_bias = bias, mle_rmse = rmse, mle_coverage = coverage) %>%
    full_join(
      bayes_summary %>% select(parameter, bayes_mean = mean, bayes_lower = lower, bayes_upper = upper, bayes_bias = bias, bayes_rmse = rmse, bayes_cover = cover),
      by = "parameter"
    ) %>%
    left_join(truth_map, by = "parameter")
  
  write_csv(comp, "output/tables/mle_bayes_comparison.csv")
  
  tab_comp <- comp %>%
    mutate(
      across(where(is.numeric), ~ round(.x, 4)),
      bayes_cover = ifelse(is.na(bayes_cover), NA, ifelse(bayes_cover, "Yes", "No"))
    )
  
  save_kable(
    latex_table(tab_comp, "Comparison of complete-data MLE and Bayesian posterior summaries.", "tab:comparison"),
    "output/tables/tab_mle_bayes_comparison.tex"
  )
  
  p_compare <- comp %>%
    filter(!is.na(truth)) %>%
    select(parameter, truth, mle_estimate, bayes_mean) %>%
    pivot_longer(c(mle_estimate, bayes_mean), names_to = "method", values_to = "estimate") %>%
    mutate(method = recode(method, mle_estimate = "MLE", bayes_mean = "Bayes")) %>%
    ggplot(aes(x = truth, y = estimate, color = method)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray35") +
    facet_wrap(~ parameter, scales = "free") +
    theme_classic(base_size = 12) +
    labs(x = "True value", y = "Estimate", color = "Method") +
    theme(strip.background = element_blank())
  
  ggsave("output/figures/mle_bayes_compare.png", p_compare, width = 12, height = 8, dpi = 300)
}

# ============================================================
# 4) SESSION INFO
# ============================================================
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
#################################################################################################################