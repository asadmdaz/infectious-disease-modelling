#################################################################################################################
# ============================================================
# Tables and figures for SEIR dynamic network paper
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(forcats)
library(kableExtra)
library(purrr)

mle_all <- read.csv("output/fit_results_all_parameter_Complete_MLE.csv", header = T)

truth_map <- tibble(
  parameter = c(
    "beta","xi","kappa","gamma",
    paste0("eta_", c("SS","SE","SI","SR","EE","EI","ER","II","IR","RR")),
    paste0("omega_", c("SS","SE","SI","SR","EE","EI","ER","II","IR","RR")),
    "p_E","p_I","s","c"
  ),
  truth = c(
    0.30, 0.05, 0.40, 0.25,
    rep(0.08, 10), rep(c(0.03,0.03,0.06,0.03,0.03,0.03,0.03,0.07,0.03,0.03), 1),
    0.60, 0.85, 0.90, 0.95
  )
)

mle_all <- mle_all %>%
  left_join(truth_map, by = "parameter")

mle_summary <- mle_all %>%
  group_by(parameter) %>%
  summarise(
    truth = first(truth),
    estimate = mean(fit_summary.estimate, na.rm = TRUE),
    se = mean(fit_summary.se, na.rm = TRUE),
    lower = mean(fit_summary.lower, na.rm = TRUE),
    upper = mean(fit_summary.upper, na.rm = TRUE),
    bias = mean(fit_summary.estimate - truth, na.rm = TRUE),
    rmse = sqrt(mean((fit_summary.estimate - truth)^2, na.rm = TRUE)),
    coverage = mean(fit_summary.lower <= truth & truth <= fit_summary.upper, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(mle_summary, "output/mle_summary.csv")

# Table 1
tab_mle <- mle_summary %>%
  mutate(
    estimate = round(estimate, 4),
    se = round(se, 4),
    lower = round(lower, 4),
    upper = round(upper, 4),
    bias = round(bias, 4),
    rmse = round(rmse, 4),
    coverage = round(coverage, 3)
  )

kable(
  tab_mle,
  format = "latex",
  booktabs = TRUE,
  caption = "Complete-data maximum likelihood estimates and simulation performance.",
  label = "tab:mle_summary"
) %>%
  kable_styling(latex_options = c("hold_position", "striped"))

mle_summary <- mle_summary %>%
  mutate(parameter = as.factor(parameter))

p_bias <- mle_summary %>%
  mutate(parameter = reorder(parameter, bias)) %>%
  ggplot(aes(x = parameter, y = bias)) +
  geom_col(fill = "#2c7fb8", width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Bias")

ggsave("output/mle_bias.png", p_bias, width = 10, height = 8, dpi = 300)

p_rmse <- mle_summary %>%
  mutate(parameter = fct_reorder(parameter, rmse)) %>%
  ggplot(aes(x = parameter, y = rmse)) +
  geom_col(fill = "#41ab5d", width = 0.7) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "RMSE")

ggsave("output/mle_rmse.png", p_rmse, width = 10, height = 8, dpi = 300)


p_cov <- mle_summary %>%
  mutate(parameter = fct_reorder(parameter, coverage)) %>%
  ggplot(aes(x = parameter, y = coverage)) +
  geom_col(fill = "#f03b20", width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Coverage")

ggsave("output/mle_coverage.png", p_cov, width = 10, height = 8, dpi = 300)



### MCMC summary ###
post <- read_csv("output/fit_results_all_parameter_bayesian_posterior_draws.csv", show_col_types = FALSE)
post <- post[,-1]

library(dplyr)
library(tidyr)
library(readr)
library(posterior)

# ------------------------------------------------------------
# Read posterior draws
# Expected columns:
# posterior_draws.beta
# posterior_draws.xi
# posterior_draws.kappa
# posterior_draws.gamma
# posterior_draws.p_E
# posterior_draws.p_I
# posterior_draws.s
# posterior_draws.c
# posterior_draws.eta_SS
# posterior_draws.omega_SS
# ...
# ------------------------------------------------------------
post <- read_csv("output/fit_results_all_parameter_bayesian_posterior_draws.csv", show_col_types = FALSE)
post <- post[,-1]

# Remove prefix if present
names(post) <- sub("^posterior_draws\\.", "", names(post))

# Keep only numeric parameter columns
post <- post %>% select(where(is.numeric))

# ------------------------------------------------------------
# Summarise posterior draws
# ------------------------------------------------------------
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
  arrange(parameter)

write_csv(bayes_summary, "output/bayes_summary1.csv")


#################################################################################################################
# ============================================================
# Bayesian MCMC summaries and plots
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(knitr)
  library(kableExtra)
  library(scales)
})

# dir.create("output", showWarnings = FALSE, recursive = TRUE)
# dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
# dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Read posterior draws
# ------------------------------------------------------------
post <- read_csv("output/fit_results_all_parameter_bayesian_posterior_draws.csv", show_col_types = FALSE)
names(post) <- sub("^posterior_draws\\.", "", names(post))

# Keep only parameter columns
post <- post %>% select(where(is.numeric))
post <- post[,-1]

# ------------------------------------------------------------
# 2. Summarise posterior draws
# ------------------------------------------------------------
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
  arrange(parameter)

write.csv(bayes_summary, "output/bayes_summary2.csv")

# Add truth and performance
truth_map <- tibble::tribble(
  ~parameter, ~truth,
  "beta", 0.30,
  "xi", 0.05,
  "kappa", 0.40,
  "gamma", 0.25,
  "p_E", 0.60,
  "p_I", 0.85,
  "s", 0.90,
  "c", 0.95
)

bayes_summary <- bayes_summary %>%
  left_join(truth_map, by = "parameter") %>%
  mutate(
    bias = mean - truth,
    abs_error = abs(mean - truth),
    cover = lower <= truth & truth <= upper
  )

write_csv(bayes_summary, "output/bayes_summary_with_truth.csv")

######################################################################
# Table for paper
bayes_summary <- read_csv("output/bayes_summary_with_truth.csv", show_col_types = FALSE)

tab_bayes <- bayes_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)),
         cover = ifelse(cover, "Yes", "No"))

kable(
  tab_bayes,
  format = "latex",
  booktabs = TRUE,
  caption = "Posterior summaries from Bayesian MCMC",
  label = "tab:bayes_summary"
) %>%
  kable_styling(latex_options = c("hold_position", "striped"))

# Truth vs estimate plot
p_truth <- bayes_summary %>%
  filter(!is.na(truth)) %>%
  ggplot(aes(x = truth, y = mean)) +
  geom_point(size = 2, alpha = 0.7, color = "#2c7fb8") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ parameter, scales = "free") +
  theme_classic(base_size = 12) +
  labs(x = "True value", y = "Posterior mean") +
  theme(strip.background = element_blank())

ggsave("output/bayes_truth_vs_mean.png", p_truth, width = 12, height = 8, dpi = 300)

# Density plots
post_long <- post %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "draw")

p_density <- post_long %>%
  ggplot(aes(x = draw, fill = parameter)) +
  geom_density(alpha = 0.25, color = NA) +
  facet_wrap(~ parameter, scales = "free", ncol = 4) +
  theme_classic(base_size = 11) +
  labs(x = "Posterior draw", y = "Density") +
  theme(
    legend.position = "none",
    strip.background = element_blank()
  )

ggsave("output/bayes_density.png", p_density, width = 14, height = 10, dpi = 300)


####################################################################
truth_map <- tibble::tribble(
  ~parameter, ~truth,
  "beta", 0.30,
  "xi", 0.05,
  "kappa", 0.40,
  "gamma", 0.25,
  "p_E", 0.60,
  "p_I", 0.85,
  "s", 0.90,
  "c", 0.95
)

bayes_summary <- read_csv("output/bayes_summary2.csv", show_col_types = FALSE)

bayes_summary = bayes_summary %>%
  left_join(truth_map, by = "parameter") %>%
  mutate(
    bias = posterior_mean - truth,
    coverage = ifelse(lower <= truth & truth <= upper, 1, 0)
  ) %>%
  select(parameter, truth, posterior_mean, posterior_median, posterior_sd, lower, upper, bias, coverage) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

latex_tab <- kable(
  bayes_summary,
  format = "latex",
  booktabs = TRUE,
  caption = "Posterior summaries for the eight main parameters based on the true simulation values.",
  label = "tab:bayes_summary"
) %>%
  kable_styling(latex_options = c("hold_position", "striped"), font_size = 9)

save_kable(latex_tab, "output/tab_bayes_summary.tex")
#----------------------------------------------------------------------------------------------------------------
#################################################################################################################