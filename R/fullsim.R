#######################################################################################################
# ============================================================
# FULL MODEL-SPECIFIC SIMULATION + BAYESIAN MCMC PIPELINE
# SEIR on dynamic status-dependent contact network
# High / Moderate / Sparse observation regimes
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

set.seed(20260607)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/sim_results", showWarnings = FALSE, recursive = TRUE)
dir.create("output/mcmc", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SETTINGS
# ============================================================
R_reps <- 500
regimes <- c("High", "Moderate", "Sparse")

epi_params <- c("beta", "xi", "kappa", "gamma")
obs_params <- c("p_E", "p_I", "s", "c")

net_pairs <- c("SS","SE","SI","SR","EE","EI","ER","II","IR","RR")
net_form_params <- paste0("eta_", net_pairs)
net_diss_params <- paste0("tau_", net_pairs)
net_params <- c(net_form_params, net_diss_params)

all_params <- c(epi_params, obs_params, net_params)

truth_map <- tibble::tribble(
  ~parameter, ~truth,
  "beta", 0.30,
  "xi", 0.05,
  "kappa", 0.40,
  "gamma", 0.25,
  "p_E", 0.60,
  "p_I", 0.85,
  "s", 0.90,
  "c", 0.95,
  "eta_SS", 0.08, "tau_SS", 0.03,
  "eta_SE", 0.08, "tau_SE", 0.03,
  "eta_SI", 0.05, "tau_SI", 0.06,
  "eta_SR", 0.08, "tau_SR", 0.03,
  "eta_EE", 0.08, "tau_EE", 0.03,
  "eta_EI", 0.08, "tau_EI", 0.03,
  "eta_ER", 0.08, "tau_ER", 0.03,
  "eta_II", 0.04, "tau_II", 0.07,
  "eta_IR", 0.08, "tau_IR", 0.03,
  "eta_RR", 0.08, "tau_RR", 0.03
)

# ============================================================
# COMPLETE-DATA LOG-LIKELIHOOD
# ============================================================
log_complete_data_lik <- function(theta, dat) {
  beta  <- theta["beta"]
  xi    <- theta["xi"]
  kappa <- theta["kappa"]
  gamma <- theta["gamma"]
  p_E   <- theta["p_E"]
  p_I   <- theta["p_I"]
  s     <- theta["s"]
  c     <- theta["c"]
  
  st <- dat$stats
  
  epi_terms <- st$N_SE_int * log(beta) - beta * st$U_SI +
    st$N_SE_ext * log(xi) - xi * st$U_S +
    st$N_EI * log(kappa) - kappa * st$U_E +
    st$N_IR * log(gamma) - gamma * st$U_I
  
  obs_terms <- st$R_E * log(p_E) + (st$M_E - st$R_E) * log(1 - p_E) +
    st$R_I * log(p_I) + (st$M_I - st$R_I) * log(1 - p_I) +
    st$R_1 * log(s) + (st$M_1 - st$R_1) * log(1 - s) +
    st$R_0 * log(c) + (st$M_0 - st$R_0) * log(1 - c)
  
  net_terms <- 0
  if (!is.null(dat$net_stats)) {
    for (nm in names(dat$net_stats)) {
      ns <- dat$net_stats[[nm]]
      eta <- theta[paste0("eta_", nm)]
      omg <- theta[paste0("tau_", nm)]
      net_terms <- net_terms +
        ns$N01 * log(eta) - eta * ns$V0 +
        ns$N10 * log(omg) - omg * ns$V1
    }
  }
  
  if (any(!is.finite(c(epi_terms, obs_terms, net_terms)))) return(-Inf)
  epi_terms + obs_terms + net_terms
}

log_prior <- function(theta) {
  lp <- 0
  lp <- lp + dgamma(theta["beta"], 2, 10, log = TRUE)
  lp <- lp + dgamma(theta["xi"], 2, 40, log = TRUE)
  lp <- lp + dgamma(theta["kappa"], 2, 5, log = TRUE)
  lp <- lp + dgamma(theta["gamma"], 2, 8, log = TRUE)
  
  lp <- lp + dbeta(theta["p_E"], 6, 4, log = TRUE)
  lp <- lp + dbeta(theta["p_I"], 8, 2, log = TRUE)
  lp <- lp + dbeta(theta["s"], 18, 2, log = TRUE)
  lp <- lp + dbeta(theta["c"], 20, 1, log = TRUE)
  
  for (nm in net_form_params) lp <- lp + dgamma(theta[nm], 2, 20, log = TRUE)
  for (nm in net_diss_params) lp <- lp + dgamma(theta[nm], 2, 20, log = TRUE)
  
  lp
}

log_post <- function(theta, dat) {
  lp <- log_prior(theta)
  ll <- log_complete_data_lik(theta, dat)
  if (!is.finite(lp) || !is.finite(ll)) return(-Inf)
  lp + ll
}

# ============================================================
# MODEL-SPECIFIC SIMULATION PLACEHOLDER
# Replace with your exact generator
# ============================================================
simulate_dataset <- function(regime, rep_id) {
  truth_params <- c(
    beta = 0.30, xi = 0.05, kappa = 0.40, gamma = 0.25,
    p_E = 0.60, p_I = 0.85, s = 0.90, c = 0.95,
    setNames(rep(0.08, length(net_pairs)), paste0("eta_", net_pairs)),
    setNames(rep(0.03, length(net_pairs)), paste0("tau_", net_pairs))
  )
  
  stats <- list(
    N_SE_int = 20, N_SE_ext = 8, N_EI = 18, N_IR = 15,
    U_SI = 65, U_S = 120, U_E = 45, U_I = 55,
    M_E = 100, R_E = 60, M_I = 100, R_I = 84,
    M_1 = 800, R_1 = 720, M_0 = 1200, R_0 = 1140
  )
  
  net_stats <- setNames(vector("list", length(net_pairs)), net_pairs)
  for (nm in net_pairs) {
    net_stats[[nm]] <- list(
      N01 = 30, N10 = 28,
      V0 = 400, V1 = 380
    )
  }
  
  list(
    regime = regime,
    replicate = rep_id,
    truth = truth_params,
    stats = stats,
    net_stats = net_stats,
    Z = NULL,
    A = NULL
  )
}

# ============================================================
# LATENT PATH UPDATES
# Replace with exact local proposals for your CTMC path
# ============================================================
update_epidemic_path <- function(Z, A, theta, dat) {
  # Suggestion:
  # - sample a node i
  # - propose a local shift in infection time or incubation time
  # - enforce S -> E -> I -> R ordering
  # - recompute sufficient stats in dat$stats
  # - accept/reject by MH using log_post
  #
  # Return updated Z
  Z
}

update_network_path <- function(Z, A, theta, dat) {
  # Suggestion:
  # - sample a dyad (i, j)
  # - propose edge event-time perturbation or local insertion/deletion
  # - preserve consistency with node states
  # - recompute dat$net_stats
  # - accept/reject by MH using log_post
  #
  # Return updated A
  A
}

# ============================================================
# POSTERIOR SUMMARIZATION
# ============================================================
summarize_draws <- function(draws, truth_map) {
  draws %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      lower = quantile(value, 0.025, na.rm = TRUE),
      upper = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(truth_map, by = "parameter") %>%
    mutate(
      bias = mean - truth,
      rmse = abs(mean - truth),
      coverage = lower <= truth & truth <= upper,
      width = upper - lower
    )
}

# ============================================================
# MCMC FITTING LOOP
# ============================================================
fit_bayesian_mcmc <- function(dat, n_iter = 12000, burn = 4000, thin = 4) {
  theta <- c(
    beta = 0.20, xi = 0.04, kappa = 0.35, gamma = 0.20,
    p_E = 0.50, p_I = 0.80, s = 0.85, c = 0.90
  )
  
  for (nm in net_pairs) {
    theta[paste0("eta_", nm)] <- 0.05
    theta[paste0("tau_", nm)] <- 0.05
  }
  
  keep_idx <- seq(burn + 1, n_iter, by = thin)
  draws <- matrix(NA_real_, nrow = length(keep_idx), ncol = length(theta))
  colnames(draws) <- names(theta)
  
  Z <- dat$Z
  A <- dat$A
  
  prop_sd <- setNames(rep(0.02, length(theta)), names(theta))
  prop_sd[c("beta", "xi", "kappa", "gamma")] <- c(0.03, 0.01, 0.03, 0.02)
  prop_sd[c("p_E", "p_I", "s", "c")] <- c(0.02, 0.02, 0.01, 0.01)
  
  keep <- 1
  for (m in 1:n_iter) {
    Z <- update_epidemic_path(Z, A, theta, dat)
    A <- update_network_path(Z, A, theta, dat)
    
    for (nm in names(theta)) {
      cur <- theta
      prop <- cur
      prop[nm] <- rnorm(1, cur[nm], prop_sd[nm])
      
      if (nm %in% c("beta", "xi", "kappa", "gamma", net_form_params, net_diss_params)) {
        if (prop[nm] <= 0) next
      }
      if (nm %in% c("p_E", "p_I", "s", "c")) {
        if (prop[nm] <= 0 || prop[nm] >= 1) next
      }
      
      lp_cur <- log_post(cur, dat)
      lp_prop <- log_post(prop, dat)
      if (log(runif(1)) < lp_prop - lp_cur) theta <- prop
    }
    
    if (m %in% keep_idx) {
      draws[keep, ] <- theta
      keep <- keep + 1
    }
  }
  
  as_tibble(draws)
}

# ============================================================
# REPLICATE DRIVER
# ============================================================
run_one_rep <- function(regime, rep_id) {
  dat <- simulate_dataset(regime, rep_id)
  draws <- fit_bayesian_mcmc(dat)
  
  summ <- summarize_draws(draws, truth_map) %>%
    mutate(regime = regime, replicate = rep_id)
  
  dir.create(file.path("output/mcmc", tolower(regime)), showWarnings = FALSE, recursive = TRUE)
  write_csv(draws, file.path("output/mcmc", tolower(regime), paste0("draws_rep_", rep_id, ".csv")))
  write_csv(summ, file.path("output/mcmc", tolower(regime), paste0("summary_rep_", rep_id, ".csv")))
  
  summ
}

# ============================================================
# BATCH RUN
# ============================================================
all_results <- map_dfr(regimes, function(rg) {
  map_dfr(1:R_reps, function(rep_id) run_one_rep(rg, rep_id))
})

write_csv(all_results, "output/tables/all_mcmc_results.csv")

# ============================================================
# REGIME-LEVEL SUMMARY
# ============================================================
summary_all <- all_results %>%
  group_by(regime, parameter) %>%
  summarise(
    truth = first(truth),
    estimate = mean(mean, na.rm = TRUE),
    sd = mean(sd, na.rm = TRUE),
    lower = mean(lower, na.rm = TRUE),
    upper = mean(upper, na.rm = TRUE),
    bias = mean(bias, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE),
    width = mean(width, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_all, "output/tables/simulation_summary_all.csv")

# ============================================================
# TABLE BUILDERS
# ============================================================
make_regime_table <- function(df, caption, label, font_size = 8) {
  tab <- df %>%
    select(parameter, regime, bias, rmse, coverage) %>%
    pivot_wider(
      names_from = regime,
      values_from = c(bias, rmse, coverage),
      names_glue = "{regime}_{.value}"
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  kbl(tab, format = "latex", booktabs = TRUE, caption = caption, label = label) %>%
    kable_styling(latex_options = c("hold_position", "striped"), font_size = font_size)
}

epi_tab <- summary_all %>% filter(parameter %in% epi_params)
net_tab <- summary_all %>% filter(parameter %in% net_params)
obs_tab <- summary_all %>% filter(parameter %in% obs_params)

save_kable(
  make_regime_table(epi_tab, "Simulation performance for epidemic parameters across observation regimes.", "tab:epi_regime", 8),
  "output/tables/tab_epi_regime.tex"
)

save_kable(
  make_regime_table(net_tab, "Simulation performance for network parameters across observation regimes.", "tab:net_regime", 7),
  "output/tables/tab_net_regime.tex"
)

save_kable(
  make_regime_table(obs_tab, "Simulation performance for observation parameters across observation regimes.", "tab:obs_regime", 8),
  "output/tables/tab_obs_regime.tex"
)

# ============================================================
# FIGURES
# ============================================================
p_bias <- summary_all %>%
  mutate(parameter = fct_reorder(parameter, bias, .fun = mean)) %>%
  ggplot(aes(x = parameter, y = bias, fill = regime)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Bias", fill = "Regime")

ggsave("output/figures/bias_by_regime.png", p_bias, width = 12, height = 9, dpi = 300)

p_rmse <- summary_all %>%
  mutate(parameter = fct_reorder(parameter, rmse, .fun = mean)) %>%
  ggplot(aes(x = parameter, y = rmse, fill = regime)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "RMSE", fill = "Regime")

ggsave("output/figures/rmse_by_regime.png", p_rmse, width = 12, height = 9, dpi = 300)

p_cov <- summary_all %>%
  mutate(parameter = fct_reorder(parameter, coverage, .fun = mean)) %>%
  ggplot(aes(x = parameter, y = coverage, fill = regime)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Coverage", fill = "Regime")

ggsave("output/figures/coverage_by_regime.png", p_cov, width = 12, height = 9, dpi = 300)

p_truth <- summary_all %>%
  ggplot(aes(x = truth, y = estimate, color = regime)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ parameter, scales = "free", ncol = 4) +
  theme_classic(base_size = 11) +
  labs(x = "True value", y = "Posterior mean", color = "Regime") +
  theme(strip.background = element_blank())

ggsave("output/figures/truth_vs_estimate_by_regime.png", p_truth, width = 15, height = 10, dpi = 300)

# ============================================================
# SESSION INFO
# ============================================================
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
#######################################################################################################