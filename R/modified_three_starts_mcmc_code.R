#######################################################################################################
# ============================================================
# PARALLEL FULL MODEL-SPECIFIC SIMULATION + BAYESIAN MCMC
# SEIR on dynamic status-dependent contact network
# Three initial values: low, medium, high
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(future.apply)
  library(future)
})

set.seed(20260607)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("mcmc_convergence_plots", showWarnings = FALSE, recursive = TRUE)

plan(multisession, workers = max(1, parallel::detectCores() - 1))

# ============================================================
# SETTINGS
# ============================================================
R_reps <- 10
regimes <- c("High", "Moderate", "Sparse")
starts <- c("low", "medium", "high")

# Required order: epidemic parameters first, then observation parameters
epi_params <- c("beta", "xi", "kappa", "gamma")
obs_params <- c("p_E", "p_I", "s", "c")
main_params <- c(epi_params, obs_params)

net_pairs <- c("SS","SE","SI","SR","EE","EI","ER","II","IR","RR")
net_form_params <- paste0("eta_", net_pairs)
net_diss_params <- paste0("tau_", net_pairs)
net_params <- c(net_form_params, net_diss_params)

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
  st <- dat$stats
  
  epi_terms <- st$N_SE_int * log(theta["beta"]) - theta["beta"] * st$U_SI +
    st$N_SE_ext * log(theta["xi"]) - theta["xi"] * st$U_S +
    st$N_EI * log(theta["kappa"]) - theta["kappa"] * st$U_E +
    st$N_IR * log(theta["gamma"]) - theta["gamma"] * st$U_I
  
  obs_terms <- st$R_E * log(theta["p_E"]) + (st$M_E - st$R_E) * log(1 - theta["p_E"]) +
    st$R_I * log(theta["p_I"]) + (st$M_I - st$R_I) * log(1 - theta["p_I"]) +
    st$R_1 * log(theta["s"]) + (st$M_1 - st$R_1) * log(1 - theta["s"]) +
    st$R_0 * log(theta["c"]) + (st$M_0 - st$R_0) * log(1 - theta["c"])
  
  net_terms <- 0
  for (nm in net_pairs) {
    ns <- dat$net_stats[[nm]]
    net_terms <- net_terms +
      ns$N01 * log(theta[paste0("eta_", nm)]) - theta[paste0("eta_", nm)] * ns$V0 +
      ns$N10 * log(theta[paste0("tau_", nm)]) - theta[paste0("tau_", nm)] * ns$V1
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
# DATA GENERATION PLACEHOLDER
# ============================================================
simulate_dataset <- function(regime, rep_id) {
  truth_params <- c(
    beta = 0.30, xi = 0.05, kappa = 0.40, gamma = 0.25,
    p_E = 0.60, p_I = 0.85, s = 0.90, c = 0.95
  )
  truth_params <- c(
    truth_params,
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
  for (nm in net_pairs) net_stats[[nm]] <- list(N01 = 30, N10 = 28, V0 = 400, V1 = 380)
  
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
# ============================================================
update_epidemic_path <- function(Z, A, theta, dat) Z
update_network_path <- function(Z, A, theta, dat) A

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
# INITIAL VALUES
# ============================================================
initial_theta <- function(start) {
  if (start == "low") {
    theta <- c(
      beta = 0.01, xi = 0.01, kappa = 0.01, gamma = 0.01,
      p_E = 0.05, p_I = 0.05, s = 0.05, c = 0.05
    )
    eta0 <- 0.01
    tau0 <- 0.01
  } else if (start == "medium") {
    theta <- c(
      beta = 0.20, xi = 0.04, kappa = 0.30, gamma = 0.20,
      p_E = 0.50, p_I = 0.80, s = 0.85, c = 0.90
    )
    eta0 <- 0.05
    tau0 <- 0.05
  } else {
    theta <- c(
      beta = 0.80, xi = 0.30, kappa = 0.80, gamma = 0.70,
      p_E = 0.95, p_I = 0.95, s = 0.98, c = 0.99
    )
    eta0 <- 0.15
    tau0 <- 0.15
  }
  
  for (nm in net_pairs) {
    theta[paste0("eta_", nm)] <- eta0
    theta[paste0("tau_", nm)] <- tau0
  }
  theta
}

# ============================================================
# MCMC FITTING
# ============================================================
fit_bayesian_mcmc <- function(dat, start = "medium", n_iter = 700, burn = 200, thin = 2) {
  theta <- initial_theta(start)
  
  iter_draws <- matrix(NA_real_, nrow = n_iter, ncol = length(theta))
  colnames(iter_draws) <- names(theta)
  
  prop_sd <- setNames(rep(0.02, length(theta)), names(theta))
  prop_sd[c("beta", "xi", "kappa", "gamma")] <- c(0.03, 0.01, 0.03, 0.02)
  
  for (m in seq_len(n_iter)) {
    Z <- update_epidemic_path(dat$Z, dat$A, theta, dat)
    A <- update_network_path(Z, dat$A, theta, dat)
    
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
    
    iter_draws[m, ] <- theta
  }
  
  list(
    iter_draws = as_tibble(iter_draws) %>%
      mutate(iteration = seq_len(n_iter)),
    post_draws = as_tibble(iter_draws[seq(burn + 1, n_iter, by = thin), , drop = FALSE])
  )
}

# ============================================================
# SINGLE REPLICATE
# ============================================================
run_one_rep <- function(regime, rep_id, start) {
  dat <- simulate_dataset(regime, rep_id)
  fit <- fit_bayesian_mcmc(dat, start = start)
  
  list(
    iter_draws = fit$iter_draws %>% mutate(regime = regime, replicate = rep_id, start = start),
    post_draws = summarize_draws(fit$post_draws, truth_map) %>%
      mutate(regime = regime, replicate = rep_id, start = start)
  )
}

# ============================================================
# PARALLEL EXECUTION
# ============================================================
grid <- expand_grid(
  regime = regimes,
  replicate = seq_len(R_reps),
  start = starts
)

all_fits <- future_lapply(seq_len(nrow(grid)), function(i) {
  run_one_rep(grid$regime[i], grid$replicate[i], grid$start[i])
}, future.seed = TRUE)

all_iter_draws <- bind_rows(lapply(all_fits, `[[`, "iter_draws"))
all_results <- bind_rows(lapply(all_fits, `[[`, "post_draws"))

write_csv(all_results, "output/tables/all_mcmc_results.csv")
write_csv(all_iter_draws, "output/tables/all_mcmc_iter_draws.csv")

# ============================================================
# CONVERGENCE PLOTS
# ============================================================
truth_main <- truth_map %>%
  filter(parameter %in% main_params) %>%
  mutate(parameter = factor(parameter, levels = main_params))

iter_plot_data <- all_iter_draws %>%
  pivot_longer(
    cols = all_of(main_params),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    regime = factor(regime, levels = c("High", "Moderate", "Sparse")),
    start = factor(start, levels = c("low", "medium", "high")),
    parameter = factor(parameter, levels = main_params)
  )

# Parsed plotmath labels:
# beta, xi, kappa and gamma are displayed as Greek symbols;
# p_E and p_I are displayed as p with E/I subscripts.
param_labels <- c(
  "beta"  = "beta",
  "xi"    = "xi",
  "kappa" = "kappa",
  "gamma" = "gamma",
  "p_E"   = "p[E]",
  "p_I"   = "p[I]",
  "s"     = "s",
  "c"     = "c"
)

# ------------------------------------------------------------
# 1. One faceted convergence plot per observation regime
# ------------------------------------------------------------
for (r in levels(iter_plot_data$regime)) {
  p <- iter_plot_data %>%
    filter(regime == r) %>%
    ggplot(aes(x = iteration, y = value, color = start, group = interaction(start, replicate))) +
    geom_line(alpha = 0.45, linewidth = 0.35) +
    geom_hline(
      data = truth_main,
      aes(yintercept = truth),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    geom_vline(xintercept = 200, linetype = "dotted", color = "gray40") +
    facet_wrap(
      ~ parameter,
      scales = "free_y",
      ncol = 2,
      labeller = as_labeller(param_labels, default = label_parsed)
    ) +
    labs(
      title = paste("MCMC convergence from three initial values:", r, "observation regime"),
      subtitle = "Low, medium, and high initial guesses are shown; dotted line marks burn-in end",
      x = "Iteration",
      y = "Parameter value",
      color = "Initial value"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 17,
        color = "#333333"
      ),
      plot.subtitle = element_text(
        size = 13,
        color = "#333333"
      ),
      strip.background = element_rect(
        fill = "grey75",
        color = "grey75"
      ),
      strip.text = element_text(
        color = "grey15",
        face = "bold",
        size = 14
      ),
      axis.title = element_text(
        size = 12,
        face = "bold",
        color = "#333333"
      ),
      axis.text = element_text(
        size = 13,
        color = "#333333"
      )
    )
  
  ggsave(
    filename = paste0("mcmc_convergence_plots/convergence_3starts_", r, ".png"),
    plot = p,
    width = 10,
    height = 10,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 2. One convergence plot per parameter across observation regimes
# ------------------------------------------------------------
for (par in main_params) {
  plot_title <- bquote(
    "MCMC convergence for" ~ .(parse(text = param_labels[[par]])[[1]]) ~
      "from three initial values"
  )
  
  p <- iter_plot_data %>%
    filter(parameter == par) %>%
    ggplot(aes(x = iteration, y = value, color = start, group = interaction(start, replicate))) +
    geom_line(alpha = 0.5, linewidth = 0.4) +
    geom_hline(
      data = truth_main %>% filter(parameter == par),
      aes(yintercept = truth),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    geom_vline(xintercept = 200, linetype = "dotted", color = "gray40") +
    facet_wrap(~ regime, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = pretty(iter_plot_data$iteration)) +
    labs(
      title = plot_title,
      subtitle = "Low, medium, and high initial guesses are shown; dotted line marks burn-in end",
      x = "Iteration",
      y = "Parameter value",
      color = "Initial value"
    ) +
    theme_classic(base_size = 12) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 17,
        color = "#333333"
      ),
      plot.subtitle = element_text(
        size = 13,
        color = "#333333"
      ),
      strip.background = element_rect(
        fill = "grey75",
        color = "grey75"
      ),
      strip.text = element_text(
        color = "grey15",
        face = "bold",
        size = 14
      ),
      axis.title = element_text(
        size = 12,
        face = "bold",
        color = "#333333"
      ),
      axis.text = element_text(
        size = 13,
        color = "#333333"
      )
    )
  
  ggsave(
    filename = paste0("mcmc_convergence_plots/convergence_3starts_", par, ".png"),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
}

message("Finished. Files written to output/tables and mcmc_convergence_plots.")
# ------------------------------------------------------------------------------
#######################################################################################################