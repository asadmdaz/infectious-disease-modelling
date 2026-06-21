#######################################################################################################
# ============================================================
# PARALLEL FULL MODEL-SPECIFIC SIMULATION + BAYESIAN MCMC
# SEIR on dynamic status-dependent contact network
# Uses future.apply for replicate-level parallelization
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
  library(future.apply)
  library(future)
})

set.seed(20260607)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/sim_results", showWarnings = FALSE, recursive = TRUE)
dir.create("output/mcmc", showWarnings = FALSE, recursive = TRUE)

# Choose backend:
# multisession works on Windows/macOS/Linux
# multicore works on Linux/macOS, not Windows
plan(multisession, workers = max(1, parallel::detectCores() - 1))

# ============================================================
# SETTINGS
# ============================================================
R_reps <- 1
regimes <- c("High", "Moderate", "Sparse")

epi_params <- c("beta", "xi", "kappa", "gamma")
obs_params <- c("p_E", "p_I", "s", "c")
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
    beta = 0.30, xi = 0.05, kappa = 0.40, gamma = 0.10,
    p_E = 0.60, p_I = 0.85, s = 0.90, c = 0.95
  )
  truth_params <- c(truth_params,
                    setNames(rep(0.08, length(net_pairs)), paste0("eta_", net_pairs)),
                    setNames(rep(0.03, length(net_pairs)), paste0("tau_", net_pairs)))
  
  stats <- list(
    N_SE_int = 20, N_SE_ext = 8, N_EI = 18, N_IR = 15,
    U_SI = 65, U_S = 120, U_E = 45, U_I = 55,
    M_E = 100, R_E = 60, M_I = 100, R_I = 84,
    M_1 = 800, R_1 = 720, M_0 = 1200, R_0 = 1140
  )
  net_stats <- setNames(vector("list", length(net_pairs)), net_pairs)
  for (nm in net_pairs) net_stats[[nm]] <- list(N01 = 30, N10 = 28, V0 = 400, V1 = 380)
  
  list(regime = regime, replicate = rep_id, truth = truth_params, stats = stats, net_stats = net_stats, Z = NULL, A = NULL)
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
# MCMC FITTING
# ============================================================
#beta = 0.01, xi = 0.01, kappa = 0.01, gamma = 0.01,
#p_E = 0.01, p_I = 0.01, s = 0.01, c = 0.01

fit_bayesian_mcmc <- function(dat, n_iter = 200, burn = 50, thin = 2) {
  theta <- c(
    beta = 0.01, xi = 0.01, kappa = 0.01, gamma = 0.01,
    p_E = 0.01, p_I = 0.01, s = 0.01, c = 0.01
  )
  for (nm in net_pairs) {
    theta[paste0("eta_", nm)] <- 0.05
    theta[paste0("tau_", nm)] <- 0.05
  }
  
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
run_one_rep <- function(regime, rep_id) {
  dat <- simulate_dataset(regime, rep_id)
  fit <- fit_bayesian_mcmc(dat)
  
  list(
    iter_draws = fit$iter_draws %>% mutate(regime = regime, replicate = rep_id),
    post_draws = summarize_draws(fit$post_draws, truth_map) %>%
      mutate(regime = regime, replicate = rep_id)
  )
}

# ============================================================
# PARALLEL EXECUTION
# ============================================================
# FIX: Create the grid dataframe combining regimes and replicates
grid <- expand_grid(
  regime = regimes,
  replicate = seq_len(R_reps)
)

all_fits <- future_lapply(seq_len(nrow(grid)), function(i) {
  rg <- grid$regime[i]
  rep_id <- grid$replicate[i]
  run_one_rep(rg, rep_id)
}, future.seed = TRUE)

all_iter_draws <- bind_rows(lapply(all_fits, `[[`, "iter_draws"))
all_results <- bind_rows(lapply(all_fits, `[[`, "post_draws"))

write_csv(all_results, "output/tables/all_mcmc_results.csv")
write_csv(all_iter_draws, "output/tables/all_mcmc_iter_draws.csv")


#xxxxxxxxxxxxxxxxxxxxx
R_reps <- 1
regimes <- c("High", "Moderate", "Sparse")

# FIX: Define burn globally so your plotting code can read it
burn <- 50  

epi_params <- c("beta", "xi", "kappa", "gamma")
obs_params <- c("p_E", "p_I", "s", "c")
# Define the parameters you actually want to plot
main_params <- c("beta", "xi", "kappa", "gamma", "p_E", "p_I", "s", "c")

# Define the missing labels object
param_labels <- c(
  "beta"  = "Beta (Transmission rate)",
  "xi"    = "Xi (External infection)",
  "kappa" = "Kappa (E to I rate)",
  "gamma" = "Gamma (Recovery rate)",
  "p_E"   = "p_E (E detection prob)",
  "p_I"   = "p_I (I detection prob)",
  "s"     = "s (True positive rate)",
  "c"     = "c (True negative rate)"
)

iter_plot_data <- all_iter_draws %>%
  # Crucial Step: Transform parameters from column names into row values
  pivot_longer(
    cols = -c(iteration, regime, replicate), 
    names_to = "parameter", 
    values_to = "value"
  ) %>%
  filter(parameter %in% main_params) %>%
  mutate(
    regime = factor(regime, levels = c("High", "Moderate", "Sparse")),
    parameter = factor(parameter, levels = main_params)
  )

# ============================================================
# PARALLEL FULL MODEL-SPECIFIC SIMULATION + BAYESIAN MCMC
# SEIR on dynamic status-dependent contact network
# Uses future.apply for replicate-level parallelization
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
  library(future.apply)
  library(future)
})

set.seed(20260607)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/sim_results", showWarnings = FALSE, recursive = TRUE)
dir.create("output/mcmc", showWarnings = FALSE, recursive = TRUE)

# Choose backend:
# multisession works on Windows/macOS/Linux
# multicore works on Linux/macOS, not Windows
plan(multisession, workers = max(1, parallel::detectCores() - 1))

# ============================================================
# SETTINGS
# ============================================================
R_reps <- 1
regimes <- c("High", "Moderate", "Sparse")

epi_params <- c("beta", "xi", "kappa", "gamma")
obs_params <- c("p_E", "p_I", "s", "c")
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
    beta = 0.30, xi = 0.05, kappa = 0.40, gamma = 0.10,
    p_E = 0.60, p_I = 0.85, s = 0.90, c = 0.95
  )
  truth_params <- c(truth_params,
                    setNames(rep(0.08, length(net_pairs)), paste0("eta_", net_pairs)),
                    setNames(rep(0.03, length(net_pairs)), paste0("tau_", net_pairs)))
  
  stats <- list(
    N_SE_int = 20, N_SE_ext = 8, N_EI = 18, N_IR = 15,
    U_SI = 65, U_S = 120, U_E = 45, U_I = 55,
    M_E = 100, R_E = 60, M_I = 100, R_I = 84,
    M_1 = 800, R_1 = 720, M_0 = 1200, R_0 = 1140
  )
  net_stats <- setNames(vector("list", length(net_pairs)), net_pairs)
  for (nm in net_pairs) net_stats[[nm]] <- list(N01 = 30, N10 = 28, V0 = 400, V1 = 380)
  
  list(regime = regime, replicate = rep_id, truth = truth_params, stats = stats, net_stats = net_stats, Z = NULL, A = NULL)
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
# MCMC FITTING
# ============================================================
#beta = 0.01, xi = 0.01, kappa = 0.01, gamma = 0.01,
#p_E = 0.01, p_I = 0.01, s = 0.01, c = 0.01

fit_bayesian_mcmc <- function(dat, n_iter = 200, burn = 50, thin = 2) {
  theta <- c(
    beta = 0.01, xi = 0.01, kappa = 0.01, gamma = 0.01,
    p_E = 0.01, p_I = 0.01, s = 0.01, c = 0.01
  )
  for (nm in net_pairs) {
    theta[paste0("eta_", nm)] <- 0.05
    theta[paste0("tau_", nm)] <- 0.05
  }
  
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
run_one_rep <- function(regime, rep_id) {
  dat <- simulate_dataset(regime, rep_id)
  fit <- fit_bayesian_mcmc(dat)
  
  list(
    iter_draws = fit$iter_draws %>% mutate(regime = regime, replicate = rep_id),
    post_draws = summarize_draws(fit$post_draws, truth_map) %>%
      mutate(regime = regime, replicate = rep_id)
  )
}

# ============================================================
# PARALLEL EXECUTION
# ============================================================
# FIX: Create the grid dataframe combining regimes and replicates
grid <- expand_grid(
  regime = regimes,
  replicate = seq_len(R_reps)
)

all_fits <- future_lapply(seq_len(nrow(grid)), function(i) {
  rg <- grid$regime[i]
  rep_id <- grid$replicate[i]
  run_one_rep(rg, rep_id)
}, future.seed = TRUE)

all_iter_draws <- bind_rows(lapply(all_fits, `[[`, "iter_draws"))
all_results <- bind_rows(lapply(all_fits, `[[`, "post_draws"))

write_csv(all_results, "output/tables/all_mcmc_results.csv")
write_csv(all_iter_draws, "output/tables/all_mcmc_iter_draws.csv")



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Convergence plots

library(dplyr)
library(tidyr)
library(ggplot2)

R_reps <- 1
regimes <- c("High", "Moderate", "Sparse")
burn <- 50  

# 1. Define the explicit order
main_params <- c("beta", "kappa", "gamma", "xi", "p_E", "p_I", "s", "c")

param_labels <- c(
  "beta"  = "beta ~ '(Transmission rate)'",
  "kappa" = "kappa ~ '(E to I rate)'",
  "gamma" = "gamma ~ '(Recovery rate)'",
  "xi"    = "xi ~ '(External infection)'",
  "p_E"   = "p[E] ~ '(E detection prob)'",
  "p_I"   = "p[I] ~ '(I detection prob)'",
  "s"     = "s ~ '(True positive rate)'",
  "c"     = "c ~ '(True negative rate)'"
)

# 2. Prepare the plot data with strict factoring
iter_plot_data <- all_iter_draws %>%
  pivot_longer(
    cols = -c(iteration, regime, replicate), 
    names_to = "parameter", 
    values_to = "value"
  ) %>%
  filter(parameter %in% main_params) %>%
  mutate(
    regime = factor(regime, levels = c("High", "Moderate", "Sparse")),
    parameter = factor(parameter, levels = main_params) # Enforced here
  )

# 3. CRITICAL FIX: Explicitly factorize the truth_map dataframe as well
cleaned_truth_map <- truth_map %>% 
  filter(parameter %in% main_params) %>%
  mutate(parameter = factor(parameter, levels = main_params)) # Enforced here too

# 4. Run the plotting loop
for (r in levels(iter_plot_data$regime)) {
  p <- iter_plot_data %>%
    filter(regime == r) %>%
    ggplot(aes(x = iteration, y = value, group = replicate)) +
    geom_line(alpha = 0.5, linewidth = 0.5, col="#41B7C4") + 
    
    # Updated to use the explicitly factored truth map data
    geom_hline(data = cleaned_truth_map,
               aes(yintercept = truth), linetype = "dashed", color = "red", linewidth = 0.5) +
    
    facet_wrap(~ parameter, scales = "free_y", ncol = 2, 
               labeller = as_labeller(param_labels, default = label_parsed)) +
    
    geom_vline(xintercept = burn, linetype = "dotted", color = "gray30") +
    labs(
      title = paste("MCMC convergence from initial values:", r, "observation regime"),
      subtitle = "All iterations are shown; dotted line marks the end of burn-in.",
      x = "Iteration",
      y = "Parameter value"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "gray90", color = "gray30", linewidth = 0.5),
      strip.text = element_text(face = "bold", color = "black")
    )
  
  ggsave(
    filename = paste0("mcmc_convergence_plots/convergence_iter_", r, ".png"),
    plot = p,
    width = 9,
    height = 10,
    dpi = 300
  )
}
#xxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx






# ============================================================
# REGIME SUMMARY
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

iter_plot_data <- all_iter_draws %>%
  pivot_longer(
    cols = all_of(main_params),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    regime = factor(regime, levels = c("High", "Moderate", "Sparse")),
    parameter = factor(parameter, levels = main_params)
  )

for (r in levels(iter_plot_data$regime)) {
  p <- iter_plot_data %>%
    filter(regime == r) %>%
    # If you want to drop burn-in iterations as discussed previously, add:
    # filter(iteration > burn) %>% 
    ggplot(aes(x = iteration, y = value, group = interaction(regime, replicate))) +
    geom_line(alpha = 0.12, linewidth = 0.25) +
    geom_vline(xintercept = 5, linetype = "dotted", color = "gray40") +
    geom_hline(
      data = truth_map %>% filter(parameter %in% main_params),
      aes(yintercept = truth),
      # FIX: Removed 'inherit.aes = FALSE' to stop the warning
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2, labeller = labeller(parameter = param_labels)) +
    labs(
      title = paste("MCMC convergence from initial values:", r, "observation regime"),
      subtitle = "All iterations are shown; dotted line marks the end of burn-in.",
      x = "Iteration",
      y = "Parameter value"
    ) +
    theme_classic(base_size = 12)
  
  ggsave(
    filename = paste0("mcmc_convergence_plots/convergence_iter_", r, ".png"),
    plot = p,
    width = 9,
    height = 10,
    dpi = 300
  )
}


# ============================================================
# TABLES
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

save_kable(
  make_regime_table(summary_all %>% filter(parameter %in% epi_params),
                    "Simulation performance for epidemic parameters across observation regimes.",
                    "tab:epi_regime", 8),
  "output/tables/tab_epi_regime.tex"
)

save_kable(
  make_regime_table(summary_all %>% filter(parameter %in% net_params),
                    "Simulation performance for network parameters across observation regimes.",
                    "tab:net_regime", 7),
  "output/tables/tab_net_regime.tex"
)

save_kable(
  make_regime_table(summary_all %>% filter(parameter %in% obs_params),
                    "Simulation performance for observation parameters across observation regimes.",
                    "tab:obs_regime", 8),
  "output/tables/tab_obs_regime.tex"
)


# Read the saved summary file
library(dplyr)
library(readr)
library(knitr)
library(kableExtra)

library(dplyr)
library(readr)
library(knitr)
library(kableExtra)

summary_all <- read_csv("output/tables/simulation_summary_all.csv", show_col_types = FALSE)

main_params <- c("beta", "xi", "kappa", "gamma", "p_E", "p_I", "s", "c")

posterior_table <- summary_all %>%
  filter(parameter %in% main_params) %>%
  mutate(
    Parameter = case_when(
      parameter == "beta" ~ "$\\beta$",
      parameter == "xi" ~ "$\\xi$",
      parameter == "kappa" ~ "$\\kappa$",
      parameter == "gamma" ~ "$\\gamma$",
      parameter == "p_E" ~ "$p_E$",
      parameter == "p_I" ~ "$p_I$",
      parameter == "s" ~ "$s$",
      parameter == "c" ~ "$c$",
      TRUE ~ parameter
    ),
    Truth = round(truth, 4),
    Mean = round(estimate, 4),
    SD = round(sd, 4),
    MSE = round((estimate - truth)^2, 4),
    Lower = round(lower, 4),
    Upper = round(upper, 4),
    Bias = round(bias, 4),
    `Abs Error` = round(abs(estimate - truth), 4),
    Coverage = ifelse(coverage >= 0.95, "Yes", "No")
  ) %>%
  select(regime, Parameter, Truth, Mean, SD, MSE, Lower, Upper, Truth, Bias, `Abs Error`, Coverage) %>%
  arrange(factor(regime, levels = c("High", "Moderate", "Sparse")), Parameter)

kable(
  posterior_table,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Posterior summaries from Bayesian MCMC",
  label = "tab:posterior_summaries"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    font_size = 8
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

writeLines(capture.output(sessionInfo()), "output/session_info.txt")

# ============================================================
# Convergence-style plots from replicate-level MCMC summaries
# For regimes: High, Moderate, Sparse
# Parameters: beta, xi, kappa, gamma, p_E, p_I, s, c
# ============================================================

################################################################################
## Better versions with greek symbol and thinner lines
# ============================================================
# Convergence-style plots from replicate-level MCMC summaries
# For regimes: High, Moderate, Sparse
# Parameters: beta, xi, kappa, gamma, p_E, p_I, s, c
# ============================================================

# ============================================================
# MCMC stability / convergence-style plots
# Mathematical parameter labels for JASA-style figures
# Parameters: beta, xi, kappa, gamma, p_E, p_I, s, c
# Regimes: High, Moderate, Sparse
# ============================================================

# ============================================================
# MCMC stability / convergence-style plots
# JASA-style figures with large mathematical / Greek symbols
# Parameters: beta, xi, kappa, gamma, p_E, p_I, s, c
# Regimes: High, Moderate, Sparse
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)
library(purrr)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

mcmc_results <- read_csv("output/tables/all_mcmc_results.csv")

# ------------------------------------------------------------
# Parameters and output folder
# ------------------------------------------------------------

main_params <- c("beta", "xi", "kappa", "gamma", "p_E", "p_I", "s", "c")

dir.create("mcmc_convergence_plots", showWarnings = FALSE)

plot_data <- mcmc_results %>%
  filter(parameter %in% main_params) %>%
  mutate(
    regime = factor(regime, levels = c("High", "Moderate", "Sparse")),
    parameter = factor(parameter, levels = main_params)
  )

# ------------------------------------------------------------
# Mathematical labels
# ------------------------------------------------------------

# Facet labels
param_labels <- c(
  beta  = "beta",
  xi    = "xi",
  kappa = "kappa",
  gamma = "gamma",
  p_E   = "p[E]",
  p_I   = "p[I]",
  s     = "s",
  c     = "c"
)

# Plot title labels
param_title_labels <- list(
  beta  = expression(beta),
  xi    = expression(xi),
  kappa = expression(kappa),
  gamma = expression(gamma),
  p_E   = expression(p[E]),
  p_I   = expression(p[I]),
  s     = expression(s),
  c     = expression(c)
)

# ------------------------------------------------------------
# Colour palettes
# ------------------------------------------------------------

#regime_pal <- c(
#  High = "#1B9E77",
#  Moderate = "#D95F02",
#  Sparse = "#7570B3"
#)

#param_pal <- c(
#  beta  = "#66C2A5",
#  xi    = "#FC8D62",
#  kappa = "#8DA0CB",
#  gamma = "#E78AC3",
#  p_E   = "#A6D854",
#  p_I   = "#FFD92F",
#  s     = "#E5C494",
#  c     = "#B3B3B3"
#)


# ------------------------------------------------------------
# 1. Faceted stability plots for each regime
#    Large mathematical symbols in facet strips
# ------------------------------------------------------------

for (r in levels(plot_data$regime)) {
  
  p <- plot_data %>%
    filter(regime == r) %>%
    ggplot(aes(
      x = replicate,
      y = mean,
      color = parameter,
      group = parameter
    )) +
    geom_line(linewidth = 0.45, alpha = 0.85, col="#009593") +
    geom_point(size = 1.4, alpha = 0.85, col="#009593") +
    geom_hline(
      aes(yintercept = truth),
      linetype = "dashed",
      color = "black",
      linewidth = 0.55
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free_y",
      ncol = 2,
      labeller = as_labeller(param_labels, label_parsed)
    ) +
    scale_color_manual(values = param_pal, guide = "none") +
    labs(
      title = paste("MCMC stability across replicates:", r, "observation regime"),
      subtitle = "Coloured traces show posterior means; dashed black line shows the true value",
      x = "Simulation replicate",
      y = "Posterior mean"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        color = "#CF597E"
      ),
      plot.subtitle = element_text(
        size = 13,
        color = "#E6886A"
      ),
      strip.background = element_rect(
        fill = "#023FA5",
        color = "#41B7C4"
      ),
      strip.text = element_text(
        color = "white",
        face = "bold",
        size = 18
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "#333333"
      ),
      axis.text = element_text(
        size = 13,
        color = "#333333"
      )
    )
  
  ggsave(
    filename = paste0("mcmc_convergence_plots/convergence_", r, ".png"),
    plot = p,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 2. One figure per parameter, comparing regimes
#    Very large mathematical symbol as plot title
# ------------------------------------------------------------

for (par in main_params) {
  
  p <- plot_data %>%
    filter(parameter == par) %>%
    ggplot(aes(
      x = replicate,
      y = mean,
      color = regime,
      group = regime
    )) +
    geom_line(linewidth = 0.45, alpha = 0.85, col="#009593") +
    geom_point(size = 1.4, alpha = 0.85, col="#009593") +
    geom_hline(
      aes(yintercept = truth),
      linetype = "dashed",
      color = "black",
      linewidth = 0.55
    ) +
    facet_wrap(
      ~ regime,
      scales = "free_y",
      ncol = 1
    ) +
    scale_color_manual(values = regime_pal) +
    labs(
      subtitle = "Each colour corresponds to an observation regime; dashed black line shows the true value",
      x = "Simulation replicate",
      y = "Posterior mean",
      colour = "Regime"
    ) +
    ggtitle(param_title_labels[[par]]) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        size = 28,
        color = "#CF597E",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 13,
        color = "#E6886A",
        hjust = 0.5
      ),
      strip.background = element_rect(
        fill = "#023FA5",
        color = "#41B7C4"
      ),
      strip.text = element_text(
        face = "bold",
        color = "white",
        size = 16
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "#333333"
      ),
      axis.text = element_text(
        size = 13,
        color = "#333333"
      )
    )
  
  ggsave(
    filename = paste0("mcmc_convergence_plots/convergence_", par, ".png"),
    plot = p,
    width = 8.5,
    height = 7,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 3. Posterior interval stability plots for each regime
#    Large mathematical symbols in facet strips
# ------------------------------------------------------------

for (r in levels(plot_data$regime)) {
  
  p <- plot_data %>%
    filter(regime == r) %>%
    ggplot(aes(
      x = replicate,
      y = mean,
      color = parameter,
      fill = parameter,
      group = parameter
    )) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.36,
      color = NA
    ) +
    geom_line(linewidth = 0.45, alpha = 0.85, col="#009593") +
    geom_point(size = 1.4, alpha = 0.85, col="#009593") +
    geom_hline(
      aes(yintercept = truth),
      linetype = "dashed",
      color = "black",
      linewidth = 0.55
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free_y",
      ncol = 2,
      labeller = as_labeller(param_labels, label_parsed)
    ) +
    scale_color_manual(values = param_pal, guide = "none") +
    scale_fill_manual(values = param_pal, guide = "none") +
    labs(
      title = paste("Posterior interval stability:", r, "observation regime"),
      subtitle = "Shaded bands show posterior intervals; dashed black line shows the true value",
      x = "Simulation replicate",
      y = "Posterior mean with interval"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        color = "#CF597E"
      ),
      plot.subtitle = element_text(
        size = 13,
        color = "#E6886A"
      ),
      strip.background = element_rect(
        fill = "#023FA5",
        color = "#41B7C4"
      ),
      strip.text = element_text(
        color = "white",
        face = "bold",
        size = 18
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "#333333"
      ),
      axis.text = element_text(
        size = 13,
        color = "#333333"
      )
    )
  
  ggsave(
    filename = paste0("mcmc_convergence_plots/interval_stability_", r, ".png"),
    plot = p,
    width = 9,
    height = 10,
    dpi = 300
  )
}
# ------------------------------------------------------------------------------
#######################################################################################################