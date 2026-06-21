#################################################################################################################
# ============================================================
# Simulation study: epidemic-network model
# ============================================================

# Packages ----------------------------------------------------
required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "tibble",
  "readr", "ggplot2", "stringr", "knitr"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_pkgs, installed)
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(ggplot2)
library(stringr)
library(knitr)

# Options -----------------------------------------------------
set.seed(12345)
options(stringsAsFactors = FALSE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

##############################################################
# ============================================================
# Joint SEIR + dynamic network simulation generator
# Partially observed epidemic processes on dynamic networks
# ============================================================

# -------------------------------
# 1. Setup
# -------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
})

set.seed(20260201)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 2. Model parameters
# -------------------------------
params <- list(
  N = 100,
  T = 30,
  beta = 0.30,   # internal infection
  xi = 0.05,     # external infection
  kappa = 0.40,  # E -> I
  gamma = 0.25   # I -> R
)

# Network parameters by dyad disease-state pair
state_levels <- c("S", "E", "I", "R")

make_pair_name <- function(a, b) {
  paste(sort(c(a, b)), collapse = "")
}

pair_names <- unique(
  sapply(state_levels, function(a) sapply(state_levels, function(b) make_pair_name(a, b)))
)

set_names <- function(x, nm) { names(x) <- nm; x }

eta <- setNames(rep(0.08, length(pair_names)), pair_names)
omega <- setNames(rep(0.03, length(pair_names)), pair_names)

# Example heterogeneity: lower formation among SI dyads
eta["SI"] <- 0.05
eta["II"] <- 0.04
omega["SI"] <- 0.06
omega["II"] <- 0.07

obs <- list(
  p_E = 0.60,
  p_I = 0.85,
  s = 0.90,
  c = 0.95
)

# Observation schedule
symptom_times <- seq(2, params$T, by = 2)
contact_times  <- seq(1, params$T, by = 3)

# -------------------------------
# 3. Initialization
# -------------------------------
N <- params$N
Tmax <- params$T

nodes <- 1:N
dyads <- combn(nodes, 2, simplify = FALSE)

init_states <- rep("S", N)
init_states[sample(nodes, 3)] <- "I"  # seed 3 infectious
init_states[sample(which(init_states == "S"), 5)] <- "E"

init_network <- matrix(0L, N, N)
diag(init_network) <- 0L

# -------------------------------
# 4. Utility functions
# -------------------------------
pair_key <- function(xi, xj) make_pair_name(xi, xj)

infectious_contacts <- function(i, states, A) {
  sum(A[i, ] == 1L & states == "I")
}

internal_infection_rate <- function(i, states, A, beta) {
  beta * infectious_contacts(i, states, A)
}

network_rate <- function(i, j, states, A, eta, omega) {
  p <- pair_key(states[i], states[j])
  if (A[i, j] == 0L) eta[p] else omega[p]
}

total_hazard <- function(states, A, params, eta, omega) {
  beta <- params$beta
  xi <- params$xi
  kappa <- params$kappa
  gamma <- params$gamma
  
  haz <- 0
  
  # Epidemic hazards
  for (i in 1:length(states)) {
    if (states[i] == "S") haz <- haz + internal_infection_rate(i, states, A, beta) + xi
    if (states[i] == "E") haz <- haz + kappa
    if (states[i] == "I") haz <- haz + gamma
  }
  
  # Network hazards
  for (ij in dyads) {
    i <- ij[1]; j <- ij[2]
    haz <- haz + network_rate(i, j, states, A, eta, omega)
  }
  
  haz
}

simulate_next_event <- function(t, states, A, params, eta, omega) {
  N <- length(states)
  beta <- params$beta
  xi <- params$xi
  kappa <- params$kappa
  gamma <- params$gamma
  
  event_list <- list()
  rates <- numeric(0)
  
  # Epidemic events
  for (i in 1:N) {
    if (states[i] == "S") {
      r_int <- internal_infection_rate(i, states, A, beta)
      r_ext <- xi
      if (r_int > 0) {
        event_list[[length(event_list) + 1]] <- list(type = "SE_int", i = i)
        rates <- c(rates, r_int)
      }
      if (r_ext > 0) {
        event_list[[length(event_list) + 1]] <- list(type = "SE_ext", i = i)
        rates <- c(rates, r_ext)
      }
    } else if (states[i] == "E") {
      event_list[[length(event_list) + 1]] <- list(type = "EI", i = i)
      rates <- c(rates, kappa)
    } else if (states[i] == "I") {
      event_list[[length(event_list) + 1]] <- list(type = "IR", i = i)
      rates <- c(rates, gamma)
    }
  }
  
  # Network events
  for (ij in dyads) {
    i <- ij[1]; j <- ij[2]
    p <- pair_key(states[i], states[j])
    if (A[i, j] == 0L) {
      event_list[[length(event_list) + 1]] <- list(type = "A01", i = i, j = j, pair = p)
      rates <- c(rates, eta[p])
    } else {
      event_list[[length(event_list) + 1]] <- list(type = "A10", i = i, j = j, pair = p)
      rates <- c(rates, omega[p])
    }
  }
  
  total <- sum(rates)
  if (total <= 0) return(NULL)
  
  dt <- rexp(1, rate = total)
  which_event <- sample(seq_along(rates), size = 1, prob = rates)
  list(t_next = t + dt, event = event_list[[which_event]], total_hazard = total)
}

apply_event <- function(states, A, event) {
  if (event$type == "SE_int" || event$type == "SE_ext") {
    states[event$i] <- "E"
  } else if (event$type == "EI") {
    states[event$i] <- "I"
  } else if (event$type == "IR") {
    states[event$i] <- "R"
  } else if (event$type == "A01") {
    A[event$i, event$j] <- 1L
    A[event$j, event$i] <- 1L
  } else if (event$type == "A10") {
    A[event$i, event$j] <- 0L
    A[event$j, event$i] <- 0L
  }
  list(states = states, A = A)
}

record_state <- function(t, states, A, latent_epi, latent_net) {
  epi_row <- tibble(
    time = t,
    id = 1:length(states),
    state = states
  )
  
  net_row <- expand.grid(i = 1:(nrow(A) - 1), j = 2:nrow(A)) %>%
    as_tibble() %>%
    filter(i < j) %>%
    mutate(
      time = t,
      Aij = map2_int(i, j, ~ A[.x, .y])
    )
  
  list(epi = bind_rows(latent_epi, epi_row),
       net = bind_rows(latent_net, net_row))
}

# -------------------------------
# 5. Simulate latent trajectory
# -------------------------------
simulate_latent <- function(params, eta, omega, init_states, init_network) {
  t <- 0
  states <- init_states
  A <- init_network
  
  latent_epi <- tibble(time = numeric(), id = integer(), state = character())
  latent_net <- tibble(time = numeric(), i = integer(), j = integer(), Aij = integer())
  event_log <- tibble(time = numeric(), type = character(), i = integer(), j = integer(), pair = character())
  
  latent <- record_state(t, states, A, latent_epi, latent_net)
  latent_epi <- latent$epi
  latent_net <- latent$net
  
  while (t < params$T) {
    nxt <- simulate_next_event(t, states, A, params, eta, omega)
    if (is.null(nxt)) break
    if (nxt$t_next > params$T) break
    
    t <- nxt$t_next
    upd <- apply_event(states, A, nxt$event)
    states <- upd$states
    A <- upd$A
    
    latent <- record_state(t, states, A, latent_epi, latent_net)
    latent_epi <- latent$epi
    latent_net <- latent$net
    
    event_log <- bind_rows(event_log, tibble(
      time = t,
      type = nxt$event$type,
      i = ifelse(is.null(nxt$event$i), NA_integer_, nxt$event$i),
      j = ifelse(is.null(nxt$event$j), NA_integer_, nxt$event$j),
      pair = ifelse(is.null(nxt$event$pair), NA_character_, nxt$event$pair)
    ))
  }
  
  list(epi = latent_epi, net = latent_net, events = event_log)
}

latent <- simulate_latent(params, eta, omega, init_states, init_network)

# -------------------------------
# 6. Interpolate latent states at observation times
# -------------------------------
last_state_before <- function(times, values, t0) {
  idx <- max(which(times <= t0))
  values[idx]
}

get_states_at <- function(latent_epi, t0) {
  out <- latent_epi %>%
    group_by(id) %>%
    arrange(time) %>%
    summarise(state = last_state_before(time, state, t0), .groups = "drop")
  out$state
}

get_network_at <- function(latent_net, t0, N) {
  A <- matrix(0L, N, N)
  diag(A) <- 0L
  
  sub <- latent_net %>% filter(time <= t0)
  if (nrow(sub) == 0) return(A)
  
  latest <- sub %>%
    group_by(i, j) %>%
    arrange(time) %>%
    summarise(Aij = dplyr::last(Aij), .groups = "drop")
  
  for (k in seq_len(nrow(latest))) {
    i <- latest$i[k]; j <- latest$j[k]
    A[i, j] <- latest$Aij[k]
    A[j, i] <- latest$Aij[k]
  }
  A
}

# -------------------------------
# 7. Generate observations
# -------------------------------
simulate_observations <- function(latent_epi, latent_net, params, obs, symptom_times, contact_times, N) {
  symp_obs <- map_dfr(symptom_times, function(t0) {
    states_t <- get_states_at(latent_epi, t0)
    tibble(
      time = t0,
      id = 1:N,
      state = states_t,
      Y = map_int(states_t, function(st) {
        if (st == "E") rbinom(1, 1, obs$p_E)
        else if (st == "I") rbinom(1, 1, obs$p_I)
        else rbinom(1, 1, 0.05)
      })
    )
  })
  
  cont_obs <- map_dfr(contact_times, function(t0) {
    A_t <- get_network_at(latent_net, t0, N)
    pairs <- combn(1:N, 2, simplify = FALSE)
    tibble(
      time = t0,
      i = map_int(pairs, 1),
      j = map_int(pairs, 2),
      A = map_int(pairs, function(v) A_t[v[1], v[2]]),
      B = map_int(seq_along(pairs), function(k) {
        a <- A_t[pairs[[k]][1], pairs[[k]][2]]
        if (a == 1L) rbinom(1, 1, obs$s) else rbinom(1, 1, 1 - obs$c)
      })
    )
  })
  
  list(symptoms = symp_obs, contacts = cont_obs)
}

observed <- simulate_observations(
  latent$epi, latent$net, params, obs, symptom_times, contact_times, N
)

# -------------------------------
# 8. Save outputs
# -------------------------------
write_csv(latent$epi, "output/latent_epidemic_states.csv")
write_csv(latent$net, "output/latent_network_states.csv")
write_csv(latent$events, "output/latent_event_log.csv")
write_csv(observed$symptoms, "output/observed_symptoms.csv")
write_csv(observed$contacts, "output/observed_contacts.csv")

# Compact simulation metadata
sim_meta <- tibble(
  N = params$N,
  T = params$T,
  beta = params$beta,
  xi = params$xi,
  kappa = params$kappa,
  gamma = params$gamma,
  p_E = obs$p_E,
  p_I = obs$p_I,
  s = obs$s,
  c = obs$c,
  n_events = nrow(latent$events),
  n_epi_rows = nrow(latent$epi),
  n_net_rows = nrow(latent$net)
)

write_csv(sim_meta, "output/simulation_metadata.csv")

# -------------------------------
# 9. Reproducibility info
# -------------------------------
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
##############################################################






# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

bias_fun <- function(est, truth) mean(est - truth, na.rm = TRUE)
rmse_fun <- function(est, truth) sqrt(mean((est - truth)^2, na.rm = TRUE))
cover_fun <- function(lower, upper, truth) mean(lower <= truth & truth <= upper, na.rm = TRUE)
mae_fun <- function(x, y) mean(abs(x - y), na.rm = TRUE)

summarise_param <- function(df, param, truth_col = "truth") {
  truth <- unique(df[[truth_col]])
  stopifnot(length(truth) == 1)
  
  est_col <- "estimate"
  low_col <- "lower"
  upp_col <- "upper"
  
  tibble(
    parameter = param,
    truth = truth,
    bias = bias_fun(df[[est_col]], truth),
    rmse = rmse_fun(df[[est_col]], truth),
    coverage = cover_fun(df[[low_col]], df[[upp_col]], truth),
    mean_width = mean(df[[upp_col]] - df[[low_col]], na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# Expected input format
# ------------------------------------------------------------
# The script assumes you have a replicate-level results table
# called sim_results with the following columns:
#
# regime   : "High", "Moderate", or "Sparse"
# parameter: parameter name
# truth    : true parameter value
# estimate : point estimate from each replicate
# lower    : lower endpoint of 95% credible/confidence interval
# upper    : upper endpoint of 95% credible/confidence interval
#
# Example rows:
# regime parameter truth estimate lower upper
# High   beta      0.30  0.31     0.28  0.35
#
# If your output object is named differently, adapt below.
# ------------------------------------------------------------

# Placeholder: load replicate-level simulation output ----------
# sim_results <- read_csv("simulation_replicate_results.csv")

# Example structure check
# glimpse(sim_results)

# ------------------------------------------------------------
# Table 1: Epidemic parameters
# ------------------------------------------------------------
epi_params <- c("beta", "xi", "kappa", "gamma")

epi_table <- sim_results %>%
  filter(parameter %in% epi_params) %>%
  group_by(regime, parameter) %>%
  group_modify(~ summarise_param(.x, unique(.x$parameter))) %>%
  ungroup() %>%
  select(regime, parameter, truth, bias, rmse, coverage, mean_width) %>%
  arrange(parameter, regime)

write_csv(epi_table, "output/table_epidemic_recovery.csv")

epi_table_fmt <- epi_table %>%
  mutate(
    bias = round(bias, 4),
    rmse = round(rmse, 4),
    coverage = round(coverage, 3),
    mean_width = round(mean_width, 4)
  )

kable(epi_table_fmt, caption = "Recovery of epidemic parameters across simulation regimes.")

# ------------------------------------------------------------
# Table 2: Network parameters
# ------------------------------------------------------------
net_params <- c("eta_SS", "eta_SE", "eta_SI", "eta_SR",
                "omega_SS", "omega_SE", "omega_SI", "omega_SR")

net_table <- sim_results %>%
  filter(parameter %in% net_params) %>%
  group_by(parameter) %>%
  group_modify(~ summarise_param(.x, unique(.x$parameter))) %>%
  ungroup() %>%
  select(parameter, truth, bias, rmse, coverage, mean_width) %>%
  arrange(parameter)

write_csv(net_table, "output/table_network_recovery.csv")

net_table_fmt <- net_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

kable(net_table_fmt, caption = "Recovery of dynamic network parameters under partial observation.")

# ------------------------------------------------------------
# Table 3: Observation parameters
# ------------------------------------------------------------
obs_params <- c("p_E", "p_I", "s", "c")

obs_table <- sim_results %>%
  filter(parameter %in% obs_params) %>%
  group_by(parameter) %>%
  group_modify(~ summarise_param(.x, unique(.x$parameter))) %>%
  ungroup() %>%
  select(parameter, truth, bias, rmse, coverage, mean_width) %>%
  arrange(parameter)

write_csv(obs_table, "output/table_observation_recovery.csv")

obs_table_fmt <- obs_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

kable(obs_table_fmt, caption = "Estimation performance for observation-model parameters.")

# ------------------------------------------------------------
# Table 4: Latent-path recovery
# ------------------------------------------------------------
# Expected columns:
# regime, disease_accuracy, edge_accuracy, infection_time_mae, edge_change_mae
latent_summary <- sim_results %>%
  filter(parameter %in% c("latent_disease_accuracy", "latent_edge_accuracy",
                          "infection_time_mae", "edge_change_mae")) %>%
  select(regime, parameter, estimate) %>%
  pivot_wider(names_from = parameter, values_from = estimate) %>%
  rename(
    disease_accuracy = latent_disease_accuracy,
    edge_accuracy = latent_edge_accuracy,
    infection_time_mae = infection_time_mae,
    edge_change_mae = edge_change_mae
  ) %>%
  arrange(regime)

write_csv(latent_summary, "output/table_latent_recovery.csv")

latent_summary_fmt <- latent_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

kable(latent_summary_fmt, caption = "Recovery of latent epidemic and network histories under different observation regimes.")

# ------------------------------------------------------------
# Table 5: Sensitivity analysis
# ------------------------------------------------------------
# Expected rows for scenarios such as:
# "baseline", "diffuse_priors", "informative_priors", "low_specificity", "low_sensitivity"
sensitivity_table <- sim_results %>%
  filter(parameter %in% c("beta_rmse", "xi_rmse", "eta_rmse", "coverage")) %>%
  select(regime, parameter, estimate) %>%
  pivot_wider(names_from = parameter, values_from = estimate) %>%
  arrange(regime)

write_csv(sensitivity_table, "output/table_sensitivity.csv")

sensitivity_table_fmt <- sensitivity_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

kable(sensitivity_table_fmt, caption = "Sensitivity analysis for prior specification and contact misclassification.")

# ------------------------------------------------------------
# Optional publication-quality plots
# ------------------------------------------------------------

plot_param_recovery <- function(df, param_group, ylab) {
  df %>%
    filter(parameter %in% param_group) %>%
    ggplot(aes(x = parameter, y = estimate, fill = regime)) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.alpha = 0.3) +
    geom_hline(aes(yintercept = truth), linetype = "dashed", linewidth = 0.6) +
    facet_wrap(~ regime, scales = "free_y") +
    labs(x = NULL, y = ylab, fill = "Regime") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
}

# Example if replicate-level draws are available
# p1 <- plot_param_recovery(sim_results, epi_params, "Posterior estimate")
# ggsave("output/epi_recovery.png", p1, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------
# Session information for reproducibility
# ------------------------------------------------------------
sink("output/session_info.txt")
sessionInfo()
sink()
# ------------------------------------------------------------------------------
#################################################################################################################