#################################################################################################################
# ============================================================
# Bayesian inference for SEIR + dynamic network model
# Complete-data latent path posterior sampler
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
})

set.seed(20260603)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------
rgamma_post <- function(a0, b0, n_events, exposure) {
  rgamma(1, shape = a0 + n_events, rate = b0 + exposure)
}

rbeta_post <- function(a0, b0, successes, trials) {
  rbeta(1, shape1 = a0 + successes, shape2 = b0 + trials - successes)
}

ci_quant <- function(x, probs = c(0.025, 0.975)) {
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}

state_time <- function(df, state_name, T_end) {
  df <- df %>% arrange(id, time)
  total <- 0
  for (i in unique(df$id)) {
    d <- df %>% filter(id == i) %>% arrange(time)
    if (nrow(d) < 1) next
    tt <- d$time
    st <- d$state
    if (length(tt) >= 2) {
      for (k in seq_len(length(tt) - 1)) {
        if (st[k] == state_name) total <- total + (tt[k + 1] - tt[k])
      }
    }
    if (st[length(st)] == state_name && max(tt) < T_end) {
      total <- total + (T_end - max(tt))
    }
  }
  total
}

pair_key <- function(a, b) paste(sort(c(a, b)), collapse = "")

node_state_at <- function(epi, t0) {
  epi_t <- epi %>%
    filter(time <= t0) %>%
    arrange(id, time) %>%
    group_by(id) %>%
    summarise(state = dplyr::last(state), .groups = "drop")
  epi_t$state[order(epi_t$id)]
}

net_at <- function(net, t0, N) {
  A <- matrix(0L, N, N)
  d <- net %>% filter(time <= t0)
  if (nrow(d) == 0) return(A)
  latest <- d %>%
    arrange(i, j, time) %>%
    group_by(i, j) %>%
    summarise(Aij = dplyr::last(Aij), .groups = "drop")
  for (k in seq_len(nrow(latest))) {
    i <- latest$i[k]; j <- latest$j[k]
    A[i, j] <- latest$Aij[k]
    A[j, i] <- latest$Aij[k]
  }
  A
}

# ------------------------------------------------------------
# Posterior sampler for one replicate
# ------------------------------------------------------------
bayes_fit_seir_dynnet <- function(epi_file,
                                  event_file,
                                  net_file,
                                  symp_file,
                                  cont_file,
                                  meta_file,
                                  n_iter = 4000,
                                  burn = 1000,
                                  thin = 5) {
  
  epi    <- read_csv(epi_file,   show_col_types = FALSE)
  events <- read_csv(event_file, show_col_types = FALSE)
  net    <- read_csv(net_file,   show_col_types = FALSE)
  symp   <- read_csv(symp_file,  show_col_types = FALSE)
  cont   <- read_csv(cont_file,  show_col_types = FALSE)
  meta   <- read_csv(meta_file,  show_col_types = FALSE)
  
  T_end <- meta$T[1]
  N <- max(epi$id)
  
  # ----------------------------------------------------------
  # Sufficient statistics from complete latent data
  # Epidemic
  # ----------------------------------------------------------
  N_SE_int <- sum(events$type == "SE_int", na.rm = TRUE)
  N_SE_ext <- sum(events$type == "SE_ext", na.rm = TRUE)
  N_EI     <- sum(events$type == "EI", na.rm = TRUE)
  N_IR     <- sum(events$type == "IR", na.rm = TRUE)
  
  U_S <- state_time(epi, "S", T_end)
  U_E <- state_time(epi, "E", T_end)
  U_I <- state_time(epi, "I", T_end)
  
  # Observation model
  M_E <- sum(symp$state == "E", na.rm = TRUE)
  R_E <- sum(symp$state == "E" & symp$Y == 1, na.rm = TRUE)
  M_I <- sum(symp$state == "I", na.rm = TRUE)
  R_I <- sum(symp$state == "I" & symp$Y == 1, na.rm = TRUE)
  
  M1 <- sum(cont$A == 1, na.rm = TRUE)
  R1 <- sum(cont$A == 1 & cont$B == 1, na.rm = TRUE)
  M0 <- sum(cont$A == 0, na.rm = TRUE)
  R0 <- sum(cont$A == 0 & cont$B == 0, na.rm = TRUE)
  
  # ----------------------------------------------------------
  # Network sufficient statistics
  # ----------------------------------------------------------
  pair_names <- c("SS","SE","SI","SR","EE","EI","ER","II","IR","RR")
  form_counts <- setNames(as.list(rep(0, length(pair_names))), pair_names)
  diss_counts <- setNames(as.list(rep(0, length(pair_names))), pair_names)
  exp0 <- setNames(as.list(rep(0, length(pair_names))), pair_names)
  exp1 <- setNames(as.list(rep(0, length(pair_names))), pair_names)
  
  ev <- events %>% filter(type %in% c("A01", "A10"))
  if (nrow(ev) > 0) {
    for (k in seq_len(nrow(ev))) {
      t0 <- ev$time[k]
      st <- node_state_at(epi, t0)
      i <- ev$i[k]; j <- ev$j[k]
      if (is.na(i) || is.na(j)) next
      p <- pair_key(st[i], st[j])
      if (p %in% pair_names) {
        if (ev$type[k] == "A01") form_counts[[p]] <- form_counts[[p]] + 1
        if (ev$type[k] == "A10") diss_counts[[p]] <- diss_counts[[p]] + 1
      }
    }
  }
  
  times_all <- sort(unique(c(0, epi$time, net$time, T_end)))
  if (length(times_all) < 2) times_all <- c(0, T_end)
  
  for (m in seq_len(length(times_all) - 1)) {
    tmid <- (times_all[m] + times_all[m + 1]) / 2
    dt <- times_all[m + 1] - times_all[m]
    if (dt <= 0) next
    st <- node_state_at(epi, tmid)
    A <- net_at(net, tmid, N)
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        p <- pair_key(st[i], st[j])
        if (!(p %in% pair_names)) next
        if (A[i, j] == 0L) exp0[[p]] <- exp0[[p]] + dt else exp1[[p]] <- exp1[[p]] + dt
      }
    }
  }
  
  # ----------------------------------------------------------
  # Priors
  # ----------------------------------------------------------
  pri <- list(
    beta = c(1, 1),
    xi = c(1, 1),
    kappa = c(1, 1),
    gamma = c(1, 1),
    p_E = c(1, 1),
    p_I = c(1, 1),
    s = c(1, 1),
    c = c(1, 1)
  )
  pri_net <- list(
    eta = c(1, 1),
    omega = c(1, 1)
  )
  
  # ----------------------------------------------------------
  # MCMC storage
  # ----------------------------------------------------------
  keep_idx <- seq(burn + 1, n_iter, by = thin)
  draws <- vector("list", length(keep_idx))
  
  # Initialize at complete-data MLEs
  beta  <- ifelse(U_S > 0, N_SE_int / U_S, 0.1)
  xi    <- ifelse(U_S > 0, N_SE_ext / U_S, 0.05)
  kappa <- ifelse(U_E > 0, N_EI / U_E, 0.2)
  gamma <- ifelse(U_I > 0, N_IR / U_I, 0.2)
  pE <- ifelse(M_E > 0, R_E / M_E, 0.6)
  pI <- ifelse(M_I > 0, R_I / M_I, 0.8)
  s_obs <- ifelse(M1 > 0, R1 / M1, 0.9)
  c_obs <- ifelse(M0 > 0, R0 / M0, 0.95)
  
  eta <- setNames(rep(0.05, length(pair_names)), pair_names)
  omega <- setNames(rep(0.05, length(pair_names)), pair_names)
  
  # ----------------------------------------------------------
  # Gibbs sampler under complete-data conjugacy
  # ----------------------------------------------------------
  idx <- 0
  for (iter in 1:n_iter) {
    # Epidemic parameters
    beta  <- rgamma_post(pri$beta[1],  pri$beta[2],  N_SE_int, U_S)
    xi    <- rgamma_post(pri$xi[1],    pri$xi[2],    N_SE_ext, U_S)
    kappa <- rgamma_post(pri$kappa[1],  pri$kappa[2], N_EI,     U_E)
    gamma <- rgamma_post(pri$gamma[1],  pri$gamma[2], N_IR,     U_I)
    
    # Observation parameters
    pE    <- rbeta_post(pri$p_E[1], pri$p_E[2], R_E, M_E)
    pI    <- rbeta_post(pri$p_I[1], pri$p_I[2], R_I, M_I)
    s_obs <- rbeta_post(pri$s[1],   pri$s[2],   R1,  M1)
    c_obs <- rbeta_post(pri$c[1],   pri$c[2],   R0,  M0)
    
    # Network parameters
    for (p in pair_names) {
      eta[p]   <- rgamma_post(pri_net$eta[1],   pri_net$eta[2],   form_counts[[p]], exp0[[p]])
      omega[p] <- rgamma_post(pri_net$omega[1], pri_net$omega[2], diss_counts[[p]], exp1[[p]])
    }
    
    if (iter %in% keep_idx) {
      idx <- idx + 1
      row <- tibble(
        beta = beta, xi = xi, kappa = kappa, gamma = gamma,
        p_E = pE, p_I = pI, s = s_obs, c = c_obs
      )
      for (p in pair_names) {
        row[[paste0("eta_", p)]] <- eta[p]
        row[[paste0("omega_", p)]] <- omega[p]
      }
      draws[[idx]] <- row
    }
  }
  
  post <- bind_rows(draws)
  
  fit_summary <- post %>%
    summarise(across(everything(),
                     list(
                       estimate = ~mean(.x, na.rm = TRUE),
                       lower = ~quantile(.x, 0.025, na.rm = TRUE),
                       upper = ~quantile(.x, 0.975, na.rm = TRUE)
                     ),
                     .names = "{.col}_{.fn}")) %>%
    tidyr::pivot_longer(everything(),
                        names_to = c("parameter", ".value"),
                        names_pattern = "^(.*)_(estimate|lower|upper)$") %>%
    select(parameter, estimate, lower, upper)
  
  list(
    posterior_draws = post,
    fit_summary = fit_summary,
    sufficient_stats = list(
      N_SE_int = N_SE_int, N_SE_ext = N_SE_ext, N_EI = N_EI, N_IR = N_IR,
      U_S = U_S, U_E = U_E, U_I = U_I,
      M_E = M_E, R_E = R_E, M_I = M_I, R_I = R_I,
      M1 = M1, R1 = R1, M0 = M0, R0 = R0,
      form_counts = form_counts,
      diss_counts = diss_counts,
      exp0 = exp0,
      exp1 = exp1
    )
  )
}

# Example use:
bayes_fit_results <- bayes_fit_seir_dynnet(
  epi_file = "output/latent_epidemic_states.csv",
  event_file = "output/latent_event_log.csv",
  net_file = "output/latent_network_states.csv",
  symp_file = "output/observed_symptoms.csv",
  cont_file = "output/observed_contacts.csv",
  meta_file = "output/simulation_metadata.csv", 
  n_iter = 4000, burn = 1000, thin = 5
)

write.csv(res[1], "output/fit_results_all_parameter_bayesian_posterior_draws.csv")
write.csv(res[2], "output/fit_results_all_parameter_bayesian_fit_summary.csv")
write.csv(res[3], "output/fit_results_all_parameter_bayesian_sufficient_stats.csv")
#----------------------------------------------------------------------------------
#################################################################################################################
