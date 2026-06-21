#################################################################################################################
# ============================================================
# Complete-data likelihood fit for one simulated replicate
# Includes epidemic, observation, and full dynamic-network MLEs
# ============================================================

library(dplyr)
library(readr)
library(tibble)

fit_seir_dynnet_likelihood_full <- function(epi_file,
                                            event_file,
                                            net_file,
                                            symp_file,
                                            cont_file,
                                            T_end) {
  
  epi   <- read_csv(epi_file,   show_col_types = FALSE)
  events <- read_csv(event_file, show_col_types = FALSE)
  net   <- read_csv(net_file,   show_col_types = FALSE)
  symp  <- read_csv(symp_file,  show_col_types = FALSE)
  cont  <- read_csv(cont_file,  show_col_types = FALSE)
  
  if (!all(c("time", "id", "state") %in% names(epi))) stop("epi file must contain time, id, state")
  if (!all(c("time", "type") %in% names(events))) stop("event file must contain time, type")
  if (!all(c("time", "i", "j", "Aij") %in% names(net))) stop("net file must contain time, i, j, Aij")
  if (!all(c("time", "id", "state", "Y") %in% names(symp))) stop("symptom file must contain time, id, state, Y")
  if (!all(c("time", "i", "j", "A", "B") %in% names(cont))) stop("contact file must contain time, i, j, A, B")
  
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
  
  # ----------------------------
  # Epidemic sufficient stats
  # ----------------------------
  N_SE_int <- sum(events$type == "SE_int", na.rm = TRUE)
  N_SE_ext <- sum(events$type == "SE_ext", na.rm = TRUE)
  N_EI     <- sum(events$type == "EI", na.rm = TRUE)
  N_IR     <- sum(events$type == "IR", na.rm = TRUE)
  
  U_S <- state_time(epi, "S", T_end)
  U_E <- state_time(epi, "E", T_end)
  U_I <- state_time(epi, "I", T_end)
  
  beta_hat  <- ifelse(U_S > 0, N_SE_int / U_S, NA_real_)
  xi_hat    <- ifelse(U_S > 0, N_SE_ext / U_S, NA_real_)
  kappa_hat <- ifelse(U_E > 0, N_EI / U_E, NA_real_)
  gamma_hat <- ifelse(U_I > 0, N_IR / U_I, NA_real_)
  
  beta_se  <- ifelse(U_S > 0, sqrt(N_SE_int) / U_S, NA_real_)
  xi_se    <- ifelse(U_S > 0, sqrt(N_SE_ext) / U_S, NA_real_)
  kappa_se <- ifelse(U_E > 0, sqrt(N_EI) / U_E, NA_real_)
  gamma_se <- ifelse(U_I > 0, sqrt(N_IR) / U_I, NA_real_)
  
  # ----------------------------
  # Observation sufficient stats
  # ----------------------------
  M_E <- sum(symp$state == "E", na.rm = TRUE)
  R_E <- sum(symp$state == "E" & symp$Y == 1, na.rm = TRUE)
  M_I <- sum(symp$state == "I", na.rm = TRUE)
  R_I <- sum(symp$state == "I" & symp$Y == 1, na.rm = TRUE)
  
  M1 <- sum(cont$A == 1, na.rm = TRUE)
  R1 <- sum(cont$A == 1 & cont$B == 1, na.rm = TRUE)
  M0 <- sum(cont$A == 0, na.rm = TRUE)
  R0 <- sum(cont$A == 0 & cont$B == 0, na.rm = TRUE)
  
  pE_hat <- ifelse(M_E > 0, R_E / M_E, NA_real_)
  pI_hat <- ifelse(M_I > 0, R_I / M_I, NA_real_)
  s_hat   <- ifelse(M1 > 0, R1 / M1, NA_real_)
  c_hat    <- ifelse(M0 > 0, R0 / M0, NA_real_)
  
  pE_se <- ifelse(M_E > 0, sqrt(pE_hat * (1 - pE_hat) / M_E), NA_real_)
  pI_se <- ifelse(M_I > 0, sqrt(pI_hat * (1 - pI_hat) / M_I), NA_real_)
  s_se  <- ifelse(M1 > 0, sqrt(s_hat  * (1 - s_hat) / M1), NA_real_)
  c_se  <- ifelse(M0 > 0, sqrt(c_hat  * (1 - c_hat) / M0), NA_real_)
  
  # ----------------------------
  # Full dynamic-network MLEs
  # Need dyad-state exposure time and event counts by pair type
  # ----------------------------
  
  # Helper: state at time t for each node
  node_state_at <- function(t0) {
    epi_t <- epi %>% filter(time <= t0) %>%
      arrange(id, time) %>%
      group_by(id) %>%
      summarise(state = dplyr::last(state), .groups = "drop")
    epi_t$state[order(epi_t$id)]
  }
  
  # Helper: network at time t
  net_at <- function(t0, N) {
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
  
  N <- max(epi$id)
  pair_names <- c("SS","SE","SI","SR","EE","EI","ER","II","IR","RR")
  
  # Event counts
  form_counts <- setNames(as.list(rep(0, length(pair_names))), pair_names)
  diss_counts <- setNames(as.list(rep(0, length(pair_names))), pair_names)
  
  ev <- events %>% filter(type %in% c("A01", "A10"))
  if (nrow(ev) > 0) {
    for (k in seq_len(nrow(ev))) {
      t0 <- ev$time[k]
      st <- node_state_at(t0)
      i <- ev$i[k]; j <- ev$j[k]
      if (is.na(i) || is.na(j)) next
      p <- pair_key(st[i], st[j])
      if (!(p %in% pair_names)) next
      if (ev$type[k] == "A01") form_counts[[p]] <- form_counts[[p]] + 1
      if (ev$type[k] == "A10") diss_counts[[p]] <- diss_counts[[p]] + 1
    }
  }
  
  # Integrated risk sets for each dyad pair type
  # For each event interval, accumulate time spent in each dyad state and whether edge was absent/present.
  # This computes exact exposure from the latent piecewise-constant paths.
  exposure0 <- setNames(as.list(rep(0, length(pair_names))), pair_names)
  exposure1 <- setNames(as.list(rep(0, length(pair_names))), pair_names)
  
  times_all <- sort(unique(c(0, epi$time, net$time, T_end)))
  times_all <- times_all[times_all >= 0 & times_all <= T_end]
  if (length(times_all) < 2) times_all <- c(0, T_end)
  
  for (m in seq_len(length(times_all) - 1)) {
    tmid <- (times_all[m] + times_all[m + 1]) / 2
    dt <- times_all[m + 1] - times_all[m]
    if (dt <= 0) next
    
    st <- node_state_at(tmid)
    A  <- net_at(tmid, N)
    
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        p <- pair_key(st[i], st[j])
        if (!(p %in% pair_names)) next
        if (A[i, j] == 0L) {
          exposure0[[p]] <- exposure0[[p]] + dt
        } else {
          exposure1[[p]] <- exposure1[[p]] + dt
        }
      }
    }
  }
  
  eta_names <- pair_names
  omega_names <- pair_names
  
  eta_hat <- sapply(eta_names, function(p) {
    if (exposure0[[p]] > 0) form_counts[[p]] / exposure0[[p]] else NA_real_
  })
  omega_hat <- sapply(omega_names, function(p) {
    if (exposure1[[p]] > 0) diss_counts[[p]] / exposure1[[p]] else NA_real_
  })
  
  eta_se <- sapply(eta_names, function(p) {
    if (exposure0[[p]] > 0) sqrt(form_counts[[p]]) / exposure0[[p]] else NA_real_
  })
  omega_se <- sapply(omega_names, function(p) {
    if (exposure1[[p]] > 0) sqrt(diss_counts[[p]]) / exposure1[[p]] else NA_real_
  })
  
  z <- qnorm(0.975)
  
  fit_summary <- bind_rows(
    tibble(parameter = c("beta","xi","kappa","gamma","p_E","p_I","s","c"),
           estimate = c(beta_hat, xi_hat, kappa_hat, gamma_hat, pE_hat, pI_hat, s_hat, c_hat),
           se = c(beta_se, xi_se, kappa_se, gamma_se, pE_se, pI_se, s_se, c_se),
           lower = estimate - z * se,
           upper = estimate + z * se),
    tibble(parameter = paste0("eta_", names(eta_hat)),
           estimate = as.numeric(eta_hat),
           se = as.numeric(eta_se),
           lower = estimate - z * se,
           upper = estimate + z * se),
    tibble(parameter = paste0("omega_", names(omega_hat)),
           estimate = as.numeric(omega_hat),
           se = as.numeric(omega_se),
           lower = estimate - z * se,
           upper = estimate + z * se)
  )
  
  list(
    fit_summary = fit_summary,
    sufficient_stats = list(
      N_SE_int = N_SE_int, N_SE_ext = N_SE_ext, N_EI = N_EI, N_IR = N_IR,
      U_S = U_S, U_E = U_E, U_I = U_I,
      M_E = M_E, R_E = R_E, M_I = M_I, R_I = R_I,
      M1 = M1, R1 = R1, M0 = M0, R0 = R0,
      form_counts = form_counts,
      diss_counts = diss_counts,
      exposure0 = exposure0,
      exposure1 = exposure1
    )
  )
}

#Example use:
mle_fit_all_parameters <- fit_seir_dynnet_likelihood_full(
  epi_file   = "output/latent_epidemic_states.csv",
  event_file = "output/latent_event_log.csv",
  net_file   = "output/latent_network_states.csv",
  symp_file  = "output/observed_symptoms.csv",
  cont_file  = "output/observed_contacts.csv",
  T_end      = 30
)
print(res$fit_summary, n=28)
write.csv(res, "output/fit_results_all_parameters.csv")
# -------------------------------------------------------------------------------
#################################################################################################################