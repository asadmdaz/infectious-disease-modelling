#################################################################################################################
# ------------------------------------------------------------
# Likelihood fit for one simulated replicate
# Input: latent epidemic/network CSVs + observed data CSVs
# Output: parameter estimates and 95% Wald intervals
# ------------------------------------------------------------

library(dplyr)
library(readr)
library(tibble)

fit_seir_dynnet_likelihood <- function(epi_file,
                                       event_file,
                                       net_file,
                                       symp_file,
                                       cont_file,
                                       T_end) {
  
  epi <- read_csv(epi_file, show_col_types = FALSE)
  events <- read_csv(event_file, show_col_types = FALSE)
  net <- read_csv(net_file, show_col_types = FALSE)
  symp <- read_csv(symp_file, show_col_types = FALSE)
  cont <- read_csv(cont_file, show_col_types = FALSE)
  
  # ----------------------------
  # Complete-data sufficient stats
  # ----------------------------
  N_SE_int <- sum(events$type == "SE_int", na.rm = TRUE)
  N_SE_ext <- sum(events$type == "SE_ext", na.rm = TRUE)
  N_EI <- sum(events$type == "EI", na.rm = TRUE)
  N_IR <- sum(events$type == "IR", na.rm = TRUE)
  
  # Integrated state times from latent epidemic history
  state_time <- function(df, state_name, T_end) {
    df <- df %>% arrange(id, time)
    total <- 0
    for (i in unique(df$id)) {
      d <- df %>% filter(id == i)
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
  
  U_S <- state_time(epi, "S", T_end)
  U_E <- state_time(epi, "E", T_end)
  U_I <- state_time(epi, "I", T_end)
  
  # Observational sufficient stats
  M_E <- sum(symp$state == "E", na.rm = TRUE)
  R_E <- sum(symp$state == "E" & symp$Y == 1, na.rm = TRUE)
  M_I <- sum(symp$state == "I", na.rm = TRUE)
  R_I <- sum(symp$state == "I" & symp$Y == 1, na.rm = TRUE)
  
  M1 <- sum(cont$A == 1, na.rm = TRUE)
  R1 <- sum(cont$A == 1 & cont$B == 1, na.rm = TRUE)
  M0 <- sum(cont$A == 0, na.rm = TRUE)
  R0 <- sum(cont$A == 0 & cont$B == 0, na.rm = TRUE)
  
  # ----------------------------
  # Closed-form complete-data MLEs
  # ----------------------------
  beta_hat  <- ifelse(U_S > 0, N_SE_int / U_S, NA_real_)
  xi_hat    <- ifelse(U_S > 0, N_SE_ext / U_S, NA_real_)
  kappa_hat <- ifelse(U_E > 0, N_EI / U_E, NA_real_)
  gamma_hat <- ifelse(U_I > 0, N_IR / U_I, NA_real_)
  
  pE_hat <- ifelse(M_E > 0, R_E / M_E, NA_real_)
  pI_hat <- ifelse(M_I > 0, R_I / M_I, NA_real_)
  s_hat  <- ifelse(M1 > 0, R1 / M1, NA_real_)
  c_hat  <- ifelse(M0 > 0, R0 / M0, NA_real_)
  
  # Wald SEs from binomial/Poisson approximations
  beta_se  <- ifelse(U_S > 0, sqrt(N_SE_int) / U_S, NA_real_)
  xi_se    <- ifelse(U_S > 0, sqrt(N_SE_ext) / U_S, NA_real_)
  kappa_se <- ifelse(U_E > 0, sqrt(N_EI) / U_E, NA_real_)
  gamma_se <- ifelse(U_I > 0, sqrt(N_IR) / U_I, NA_real_)
  
  pE_se <- ifelse(M_E > 0, sqrt(pE_hat * (1 - pE_hat) / M_E), NA_real_)
  pI_se <- ifelse(M_I > 0, sqrt(pI_hat * (1 - pI_hat) / M_I), NA_real_)
  s_se  <- ifelse(M1 > 0, sqrt(s_hat  * (1 - s_hat) / M1), NA_real_)
  c_se  <- ifelse(M0 > 0, sqrt(c_hat  * (1 - c_hat) / M0), NA_real_)
  
  z <- qnorm(0.975)
  
  fit_summary <- tibble(
    parameter = c("beta", "xi", "kappa", "gamma", "p_E", "p_I", "s", "c"),
    estimate  = c(beta_hat, xi_hat, kappa_hat, gamma_hat, pE_hat, pI_hat, s_hat, c_hat),
    se        = c(beta_se, xi_se, kappa_se, gamma_se, pE_se, pI_se, s_se, c_se),
    lower     = estimate - z * se,
    upper     = estimate + z * se
  )
  
  list(
    fit_summary = fit_summary,
    sufficient_stats = list(
      N_SE_int = N_SE_int, N_SE_ext = N_SE_ext, N_EI = N_EI, N_IR = N_IR,
      U_S = U_S, U_E = U_E, U_I = U_I,
      M_E = M_E, R_E = R_E, M_I = M_I, R_I = R_I,
      M1 = M1, R1 = R1, M0 = M0, R0 = R0
    )
  )
}

# Example use
res <- fit_seir_dynnet_likelihood(
  epi_file   = "output/latent_epidemic_states.csv",
  event_file = "output/latent_event_log.csv",
  net_file   = "output/latent_network_states.csv",
  symp_file  = "output/observed_symptoms.csv",
  cont_file  = "output/observed_contacts.csv",
  T_end      = 30
)
print(res$fit_summary)
write.csv(res, "output/fit_results_parameters.csv")
#################################################################################################################