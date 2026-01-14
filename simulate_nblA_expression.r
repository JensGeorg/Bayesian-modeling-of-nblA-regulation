library(cmdstanr)
library(posterior)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# 1. LOAD PARAMETERS FROM FIT


  draws_df <- as_draws_df(fit$draws())
  
  # Helper to get mean of a parameter
  get_param <- function(p) {
    if (p %in% names(draws_df)) return(mean(draws_df[[p]]))
    return(NA)
  }
  
  # Create parameter list
  parms <- list(
    # Kinetics
    Kd_sRNA             = get_param("Kd_sRNA"),
    decay_R_free        = get_param("decay_R_free"),
    decay_complex_sRNA  = get_param("decay_complex_sRNA"),
    decay_S_free        = get_param("decay_S_free"),
    decay_asR           = get_param("decay_asR"),
    
    # Synthesis
    base_syn            = get_param("base_syn"),
    base_syn_S          = get_param("base_syn_S"),
    synthesis_leak      = get_param("synthesis_leak"),
    
    # Interference
    term_rate_asR_on_RNA = get_param("term_rate_asR_on_RNA"),
    term_rate_RNA_on_asR = get_param("term_rate_RNA_on_asR"),
    
    # Forcing Data Points (Parameters estimated by Stan)
    f_RNA_9             = get_param("f_RNA_9"),
    f_RNA_24            = get_param("f_RNA_24"),
    f_asRNA_9           = get_param("f_asRNA_9"),
    f_asRNA_24          = get_param("f_asRNA_24"),
    f_sRNA_9            = get_param("f_sRNA_9"),
    f_sRNA_24           = get_param("f_sRNA_24")
  )



#2. DEFINE HELPER FUNCTIONS From Stan Model

# Softplus function (log1p_exp) for smooth clamping
log1p_exp <- function(x) {
  # Avoid overflow for large x
  ifelse(x > 40, x, log1p(exp(x)))
}

# Interpolation function
interp_forcing <- function(t, v0, v9, v24) {
  if (t <= 0) return(v0)
  else if (t <= 9) return(v0 + (v9 - v0) * t / 9.0)
  else if (t <= 24) return(v9 + (v24 - v9) * (t - 9.0) / 15.0)
  else return(v24)
}

# QSSA Complex Solver
solve_complex <- function(R_tot, S_tot, Kd) {
  Rt <- pmax(0, R_tot)
  St <- pmax(0, S_tot)
  sum_all <- Rt + St + Kd
  discrim <- sum_all^2 - 4.0 * Rt * St
  
  C_raw <- (sum_all - sqrt(pmax(discrim, 1e-12))) / 2.0
  
  # Apply soft clamp logic from Stan (source: 9)
  return(log1p_exp(10.0 * C_raw) / 10.0)
}

# Analytic Steady State Solver (for Initial Conditions)
solve_steady_state_analytic <- function(syn_R, syn_S, Kd, k_R, k_S, k_C) {
  # Coefficients for Quadratic: a*C^2 + b*C + c = 0
  a <- k_C^2
  b <- -1.0 * (k_C * (syn_R + syn_S) + k_R * k_S * Kd)
  c <- syn_R * syn_S
  
  discriminant <- b^2 - 4.0 * a * c
  
  # Solve Quadratic
  C <- (-b - sqrt(pmax(0.0, discriminant))) / (2.0 * a)
  C <- pmax(0.0, C)
  
  R_free <- pmax(0.0, (syn_R - k_C * C) / k_R)
  S_free <- pmax(0.0, (syn_S - k_C * C) / k_S)
  
  return(c(R_tot = R_free + C, S_tot = S_free + C, asR_tot = 0)) # asR handled separately
}


# 3. DEFINE ODE MODEL


rna_ode_qssa <- function(t, y, parms) {
 
  R_tot <- y["R_tot"]
  S_tot <- y["S_tot"]
  asR_tot <- y["asR_tot"]
  

  is_sRNA_KO <- parms$is_sRNA_KO
  is_asRNA_KO <- parms$is_asRNA_KO
  is_heterocyst <- parms$is_heterocyst
  is_st9_0 <- parms$is_st9_0
  

  f_R_val <- interp_forcing(t, 1.0, parms$f_RNA_9, parms$f_RNA_24)
  f_asR_val <- interp_forcing(t, 102.222, parms$f_asRNA_9, parms$f_asRNA_24)
  f_S_val <- interp_forcing(t, 1.0, parms$f_sRNA_9, parms$f_sRNA_24)
  
  # Heterocyst logic 
  if (is_heterocyst) {
    f_R_val <- interp_forcing(t, 1.0, parms$f_RNA_9* 3.0, parms$f_RNA_24 * 3.0)
  }
  if(is_st9_0){
	 f_S_val <- interp_forcing(t, 1.0, 0, parms$f_sRNA_24)
  }
  
  #Synthesis Rates
  syn_R_raw <- parms$base_syn * f_R_val
  syn_asR_raw <- parms$base_syn * f_asR_val
  syn_S_raw <- parms$base_syn_S * f_S_val
  
  if (is_sRNA_KO) syn_S_raw <- 0
  if (is_asRNA_KO) syn_asR_raw <- 0
  
  # Interference (TI)
  # R synthesis (blocked by asR)
  diff_ode <- syn_R_raw - parms$term_rate_asR_on_RNA * syn_asR_raw
  # Smooth approximation from Stan
  syn_R_eff_main <- log1p_exp(5.0 * diff_ode) / 5.0
  # Add Leak
  syn_R_eff_total <- syn_R_eff_main + parms$synthesis_leak
  
  # asR synthesis (blocked by R)
  diff_ode_as <- syn_asR_raw - parms$term_rate_RNA_on_asR * syn_R_raw
  syn_asR_eff <- log1p_exp(5.0 * diff_ode_as) / 5.0
  
  # Solve Complex (QSSA)
  C_sRNA <- solve_complex(R_tot, S_tot, parms$Kd_sRNA)
  R_free <- pmax(0, R_tot - C_sRNA)
  S_free <- pmax(0, S_tot - C_sRNA)
  
  # Derivatives 
  dR <- syn_R_eff_total - parms$decay_R_free * R_free - parms$decay_complex_sRNA * C_sRNA
  dS <- syn_S_raw - parms$decay_S_free * S_free - parms$decay_complex_sRNA * C_sRNA
  dasR <- syn_asR_eff - parms$decay_asR * asR_tot
  
  return(list(c(dR, dS, dasR)))
}

# 4. RUN SIMULATIONS

# Simulation Settings
times <- seq(0, 15, by = 0.1)
scenarios <- c("Wild Type", "sRNA Knockout", "asRNA Knockout", "Heterocyst")
results_list <- list()

# 1. Calculate t=0 Initial Conditions for WT (Used as baseline)
#    Recalculate these for each scenario to ensure steady state at t=0
calculate_init <- function(p, is_s_ko=F, is_as_ko=F) {
  
  # t=0 Synthesis
  syn_R_raw_0 <- p$base_syn * 1.0
  syn_asR_raw_0 <- p$base_syn * 102.222
  syn_S_raw_0 <- p$base_syn_S * 1.0
  
  if (is_s_ko) syn_S_raw_0 <- 0
  if (is_as_ko) syn_asR_raw_0 <- 0
  
  # TI at t=0
  diff_term <- syn_R_raw_0 - p$term_rate_asR_on_RNA * syn_asR_raw_0
  syn_R_main_0 <- log1p_exp(5.0 * diff_term) / 5.0
  syn_R_total_0 <- syn_R_main_0 + p$synthesis_leak
  
  syn_asR_eff_0 <- pmax(1e-6, syn_asR_raw_0 - p$term_rate_RNA_on_asR * syn_R_raw_0)
  
  # Steady State solution
  ss_vals <- solve_steady_state_analytic(
    syn_R_total_0, syn_S_raw_0, p$Kd_sRNA,
    p$decay_R_free, p$decay_S_free, p$decay_complex_sRNA
  )
  
  # asR steady state
  asR_ss <- syn_asR_eff_0 / p$decay_asR
  
  return(c(R_tot = ss_vals[["R_tot"]], S_tot = ss_vals[["S_tot"]], asR_tot = asR_ss))
}


for (scen in scenarios) {
  cat("Simulating:", scen, "\n")
  
  # Set Scenario Flags
  p_local <- parms
  p_local$is_sRNA_KO <- FALSE
  p_local$is_asRNA_KO <- FALSE
  p_local$is_heterocyst <- FALSE
  p_local$is_st9_0 <- FALSE
  
  
  if (scen == "sRNA Knockout") p_local$is_sRNA_KO <- TRUE
  if (scen == "asRNA Knockout") p_local$is_asRNA_KO <- TRUE
  if (scen == "Heterocyst") p_local$is_heterocyst <- TRUE
  
  # Calculate Initial Conditions specific to this scenario
  # Assuming the system starts at steady state for that condition
  y0 <- calculate_init(p_local, p_local$is_sRNA_KO, p_local$is_asRNA_KO)
  
  # Run ODE
  out <- ode(y = y0, times = times, func = rna_ode_qssa, parms = p_local)
  
  # Process Output
  df <- as.data.frame(out)
  df$Scenario <- scen
  
  df$C_sRNA <- mapply(solve_complex, df$R_tot, df$S_tot, MoreArgs = list(Kd = p_local$Kd_sRNA))
  df$R_free <- pmax(0, df$R_tot - df$C_sRNA)
  df$S_free <- pmax(0, df$S_tot - df$C_sRNA)
  
  results_list[[scen]] <- df
}

combined_results <- bind_rows(results_list)

# 5. PLOTTING

# A. Total RNA Plot 
plot_rna <- ggplot(combined_results, aes(x = time, y = R_tot, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "nblA Total RNA Dynamics (QSSA Model)",
    subtitle = "Includes sRNA complex formation and Interference",
    x = "Time [h]",
    y = "Total RNA Concentration [au]"
  ) +
  theme_bw(base_size = 14) +
  ggsci::scale_color_npg() +
  theme(legend.position = "bottom")

# B. Faceted View (RNA, sRNA, asRNA)
long_df <- combined_results %>%
  select(time, Scenario, R_tot, S_tot, asR_tot) %>%
  pivot_longer(cols = c(R_tot, S_tot, asR_tot), names_to = "Species", values_to = "Concentration") %>%
  mutate(Species = factor(Species, levels = c("R_tot", "S_tot", "asR_tot"), 
                          labels = c("mRNA (nblA)", "sRNA", "asRNA")))

plot_faceted <- ggplot(long_df, aes(x = time, y = Concentration, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Species, scales = "free_y") +
  labs(title = "Dynamics of all Species") +
  theme_bw() +
  ggsci::scale_color_npg() +
  theme(legend.position = "bottom")

# Display
# print(plot_rna)
# print(plot_faceted) 

# Save
ggsave("nbla_simulation_qssa.pdf", plot_rna, width = 8, height = 6)
ggsave("nbla_simulation_qssa_all_species.pdf", plot_faceted, width = 8, height = 6)
ggsave("nbla_simulation_qssa.png", plot_rna, width = 8, height = 6)
ggsave("nbla_simulation_qssa_all_species.png", plot_faceted, width = 8, height = 6)

# A. Total RNA Plot (Matches your request)
plot_rna <- ggplot(combined_results, aes(x = time, y = R_tot, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "nblA Total RNA Dynamics (QSSA Model)",
    subtitle = "Includes sRNA complex formation and Interference",
    x = "Time [h]",
    y = "Total RNA Concentration [au]"
  ) +
  xlim(0,4) +
  ylim(0,1500) +
  theme_bw(base_size = 14) +
  ggsci::scale_color_npg() +
  theme(legend.position = "bottom")

# B. Faceted View (RNA, sRNA, asRNA)
long_df <- combined_results %>%
  select(time, Scenario, R_tot, S_tot, asR_tot) %>%
  pivot_longer(cols = c(R_tot, S_tot, asR_tot), names_to = "Species", values_to = "Concentration") %>%
  mutate(Species = factor(Species, levels = c("R_tot", "S_tot", "asR_tot"), 
                          labels = c("mRNA (nblA)", "sRNA", "asRNA")))

plot_faceted <- ggplot(long_df, aes(x = time, y = Concentration, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Species, scales = "free_y") +
  labs(title = "Dynamics of all Species") +
  theme_bw() +
  ggsci::scale_color_npg() +
  theme(legend.position = "bottom")

# Display
# print(plot_rna)
# print(plot_faceted) 

# Save
ggsave("nbla_simulation_qssa_0to4.pdf", plot_rna, width = 8, height = 6)
ggsave("nbla_simulation_qssa_all_species_0to4.pdf", plot_faceted, width = 8, height = 6)
ggsave("nbla_simulation_qssa_0to4.png", plot_rna, width = 8, height = 6)
ggsave("nbla_simulation_qssa_all_species_0to4.png", plot_faceted, width = 8, height = 6)

