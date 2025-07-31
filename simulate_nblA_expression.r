library(rstan)
library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

params_draws <- rstan::extract(fit)

params_mean <- list()

for (p_name in names(params_draws)) {
  p_draws <- params_draws[[p_name]]
  if (is.matrix(p_draws)) {
    params_mean[[p_name]] <- apply(p_draws, 2, mean)
  } else {
    params_mean[[p_name]] <- mean(p_draws)
  }
}

stan_data_ts <- list(
  rna_ts = c(0, 9, 24),
  asrna_ts = c(0, 9, 24),
  srna_ts = c(0, 6, 8, 12, 24)
)

rna_hetero_interpolator <- approxfun(c(0, 24), c(params_mean$f_RNA_vals[1], params_mean$f_RNA_vals[3] * 3), rule = 2)
rna_interpolator <- approxfun(stan_data_ts$rna_ts, params_mean$f_RNA_vals, rule = 2)
asrna_interpolator <- approxfun(stan_data_ts$asrna_ts, params_mean$f_asRNA_vals, rule = 2)
srna_interpolator <- approxfun(stan_data_ts$srna_ts, params_mean$f_sRNA_vals, rule = 2)

params_mean$k_on_sRNA <- params_mean$k_on_sRNA_gq

rna_dynamic_model_scenarios <- function(t, y, parms) {
  names(y) <- c("RNA_free", "RNA_complex_sRNA", "asRNA_free", "sRNA_free")
  y[which(y < 0)] <- 0

  RNA_hetero_synth_t_fitted <- parms$base_RNA_synth * parms$rna_hetero_interpolator(t)
  RNA_synth_t_fitted <- parms$base_RNA_synth * parms$rna_interpolator(t)
  asRNA_synth_t_fitted <- parms$base_RNA_synth * parms$asrna_interpolator(t)
  sRNA_synth_t_fitted <- parms$base_RNA_synth_sRNA * parms$srna_interpolator(t)
  
  asRNA_synth_t <- if (parms$scenario == "asRNA Knockout") 0 else asRNA_synth_t_fitted
  sRNA_synth_t <- if (parms$scenario == "sRNA Knockout") 0 else sRNA_synth_t_fitted
  
  if (parms$scenario == "heterocyst") {
    RNA_synth_t_fitted <- RNA_hetero_synth_t_fitted
  }
  
  RNA_synth_t <- RNA_synth_t_fitted
  
  effective_RNA_synth <- max(0, RNA_synth_t - asRNA_synth_t * parms$termination_rate_RNA)
  effective_asRNA_synth <- max(0, asRNA_synth_t - RNA_synth_t * parms$termination_rate_aRNA)

  dydt_RNA_free <- effective_RNA_synth - parms$k_on_sRNA * y["sRNA_free"] * y["RNA_free"] - parms$RNA_decay * y["RNA_free"] + parms$k_off_sRNA * y["RNA_complex_sRNA"]
  dydt_RNA_complex_sRNA <- parms$k_on_sRNA * y["sRNA_free"] * y["RNA_free"] - parms$k_off_sRNA * y["RNA_complex_sRNA"] - parms$complex_decay_sRNA * y["RNA_complex_sRNA"]
  dydt_asRNA_free <- effective_asRNA_synth - parms$asRNA_decay * y["asRNA_free"]
  dydt_sRNA_free <- sRNA_synth_t - parms$sRNA_decay * y["sRNA_free"] - parms$k_on_sRNA * y["sRNA_free"] * y["RNA_free"] + parms$k_off_sRNA * y["RNA_complex_sRNA"]

  return(list(c(dydt_RNA_free, dydt_RNA_complex_sRNA, dydt_asRNA_free, dydt_sRNA_free)))
}

scenarios_to_run <- c("Wild Type", "sRNA Knockout", "asRNA Knockout", "heterocyst")
all_scenario_results <- list()

basal_asRNA_synth <- params_mean$base_RNA_synth * params_mean$f_asRNA_vals_out[1]
basal_sRNA_synth <- params_mean$base_RNA_synth_sRNA * params_mean$f_sRNA_vals_out[1]

sim_params <- params_mean
sim_params$RNA_synth_t <- params_mean$base_RNA_synth * params_mean$f_RNA_vals_out[1]
sim_params$asRNA_synth_t <- params_mean$base_RNA_synth * params_mean$f_asRNA_vals_out[1]
sim_params$sRNA_synth_t <- params_mean$base_RNA_synth_sRNA * params_mean$f_sRNA_vals_out[1]

y_initial <- c(RNA_free = 0, RNA_complex_sRNA = 0, asRNA_free = 0, sRNA_free = 0)

names(y_initial) <- c("RNA_free", "RNA_complex_sRNA", "asRNA_free", "sRNA_free")
times <- seq(from = 0, to = 30, by = 0.1)

for (current_scenario in scenarios_to_run) {
  cat("Running simulation for scenario:", current_scenario, "\n")
  
  sim_params <- params_mean
  sim_params$rna_interpolator <- rna_interpolator
  sim_params$rna_hetero_interpolator <- rna_hetero_interpolator
  sim_params$asrna_interpolator <- asrna_interpolator
  sim_params$srna_interpolator <- srna_interpolator
  sim_params$scenario <- current_scenario
  
  simulation_output <- deSolve::ode(
    y = y_initial,
    times = times,
    func = rna_dynamic_model_scenarios,
    parms = sim_params
  )
  
  results_df <- as.data.frame(simulation_output)
  results_df$scenario <- current_scenario
  all_scenario_results[[current_scenario]] <- results_df
}

combined_results <- bind_rows(all_scenario_results)

plot_data <- combined_results %>%
  mutate(Total_RNA = RNA_free + RNA_complex_sRNA)

sRNA_comparison_plot <- ggplot(plot_data, aes(x = time, y = Total_RNA, color = scenario)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "nblA concentration simulation",
    x = "Time [h]",
    y = "Total observed nblA concentration",
    color = "Scenario"
  ) +
  theme_bw(base_size = 14) +
  ggsci::scale_color_npg() +
  theme(legend.position = "bottom")


print(sRNA_comparison_plot)
