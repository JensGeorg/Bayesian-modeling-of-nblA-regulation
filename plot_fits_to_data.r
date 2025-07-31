library(ggplot2)
library(dplyr)
library(rstan)
library(tidyr)

obs_times <- c(0, 6, 8, 12, 24)
obs_RNA <- c(1, 11.23555901, 12.55334557, 14.32040113, 6.453134074)
obs_asRNA <- c(1, 1.180992661, 0.590496331, 0.726986259, 0.946057647)
obs_sRNA <- c(1, 0.271683716, 0.332171454, 0.320856474, 0.77916458)
n_obs <- length(obs_times)

decay_obs_times <- c(0, 24)
obs_RNA_decay <- c(13.73, 4.4)
n_decay_obs <- length(decay_obs_times)

nb_times_srna_ko <- c(0, 2, 4, 6, 8)
obs_RNA_nb_srna_ko <- c(1, 5.76018571, 9.76801526, 38.1105133, 46.8838844)
n_obs_nb_srna_ko <- length(nb_times_srna_ko)

nb_times_wt1 <- c(0, 3, 6, 9, 12, 24)
obs_RNA_nb_wt1 <- c(1, 9, 14, 19, 18, 7)
n_obs_nb_wt1 <- length(nb_times_wt1)

nb_times_wt2 <- c(0, 8, 24)
obs_RNA_nb_wt2 <- c(1.0615e-11, 0.80544218, 1.30858618)
n_obs_nb_wt2 <- length(nb_times_wt2)

nb_times_asrna_ko <- c(0, 8, 24)
obs_RNA_nb_asrna_ko <- c(1.82039604, 4.46959099, 3.34963486)
n_obs_nb_asrna_ko <- length(nb_times_asrna_ko)

k_off_lit_estimate <- 18000
k_off_lower_bound <- 1080
k_off_upper_bound <- 36000

k_off_log_mean_lit_val <- log(k_off_lit_estimate)
k_off_log_sd_lit_val <- (log(k_off_upper_bound) - log(k_off_lower_bound)) / (2 * qnorm(0.975))

t0 <- 0

y0_data <- c(10.96, 3, 117.91, 62.84)

fix_y0 <- 0
sRNA_strategy <- 2

rna_ts <- c(0, 9, 24)
rna_vals_data <- c(1, 10.88888889, 15.77777778)
n_RNA <- length(rna_ts)

asrna_ts <- c(0, 9, 24)
asrna_vals_data <- c(1, 1.180992661, 0.946057647) * 102.2222222
n_asRNA <- length(asrna_ts)

srna_ts <- c(0, 6, 8, 12, 24)
n_sRNA <- length(srna_ts)

L <- n_RNA + n_asRNA + n_sRNA

kd_estimate_nls <- 117.4665
ci_lower_nls <- 96.13608
ci_upper_nls <- 143.16898

z_score_95 <- qnorm(0.975)

mu_log_Kd_prior_val <- mean(c(log(ci_lower_nls), log(ci_upper_nls)))
sigma_log_Kd_prior_val <- (log(ci_upper_nls) - log(ci_lower_nls)) / (2 * z_score_95)


stan_data <- list(
  n_obs = n_obs,
  obs_times = obs_times,
  obs_RNA = obs_RNA,
  obs_asRNA = obs_asRNA,
  obs_sRNA = obs_sRNA,
  n_decay_obs = n_decay_obs,
  decay_obs_times = decay_obs_times,
  obs_RNA_decay = obs_RNA_decay,
  n_obs_nb_wt1 = n_obs_nb_wt1,
  nb_times_wt1 = nb_times_wt1,
  obs_RNA_nb_wt1 = obs_RNA_nb_wt1,
  n_obs_nb_srna_ko = n_obs_nb_srna_ko,
  nb_times_srna_ko = nb_times_srna_ko,
  obs_RNA_nb_srna_ko = obs_RNA_nb_srna_ko,
  n_obs_nb_wt2 = n_obs_nb_wt2,
  nb_times_wt2 = nb_times_wt2,
  obs_RNA_nb_wt2 = obs_RNA_nb_wt2,
  n_obs_nb_asrna_ko = n_obs_nb_asrna_ko,
  nb_times_asrna_ko = nb_times_asrna_ko,
  obs_RNA_nb_asrna_ko = obs_RNA_nb_asrna_ko,
  fix_y0 = fix_y0,
  y0_data = y0_data,
  t0 = t0,
  sRNA_strategy = sRNA_strategy,
  n_RNA = n_RNA,
  rna_ts = rna_ts,
  n_asRNA = n_asRNA,
  asrna_ts = asrna_ts,
  n_sRNA = n_sRNA,
  srna_ts = srna_ts,
  L = L,
  rna_vals_data = rna_vals_data,
  asrna_vals_data = asrna_vals_data,
  kd_srna_log_mean_prior = mu_log_Kd_prior_val,
  kd_srna_log_sd_prior = sigma_log_Kd_prior_val,
  k_off_log_mean_lit = k_off_log_mean_lit_val,
  k_off_log_sd_lit = k_off_log_sd_lit_val
)

posterior_draws <- as.matrix(fit)

get_stats <- function(draws) {
  median_val <- median(draws)
  lower_ci <- quantile(draws, 0.025)
  upper_ci <- quantile(draws, 0.975)
  return(c(median = median_val, lower = lower_ci, upper = upper_ci))
}

preds_list <- list()

rna_pred_samples <- posterior_draws[, grep("^RNA_pred_out", colnames(posterior_draws)), drop = FALSE]
rna_preds_df <- as.data.frame(t(apply(rna_pred_samples, 2, get_stats)))
rna_preds_df$time <- obs_times
rna_preds_df$type <- "Microarray RNA"
preds_list[[length(preds_list) + 1]] <- rna_preds_df

asrna_pred_samples <- posterior_draws[, grep("^asRNA_pred_out", colnames(posterior_draws)), drop = FALSE]
asrna_preds_df <- as.data.frame(t(apply(asrna_pred_samples, 2, get_stats)))
asrna_preds_df$time <- obs_times
asrna_preds_df$type <- "Microarray asRNA"
preds_list[[length(preds_list) + 1]] <- asrna_preds_df

sRNA_pred_samples <- posterior_draws[, grep("^sRNA_pred_out", colnames(posterior_draws)), drop = FALSE]
sRNA_preds_df <- as.data.frame(t(apply(sRNA_pred_samples, 2, get_stats)))
sRNA_preds_df$time <- obs_times
sRNA_preds_df$type <- "Microarray sRNA"
preds_list[[length(preds_list) + 1]] <- sRNA_preds_df

decay_pred_samples <- posterior_draws[, grep("^RNA_decay_pred_out", colnames(posterior_draws)), drop = FALSE]
decay_preds_df <- as.data.frame(t(apply(decay_pred_samples, 2, get_stats)))
decay_preds_df$time <- decay_obs_times
decay_preds_df$type <- "RNA Decay Rate"
decay_preds_df <- decay_preds_df[decay_preds_df$median >= 0, ]
preds_list[[length(preds_list) + 1]] <- decay_preds_df

nb_wt1_pred_samples <- posterior_draws[, grep("^RNA_pred_nb_wt1_out", colnames(posterior_draws)), drop = FALSE]
nb_wt1_preds_df <- as.data.frame(t(apply(nb_wt1_pred_samples, 2, get_stats)))
nb_wt1_preds_df$time <- nb_times_wt1
nb_wt1_preds_df$type <- "Northern Blot WT1"
preds_list[[length(preds_list) + 1]] <- nb_wt1_preds_df

nb_srna_ko_pred_samples <- posterior_draws[, grep("^RNA_pred_nb_srna_ko_out", colnames(posterior_draws)), drop = FALSE]
nb_srna_ko_preds_df <- as.data.frame(t(apply(nb_srna_ko_pred_samples, 2, get_stats)))
nb_srna_ko_preds_df$time <- nb_times_srna_ko
nb_srna_ko_preds_df$type <- "Northern Blot sRNA KO"
preds_list[[length(preds_list) + 1]] <- nb_srna_ko_preds_df

nb_wt2_pred_samples <- posterior_draws[, grep("^RNA_pred_nb_wt2_out", colnames(posterior_draws)), drop = FALSE]
nb_wt2_preds_df <- as.data.frame(t(apply(nb_wt2_pred_samples, 2, get_stats)))
nb_wt2_preds_df$time <- nb_times_wt2
nb_wt2_preds_df$type <- "Northern Blot WT2"
preds_list[[length(preds_list) + 1]] <- nb_wt2_preds_df

nb_asrna_ko_pred_samples <- posterior_draws[, grep("^RNA_pred_nb_asrna_ko_out", colnames(posterior_draws)), drop = FALSE]
nb_asrna_ko_preds_df <- as.data.frame(t(apply(nb_asrna_ko_pred_samples, 2, get_stats)))
nb_asrna_ko_preds_df$time <- nb_times_asrna_ko
nb_asrna_ko_preds_df$type <- "Northern Blot asRNA KO"
preds_list[[length(preds_list) + 1]] <- nb_asrna_ko_preds_df

all_predictions_df <- do.call(rbind, preds_list)
colnames(all_predictions_df) <- c("median", "lower", "upper", "time", "type")

df_obs_microarray <- tibble(
  time = rep(obs_times, 3),
  type = c(rep("Microarray RNA", n_obs),
           rep("Microarray asRNA", n_obs),
           rep("Microarray sRNA", n_obs)),
  observed_value = c(obs_RNA, obs_asRNA, obs_sRNA)
)

df_obs_decay <- tibble(
  time = decay_obs_times,
  type = "RNA Decay Rate",
  observed_value = obs_RNA_decay
)

df_obs_nb_wt1 <- tibble(
  time = nb_times_wt1,
  type = "Northern Blot WT1",
  observed_value = obs_RNA_nb_wt1
)

df_obs_nb_srna_ko <- tibble(
  time = nb_times_srna_ko,
  type = "Northern Blot sRNA KO",
  observed_value = obs_RNA_nb_srna_ko
)

df_obs_nb_wt2 <- tibble(
  time = nb_times_wt2,
  type = "Northern Blot WT2",
  observed_value = obs_RNA_nb_wt2
)

df_obs_nb_asrna_ko <- tibble(
  time = nb_times_asrna_ko,
  type = "Northern Blot asRNA KO",
  observed_value = obs_RNA_nb_asrna_ko
)

all_observations <- bind_rows(
  df_obs_microarray,
  df_obs_decay,
  df_obs_nb_wt1,
  df_obs_nb_srna_ko,
  df_obs_nb_wt2,
  df_obs_nb_asrna_ko
)

facet_limits <- all_predictions_df %>%
  group_by(type) %>%
  summarise(
    min_val = 0,
    max_val = max(2, max(upper, na.rm = TRUE))
  )

dummy_data_for_limits <- facet_limits %>%
  select(type, min_val, max_val) %>%
  tidyr::pivot_longer(cols = c(min_val, max_val), names_to = "limit_type", values_to = "dummy_y") %>%
  mutate(time = 0)

plot_timeseries_custom_range <- ggplot() +
  geom_blank(data = dummy_data_for_limits, aes(x = time, y = dummy_y, group = type)) +
  geom_ribbon(data = all_predictions_df, aes(x = time, ymin = lower, ymax = upper, fill = type),
              alpha = 0.2, show.legend = FALSE) +
  geom_line(data = all_predictions_df, aes(x = time, y = median, color = type), linewidth = 1) +
  geom_point(data = all_observations, aes(x = time, y = observed_value, color = type), size = 2) +
  facet_wrap(~ type, scales = "free_y", ncol = 3) +
  labs(title = "Model Predictions vs. Observed Data with 95% Credible Intervals",
       x = "Time (hours)",
       y = "Value") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  scale_color_brewer(palette = "Set1")


print(plot_timeseries_custom_range)
