
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)


# 1. Fit Data

if (!exists("fit")) {
  # Replace with your actual fit object name/path if different
  fit <- load("fit.Rdata")
}
write.table(fit$summary(), file="out.txt", sep="\t")



#  A. Microarray Data
# t=0 is excluded. Scaled to the model at t=6.
obs_data_microarray <- data.frame(
  Time = c(6, 8, 12, 24),
  RNA_rel = c(11.24, 12.55, 14.32, 6.45),
  sRNA_rel = c(0.272, 0.332, 0.321, 0.779),
  asRNA_rel = c(1.181, 0.590, 0.727, 0.946)
)

# B. RNA-seq Data (Fold-change relative to t=0)
# Derived from raw counts 
obs_data_rnaseq <- data.frame(
  Time = c(0, 9),
  RNA_rel = c(1, 46.38),
  sRNA_rel = c(1, 0.24),
  asRNA_rel = c(1, 0.88)
)

# C. Additional Constraints  ---
# NB1 sRNA KO: Raw quantification
data_nb1_srna_ko <- data.frame(
  Time = c(0, 2, 4, 6, 8),
  Value_Raw = c(1, 5.76, 9.77, 38.1, 46.9),
  Species = "sRNA KO (NB1)"
)

# NB1 WT: Raw quantification
data_nb1_wt <- data.frame(
  Time = c(0, 3, 6, 9, 12, 24),
  Value_Raw = c(1, 9, 14, 19, 18, 7),
  Species = "WT (NB1)"
)

# qPCR asRNA KO: Raw quantification (t0 vs t9)
data_qpcr_asrna_ko <- data.frame(
  Time = c(0, 9),
  Value_Raw = c(1.82, 4.47),
  Species = "asRNA KO (qPCR)"
)

# Decay: Absolute rates
data_decay <- data.frame(
  Time = c(0, 9),
  Value_Raw = c(13.73, 4.4),
  Species = "Effective Decay Rate"
)

# Model Time Vector
model_times <- c(0.2, 0.5, 1, 2, 3, 4, 6, 8, 9, 12, 16, 24)


# 2. HELPER FUNCTIONS

out <- read.table("out.txt", header = TRUE, stringsAsFactors = FALSE)

extract_trajectory <- function(var_pattern, time_vec, t0_var_name) {
  t0_row <- out[out$variable == t0_var_name, ]
  t0_df <- data.frame(Time = 0, Mean = t0_row$mean, Lower = t0_row$q5, Upper = t0_row$q95)
  rows <- out[grep(paste0("^", var_pattern, "\\["), out$variable), ]
  rows$idx <- as.numeric(gsub(paste0("^", var_pattern, "\\[([0-9]+)\\]"), "\\1", rows$variable))
  rows <- rows[order(rows$idx), ]
  rows <- rows[1:length(time_vec), ]
  traj_df <- data.frame(Time = time_vec, Mean = rows$mean, Lower = rows$q5, Upper = rows$q95)
  return(rbind(t0_df, traj_df))
}

# Normalize model to 1 at t=0
normalize_to_t0 <- function(df) {
  ref_mean <- df$Mean[df$Time == 0]
  df$Mean  <- df$Mean / ref_mean
  df$Lower <- df$Lower / ref_mean
  df$Upper <- df$Upper / ref_mean
  return(df)
}

# 3. PROCESS MODEL TRAJECTORIES

# Plot 1: Main Species (Normalized to t=0)
# RNA
df_rna <- extract_trajectory("RNA_wt_trajectory", model_times, "RNA_wt_t0")
df_rna <- normalize_to_t0(df_rna); df_rna$Species <- "nblA mRNA"

# sRNA
sRNA_rows <- out[grep("y_hat\\[[0-9]+,2\\]", out$variable), ]
sRNA_rows$idx <- as.numeric(gsub("y_hat\\[([0-9]+),2\\]", "\\1", sRNA_rows$variable))
sRNA_rows <- sRNA_rows[order(sRNA_rows$idx), ]
t0_sRNA <- out[out$variable == "S0_wt", ]
df_srna <- rbind(data.frame(Time=0, Mean=t0_sRNA$mean, Lower=t0_sRNA$q5, Upper=t0_sRNA$q95),
                 data.frame(Time=model_times, Mean=sRNA_rows$mean, Lower=sRNA_rows$q5, Upper=sRNA_rows$q95))
df_srna <- normalize_to_t0(df_srna); df_srna$Species <- "NsrR1"

# asRNA
asRNA_rows <- out[grep("y_hat\\[[0-9]+,3\\]", out$variable), ]
asRNA_rows$idx <- as.numeric(gsub("y_hat\\[([0-9]+),3\\]", "\\1", asRNA_rows$variable))
asRNA_rows <- asRNA_rows[order(asRNA_rows$idx), ]
t0_asRNA <- out[out$variable == "asR0_tot", ]
if(nrow(t0_asRNA)==0) t0_asRNA <- out[out$variable == "y0_wt[3]", ]
df_asrna <- rbind(data.frame(Time=0, Mean=t0_asRNA$mean, Lower=t0_asRNA$q5, Upper=t0_asRNA$q95),
                  data.frame(Time=model_times, Mean=asRNA_rows$mean, Lower=asRNA_rows$q5, Upper=asRNA_rows$q95))
df_asrna <- normalize_to_t0(df_asrna); df_asrna$Species <- "as_nblA"

df_main_model <- rbind(df_asrna, df_rna, df_srna)

# Plot 2: Validation Constraints (Absolute or Scaled) ---
df_srna_ko <- extract_trajectory("RNA_srna_ko_trajectory", model_times, "RNA_srna_ko_t0")
df_srna_ko$Species <- "sRNA KO (NB1)"

df_wt_nb1 <- extract_trajectory("RNA_wt_trajectory", model_times, "RNA_wt_t0")
df_wt_nb1$Species <- "WT (NB1)"

df_asrna_ko <- extract_trajectory("RNA_asrna_ko_trajectory", model_times, "RNA_asrna_ko_t0")
df_asrna_ko$Species <- "asRNA KO (qPCR)"

df_decay <- extract_trajectory("decay_eff_wt", model_times, "decay_eff_wt_t0")
df_decay$Species <- "Effective Decay Rate"

df_val_model <- rbind(df_srna_ko, df_wt_nb1, df_asrna_ko, df_decay)


# 4. SCALE DATA TO MATCH MODEL

# Helper to scale data based on a specific time anchor
scale_to_model <- function(obs_df, model_df, anchor_time, val_col, species_name) {
  # Get Model value at anchor
  if(anchor_time == 0) {
    model_val <- model_df$Mean[model_df$Time == 0 & model_df$Species == species_name]
  } else {
    traj <- model_df[model_df$Species == species_name, ]
    model_val <- approx(traj$Time, traj$Mean, xout = anchor_time)$y
  }
  
  # Get Data value at anchor
  data_val <- obs_df[[val_col]][obs_df$Time == anchor_time]
  
  # Calculate factor
  k <- model_val / data_val
  
  obs_df$Value_Scaled <- obs_df[[val_col]] * k
  return(obs_df)
}

#Plot 1 Scaling 
# 1. Microarray: Scale to Model at t=6 (First reliable point)
obs_micro_rna   <- data.frame(Time=obs_data_microarray$Time, Val=obs_data_microarray$RNA_rel) %>%
  scale_to_model(df_main_model, 6, "Val", "nblA mRNA")
obs_micro_srna  <- data.frame(Time=obs_data_microarray$Time, Val=obs_data_microarray$sRNA_rel) %>%
  scale_to_model(df_main_model, 6, "Val", "NsrR1")
obs_micro_asrna <- data.frame(Time=obs_data_microarray$Time, Val=obs_data_microarray$asRNA_rel) %>%
  scale_to_model(df_main_model, 6, "Val", "as_nblA")

obs_micro_combined <- rbind(
  data.frame(Time=obs_micro_rna$Time, Value=obs_micro_rna$Value_Scaled, Species="nblA mRNA"),
  data.frame(Time=obs_micro_srna$Time, Value=obs_micro_srna$Value_Scaled, Species="NsrR1"),
  data.frame(Time=obs_micro_asrna$Time, Value=obs_micro_asrna$Value_Scaled, Species="as_nblA")
)
obs_micro_combined$DataType <- "Microarray"

# 2. RNA-seq: No scaling needed (Already relative to t=0, model is relative to t=0)
obs_rnaseq_combined <- data.frame(
  Time = rep(obs_data_rnaseq$Time, 3),
  Value = c(obs_data_rnaseq$asRNA_rel, obs_data_rnaseq$RNA_rel, obs_data_rnaseq$sRNA_rel),
  Species = rep(c("as_nblA", "nblA mRNA", "NsrR1"), each=2),
  DataType = "RNA-seq"
)

obs_p1_final <- rbind(obs_micro_combined, obs_rnaseq_combined)

# Plot 2 Scaling 
# 1. NB1 sRNA KO: Scale to t=2 (Constraint definition)
nb1_srna_scaled <- scale_to_model(data_nb1_srna_ko, df_val_model, 2, "Value_Raw", "sRNA KO (NB1)")

# 2. NB1 WT: Scale to t=3 (Constraint definition)
nb1_wt_scaled <- scale_to_model(data_nb1_wt, df_val_model, 3, "Value_Raw", "WT (NB1)")

# 3. qPCR: Scale to t=0 (Fold change constraint)
qpcr_scaled <- scale_to_model(data_qpcr_asrna_ko, df_val_model, 0, "Value_Raw", "asRNA KO (qPCR)")

# 4. Decay: No scaling
decay_data <- data_decay
decay_data$Value_Scaled <- decay_data$Value_Raw

# Combine and FILTER t=0 for NB1
obs_p2_final <- rbind(
  data.frame(Time=nb1_srna_scaled$Time, Value=nb1_srna_scaled$Value_Scaled, Species="sRNA KO (NB1)"),
  data.frame(Time=nb1_wt_scaled$Time, Value=nb1_wt_scaled$Value_Scaled, Species="WT (NB1)"),
  data.frame(Time=qpcr_scaled$Time, Value=qpcr_scaled$Value_Scaled, Species="asRNA KO (qPCR)"),
  data.frame(Time=decay_data$Time, Value=decay_data$Value_Scaled, Species="Effective Decay Rate")
)

obs_p2_final <- obs_p2_final %>%
  filter(!(Time == 0 & Species %in% c("sRNA KO (NB1)", "WT (NB1)")))


# 5. PLOTTING

# Factor Levels
levels_p1 <- c("as_nblA", "nblA mRNA", "NsrR1")
df_main_model$Species <- factor(df_main_model$Species, levels = levels_p1)
obs_p1_final$Species <- factor(obs_p1_final$Species, levels = levels_p1)

levels_p2 <- c("sRNA KO (NB1)", "WT (NB1)", "asRNA KO (qPCR)", "Effective Decay Rate")
df_val_model$Species <- factor(df_val_model$Species, levels = levels_p2)
obs_p2_final$Species <- factor(obs_p2_final$Species, levels = levels_p2)

theme_clean <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "top"
  )

# --- Plot 1 ---
p1 <- ggplot() +
  geom_ribbon(data = df_main_model, aes(x = Time, ymin = Lower, ymax = Upper), fill = "#ffcccc", alpha = 0.6) +
  geom_line(data = df_main_model, aes(x = Time, y = Mean), size = 1, color = "black") +
  geom_point(data = filter(obs_p1_final, DataType == "Microarray"), aes(x = Time, y = Value, shape = DataType), size = 2.5) +
  geom_point(data = filter(obs_p1_final, DataType == "RNA-seq"), aes(x = Time, y = Value, shape = DataType), size = 3, color = "red", fill = "red") +
  scale_shape_manual(values = c("Microarray" = 16, "RNA-seq" = 15)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25), limits = c(0, 25)) +  
  facet_wrap(~Species, ncol = 1, scales = "free_y") +
  expand_limits(y = 0) + 
  labs(x = expression(bold(paste("Time after ", NH[4]^"+", " removal (h)"))), y = "Relative concentration (Norm. to t=0)") +
  theme_clean

# --- Plot 2 ---
p2 <- ggplot() +
  geom_ribbon(data = df_val_model, aes(x = Time, ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.2) +
  geom_line(data = df_val_model, aes(x = Time, y = Mean), color = "steelblue", size = 1) +
  geom_point(data = obs_p2_final, aes(x = Time, y = Value), color = "black", size = 2.5) +
  facet_wrap(~Species, scales = "free_y", ncol = 2) +
  expand_limits(y = 0) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24), limits = c(0, 24.5)) +
  labs(x = "Time (h)", y = "Concentration / Rate (Scaled)") +
  theme_clean

# Save plots
# ggsave("trajectories_part1.pdf", plot = p1, width = 4, height = 8, dpi = 300)
# ggsave("trajectories_part2.pdf", plot = p2, width = 8, height = 6, dpi = 300)
# ggsave("trajectories_part1.png", plot = p1, width = 4, height = 8, dpi = 300)
# ggsave("trajectories_part2.png", plot = p2, width = 8, height = 6, dpi = 300)
