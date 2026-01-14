library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid) 

# 1. LOAD MODEL FIT

if (!exists("fit")) {
  # Replace with your actual fit object name/path if different
 load("fit.Rdata")
}

# Extract draws
draws_df <- as_draws_df(fit$draws())

# --- DATA CONSTANTS ---
f_RNA_9_data    <- 10.88888889
f_asRNA_9_data  <- 120.7
f_asRNA_24_data <- 96.7

# PARAMETER BOUNDS (from Stan model)
bounds_list <- list(
  Kd_sRNA              = c(10, 5000),
  decay_R_free         = c(0.1, 5),
  decay_complex_sRNA   = c(5, 50),
  decay_S_free         = c(5, 100),
  decay_asR            = c(10, 100),
  base_syn             = c(10, 50000),
  base_syn_S           = c(500, 50000),
  term_rate_asR_on_RNA = c(0.001, 0.1),
  term_rate_RNA_on_asR = c(0, 1),
  synthesis_leak       = c(1, 200),
  f_RNA_9              = c(1, 50),
  f_RNA_24             = c(1, 20),
  f_asRNA_9            = c(50, 300),
  f_asRNA_24           = c(50, 300),
  f_sRNA_9             = c(0, 10),
  f_sRNA_24            = c(0, 10),
  sigma_rnaseq         = c(0, Inf),
  sigma_drnaseq         = c(0, Inf),
  sigma_microarray     = c(0, Inf),
  sigma_qpcr           = c(0, Inf),
  sigma_nb1            = c(0, Inf),
  sigma_decay          = c(0, Inf),
  sigma_reporter       = c(0, Inf)
)

# 2.PRIORS & BOUNDS DATAFRAME

# Initialize empty dataframe
prior_data <- data.frame(
  Parameter = character(),
  Lower_Bound = numeric(),
  Upper_Bound = numeric(),
  Prior_Distribution = character(),
  Prior_Parameters = character(),
  stringsAsFactors = FALSE
)

# Helper function to add rows
add_prior_row <- function(param, lower, upper, dist_name, dist_params) {
  new_row <- data.frame(
    Parameter = param,
    Lower_Bound = lower,
    Upper_Bound = upper,
    Prior_Distribution = dist_name,
    Prior_Parameters = dist_params,
    stringsAsFactors = FALSE
  )
  return(rbind(prior_data, new_row))
}

#Populate Manual Prior Info

# Kinetics
prior_data <- add_prior_row("Kd_sRNA", 10, 5000, "lognormal", "log(600), 0.5")
prior_data <- add_prior_row("decay_R_free", 0.1, 5, "lognormal", "log(1.5), 0.2")
prior_data <- add_prior_row("decay_complex_sRNA", 5, 50, "lognormal", "log(19), 0.2")
prior_data <- add_prior_row("decay_S_free", 5, 100, "lognormal", "log(10), 0.3")
prior_data <- add_prior_row("decay_asR", 10, 100, "lognormal", "log(60), 0.8")

# Synthesis
prior_data <- add_prior_row("base_syn", 10, 50000, "lognormal", "log(25000), 0.5")
prior_data <- add_prior_row("base_syn_S", 500, 50000, "lognormal", "log(15000), 0.5")

# Interference
prior_data <- add_prior_row("term_rate_asR_on_RNA", 0.001, 0.1, "uniform", "min=0.001, max=0.1")
prior_data <- add_prior_row("term_rate_RNA_on_asR", 0, 1, "beta", "shape1=2, shape2=2")

# Leak
prior_data <- add_prior_row("synthesis_leak", 1, 200, "normal (dynamic)", "mean=expected_leak, sd=0.1*expected_leak")

# Forcing Functions
prior_data <- add_prior_row("f_sRNA_9", 0, 10, "lognormal", "log(0.3), 0.2")
prior_data <- add_prior_row("f_sRNA_24", 0, 10, "lognormal", "log(0.8), 0.2")
prior_data <- add_prior_row("f_RNA_24", 1, 20, "lognormal", "log(5.5), 1")

# Dynamic Forcing Priors
prior_data <- add_prior_row("f_RNA_9", 1, 50, "lognormal", paste0("log(", round(f_RNA_9_data, 2), "), sigma_reporter"))
prior_data <- add_prior_row("f_asRNA_9", 50, 300, "lognormal", paste0("log(", round(f_asRNA_9_data, 2), "), sigma_reporter"))
prior_data <- add_prior_row("f_asRNA_24", 50, 300, "lognormal", paste0("log(", round(f_asRNA_24_data, 2), "), sigma_reporter"))

# Sigmas
prior_data <- add_prior_row("sigma_rnaseq", 0, Inf, "exponential", "rate=5")
prior_data <- add_prior_row("sigma_microarray", 0, Inf, "exponential", "rate=20")
prior_data <- add_prior_row("sigma_qpcr", 0, Inf, "exponential", "rate=10")
prior_data <- add_prior_row("sigma_nb1", 0, Inf, "exponential", "rate=5")
prior_data <- add_prior_row("sigma_decay", 0, Inf, "exponential", "rate=5")
prior_data <- add_prior_row("sigma_reporter", 0, Inf, "exponential", "rate=1")



# 3. CALCULATE POSTERIOR STATISTICS 


# Initialize columns for stats
prior_data$Posterior_Mean <- NA
prior_data$Posterior_Median <- NA
prior_data$Posterior_Q5 <- NA
prior_data$Posterior_Q95 <- NA

for (i in 1:nrow(prior_data)) {
  param_name <- prior_data$Parameter[i]
  
  if (param_name %in% names(draws_df)) {
    samples <- draws_df[[param_name]]
    
    prior_data$Posterior_Mean[i]   <- mean(samples)
    prior_data$Posterior_Median[i] <- median(samples)
    prior_data$Posterior_Q5[i]     <- quantile(samples, 0.05)
    prior_data$Posterior_Q95[i]    <- quantile(samples, 0.95)
  }
}

final_table <- prior_data %>%
  select(Parameter, 
         Posterior_Mean, Posterior_Median, Posterior_Q5, Posterior_Q95,
         Lower_Bound, Upper_Bound, 
         Prior_Distribution, Prior_Parameters)

write.csv(final_table, "parameter_posteriors.csv", row.names = FALSE)


# 4. PLOTTING FUNCTION

plot_param_posterior <- function(param_name, 
                                 prior_fun = NULL, 
                                 prior_args = list(), 
                                 x_label = NULL,
                                 xlims = NULL) {
  
  if (!param_name %in% names(draws_df)) {
    warning(paste("Parameter", param_name, "not found. Skipping."))
    return(NULL)
  }
  
  samples <- draws_df[[param_name]]
  median_val <- median(samples)
  
  if (is.null(xlims)) {
    q_range <- quantile(samples, c(0.005, 0.995))
    bw <- diff(q_range)
    xlims <- c(q_range[1] - bw*0.5, q_range[2] + bw*0.5)
    
    if (param_name %in% names(bounds_list)) {
      b <- bounds_list[[param_name]]
      if (abs(q_range[1] - b[1]) < bw) xlims[1] <- min(xlims[1], b[1] - bw*0.1)
      if (abs(q_range[2] - b[2]) < bw) xlims[2] <- max(xlims[2], b[2] + bw*0.1)
    }
    if (mean(samples > 0) == 1) xlims[1] <- max(0, xlims[1])
  }
  
  # --- Plot ---
  p <- ggplot(data.frame(x = samples), aes(x = x)) +
    geom_density(fill = "steelblue", alpha = 0.4, color = NA) +
    geom_vline(xintercept = median_val, color = "black", linewidth = 1) +
    labs(title = (if (is.null(x_label)) param_name else x_label),
         x = "Value", y = "Density") +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold"))
  
  # --- Add Bounds ---
  if (param_name %in% names(bounds_list)) {
    b <- bounds_list[[param_name]]
    view_width <- xlims[2] - xlims[1]
    
    if (b[1] >= (xlims[1] - view_width)) {
      p <- p + geom_vline(xintercept = b[1], linetype = "dotted", size = 1, color = "black")
    }
    if (b[2] <= (xlims[2] + view_width)) {
      p <- p + geom_vline(xintercept = b[2], linetype = "dotted", size = 1, color = "black")
    }
  }
  
  # --- Add Prior ---
  if (!is.null(prior_fun)) {
    x_seq <- seq(xlims[1], xlims[2], length.out = 300)
    y_prior <- do.call(prior_fun, c(list(x = x_seq), prior_args))
    prior_df <- data.frame(x = x_seq, y = y_prior)
    
    p <- p + geom_line(data = prior_df, aes(x = x, y = y), 
                       color = "red", linetype = "dashed", linewidth = 1)
  }
  
  p <- p + coord_cartesian(xlim = xlims)
  return(p)
}



plot_list <- list()
idx <- 1

sigma_rep_med <- median(draws_df$sigma_reporter)

# --- 1. KINETICS ---
plot_list[[idx]] <- plot_param_posterior("Kd_sRNA", dlnorm, list(meanlog = log(600), sdlog = 0.5), "Kd (sRNA)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("decay_R_free", dlnorm, list(meanlog = log(1.5), sdlog = 0.2), "Decay (Free RNA)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("decay_complex_sRNA", dlnorm, list(meanlog = log(19), sdlog = 0.2), "Decay (Complex)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("decay_S_free", dlnorm, list(meanlog = log(10), sdlog = 0.3), "Decay (Free sRNA)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("decay_asR", dlnorm, list(meanlog = log(60), sdlog = 0.8), "Decay (asRNA)"); idx <- idx + 1

# --- 2. SYNTHESIS ---
plot_list[[idx]] <- plot_param_posterior("base_syn", dlnorm, list(meanlog = log(25000), sdlog = 0.5), "Basal Syn (RNA)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("base_syn_S", dlnorm, list(meanlog = log(15000), sdlog = 0.5), "Basal Syn (sRNA)"); idx <- idx + 1

# --- 3. INTERFERENCE ---
plot_list[[idx]] <- plot_param_posterior("term_rate_asR_on_RNA", dunif, list(min = 0.001, max = 0.1), "TI Rate (asR on RNA)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("term_rate_RNA_on_asR", dbeta, list(shape1 = 2, shape2 = 2), "TI Rate (RNA on asR)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("synthesis_leak", NULL, list(), "Synthesis Leak (Dynamic Prior)"); idx <- idx + 1

# --- 4. FORCING FUNCTIONS ---
plot_list[[idx]] <- plot_param_posterior("f_sRNA_9", dlnorm, list(meanlog = log(0.3), sdlog = 0.2), "f sRNA (t9)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("f_sRNA_24", dlnorm, list(meanlog = log(0.8), sdlog = 0.2), "f sRNA (t24)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("f_RNA_24", dlnorm, list(meanlog = log(5.5), sdlog = 1), "f RNA (t24)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("f_RNA_9", dlnorm, list(meanlog = log(f_RNA_9_data), sdlog = sigma_rep_med), "f RNA (t9)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("f_asRNA_9", dlnorm, list(meanlog = log(f_asRNA_9_data), sdlog = sigma_rep_med), "f asRNA (t9)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("f_asRNA_24", dlnorm, list(meanlog = log(f_asRNA_24_data), sdlog = sigma_rep_med), "f asRNA (t24)"); idx <- idx + 1

# --- 5. SIGMAS ---
plot_list[[idx]] <- plot_param_posterior("sigma_rnaseq", dexp, list(rate = 5), "Sigma (RNA-seq)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("sigma_microarray", dexp, list(rate = 20), "Sigma (Microarray)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("sigma_qpcr", dexp, list(rate = 10), "Sigma (qPCR)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("sigma_nb1", dexp, list(rate = 5), "Sigma (NB1)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("sigma_decay", dexp, list(rate = 5), "Sigma (Decay)"); idx <- idx + 1
plot_list[[idx]] <- plot_param_posterior("sigma_reporter", dexp, list(rate = 1), "Sigma (Reporter)"); idx <- idx + 1

plot_list <- plot_list[!sapply(plot_list, is.null)]

final_plot <- grid.arrange(grobs = plot_list, ncol = 5, 
                           top = textGrob("Posterior (Blue) vs Prior (Red Dashed) & Bounds (Black Dotted)", 
                                          gp = gpar(fontsize = 15, fontface = "bold")))

ggsave("parameter_posteriors.pdf", plot = final_plot, width = 20, height = 15, dpi = 300)
