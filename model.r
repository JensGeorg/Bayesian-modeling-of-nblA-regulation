



options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_model_code <- "
functions {
  real linear_interp(real t, int n, real[] ts, real[] ys) {
    if (n == 0)
      return 0.0;
    if (t <= ts[1])
      return ys[1];
    if (t >= ts[n])
      return ys[n];
    for (i in 1:(n - 1)) {
      if (t >= ts[i] && t < ts[i + 1]) {
        if (ts[i + 1] == ts[i])
          return ys[i];
        real slope = (ys[i + 1] - ys[i]) / (ts[i + 1] - ts[i]);
        return ys[i] + slope * (t - ts[i]);
      }
    }
    return ys[n];
  }

  real[] rna_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    int K_plus_1 = size(theta);
    int K = K_plus_1 - 1;
    real flag_real = (K > 0 && K_plus_1 <= size(theta)) ? theta[K_plus_1] : 0.0;
    int condition_flag = 0;
    if (fabs(flag_real - 1.0) < 1e-9) {
      condition_flag = 1;
    } else if (fabs(flag_real - 2.0) < 1e-9) {
      condition_flag = 2;
    }

    int n_RNA = x_i[1];
    int n_asRNA = x_i[2];
    int n_sRNA = x_i[3];

    real base_RNA_synth = K >= 1 ? theta[1] : 0.0;
    real RNA_decay = K >= 2 ? theta[2] : 0.0;
    real k_on_sRNA = K >= 3 ? theta[3] : 0.0;
    real k_off_sRNA = K >= 4 ? theta[4] : 0.0;
    real termination_rate_RNA = K >= 5 ? theta[5] : 0.0;
    real termination_rate_aRNA = K >= 6 ? theta[6] : 0.0;
    real base_RNA_synth_sRNA = K >= 7 ? theta[7] : 0.0;
    real sRNA_decay = K >= 8 ? theta[8] : 0.0;
    real asRNA_decay = K >= 9 ? theta[9] : 0.0;
    real complex_decay_sRNA = K >= 10 ? theta[10] : 0.0;

    int offset_asRNA_ts = n_RNA;
    int offset_sRNA_ts = n_RNA + n_asRNA;
    int offset_f_RNA = 10;
    int offset_f_asRNA = 10 + n_RNA;
    int offset_f_sRNA = 10 + n_RNA + n_asRNA;

    real f_RNA = 0.0;
    real f_asRNA = 0.0;
    real f_sRNA = 0.0;

    if (n_RNA > 0) {
      int start_idx = offset_f_RNA + 1;
      int end_idx = offset_f_RNA + n_RNA;
      if (start_idx <= K && end_idx <= K && start_idx <= end_idx) {
        f_RNA = linear_interp(t, n_RNA, x_r[1 : n_RNA], theta[start_idx : end_idx]);
      }
    }

    if (n_asRNA > 0) {
      int start_idx = offset_f_asRNA + 1;
      int end_idx = offset_f_asRNA + n_asRNA;
      if (start_idx <= K && end_idx <= K && start_idx <= end_idx) {
        f_asRNA = linear_interp(t, n_asRNA, x_r[(offset_asRNA_ts + 1) : (offset_asRNA_ts + n_asRNA)], theta[start_idx : end_idx]);
      }
    }

    if (n_sRNA > 0) {
      int start_idx = offset_f_sRNA + 1;
      int end_idx = offset_f_sRNA + n_sRNA;
      if (start_idx <= K && end_idx <= K && start_idx <= end_idx) {
        f_sRNA = linear_interp(t, n_sRNA, x_r[(offset_sRNA_ts + 1) : (offset_sRNA_ts + n_sRNA)], theta[start_idx : end_idx]);
      }
    }

    real RNA_synth_t = base_RNA_synth * f_RNA;
    real asRNA_synth_t = base_RNA_synth * f_asRNA;
    real sRNA_synth_t = base_RNA_synth_sRNA * f_sRNA;

    if (condition_flag == 1) {
      sRNA_synth_t = 0.0;
    }
    if (condition_flag == 2) {
      asRNA_synth_t = 0.0;
    }

    real effective_RNA_synth = fmax(0, RNA_synth_t - asRNA_synth_t * termination_rate_RNA);
    real effective_asRNA_synth = fmax(0, asRNA_synth_t - RNA_synth_t * termination_rate_aRNA);

    real dydt[4];
    dydt[1] = effective_RNA_synth - k_on_sRNA * y[4] * y[1] - RNA_decay * y[1] + k_off_sRNA * y[2];
    dydt[2] = k_on_sRNA * y[4] * y[1] - k_off_sRNA * y[2] - complex_decay_sRNA * y[2];
    dydt[3] = effective_asRNA_synth - asRNA_decay * y[3];
    dydt[4] = sRNA_synth_t - sRNA_decay * y[4] - k_on_sRNA * y[4] * y[1] + k_off_sRNA * y[2];
    return dydt;
  }
}

data {
  int<lower=0> n_obs;
  real obs_times[n_obs];
  real obs_RNA[n_obs];
  real obs_asRNA[n_obs];
  real obs_sRNA[n_obs];
  int<lower=0> n_decay_obs;
  real decay_obs_times[n_decay_obs];
  real obs_RNA_decay[n_decay_obs];
  int<lower=0> n_obs_nb_wt1;
  real nb_times_wt1[n_obs_nb_wt1];
  real obs_RNA_nb_wt1[n_obs_nb_wt1];
  int<lower=0> n_obs_nb_srna_ko;
  real nb_times_srna_ko[n_obs_nb_srna_ko];
  real obs_RNA_nb_srna_ko[n_obs_nb_srna_ko];
  int<lower=0> n_obs_nb_wt2;
  real nb_times_wt2[n_obs_nb_wt2];
  real obs_RNA_nb_wt2[n_obs_nb_wt2];
  int<lower=0> n_obs_nb_asrna_ko;
  real nb_times_asrna_ko[n_obs_nb_asrna_ko];
  real obs_RNA_nb_asrna_ko[n_obs_nb_asrna_ko];
  real<lower=0> t0;
  int<lower=0, upper=1> fix_y0;
  real y0_data[4];
  int<lower=0> n_RNA;
  real rna_ts[n_RNA];
  int<lower=0> n_asRNA;
  real asrna_ts[n_asRNA];
  int<lower=0> n_sRNA;
  real srna_ts[n_sRNA];
  int<lower=0> L;
  real rna_vals_data[n_RNA];
  real asrna_vals_data[n_asRNA];
  int<lower=1, upper=2> sRNA_strategy;
  real kd_srna_log_mean_prior;
  real<lower=0> kd_srna_log_sd_prior;
  real k_off_log_mean_lit;
  real<lower=0> k_off_log_sd_lit;
}

transformed data {
  int x_i[3];
  x_i[1] = n_RNA;
  x_i[2] = n_asRNA;
  x_i[3] = n_sRNA;
  int tot_ts = n_RNA + n_asRNA + n_sRNA;
  real x_r[tot_ts];
  int current_pos = 1;

  if (n_RNA > 0) {
    for (i in 1:n_RNA) {
      if (current_pos <= tot_ts)
        x_r[current_pos] = rna_ts[i];
      current_pos += 1;
    }
  }

  if (n_asRNA > 0) {
    for (i in 1:n_asRNA) {
      if (current_pos <= tot_ts)
        x_r[current_pos] = asrna_ts[i];
      current_pos += 1;
    }
  }

  if (n_sRNA > 0) {
    for (i in 1:n_sRNA) {
      if (current_pos <= tot_ts)
        x_r[current_pos] = srna_ts[i];
      current_pos += 1;
    }
  }

  int decay_obs_indices[n_decay_obs];
  if (n_decay_obs > 0 && n_obs > 0) {
    for (i in 1:n_decay_obs) {
      int found_idx = 0;
      for (j in 1:n_obs) {
        if (fabs(obs_times[j] - decay_obs_times[i]) < 1e-9) {
          found_idx = j;
          break;
        }
      }
      decay_obs_indices[i] = found_idx;
    }
  }

  int max_steps = 100000;
  real rel_tol = 1e-5;
  real abs_tol = 1e-5;
  int K = 10 + L;

  int N_total_obs = 0;
  if (n_obs > 0) {
    N_total_obs += n_obs;
    N_total_obs += n_obs;
    N_total_obs += n_obs;
  }
  if (n_obs_nb_wt1 > 0) {
    N_total_obs += 1;
    if (n_obs_nb_wt1 > 1) {
      N_total_obs += (n_obs_nb_wt1 - 1);
    }
  }
  if (n_obs_nb_srna_ko > 0) {
    N_total_obs += 1;
    if (n_obs_nb_srna_ko > 1) {
      N_total_obs += (n_obs_nb_srna_ko - 1);
    }
  }
  if (n_obs_nb_wt2 > 0) {
    N_total_obs += 1;
    if (n_obs_nb_wt2 > 1) {
      N_total_obs += (n_obs_nb_wt2 - 1);
    }
  }
  if (n_obs_nb_asrna_ko > 0) {
    N_total_obs += 1;
    if (n_obs_nb_asrna_ko > 1) {
      N_total_obs += (n_obs_nb_asrna_ko - 1);
    }
  }
}

parameters {
  real<lower=0> param_log_RNA_free0;
  real<lower=0> param_log_sRNA_free0;
  real log_k_off_sRNA;
  real log_RNA_free0_sKO;
  real log_RNA_free0_asKO;
  real log_asRNA_free0_param;
  real<lower=0> base_RNA_synth;
  real<lower=0> RNA_decay;
  real<lower=0> Kd_sRNA;
  real logit_termination_rate_RNA;
  real<lower=0, upper=1> termination_rate_aRNA;
  real<lower=0> base_RNA_synth_sRNA;
  real<lower=0> sRNA_decay;
  real<lower=0> asRNA_decay;
  real<lower=0> complex_decay_sRNA;
  vector<lower=0>[n_RNA] f_RNA_vals;
  vector<lower=0>[n_asRNA] f_asRNA_vals;
  real<lower=0> f_sRNA_base;
  vector[n_sRNA > 0 ? n_sRNA - 1 : 0] d_sRNA_diff_raw;
  real<lower=0> sigma_intensity;
  real<lower=0> sigma_decay_t0;
  real<lower=0> sigma_decay_t24;
  real<lower=0> scale_nb_set1_2;
  real<lower=0> scale_nb_set3;
  real<lower=0> sigma_nb;
  real<lower=0> sigma_nb_t0;
}

transformed parameters {
  real y0[4];
  real RNA_free0_sKO = exp(log_RNA_free0_sKO);
  real RNA_free0_asKO = exp(log_RNA_free0_asKO);
  real local_RNA_free0_sKO = (fix_y0 == 0) ? RNA_free0_sKO : -1.0;
  real local_RNA_free0_asKO = (fix_y0 == 0) ? RNA_free0_asKO : -1.0;

  if (fix_y0 == 1) {
    y0 = y0_data;
    local_RNA_free0_sKO = y0[1];
    local_RNA_free0_asKO = y0[1];
  } else {
    real RNA_free0_calc;
    real sRNA_free0_calc;
    real RNA_complex_sRNA0_calc;
    real asRNA_free0_calc;
    real total_RNA_unit_t0_calc;

    RNA_free0_calc = exp(param_log_RNA_free0);
    sRNA_free0_calc = exp(param_log_sRNA_free0);
    RNA_complex_sRNA0_calc = (RNA_free0_calc * sRNA_free0_calc) / fmax(1e-9, Kd_sRNA);

    total_RNA_unit_t0_calc = RNA_free0_calc + RNA_complex_sRNA0_calc;
    asRNA_free0_calc = exp(log_asRNA_free0_param);

    y0[1] = RNA_free0_calc;
    y0[2] = RNA_complex_sRNA0_calc;
    y0[3] = asRNA_free0_calc;
    y0[4] = sRNA_free0_calc;
  }

  real termination_rate_RNA = inv_logit(logit_termination_rate_RNA);
  real k_off_sRNA = exp(log_k_off_sRNA);
  real k_on_sRNA_derived = k_off_sRNA / Kd_sRNA;

  real f_sRNA_vals[n_sRNA];
  if (n_sRNA > 0) {
    if (sRNA_strategy == 2) {
      if (n_sRNA == 1) {
        f_sRNA_vals[1] = f_sRNA_base;
      } else {
        vector[n_sRNA - 1] d_sRNA_log_diff = d_sRNA_diff_raw;
        f_sRNA_vals[1] = f_sRNA_base;
        for (i in 2:n_sRNA) {
          f_sRNA_vals[i] = fmax(0.0, f_sRNA_vals[i - 1] * exp(d_sRNA_log_diff[i - 1]));
        }
      }
    } else {
      for (i in 1:n_sRNA) {
        f_sRNA_vals[i] = f_sRNA_base;
      }
    }
  }

  real theta_base[K];
  theta_base[1] = base_RNA_synth;
  theta_base[2] = RNA_decay;
  theta_base[3] = k_on_sRNA_derived;
  theta_base[4] = k_off_sRNA;
  theta_base[5] = termination_rate_RNA;
  theta_base[6] = termination_rate_aRNA;
  theta_base[7] = base_RNA_synth_sRNA;
  theta_base[8] = sRNA_decay;
  theta_base[9] = asRNA_decay;
  theta_base[10] = complex_decay_sRNA;

  if (n_RNA > 0) {
    if (10 + n_RNA <= K)
      theta_base[(10 + 1) : (10 + n_RNA)] = to_array_1d(f_RNA_vals);
  }

  if (n_asRNA > 0) {
    int start_idx = 10 + n_RNA + 1;
    int end_idx = 10 + n_RNA + n_asRNA;
    if (start_idx <= K && end_idx <= K && start_idx <= end_idx)
      theta_base[start_idx : end_idx] = to_array_1d(f_asRNA_vals);
  }

  if (n_sRNA > 0) {
    int start_idx = 10 + n_RNA + n_asRNA + 1;
    if (start_idx <= K && (10 + L) <= K && start_idx <= (10 + L))
      theta_base[start_idx : K] = f_sRNA_vals;
  }

  real theta_with_flag[K + 1];
  if (K > 0)
    theta_with_flag[1:K] = theta_base;

  real y_ode[n_obs, 4];
  real y_ode_wt_nb1[n_obs_nb_wt1, 4];
  real y_ode_srna_ko[n_obs_nb_srna_ko, 4];
  real y_ode_wt_nb2[n_obs_nb_wt2, 4];
  real y_ode_asrna_ko[n_obs_nb_asrna_ko, 4];
  real RNA_pred[n_obs];
  real asRNA_pred[n_obs];
  real sRNA_pred[n_obs];
  real RNA_pred_nb_wt1[n_obs_nb_wt1];
  real RNA_pred_nb_srna_ko[n_obs_nb_srna_ko];
  real RNA_pred_nb_wt2[n_obs_nb_wt2];
  real RNA_pred_nb_asrna_ko[n_obs_nb_asrna_ko];
  real RNA_decay_pred[n_decay_obs];

  if (n_obs > 0) {
    theta_with_flag[K + 1] = 0.0;
    if (obs_times[1] == t0) {
      y_ode[1] = y0;
      if (n_obs > 1) {
        y_ode[2:n_obs] = integrate_ode_bdf(rna_ode, y0, t0, obs_times[2:n_obs], theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
      }
    } else {
      y_ode = integrate_ode_bdf(rna_ode, y0, t0, obs_times, theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
    }

    real baseline_RNA = y0[1] + y0[2] + 1e-12;
    real baseline_asRNA = y0[3] + 1e-12;
    real baseline_sRNA = y0[4] + y0[2] + 1e-12;

    for (i in 1:n_obs) {
      RNA_pred[i] = (y_ode[i, 1] + y_ode[i, 2]) / baseline_RNA;
      asRNA_pred[i] = y_ode[i, 3] / baseline_asRNA;
      sRNA_pred[i] = (y_ode[i, 4] + y_ode[i, 2]) / baseline_sRNA;
    }
  }

  if (n_obs_nb_wt1 > 0) {
    theta_with_flag[K + 1] = 0.0;
    if (nb_times_wt1[1] == t0) {
      y_ode_wt_nb1[1] = y0;
      if (n_obs_nb_wt1 > 1) {
        y_ode_wt_nb1[2:n_obs_nb_wt1] = integrate_ode_bdf(rna_ode, y0, t0, nb_times_wt1[2:n_obs_nb_wt1], theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
      }
    } else {
      y_ode_wt_nb1 = integrate_ode_bdf(rna_ode, y0, t0, nb_times_wt1, theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
    }
    for (i in 1:n_obs_nb_wt1) {
      RNA_pred_nb_wt1[i] = scale_nb_set1_2 * (y_ode_wt_nb1[i, 1] + y_ode_wt_nb1[i, 2]);
    }
  }

  if (n_obs_nb_srna_ko > 0) {
    real y0_srna_ko[4];
    y0_srna_ko[1] = local_RNA_free0_sKO;
    y0_srna_ko[2] = 0.0;
    y0_srna_ko[3] = y0[3];
    y0_srna_ko[4] = 0.0;
    theta_with_flag[K + 1] = 1.0;
    if (nb_times_srna_ko[1] == t0) {
      y_ode_srna_ko[1] = y0_srna_ko;
      if (n_obs_nb_srna_ko > 1) {
        y_ode_srna_ko[2:n_obs_nb_srna_ko] = integrate_ode_bdf(rna_ode, y0_srna_ko, t0, nb_times_srna_ko[2:n_obs_nb_srna_ko], theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
      }
    } else {
      y_ode_srna_ko = integrate_ode_bdf(rna_ode, y0_srna_ko, t0, nb_times_srna_ko, theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
    }
    for (i in 1:n_obs_nb_srna_ko) {
      RNA_pred_nb_srna_ko[i] = scale_nb_set1_2 * (y_ode_srna_ko[i, 1] + y_ode_srna_ko[i, 2]);
    }
  }

  if (n_obs_nb_wt2 > 0) {
    theta_with_flag[K + 1] = 0.0;
    if (nb_times_wt2[1] == t0) {
      y_ode_wt_nb2[1] = y0;
      if (n_obs_nb_wt2 > 1) {
        y_ode_wt_nb2[2:n_obs_nb_wt2] = integrate_ode_bdf(rna_ode, y0, t0, nb_times_wt2[2:n_obs_nb_wt2], theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
      }
    } else {
      y_ode_wt_nb2 = integrate_ode_bdf(rna_ode, y0, t0, nb_times_wt2, theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
    }
    for (i in 1:n_obs_nb_wt2) {
      RNA_pred_nb_wt2[i] = scale_nb_set3 * (y_ode_wt_nb2[i, 1] + y_ode_wt_nb2[i, 2]);
    }
  }

  if (n_obs_nb_asrna_ko > 0) {
    real y0_asrna_ko[4];
    y0_asrna_ko[1] = local_RNA_free0_asKO;
    y0_asrna_ko[2] = y0[2];
    y0_asrna_ko[3] = 0.0;
    y0_asrna_ko[4] = y0[4];
    theta_with_flag[K + 1] = 2.0;
    if (nb_times_asrna_ko[1] == t0) {
      y_ode_asrna_ko[1] = y0_asrna_ko;
      if (n_obs_nb_asrna_ko > 1) {
        y_ode_asrna_ko[2:n_obs_nb_asrna_ko] = integrate_ode_bdf(rna_ode, y0_asrna_ko, t0, nb_times_asrna_ko[2:n_obs_nb_asrna_ko], theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
      }
    } else {
      y_ode_asrna_ko = integrate_ode_bdf(rna_ode, y0_asrna_ko, t0, nb_times_asrna_ko, theta_with_flag, x_r, x_i, rel_tol, abs_tol, max_steps);
    }
    for (i in 1:n_obs_nb_asrna_ko) {
      RNA_pred_nb_asrna_ko[i] = scale_nb_set3 * (y_ode_asrna_ko[i, 1] + y_ode_asrna_ko[i, 2]);
    }
  }

  if (n_decay_obs > 0 && n_obs > 0) {
    for (i in 1:n_decay_obs) {
      if (decay_obs_indices[i] > 0) {
        int idx = decay_obs_indices[i];
        real sum_states = y_ode[idx, 1] + y_ode[idx, 2] + 1e-9;
        real w1_dynamic = y_ode[idx, 1] / sum_states;
        real w2_dynamic = y_ode[idx, 2] / sum_states;
        RNA_decay_pred[i] = fmax(0, RNA_decay * w1_dynamic + complex_decay_sRNA * w2_dynamic);
      } else {
        RNA_decay_pred[i] = -1.0;
      }
    }
  }
}

model {
  if (fix_y0 == 0) {
    param_log_RNA_free0 ~ normal(log(1.0), 1.5);
    param_log_sRNA_free0 ~ normal(log(6.0), 1.5);
    log_RNA_free0_sKO ~ normal(log(1.0), 1.5);
    log_RNA_free0_asKO ~ normal(log(3.5), 1.5);
    log_asRNA_free0_param ~ normal(log(14.28), 1.0);
  }

  log_k_off_sRNA ~ normal(k_off_log_mean_lit, k_off_log_sd_lit);
  log(base_RNA_synth) ~ normal(log(25.0), 1.0);
  log(RNA_decay) ~ normal(log(2), 0.5);
  Kd_sRNA ~ lognormal(kd_srna_log_mean_prior, kd_srna_log_sd_prior);
  logit_termination_rate_RNA ~ normal(-2.2, 0.5);
  termination_rate_aRNA ~ beta(1, 1);
  log(base_RNA_synth_sRNA) ~ normal(log(100.0), 1.5);
  log(sRNA_decay) ~ normal(log(15.0), 1);
  log(asRNA_decay) ~ normal(log(57.0), 1);
  log(complex_decay_sRNA) ~ normal(log(50.0), 1);

  if (n_RNA > 0)
    f_RNA_vals ~ normal(rna_vals_data, 1.5);
  if (n_asRNA > 0)
    f_asRNA_vals ~ normal(asrna_vals_data, 10.0);
  if (n_sRNA > 0) {
    if (sRNA_strategy == 2) {
      log(f_sRNA_base) ~ normal(log(1.0), 1.5);
      if (n_sRNA > 1)
        d_sRNA_diff_raw ~ normal(0, 3);
    } else {
      log(f_sRNA_base) ~ normal(log(1.0), 1.5);
    }
  }

  sigma_intensity ~ normal(0, 5.0);
  sigma_decay_t0 ~ normal(0, 2.0);
  sigma_decay_t24 ~ normal(0, 1.5);
  log(scale_nb_set1_2) ~ normal(0, 1.5);
  log(scale_nb_set3) ~ normal(0, 1.5);
  sigma_nb ~ normal(0, 1.0);
  sigma_nb_t0 ~ normal(0, 2.5);

  if (n_obs > 0) {
    target += normal_lpdf(to_vector(obs_RNA) | to_vector(RNA_pred), sigma_intensity);
    target += normal_lpdf(to_vector(obs_asRNA) | to_vector(asRNA_pred), sigma_intensity);
    target += normal_lpdf(to_vector(obs_sRNA) | to_vector(sRNA_pred), sigma_intensity);
  }

  if (n_decay_obs > 0) {
    for (j in 1:n_decay_obs) {
      if (decay_obs_indices[j] > 0) {
        if (decay_obs_times[j] == 0) {
          target += normal_lpdf(obs_RNA_decay[j] | RNA_decay_pred[j], sigma_decay_t0);
        } else if (decay_obs_times[j] == 24) {
          target += normal_lpdf(obs_RNA_decay[j] | RNA_decay_pred[j], sigma_decay_t24);
        }
      }
    }
  }

  if (n_obs_nb_wt1 > 0) {
    target += normal_lpdf(obs_RNA_nb_wt1[1] | RNA_pred_nb_wt1[1], sigma_nb_t0);
    if (n_obs_nb_wt1 > 1) {
      target += normal_lpdf(to_vector(obs_RNA_nb_wt1[2:n_obs_nb_wt1]) | to_vector(RNA_pred_nb_wt1[2:n_obs_nb_wt1]), sigma_nb);
    }
  }

  if (n_obs_nb_srna_ko > 0) {
    target += normal_lpdf(obs_RNA_nb_srna_ko[1] | RNA_pred_nb_srna_ko[1], sigma_nb_t0);
    if (n_obs_nb_srna_ko > 1) {
      target += normal_lpdf(to_vector(obs_RNA_nb_srna_ko[2:n_obs_nb_srna_ko]) | to_vector(RNA_pred_nb_srna_ko[2:n_obs_nb_srna_ko]), sigma_nb);
    }
  }

  if (n_obs_nb_wt2 > 0) {
    target += normal_lpdf(obs_RNA_nb_wt2[1] | RNA_pred_nb_wt2[1], sigma_nb_t0);
    if (n_obs_nb_wt2 > 1) {
      target += normal_lpdf(to_vector(obs_RNA_nb_wt2[2:n_obs_nb_wt2]) | to_vector(RNA_pred_nb_wt2[2:n_obs_nb_wt2]), sigma_nb);
    }
  }

  if (n_obs_nb_asrna_ko > 0) {
    target += normal_lpdf(obs_RNA_nb_asrna_ko[1] | RNA_pred_nb_asrna_ko[1], sigma_nb_t0);
    if (n_obs_nb_asrna_ko > 1) {
      target += normal_lpdf(to_vector(obs_RNA_nb_asrna_ko[2:n_obs_nb_asrna_ko]) | to_vector(RNA_pred_nb_asrna_ko[2:n_obs_nb_asrna_ko]), sigma_nb);
    }
  }
}

generated quantities {
  vector[N_total_obs] log_lik;
  int current_log_lik_idx = 1;
  real RNA_pred_out[n_obs];
  real asRNA_pred_out[n_obs];
  real sRNA_pred_out[n_obs];
  real RNA_decay_pred_out[n_decay_obs];
  real RNA_pred_nb_wt1_out[n_obs_nb_wt1];
  real RNA_pred_nb_srna_ko_out[n_obs_nb_srna_ko];
  real RNA_pred_nb_wt2_out[n_obs_nb_wt2];
  real RNA_pred_nb_asrna_ko_out[n_obs_nb_asrna_ko];
  real sRNA_pred_nb_wt1_out[n_obs_nb_wt1];
  real asRNA_pred_nb_wt1_out[n_obs_nb_wt1];
  real sRNA_pred_nb_srna_ko_out[n_obs_nb_srna_ko];
  real asRNA_pred_nb_srna_ko_out[n_obs_nb_srna_ko];
  real sRNA_pred_nb_wt2_out[n_obs_nb_wt2];
  real asRNA_pred_nb_wt2_out[n_obs_nb_wt2];
  real sRNA_pred_nb_asrna_ko_out[n_obs_nb_asrna_ko];
  real asRNA_pred_nb_asrna_ko_out[n_obs_nb_asrna_ko];
  vector[n_RNA] f_RNA_vals_out;
  vector[n_asRNA] f_asRNA_vals_out;
  real f_sRNA_vals_out[n_sRNA];
  real y0_out[4];
  real Kd_sRNA_ArbUnits_calc = Kd_sRNA;
  real k_on_sRNA_gq = k_on_sRNA_derived;

  real ratio_total_sRNA_to_total_RNA_at_t0;
  if (fix_y0 == 0) {
    real total_RNA_unit_t0_gq = y0[1] + y0[2];
    if (total_RNA_unit_t0_gq > 1e-9) {
      ratio_total_sRNA_to_total_RNA_at_t0 = (y0[4] + y0[2]) / total_RNA_unit_t0_gq;
    } else {
      ratio_total_sRNA_to_total_RNA_at_t0 = -1.0;
    }
  } else {
    ratio_total_sRNA_to_total_RNA_at_t0 = -1.0;
  }

  if (n_obs > 0) {
    for (i in 1:n_obs) {
      log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA[i] | RNA_pred[i], sigma_intensity);
      current_log_lik_idx += 1;
    }
    for (i in 1:n_obs) {
      log_lik[current_log_lik_idx] = normal_lpdf(obs_asRNA[i] | asRNA_pred[i], sigma_intensity);
      current_log_lik_idx += 1;
    }
    for (i in 1:n_obs) {
      log_lik[current_log_lik_idx] = normal_lpdf(obs_sRNA[i] | sRNA_pred[i], sigma_intensity);
      current_log_lik_idx += 1;
    }
  }

  if (n_obs_nb_wt1 > 0) {
    log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA_nb_wt1[1] | RNA_pred_nb_wt1[1], sigma_nb_t0);
    current_log_lik_idx += 1;
    if (n_obs_nb_wt1 > 1) {
      for (i in 2:n_obs_nb_wt1) {
        log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA_nb_wt1[i] | RNA_pred_nb_wt1[i], sigma_nb);
        current_log_lik_idx += 1;
      }
    }
  }

  if (n_obs_nb_srna_ko > 0) {
    log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA_nb_srna_ko[1] | RNA_pred_nb_srna_ko[1], sigma_nb_t0);
    current_log_lik_idx += 1;
    if (n_obs_nb_srna_ko > 1) {
      for (i in 2:n_obs_nb_srna_ko) {
        log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA_nb_srna_ko[i] | RNA_pred_nb_srna_ko[i], sigma_nb);
        current_log_lik_idx += 1;
      }
    }
  }

  if (n_obs_nb_wt2 > 0) {
    log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA_nb_wt2[1] | RNA_pred_nb_wt2[1], sigma_nb_t0);
    current_log_lik_idx += 1;
    if (n_obs_nb_wt2 > 1) {
      for (i in 2:n_obs_nb_wt2) {
        log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA_nb_wt2[i] | RNA_pred_nb_wt2[i], sigma_nb);
        current_log_lik_idx += 1;
      }
    }
  }

  if (n_obs_nb_asrna_ko > 0) {
    log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA_nb_asrna_ko[1] | RNA_pred_nb_asrna_ko[1], sigma_nb_t0);
    current_log_lik_idx += 1;
    if (n_obs_nb_asrna_ko > 1) {
      for (i in 2:n_obs_nb_asrna_ko) {
        log_lik[current_log_lik_idx] = normal_lpdf(obs_RNA_nb_asrna_ko[i] | RNA_pred_nb_asrna_ko[i], sigma_nb);
        current_log_lik_idx += 1;
      }
    }
  }

  if (n_obs > 0) {
    for (i in 1:n_obs) {
      RNA_pred_out[i] = RNA_pred[i];
      asRNA_pred_out[i] = asRNA_pred[i];
      sRNA_pred_out[i] = sRNA_pred[i];
    }
  }
  if (n_decay_obs > 0) {
    for (i in 1:n_decay_obs) {
      RNA_decay_pred_out[i] = (decay_obs_indices[i] > 0) ? RNA_decay_pred[i] : -1.0;
    }
  }
  if (n_obs_nb_wt1 > 0) {
    RNA_pred_nb_wt1_out = RNA_pred_nb_wt1;
    for (i in 1:n_obs_nb_wt1) {
      sRNA_pred_nb_wt1_out[i] = scale_nb_set1_2 * (y_ode_wt_nb1[i, 4] + y_ode_wt_nb1[i, 2]);
      asRNA_pred_nb_wt1_out[i] = scale_nb_set1_2 * y_ode_wt_nb1[i, 3];
    }
  }
  if (n_obs_nb_srna_ko > 0) {
    RNA_pred_nb_srna_ko_out = RNA_pred_nb_srna_ko;
    for (i in 1:n_obs_nb_srna_ko) {
      sRNA_pred_nb_srna_ko_out[i] = scale_nb_set1_2 * (y_ode_srna_ko[i, 4] + y_ode_srna_ko[i, 2]);
      asRNA_pred_nb_srna_ko_out[i] = scale_nb_set1_2 * y_ode_srna_ko[i, 3];
    }
  }
  if (n_obs_nb_wt2 > 0) {
    RNA_pred_nb_wt2_out = RNA_pred_nb_wt2;
    for (i in 1:n_obs_nb_wt2) {
      sRNA_pred_nb_wt2_out[i] = scale_nb_set3 * (y_ode_wt_nb2[i, 4] + y_ode_wt_nb2[i, 2]);
      asRNA_pred_nb_wt2_out[i] = scale_nb_set3 * y_ode_wt_nb2[i, 3];
    }
  }
  if (n_obs_nb_asrna_ko > 0) {
    RNA_pred_nb_asrna_ko_out = RNA_pred_nb_asrna_ko;
    for (i in 1:n_obs_nb_asrna_ko) {
      sRNA_pred_nb_asrna_ko_out[i] = scale_nb_set3 * (y_ode_asrna_ko[i, 4] + y_ode_asrna_ko[i, 2]);
      asRNA_pred_nb_asrna_ko_out[i] = scale_nb_set3 * y_ode_asrna_ko[i, 3];
    }
  }

  if (n_RNA > 0)
    f_RNA_vals_out = f_RNA_vals;
  if (n_asRNA > 0)
    f_asRNA_vals_out = f_asRNA_vals;
  if (n_sRNA > 0) {
    for (i in 1:n_sRNA)
      f_sRNA_vals_out[i] = f_sRNA_vals[i];
  }
  y0_out = y0;
  real RNA_free0_sKO_out = (fix_y0 == 0) ? RNA_free0_sKO : -1.0;
  real RNA_free0_asKO_out = (fix_y0 == 0) ? RNA_free0_asKO : -1.0;
  real k_off_sRNA_gq = k_off_sRNA;
  real Kd_sRNA_gq = Kd_sRNA;
}

"



stan_model_fit <- stan_model(model_code = stan_model_code)
