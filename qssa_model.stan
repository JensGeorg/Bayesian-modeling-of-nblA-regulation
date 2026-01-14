functions {
  // Linear interpolation for synthesis rates
  real interp_forcing(real t, real v0, real v9, real v24) {
    if (t <= 0) return v0;
    else if (t <= 9) return v0 + (v9 - v0) * t / 9.0;
    else if (t <= 24) return v9 + (v24 - v9) * (t - 9.0) / 15.0;
    else return v24;
  }
  
  // QSSA solver (sRNA complex only)
  real solve_complex(real R_tot, real S_tot, real Kd) {
    real Rt = fmax(0.0, R_tot);
    real St = fmax(0.0, S_tot);
    real sum_all = Rt + St + Kd;
    real discrim = square(sum_all) - 4.0 * Rt * St;
	real C_raw = (sum_all - sqrt(fmax(discrim, 1e-12))) / 2.0;
	return log1p_exp(10.0 * C_raw) / 10.0;
  }
  
  // Effective decay calc
  real calc_eff_decay(real R_tot, real S_tot, 
                      real Kd_sRNA,
                      real k_free, real k_complex_sRNA) {
    real C_sRNA = solve_complex(R_tot, S_tot, Kd_sRNA);
    real R_free = fmax(0.0, R_tot - C_sRNA);
    real total_decay = k_free * R_free + k_complex_sRNA * C_sRNA;
    return total_decay / (R_tot + 1e-9);
  }
  
  // Calculate steady-state RNA given synthesis, decay parameters, and sRNA
  // At steady state: syn_total = decay_R_free * R_free + decay_complex * C
  // syn_total = syn_eff_main + leak
  real calc_R_steady_state(real syn_R_total, real S_tot, real Kd_sRNA,
                           real decay_R_free, real decay_complex_sRNA) {
    // Iterative solver for R_tot at steady state
    real R_tot = syn_R_total / decay_R_free;  // Initial guess (no complex)
    
    for (iter in 1:20) {
      real C = solve_complex(R_tot, S_tot, Kd_sRNA);
      real R_free = fmax(0.0, R_tot - C);
      real R_tot_new = (syn_R_total - (decay_complex_sRNA - decay_R_free) * C) / decay_R_free;
      R_tot_new = fmax(1e-6, R_tot_new);
      
      // Check convergence
      if (abs(R_tot_new - R_tot) < 1e-6 * R_tot) break;
      R_tot = R_tot_new;
    }
    return fmax(1e-6, R_tot);
  }
  // EXACT STEADY STATE SOLVER (Analytic Solution)
  // Solves for Complex (C) directly, then infers Totals.
  // Returns vector: [R_tot, S_tot, C]
  vector solve_steady_state_analytic(real syn_R, real syn_S, real Kd, 
                                     real k_R, real k_S, real k_C) {
    // Coefficients for Quadratic: a*C^2 + b*C + c = 0
    // Derived from: (Syn_R - k_C*C)/k_R * (Syn_S - k_C*C)/k_S = Kd * C
    
    real a = square(k_C);
    real b = -1.0 * (k_C * (syn_R + syn_S) + k_R * k_S * Kd);
    real c = syn_R * syn_S;
    

    real discriminant = square(b) - 4.0 * a * c;
    real C = (-b - sqrt(fmax(0.0, discriminant))) / (2.0 * a);
    
    C = fmax(0.0, C);
    
    // Calculate Free species
    // Constraint: Free species must be non-negative. 
    // If Syn is very low, C might absorb all flux.
    real R_free = fmax(0.0, (syn_R - k_C * C) / k_R);
    real S_free = fmax(0.0, (syn_S - k_C * C) / k_S);
    
    return [R_free + C, S_free + C, C]';
  }

  // ODE system
  vector rna_ode(real t, vector y, 
                 real Kd_sRNA,
                 real decay_R_free, real decay_complex_sRNA,
                 real decay_S_free, real decay_asR,
                 real base_syn, real base_syn_S,
                 real term_rate_asR_on_RNA, real term_rate_RNA_on_asR,
                 real synthesis_leak, 
                 real f_RNA_9, real f_RNA_24,
                 real f_asRNA_9, real f_asRNA_24,
                 real f_sRNA_9, real f_sRNA_24,
                 int is_sRNA_KO, int is_asRNA_KO) {
    
    real R_tot = y[1];
    real S_tot = y[2];
    real asR_tot = y[3];

    // synthesis rate functions
    real f_R = interp_forcing(t, 1.0, f_RNA_9, f_RNA_24);
    real f_asR = interp_forcing(t, 102.222, f_asRNA_9, f_asRNA_24);
    real f_S = interp_forcing(t, 1.0, f_sRNA_9, f_sRNA_24);

    real syn_R_raw = base_syn * f_R;
    real syn_asR_raw = base_syn * f_asR;
    real syn_S_raw = base_syn_S * f_S;
    
    if (is_sRNA_KO == 1) syn_S_raw = 0;
    if (is_asRNA_KO == 1) syn_asR_raw = 0;

    // TI (Transcriptional Interference)
    // 1. Calculate main synthesis (can go to 0.0)
	real diff_ode = syn_R_raw - term_rate_asR_on_RNA * syn_asR_raw;
	real syn_R_eff_main = log1p_exp(1.0 * diff_ode) / 1.0;
    
    // 2. Add Leak (always present) -> Total Effective Synthesis acounts for existing premature transcripts
    real syn_R_eff_total = syn_R_eff_main + synthesis_leak;
	 
    
    // asRNA interference
	real diff_ode_as = syn_asR_raw - term_rate_RNA_on_asR * syn_R_raw;
	real syn_asR_eff = log1p_exp(5.0 * diff_ode_as) / 5.0;
    
    // Complex
    real C_sRNA = solve_complex(R_tot, S_tot, Kd_sRNA);
    real R_free = fmax(0.0, R_tot - C_sRNA);
    real S_free = fmax(0.0, S_tot - C_sRNA);

    // Derivatives
    real dR = syn_R_eff_total - decay_R_free * R_free - decay_complex_sRNA * C_sRNA;
    real dS = syn_S_raw - decay_S_free * S_free - decay_complex_sRNA * C_sRNA;
    real dasR = syn_asR_eff - decay_asR * asR_tot;

    return [dR, dS, dasR]';
  }
}

data {
  int<lower=1> n_solver;
  array[n_solver] real solver_times;
  real t0;
  
  // Data
  real obs_asRNA_RNA_ratio_t9;
  real obs_asRNA_RNA_ratio_t0;
  real obs_sRNA_RNA_ratio_t0;
  
  int<lower=0> n_microarray;
  array[n_microarray] int microarray_time_idx;
  array[n_microarray] real obs_RNA_microarray;
  array[n_microarray] real obs_sRNA_microarray;
  array[n_microarray] real obs_asRNA_microarray;
  
  real obs_asrna_ko_fc_t9_t0;
  real obs_asrna_ko_wt_ratio_t9;
  
  int<lower=1> n_nb1_srna_ko;
  array[n_nb1_srna_ko] int nb1_srna_ko_time_idx;
  array[n_nb1_srna_ko] real obs_RNA_srna_ko_fc_from_t2;
  int srna_ko_t2_idx;
  
  int<lower=1> n_nb1_wt;
  array[n_nb1_wt] int nb1_wt_time_idx;
  array[n_nb1_wt] real obs_RNA_nb1_wt_fc_from_t3;
  int wt_t3_idx;
  
  int<lower=1> n_decay_obs;
  array[n_decay_obs] int decay_time_idx;
  array[n_decay_obs] real obs_RNA_decay;
  
  real f_RNA_9_data;
  real f_asRNA_9_data;
  real f_asRNA_24_data;
  
  real obs_ratio_wt_srna_ko_t6; 
  real obs_ratio_asrna_ko_wt_t8;
  real obs_ratio_asrna_ko_wt_t24;
}

parameters {
  // Kinetic Parameters
  real<lower=10, upper=5000> Kd_sRNA;
  real<lower=0.1, upper=5> decay_R_free;
  real<lower=5, upper=50> decay_complex_sRNA;
  real<lower=5, upper=100> decay_S_free;
  real<lower=10, upper=100> decay_asR; 
  
  real<lower=10, upper=50000> base_syn;
  real<lower=500, upper=50000> base_syn_S;

  // Interference parameters
  // asRNA synthesis is 102x RNA synthesis at t=0
  real<lower=0.001, upper=0.1> term_rate_asR_on_RNA;  
  real<lower=0, upper=1> term_rate_RNA_on_asR;
  
  real<lower=1, upper=200> synthesis_leak;
  
  // Synthesis rates
  real<lower=1, upper=50> f_RNA_9;
  real<lower=1, upper=20> f_RNA_24;
  real<lower=50, upper=300> f_asRNA_9;
  real<lower=50, upper=300> f_asRNA_24;
  real<lower=0, upper=10> f_sRNA_9;
  real<lower=0, upper=10> f_sRNA_24;
  
  // Errors
  real<lower=0> sigma_rnaseq;
  real<lower=0> sigma_microarray;
  real<lower=0> sigma_qpcr;
  real<lower=0> sigma_nb1;
  real<lower=0> sigma_decay;
  real<lower=0> sigma_reporter;
}

transformed parameters {
  // 1. Raw synthesis rates at t=0
  real syn_R_raw_0 = base_syn;             // f_R(0)=1
  real syn_asR_raw_0 = base_syn * 102.222; // f_asR(0)=102.22
  real syn_S_raw_0 = base_syn_S;           // f_S(0)=1

  // 2. WT Effective Synthesis (TI + Leak)
  real diff_term = syn_R_raw_0 - term_rate_asR_on_RNA * syn_asR_raw_0;
  real syn_R_main_wt_0 = log1p_exp(5.0 * diff_term) / 5.0;
  real syn_R_total_wt_0 = syn_R_main_wt_0 + synthesis_leak;
  
  real syn_asR_eff_0 = fmax(1e-6, syn_asR_raw_0 - term_rate_RNA_on_asR * syn_R_raw_0);

  // 3. asRNA KO Effective Synthesis (No TI, just Leak + Raw)
  real syn_R_total_asrna_ko_0 = syn_R_raw_0 + synthesis_leak;

  // 4. asRNA steady state
  real asR0_tot = syn_asR_eff_0 / decay_asR;

  // WT Steady State
  vector[3] ss_wt = solve_steady_state_analytic(
      syn_R_total_wt_0, 
      syn_S_raw_0, 
      Kd_sRNA, 
      decay_R_free, decay_S_free, decay_complex_sRNA
  );
  real R0_wt = fmax(1e-6, ss_wt[1]);
  real S0_wt = fmax(1e-6, ss_wt[2]);
  
  real S0_tot = S0_wt;
  
  // sRNA KO Steady State (Simple: sRNA=0, so R = Syn/k_R)
  real R0_srna_ko = syn_R_total_wt_0 / decay_R_free;
  R0_srna_ko = fmax(1e-6, R0_srna_ko);

  // asRNA KO Steady State (Uses different RNA synthesis rate)
  vector[3] ss_asrna_ko = solve_steady_state_analytic(
      syn_R_total_asrna_ko_0, 
      syn_S_raw_0, 
      Kd_sRNA, 
      decay_R_free, decay_S_free, decay_complex_sRNA
  );
  real R0_asrna_ko = fmax(1e-6, ss_asrna_ko[1]);
  real S0_asrna_ko = fmax(1e-6, ss_asrna_ko[2]);

  real C_wt_t0_tp = solve_complex(R0_wt, S0_wt, Kd_sRNA);
  real R_free_wt_t0_tp = fmax(0.0, R0_wt - C_wt_t0_tp);
  
  real expected_leak = decay_R_free * R_free_wt_t0_tp + decay_complex_sRNA * C_wt_t0_tp;

  vector[3] y0_wt = [R0_wt, S0_wt, asR0_tot]';
  vector[3] y0_srna_ko = [R0_srna_ko, 0.0, asR0_tot]';
  vector[3] y0_asrna_ko = [R0_asrna_ko, S0_asrna_ko, 0.0]';

  array[n_solver] vector[3] y_hat = ode_rk45_tol(
    rna_ode, y0_wt, t0, solver_times,
    1e-5, 1e-5, 5000,
    Kd_sRNA, decay_R_free, decay_complex_sRNA, decay_S_free, decay_asR,
    base_syn, base_syn_S, term_rate_asR_on_RNA, term_rate_RNA_on_asR,
    synthesis_leak, f_RNA_9, f_RNA_24, f_asRNA_9, f_asRNA_24, f_sRNA_9, f_sRNA_24,
    0, 0);

  array[n_solver] vector[3] y_hat_asrna_ko = ode_rk45_tol(
    rna_ode, y0_asrna_ko, t0, solver_times,
    1e-5, 1e-5, 5000,
    Kd_sRNA, decay_R_free, decay_complex_sRNA, decay_S_free, decay_asR,
    base_syn, base_syn_S, term_rate_asR_on_RNA, term_rate_RNA_on_asR,
    synthesis_leak, f_RNA_9, f_RNA_24, f_asRNA_9, f_asRNA_24, f_sRNA_9, f_sRNA_24,
    0, 1);

  array[n_solver] vector[3] y_hat_srna_ko = ode_rk45_tol(
    rna_ode, y0_srna_ko, t0, solver_times,
    1e-5, 1e-5, 5000,
    Kd_sRNA, decay_R_free, decay_complex_sRNA, decay_S_free, decay_asR,
    base_syn, base_syn_S, term_rate_asR_on_RNA, term_rate_RNA_on_asR,
    synthesis_leak, f_RNA_9, f_RNA_24, f_asRNA_9, f_asRNA_24, f_sRNA_9, f_sRNA_24,
    1, 0);
}

model {
  //PRIORS
  Kd_sRNA ~ lognormal(log(600), 0.5);
  decay_R_free ~ lognormal(log(1.5), 0.2); 
  decay_complex_sRNA ~ lognormal(log(19), 0.2);
  decay_S_free ~ lognormal(log(10), 0.3);
  decay_asR ~ lognormal(log(60), 0.8); 
  
  // SYNTHESIS: HIGH FLUX 
  base_syn ~ lognormal(log(25000), 0.5); 
  base_syn_S ~ lognormal(log(15000), 0.5);

  //INTERFERENCE 
  term_rate_RNA_on_asR ~ beta(2, 2);
  
  // LEAK - constrained by steady-state balance at t=0
  // At steady state: synthesis_leak = decay_R_free * R_free + decay_complex * C
  synthesis_leak ~ normal(expected_leak, fmax(expected_leak * 0.1, 10)); 
  
  //SIGMAS & synthesis rates
  f_sRNA_9 ~ lognormal(log(0.3), 0.2);   
  f_sRNA_24 ~ lognormal(log(0.8), 0.2); 
  
  sigma_rnaseq ~ exponential(5);
  sigma_microarray ~ exponential(20);
  sigma_qpcr ~ exponential(10);
  sigma_nb1 ~ exponential(5);
  sigma_decay ~ exponential(5);
  sigma_reporter ~ exponential(2);
  
  f_RNA_9 ~ lognormal(log(f_RNA_9_data), sigma_reporter);
  f_RNA_24 ~ lognormal(log(5.5), 1); 
  f_asRNA_9 ~ lognormal(log(f_asRNA_9_data), sigma_reporter);
  f_asRNA_24 ~ lognormal(log(f_asRNA_24_data), sigma_reporter);
  
  // LIKELIHOODS
  
  // 1. t=0 Ratios
  real pred_asRNA_ratio_t0 = asR0_tot / R0_wt;
  obs_asRNA_RNA_ratio_t0 ~ lognormal(log(pred_asRNA_ratio_t0), sigma_rnaseq);

  real pred_sRNA_ratio_t0 = S0_tot / R0_wt;
  obs_sRNA_RNA_ratio_t0 ~ lognormal(log(pred_sRNA_ratio_t0), sigma_rnaseq);

  // asRNA_KO / WT ratio at t=0 (qPCR data)
  // Data shows: WT(t=0) ≈ 0 (1e-11), asRNA_KO(t=0) = 1.82
  // At t=8: ratio = 4.47/0.805 = 5.55x
  // At t=24: ratio = 3.35/1.31 = 2.56x  
  // The ratio DECREASES over time, so at t=0 it should be HIGHER than 5.55x
  // Since WT is essentially at detection limit at t=0, we expect a large ratio
  {
    real ratio_ko_wt_t0 = R0_asrna_ko / fmax(R0_wt, 1e-6);
    ratio_ko_wt_t0 ~ lognormal(log(20), 0.7);
  }
  
  // WT RNA fold-change from t=0 to t=9 (RNA-seq) 
  // RNA-seq data: RNA_t0 = 19.9, RNA_t9 = 923
  // Fold change = 923 / 19.9 = 46.4x
  // This constrains that WT at t=0 should be much lower than at t=9
  // Using wider uncertainty to allow flexibility with NB1 data
  {
    int t9_idx = 0;
    for (i in 1:n_solver) {
      if (abs(solver_times[i] - 9.0) < 0.1) t9_idx = i;
    }
    if (t9_idx > 0) {
      real R_wt_t9 = fmax(y_hat[t9_idx][1], 1e-6);
      real fc_t9_t0 = R_wt_t9 / fmax(R0_wt, 1e-6);
      fc_t9_t0 ~ lognormal(log(46.4), 0.5);
    }
  }

  // Ratio at t=9
  {
    int t9_idx = 0;
    for (i in 1:n_solver) {
      if (abs(solver_times[i] - 9.0) < 0.1) t9_idx = i;
    }
    if (t9_idx > 0) {
      real pred_ratio_t9 = fmax(y_hat[t9_idx][3], 1e-6) / fmax(y_hat[t9_idx][1], 1e-6);
      obs_asRNA_RNA_ratio_t9 ~ lognormal(log(pred_ratio_t9), sigma_rnaseq);
    }
  }
  
  // Microarray Data
  if (n_microarray > 0) {
    vector[3] y_ref = y_hat[microarray_time_idx[1]];
    real RNA_ref = fmax(y_ref[1], 1e-6);
    real sRNA_ref = fmax(y_ref[2], 1e-6);
    real asRNA_ref = fmax(y_ref[3], 1e-6);
    for (i in 1:n_microarray) {
      vector[3] y_val = y_hat[microarray_time_idx[i]];
      obs_RNA_microarray[i] / obs_RNA_microarray[1] ~ 
        lognormal(log(fmax(y_val[1], 1e-6) / RNA_ref), sigma_microarray);
      obs_sRNA_microarray[i] / obs_sRNA_microarray[1] ~ 
        lognormal(log(fmax(y_val[2], 1e-6) / sRNA_ref), sigma_microarray);
      obs_asRNA_microarray[i] / obs_asRNA_microarray[1] ~ 
        lognormal(log(fmax(y_val[3], 1e-6) / asRNA_ref), sigma_microarray);
    }
  }
  
  // asRNA Knockout qPCR
  {
    int t9_idx = 0;
    for (i in 1:n_solver) {
      if (abs(solver_times[i] - 9.0) < 0.1) t9_idx = i;
    }
    if (t9_idx > 0) {
      real R_wt_t9 = fmax(y_hat[t9_idx][1], 1e-6);
      real R_ko_t9 = fmax(y_hat_asrna_ko[t9_idx][1], 1e-6);
      
      // Fold-Change t9/t0 for asRNA KO
      obs_asrna_ko_fc_t9_t0 ~ lognormal(log(R_ko_t9 / R0_asrna_ko), sigma_qpcr);
      // KO vs WT Ratio at t9
      obs_asrna_ko_wt_ratio_t9 ~ lognormal(log(R_ko_t9 / R_wt_t9), sigma_qpcr);
    }
  }
  
  // sRNA Knockout - relative to t=2
  {
    real R_srna_ko_t2 = fmax(y_hat_srna_ko[srna_ko_t2_idx][1], 1e-6);
    for (i in 1:n_nb1_srna_ko) {
      real R_srna_ko = fmax(y_hat_srna_ko[nb1_srna_ko_time_idx[i]][1], 1e-6);
      obs_RNA_srna_ko_fc_from_t2[i] ~ lognormal(log(R_srna_ko / R_srna_ko_t2), sigma_nb1);
    }
  }
  
  // WT NB1 Data - relative to t=3
  {
    real R_wt_t3 = fmax(y_hat[wt_t3_idx][1], 1e-6);
    for (i in 1:n_nb1_wt) {
      real R_wt = fmax(y_hat[nb1_wt_time_idx[i]][1], 1e-6);
      obs_RNA_nb1_wt_fc_from_t3[i] ~ lognormal(log(R_wt / R_wt_t3), sigma_nb1);
    }
  }
  
  // Decay Data
  for (i in 1:n_decay_obs) {
    vector[3] y_val = decay_time_idx[i] == 0 ? y0_wt : y_hat[decay_time_idx[i]];
    real eff_decay = calc_eff_decay(fmax(y_val[1], 1e-6), fmax(y_val[2], 1e-6), 
                                     Kd_sRNA,
                                     decay_R_free, decay_complex_sRNA);
	obs_RNA_decay[i] ~ lognormal(log(fmax(eff_decay, 1e-9)), sigma_decay);
  }
  
  // WT vs sRNA KO at t=6 
  {
    int t6_idx = 0;
    for (i in 1:n_solver) {
      if (abs(solver_times[i] - 6.0) < 0.1) t6_idx = i;
    }
    if (t6_idx > 0) {
       real wt_val = fmax(y_hat[t6_idx][1], 1e-6);
       real srna_ko_val = fmax(y_hat_srna_ko[t6_idx][1], 1e-6);
       // We fit the ratio WT / sRNA_KO
       obs_ratio_wt_srna_ko_t6 ~ lognormal(log(wt_val / srna_ko_val), sigma_nb1);
    }
  }

  // asRNA KO vs WT at t=8 and t=24 
  {
    int t8_idx = 0;
    int t24_idx = 0;
    for (i in 1:n_solver) {
      if (abs(solver_times[i] - 8.0) < 0.1) t8_idx = i;
      if (abs(solver_times[i] - 24.0) < 0.1) t24_idx = i;
    }
    
    if (t8_idx > 0) {
       real wt_8 = fmax(y_hat[t8_idx][1], 1e-6);
       real ko_8 = fmax(y_hat_asrna_ko[t8_idx][1], 1e-6);
       obs_ratio_asrna_ko_wt_t8 ~ lognormal(log(ko_8 / wt_8), sigma_qpcr);
    }
    
    if (t24_idx > 0) {
       real wt_24 = fmax(y_hat[t24_idx][1], 1e-6);
       real ko_24 = fmax(y_hat_asrna_ko[t24_idx][1], 1e-6);
       obs_ratio_asrna_ko_wt_t24 ~ lognormal(log(ko_24 / wt_24), sigma_qpcr);
    }
  }
}

generated quantities {
  // Trajectories
  array[n_solver] real RNA_wt_trajectory;
  array[n_solver] real RNA_asrna_ko_trajectory;
  array[n_solver] real RNA_srna_ko_trajectory;
  array[n_solver] real decay_eff_wt;
  
  // Initial values
  real RNA_wt_t0 = R0_wt;
  real RNA_asrna_ko_t0 = R0_asrna_ko;
  real RNA_srna_ko_t0 = R0_srna_ko;
  
  for (i in 1:n_solver) {
    RNA_wt_trajectory[i] = fmax(y_hat[i][1], 1e-6);
    RNA_asrna_ko_trajectory[i] = fmax(y_hat_asrna_ko[i][1], 1e-6);
    RNA_srna_ko_trajectory[i] = fmax(y_hat_srna_ko[i][1], 1e-6);
    
    decay_eff_wt[i] = calc_eff_decay(
        fmax(y_hat[i][1], 1e-6),
        fmax(y_hat[i][2], 1e-6),
        Kd_sRNA,
        decay_R_free, 
        decay_complex_sRNA
    );
  }

  real decay_eff_wt_t0 = calc_eff_decay(R0_wt, S0_wt, Kd_sRNA,
                                         decay_R_free, decay_complex_sRNA);
                                         
  
  // Check steady-state 
  // For WT at t=0:
  real C_wt_t0 = solve_complex(R0_wt, S0_wt, Kd_sRNA);
  real R_free_wt_t0 = fmax(0.0, R0_wt - C_wt_t0);
  real dR_dt_wt_t0 = syn_R_total_wt_0 - decay_R_free * R_free_wt_t0 - decay_complex_sRNA * C_wt_t0;
  
  // For asRNA KO at t=0:
  real C_asrna_ko_t0 = solve_complex(R0_asrna_ko, S0_asrna_ko, Kd_sRNA);
  real R_free_asrna_ko_t0 = fmax(0.0, R0_asrna_ko - C_asrna_ko_t0);
  real dR_dt_asrna_ko_t0 = syn_R_total_asrna_ko_0 - decay_R_free * R_free_asrna_ko_t0 - decay_complex_sRNA * C_asrna_ko_t0;
  
  
  // 1. Synthesis rates at key timepoints
  real syn_R_raw_t0 = base_syn * 1.0;
  real syn_R_raw_t9 = base_syn * f_RNA_9;
  real syn_asR_raw_t0 = base_syn * 102.222;
  real syn_asR_raw_t9 = base_syn * f_asRNA_9;
  
  // 2. Interference at t=0 and t=9
  real interference_t0 = term_rate_asR_on_RNA * syn_asR_raw_t0;
  real interference_t9 = term_rate_asR_on_RNA * syn_asR_raw_t9;
  real blocking_fraction_t0 = interference_t0 / syn_R_raw_t0;
  real blocking_fraction_t9 = interference_t9 / syn_R_raw_t9;
  
  // 3. Effective synthesis at t=0 and t=9
  real syn_R_eff_main_t0 = fmax(0.0, syn_R_raw_t0 - interference_t0);
  real syn_R_eff_main_t9 = fmax(0.0, syn_R_raw_t9 - interference_t9);
  real syn_R_total_t0 = syn_R_eff_main_t0 + synthesis_leak;
  real syn_R_total_t9 = syn_R_eff_main_t9 + synthesis_leak;
  
  // 4. Synthesis fold-changes
  real FC_raw_synthesis = syn_R_raw_t9 / syn_R_raw_t0;
  real FC_eff_synthesis = syn_R_total_t9 / syn_R_total_t0;
  real leak_fraction_t0 = synthesis_leak / syn_R_total_t0;
  real leak_fraction_t9 = synthesis_leak / syn_R_total_t9;
  
  // 5. Decay rate analysis 
  // Get decay at t=9
  int t9_idx_decay = 0;
  for (i in 1:n_solver) {
    if (abs(solver_times[i] - 9.0) < 0.1) t9_idx_decay = i;
  }
  real decay_eff_wt_t9 = t9_idx_decay > 0 ? 
    calc_eff_decay(fmax(y_hat[t9_idx_decay][1], 1e-6),
                   fmax(y_hat[t9_idx_decay][2], 1e-6),
                   Kd_sRNA, decay_R_free, decay_complex_sRNA) : 0;
  
  real FC_decay_rate = decay_eff_wt_t9 / decay_eff_wt_t0;
  
  // 6. Complex formation at t=0 and t=9
  real fraction_bound_t0 = C_wt_t0 / fmax(R0_wt, 1e-6);
  real C_wt_t9 = t9_idx_decay > 0 ? 
    solve_complex(fmax(y_hat[t9_idx_decay][1], 1e-6),
                  fmax(y_hat[t9_idx_decay][2], 1e-6),
                  Kd_sRNA) : 0;
  real fraction_bound_t9 = t9_idx_decay > 0 ?
    C_wt_t9 / fmax(y_hat[t9_idx_decay][1], 1e-6) : 0;
  
  // 7. Steady-state predictions (synthesis / decay)
  real R_predicted_ss_t0 = syn_R_total_t0 / decay_eff_wt_t0;
  real R_predicted_ss_t9 = syn_R_total_t9 / decay_eff_wt_t9;
  real FC_predicted_ss = R_predicted_ss_t9 / R_predicted_ss_t0;
  
   
  // 8. sRNA levels at t=0 vs t=9
  real S_t0 = S0_wt;
  real S_t9 = t9_idx_decay > 0 ? y_hat[t9_idx_decay][2] : 0;
  real sRNA_FC = S_t9 / fmax(S_t0, 1e-6);
}
