data {
  int<lower=0> L;
  int<lower=0> N_species;
  
  array[L] int<lower=0> rep_species;
  vector<lower=0>[L] y_base;
  array[L] int<lower=0> y_trans;  
  vector[L] area_log;

}


transformed data {
  vector[L] y_base_log = log(y_base);
}


parameters {
  // real<lower=0, upper=1> mu;
  // real<lower=0> kappa_inv;
  // vector<lower=0, upper=1>[N_species] prop;
  real rate_global;
  real<lower=0> sigma_raw;
  vector[N_species] rate_contrast;
}


transformed parameters {
  // real<lower=0> kappa = inv(kappa_inv);
  real sigma_global = sigma_raw * 0.1;
  vector[N_species] rate_log = rate_global + rate_contrast * sigma_global;
  vector[L] y_hat_log = y_base_log + rate_log[rep_species] + area_log;
}

model {
  //// Priors not necessary for regularization
  // mu ~ normal(0.01, 0.05);
  // kappa_inv ~ normal(0, 0.1);
  
  // Prior
  sigma_raw ~ std_normal(); // sigma_global ~ normal(0, 0.5);

  // Global level, non-centered
  rate_contrast ~ std_normal();
  //// Equivalent to:
  // rate_log ~ normal(mu_global, sigma_global);
  
  // Species level
  y_trans ~ poisson_log(y_hat_log);
}
