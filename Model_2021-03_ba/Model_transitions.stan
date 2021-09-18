data {
  int<lower=0> L;
  int<lower=0> N_species;
  
  array[L] int<lower=0> rep_species;
  array[L] int<lower=0> y_base;
  array[L] int<lower=0> y_trans;
}

parameters {
  // real<lower=0, upper=1> mu;
  // real<lower=0> kappa_inv;
  // vector<lower=0, upper=1>[N_species] prop;
  real prop_global;
  real<lower=0> sigma_raw;
  vector[N_species] prop_contrast;
}

transformed parameters {
  // real<lower=0> kappa = inv(kappa_inv);
  real sigma_global = sigma_raw * 0.01;
  vector[N_species] prop_logit = prop_global + prop_contrast * sigma_global;
}

model {
  //// Priors not necessary for regularization
  // mu ~ normal(0.01, 0.05);
  // kappa_inv ~ normal(0, 0.1);
  
  // Prior
  sigma_raw ~ std_normal(); // sigma_global ~ normal(0, 0.5);

  // Global level, non-centered
  prop_contrast ~ std_normal();
  //// Equivalent to:
  // prop_logit ~ normal(mu_global, sigma_global);
  
  // Species level
  y_trans ~ binomial_logit(y_base, prop_logit[rep_species]);
}
