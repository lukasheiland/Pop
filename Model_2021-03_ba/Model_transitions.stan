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
  real mu_global;
  real<lower=0> sigma_global;
  vector[N_species] prop_logit;
}

transformed parameters {
  // real<lower=0> kappa = inv(kappa_inv);
}

model {
  //// Priors not necessary for regularization
  // mu ~ normal(0.01, 0.05);
  // kappa_inv ~ normal(0, 0.1);
  
  // Prior
  sigma_global ~ normal(0, 0.001);
  
  // Global level
  prop_logit ~ normal(mu_global, sigma_global);
  
  // Species level
  y_trans ~ binomial_logit(y_base, prop_logit[rep_species]);
}
