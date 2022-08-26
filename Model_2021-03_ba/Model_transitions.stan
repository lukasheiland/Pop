functions {

}


data {
  
  int<lower=0> L;
  int<lower=0> N_species;
  
  array[L] int<lower=0> rep_species;
  vector<lower=0>[L] y_base;
  array[L] int<lower=0> y_trans;  
  vector[L] area_log;
  
  real prior_rate;

}


transformed data {
  
  vector[L] y_base_log = log(y_base);
  
}


parameters {
  
  vector[N_species] rate_log;
  
}


transformed parameters {
  
  vector<lower=0>[L] y_hat = exp(y_base_log + rate_log[rep_species] + area_log);

}


model {

  //// Priors
  rate_log ~ normal(prior_rate, 3);
  
  //// Model
  y_trans ~ poisson(y_hat);

}

generated quantities {

  // array[L] int y_sim = neg_binomial_2_rng(y_hat, phi_trans_rep);
  array[L] int y_sim = poisson_rng(y_hat);

}
