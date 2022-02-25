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
  //// For hierarchical version
  // real rate_global;
  // real<lower=0> sigma_raw;
  // vector[N_species] rate_contrast;
  
  vector[N_species] rate_log;
}


transformed parameters {
  
  //// For hierarchical version
  // vector[N_species] rate_log = rate_global + rate_contrast * sigma_raw; // equivalent to rate_log ~ normal(mu_global, sigma_raw);
  
  vector[L] y_hat_log = y_base_log + rate_log[rep_species] + area_log; 
}


model {  
  
  //// For hierarchical version
  /// Hyperpriors
  // sigma_raw ~ normal(0, 0.1);
  // rate_contrast ~ normal(0, 0.1); /// Global level, non-centered

  y_trans ~ poisson_log(y_hat_log);
}
