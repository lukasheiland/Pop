functions {
  // Implementation of negbinomial probability density with zero inflation
  real neg_binomial_0_lpmf(int y, real y_hat, real phi_obs, real theta) {
  
   real t; // target
   // From a process point of view this is just a negbinom model, with some values multiplied by zero.
   // rbinom(100, 1, prob = 0.2) * rnbinom(100, size = 1100, mu = 10)
  
   if (y == 0) {
     // Joint Likelihood of 0 coming from probability theta or negbinonial
   	t = log_sum_exp(bernoulli_lpmf(1 | theta),
                        bernoulli_lpmf(0 | theta) + neg_binomial_2_lpmf(y | y_hat, phi_obs));
   } else {
  // Joint Likelihood of 0 coming from probability theta_rep or negbinonial
   	t = bernoulli_lpmf(0 | theta) +  // log1m(theta) synonymous to bernoulli_lpmf(0 | theta_rep)?
   		neg_binomial_2_lpmf(y | y_hat, phi_obs);
   }
   return t; // target wich will get summed up at each run
  }
   
}


data {
  int<lower=0> L;
  int<lower=0> N_species;
  
  array[L] int<lower=0> rep_species;
  vector<lower=0>[L] y_base;
  /// array[L] int<lower=0> y_base;
  array[L] int<lower=0> y_trans;  
  vector[L] area_log;
  
  real prior_rate;

}


transformed data {
  //// Version with base population as data
  vector[L] y_base_log = log(y_base);
}


parameters {
  //// Hierarchical version
  // real rate_global;
  // real<lower=0> sigma_raw;
  // vector[N_species] rate_contrast;
  
  //// ZI-version
  // real<lower=0,upper=1> theta;
  
  //// Version with latent base pop
  // vector<upper=17>[L] y_base_log; // max is log(15)
  
  vector[N_species] rate_log;
  
  // vector<lower=0>[N_species] phi_base_inv;
  // vector<lower=0>[N_species] phi_trans_inv;
  
}


transformed parameters {
  
  //// Hierarchical version
  // vector[N_species] rate_log = rate_global + rate_contrast * sigma_raw; // equivalent to rate_log ~ normal(mu_global, sigma_raw);
  
  //// ZI-version
  // vector[L] theta_rep = rep_vector(theta, L);
  
  // vector<lower=0>[N_species] phi_base = inv(phi_base_inv);
  // vector<lower=0>[N_species] phi_trans = inv(phi_trans_inv);
  
  vector<lower=0>[L] y_hat = exp(y_base_log + rate_log[rep_species] + area_log);
  
  // vector<lower=0>[L] phi_base_rep = phi_base[rep_species];
  // vector<lower=0>[L] phi_trans_rep = phi_trans[rep_species];

}


model {

  //// Priors
  // theta ~ beta(1, 20);
  
  // phi_base_inv ~ normal(0, 5);
  // phi_trans_inv ~ normal(0, 5);
  
  /// y_base_log ~ normal(10, 10);
  
  rate_log ~ normal(prior_rate, 3);
  
  //// Hierarchical version
  /// Hyperpriors
  // sigma_raw ~ normal(0, 0.1);
  // rate_contrast ~ normal(0, 0.1); /// Global level, non-centered
  
  //// ZI-version
  // for (l in 1:L) {
  //  y_trans[l] ~ neg_binomial_0(y_hat[l], phi_rep[l], theta);
  //}
  
  //// Version with latent base pop
  /// y_base ~ poisson(exp(y_base_log));
  // y_base ~ neg_binomial_2(exp(y_base_log), phi_base_rep);
  
  y_trans ~ poisson(y_hat);
  // y_trans ~ neg_binomial_2(y_hat, phi_trans_rep);
  
}

generated quantities {

  // array[L] int y_sim = neg_binomial_2_rng(y_hat, phi_trans_rep);
  array[L] int y_sim = poisson_rng(y_hat);

}
