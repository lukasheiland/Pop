
data {

  int<lower=1> N;
  array[N] int<lower=0> y;
  
  int<lower=1> N_offset;
  array[N] int<lower=1> rep_offset;
  vector[N] offset;
  
  vector[N] ba;
  // vector[N] l_smooth;
  // vector[N] ba_sum;
  // vector[N] offset_scaled;
}


parameters {
  
  real k_log;
  // real l_log;
  real r_log;
  
  // real o_log; // slope for offset
  vector[N_offset] o_log;
  real<lower=0> sigma_o;

  vector<lower=0, upper=1>[N_offset] theta;
  vector<lower=0>[N_offset] phi_inv_sqrt;  
}

transformed parameters {

      vector<lower=0>[N_offset] phi = inv_square(phi_inv_sqrt);
      
      //// k
      vector<lower=0>[N] y_hat_ha = exp(k_log + o_log[rep_offset]) + exp(r_log) * ba;
      
      vector<lower=0>[N] y_hat =  y_hat_ha .* offset;
}


model {

  //// Priors
  theta ~ beta(2, 2);
  phi_inv_sqrt ~ std_normal();
  
  k_log ~ normal(0, 2);
  // l_log ~ normal(0, 2);
  r_log ~ normal(0, 2);
  // s_log ~ normal(0, 2);

  sigma_o ~ std_normal();
  o_log ~ normal(0, sigma_o);
  
 for (n in 1:N) {
   if (y[n] == 0)
     target += log_sum_exp(bernoulli_lpmf(1 | theta[rep_offset[n]]),
                           bernoulli_lpmf(0 | theta[rep_offset[n]])
                             + poisson_lpmf(y[n] | y_hat[n]));
   else
     target += bernoulli_lpmf(0 | theta[rep_offset[n]])
                 + neg_binomial_2_lpmf(y[n] | y_hat[n], phi[rep_offset[n]]);
  }
  
  //// version without zi
  // y ~ neg_binomial_2(y_hat, phi[rep_offset]);
  
}


generated quantities {

  array[N] int y_sim;
  
  for (n in 1:N) {
    y_sim[n] = neg_binomial_2_rng(y_hat[n], phi[rep_offset[n]]) * !bernoulli_rng(theta[rep_offset[n]]);
  }

  //// version without zi
 //  y_sim = neg_binomial_2_rng(y_hat, phi[rep_offset]);
  
}
