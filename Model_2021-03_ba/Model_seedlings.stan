
data {

  int<lower=1> N;
  array[N] int<lower=0> y;
  
  int<lower=1> N_offset;
  array[N] int<lower=1> rep_offset;
  vector[N] offset;
  
  vector[N] ba;
  vector<lower=0>[N] l_smooth;
  // vector[N] offset_scaled;
  // vector[N] ba_sum;
  
}


parameters {
  
  real l_log;
  real r_log;

  // real theta_logit;
  vector[N_offset] theta_logit;
  real m_logit;
  
  vector<lower=0>[N_offset] phi_inv_sqrt;
}

transformed parameters {

  vector<lower=0>[N_offset] phi = inv_square(phi_inv_sqrt);

  vector<lower=0>[N] y_hat_ha = exp(l_log) * l_smooth + exp(r_log) * ba;
  vector<lower=0>[N] y_hat =  y_hat_ha .* offset;
  
  real m_logit_scaled = m_logit * 1e-4; // account for the fact, that y_hat is on the hectare scale, so that there are huge numbers. 1e-4 converts to the m^2 scale
  
  vector<lower=0, upper=1>[N] prob_0 = inv_logit(theta_logit[rep_offset] + m_logit_scaled * y_hat_ha);
}


model {

  //// Priors
  theta_logit ~ normal(0, 2); // beta(2, 4);
  m_logit ~ std_normal();
  
  phi_inv_sqrt ~ std_normal();
  
  l_log ~ normal(0, 3);
  r_log ~ normal(0, 3);
  

  //// Model
  for (n in 1:N) {
   if (y[n] == 0)
      1 ~ bernoulli(prob_0[n]);
    else {
      0 ~ bernoulli(prob_0[n]);
      y[n] ~ neg_binomial_2(y_hat[n], phi[rep_offset[n]]) T[1, ];
    }
   }
   
}


generated quantities {

  // Beware! Quick and dirty, not-completely correct implementation.
  array[N] int y_sim;
  
  for (n in 1:N) {
  	if(y[n] == 0) {
	  	y_sim[n] = 0;
  	} else {
  		y_sim[n] = neg_binomial_2_rng(y_hat[n], phi[rep_offset[n]]); // 
  		while(y_sim[n] == 0) y_sim[n] = neg_binomial_2_rng(y_hat[n], phi[rep_offset[n]]); // 
  	}
  }

}
