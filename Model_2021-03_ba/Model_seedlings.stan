
data {

  int<lower=1> N;
  int<lower=1> N_plots;
  array[N] int<lower=0> y;
  
  int<lower=1> N_offset;
  array[N] int<lower=1> rep_offset;
  array[N] int<lower=1> rep_species;
  array[N] int<lower=1> rep_plot;
  vector[N] offset_data;
  
  vector[N] ba;
  vector<lower=0>[N] l_smooth;
  
}


parameters {
  
  vector[2] k_log;
  vector[2] r_log;

  // real theta_logit;
  vector[N_offset] theta_logit;
  real m_logit;
  
  vector<lower=0>[N_offset] phi_inv_sqrt;
}

transformed parameters {

  vector<lower=0>[N_offset] phi = inv_square(phi_inv_sqrt);
  vector<lower=0>[N_plots] y_hat_ha_total = rep_vector(0, N_plots);
  vector[N] y_hat_ha_total_centered;

  vector<lower=0>[N] y_hat_ha = exp(k_log[rep_species]) + exp(r_log[rep_species]) .* ba;
  vector<lower=0>[N] y_hat =  y_hat_ha .* offset_data;
  
  for (n in 1:N) {
    y_hat_ha_total[rep_plot[n]] += y_hat_ha[n];
  }
  
  real center_y_hat_ha_total = mean(y_hat_ha_total);
  
  for (n in 1:N) {
    y_hat_ha_total_centered[n] = y_hat_ha_total[rep_plot[n]] - center_y_hat_ha_total;
  }

  real m_logit_scaled = m_logit * 1e-4; // account for the fact, that y_hat is on the hectare scale, so that there are huge numbers. 1e-4 converts to the m^2 scale

  
  vector<lower=0, upper=1>[N] prob_0 = inv_logit(theta_logit[rep_offset] + m_logit_scaled * y_hat_ha_total_centered);

}


model {

  //// Priors
  theta_logit ~ normal(0, 2); // beta(2, 4);
  m_logit ~ std_normal();
  
  phi_inv_sqrt ~ std_normal();
  
  k_log ~ normal(4, 5);
  // l_log ~ normal(0, 3);
  r_log ~ normal(3, 5);
  

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

  // Beware! Quick and dirty, incorrect implementation of y_sim, only for rough residuals
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
