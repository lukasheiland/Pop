functions {
  
  //// Difference equations
  vector[] simulate(vector initialstate, int time_max, int[] times,
                  vector g, vector r, vector s, vector l,
                  vector b, vector c_j, vector c_b, vector h,
                  real ba_a_avg, real ba_a_upper,
                  int N_spec, int N_pops,
                  // matrix u, // a matrix[N_pops, time_max-1]
                  int[] i_j, int[] i_a, int[] i_b) {
    
    
    // State matrix with [species, times].
    vector[N_pops] State[time_max];
    State[1,] = initialstate;
    // print("State 1 before: ", State[1,]);

    for (t in 2:time_max) {
      // Structure of state[N_pops]: stage/species
      // State has times rows (array row-major),  u has times cols (matrix: col-major access)
      
      vector[N_spec] J = State[t-1, i_j];
      vector[N_spec] A = State[t-1, i_a];
      vector[N_spec] B = State[t-1, i_b];

      vector[N_spec] BA = A*ba_a_avg + B;
      real BA_sum = sum(BA);

      State[t, i_j]  =  (r.*BA + l + (J - g.*J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      State[t, i_a]  =  (g.*J + (A - h.*A)) ./ (1 + c_b*BA_sum);
      State[t, i_b]  =  (1+b).*((h.*A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
      
    }
    
    // print("State 1 after: ", State[1,]);

    return State[times,];
  }
  
}



data {
  
  //// N — number of observations/groups
  // full rectangular structure, i.e. same numner of times, species, pops, plots within locs
  int<lower=0> N_locs; // overall number of locations. This is the major running variable to iterate over the model vectors.
  int<lower=0> N_plots;
  int<lower=0> N_times;
  int<lower=0> N_species;
  int<lower=0> N_pops; // populations are species*stages
  
  int<lower=0> N_beta;
  
  //// actual data
  // int<lower=1> species[N_species]; // not needed?
  // int<lower=1> pops[N_pops];

  // obsmethod — factor (1, 2)
  int<lower=1> rep_obsmethod2pops[N_pops]; // 1 1 1 1 2 2 2 2 2 2 2 2

  int<lower=1> i_j[N_species]; // e.g. 1:4
  int<lower=1> i_a[N_species]; // e.g. 5:9, etc.
  int<lower=1> i_b[N_species];

  real dbh_lower_a; // 100
  real dbh_lower_b; // 200
  
  int times[N_locs, N_times]; // assumes start at 1. (Even though the number of times is the same, the times can deviate)
  int time_max[N_locs];
  int timespan_max; // max(time_max) - time_globalmin

  matrix[N_locs, N_beta] X; // design matrix
  
  // The response.
  vector[N_pops] y0_loc [N_locs];
  vector[N_pops] y [N_locs, N_times, N_plots];
}


transformed data {
  // Shift all times to start at 1.
  real ba_a_upper = pi() * (dbh_lower_b/2)^2 * 1e-6; // # pi*r^2, mm^2 to m^2
  real ba_a_avg = pi() * ((dbh_lower_a + dbh_lower_b)/2/2)^2 * 1e-6;
  
  //// Data for separate fitting of the initial state
  // vector[N_pops] y0 [N_locs, N_plots] = y[ , , 1, ];
}


parameters {
  //// Level 1, (global, species): species-specific rates, indepndent of environment.
  // … dependent on environment.
  
  // matrix[N_beta, N_species] Beta_g; // (J-, A+) transition rate from J to A
  // matrix[N_beta, N_species] Beta_r; // // (J+) flow into the system, dependent on env
  // matrix[N_beta, N_species] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A
  matrix[N_beta, N_species] Beta_l;


  // … independent of environment
  vector<lower=0>[N_species] b;
  // vector<lower=0>[N_species] c_a;
  vector<lower=0>[N_species] c_b;
  vector<lower=0>[N_species] c_j;
  vector<lower=0,upper=1>[N_species] g;
  vector<lower=0,upper=1>[N_species] h; // (A-, B+), here still unconstrained on log scale, flow out of the system, independent of environment ("intercept only")
  vector<lower=0>[N_species] l;
  vector<lower=0>[N_species] r;
  vector<lower=0>[N_species] s;
  
  // vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  // vector<lower=0>[2] sigma_obs; // observation error
  vector<lower=0>[2] alpha_obs; // observation error
  
  // matrix[N_pops, timespan_max] u[N_locs];

  vector<lower=0>[N_pops] state_init[N_locs];
}


transformed parameters {
  
  vector[N_pops] y_hat[N_locs, N_times];
  
  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  // matrix[N_locs, N_species] g_logit = X * Beta_g;
  // matrix[N_locs, N_species] r_log = X * Beta_r;
  // matrix[N_locs, N_species] s_log = X * Beta_s;
  matrix[N_locs, N_species] l_log = X * Beta_l;

  
  for(loc in 1:N_locs) {
    
    y_hat[loc, ] = simulate(state_init[loc], time_max[loc], times[loc, ],
                              // exp(g_log[loc, ]'), exp(r_log[loc, ]'), exp(s_log[loc, ]'),
                              g, r, s,  exp(l_log[loc, ]'),
                              b, c_j, c_b, h,
                              ba_a_avg, ba_a_upper,
                              N_species, N_pops,
                              // u[loc],
                              i_j, i_a, i_b);
    }

}


model {
  
  //---------- PRIORS ---------------------------------
  
  // Beta_r[1,] ~ normal(2, 2); // intercept
  // 1 ./ shape_par ~ normal(0, 100); 
  // sigma_process ~ normal(0, 0.01);
  // sigma_obs ~ normal(0, [0.5, 0.1]); // for observations from predictions
  
  // r ~ normal(20, 4);
  
  // Beta_g[1,] ~ normal(-12, 1);
  // to_vector(Beta_g[2:N_beta,]) ~ std_normal();
  
  // h ~ normal(0.4, 1);
  // g ~ gamma(5, 5/0.007);
  // b ~ gamma(10, 10/0.4);
  
  //---------- MODEL ---------------------------------
  
  for(loc in 1:N_locs) {

    // Some priors
    // normal: // state_init[loc] ~ normal(y0_loc[loc], sigma_obs[rep_obsmethod2pops]);
    state_init[loc] ~ gamma(alpha_obs[rep_obsmethod2pops], alpha_obs[rep_obsmethod2pops] ./ y0_loc[loc]);
    
    // to_vector(u[loc]) ~ normal(0, 0.1);
    
    for (t in 1:N_times) {
      
      //// Debugging alternatives
      // print(y_hat_log[loc, p, ]);
      // y0_log[loc, p, ] ~ normal(state_init_log[loc, p, ], sigma_obs[rep_obsmethod2pops]);
        
  
      // for (t in 2:N_times) { // in case of separate y0_log fitting
      for (p in 1:N_plots) {
        
        // y_log[loc, t, p, ] ~ normal(y_hat_log[loc, t], sigma_obs[rep_obsmethod2pops]);
        y[loc, t, p, ] ~ gamma(alpha_obs[rep_obsmethod2pops], alpha_obs[rep_obsmethod2pops] ./ y_hat[loc, t]);

      }
    }
  }
}

generated quantities {

  real y_sim [N_locs, N_times, N_plots, N_pops];
  vector[N_pops] y_hat_rep[N_locs, N_times, N_plots];

  for(loc in 1:N_locs) {

    for (t in 1:N_times) {

      for (p in 1:N_plots) {

        y_hat_rep[loc, t, p,] = y_hat[loc, t,];
        
        // y_sim[loc, t, p, ] = normal_rng(y_hat[loc, t], sigma_obs[rep_obsmethod2pops]);
        y_sim[loc, t, p, ] = gamma_rng(alpha_obs[rep_obsmethod2pops], alpha_obs[rep_obsmethod2pops]./y_hat[loc, t]);



      }
    }
  }

}
