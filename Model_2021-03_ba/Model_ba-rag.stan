functions {
  
  //// Difference equations
  matrix simulate(vector initialstate, int time_max,
                  vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                  vector ba_a_avg, real ba_a_upper,
                  int N_spec, int N_pops,
                  int[] i_j, int[] i_a, int[] i_b) {
    
    // State matrix with [species, times]. Times columns is sensible in order to have col-major access in matrices and for to_vector() later on
    matrix[N_pops, time_max] State;
    State[,1] = initialstate;
    
    for (t in 2:time_max) {
      // Structure of state[N_pops]: stage/species
      
      vector[N_spec] J = State[i_j, t-1];
      vector[N_spec] A = State[i_a, t-1];
      vector[N_spec] B = State[i_b, t-1];

      vector[N_spec] BA = A .* ba_a_avg + B;
      real BA_sum = sum(BA);
      
      /// Ricker model: assign to the the log states
      // Note: log1p(expm1(a) + exp(b)) == log(exp(a) + exp(b)); It is important to expm1() some state (here J), because the rates are positive anyway
      State[i_j, t]  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      State[i_a, t]  =  (g .* J + (A - h .*A )) ./ (1 + c_a*BA_sum);
      State[i_b, t]  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
    
    }
    
    return State;
  }
  
  
  //// ODE integration and data assignment to one long vector is wrapped into a function here because
  // - the transformed parameters block does not allow declaring integers (necessary for ragged data structure indexing),
  // - the model block does not alllow for declaring variables within a loop.
  vector unpack(vector state_init, int[] time_max, int[] times, int[] species, int[] pops,
                matrix g, matrix r, matrix s,  matrix l, // env-dependent on log scale
                vector b, vector c_j, vector c_b, vector h, // env-independent, on state model scale
                vector ba_a_avg, real ba_a_upper,
                int[] n_species, int[] n_pops, int[] n_reobs, int[] n_yhat,
                int L_yhat, int N_locs) {

    int pos_species = 1;
    int pos_pops = 1;
    int pos_times = 1; // segmenting times[L_y]
    int pos_yhat = 1;
    vector[L_yhat] y_hat;

    for (loc in 1:N_locs) {
      int n_s = n_species[loc];
      int n_p = n_pops[loc];
      int n_t = n_reobs[loc];
      
      // Slicing the global parameter vectors down to local level present species.
      int s_species[n_s] = segment(species, pos_species, n_s);
      
      // Getting indices of stages for the state vector
      int i_j[n_s] = segment(pops, pos_pops, n_s);
      int i_a[n_s] = segment(pops, pos_pops+n_s, n_s);
      int i_b[n_s] = segment(pops, pos_pops+n_s+n_s, n_s);

      //// Returns an matrix State[p, t]
      // print("r", exp(r[loc, s_species])'); print("h", h[loc, s_species]);
      matrix[n_p, time_max[loc]] States =
                        simulate(segment(state_init, pos_pops, n_p), // segment has length n_pops. Structure state_init: locs/stage/species 
                                 time_max[loc],
                                 exp(g[loc, s_species]'), exp(r[loc, s_species]'), exp(s[loc, s_species]'), exp(l[loc, s_species]'), // link function: r, loc, g, m are still on the log scale
                                 b[s_species], c_j[s_species], c_b[s_species], h[s_species],  // environmentally-independent parameters have been transformed through lognormal
                                 ba_a_avg, ba_a_upper,
                                 n_s, n_p,
                                 i_j, i_a, i_b);
      
  
      // Flattening the matrix into a vector for the location and append it to y_hat local vector yhat[m], and then into function-global vector y_hat[L_yhat].
      // to_vector converts matrix to a column vector in column-major order.
      y_hat[pos_yhat:(pos_yhat - 1 + n_yhat[loc])] =
                       to_vector(States[ , segment(times, pos_times, n_t)]); // only select columns with times in the data
      
      
      pos_species = pos_species + n_s;
      pos_pops = pos_pops + n_p;
      pos_times = pos_times + n_t;
      pos_yhat = pos_yhat + n_yhat[loc];
    }

  return y_hat; // Structure: locations/resurveys/pops(==stages/species)
  }

}

data {
  
  //// On assumed data structure
  // The most comprehensive data set is L_y (resp. L_y0) with grouping.
  // Everything is subset from the master subsets with these groupings (*_reobs, *_y0) and thus consistently sorted.
  // Ratio: "resurveys" includes "pops" due of the return structure in ode_*(); "pops" includes "plots" because of the loc == population assumption (randowm effect from location/resurvey/pop to .../plots).
  // (With population interactions and ODE integration, n_species, n_stages, and n_times have to be constant within the process model level, i.e. here location.)
  // Factor "pops" is structured stages/species.
  
  //// N — number of observations/groups
  int<lower=0> N_locs; // overall number of locations. This is the major running variable to iterate over the model vectors.
  int<lower=0> N_species; // overall number of unique species across plot and locations (not nested!)
  int<lower=0> N_beta;
  
  int<lower=0> L_times; // locations/reobs
  int<lower=0> L_init; // locations/pops
  int<lower=0> L_yhat; // locations/resurveys/pops
  int<lower=0> L_y0; // locations/pops/plots
  int<lower=0> L_y; // locations/resurveys/pops/plots
  int<lower=0> L_species; // locations/species — length of the vector species



  //// n - number of levels within locs
  int<lower=0> n_species[N_locs]; // n species within loc; this is the length of initial values!
  int<lower=0> n_pops[N_locs]; // n pops (species*stages) within loc; this is the length of initial values!
  int<lower=0> n_reobs[N_locs]; // n solution times within locs
  int<lower=0> n_yhat[N_locs]; // n pops*reobs within locs

  //// rep - repeat indices within groups for broadcasting to the subsequent hierarchical levels
  int<lower=1> rep_init2y0[L_y0]; // repeat predictions on level "locations/pops" n_plots times to "locations/pops/plots"
  int<lower=1> rep_yhat2y[L_y]; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  int<lower=1> rep_locs2times[L_times]; // repeat locations n_times per location to "locations/pops/resurveys/plots"

  
  //// actual data
  int<lower=1> species[L_species]; // locations/species
  int<lower=1> pops[L_init]; // locations/pops: just numbering, e.g 1:12 for indexing

  // obsmethod — factor (1, 2)
  int<lower=1> rep_obsmethod2y0[L_y0]; // repeat predictions on level "locations/pops" n_plots times to "locations/pops/plots"
  int<lower=1> rep_obsmethod2y[L_y]; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  
  
  /// FIX by getting data block from ba
  
  real dbh_lower_a;
  real dbh_lower_b;
  
  int time_init[N_locs];
  int time_max_data[N_locs];
  int times_data[L_times]; // locations/resurveys // old: int times[L_y]; // locations/plots/pops/resurveys

  matrix[N_locs, N_beta] X; // design matrix
  
  vector[L_y0] y0_log;
  vector[L_y] y_log;
  
}

transformed data {

  for (loc in 1:N_locs) time_max[loc] = time_max_data[loc] - time_init[loc] + 1; // no vector operations for arrays
  for (t in 1:L_times) times[t] = times_data[t] - time_init[rep_locs2times[t]] + 1;
  
}


parameters {
  //// Level 1, (global, species): species-specific rates, indepndent of environment.
  // … dependent on environment.
  matrix[N_beta, N_species] Beta_g; // (J-, A+) transition rate from J to A
  // matrix[N_beta, N_species] Beta_m_j; // (J-) density independent mortalitty
  matrix[N_beta, N_species] Beta_r; // // (J+) flow into the system, dependent on env
  matrix[N_beta, N_species] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A
  matrix[N_beta, N_species] Beta_l; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A


  // … independent of environment
  vector<lower=0>[N_species] b;
  vector<lower=0>[N_species] c_j;
  // vector<lower=0>[N_species] c_a;
  vector<lower=0>[N_species] c_b;
  vector<lower=0>[N_species] h; // (A-, B+), here still unconstrained on log scale, flow out of the system, independent of environment ("intercept only")
  // vector<lower=0>[N_species] m_a;
  
  vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  vector<lower=0>[2] sigma_obs; // lognormal error for observations from predictions

  vector[L_init] state_init_log; // log scale
    // remember to change getTrueInits() when choosing a rectangular version as init parameters
}


transformed parameters {

  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  matrix<lower=0, upper = 1>[N_locs, N_species] g_log = X * Beta_g;
  // matrix[N_locs, N_species] m_j_log = X * Beta_m_j;
  matrix[N_locs, N_species] r_log = X * Beta_r;
  matrix[N_locs, N_species] s_log = X * Beta_s;
  matrix[N_locs, N_species] l_log = X * Beta_l;


}

model {

  //---------- PRIORS ---------------------------------
  
  // Beta_r[1,] ~ normal(2, 2); // intercept
  // sigma_process ~ normal(0, 0.01);
  // sigma_obs ~ normal(0, [2, 0.2]); // for observations from predictions

  //---------- MODEL ---------------------------------


  // Separately fit initial state to feed to integrator.
  y0_log ~ normal(state_init_log[rep_init2y0], sigma_obs[rep_obsmethod2y0]); //  y0 ~ neg_binomial_0(state_init[rep_init2y0], theta, sigma_obs);
  

  // Level 2 (locs) to level 3 y.
  vector[L_yhat] y_hat_log = unpack(state_init_log, time_max, times, species, pops,
                                    g_log, r_log, s_log, l_log, // rates matrix[N_locs, N_species]; will have to be transformed
                                    b, c_j, c_b, h, // rates vector[N_species]
                                    ba_a_avg, ba_a_upper,
                                    n_species, n_pops, n_reobs, n_yhat,
                                    L_yhat, N_locs);
  
  
  // Fit predictions to data. (level: location/plot/resurvey/species)
  y_log ~ normal(y_hat_log[rep_yhat2y], sigma_obs[rep_obsmethod2y]); //   y ~ neg_binomial_0(y_hat[rep_yhat2y], theta, sigma_obs);
}


generated quantities {

  real y_log_sim[L_y];
  vector[L_y] y_hat_log_rep;
  
  vector[L_yhat] y_hat_log = unpack(state_init_log, time_max, times, species, pops,
                                    g_log, r_log, s_log, l_log, // rates matrix[N_locs, N_species]; will have to be transformed
                                    b, c_j, c_b, h, // rates vector[N_species]
                                    ba_a_avg, ba_a_upper,
                                    n_species, n_pops, n_reobs, n_yhat,
                                    L_yhat, N_locs);
  
  y_hat_log_rep = y_hat_log[rep_yhat2y];
  y_log_sim = normal_rng(y_hat_log_rep, sigma_obs[rep_obsmethod2y]);

}

