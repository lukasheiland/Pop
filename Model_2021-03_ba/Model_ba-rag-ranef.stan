functions {
  
  //// Difference equations
  matrix simulate(vector initialstate, int time_max,
                  vector g, vector m_j, vector r, vector s,
                  vector b, vector c_j, vector c_a, vector c_b, vector h, vector m_a,
                  vector[] u,
                  real ba_a_avg, real ba_a_upper,
                  int N_spec, int N_pops,
                  int[] i_j, int[] i_a, int[] i_b) {
    
    // State matrix with [species, times]. Times columns is sensible in order to have col-major access in matrices and for to_vector()
    matrix[N_pops, time_max] State;
    State[,1] = initialstate;
    
    for (t in 2:time_max) {
      // Structure of state[N_pops]: stage/species
      
      vector[N_spec] J_log = State[i_j, t-1];
      vector[N_spec] A_log = State[i_a, t-1];

      vector[N_spec] J = exp(J_log);
      vector[N_spec] A = exp(A_log);
      vector[N_spec] B = exp(State[i_b, t-1]);
      
      real BA = sum(A*ba_a_avg + B);
      
      /// Ricker model: fill the log states
      // Note: log1p(expm1(a) + exp(b)) == log(exp(a) + exp(b))
      State[i_j, t]  =  log1p(expm1(J_log) + r) - c_j*sum(J) - s*BA - m_j - g + u[1, ]; // equivalent to … =  J + (r - (c_j*sum(J) + s*AB + g + m_j).*J + u[1, i_j]) * dt;
      State[i_a, t]  =  log1p(A + expm1(J_log - g)) - c_a*BA - m_a - h + u[2, ];
      State[i_b, t]  =  log1p(B + expm1(A_log - h) * ba_a_upper) + b - c_b*sum(B) + u[3, ]; // exp(A_log - h) * ba_a_upper is the input of new basal area through ingrowth of count A*exp(-h) from A
    
    }
    
    return State;
  }
  
  
  //// ODE integration and data assignment to one long vector is wrapped into a function here because
  // - the transformed parameters block does not allow declaring integers (necessary for ragged data structure indexing),
  // - the model block does not alllow for declaring variables within a loop.
  vector unpack(vector state_init, int[] time_max, int[] times, int[] species, int[] pops,
                matrix g,  matrix m_j, matrix r, matrix s,  // env-dependent on log scale
                vector[] b, vector[] c_j, vector[] c_a, vector[] c_b, vector[] h, vector[] m_a, // env-independent, on state model scale
                vector[] u,
                real ba_a_avg, real ba_a_upper,
                int[] n_species, int[] n_pops, int[] n_reobs, int[] n_yhat,
                int N_yhat, int N_locs) {

    int pos_species = 1;
    int pos_pops = 1;
    int pos_times = 1; // segmenting times[N_y]
    int pos_yhat = 1;
    vector[N_yhat] y_hat;

    for (l in 1:N_locs) {
      int n_s = n_species[l];
      int n_p = n_pops[l];
      int n_t = n_reobs[l];
      
      // Slicing the global parameter vectors down to local level present species.
      int s_species[n_s] = segment(species, pos_species, n_s);
      
      // Getting indices of stages for the state vector
      int i_j[n_s] = segment(pops, pos_pops, n_s);
      int i_a[n_s] = segment(pops, pos_pops+n_s, n_s);
      int i_b[n_s] = segment(pops, pos_pops+n_s+n_s, n_s);

      //// Returns an matrix State[p, t]
      // print("r", exp(r[l, s_species])'); print("h", h[l, s_species]);
      matrix[n_p, time_max[l]] States =
                        simulate(segment(state_init, pos_pops, n_p), // segment has length n_pops. Structure state_init: locs/stage/species 
                                 time_max[l],
                                 exp(g[l, s_species]'), exp(m_j[l, s_species]'), exp(r[l, s_species]'), exp(s[l, s_species]'), // link function: r, l, g, m are still on the log scale
                                 b[l, s_species], c_j[l, s_species], c_a[l, s_species], c_b[l, s_species], h[l, s_species], m_a[l, s_species],  // environmentally-independent parameters have been transformed through lognormal
                                 u[,s_species], // process error
                                 ba_a_avg, ba_a_upper,
                                 n_s, n_p,
                                 i_j, i_a, i_b);
      
  
      // Flattening the matrix into a vector for the location and append it to y_hat local vector yhat[m], and then into function-global vector y_hat[N_yhat].
      // to_vector converts matrix to a column vector in column-major order.
      y_hat[pos_yhat:(pos_yhat - 1 + n_yhat[l])] =
                       to_vector(States[ , segment(times, pos_times, n_t)]); // only select columns with times in the data
      
      
      // Debugging loop
      // int containsinf = 0;
      // 
      // for (i in pos_yhat:(pos_yhat - 1 + n_yhat[l])) {
      //   if (y_hat[i] == negative_infinity() ) containsinf = 1;
      // }
      // 
      // if(containsinf) {
      //   print("states:", y_hat[pos_yhat:(pos_yhat - 1 + n_yhat[l])]);
      // }
      
      
      pos_species = pos_species + n_s;
      pos_pops = pos_pops + n_p;
      pos_times = pos_times + n_t;
      pos_yhat = pos_yhat + n_yhat[l];
    }

  return y_hat; // Structure: locations/resurveys/pops(==stages/species)
  }

}

data {
  
  //// On assumed data structure
  // The most comprehensive data set is N_y (resp. N_y0) with grouping.
  // Everything is subset from the master subsets with these groupings (*_reobs, *_y0) and thus consistently sorted.
  // Ratio: "resurveys" includes "pops" due of the return structure in ode_*(); "pops" includes "plots" because of the loc == population assumption (randowm effect from location/resurvey/pop to .../plots).
  // (With population interactions and ODE integration, n_species, n_stages, and n_times have to be constant within the process model level, i.e. here location.)
  // Factor "pops" is structured stages/species.
  
  //// N — number of observations/groups
  int<lower=0> N_locs; // overall number of locations. This is the major running variable to iterate over the model vectors.
  int<lower=0> N_times; // locations/reobs
  int<lower=0> N_init; // locations/pops
  int<lower=0> N_yhat; // locations/resurveys/pops
  int<lower=0> N_y0; // locations/pops/plots
  int<lower=0> N_y; // locations/resurveys/pops/plots
  int<lower=0> N_species; // locations/species — length of the vector species


  int<lower=0> N_totalspecies; // overall number of unique species across plot and locations (not nested!)
  int<lower=0> N_beta;

  //// n - number of levels within locs
  int<lower=0> n_species[N_locs]; // n species within loc; this is the length of initial values!
  int<lower=0> n_pops[N_locs]; // n pops (species*stages) within loc; this is the length of initial values!
  int<lower=0> n_reobs[N_locs]; // n solution times within locs
  int<lower=0> n_yhat[N_locs]; // n pops*reobs within locs

  //// rep - repeat indices within groups for broadcasting to the subsequent hierarchical levels
  int<lower=1> rep_init2y0[N_y0]; // repeat predictions on level "locations/pops" n_plots times to "locations/pops/plots"
  int<lower=1> rep_yhat2y[N_y]; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  int<lower=1> rep_locs2times[N_times]; // repeat locations n_times per location to "locations/pops/resurveys/plots"

  
  //// actual data
  int<lower=1> species[N_species]; // locations/species
  int<lower=1> pops[N_init]; // locations/pops: just numbering, e.g 1:12 for indexing

  // obsmethod — factor (1, 2)
  int<lower=1> obsmethod_y0[N_y0]; // repeat predictions on level "locations/pops" n_plots times to "locations/pops/plots"
  int<lower=1> obsmethod_y[N_y]; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  
  real dbh_lower_a;
  real dbh_lower_b;
  
  int time_init[N_locs];
  int time_max_data[N_locs];
  int times_data[N_times]; // locations/resurveys // old: int times[N_y]; // locations/plots/pops/resurveys

  matrix[N_locs, N_beta] X; // design matrix
  
  vector[N_y0] y0_log;
  vector[N_y] y_log;
  
}

transformed data {
  // Shift all times to start at 1.
  int time_max[N_locs];
  int times[N_times];
  real ba_a_upper = pi() * (dbh_lower_b/2)^2 * 1e-6; // # pi*r^2, mm^2 to m^2
  real ba_a_avg = pi() * ((dbh_lower_a + dbh_lower_b)/2/2)^2 * 1e-6;

  for (l in 1:N_locs) time_max[l] = time_max_data[l] - time_init[l] + 1; // no vector operations for arrays
  for (t in 1:N_times) times[t] = times_data[t] - time_init[rep_locs2times[t]] + 1;
  
}


parameters {
  //// Level 1, (global, species): species-specific rates, indepndent of environment.
  // … dependent on environment.
  matrix[N_beta, N_totalspecies] Beta_g; // (J-, A+) transition rate from J to A
  matrix[N_beta, N_totalspecies] Beta_m_j; // (J-) density independent mortalitty
  matrix[N_beta, N_totalspecies] Beta_r; // // (J+) flow into the system, dependent on env
  matrix[N_beta, N_totalspecies] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A

  // … independent of environment
  vector<lower=0>[N_totalspecies] b;
  vector<lower=0>[N_totalspecies] c_j;
  vector<lower=0>[N_totalspecies] c_a;
  vector<lower=0>[N_totalspecies] c_b;
  vector<lower=0>[N_totalspecies] h; // (A-, B+), here still unconstrained on log scale, flow out of the system, independent of environment ("intercept only")
  vector<lower=0>[N_totalspecies] m_a;
  
  vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  vector<lower=0>[2] sigma_obs; // lognormal error for observations from predictions
  real<lower=0> shape_par; // normal error for time plot (plots) from locations
  vector[N_totalspecies] u[3]; // normal error for time plot (plots) from locations


  //// Level 2 (locations):
  // "random effect" h rate
  vector<lower=0>[N_totalspecies] b_loc[N_locs];
  vector<lower=0>[N_totalspecies] c_j_loc[N_locs];
  vector<lower=0>[N_totalspecies] c_a_loc[N_locs];
  vector<lower=0>[N_totalspecies] c_b_loc[N_locs];
  vector<lower=0>[N_totalspecies] h_loc[N_locs];
  vector<lower=0>[N_totalspecies] m_a_loc[N_locs];

  
  vector[N_init] state_init_log; // log scale
}


transformed parameters {

  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  matrix[N_locs, N_totalspecies] g_log = X * Beta_g;
  matrix[N_locs, N_totalspecies] m_j_log = X * Beta_m_j;
  matrix[N_locs, N_totalspecies] r_log = X * Beta_r;
  matrix[N_locs, N_totalspecies] s_log = X * Beta_s;

}

model {

  //---------- PRIORS ---------------------------------
  
  // Beta_r[1,] ~ normal(2, 2); // intercept
  // 1 ./ shape_par ~ normal(0, 100); 
  // sigma_process ~ normal(0, 0.01);
  sigma_obs ~ gamma(10, 10/2); // for observations from predictionsions

  //---------- MODEL ---------------------------------
  // debugging: print("shape_par ", shape_par); print("r_loc ", r_loc[l,]);
  
  //// Level 1 (species) to 2 (locs). Transition rate from A to B, directly, without transformation through env, from level 1 to 2.
  for (l in 1:N_locs) {
    // r[l,]    ~ lognormal(r_log, shape_par); // lognormal acts like exp (log-link) on the logmedian
    b_loc[l,]    ~ gamma(shape_par, shape_par./b);
    c_j_loc[l,]  ~ gamma(shape_par, shape_par./c_j);
    c_a_loc[l,]  ~ gamma(shape_par, shape_par./c_a);
    c_b_loc[l,]  ~ gamma(shape_par, shape_par./c_b);
    h_loc[l,]    ~ gamma(shape_par, shape_par./h);
    m_a_loc[l,]  ~ gamma(shape_par, shape_par./m_a);

  }
  
  //// Level 2 (locs) to 3 (plots).
  // Separately fit initial state to feed to integrator.
  y0_log ~ normal(state_init_log[rep_init2y0], sigma_obs[obsmethod_y0]); //  y0 ~ neg_binomial_0(state_init[rep_init2y0], theta, sigma_obs);
  
  // Process error
  for (stage in 1:3) {
      u[stage,] ~ normal(0, sigma_process[stage]);
  }
  
  // Level 2 (locs).
  vector[N_yhat] y_hat_log = unpack(state_init_log, time_max, times, species, pops,
                                    g_log, m_j_log, r_log, s_log, // rates matrix[N_locs, N_totalspecies]; will have to be transformed
                                    b_loc, c_j_loc, c_a_loc, c_b_loc, h_loc, m_a_loc, // rates vector[N_totalspecies][N_locs]
                                    u, // process error
                                    ba_a_avg, ba_a_upper,
                                    n_species, n_pops, n_reobs, n_yhat,
                                    N_yhat, N_locs);
  
  
  // Fit predictions to data. (level: location/plot/resurvey/species)
  y_log ~ normal(y_hat_log[rep_yhat2y], sigma_obs[obsmethod_y]); //   y ~ neg_binomial_0(y_hat[rep_yhat2y], theta, sigma_obs);
}
