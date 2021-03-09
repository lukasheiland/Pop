functions {
  
  //// Difference equations
  matrix simulate(vector initialstate, int time_max,
                  vector r, vector s, vector g, vector m_j,
                  vector h, vector c_j, vector c_ab,
                  vector[] u,
                  int N_spec, int N_pops,
                  int[] i_j, int[] i_a, int[] i_b) {
    
    
    // int i_j[N_spec] = {1:N_spec};
    // int i_a[N_spec] = N_spec+(1:N_spec);
    // int i_b[N_spec] = (N_spec+N_spec+1):N_pops;
    
    // State matrix with [species, times]. Times columns is sensible in order to have col-major access in matrices and for to_vector()
    matrix[N_pops, time_max] State;
    State[,1] = initialstate;
    
    for (t in 2:time_max) {
      // Structure of state[N_pops]: stage/species

      vector[N_spec] J = State[i_j, t-1];
      vector[N_spec] A = State[i_a, t-1];
      vector[N_spec] B = State[i_b, t-1];
      
      real AB = sum(A + B);
      
      // Euler method with assumed constant time interval 1.
      vector[N_spec] J1  =  J + r - (c_j*sum(J) + s*AB + g + m_j).*J + u[1,]; // equivalent to … =  J + (r - (c_j*sum(J) + s*AB + g + m_j).*J + u[1, i_j]) * dt;
      vector[N_spec] A1  =  A + g .* J - (c_ab*AB + h).*A + u[2,];
      vector[N_spec] B1  =  B + h .* A - (c_ab*AB).*B + u[3,];
      
      // (supposedly) efficient vectorized way to set x<=0 to eps
      State[i_j, t] = (fabs(J1) + J1)*0.5 + machine_precision();
      State[i_a, t] = (fabs(A1) + A1)*0.5 + machine_precision();
      State[i_b, t] = (fabs(B1) + B1)*0.5 + machine_precision();
    }
    
    return State;
  }
  
  
  //// ODE integration and data assignment to one long vector is wrapped into a function here because
  // - the transformed parameters block does not allow declaring integers (necessary for ragged data structure indexing),
  // - the model block does not alllow for declaring variables within a loop.
  vector unpack(vector state_init, int[] time_max, int[] times, int[] species, int[] pops,
                matrix r, matrix s, matrix g, matrix m_j, // env-dependent
                vector[] h, vector[] c_j, vector[] c_ab, // env-independent
                vector[] u, 
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
                                 exp(r[l, s_species]'), exp(s[l, s_species]'), exp(g[l, s_species]'), exp(m_j[l, s_species]'), // link function: r, l, g, m are still on the log scale
                                 h[l, s_species], c_j[l, s_species], c_ab[l, s_species],  // environmentally-independent parameters have been transformed through lognormal
                                 u[,s_species], // process error
                                 n_s, n_p,
                                 i_j, i_a, i_b);
      
  
      
      // Flattening the matrix into a vector for the location and append it to y_hat local vector yhat[m], and then into function-global vector y_hat[N_yhat].
      // to_vector converts matrix to a column vector in column-major order.

      y_hat[pos_yhat:(pos_yhat - 1 + n_yhat[l])] =
                       to_vector(States[ , segment(times, pos_times, n_t)]); // only select columns with times in the data
      
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
  
  
  int time_init[N_locs];
  int time_max_data[N_locs];
  int times_data[N_times]; // locations/resurveys // old: int times[N_y]; // locations/plots/pops/resurveys

  matrix[N_locs, N_beta] X; // design matrix
  
  int<lower=0> y0[N_y0];
  int<lower=0> y[N_y];
  
}

transformed data {
  // Shift all times to start at 1.
  int time_max[N_locs];
  for (l in 1:N_locs) time_max[l] = time_max_data[l] - time_init[l] + 1; // no vector operations for arrays
  
  int times[N_times];
  for (t in 1:N_times) times[t] = times_data[t] - time_init[rep_locs2times[t]] + 1;
}


parameters {
  //// Level 1, (global, species): species-specific rates, indepndent of environment.
  // … dependent on environment.
  matrix[N_beta, N_totalspecies] Beta_r; // // (J+) flow into the system, dependent on env
  matrix[N_beta, N_totalspecies] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A
  matrix[N_beta, N_totalspecies] Beta_g; // (J-, A+) transition rate from J to A
  matrix[N_beta, N_totalspecies] Beta_m_j; // (J-) density independent mortalitty

  // … independent of environment.
  vector[N_totalspecies] h_log; // (A-, B+), here still unconstrained on log scale, flow out of the system, independent of environment ("intercept only")
  vector[N_totalspecies] c_j_log;
  vector[N_totalspecies] c_ab_log;
  
  vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  vector<lower=0>[2] phi_obs; // lognormal error for observations from predictions
  real<lower=0> shape_par; // normal error for time plot (plots) from locations
  vector[N_totalspecies] u[3]; // normal error for time plot (plots) from locations


  //// Level 2 (locations):
  // "random effect" h rate
  vector<lower=0>[N_totalspecies] h[N_locs];
  vector<lower=0>[N_totalspecies] c_j[N_locs];
  vector<lower=0>[N_totalspecies] c_ab[N_locs];
  
  vector[N_init] state_init_log; // log scale
}


transformed parameters {

  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  matrix[N_locs, N_totalspecies] r_log = X * Beta_r;
  matrix[N_locs, N_totalspecies] s_log = X * Beta_s;
  matrix[N_locs, N_totalspecies] g_log = X * Beta_g;
  matrix[N_locs, N_totalspecies] m_j_log = X * Beta_m_j;

  vector[N_init] state_init = exp(state_init_log);
}

model {

  //---------- PRIORS ---------------------------------
    
  // 1 ./ shape_par ~ normal(0, 1); // for time plot (plots) from locations; Gaussian has the least heavy tails, even less than expontial.
  1 ./ phi_obs ~ normal(0, 0.001); // for observations from predictionsions
    // On prior choice for the overdispersion parameter phi in neg_binomial_2: https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#story-when-the-generic-prior-fails-the-case-of-the-negative-binomial
  
  
  //---------- MODEL ---------------------------------
  // debugging: print("shape_par ", shape_par); print("r_loc ", r_loc[l,]);
  
  //// Level 1 (species) to 2 (locs). Transition rate from A to B, directly, without transformation through env, from level 1 to 2.
  for (l in 1:N_locs) {
    // r[l,]    ~ lognormal(r_log, shape_par); // lognormal acts like exp (log-link) on the logmedian
    h[l,]    ~ gamma(shape_par, shape_par./exp(h_log));
    c_j[l,]  ~ gamma(shape_par, shape_par./exp(c_j_log));
    c_ab[l,] ~ gamma(shape_par, shape_par./exp(c_ab_log));
  }
  
  //// Level 2 (locs) to 3 (plots).
  // Separately fit initial state to feed to integrator.
  y0 ~ neg_binomial_2(state_init[rep_init2y0], phi_obs[obsmethod_y0]); //  y0 ~ neg_binomial_0(state_init[rep_init2y0], theta, phi_obs);
  
  // Process error
  for (stage in 1:3) {
      u[stage,] ~ normal(0, sigma_process[stage]);
  }
  
  // Level 2 (locs).
  vector[N_yhat] y_hat = unpack(state_init, time_max, times, species, pops,
                              r_log, s_log, g_log, m_j_log, // rates matrix[N_locs, N_totalspecies]; will have to be transformed
                              h, c_j, c_ab, // rates vector[N_totalspecies][N_locs]
                              u, // process error
                              n_species, n_pops, n_reobs, n_yhat,
                              N_yhat, N_locs);
  
  
  // Fit predictions to data. (level: location/plot/resurvey/species)
  y ~ neg_binomial_2(y_hat[rep_yhat2y], phi_obs[obsmethod_y]); //   y ~ neg_binomial_0(y_hat[rep_yhat2y], theta, phi_obs);
}
