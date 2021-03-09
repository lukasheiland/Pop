functions {
  
  //// ODE system 
  vector ds_dt(real time, vector state, vector r, vector s, vector g, vector m_j,
               vector h, vector c_j, vector c_ab, int N_spec) {
    
    // Structure of state[N_pops]: stage/species
    
    vector[N_spec] J = segment(state, 1, N_spec);
    vector[N_spec] A = segment(state, N_spec+1, N_spec);
    vector[N_spec] B = segment(state, N_spec+N_spec+1, N_spec);
    real AB = sum(A + B);
    
    vector[N_spec] dJ = r - (c_j*sum(J) + s*AB + g + m_j).*J;
    vector[N_spec] dA = g .* J - (c_ab*AB + h).*A;
    vector[N_spec] dB = h .* A - (c_ab*AB).*B;
    
    return append_row(dJ, append_row(dA, dB));
  }
  
  
  //// ODE integration and data assignment to one long vector is wrapped into a function here because
  // - the transformed parameters block does not allow declaring integers (necessary for ragged data structure indexing),
  // - the model block does not alllow for declaring variables within a loop.
  vector simulate(vector state_init, real[] time_init, real[] times, int[] species,
                  matrix r, matrix s, matrix g, matrix m_j, // env-dependent
                  vector[] h, vector[] c_j, vector[] c_ab, // env-independent
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

      // Returns an array of vectors, accessed vie State[i]
      // ode_rk45_tol(function ode, vector initial_state, real initial_time, real[] times, real rel_tol, real abs_tol, int max_num_steps, ...)
      // print("r", exp(r[l, s_species])'); print("h", h[l, s_species]);
      vector[n_p] States[n_t] =
        ode_rk45_tol(ds_dt,
                     segment(state_init, pos_pops, n_p), // segment has length n_pops. Structure state_init: locs/stage/species 
                     time_init[l],
                     segment(times, pos_times, n_t),
                     1e-5, 1e-3, 10000, // integrator tolerance settings: real rel_tol, real abs_tol, int max_num_steps
                     exp(r[l, s_species]'), exp(s[l, s_species]'), exp(g[l, s_species]'), exp(m_j[l, s_species]'), // link function: r, l, g, m are still on the log scale
                     h[l, s_species], c_j[l, s_species], c_ab[l, s_species],  // environmentally-independent parameters have been transformed through lognormal
                     n_s);
      
      // Flattening array into local vector s_yhat[m], and then into function-global vector y_hat[N_yhat].
      vector[n_p] s_yhat = States[1];
      for (t in 2:n_t) s_yhat = append_row(s_yhat, States[t]);

      y_hat[pos_yhat:(pos_yhat - 1 + n_yhat[l])] = s_yhat; // Location wise casting to array. Structure within s_hat: resurveys/stages/species
      
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
  int<lower=0> rep_init2y0[N_y0]; // repeat predictions on level "locations/pops" n_plots times to "locations/pops/plots"
  int<lower=0> rep_yhat2y[N_y]; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  
  
  //// actual data
  int species[N_species]; // locations/species
  real time_init[N_locs];
  real times[N_times]; // locations/plots/pops/resurveys
  
  matrix[N_locs, N_beta] X; // design matrix
  
  int<lower=0> y0[N_y0];
  int<lower=0> y[N_y];
  
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
  
  real<lower=0> phi_obs; // negbinomial error for observations from predictions
  real<lower=0> shape_par; // gamma error for time plot (plots) from locations
  
  
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
  
  //// Level 2 (locs).
  vector[N_yhat] y_hat = simulate(state_init, time_init, times, species,
                                  r_log, s_log, g_log, m_j_log, // rates matrix[N_locs, N_totalspecies]; will have to be transformed
                                  h, c_j, c_ab, // rates vector[N_totalspecies][N_locs]
                                  n_species, n_pops, n_reobs, n_yhat,
                                  N_yhat, N_locs);
                                     
}

model {

  //---------- PRIORS ---------------------------------
    
  1/shape_par ~ normal(0, 1); // for time plot (plots) from locations; Gaussian has the least heavy tails, even less than expontial.
  1/phi_obs ~ normal(0, 0.001);  // for observations from predictionsions
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
  y0 ~ neg_binomial_2(state_init[rep_init2y0], phi_obs); //  y0 ~ neg_binomial_0(state_init[rep_init2y0], theta, phi_obs);
  
  // Fit predictions to data. (level: location/plot/resurvey/species)
  y ~ neg_binomial_2(y_hat[rep_yhat2y], phi_obs); //   y ~ neg_binomial_0(y_hat[rep_yhat2y], theta, phi_obs);
}
