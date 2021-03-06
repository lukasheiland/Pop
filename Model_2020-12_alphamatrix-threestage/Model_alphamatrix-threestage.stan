functions {
  
  //// Alpha matrix (local) from global Alpha matrix and species vector
  // matrix constructAlpha(int species, matrix Alpha) {
  //   return Alpha[species, species];
  // }
  
  
  //// ODE system
  vector ds_dt(real time, vector state, vector r, vector s, vector g, vector h, matrix Alpha_j, matrix Alpha_ab, int N_species) {
    
    // Structure of state[N_pops]: species/stage
    
    vector[N_species] J = segment(state, 1, N_species);
    vector[N_species] A = segment(state, N_species+1, N_species);
    vector[N_species] B = segment(state, N_species+N_species+1, N_species);
    vector[N_species] AB = A + B;

    
    vector[N_species] dJ = r - (Alpha_j * J  + s*sum(AB)).*J - g .* J;
    vector[N_species] dA = g .* J - (Alpha_ab * AB).*A - h .* A;
    vector[N_species] dB = h .* A - (Alpha_ab * AB).*B;
    
    return append_row(dJ, append_row(dA, dB));
  }
  
  
  //// ODE integration and data assignment to one long vector is wrapped into a function here because
  // - the transformed parameters block does not allow declaring integers (necessary for ragged data structure indexing),
  // - the model block does not alllow for declaring variables within a loop.
  vector simulate(vector state_init, real[] time_init, real[] times, int[] species,
                  matrix r, matrix s, matrix g, vector[] h, matrix Alpha_j, matrix Alpha_ab,
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
      
      // Slicing the global matrix down to local level present species.
      int s_species[n_s] = segment(species, pos_species, n_s);
      matrix[n_s, n_s] Alpha_j_loc = Alpha_j[s_species, s_species];
      matrix[n_s, n_s] Alpha_ab_loc = Alpha_ab[s_species, s_species];

      // Returns an array of vectors, accessed vie State[i]
      // ode_rk45_tol(function ode, vector initial_state, real initial_time, real[] times, real rel_tol, real abs_tol, int max_num_steps, ...)
      // print("r", exp(r[l, s_species])'); print("h", h[l, s_species]);
      vector[n_p] States[n_t] =
        ode_rk45_tol(ds_dt,
                     segment(state_init, pos_pops, n_p), // segment has length n_pops. Structure state_init: locs/stage/species 
                     time_init[l],
                     segment(times, pos_times, n_t),
                     1e-5, 1e-3, 10000, // integrator tolerance settings: real rel_tol, real abs_tol, int max_num_steps
                     exp(r[l, s_species]'), exp(s[l, s_species]'), exp(g[l, s_species]'), h[l, s_species], // link function: r, l, g are still on the log scale
                     Alpha_j_loc, Alpha_ab_loc,
                     n_s);
      
      // Flattening array into local vector s_yhat[m], and then into function-global vector y_hat[N_yhat].
      vector[n_p] s_yhat = States[1];
      for (t in 2:n_t) s_yhat = append_row(s_yhat, States[t]);

      y_hat[pos_yhat:(pos_yhat - 1 + n_yhat[l])] = s_yhat; // Location wise casting to array. Structure within s_hat: resurveys/stages/species
      
      //// Note for vectorized version of negative replacements. Highly // y_hat = (fabs(y_hat)+y_hat)/2;

      pos_species = pos_species + n_s;
      pos_pops = pos_species + n_p;
      pos_times = pos_times + n_t;
      pos_yhat = pos_yhat + n_yhat[l];
    }

  return y_hat; // Structure: locations/resurveys/pops(==stages/species)
  }

  
  //// Implementation of the negative binomial probability mass — including definition of support 0 (R style)
  // real neg_binomial_0_lpmf(int[] y, vector mu, real theta, real phi) {
  //   real t;
  //   for (i in 1:size(mu)) {
  //     if (y[i] == 0) {
  //       t += log_sum_exp(bernoulli_lpmf(1 | theta),
  //                        bernoulli_lpmf(0 | theta) + poisson_lpmf(y[i] | mu[i]));
  //     }
  //     else {
  //       t += bernoulli_lpmf(0 | theta) + poisson_lpmf(y[i] | mu[i]);
  //     }
  //   }
  //   return t;
  // }
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
  real times[N_y]; // locations/plots/pops/resurveys
  
  matrix[N_locs, N_beta] X; // design matrix
  
  int<lower=0> y0[N_y0];
  int<lower=0> y[N_y];
  
}


parameters {
  /// Global species-specific rates, indepndent of environment.
  // … depndent on environment.
  matrix[N_beta, N_totalspecies] Beta_r; // // (J+) flow into the system, dependent on env
  matrix[N_beta, N_totalspecies] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A
  matrix[N_beta, N_totalspecies] Beta_g; // (J-, A+) transition rate from J to A
  // … indepndent of environment.
  vector[N_totalspecies] h_log; // (A-, B+), here still unconstrained on log scale, flow out of the system, independent of environment ("intercept only")
  
  // Per plot inital state
  vector<lower=0>[N_init] state_init;
  
  //// Version with random rates. Per time series (plot) level parameters:
  // vector<lower=0>[N_totalspecies] r[N_locs];
  // vector<lower=0>[N_totalspecies] s[N_locs];
  // vector<lower=0>[N_totalspecies] g[N_locs];
  
  vector<lower=0>[N_totalspecies] h[N_locs];

  // Global within-size class interaction matrices, components of A
  matrix<lower=0>[N_totalspecies, N_totalspecies] Alpha_j; // (-) juvenile competition matrix
  matrix<lower=0>[N_totalspecies, N_totalspecies] Alpha_ab; // (-) adult competition matrix
  
  // real<lower=0, upper=1> theta; // zero inflation prob
  
  real<lower=0> phi_obs; // lognormal error for observations from predictions
  real<lower=0> sigma_par; // normal error for time plot (plots) from locations
}


transformed parameters {

  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  matrix[N_locs, N_totalspecies] r_log = X * Beta_r;
  matrix[N_locs, N_totalspecies] s_log = X * Beta_s;
  matrix[N_locs, N_totalspecies] g_log = X * Beta_g;
  
  //// Level 2 (locs).
  vector[N_yhat] y_hat = simulate(state_init, time_init, times, species,
                                  r_log, s_log, g_log, h, // model rates vector[N_totalspecies][N_locs]
                                  Alpha_j, Alpha_ab,
                                  n_species, n_pops, n_reobs, n_yhat,
                                  N_yhat, N_locs);
  
  
  //// Version with random rates. Level 2 (locs).
  // vector[N_yhat] y_hat = simulate(state_init, time_init, times, species,
  //                                 r, s, g, h, // model rates vector[N_totalspecies][N_locs]
  //                                 Alpha_j, Alpha_ab,
  //                                 n_species, n_pops, n_reobs, n_yhat,
  //                                 N_yhat, N_locs);
}

model {

  //---------- PRIORS ---------------------------------
    
  sigma_par ~ normal(0, 0.01); // for time plot (plots) from locations; Gaussian has the least heavy tails, even less than expontial.
  1/phi_obs ~ normal(0, 1); // for observations from predictionsions
    // On prior choice for the overdispersion in negative binomial 2: https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#story-when-the-generic-prior-fails-the-case-of-the-negative-binomial
  
  
  //---------- MODEL ---------------------------------
  //// debugging printout for checking initial value acceptance
  // print("h ", h); print("h_log", h_log);  print("state_init ", state_init); print("y_hat ", y_hat);
  
  for (l in 1:N_locs) { // debugging: print("sigma_par ", sigma_par); print("r_loc ", r_loc[l,]);
    
    //// Version with random rates:
    //// Level 2 (locs). Population rates dependent on environment.
    // r[l,] ~ lognormal(r_log[l,], sigma_par);
    // s[l,] ~ lognormal(s_log[l,], sigma_par);
    // g[l,] ~ lognormal(g_log[l,], sigma_par);

    //// Level 1 (species) to 2 (locs). Transition rate from A to B, directly, without transformation through env, from level 1 to 2.
    h[l,] ~ lognormal(h_log, sigma_par); // lognormal acts like exp on the logmedian
  }
  
  //// Level 2 (locs). Separately fit initial state to feed to integrator.
  y0 ~ neg_binomial_2(state_init[rep_init2y0], phi_obs); //  y0 ~ neg_binomial_0(state_init[rep_init2y0], theta, phi_obs);
  
  //// Level 2 (locs) to 3 (plots). Fit predictions on level location/plot/resurvey/species to data.
  y ~ neg_binomial_2(y_hat[rep_yhat2y], phi_obs); //   y ~ neg_binomial_0(y_hat[rep_yhat2y], theta, phi_obs);
}
