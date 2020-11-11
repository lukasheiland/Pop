// A generic implementation of the generalized Lotka-Volterra system.
functions {

  vector ds_dt(real time, vector state, vector r, matrix A, int n_pops) {

    vector[n_pops] ds = r + (A * state);
    for (i in 1:n_pops) if (ds[i] < 0 || is_nan(ds[i]) ) ds[i] = 0;
    return ds;
  }
}

data {
  // N
  int<lower=0> N_locs;
  int<lower=0> N_seriesperloc;
  int<lower=0> N_obs; // n solution times, not including t0, for now N_obs have to be the same in every series :(
  int<lower=0> N_species; // for now N_species have to be the same in every series :(
  int<lower=0> N_pops;
  int<lower=0> N_beta;

  // indices for matrix layout
  int<lower=0> I_g[N_species, 3];
  int<lower=0> I_s[N_species*N_species, 3];
  int<lower=0> I_0[N_species*N_species-N_species, 3];
  int<lower=0> I_lg[N_species, 3];
  // int<lower=0> I_l[N_species, 3]; // These are just part of the A_a matrix.
  // indices mapping of life stage matrices to full matrix
  int<lower=0> Map_j[N_species*N_species, 4];
  int<lower=0> Map_a[N_species*N_species, 4];
  // logical: which are adults, which juveniles
  int<lower=1> i_j[N_species];
  int<lower=1> i_a[N_species];

  // actual data
  real time_init[N_locs, N_seriesperloc];
  real times[N_locs, N_seriesperloc, N_obs]; // arrays are row major!
  matrix[N_locs, N_beta] X; // design matrix
  real<lower=0> y_init[N_locs, N_seriesperloc, N_pops];
  vector<lower=0>[N_pops] Y[N_locs, N_seriesperloc, N_obs]; // two-dimensional array of vectors[N_pops]
}


parameters {
  // Global species-specific components of A, dependent on environment, log scale!
  matrix[N_beta, N_species] Beta_s; // (-), here still unconstrained on log scale, shading affectedness of juveniles from A
  matrix[N_beta, N_species] Beta_g; // (+) transition rate from J to A
  
  // Global system border flows, components of f
  matrix[N_beta, N_species] Beta_r; // // (+) flow into the system, dependent on env
  vector[N_species] o; // (-), here still positive, flow out of the system, independent of environment ("intercept only")

  // Global within-size class interaction matrices, components of A
  matrix<upper=0>[N_species, N_species] A_j; // (-) juvenile competition matrix
  matrix<upper=0>[N_species, N_species] A_a; // (-) adult competition matrix

  // Per series inital state
  vector<lower=0>[N_pops] state_init[N_locs, N_seriesperloc];

  real<lower=0> sigma_obs; // lognormal error for observations from predictions
  real<lower=0> sigma_par; // lognormal error for time series (plots) from locations
}


transformed parameters {
  // Location level parameters (unconstrained because on log scale)
  matrix[N_locs, N_species] s_loc;
  matrix[N_locs, N_species] g_loc;
  matrix[N_locs, N_species] r_loc;

  // Series (plot) level parameters:
  vector<lower=0>[N_species] s_series[N_locs, N_seriesperloc]; // (-), here still positive
  vector<lower=0>[N_species] g_series[N_locs, N_seriesperloc];
  vector<lower=0>[N_species] r_series[N_locs, N_seriesperloc];
  vector<lower=0>[N_species] o_series[N_locs, N_seriesperloc]; // (-), here still positive
  
  matrix[N_pops, N_pops] A_series[N_locs, N_seriesperloc]; // combined full matrix at the location level
  vector[N_pops] f_series[N_locs, N_seriesperloc]; // combined vector of r and o at the location level
  // vector<upper=0>[N_pops] l = diagonal(A); // juvenile self-competition
  vector<lower=0>[N_pops] State[N_locs,N_seriesperloc,N_obs];

  //// Level 1 to 2. Environmental effects, on the log scale!
  r_loc = X * Beta_r;
  s_loc = X * Beta_s;
  g_loc = X * Beta_g;
  
  //// Level 2 to 3. Full matrix on series level is broadcasted from species and location level.
  // variables *_series are just matrices/vectors inside arrays, which are accesed e.g. througn matrix[l, z, n, m]
  for (l in 1:N_locs) {
    for (z in 1:N_seriesperloc) {
    
        for (n in 1:(N_species*N_species)) {
          //// Level 3. adult-juvenile competition on loc level
          A_series[l, z, I_s[n,2], I_s[n,3]] = - s_series[l, z, I_s[n,1]]; // (-)!
          //// Level 1 to 3. Within size-class competition is just broadcast from the global parameters
          A_series[l, z, Map_j[n,3], Map_j[n,4]] = A_j[Map_j[n,1], Map_j[n,2]];
          A_series[l, z, Map_a[n,3], Map_a[n,4]] = A_a[Map_a[n,1], Map_a[n,2]];
        }
        
        for (n in 1:(N_species*N_species-N_species)) {
          //// Level 3. Set competition of juveniles on adults to 0.
          A_series[l, z, I_0[n,2], I_0[n,3]] = 0.0;
        }
        
        for (n in 1:N_species) {
          //// Level 3.
          A_series[l, z, I_g[n,2], I_g[n,3]] = g_series[l, z, n];
          A_series[l, z, I_lg[n,2], I_lg[n,3]] = A_series[l, z, I_lg[n,2], I_lg[n,3]] - g_series[l, z, n]; // Juveniles act negatively upon themselves through both competition and transition.
        }
 
        f_series[l, z, i_j] = to_vector(r_series[l, z,]);
        f_series[l, z, i_a] = - to_vector(o_series[l, z,]); // (-)!
 
        // Integrate ode.
        // returns an array of column vectors, 'columns' accessed vie State[i]
        State[l, z,] = ode_rk45_tol(ds_dt,
                                state_init[l, z], time_init[l, z], times[l, z,],
                                1e-5, 1e-3, 10000, // integrator tolerance settings: real rel_tol, real abs_tol, int max_num_steps
                                f_series[l, z], A_series[l, z], // model parameters
                                N_pops);
    } // loop: z in N_series
  } // loop: l in N_locs
}


model {
  
  //---------- PRIORS ---------------------------------
  
  //// Level 1. log-scale
  o ~ normal(-2, 0.02);
  
  // sigma ~ cauchy(0, 0.000001);
  sigma_obs ~ exponential(50000); // for observations from predictionsions
  sigma_par ~ exponential(100000); // for time series (plots) from locations


  //---------- MODEL ---------------------------------
  
  //// Level 2 to 3. Sample series level (plot level) parameters from location parameters.
  for (l in 1:N_locs) {
    for (z in 1:N_seriesperloc) {
      r_series[l, z] ~ lognormal(r_loc[l], sigma_par);
      s_series[l, z] ~ lognormal(s_loc[l], sigma_par); // (-), here still positive
      g_series[l, z] ~ lognormal(g_loc[l], sigma_par);
      //// Directly from level 1 to 3.
      o_series[l, z] ~ lognormal(o, sigma_obs); // (-), here still positive
      
      //// Level 3. Fit time series level predictions to data.
      y_init[l, z, ] ~ lognormal(log(state_init[l, z]), sigma_obs); // Separately fitting initial state to feed to integrator.

      for (t in 1:N_obs) {
        Y[l,z,t] ~ lognormal(log(State[l, z, t]), sigma_obs); // Lognormal model of states.
      } // loop: t in N_obs
      
    } // loop: z in N_seriesperloc
  } // loop: l in N_locs
}

