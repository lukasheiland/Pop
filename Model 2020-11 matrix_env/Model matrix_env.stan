// A generic implementation of the generalized Lotka-Volterra system.
functions {

  vector ds_dt(real time, vector state, vector r, matrix A, int n_pops) {

    vector[n_pops] ds = r + (A * state);
    // for(i in 1:n_pops) if(ds[i] < 0) ds[i] = 0;
    return ds;
  }
}

data {
  // N
  int<lower=0> N_locs;
  int<lower=0> N_seriesperloc;
  int<lower=0> N_obs; // not including t0, for now N_obs have to be the same in every series :(
  int<lower=0> N_species; // for now N_species have to be the same in every series :(
  int<lower=0> N_pops;

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
  real<lower=0> y_init[N_locs, N_seriesperloc, N_pops];
  vector<lower=0>[N_pops] Y[N_locs, N_seriesperloc, N_obs]; // two-dimensional array of vectors[N_pops]
}

parameters {
  // Global species specific parameters
  // components of A, i.e. population interactions
  vector<lower=0>[N_species] g; // (+) transition rate from J to A
  vector<upper=0>[N_species] s; // (-) shading affectedness of juveniles from A

  matrix<upper=0>[N_species, N_species] A_j; // (-) juvenile competition matrix
  matrix<upper=0>[N_species, N_species] A_a; // (-) adult competition matrix

  // system border flows
  vector<lower=0>[N_species] r; // (+) flow into the system
  vector<upper=0>[N_species] o; // (-) flow out of the system

  // inital state
  vector<lower=0>[N_pops] state_init[N_locs, N_seriesperloc];

  // lognormal error at observation state
  real<lower=0> sigma; // vector<lower=0>[N_species] sigma;
  // lognotmal error for location parameters
  real<lower=0> theta;
}

transformed parameters {
  // Location level parameters
  vector[N_species] g_loc[N_locs];
  vector[N_species] s_loc[N_locs];
  vector[N_species] r_loc[N_locs];
  vector[N_species] o_loc[N_locs];

  matrix[N_pops, N_pops] A_loc[N_locs]; // combined full matrix at the location level
  vector[N_pops] f_loc[N_locs]; // combined vector of r and o at the location level
  // vector<upper=0>[N_pops] l = diagonal(A); // juvenile self-competition
  vector<lower=0>[N_pops] State[N_locs,N_seriesperloc,N_obs];

  for(l in 1:N_locs) {
    // variables *_loc are just matrices/vectors inside arrays, which are accesed e.g. througn matrix[l, n, m]

    // broadcasting of parameters into the interaction matrix A
    for (n in 1:(N_species*N_species)) {
      A_loc[l,  I_s[n,2], I_s[n,3]] = s_loc[l,  I_s[n,1]]; // adult-juvenile competition on loc level

      // Within size-class competition is just broadcast from the global parameters
      A_loc[l,  Map_j[n,3], Map_j[n,4]] = A_j[Map_j[n,1], Map_j[n,2]];
      A_loc[l,  Map_a[n,3], Map_a[n,4]] = A_a[Map_a[n,1], Map_a[n,2]];
    }
    for (n in 1:(N_species*N_species-N_species)) {
      A_loc[l,  I_0[n,2], I_0[n,3]] = 0.0;
    }
    for (n in 1:N_species) {
      A_loc[l,  I_g[n,2], I_g[n,3]] = g_loc[l,  n];
      A_loc[l,  I_lg[n,2], I_lg[n,3]] = A_loc[l,  I_lg[n,2], I_lg[n,3]] - g_loc[l,  n]; // Juveniles act negatively upon themselves through both competition and transition.
    }

    f_loc[l, i_j] = r_loc[l];
    f_loc[l, i_a] = o_loc[l];

    // Integrate ode.
    // returns an array of column vectors, 'columns' accessed vie State[i]
    for(z in 1:N_seriesperloc) {
        State[l, z,] = ode_rk45_tol(ds_dt,
                                state_init[l, z], time_init[l, z], times[l, z,],
                                1e-5, 1e-3, 10000, // old integrator tolerances settings: real rel_tol, real abs_tol, int max_num_steps
                                f_loc[l], A_loc[l], // model parameters
                                N_pops);
    } // loop: z in N_series
  } // loop: l in N_locs
}


model {
  // Priors
  // for (s in 1:N_pops) nA[, s] ~ exponential(10); // negative decaying prior, because all-competitive is assumed, matrices col-major!
  // sigma ~ cauchy(0, 0.000001); // sigma is constrained to be > 0, thus half-cauchy
  sigma ~ exponential(50000); // sigma is constrained to be > 0, thus half-cauchy
  theta ~ exponential(10000);


  // Model
  for (l in 1:N_locs) {
    // sample location level parameters
    g_loc[l] ~ lognormal(log(g), theta);
    s_loc[l] ~ lognormal(log(s), theta);
    r_loc[l] ~ lognormal(log(r), theta);
    o_loc[l] ~ lognormal(log(o), theta);

    for (z in 1:N_seriesperloc) {
      y_init[l, z, ] ~ lognormal(log(state_init[l, z]), sigma); // Separately fitting initial state to feed to integrator.

      for (t in 1:N_obs) {
        Y[l,z,t] ~ lognormal(log(State[l, z, t]), sigma); // Lognormal model of states.
      }
    } // loop: z in N_seriesperloc
  } // loop: l in N_locs
}

