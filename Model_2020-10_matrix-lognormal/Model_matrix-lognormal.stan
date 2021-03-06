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
  int<lower=0> N_series;
  int<lower=0> N_obs; // not including t0, for now N_obs have to be the same in every series :(
  int<lower=0> N_species; // for now N_species have to be the same in every series :(
  int<lower=0> N_pops;

  // indices for matrix layout
  int<lower=0> I_g[N_species, 3];
  int<lower=0> I_s[N_species*N_species, 3];
  int<lower=0> I_0[N_species*N_species-N_species, 3];
  int<lower=0> I_lg[N_species, 3];
  // int<lower=0> I_l[N_species, 3]; // These are just part of the A_a matrix.
  int<lower=0> Map_j[N_species*N_species, 4];
  int<lower=0> Map_a[N_species*N_species, 4];

  int<lower=1> i_j[N_species];
  int<lower=1> i_a[N_species];

  real time_init[N_series];
  real times[N_series, N_obs]; // arrays are row major!
  real<lower=0> y_init[N_series, N_pops];
  vector<lower=0>[N_pops] Y[N_series, N_obs]; // two-dimensional array of vectors[N_pops]
}

parameters {
  // components of A, i.e. population interactions
  vector<lower=0>[N_species] g; // (+) transition rate from J to A
  vector<upper=0>[N_species] s; // (-) shading affectedness of juveniles from A

  matrix<upper=0>[N_species, N_species] A_j; // (-) juvenile competition matrix
  matrix<upper=0>[N_species, N_species] A_a; // (-) adult competition matrix

  // system border flows
  vector<lower=0>[N_species] r; // (+) flow into the system
  vector<upper=0>[N_species] o; // (-) flow out of the system

  // inital state
  vector<lower=0>[N_pops] state_init[N_series];

  // lognormal error at observation state
  real<lower=0> sigma; // vector<lower=0>[N_species] sigma;
}

transformed parameters {

  matrix[N_pops, N_pops] A;
  vector[N_pops] f;
  // vector<upper=0>[N_pops] l = diagonal(A); // juvenile self-competition
  vector<lower=0>[N_pops] State[N_series,N_obs];

  // // broadcasting of parameters into the interaction matrix A
  for (n in 1:(N_species*N_species)) {
    A[I_s[n,2], I_s[n,3]] = s[I_s[n,1]];
    A[Map_j[n,3], Map_j[n,4]] = A_j[Map_j[n,1], Map_j[n,2]];
    A[Map_a[n,3], Map_a[n,4]] = A_a[Map_a[n,1], Map_a[n,2]];
  }
  for (n in 1:(N_species*N_species-N_species)) {
    A[I_0[n,2], I_0[n,3]] = 0.0;
  }
  for (n in 1:N_species) {
    A[I_g[n,2], I_g[n,3]] = g[n];
    A[I_lg[n,2], I_lg[n,3]] = A[I_lg[n,2], I_lg[n,3]] - g[n]; // Juveniles act negatively upon themselves through competition and transition.
  }

  f[i_j] = r;
  f[i_a] = o;

  // Integrate ode.
  // returns an array of column vectors, 'columns' accessed vie State[i]
  for(z in 1:N_series) {
      State[z,] = ode_rk45_tol(ds_dt,
                              state_init[z], time_init[z], times[z,],
                              1e-5, 1e-3, 10000, // old integrator tolerances settings: real rel_tol, real abs_tol, int max_num_steps
                              f, A, // model parameters
                              N_pops);
  }

}


model {
  // Priors
  // for (s in 1:N_pops) nA[, s] ~ exponential(10); // negative decaying prior, because all-competitive is assumed, matrices col-major!
  // sigma ~ cauchy(0, 0.000001); // sigma is constrained to be > 0, thus half-cauchy
  sigma ~ exponential(50000); // sigma is constrained to be > 0, thus half-cauchy


  // Model
  for(z in 1:N_series) {
    y_init[z, ] ~ lognormal(log(state_init[z]), sigma); // Separately fitting initial state to feed to integrator.

    for (t in 1:N_obs) {
      Y[z,t] ~ lognormal(log(State[z,t]), sigma); // Lognormal model of states.
    }

  }
}

