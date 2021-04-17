functions {
  
  //// Difference equations
  vector[] simulate(vector initialstate, int time_max, int[] times,
                  vector g, vector m_j, vector r, vector s,
                  vector b, vector c_j, vector c_a, vector c_b, vector h, vector m_a,
                  real ba_a_avg, real ba_a_upper,
                  int N_spec, int N_pops,
                  matrix u,
                  int[] i_j, int[] i_a, int[] i_b) {
    
    
    // State matrix with [species, times].
    vector[N_pops] State[time_max];
    State[1, ] = initialstate;

    for (t in 2:time_max) {
      // Structure of state[N_pops]: stage/species
      // State has times rows (array row-major),  u has times cols (matrix: col-major access)
      
      vector[N_spec] J_log = State[t-1, i_j];
      vector[N_spec] A_log = State[t-1, i_a];

      vector[N_spec] J = exp(J_log);
      vector[N_spec] A = exp(A_log);
      vector[N_spec] B = exp(State[t-1, i_b]);
      
      real BA = sum(A*ba_a_avg + B);
      
      /// Ricker model: assign to the the log states
      // Note: log1p(expm1(a) + exp(b)) == log(exp(a) + exp(b))
      State[t, i_j]  =  log1p(expm1(J_log) + r) - c_j*sum(J) - s*BA - m_j - g + u[i_j, t-1]; // equivalent to … =  J + (r - (c_j*sum(J) + s*AB + g + m_j).*J * dt;
      State[t, i_a]  =  log1p(A + expm1(J_log - g)) - c_a*BA - m_a - h + + u[i_a, t-1];
      State[t, i_b]  =  log1p(B + expm1(A_log - h) * ba_a_upper) + b - c_b*sum(B) + u[i_b, t-1]; // exp(A_log - h) * ba_a_upper is the input of new basal area through ingrowth of count A*exp(-h) from A
      
    }
    
    // print("State: ", State);
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
  int<lower=1> rep_obsmethod2pops[N_pops];

  int<lower=1> i_j[N_species]; // e.g. 1:4
  int<lower=1> i_a[N_species]; // e.g. 5:9, etc.
  int<lower=1> i_b[N_species];

  real dbh_lower_a;
  real dbh_lower_b;
  
  int times[N_locs, N_times]; // assumes start at 1. (Even though the number of times is the same, the times can deviate)
  int time_max[N_locs];
  int timespan_max; // max(time_max) - time_globalmin

  matrix[N_locs, N_beta] X; // design matrix
  
  // The response.
  vector[N_pops] y_log [N_locs, N_plots, N_times];
}


transformed data {
  // Shift all times to start at 1.
  real ba_a_upper = pi() * (dbh_lower_b/2)^2 * 1e-6; // # pi*r^2, mm^2 to m^2
  real ba_a_avg = pi() * ((dbh_lower_a + dbh_lower_b)/2/2)^2 * 1e-6;
  
  //// Data for separate fitting of the initial state
  // vector[N_pops] y0_log [N_locs, N_plots] = y_log[ , , 1, ];
}


parameters {
  //// Level 1, (global, species): species-specific rates, indepndent of environment.
  // … dependent on environment.
  matrix[N_beta, N_species] Beta_g; // (J-, A+) transition rate from J to A
  matrix[N_beta, N_species] Beta_m_j; // (J-) density independent mortalitty
  matrix[N_beta, N_species] Beta_r; // // (J+) flow into the system, dependent on env
  matrix[N_beta, N_species] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A

  // … independent of environment
  vector<lower=0>[N_species] b;
  vector<lower=0>[N_species] c_j;
  vector<lower=0>[N_species] c_a;
  vector<lower=0>[N_species] c_b;
  vector<lower=0>[N_species] h; // (A-, B+), here still unconstrained on log scale, flow out of the system, independent of environment ("intercept only")
  vector<lower=0>[N_species] m_a;
  
  // vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  vector<lower=0>[2] sigma_obs; // lognormal error for observations from predictions

  matrix[N_pops, timespan_max] u[N_locs];

  vector[N_pops] state_init_log[N_locs, N_plots];
}


transformed parameters {
  
  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  matrix[N_locs, N_species] g_log = X * Beta_g;
  matrix[N_locs, N_species] m_j_log = X * Beta_m_j;
  matrix[N_locs, N_species] r_log = X * Beta_r;
  matrix[N_locs, N_species] s_log = X * Beta_s;
  
}


model {
  
  vector[N_pops] y_hat_log[N_locs, N_plots, N_times];
  
  //---------- PRIORS ---------------------------------
  
  // Beta_r[1,] ~ normal(2, 2); // intercept
  // 1 ./ shape_par ~ normal(0, 100); 
  // sigma_process ~ normal(0, 0.01);
  sigma_obs ~ normal(0, [2, 0.2]); // for observations from predictions
  Beta_g[1,] ~ normal(-12, 1);
  to_vector(Beta_g[2:N_beta,]) ~ std_normal();
  
  h ~ normal(0.5, 0.2);
  
  //---------- MODEL ---------------------------------
  
  for(l in 1:N_locs) {
    
    target += std_normal_lpdf(to_vector(u[l])); // to_vector(u) ~ std_normal();
    
    for (p in 1:N_plots) {
      
      //// Debugging alternatives
      // print(y_hat_log[l, p, ]);
      // y0_log[l, p, ] ~ normal(state_init_log[l, p, ], sigma_obs[rep_obsmethod2pops]);
        
      y_hat_log[l, p, ] = simulate(state_init_log[l, p, ], time_max[l], times[l, ],
                                     exp(g_log[l, ]'), exp(m_j_log[l, ]'), exp(r_log[l, ]'), exp(s_log[l, ]'),
                                     b, c_j, c_a, c_b, h, m_a,
                                     ba_a_avg, ba_a_upper,
                                     N_species, N_pops,
                                     u[l],
                                     i_j, i_a, i_b);
  
      // for (t in 2:N_times) { // in case of separate y0_log fitting
      for (t in 1:N_times) {
        
        y_log[l, p, t] ~ normal(y_hat_log[l, p, t], sigma_obs[rep_obsmethod2pops]);
      
      }
    }
  }
}
