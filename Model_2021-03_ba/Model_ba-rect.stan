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

      State[t, i_j]  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      State[t, i_a]  =  (g .* J + (A - h .*A )) ./ (1 + c_b*BA_sum);
      State[t, i_b]  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
      
    }
    
    // print("State 1 after: ", State[1,]);

    return State[times,];
  }
  
    //// Gets the specific Jacobian, given states and parameters
  matrix jacobian(real J1, real J2,   real A1, real A2,   real B1, real B2,
                real b1, real b2,   real c_b1, real c_b2,   real c_j1, real c_j2,   real g1, real g2,   real h1, real h2,   real l1, real l2,     real r1, real r2,   real s1, real s2,
                real p, real q) {
      
      matrix[6, 6] J;
      
      J = [
              [- (g1 - 1)/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1) - (c_j1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                             -(c_j1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2, (p*r1)/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1) - (p*s1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                         -(p*s1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2, r1/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1) - (s1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                     -(s1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2],
              [                                                            -(c_j2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2, - (g2 - 1)/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1) - (c_j2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                         -(p*s2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2, (p*r2)/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1) - (p*s2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                     -(s2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2, r2/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1) - (s2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2],
              [                                                                                                                    g1/(c_b1*(B1 + B2 + A1*p + A2*p) + 1),                                                                                                                                                       0,                                          - (h1 - 1)/(c_b1*(B1 + B2 + A1*p + A2*p) + 1) - (p*c_b1*(A1 - A1*h1 + J1*g1))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                       -(p*c_b1*(A1 - A1*h1 + J1*g1))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                   -(c_b1*(A1 - A1*h1 + J1*g1))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                   -(c_b1*(A1 - A1*h1 + J1*g1))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2],
              [                                                                                                                                                      0,                                                                                                                     g2/(c_b2*(B1 + B2 + A1*p + A2*p) + 1),                                                                                       -(p*c_b2*(A2 - A2*h2 + J2*g2))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                          - (h2 - 1)/(c_b2*(B1 + B2 + A1*p + A2*p) + 1) - (p*c_b2*(A2 - A2*h2 + J2*g2))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                   -(c_b2*(A2 - A2*h2 + J2*g2))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                   -(c_b2*(A2 - A2*h2 + J2*g2))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2],
              [                                                                                                                                                      0,                                                                                                                                                       0,                     (h1*q + b1*h1*q)/(c_b1*(B1 + B2 + A1*p + A2*p) + 1) - (p*c_b1*(B1 + B1*b1 + A1*h1*q + A1*b1*h1*q))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                        -(p*c_b1*(B1 + B1*b1 + A1*h1*q + A1*b1*h1*q))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                         (b1 + 1)/(c_b1*(B1 + B2 + A1*p + A2*p) + 1) - (c_b1*(B1 + B1*b1 + A1*h1*q + A1*b1*h1*q))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                    -(c_b1*(B1 + B1*b1 + A1*h1*q + A1*b1*h1*q))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2],
              [                                                                                                                                                      0,                                                                                                                                                       0,                                                                        -(p*c_b2*(B2 + B2*b2 + A2*h2*q + A2*b2*h2*q))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                     (h2*q + b2*h2*q)/(c_b2*(B1 + B2 + A1*p + A2*p) + 1) - (p*c_b2*(B2 + B2*b2 + A2*h2*q + A2*b2*h2*q))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                    -(c_b2*(B2 + B2*b2 + A2*h2*q + A2*b2*h2*q))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                         (b2 + 1)/(c_b2*(B1 + B2 + A1*p + A2*p) + 1) - (c_b2*(B2 + B2*b2 + A2*h2*q + A2*b2*h2*q))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2]
          ];
      
      return J;
  }
  
  
  //// Returns the two norm of a matrix.
  real norm(matrix M) {
    return max(singular_values(M)); // "the maximum amplification is given by the maximum singular value" hence this is the 2-norm of a matrix
    // return max(fabs(eigenvalues_sym(M))); # not symmetric
  }
  
  
  //// Difference equations simulated up to the fix point given a maximum tolerance over all states.
  // Expects a state vector[N_pops]
  // returns a state vector of the form [J1, …, A1, …, B1, …, BA1, …, iterations]
  vector iterateFix(vector state_0,
                     vector g, vector r, vector s, vector l,
                     vector b, vector c_j, vector c_b, vector h,
                     real ba_a_avg, real ba_a_upper,
                     int N_spec, int N_pops,
                     int[] i_j, int[] i_a, int[] i_b,
                     real tolerance_fix) {
                       
    vector[N_pops+N_spec] s_0 = append_row(state_0, [0, 0]'); // two additional states for BA
    
    // initialize while loop variables
    vector[N_pops+N_spec] s_1 = s_0 - tolerance_fix - 1; 
    int i = 0;
    int notconvergent = 1;
    vector[N_spec] eps_ba = [1.0, 1.0]'; // tolerance_fix is set to <upper=0.5>
    
    while ( (notconvergent && max(eps_ba) > tolerance_fix) || i < 5000 ) { // 
      
      s_0 = s_1;
      
      vector[N_spec] J = s_0[i_j];
      vector[N_spec] A = s_0[i_a];
      vector[N_spec] B = s_0[i_b];

      vector[N_spec] BA = A*ba_a_avg + B;
      vector[N_spec] BA_1;
      real BA_sum = sum(BA);
      

      s_1[i_j]  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      s_1[i_a]  =  (g .* J + (A - h .*A )) ./ (1 + c_b*BA_sum);
      s_1[i_b]  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
      
      BA_1 =  s_1[i_a]*ba_a_avg + s_1[i_b]; // New BA as additional state.
      s_1[(N_pops+1):] = BA_1;
      
      notconvergent = (1 <= norm(jacobian(s_1[i_j[1]], s_1[i_j[2]], s_1[i_a[1]], s_1[i_a[2]], s_1[i_b[1]], s_1[i_b[2]], b[1], b[2], c_b[1], c_b[2], c_j[1], c_j[2], g[1], g[2], h[1], h[2], l[1], l[2], r[1], r[2], s[1], s[2], ba_a_avg, ba_a_upper)) );
      eps_ba = fabs((BA_1 - BA)./BA_1);
      
      i += 1;
    }
    
    return append_row(s_1, append_row(eps_ba, i)); // int i gets cast to real
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


  
  int times[N_locs, N_times]; // assumes start at 1. (Even though the number of times is the same, the times can deviate)
  int time_max[N_locs];
  int timespan_max; // max(time_max) - time_globalmin

  matrix[N_locs, N_beta] X; // design matrix
  matrix[N_locs, N_species] L_loc; // design matrix

  
  // The response.
  vector[N_pops] y0_loc [N_locs];
  vector[N_pops] y [N_locs, N_times, N_plots];
  
  
  // Settings
  real<upper=0.5> tolerance_fix;
  real dbh_lower_a; // 100
  real dbh_lower_b; // 200
  
}


transformed data {
  
  int N_genstates = N_pops + N_species;
  
  // Shift all times to start at 1.
  real ba_a_upper = pi() * (dbh_lower_b/2)^2 * 1e-6; // / pi*r^2, mm^2 to m^2
  real ba_a_avg = pi() * ((dbh_lower_a + dbh_lower_b)/2/2)^2 * 1e-6;
  
  vector[N_pops] y0_loc_log [N_locs] = log(y0_loc);
  
  // priors
  vector[N_species] mu_b_log = [-1.5, -1.5]';
  // vector[N_species] mu_c_a_log = [, ]';
  vector[N_species] mu_c_b_log = [-2, -2]';
  vector[N_species] mu_c_j_log = [-1.5, -1.5]';
  vector[N_species] mu_g_logit = [0.7, -0.2]';
  vector[N_species] mu_h_logit = [-0.5, -0.01]';
  // vector[N_species] mu_l_log = [, ]';
  vector[N_species] mu_r_log = [2, 2]';
  vector[N_species] mu_s_log = [-1.5, -1.5]';
  
  vector[N_species] sigma_b_log = [1, 1]';
  // vector[N_species] sigma_c_a_log = [1, 1]';
  vector[N_species] sigma_c_b_log = [2, 2]';
  vector[N_species] sigma_c_j_log = [1, 1]';
  vector[N_species] sigma_g_logit = [1, 1]';
  vector[N_species] sigma_h_logit = [0.5, 0.5]';
  // vector[N_species] sigma_l_log = [1, 1]';
  vector[N_species] sigma_r_log = [2, 2]';
  vector[N_species] sigma_s_log = [2, 2]';
  
  //// Data for separate fitting of the initial state
  // vector[N_pops] y0 [N_locs, N_plots] = y[ , , 1, ];
}


parameters {
  //// Level 1, (global, species): species-specific rates, indepndent of environment.
  // … dependent on environment.
  
  // matrix[N_beta, N_species] Beta_g; // (J-, A+) transition rate from J to A
  // matrix[N_beta, N_species] Beta_r; // // (J+) flow into the system, dependent on env
  // matrix[N_beta, N_species] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A
  // matrix[N_beta, N_species] Beta_l;


  // … independent of environment
  vector[N_species] b_log;
  // vector[N_species] c_a_log;
  vector[N_species] c_b_log;
  vector[N_species] c_j_log;
  vector[N_species] g_logit;
  vector[N_species] h_logit; // (A-, B+), here still unconstrained on log scale, flow out of the system, independent of environment ("intercept only")
  // vector[N_species] l_log;
  vector[N_species] r_log;
  vector[N_species] s_log;
  
  // vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  // vector<lower=0>[2] sigma_obs; // observation error
  vector<lower=0>[2] alpha_obs_inv; // observation error
  
  // matrix[N_pops, timespan_max] u[N_locs];

  vector[N_pops] state_init_log[N_locs];
  // vector[N_pops] state_init_log_raw[N_locs];
}


transformed parameters {
  
  vector[N_pops] y_hat[N_locs, N_times];
  // vector[N_pops] state_init_log[N_locs];

  
  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  // matrix[N_locs, N_species] g_logit = X * Beta_g;
  // matrix[N_locs, N_species] r_log = X * Beta_r;
  // matrix[N_locs, N_species] s_log = X * Beta_s;
  // matrix[N_locs, N_species] l_log = X * Beta_l;

  vector<lower=0>[2] alpha_obs = inv(alpha_obs_inv);

  for(loc in 1:N_locs) {
    
    // state_init_log[loc] = y0_loc_log[loc] + sigma_obs[rep_obsmethod2pops] .* state_init_log_raw[loc];
    
    y_hat[loc, ] = simulate(exp(state_init_log[loc]), time_max[loc], times[loc, ],
                              // exp(g_log[loc, ]'), exp(r_log[loc, ]'), exp(s_log[loc, ]'),
                              inv_logit(g_logit), exp(r_log), exp(s_log),
                              L_loc[loc, ]', // exp(l_log[loc, ]'),
                              exp(b_log), exp(c_j_log), exp(c_b_log), inv_logit(h_logit),
                              ba_a_avg, ba_a_upper,
                              N_species, N_pops,
                              // u[loc],
                              i_j, i_a, i_b);
    }

}


model {
  
  //---------- PRIORS ---------------------------------
  
  
  // sigma_process ~ normal(0, 0.01);
  // sigma_obs ~ normal(0, [0.5, 0.1]); // for observations from predictions
  alpha_obs_inv ~ normal(0, [0.1, 0.01]);
  
  // to_vector(Beta_g[2:N_beta,]) ~ std_normal();
  
  b_log   ~ normal(mu_b_log, sigma_b_log);
  c_b_log ~ normal(mu_c_b_log, sigma_c_b_log);
  c_j_log ~ normal(mu_c_j_log, sigma_c_j_log);
  g_logit ~ normal(mu_g_logit, sigma_g_logit);
  h_logit ~ normal(mu_h_logit, sigma_h_logit);
  r_log   ~ normal(mu_r_log, sigma_r_log);
  s_log   ~ normal(mu_s_log, sigma_s_log);
  
  // b_log_2nd ~ normal(mu_b_log, sigma_b_log);
  // c_b_log_2nd ~ normal(mu_c_b_log, sigma_c_b_log);
  // c_j_log_2nd ~ normal(mu_c_j_log, sigma_c_j_log);
  // g_logit_2nd ~ normal(mu_g_logit, sigma_g_logit);
  // h_logit_2nd ~ normal(mu_h_logit, sigma_h_logit);
  // r_log_2nd ~ normal(mu_r_log, sigma_r_log);
  // s_log_2nd ~ normal(mu_s_log, sigma_s_log);
  // 
  // Beta_b[1,] = b_log;
  // Beta_g[1,] = g_log;
  // Beta_h[1,] = h_log;
  // Beta_r[1,] = r_log;
  // Beta_s[1,] = s_log;
  // 
  // Beta_b[2,] ~ std_normal(); // or better smaller
  // Beta_g[2,] ~ std_normal();
  // Beta_h[2,] ~ std_normal();
  // Beta_r[2,] ~ std_normal();
  // Beta_s[2,] ~ std_normal();
  // 
  // Beta_b[3,] = b_log_2nd;
  // Beta_g[3,] = g_logit_2nd;
  // Beta_h[3,] = h_logit_2nd;
  // Beta_r[3,] = r_log_2nd;
  // Beta_s[3,] = s_log_2nd;

  
  //---------- MODEL ---------------------------------
  
  for(loc in 1:N_locs) {

    // Some priors
    // normal: // state_init[loc] ~ normal(y0_loc[loc], sigma_obs[rep_obsmethod2pops]);
    // state_init_log_raw[loc] ~ std_normal();
    //  state_init_log[loc] ~ normal(y0_loc_log[loc], 0.1);
    
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
  array[N_locs, N_times, N_plots, N_pops] real y_sim ;
  vector[N_pops] y_hat_rep[N_locs, N_times, N_plots]; // mainly for recovery test purposes
  
  array[N_locs] matrix[N_pops, N_pops] Jacobian_3; // N_stages * N_species. because each has it's own function
  array[N_locs] matrix[N_pops, N_pops] Jacobian_fix; // N_stages * N_species. because each has it's own function
  array[N_locs] real rho_3;
  array[N_locs] real rho_fix;
  array[N_locs] int converged; // tolerance has been reached
  array[N_locs] int convergent; // rho < 1
  array[N_locs] real iterations_fix;
  array[N_locs] vector[N_genstates+N_species+1] state_fix; // state_fix is a vector [J1, …, A1, …, B1, …, BA1, …, eps_ba1, …, iterations]
  array[N_locs] int dominant_fix;
  array[N_locs] int major_fix;

  // array[N_locs] matrix[N_pops, N_pops] Jacobian_c_b_0;
  // array[N_locs] int converged_c_b_0;
  array[N_locs] vector[N_genstates+N_species+1] state_fix_c_b_0;
  array[N_locs] int dominant_fix_c_b_0;
  array[N_locs] int major_fix_c_b_0;
  
  array[N_locs] vector[N_pops+N_species+1] state_fix_g_half;
  array[N_locs] vector[N_pops+N_species] contribution_g_half;
  
  array[N_locs] vector[N_genstates+N_species+1] state_fix_r_0;
  array[N_locs] vector[N_genstates] contribution_r;
  
  // array[N_locs] matrix[N_pops, N_pops] Jacobian_s_0;
  // array[N_locs] int converged_s_0;
  array[N_locs] vector[N_genstates+N_species+1] state_fix_s_0;
  array[N_locs] int dominant_fix_s_0;
  array[N_locs] int major_fix_s_0;
  array[N_locs] vector[N_genstates] contribution_s;
  

  for(loc in 1:N_locs) {
    

    Jacobian_3[loc] = jacobian(y_hat[loc, N_times, i_j[1]], y_hat[loc, N_times, i_j[2]], y_hat[loc, N_times, i_a[1]], y_hat[loc, N_times, i_a[2]], y_hat[loc, N_times, i_b[1]], y_hat[loc, N_times, i_b[2]],
                      exp(b_log[1]), exp(b_log[2]),   exp(c_b_log[1]), exp(c_b_log[2]),   exp(c_j_log[1]), exp(c_j_log[2]),   inv_logit(g_logit[1]), inv_logit(g_logit[2]),   inv_logit(h_logit[1]), inv_logit(h_logit[2]),  L_loc[loc, 1],  L_loc[loc, 2],   exp(r_log[1]), exp(r_log[2]),   exp(s_log[1]), exp(s_log[2]),
                      ba_a_avg, ba_a_upper);
    rho_3[loc] = norm(Jacobian_3[loc]);
    
    
    //// fix point, given parameters
    state_fix[loc] = iterateFix(y_hat[loc, N_times], // use the third time as initial value
                                   inv_logit(g_logit), exp(r_log), exp(s_log),
                                   L_loc[loc, ]', // exp(l_log[loc, ]'),
                                   exp(b_log), exp(c_j_log), exp(c_b_log), inv_logit(h_logit),
                                   ba_a_avg, ba_a_upper,
                                   N_species, N_pops,
                                   i_j, i_a, i_b,
                                   tolerance_fix);
                                   
    iterations_fix[loc] = state_fix[loc, N_genstates+N_species+1];
    converged[loc] = iterations_fix[loc] < 5000;

    Jacobian_fix[loc] = jacobian(state_fix[loc, i_j[1]], state_fix[loc, i_j[2]], state_fix[loc, i_a[1]], state_fix[loc, i_a[2]], state_fix[loc, i_b[1]], state_fix[loc, i_b[2]],
                      exp(b_log[1]), exp(b_log[2]), exp(c_b_log[1]), exp(c_b_log[2]), exp(c_j_log[1]), exp(c_j_log[2]), inv_logit(g_logit[1]), inv_logit(g_logit[2]),   inv_logit(h_logit[1]), inv_logit(h_logit[2]),  L_loc[loc, 1],  L_loc[loc, 2],   exp(r_log[1]), exp(r_log[2]),   exp(s_log[1]), exp(s_log[2]),
                      ba_a_avg, ba_a_upper);
    
    rho_fix[loc] = norm(Jacobian_fix[loc]);
    convergent[loc] = rho_fix[loc] < 1;
    
    
    if (converged[loc]) { # && convergent[loc]
      
      dominant_fix[loc] = state_fix[loc, N_pops+1]/state_fix[loc, N_genstates] > 3; // BA_1 > 75%
      major_fix[loc] = state_fix[loc, N_pops+1] > state_fix[loc, N_genstates]; // BA_1 > 50%
      
      
      //// ... given c_b == 0
      state_fix_c_b_0[loc] = iterateFix(state_fix[loc, 1:N_pops], // use the fixed point as initial value
                                    inv_logit(g_logit), exp(r_log), exp(s_log),
                                    L_loc[loc, ]', // exp(l_log[loc, ]'),
                                    exp(b_log), exp(c_j_log), [0.0, 0.0]', inv_logit(h_logit),
                                    ba_a_avg, ba_a_upper,
                                    N_species, N_pops,
                                    i_j, i_a, i_b,
                                    tolerance_fix);
                                    
      dominant_fix_c_b_0[loc] = state_fix_c_b_0[loc, N_pops+1]/state_fix_c_b_0[loc, N_genstates] > 3; // BA_1 > 75%
      major_fix_c_b_0[loc] = state_fix_c_b_0[loc, N_pops+1] > state_fix_c_b_0[loc, N_genstates]; // BA_1 > 50%
      
      //// ... given g == 0.5*g
      state_fix_g_half[loc] = iterateFix(state_fix[loc, 1:N_pops], // use the fixed point as initial value
                                                 inv_logit(g_logit)*0.5, exp(r_log), exp(s_log),
                                                 L_loc[loc, ]', // exp(l_log[loc, ]'),
                                                 exp(b_log), exp(c_j_log), exp(c_b_log), inv_logit(h_logit),
                                                 ba_a_avg, ba_a_upper,
                                                 N_species, N_pops,
                                                 i_j, i_a, i_b,
                                                 tolerance_fix);
      
      contribution_g_half[loc] = state_fix[loc, 1:N_genstates] - state_fix_g_half[loc, 1:N_genstates];

      
      //// ... given r == 0
      state_fix_r_0[loc] = iterateFix(state_fix[loc, 1:N_pops], // use the fixed point as initial value
                                                 inv_logit(g_logit), [0.0, 0.0]', exp(s_log),
                                                 L_loc[loc, ]', // exp(l_log[loc, ]'),
                                                 exp(b_log), exp(c_j_log), exp(c_b_log), inv_logit(h_logit),
                                                 ba_a_avg, ba_a_upper,
                                                 N_species, N_pops,
                                                 i_j, i_a, i_b,
                                                 tolerance_fix);
      
      contribution_r[loc] = state_fix[loc, 1:N_genstates] - state_fix_r_0[loc, 1:N_genstates];


      //// ... given s == 0
      state_fix_s_0[loc] = iterateFix(state_fix[loc, 1:N_pops], // use the fixed point as initial value
                                                         inv_logit(g_logit), exp(r_log), [0.0, 0.0]',
                                                         L_loc[loc, ]', // exp(l_log[loc, ]'),
                                                         exp(b_log), exp(c_j_log), exp(c_b_log), inv_logit(h_logit),
                                                         ba_a_avg, ba_a_upper,
                                                         N_species, N_pops,
                                                         i_j, i_a, i_b,
                                                         tolerance_fix);
      
      dominant_fix_s_0[loc] = state_fix_s_0[loc, N_pops+1]/state_fix_s_0[loc, N_genstates] > 3; // BA_1 > 75%
      major_fix_s_0[loc] = state_fix_s_0[loc, N_pops+1] > state_fix_s_0[loc, N_genstates]; // BA_1 > 50%
      contribution_s[loc] = state_fix[loc, 1:N_genstates] - state_fix_s_0[loc, 1:N_genstates];

    } else {
      
      print("System has not converged.");
    
    }
    
    // still nested in locs, the simulated states
    for (t in 1:N_times) {

      for (p in 1:N_plots) {

        y_hat_rep[loc, t, p,] = y_hat[loc, t,];
        
        // y_sim[loc, t, p, ] = normal_rng(y_hat[loc, t], sigma_obs[rep_obsmethod2pops]);
        y_sim[loc, t, p, ] = gamma_rng(alpha_obs[rep_obsmethod2pops], alpha_obs[rep_obsmethod2pops]./y_hat[loc, t]);

      }
    }
  }

}
