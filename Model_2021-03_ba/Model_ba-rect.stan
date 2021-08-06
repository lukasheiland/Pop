functions {
  
  //// Difference equations
  vector[] simulate(vector initialstate, int time_max, int[] times,
                  vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                  vector ba_a_avg, real ba_a_upper,
                  int N_spec, int N_pops,
                  // matrix u, // a matrix[N_pops, time_max-1]
                  int[] i_j, int[] i_a, int[] i_b) {
    
    
    // State array
    array[time_max] vector[N_pops] State;
    State[1,] = initialstate;
    // print("State 1 before: ", State[1,]);

    for (t in 2:time_max) {
      // Structure of state[N_pops]: stage/species
      // State has times rows (array row-major),  u has times cols (matrix: col-major access)
      
      vector[N_spec] J = State[t-1, i_j];
      vector[N_spec] A = State[t-1, i_a];
      vector[N_spec] B = State[t-1, i_b];

      vector[N_spec] BA = A .* ba_a_avg + B;
      real BA_sum = sum(BA);

      State[t, i_j]  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      State[t, i_a]  =  (g .* J + (A - h .*A )) ./ (1 + c_a*BA_sum);
      State[t, i_b]  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
      
    }
    
    // print("State 1 after: ", State[1,]);

    return State[times,];
  }
  
  //// Gets the specific Jacobian, given states and parameters
  // outdated: uses c_b not c_a for A
  // matrix jacobian(real J1, real J2, real A1, real A2, real B1, real B2,
  //               real b1, real b2, real c_b1, real c_b2, real c_j1, real c_j2, real g1, real g2, real h1, real h2, real l1, real l2, real r1, real r2, real s1, real s2,
  //               real p, real q) {
  //     
  //     matrix[6, 6] J;
  //     
  //     J = [
  //             [- (g1 - 1)/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1) - (c_j1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                             -(c_j1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2, (p*r1)/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1) - (p*s1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                         -(p*s1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2, r1/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1) - (s1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                     -(s1*(J1 + l1 - J1*g1 + B1*r1 + A1*p*r1))/(c_j1*(J1 + J2) + s1*(B1 + B2 + A1*p + A2*p) + 1)^2],
  //             [                                                            -(c_j2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2, - (g2 - 1)/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1) - (c_j2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                         -(p*s2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2, (p*r2)/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1) - (p*s2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                     -(s2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2, r2/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1) - (s2*(J2 + l2 - J2*g2 + B2*r2 + A2*p*r2))/(c_j2*(J1 + J2) + s2*(B1 + B2 + A1*p + A2*p) + 1)^2],
  //             [                                                                                                                    g1/(c_b1*(B1 + B2 + A1*p + A2*p) + 1),                                                                                                                                                       0,                                          - (h1 - 1)/(c_b1*(B1 + B2 + A1*p + A2*p) + 1) - (p*c_b1*(A1 - A1*h1 + J1*g1))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                       -(p*c_b1*(A1 - A1*h1 + J1*g1))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                   -(c_b1*(A1 - A1*h1 + J1*g1))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                   -(c_b1*(A1 - A1*h1 + J1*g1))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2],
  //             [                                                                                                                                                      0,                                                                                                                     g2/(c_b2*(B1 + B2 + A1*p + A2*p) + 1),                                                                                       -(p*c_b2*(A2 - A2*h2 + J2*g2))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                          - (h2 - 1)/(c_b2*(B1 + B2 + A1*p + A2*p) + 1) - (p*c_b2*(A2 - A2*h2 + J2*g2))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                   -(c_b2*(A2 - A2*h2 + J2*g2))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                                   -(c_b2*(A2 - A2*h2 + J2*g2))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2],
  //             [                                                                                                                                                      0,                                                                                                                                                       0,                     (h1*q + b1*h1*q)/(c_b1*(B1 + B2 + A1*p + A2*p) + 1) - (p*c_b1*(B1 + B1*b1 + A1*h1*q + A1*b1*h1*q))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                        -(p*c_b1*(B1 + B1*b1 + A1*h1*q + A1*b1*h1*q))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                         (b1 + 1)/(c_b1*(B1 + B2 + A1*p + A2*p) + 1) - (c_b1*(B1 + B1*b1 + A1*h1*q + A1*b1*h1*q))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                    -(c_b1*(B1 + B1*b1 + A1*h1*q + A1*b1*h1*q))/(c_b1*(B1 + B2 + A1*p + A2*p) + 1)^2],
  //             [                                                                                                                                                      0,                                                                                                                                                       0,                                                                        -(p*c_b2*(B2 + B2*b2 + A2*h2*q + A2*b2*h2*q))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                     (h2*q + b2*h2*q)/(c_b2*(B1 + B2 + A1*p + A2*p) + 1) - (p*c_b2*(B2 + B2*b2 + A2*h2*q + A2*b2*h2*q))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                                                                    -(c_b2*(B2 + B2*b2 + A2*h2*q + A2*b2*h2*q))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2,                         (b2 + 1)/(c_b2*(B1 + B2 + A1*p + A2*p) + 1) - (c_b2*(B2 + B2*b2 + A2*h2*q + A2*b2*h2*q))/(c_b2*(B1 + B2 + A1*p + A2*p) + 1)^2]
  //         ];
  //     
  //     return J;
  // }
  
  
  //// Transforms the vertex form into normal polynomial.
  vector[] transformToNormal(vector[] V) {
    
    // P is an: array[2] vector[N_beta], where array[1] is the mu, array[2] sigma;
    array[2] vector[5] P = V;
    
    // assumes vertex form f(x,y) == a_1*(x−p)^2 + a_2*(y−q)^2 + z
    // and input parameters vector[z, p, a_1, q, a_2]
    real z = V[1, 1];
    real p = V[1, 2];
    real a_1 = V[1, 3];
    real q = V[1, 4];
    real a_2 = V[1, 5]; // a == a
    
    // and output polynomial parameters vector[c, b_1, a_1, b_2, a_2]
    // as in f(x, y) = a_1*x^2 + b_1*x + a_2*y^2 + b_2*y + c
    P[1, 1] = a_2*q^2 + a_1*p^2 + z; // replace z with with c
    P[1, 2] = -2*a_1*p; // replace p and q with b1 and b2
    P[1, 4] = -2*a_2*q;
    
    return P;
  }
  
  //// Transforms polynomial form into vertex form
  vector[] transformToVertex(matrix P) {
    
    // P is a matrix[N_beta, N_species]
    array[2] vector[5] V;
    
    for (i in 1:2) {
      
      real c = P[1, i];
      real b_1 = P[2, i];
      real a_1 = P[3, i];
      real b_2 = P[4, i];
      real a_2 = P[5, i];
      
      // assumes vertex form f(x,y) == a_1*(x−p)^2 + a_2*(y−q)^2 + z
      // output vector[z, p, a_1, q, a_2]
      V[i, 2] = -b_1 / (2*a_1); // p
      V[i, 4] = -b_2 / (2*a_2); // q
      V[i, 1] = c - a_2*V[i, 4]^2 - a_1*V[i, 2]^2; // z
      V[i, 3] = a_1;
      V[i, 5] = a_2;
    }
    return V;
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
                    vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                    vector ba_a_avg, real ba_a_upper,
                    int N_spec, int N_pops,
                    int[] i_j, int[] i_a, int[] i_b,
                    real tolerance_fix, int fixiter_max) {
                       
    vector[N_pops+N_spec] s_0 = append_row(state_0, [0, 0]'); // two additional states for BA
    
    //// initialize while loop conditions
    vector[N_pops+N_spec] s_1; 
    vector[N_spec] eps_ba = [1.0, 1.0]'; // tolerance_fix is set to <upper=0.5>, that's why it is enough to set it to one for the while loop to run
    int i = 0;
    // int notconvergent = 1;
    
    
    while ( max(eps_ba) > tolerance_fix && i < fixiter_max ) { // if notconvergent were a good criterion: (notconvergent && max(eps_ba) > tolerance_fix)
      
      vector[N_spec] J;
      vector[N_spec] A;
      vector[N_spec] B;
      
      J = s_0[i_j];
      A = s_0[i_a];
      B = s_0[i_b];
      
      vector[N_spec] BA = A .* ba_a_avg + B;
      vector[N_spec] BA_1;
      real BA_sum = sum(BA);
      
      s_1[i_j]  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      s_1[i_a]  =  (g .* J + (A - h .*A )) ./ (1 + c_a*BA_sum);
      s_1[i_b]  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
      
      BA_1 =  s_1[i_a] .* ba_a_avg + s_1[i_b]; // New BA as additional state.
      s_1[(N_pops+1):] = BA_1;
      
      // notconvergent = (1 <= norm(jacobian(s_1[i_j[1]], s_1[i_j[2]], s_1[i_a[1]], s_1[i_a[2]], s_1[i_b[1]], s_1[i_b[2]], b[1], b[2], c_b[1], c_b[2], c_j[1], c_j[2], g[1], g[2], h[1], h[2], l[1], l[2], r[1], r[2], s[1], s[2], ba_a_avg, ba_a_upper)) );
      eps_ba = fabs((BA_1 - BA) ./ BA_1);
      s_0 = s_1;
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
  
  int<lower=0> N_beta; // number of environmental effects: degree 2 polynomial + intercept == (1 + 2 * n_env) == 5
  

  // obsmethod — factor (1, 2)
  array[N_pops] int<lower=1> rep_obsmethod2pops; // 1 1 1 1 2 2 2 2 2 2 2 2

  array[N_species] int<lower=1> i_j; // e.g. 1:4
  array[N_species] int<lower=1> i_a; // e.g. 5:9, etc.
  array[N_species] int<lower=1> i_b;


  
  array[N_locs, N_times] int times; // assumes start at 1. (Even though the number of times is the same, the times can deviate)
  array[N_locs] int time_max;
  int timespan_max; // max(time_max) - time_globalmin

  matrix[N_locs, N_beta] X; // design matrix
  array[N_locs] vector<lower=0>[N_species] L_smooth;

  // The response.
  // array[N_locs] vector[N_pops] y0_loc;
  array[N_locs, N_times, N_plots] vector[N_pops] y;
  
  //// Priors. The 2 reflect the two parameters mu and sigma
  // environmentally-dependent priors are species-agnostic on purpose
  // assumes vertex form f(x,y) == a_1*(x−p)^2 + a_2*(y−q)^2 + z
  // and parameters vector[z, p, q, a_1, a_2]
  array[2] vector[N_beta] prior_Vertex_c_j;
  array[2] vector[N_beta] prior_Vertex_g;
  array[2] vector[N_beta] prior_Vertex_r;
  array[2] vector[N_beta] prior_Vertex_s;
  
  array[2] vector[N_species] prior_b_log;
  array[2] vector[N_species] prior_c_a_log;
  array[2] vector[N_species] prior_c_b_log;
  array[2] vector[N_species] prior_h_logit;
  
  array[2] vector[N_species] prior_l_log;
  
  // Settings
  real<upper=0.5> tolerance_fix;
  // real dbh_lower_a; // 100
  // real dbh_lower_b; // 200
  real ba_a_upper;
  vector[N_species] ba_a_avg;
}


transformed data {
  
  //// This is dealt with empirically in R
  // real ba_a_upper = pi() * (dbh_lower_b/2)^2 * 1e-6; // / pi*r^2, mm^2 to m^2
  // real ba_a_avg = pi() * ((dbh_lower_a + dbh_lower_b)/2/2)^2 * 1e-6;
  
  // array[N_locs] vector[N_pops] y0_loc_log = log(y0_loc);
  
  // priors
  array[2] vector[N_beta] prior_Beta_c_j = transformToNormal(prior_Vertex_c_j);
  array[2] vector[N_beta] prior_Beta_g = transformToNormal(prior_Vertex_g);
  array[2] vector[N_beta] prior_Beta_r = transformToNormal(prior_Vertex_r);
  array[2] vector[N_beta] prior_Beta_s = transformToNormal(prior_Vertex_s);
  
  //// Data for separate fitting of the initial state
  // vector[N_pops] y0 [N_locs, N_plots] = y[ , , 1, ];
  
  //// Data for generated quantities
  int N_genstates = N_pops + N_species;
}


parameters {
  //// Level 1, (global, species): species-specific rates, indepndent of environment.
  
  // … independent of environment
  vector[N_species] b_log;
  vector[N_species] c_a_log;
  vector[N_species] c_b_log;
  vector[N_species] h_logit;
  // vector[N_species] c_a_log;
  
  // … dependent on environment. Matrix for convenient matrix multiplication
  matrix[N_beta, N_species] Beta_c_j;
  matrix[N_beta, N_species] Beta_g; // (J-, A+) transition rate from J to A
  matrix[N_beta, N_species] Beta_r; // // (J+) flow into the system, dependent on env
  matrix[N_beta, N_species] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A
  
  // Special case l
  vector[N_species] l_log;
  matrix[N_locs, N_species] L_random; // array[N_locs] vector[N_species] L_random; // here real array is used for compatibility with to_vector
  vector<lower=0>[N_species] sigma_l;

  
  // Errors
  vector<lower=0>[2] alpha_obs_inv; // observation error
    // vector<lower=0>[2] sigma_obs; // observation error
    // vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  
  // matrix[N_pops, timespan_max] u[N_locs];
  
  array[N_locs] vector[N_pops] state_init_log;
    // vector[N_pops] state_init_log_raw[N_locs];
}


transformed parameters {
  
  array[N_locs] vector[N_species] L_loc;
  
  array[N_locs, N_times] vector[N_pops] y_hat;
  // array[N_locs] vector[N_pops] state_init_log;

  
  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  matrix[N_locs, N_species] C_j_log = X * Beta_c_j;
  matrix[N_locs, N_species] G_logit = X * Beta_g;
  matrix[N_locs, N_species] R_log = X * Beta_r;
  matrix[N_locs, N_species] S_log = X * Beta_s;

  vector<lower=0>[2] alpha_obs = inv(alpha_obs_inv);

  for(loc in 1:N_locs) {
    
    L_loc[loc, ] = exp(l_log) .* L_smooth[loc, ] + // Offset with fixed coefficient. "The smooth to real number coefficient"
                 sigma_l .* L_random[loc, ]'; // non-centered loc-level random effects
    
    
    //// For non-centered observation error
    // state_init_log[loc] = y0_loc_log[loc] + sigma_obs[rep_obsmethod2pops] .* state_init_log_raw[loc];
    
    y_hat[loc, ] = simulate(exp(state_init_log[loc]), time_max[loc], times[loc, ],
                            exp(b_log), exp(c_a_log), exp(c_b_log), exp(C_j_log[loc,]'), inv_logit(G_logit[loc,]'), inv_logit(h_logit), L_loc[loc, ], exp(R_log[loc,]'), exp(S_log[loc,]'),
                            ba_a_avg, ba_a_upper,
                            N_species, N_pops,
                            // u[loc],
                            i_j, i_a, i_b);
                              
    }

}


model {
  
  //---------- PRIORS ---------------------------------
  
  
  // sigma_pr ocess ~ normal(0, 0.01);
  // sigma_obs ~ normal(0, [0.5, 0.1]); // for observations from predictions
  
  //// Hyperpriors
  alpha_obs_inv ~ normal(0, [0.1, 0.01]); // Observation error
  
  // ... for special offset L
  to_vector(L_random) ~ std_normal(); // Random part around slope for l
  sigma_l ~ cauchy(0, 2); // Regularizing half-normal on sigma for random slope for l
  l_log ~ normal(prior_l_log[1,], prior_l_log[2,]);
  
  
  //// Priors on Parameters
  // prior_*[2, N_species]
  b_log   ~ normal(prior_b_log[1,], prior_b_log[2,]);
  c_a_log ~ normal(prior_c_a_log[1,], prior_c_a_log[2,]);
  c_b_log ~ normal(prior_c_b_log[1,], prior_c_b_log[2,]);
  h_logit ~ normal(prior_h_logit[1,], prior_c_b_log[2,]);
  
  // same priors for both species
  for (spec in 1:N_species) {
    // prior_Beta_*[2, N_beta]
    Beta_c_j[, spec] ~ normal(prior_Beta_c_j[1, ], prior_Beta_c_j[2, ]);
    Beta_g[, spec]  ~ normal(prior_Beta_g[1, ], prior_Beta_g[2, ]);
    Beta_r[, spec] ~ normal(prior_Beta_r[1, ], prior_Beta_r[2, ]);
    Beta_s[, spec] ~ normal(prior_Beta_s[1, ], prior_Beta_s[2, ]);
  }
  

  
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
  
  int fixiter_max = 5000;

  array[N_locs, N_times, N_plots, N_pops] real y_sim ;
  vector[N_pops] y_hat_rep[N_locs, N_times, N_plots]; // mainly for recovery test purposes
  
  // array[N_locs] matrix[N_pops, N_pops] Jacobian_last; // N_stages * N_species. because each has it's own function
  // array[N_locs] matrix[N_pops, N_pops] Jacobian_fix; // N_stages * N_species. because each has it's own function
  // array[N_locs] real rho_last; // spectral radius at the most recent time in the data, however this convergence criterion seems to not apply to our fixed point (sufficient but not necessary, see Notes.md)
  // array[N_locs] real rho_fix; // ... at the fixed point
  // array[N_locs] int convergent; // rho < 1
  
  // array[N_locs] matrix[N_pops, N_pops] Jacobian_s_0;
  // array[N_locs] int converged_s_0;
  
  array[N_locs] int converged; // tolerance has been reached
  array[N_locs] real iterations_fix;
  array[N_locs] vector[N_genstates+N_species+1] state_fix; // state_fix is a vector [J1, …, A1, …, B1, …, BA1, …, eps_ba1, …, iterations]
  array[N_locs] int dominant_fix;
  array[N_locs] int major_fix;

  //// Here, the idea was to set c_b to 0 but this will also lead to Inf pop growth without m
  // array[N_locs] vector[N_genstates+N_species+1] state_fix_c_b_0;
  // array[N_locs] int dominant_fix_c_b_0;
  // array[N_locs] int major_fix_c_b_0;
  
  array[N_locs] vector[N_pops+N_species+1] state_fix_g_half;
  array[N_locs] vector[N_pops+N_species] contribution_g_half;
  
  array[N_locs] vector[N_genstates+N_species+1] state_fix_r_0;
  array[N_locs] vector[N_genstates] contribution_r;
  

  array[N_locs] vector[N_genstates+N_species+1] state_fix_s_0;
  array[N_locs] int dominant_fix_s_0;
  array[N_locs] int major_fix_s_0;
  array[N_locs] vector[N_genstates] contribution_s;
  
  // Vertices of the quadratic plane with [env_1 | env_2 | parameter]
  array[N_species] vector[N_beta] Vertex_c_j = transformToVertex(Beta_c_j);
  array[N_species] vector[N_beta] Vertex_g = transformToVertex(Beta_g);
  array[N_species] vector[N_beta] Vertex_r = transformToVertex(Beta_r);
  array[N_species] vector[N_beta] Vertex_s = transformToVertex(Beta_s);

  for(loc in 1:N_locs) {
    
    // Jacobian_last[loc] = jacobian(y_hat[loc, N_times, i_j[1]], y_hat[loc, N_times, i_j[2]], y_hat[loc, N_times, i_a[1]], y_hat[loc, N_times, i_a[2]], y_hat[loc, N_times, i_b[1]], y_hat[loc, N_times, i_b[2]],
    //                   exp(b_log[1]), exp(b_log[2]),   exp(c_b_log[1]), exp(c_b_log[2]),   exp(c_j_log[1]), exp(c_j_log[2]),   inv_logit(g_logit[1]), inv_logit(g_logit[2]),   inv_logit(h_logit[1]), inv_logit(h_logit[2]),  L_loc[loc, 1],  L_loc[loc, 2],   exp(r_log[1]), exp(r_log[2]),   exp(s_log[1]), exp(s_log[2]),
    //                   ba_a_avg, ba_a_upper);
    // rho_last[loc] = norm(Jacobian_last[loc]);
    
    //// fix point, given parameters
    state_fix[loc] = iterateFix(y_hat[loc, N_times], // use the third time as initial value
                                exp(b_log), exp(c_a_log), exp(c_b_log), exp(C_j_log[loc,]'), inv_logit(G_logit[loc,]'), inv_logit(h_logit), L_loc[loc, ], exp(R_log[loc,]'), exp(S_log[loc,]'),
                                ba_a_avg, ba_a_upper,
                                N_species, N_pops,
                                i_j, i_a, i_b,
                                tolerance_fix, fixiter_max);
                                   
    iterations_fix[loc] = state_fix[loc, N_genstates+N_species+1];
    converged[loc] = iterations_fix[loc] < fixiter_max;
    
    // Jacobian_fix[loc] = jacobian(state_fix[loc, i_j[1]], state_fix[loc, i_j[2]], state_fix[loc, i_a[1]], state_fix[loc, i_a[2]], state_fix[loc, i_b[1]], state_fix[loc, i_b[2]],
    //                   exp(b_log[1]), exp(b_log[2]), exp(c_b_log[1]), exp(c_b_log[2]), exp(c_j_log[1]), exp(c_j_log[2]), inv_logit(g_logit[1]), inv_logit(g_logit[2]),   inv_logit(h_logit[1]), inv_logit(h_logit[2]),  L_loc[loc, 1],  L_loc[loc, 2],   exp(r_log[1]), exp(r_log[2]),   exp(s_log[1]), exp(s_log[2]),
    //                   ba_a_avg, ba_a_upper);
    // rho_fix[loc] = norm(Jacobian_fix[loc]);
    // convergent[loc] = rho_fix[loc] < 1;
    
    
    if (converged[loc]) { // && convergent[loc]
      
      dominant_fix[loc] = state_fix[loc, N_pops+1]/state_fix[loc, N_genstates] > 3; // BA_1 > 75%
      major_fix[loc] = state_fix[loc, N_pops+1] > state_fix[loc, N_genstates]; // BA_1 > 50%
      
      //// ... given g == 0.5*g
      state_fix_g_half[loc] = iterateFix(state_fix[loc, 1:N_pops], // use the fixed point as initial value
                                          exp(b_log), exp(c_a_log), exp(c_b_log), exp(C_j_log[loc,]'), inv_logit(G_logit[loc,]') * 0.5, inv_logit(h_logit), L_loc[loc, ], exp(R_log[loc,]'), exp(S_log[loc,]'),
                                          ba_a_avg, ba_a_upper,
                                          N_species, N_pops,
                                          i_j, i_a, i_b,
                                          tolerance_fix, fixiter_max);
      
      contribution_g_half[loc] = state_fix[loc, 1:N_genstates] - state_fix_g_half[loc, 1:N_genstates];

      
      //// ... given r == 0
      state_fix_r_0[loc] = iterateFix(state_fix[loc, 1:N_pops], // use the fixed point as initial value
                                      exp(b_log), exp(c_a_log), exp(c_b_log), exp(C_j_log[loc,]'), inv_logit(G_logit[loc,]'), inv_logit(h_logit), L_loc[loc, ], [0.0, 0.0]', exp(S_log[loc,]'),
                                      ba_a_avg, ba_a_upper,
                                      N_species, N_pops,
                                      i_j, i_a, i_b,
                                      tolerance_fix, fixiter_max);
      
      contribution_r[loc] = state_fix[loc, 1:N_genstates] - state_fix_r_0[loc, 1:N_genstates];


      //// ... given s == 0
      state_fix_s_0[loc] = iterateFix(state_fix[loc, 1:N_pops], // use the fixed point as initial value
                                      exp(b_log), exp(c_a_log), exp(c_b_log), exp(C_j_log[loc,]'), inv_logit(G_logit[loc,]'), inv_logit(h_logit), L_loc[loc, ], exp(R_log[loc,]'), [0.0, 0.0]',
                                      ba_a_avg, ba_a_upper,
                                      N_species, N_pops,
                                      i_j, i_a, i_b,
                                      tolerance_fix, fixiter_max);
      
      dominant_fix_s_0[loc] = state_fix_s_0[loc, N_pops+1]/state_fix_s_0[loc, N_genstates] > 3; // BA_1 > 75%
      major_fix_s_0[loc] = state_fix_s_0[loc, N_pops+1] > state_fix_s_0[loc, N_genstates]; // BA_1 > 50%
      contribution_s[loc] = state_fix[loc, 1:N_genstates] - state_fix_s_0[loc, 1:N_genstates];

    } else {
      
      print("System has not converged after ", fixiter_max, "iterations.");
    
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
