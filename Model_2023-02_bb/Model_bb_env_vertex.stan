functions {
  
  //// transformToNormal(): Transforms parameters ...
  // from the vertex form: vectors[N_species]
  // into N_species column vectors of normal polynomial coefficients: Beta[N_species, N_beta].
  matrix transformToNormal(vector m,
                           vector center_env1, vector center_env2,
                           vector spread_env1, vector spread_env2) {
    
    // Beta is an: matrix[N_beta, N_species]
    matrix[5, 2] Beta;
    
    // assumes vertex form f(x,y) == spread_env1*(env1 − center_env1)^2 + spread_env2*(env2 − center_env2)^2 + m
    //
    // transform to polynomial parameters vector[c, b_env1, spread_env1, b_env2, spread_env2]
    // as in f(x, y) = c + b_env1*env1 + spread_env1*env1^2 +  b_env2*env2 + spread_env2*env2^2
    
    Beta[1,] = to_row_vector(spread_env1 .* center_env1^2 + spread_env2 .* center_env2^2 + m); // intercept c
    Beta[2,] = to_row_vector(-2 * spread_env1 .* center_env1); // linear coef b_env1
    Beta[3,] = to_row_vector(spread_env1); // quadratic coef spread_env1
    Beta[4,] = to_row_vector(-2 * spread_env2 .* center_env2); // linear coef b_env2
    Beta[5,] = to_row_vector(spread_env2); // quadratic coef spread_env2
    
    return Beta;
  }


  
  //// simulate(): Difference equations of the JAB model
  matrix simulate(vector initialstate, int time_max,
                  vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                  vector ba_a_avg, real ba_a_upper,
                  int N_spec, int N_pops,
                  array[] int i_j, array[] int i_a, array[] int i_b) {
    
    // State matrix with [species, times]. Times columns is sensible in order to have col-major access in matrices and for to_vector() later on
    matrix[N_pops, time_max] State;
    State[,1] = initialstate;
    
    for (t in 2:time_max) {
      // Structure of state[N_pops]: stage/species
      
      vector[N_spec] J = State[i_j, t-1];
      vector[N_spec] A = State[i_a, t-1];
      vector[N_spec] B = State[i_b, t-1];

      vector[N_spec] BA = (A .* ba_a_avg) + B;
      real BA_sum = sum(BA);
      
      /// Model
      vector[N_spec] lim_J = (1 + c_j*sum(J) + s*BA_sum);
      State[i_j, t]  =  l + r .* BA + (J - g .* J) ./ lim_J;
      vector[N_spec] lim_A = (1 + c_a*BA_sum);
      State[i_a, t]  =  (g .* J) ./ lim_J + (A - h .* A) ./ lim_A;
      State[i_b, t]  =  (ba_a_upper * h .* A ./ lim_A) + (1 + b) .* B ./ (1 + c_b*BA_sum);
    
    }
    
    return State;
  }
  

  
  //// unpack(): simulation and data assignment to one long vector is wrapped into a function here because
  // - the transformed parameters block does not allow declaring integers (necessary for ragged data structure indexing),
  // - the model block does not alllow for declaring variables within a loop.
  vector unpack(array[] vector state_init, array[] int time_max, array[] int times,    //* vector unpack(array[] vector state_init_log, array[] int time_max, array[] int times,
                array[] vector B, array[] vector C_a, array[] vector C_b, array[] vector C_j, array[] vector G, array[] vector H, array[] vector L_loc, array[] vector R, array[] vector S,
                vector ba_a_avg, real ba_a_upper,
                array[] int n_obs, array[] int n_yhat,
                int N_species, int N_pops, int L_y, int N_locs,
                array[] int i_j, array[] int i_a, array[] int i_b) {

    int pos_times = 1; // segmenting times[L_y]
    int pos_yhat = 1;
    vector[L_y] y_hat;

    for (loc in 1:N_locs) {
      int n_o = n_obs[loc];
      int n_y = n_yhat[loc];
      
      //// Returns a matrix State[p, t]
      // print("r", R_log);
      matrix[N_pops, time_max[loc]] States =
                        
                        simulate(state_init[loc],
                                 time_max[loc],
                                 B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc], R[loc], S[loc],
                                 ba_a_avg, ba_a_upper,
                                 N_species, N_pops,
                                 i_j, i_a, i_b);
      
  
      // Flattening the matrix into a vector for the location and append it to y_hat local vector yhat[m], and then into function-global vector y_hat[L_y].
      // to_vector converts matrix to a column vector in column-major order.
      y_hat[pos_yhat:(pos_yhat - 1 + n_y)] =
                       to_vector(States[ , segment(times, pos_times, n_o)]); // only select columns with times in the data
      
      
      pos_times = pos_times + n_o;
      pos_yhat = pos_yhat + n_y;
    }

  return y_hat; // Structure: locations/observations/pops(==stages/species)
  }
    


  //// iterateFix(): Difference equations simulated up to the fix point (i.e. equilibrium) given a maximum tolerance over all states.
  // Expects a state vector[N_pops]
  // returns a state vector of the form [J1, J2, A1, A2, B1, B2, BA1, BA2, eps_BA1, eps_BA2, iterations]
  array[] vector iterateFix(vector state_0,
                            vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                            vector ba_a_avg, real ba_a_upper,
                            int N_spec, array[] int i_j, array[] int i_a, array[] int i_b,
                            real tolerance_fix, int fixiter_max, int fixiter_min, int N_fix) {
                       

    /// initialize while loop conditions
    vector[N_spec] J = state_0[i_j];
    vector[N_spec] A = state_0[i_a];
    vector[N_spec] B = state_0[i_b];
    
    vector[N_spec] J_1;
    vector[N_spec] A_1;
    vector[N_spec] B_1;
    vector[N_spec] BA_1;

    vector[N_spec] eps_ba = rep_vector(1.0, N_spec); // tolerance_fix is set to <upper=0.5>, that's why it is enough to set it to one for the while loop to run
    int i = 0;

    
    while ( i < fixiter_min || (max(eps_ba) > tolerance_fix && i < fixiter_max) ) {
            
      vector[N_spec] BA = A .* ba_a_avg + B;
      real BA_sum = sum(BA);
      
      
      vector[N_spec] lim_J = (1 + c_j*sum(J) + s*BA_sum);
      J_1  =  l + r .* BA + (J - g .* J) ./ lim_J;
      vector[N_spec] lim_A = (1 + c_a*BA_sum);
      A_1  =  (g .* J) ./ lim_J + (A - h .* A) ./ lim_A;
      B_1  =  (ba_a_upper * h .* A ./ lim_A) + (1 + b) .* B ./ (1 + c_b*BA_sum);

      
      BA_1 = A_1 .* ba_a_avg + B_1; // New BA as additional state.

      eps_ba = abs((BA_1 - BA) ./ BA_1);
      
      /// !
      J = J_1;
      A = A_1;
      B = B_1;
      
      i += 1;

    } // end while i < fixiter_max
    
    // array with 3 (states) + 1 (BA) + 1 (eps) + 1 (n_iter) + some variables (overall N_fix)
    array[N_fix] vector[N_spec] fix = {J_1, A_1, B_1, BA_1,
                                       eps_ba, rep_vector(i, N_spec) // int i gets cast to real
                                       };
                                    
    return fix;
  
  }
  
}


///////////////////////////////////////////////////////////////////////////////



data {

  //// N — number of observations/groups; L - lengths of ragged vectors
  int<lower=0> L_times; // locations/obsid
  int<lower=0> L_y;
  int<lower=0> N_locs; // overall number of locations. This is the major running variable to iterate over the model vectors.
  int<lower=0> N_species; // overall number of unique species across plot and locations (not nested!)
  int<lower=0> N_pops; // (species*stages) within loc; this is the length of initial values!
  int<lower=0> N_env;
  int<lower=0> N_beta; // no of coefficients with predictors, depends on N_env and the specified polynomial/model matrix
  int<lower=0> N_protocol;
  int<lower=0> N_protocolTax;
  int<lower=0> N_obsidPop;

  //// n - number of levels within locs for more ragged models
  array[N_locs] int<lower=0> n_obs; // n solution times within locs
  array[N_locs] int<lower=0> n_yhat; // number of pops*obs within loc

  //// i — indices of stages
  array[N_species] int<lower=1> i_j; // e.g., 1:4
  array[N_species] int<lower=1> i_a; // e.g., 5:9, etc.
  array[N_species] int<lower=1> i_b;

  //// rep - repeat indices within groups for broadcasting to the subsequent hierarchical levels
  //+ array[L_y] int<lower=1> rep_yhat2y;
  array[L_y] int<lower=1> rep_pops2y; // factor (1:6)
  array[L_y] int<lower=1> rep_protocolTax2y;
  // array[L_y] int<lower=1> rep_protocol2y; // factor (1:5)
  array[L_y] int<lower=1> rep_obsidPop2y; // factor (1:12)
    
  //// actual data
  array[N_locs] int time_max;
  array[L_times] int times; // locations/observations
  array[N_locs] vector[N_species] L_smooth_log;
  vector<lower=0>[L_y] offset_data;
  array[L_y] int y; // the response
  matrix[N_locs, N_beta] X; // model matrix for environmental effects

  //// Settings
  real<upper=0.5> tolerance_fix;
  real ba_a_upper;
  vector[N_species] ba_a_avg;
  int<lower=0,upper=1> generateposteriorq;

  //// Priors
  // Parameters for gamma prior of initial variables
  array[N_locs] vector<lower=0>[N_pops] alpha_init;
  array[N_locs] vector<lower=0>[N_pops] beta_init;
  vector<lower=0>[N_pops] upper_init; // The upper is provided for linear rescaling of the data to sample the parameter in [0, 1]

  // Non-species-specific priors
  vector[2] prior_b_log;
  vector[2] prior_c_a_log;
  vector[2] prior_c_b_log;
  vector[2] prior_c_j_log;
  vector[2] prior_g_log;
  vector[2] prior_h_log;
  vector[2] prior_l_log;
  vector[2] prior_r_log;
  vector[2] prior_s_log;

  // Species pecific priors, e.g.
  // array[2] vector[N_species] prior_g_log;
  // array[2] vector[N_species] prior_h_log;
  // array[2] vector[N_species] prior_l_log;
  // array[2] vector[N_species] prior_r_log;

}


///////////////////////////////////////////////////////////////////////////////



transformed data {
  
  // Times are all assumed to be shifted to start at 1!
  
  //// Data for generated quantities
  int N_fix = 6; // an array of vectors[N_species] { J, A, B, BA, eps, n_iter}
  
  real targetslope = 0.1;
  
  real<lower=0> prior_sigma_b_log_spread = lambert_w0(targetslope * exp(-prior_b_log[1]) / 2);
  real<lower=0> prior_sigma_c_a_log_spread = lambert_w0(targetslope * exp(-prior_c_a_log[1]) / 2);
  real<lower=0> prior_sigma_c_b_log_spread = lambert_w0(targetslope * exp(-prior_c_b_log[1]) / 2);
  real<lower=0> prior_sigma_c_j_log_spread = lambert_w0(targetslope * exp(-prior_c_j_log[1]) / 2);
  real<lower=0> prior_sigma_g_log_spread = lambert_w0(targetslope * exp(-prior_g_log[1]) / 2);
  real<lower=0> prior_sigma_h_log_spread = lambert_w0(targetslope * exp(-prior_h_log[1]) / 2);
  real<lower=0> prior_sigma_r_log_spread = lambert_w0(targetslope * exp(-prior_r_log[1]) / 2);
  real<lower=0> prior_sigma_s_log_spread = lambert_w0(targetslope * exp(-prior_s_log[1]) / 2);
  
}


///////////////////////////////////////////////////////////////////////////////



parameters {
  //// Model parameters  
  // Maxima of the parameters
  vector<upper=0>[N_species] b_log;
  vector[N_species] c_a_log;
  vector[N_species] c_b_log;
  vector[N_species] c_j_log;
  vector<upper=0>[N_species] g_log;
  vector<upper=0>[N_species] h_log;
  vector[N_species] l_log;
  vector[N_species] r_log;
  vector[N_species] s_log;
  
  // Centers of the parameters along environmental axes
  vector[N_species] b_log_center_env1;
  vector[N_species] c_a_log_center_env1;
  vector[N_species] c_b_log_center_env1;
  vector[N_species] c_j_log_center_env1;
  vector[N_species] g_log_center_env1;
  vector[N_species] h_log_center_env1;
  vector[N_species] r_log_center_env1;
  vector[N_species] s_log_center_env1;
  
  vector[N_species] b_log_center_env2;
  vector[N_species] c_a_log_center_env2;
  vector[N_species] c_b_log_center_env2;
  vector[N_species] c_j_log_center_env2;
  vector[N_species] g_log_center_env2;
  vector[N_species] h_log_center_env2;
  vector[N_species] r_log_center_env2;
  vector[N_species] s_log_center_env2;
  
  // Spread of the parameters along environmental axes
  vector<upper=0>[N_species] b_log_spread_env1_x; // bell-shaped
  vector<lower=0>[N_species] c_a_log_spread_env1_x;  // inverse bell-shaped
  vector<lower=0>[N_species] c_b_log_spread_env1_x; // inverse bell-shaped
  vector<lower=0>[N_species] c_j_log_spread_env1_x; // inverse bell-shaped
  vector<upper=0>[N_species] g_log_spread_env1_x; // bell-shaped
  vector<upper=0>[N_species] h_log_spread_env1_x; // bell-shaped
  vector<upper=0>[N_species] r_log_spread_env1_x; // bell-shaped
  vector<lower=0>[N_species] s_log_spread_env1_x; // inverse bell-shaped
  
  vector<upper=0>[N_species] b_log_spread_env2_x;
  vector<lower=0>[N_species] c_a_log_spread_env2_x;
  vector<lower=0>[N_species] c_b_log_spread_env2_x;
  vector<lower=0>[N_species] c_j_log_spread_env2_x;
  vector<upper=0>[N_species] g_log_spread_env2_x;
  vector<upper=0>[N_species] h_log_spread_env2_x;
  vector<upper=0>[N_species] r_log_spread_env2_x;
  vector<lower=0>[N_species] s_log_spread_env2_x;


  //// Dispersion
  vector<lower=0>[N_protocolTax] phi_obs_inv; // error in neg_binomial per tax and stage
  //? vector<lower=0>[N_protocolTax] phi_obs_inv_sqrt;
  
  //// Initial state
  array[N_locs] vector<lower=0, upper=1>[N_pops] state_init_raw;

}


///////////////////////////////////////////////////////////////////////////////



transformed parameters {
  
  //// Local variables
  array[N_locs] vector<lower=0>[N_species] L_loc;
  array[N_locs] vector<lower=0>[N_pops] state_init;
  
  // Spread of the parameters along environmental axes
  vector<upper=0>[N_species] b_log_spread_env1 = b_log_spread_env1_x .* [1e-1, 1e-2]';
  vector<lower=0>[N_species] c_a_log_spread_env1 = c_a_log_spread_env1_x .* [1e-1, 1e-1]';
  vector<lower=0>[N_species] c_b_log_spread_env1 = c_b_log_spread_env1_x .* [1e-1, 1e-2]';
  vector<lower=0>[N_species] c_j_log_spread_env1 = c_j_log_spread_env1_x * 1e-1;
  vector<upper=0>[N_species] g_log_spread_env1 = g_log_spread_env1_x * 1e-1;
  vector<upper=0>[N_species] h_log_spread_env1 = h_log_spread_env1_x .* [1e-1, 1e-2]';
  vector<upper=0>[N_species] r_log_spread_env1 = r_log_spread_env1_x * 1e-3;
  vector<lower=0>[N_species] s_log_spread_env1 = s_log_spread_env1_x .* [1e-1, 1e-1]';
  
  vector<upper=0>[N_species] b_log_spread_env2 = b_log_spread_env2_x .* [1e-1, 1e-1]';
  vector<lower=0>[N_species] c_a_log_spread_env2  = c_a_log_spread_env2_x .* [1e0, 1e-1]';
  vector<lower=0>[N_species] c_b_log_spread_env2 = c_b_log_spread_env2_x .* [1e-1, 1e-2]';
  vector<lower=0>[N_species] c_j_log_spread_env2 = c_j_log_spread_env2_x * 1e-1;
  vector<upper=0>[N_species] g_log_spread_env2 = g_log_spread_env2_x * 1e-1;
  vector<upper=0>[N_species] h_log_spread_env2 = h_log_spread_env2_x * 1e-1;
  vector<upper=0>[N_species] r_log_spread_env2 = r_log_spread_env2_x * 1e-3;
  vector<lower=0>[N_species] s_log_spread_env2 = s_log_spread_env2_x * 1e-1;
  
  //// Environmental effects
  matrix[N_beta, N_species] Beta_b_log = transformToNormal(b_log, b_log_center_env1, b_log_center_env2, b_log_spread_env1, b_log_spread_env2);
  matrix[N_beta, N_species] Beta_c_a_log = transformToNormal(c_a_log, c_a_log_center_env1, c_a_log_center_env2, c_a_log_spread_env1, c_a_log_spread_env2);
  matrix[N_beta, N_species] Beta_c_b_log = transformToNormal(c_b_log, c_b_log_center_env1, c_b_log_center_env2, c_b_log_spread_env1, c_b_log_spread_env2);
  matrix[N_beta, N_species] Beta_c_j_log = transformToNormal(c_j_log, c_j_log_center_env1, c_j_log_center_env2, c_j_log_spread_env1, c_j_log_spread_env2);
  matrix[N_beta, N_species] Beta_g_log = transformToNormal(g_log, g_log_center_env1, g_log_center_env2, g_log_spread_env1, g_log_spread_env2);
  matrix[N_beta, N_species] Beta_h_log = transformToNormal(h_log, h_log_center_env1, h_log_center_env2, h_log_spread_env1, h_log_spread_env2);
  matrix[N_beta, N_species] Beta_r_log = transformToNormal(r_log, r_log_center_env1, r_log_center_env2, r_log_spread_env1, r_log_spread_env2);
  matrix[N_beta, N_species] Beta_s_log = transformToNormal(s_log,
                                                           s_log_center_env1, s_log_center_env2,
                                                           s_log_spread_env1, s_log_spread_env2);

                                                           
  matrix[N_locs, N_species] B_log = X * Beta_b_log; // two matrices X[N_locs, N_beta] * Beta[N_beta, N_species] = P_log[N_locs, N_species]; the columns of X[,N_beta] times the rows of Beta[N_beta,] will be added up as the value per loc.
  matrix[N_locs, N_species] C_a_log = X * Beta_c_a_log;
  matrix[N_locs, N_species] C_b_log = X * Beta_c_b_log;
  matrix[N_locs, N_species] C_j_log = X * Beta_c_j_log;
  matrix[N_locs, N_species] G_log = X * Beta_g_log;
  matrix[N_locs, N_species] H_log = X * Beta_h_log;
  matrix[N_locs, N_species] R_log = X * Beta_r_log;
  matrix[N_locs, N_species] S_log = X * Beta_s_log;
  
  array[N_locs] vector[N_species] B;
  array[N_locs] vector[N_species] C_a;
  array[N_locs] vector[N_species] C_b;
  array[N_locs] vector[N_species] C_j;
  array[N_locs] vector[N_species] G;
  array[N_locs] vector[N_species] H;
  array[N_locs] vector[N_species] R;
  array[N_locs] vector[N_species] S;

  for(loc in 1:N_locs) {
    
    state_init[loc] = state_init_raw[loc] .* upper_init;

    L_loc[loc, ] = exp(l_log + L_smooth_log[loc, ]); /// l * L_smooth == exp(l_log + L_smooth_log)
    // L_loc[loc, ] = exp(L_log[loc, ] + L_smooth_log[loc, ]); ///** version with random L
    
    B[loc] = exp(B_log[loc]');
    C_a[loc] = exp(C_a_log[loc]');
    C_b[loc] = exp(C_b_log[loc]');
    C_j[loc] = exp(C_j_log[loc]');
    G[loc] = exp(G_log[loc]');
    H[loc] = exp(H_log[loc]');
    R[loc] = exp(R_log[loc]');
    S[loc] = exp(S_log[loc]');
    
  }
  
  vector<lower=0>[L_y] y_hat = unpack(state_init, time_max, times,
                                      B,  C_a, C_b, C_j, G, H, L_loc, R, S,
                                      ba_a_avg, ba_a_upper,
                                      n_obs, n_yhat,
                                      N_species, N_pops, L_y, N_locs, // fixed numbers
                                      i_j, i_a, i_b);
                                
  vector[L_y] y_hat_offset = y_hat .* offset_data;
  
  vector<lower=0>[N_protocolTax] phi_obs = inv(phi_obs_inv .*
                                               [1e-3, 1e-3,      1e0, 1e0, 1e0, 1e0,          1e-3, 1e-3, 1e-3, 1e-3,      1e-1, 1e-1,      1e-1, 1e-1]'
                                               //[F.J.1, o.J.1,  F.J.2, o.J.2, F.J.3, o.J.3,  F.A.1, o.A.1, F.B.1, o.B.1,  F.A.23, o.A.23,  F.B.23, o.B.23]
                                               );
  
  //? vector<lower=0>[N_protocolTax] phi_obs = inv_square(phi_obs_inv_sqrt .*
  //?                                                    [1e-2, 1e-2, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1]'
  //?                                                    //[F.J.init, o.J.init, F.J.2, o.J.2, F.J.3, o.J.3, F.A.init, o.A.init, F.B.init, o.B.init, F.A.23, o.A.23, F.B.23, o.B.23]
  //?                                                    );
  
  vector[L_y] phi_obs_rep = phi_obs[rep_protocolTax2y];
  
}


///////////////////////////////////////////////////////////////////////////////



model {

  //—————————————————————————————————————————————————————————————————————//
  // Priors       ------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  //// Hyperpriors
  phi_obs_inv ~ std_normal();
  //? phi_obs_inv_sqrt ~ std_normal();
  
  //// Priors for Parameters  
  b_log   ~ normal(prior_b_log[1], prior_b_log[2]);
  c_a_log ~ normal(prior_c_a_log[1], prior_c_a_log[2]);
  c_b_log ~ normal(prior_c_b_log[1], prior_c_b_log[2]);
  c_j_log ~ normal(prior_c_j_log[1], prior_c_j_log[2]);
  g_log   ~ normal(prior_g_log[1], prior_g_log[2]);
  h_log   ~ normal(prior_h_log[1], prior_h_log[2]);
  l_log   ~ normal(prior_l_log[1], prior_l_log[2]);
  r_log   ~ normal(prior_r_log[1], prior_r_log[2]);
  s_log   ~ normal(prior_s_log[1], prior_s_log[2]);
  
  
  //// Priors for optimum value of parameters
  // note that there are no vertex parameters for l in the following, because for l we fit an intercept only:
  b_log_center_env1   ~ std_normal();
  c_a_log_center_env1 ~ std_normal();
  c_b_log_center_env1 ~ std_normal();
  c_j_log_center_env1 ~ std_normal();
  g_log_center_env1   ~ std_normal();
  h_log_center_env1   ~ std_normal();
  r_log_center_env1   ~ std_normal();
  s_log_center_env1   ~ std_normal();

  b_log_center_env2   ~ std_normal();
  c_a_log_center_env2 ~ std_normal();
  c_b_log_center_env2 ~ std_normal();
  c_j_log_center_env2 ~ std_normal();
  g_log_center_env2   ~ std_normal();
  h_log_center_env2   ~ std_normal();
  r_log_center_env2   ~ std_normal();
  s_log_center_env2   ~ std_normal();

  //// Priors for spread of parameters
  // Caution: exponential(rate) while double_exponential(mean, scale == 1/rate)
  -b_log_spread_env1  ~  normal(0, prior_sigma_b_log_spread); // exponential(5.0); // negative!
  c_a_log_spread_env1 ~  normal(0, prior_sigma_c_a_log_spread); // exponential(5.0);
  c_b_log_spread_env1 ~  normal(0, prior_sigma_c_b_log_spread); // exponential(5.0);
  c_j_log_spread_env1 ~  normal(0, prior_sigma_c_j_log_spread); // exponential(5.0);
  -g_log_spread_env1  ~  normal(0, prior_sigma_g_log_spread); // exponential(5.0);
  -h_log_spread_env1  ~  normal(0, prior_sigma_h_log_spread); // exponential(5.0);
  -r_log_spread_env1  ~  normal(0, prior_sigma_r_log_spread); // exponential(5.0);
  s_log_spread_env1   ~  normal(0, prior_sigma_s_log_spread); // exponential(5.0);
  
  -b_log_spread_env2  ~  normal(0, prior_sigma_b_log_spread);
  c_a_log_spread_env2 ~  normal(0, prior_sigma_c_a_log_spread);
  c_b_log_spread_env2 ~  normal(0, prior_sigma_c_b_log_spread);
  c_j_log_spread_env2 ~  normal(0, prior_sigma_c_j_log_spread);
  -g_log_spread_env2  ~  normal(0, prior_sigma_g_log_spread); 
  -h_log_spread_env2  ~  normal(0, prior_sigma_h_log_spread); 
  -r_log_spread_env2  ~  normal(0, prior_sigma_r_log_spread); 
  s_log_spread_env2   ~  normal(0, prior_sigma_s_log_spread); 
  
  
  for(l in 1:N_locs) { 
    
    //// Prior for initial state
    state_init[l] ~ gamma(alpha_init[l], beta_init[l]); // state_init is just a linear transform. -> No Jacobian correction necessary.
    
  }

  
  //—————————————————————————————————————————————————————————————————————//
  // Model       -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  y ~ neg_binomial_2(y_hat_offset, phi_obs_rep);

}


///////////////////////////////////////////////////////////////////////////////



generated quantities {
  //—————————————————————————————————————————————————————————————————————//
  // Print for debugging -----------------------------------------------//
  //———————————————————————————————————————————————————————————————————//  
  print(c_j_log);
  

  //—————————————————————————————————————————————————————————————————————//
  // Prediction  -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    

  array[L_y] real<lower=0> y_sim;
  y_sim = neg_binomial_2_rng(y_hat_offset, phi_obs_rep);


  //—————————————————————————————————————————————————————————————————————//
  // Averages over locations  ------------------------------------------//
  //———————————————————————————————————————————————————————————————————//
  vector[N_pops] avg_state_init;
  for (p in 1:N_pops) avg_state_init[p] = mean(state_init[, p]);
  
  vector[N_species] avg_L_loc;

  vector[N_species] avg_b;
  vector[N_species] avg_c_a;
  vector[N_species] avg_c_b;
  vector[N_species] avg_c_j;
  vector[N_species] avg_g;
  vector[N_species] avg_h;
  vector[N_species] avg_r;
  vector[N_species] avg_s;

  for (s in 1:N_species) {
    
    avg_L_loc[s] = mean(L_loc[, s]);
    
    avg_b[s] = mean(B[, s]);
    avg_c_a[s] = mean(C_a[, s]);
    avg_c_b[s] = mean(C_b[, s]);
    avg_c_j[s] = mean(C_j[, s]);
    avg_g[s] = mean(G[, s]);
    avg_h[s] = mean(H[, s]);
    avg_r[s] = mean(R[, s]);
    avg_s[s] = mean(S[, s]);
  }


  //—————————————————————————————————————————————————————————————————————//
  // Priors                  -------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  ///// Priors -----------------------------------
  real b_log_prior = normal_rng(prior_b_log[1], prior_b_log[2]);
  real c_a_log_prior = normal_rng(prior_c_a_log[1], prior_c_a_log[2]);
  real c_b_log_prior = normal_rng(prior_c_b_log[1], prior_c_b_log[2]);
  real c_j_log_prior = normal_rng(prior_c_j_log[1], prior_c_j_log[2]);
  real<upper=0> g_log_prior = -sqrt(square(normal_rng(prior_g_log[1], prior_g_log[2])));
  real<upper=0> h_log_prior = -sqrt(square(normal_rng(prior_h_log[1], prior_h_log[2])));
  real l_log_prior = normal_rng(prior_l_log[1], prior_l_log[2]);
  // vector[N_species] l_log_prior = to_vector(normal_rng(prior_l_log[1,], prior_l_log[2,]));
  real r_log_prior = normal_rng(prior_r_log[1], prior_r_log[2]);  
  real s_log_prior = normal_rng(prior_s_log[1], prior_s_log[2]);
  
  
  //—————————————————————————————————————————————————————————————————————————//
  // Posterior quantities  -------------------------------------------------//
  //———————————————————————————————————————————————————————————————————————//
  
  
  //// Rate tests -------------------------------------
  // int greater_b = b_log[1] > b_log[2];
  // int greater_c_a = c_a_log[1] > c_a_log[2];
  // int greater_c_b = c_b_log[1] > c_b_log[2];
  // int greater_c_j = c_j_log[1] > c_j_log[2];
  // int greater_g = g_log[1] > g_log[2];
  // int greater_h = h_log[1] > h_log[2];
  // int greater_l = l_log[1] > l_log[2];
  // int greater_r = r_log[1] > r_log[2];
  // int greater_s = s_log[1] > s_log[2];
  
  
  //// Declarations of posterior quantites --------------------------
  // - as global variables (local variables in blocks won't be saved)
  // - directly initiated with zeroes or 9, so that there are never NaNs in generated quantities.

  array[N_locs] vector[N_species] J_init = rep_array(rep_vector(0.0, N_species), N_locs);
  array[N_locs] vector[N_species] A_init = J_init;
  array[N_locs] vector[N_species] B_init = J_init;
  array[N_locs] vector[N_species] J_fix = J_init;
  array[N_locs] vector[N_species] A_fix = J_init;
  array[N_locs] vector[N_species] B_fix = J_init;
  array[N_locs] vector[N_species] ba_init = J_init;
  array[N_locs] vector[N_species] ba_fix = J_init;
  
  array[N_locs] real ba_sum_init_loc = rep_array(0.0, N_locs);
  array[N_locs] real ba_sum_fix_loc = ba_sum_init_loc;
  array[N_locs] vector[N_species] ba_frac_init = J_init;
  array[N_locs] vector[N_species] ba_frac_fix = J_init;
  
  array[N_locs] vector[N_species] eps_ba_fix = J_init;
  array[N_locs] real iterations_fix = ba_sum_init_loc;

  int fixiter_max = 6000;
  int fixiter_min = 250;

  
  array[N_locs] int converged_fix = rep_array(9, N_locs); // tolerance has been reached
  

  // array[N_locs] int dominant_init = converged_fix;
  // array[N_locs] int dominant_fix = converged_fix;
  array[N_locs] int major_init = converged_fix;
  array[N_locs] int major_fix = converged_fix;
  
  
  //// Parameters after limitation -------------------------------------
  
  array[N_locs] vector[N_species] G_lim_init_log = J_init;
  array[N_locs] vector[N_species] H_lim_init_log = J_init;
  array[N_locs] vector[N_species] B_lim_init_log = J_init;
  array[N_locs] vector[N_species] G_lim_fix_log = J_init;
  array[N_locs] vector[N_species] H_lim_fix_log = J_init;
  array[N_locs] vector[N_species] B_lim_fix_log = J_init;

  
  //// Contributions with average initial value --------------------------
  //
  // array[N_locs] vector[N_species] eps_ba_fix_avg = J_init;
  // array[N_locs] real iterations_fix_avg = ba_sum_init_loc;
  //
  // array[N_locs] int converged_fix_avg = converged_fix; // tolerance has been reached
  
  
  
  //// Declarations of counterfactual posterior quantities --------------------------
  
  // 0. Fagus getting the s minimum of others
  array[N_locs] vector[N_species] ba_fix_other_s = J_init;
  vector[N_locs] ba_frac_fix_other_s = rep_vector(9.0, N_locs);
  array[N_locs] int major_fix_other_s = converged_fix;
  
  // 1. K.O. of all regeneration
  array[N_locs] vector[N_species] ba_fix_ko_b_l_r = J_init;
  // array[N_locs] vector[N_species] ba_fix_ko_b_l_r_ko = J_init;
  
  // 2. K.O. of all environmental variation
  array[N_locs] vector[N_species] ba_fix_ko_1_env_b          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_b_c_b      = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_c_a        = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_c_b        = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_c_j        = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_g          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_g_c_j_s    = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_h          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_h_c_a      = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_l          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_r          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_1_env_s          = J_init;
  
  array[N_locs] vector[N_species] ba_fix_ko_2_env_b          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_b_c_b      = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_c_a        = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_c_b        = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_c_j        = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_g          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_g_c_j_s    = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_h          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_h_c_a      = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_l          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_r          = J_init;
  array[N_locs] vector[N_species] ba_fix_ko_2_env_s          = J_init;
  
  
  vector[N_locs] ba_frac_fix_ko_1_env_b          = rep_vector(9.0, N_locs);
  vector[N_locs] ba_frac_fix_ko_1_env_b_c_b      = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_c_a        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_c_b        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_c_j        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_g          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_g_c_j_s    = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_h          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_h_c_a      = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_l          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_r          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_1_env_s          = ba_frac_fix_ko_1_env_b;
  
  vector[N_locs] ba_frac_fix_ko_2_env_b          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_b_c_b      = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_c_a        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_c_b        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_c_j        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_g          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_g_c_j_s    = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_h          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_h_c_a      = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_l          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_r          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_fix_ko_2_env_s          = ba_frac_fix_ko_1_env_b;
  

  vector[N_locs] ba_frac_diff_fix_ko_1_env_b          = ba_frac_fix_ko_1_env_b;
  // vector[N_locs] ba_frac_diff_fix_ko_1_env_b_other_s  = ba_frac_fix_ko_1_env_b; //#
  vector[N_locs] ba_frac_diff_fix_ko_1_env_b_c_b      = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_c_a        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_c_b        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_c_j        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_g          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_g_c_j_s    = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_h          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_h_c_a      = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_l          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_r          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_1_env_s          = ba_frac_fix_ko_1_env_b;
  
  vector[N_locs] ba_frac_diff_fix_ko_2_env_b          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_b_c_b      = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_c_a        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_c_b        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_c_j        = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_g          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_g_c_j_s    = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_h          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_h_c_a      = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_l          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_r          = ba_frac_fix_ko_1_env_b;
  vector[N_locs] ba_frac_diff_fix_ko_2_env_s          = ba_frac_fix_ko_1_env_b;
  
  
  array[N_locs] int major_fix_ko_1_env_b         = converged_fix;
  array[N_locs] int major_fix_ko_1_env_b_c_b     = converged_fix;
  array[N_locs] int major_fix_ko_1_env_c_a       = converged_fix;
  array[N_locs] int major_fix_ko_1_env_c_b       = converged_fix;
  array[N_locs] int major_fix_ko_1_env_c_j       = converged_fix;
  array[N_locs] int major_fix_ko_1_env_g         = converged_fix;
  array[N_locs] int major_fix_ko_1_env_g_c_j_s   = converged_fix;
  array[N_locs] int major_fix_ko_1_env_h         = converged_fix;
  array[N_locs] int major_fix_ko_1_env_h_c_a     = converged_fix;
  array[N_locs] int major_fix_ko_1_env_l         = converged_fix;
  array[N_locs] int major_fix_ko_1_env_r         = converged_fix;
  array[N_locs] int major_fix_ko_1_env_s         = converged_fix;
   
  array[N_locs] int major_fix_ko_2_env_b         = converged_fix;
  array[N_locs] int major_fix_ko_2_env_b_c_b     = converged_fix;
  array[N_locs] int major_fix_ko_2_env_c_a       = converged_fix;
  array[N_locs] int major_fix_ko_2_env_c_b       = converged_fix;
  array[N_locs] int major_fix_ko_2_env_c_j       = converged_fix;
  array[N_locs] int major_fix_ko_2_env_g         = converged_fix;
  array[N_locs] int major_fix_ko_2_env_g_c_j_s   = converged_fix;
  array[N_locs] int major_fix_ko_2_env_h         = converged_fix;
  array[N_locs] int major_fix_ko_2_env_h_c_a     = converged_fix;
  array[N_locs] int major_fix_ko_2_env_l         = converged_fix;
  array[N_locs] int major_fix_ko_2_env_r         = converged_fix;
  array[N_locs] int major_fix_ko_2_env_s         = converged_fix;

  //———————————————————————————————————————————————————————————————————//
  // Generate Posterior quantities conditioned on setting  -----------//
  //—————————————————————————————————————————————————————————————————//
  
  if (generateposteriorq) {
    
    
    array[N_locs] real J_sum_fix_loc; // additional quantity for generating d.-d. growth terms
    array[N_locs] real J_sum_init_loc; // additional quantity for generating d.-d. growth terms

  
    //// Fix point iteration -------------------------------------------
    for(loc in 1:N_locs) {
      
      J_init[loc] = state_init[loc, 1:N_species];
      A_init[loc] = state_init[loc, (N_species+1):(N_species+N_species)];
      B_init[loc] = state_init[loc, (N_pops-N_species+1):N_pops];
      
      ba_init[loc] = state_init[loc, (N_pops-N_species+1):N_pops] + // State B
                     ba_a_avg .* state_init[loc, (N_species+1):(N_species+N_species)]; // State A * ba
                     
      ba_sum_init_loc[loc] = sum(ba_init[loc]);
      J_sum_init_loc[loc] = sum(J_init[loc]);
      
      ba_frac_init[loc] = ba_init[loc] / ba_sum_init_loc[loc];
      
      //// Density-dependent growth terms at the loc level
      B_lim_init_log[loc] = log( B[loc] ./ (1 + C_b[loc] * ba_sum_init_loc[loc]) );
      G_lim_init_log[loc] = log( G[loc] ./ (1 + C_j[loc] * J_sum_init_loc[loc] + S[loc] * ba_sum_init_loc[loc]) );
      H_lim_init_log[loc] = log( H[loc] ./ (1 + C_a[loc] * ba_sum_init_loc[loc]) );

      //// Booleans at init
      // dominant_init[loc] = (ba_init[loc, 1]/ba_init[loc, 2]) > 3; // ba_1 > 75%
      major_init[loc] = ba_init[loc, 1] > ba_init[loc, 2]; // ba_1 > 50%
      
      
      //// Simulate fix point, given parameters
      array[N_fix] vector[N_species] Fix; // N_locs arrays of vectors[N_specices] { J, A, B, BA, eps, n_iter, 2 * 9 * diff_ko_parameter }
      Fix = iterateFix(state_init[loc],
                      B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc, ], R[loc], S[loc],
                      ba_a_avg, ba_a_upper,
                      N_species, i_j, i_a, i_b,
                      tolerance_fix, fixiter_max, fixiter_min, N_fix);
                                          
                                     
      iterations_fix[loc] = Fix[6, 1]; // the 6th element is the vector: [n_iter, n_iter]'
      converged_fix[loc] = iterations_fix[loc] < fixiter_max; // (i starts at 0), when fixiter_max is reached the model ran 5001 times
      
      if (converged_fix[loc]) { // && convergent[loc]
      
        //// unpack Fix
        J_fix[loc] = Fix[1];
        A_fix[loc] = Fix[2];
        B_fix[loc] = Fix[3];
        ba_fix[loc] = Fix[4];
        eps_ba_fix[loc] = Fix[5];        
        // Fix[6] is unpacked before
        
        
        //// ba calculations
        ba_sum_fix_loc[loc] = sum(ba_fix[loc]);
        J_sum_fix_loc[loc] = sum(J_fix[loc]);
        ba_frac_fix[loc] = ba_fix[loc] / ba_sum_fix_loc[loc];
        
        //// Booleans at fixpoint
        // dominant_fix[loc] = (ba_fix[loc, 1]/ba_fix[loc, 2]) > 3; // ba_1 > 75%
        major_fix[loc] = ba_fix[loc, 1] > ba_fix[loc, 2]; // ba_1 > 50%
        
        //// Density-dependent growth terms at the loc level
        B_lim_fix_log[loc] = log( B[loc] ./ (1 + C_b[loc] * ba_sum_fix_loc[loc]) );
        G_lim_fix_log[loc] = log( G[loc] ./ (1 + C_j[loc] * J_sum_fix_loc[loc] + S[loc] * ba_sum_fix_loc[loc]) ); 
        H_lim_fix_log[loc] = log( H[loc] ./ (1 + C_a[loc] * ba_sum_fix_loc[loc]) );
 
        
        //// Counterfactual fix point iteration --------------------------
        // 0. Fagus getting the other s //#
        matrix[N_beta, N_species] Beta_s_log_other_s = transformToNormal(s_log[{2, 2}], s_log_center_env1, s_log_center_env2, s_log_spread_env1, s_log_spread_env2);
        row_vector[N_species] s_log_other_s = X[loc,] * Beta_s_log_other_s;
        vector[N_species] other_s = exp(s_log_other_s');
        
        array[N_fix] vector[N_species] Fix_other_s = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], other_s, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        ba_fix_other_s[loc] = Fix_other_s[4];
        ba_frac_fix_other_s[loc] = ba_fix_other_s[loc, 1] / sum(ba_fix_other_s[loc]);
        major_fix_other_s[loc] = ba_fix_other_s[loc, 1] > ba_fix_other_s[loc, 2];
        
        // array[N_fix] vector[N_species] Fix_ko_1_env_b_other_s = iterateFix(state_init[loc], [ avg_b[1], B[loc, 2] ]', C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], other_s, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix); //#
        // vector[N_species] ba_fix_ko_1_env_b_other_s = Fix_ko_1_env_b_other_s[4];
        // real ba_frac_fix_ko_1_env_b_other_s = ba_fix_ko_1_env_b_other_s[1] / sum(ba_fix_ko_1_env_b_other_s);

        
        
        
        // 1. K.O. of all regeneration
        vector[N_species] ko_1_b = [0, B[loc, 2]]';
        vector[N_species] ko_2_b = [B[loc, 1], 0]';
        vector[N_species] ko_1_r = [0, R[loc, 2]]';
        vector[N_species] ko_2_r = [R[loc, 1], 0]';
        vector[N_species] ko_1_l = [0, L_loc[loc, 2]]';
        vector[N_species] ko_2_l = [L_loc[loc, 1], 0]';
        
        array[N_fix] vector[N_species] Fix_ko_1_b_l_r;
        array[N_fix] vector[N_species] Fix_ko_2_b_l_r;
        
        Fix_ko_1_b_l_r = iterateFix(state_init[loc],
                                    ko_1_b, C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], ko_1_l, ko_1_r, S[loc],
                                    ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        Fix_ko_2_b_l_r = iterateFix(state_init[loc],
                                    ko_2_b, C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], ko_2_l, ko_2_r, S[loc],
                                    ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        
        ba_fix_ko_b_l_r[loc] = [Fix_ko_2_b_l_r[4, 1], Fix_ko_1_b_l_r[4, 2]]'; // basal area of the species that is not knocked out, respectively [1,2]
        // ba_fix_ko_b_l_r_ko[loc] = [Fix_ko_1_b_l_r[4, 1], Fix_ko_2_b_l_r[4, 2]]'; // basal area of the species that is knocked out, respectively [1,2]
        
        
        // 2. K.O. of environmental variation
        array[N_fix] vector[N_species] Fix_ko_1_env_b = iterateFix(state_init[loc], [ avg_b[1], B[loc, 2] ]', C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_b_c_b = iterateFix(state_init[loc], [ avg_b[1], B[loc, 2] ]', C_a[loc], [ avg_c_b[1], C_b[loc, 2] ]', C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_c_a = iterateFix(state_init[loc], B[loc], [ avg_c_a[1], C_a[loc, 2] ]', C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_c_b = iterateFix(state_init[loc], B[loc], C_a[loc], [ avg_c_b[1], C_b[loc, 2] ]', C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_c_j = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], [ avg_c_j[1], C_j[loc, 2] ]', G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_g = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], [ avg_g[1], G[loc, 2] ]', H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_g_c_j_s = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], [ avg_c_j[1], C_j[loc, 2] ]', [ avg_g[1], G[loc, 2] ]', H[loc], L_loc[loc,], R[loc], [ avg_s[1], S[loc, 2] ]', ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_h = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], [ avg_h[1], H[loc, 2] ]', L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);        
        array[N_fix] vector[N_species] Fix_ko_1_env_h_c_a = iterateFix(state_init[loc], B[loc], [ avg_c_a[1], C_a[loc, 2] ]', C_b[loc], C_j[loc], G[loc], [ avg_h[1], H[loc, 2] ]', L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_l = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], [ avg_L_loc[1], L_loc[loc, 2] ]', R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_r = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], [ avg_r[1], R[loc, 2] ]', S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_1_env_s = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], [ avg_s[1], S[loc, 2] ]', ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        
        array[N_fix] vector[N_species] Fix_ko_2_env_b = iterateFix(state_init[loc], [ B[loc, 1], avg_b[2] ]', C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_b_c_b = iterateFix(state_init[loc], [ B[loc, 1], avg_b[2] ]', C_a[loc], [ C_b[loc, 1], avg_c_b[2] ]', C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_c_a = iterateFix(state_init[loc], B[loc], [ C_a[loc, 1], avg_c_a[2] ]', C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_c_b = iterateFix(state_init[loc], B[loc], C_a[loc], [ C_b[loc, 1], avg_c_b[2] ]', C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_c_j = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], [ C_j[loc, 1], avg_c_j[2] ]', G[loc], H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_g = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], [ G[loc, 1], avg_g[2] ]', H[loc], L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_g_c_j_s = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], [ C_j[loc, 1], avg_c_j[2] ]', G[loc], H[loc], L_loc[loc,], R[loc], [ S[loc, 1], avg_s[2] ]', ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_h = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], [ G[loc, 1], avg_g[2] ]', [ H[loc, 1], avg_h[2] ]', L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);        
        array[N_fix] vector[N_species] Fix_ko_2_env_h_c_a = iterateFix(state_init[loc], B[loc], [ C_a[loc, 1], avg_c_a[2] ]', C_b[loc], C_j[loc], G[loc], [ H[loc, 1], avg_h[2] ]', L_loc[loc,], R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_l = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], [ L_loc[loc, 1], avg_L_loc[2] ]', R[loc], S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_r = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], [ R[loc, 1], avg_r[2] ]', S[loc], ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        array[N_fix] vector[N_species] Fix_ko_2_env_s = iterateFix(state_init[loc], B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc,], R[loc], [ S[loc, 1], avg_s[2] ]', ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, 50, N_fix);
        
        ba_fix_ko_1_env_b[loc]          = Fix_ko_1_env_b[4];
        ba_fix_ko_1_env_b_c_b[loc]      = Fix_ko_1_env_b_c_b[4];
        ba_fix_ko_1_env_c_a[loc]        = Fix_ko_1_env_c_a[4];
        ba_fix_ko_1_env_c_b[loc]        = Fix_ko_1_env_c_b[4];
        ba_fix_ko_1_env_c_j[loc]        = Fix_ko_1_env_c_j[4];
        ba_fix_ko_1_env_g[loc]          = Fix_ko_1_env_g[4];
        ba_fix_ko_1_env_g_c_j_s[loc]    = Fix_ko_1_env_g_c_j_s[4];
        ba_fix_ko_1_env_h[loc]          = Fix_ko_1_env_h[4];
        ba_fix_ko_1_env_h_c_a[loc]      = Fix_ko_1_env_h_c_a[4];
        ba_fix_ko_1_env_l[loc]          = Fix_ko_1_env_l[4];
        ba_fix_ko_1_env_r[loc]          = Fix_ko_1_env_r[4];
        ba_fix_ko_1_env_s[loc]          = Fix_ko_1_env_s[4];
             
        ba_fix_ko_2_env_b[loc]          = Fix_ko_2_env_b[4];
        ba_fix_ko_2_env_b_c_b[loc]      = Fix_ko_2_env_b_c_b[4];
        ba_fix_ko_2_env_c_a[loc]        = Fix_ko_2_env_c_a[4];
        ba_fix_ko_2_env_c_b[loc]        = Fix_ko_2_env_c_b[4];
        ba_fix_ko_2_env_c_j[loc]        = Fix_ko_2_env_c_j[4];
        ba_fix_ko_2_env_g[loc]          = Fix_ko_2_env_g[4];
        ba_fix_ko_2_env_g_c_j_s[loc]    = Fix_ko_2_env_g_c_j_s[4];
        ba_fix_ko_2_env_h[loc]          = Fix_ko_2_env_h[4];
        ba_fix_ko_2_env_h_c_a[loc]      = Fix_ko_2_env_h_c_a[4];
        ba_fix_ko_2_env_l[loc]          = Fix_ko_2_env_l[4];
        ba_fix_ko_2_env_r[loc]          = Fix_ko_2_env_r[4];
        ba_fix_ko_2_env_s[loc]          = Fix_ko_2_env_s[4];
        
        ba_frac_fix_ko_1_env_b[loc]          = ba_fix_ko_1_env_b[loc, 1]          / sum(ba_fix_ko_1_env_b[loc]);
        ba_frac_fix_ko_1_env_b_c_b[loc]      = ba_fix_ko_1_env_b_c_b[loc, 1]      / sum(ba_fix_ko_1_env_b_c_b[loc]);
        ba_frac_fix_ko_1_env_c_a[loc]        = ba_fix_ko_1_env_c_a[loc, 1]        / sum(ba_fix_ko_1_env_c_a[loc]);
        ba_frac_fix_ko_1_env_c_b[loc]        = ba_fix_ko_1_env_c_b[loc, 1]        / sum(ba_fix_ko_1_env_c_b[loc]);
        ba_frac_fix_ko_1_env_c_j[loc]        = ba_fix_ko_1_env_c_j[loc, 1]        / sum(ba_fix_ko_1_env_c_j[loc]);
        ba_frac_fix_ko_1_env_g[loc]          = ba_fix_ko_1_env_g[loc, 1]          / sum(ba_fix_ko_1_env_g[loc]);
        ba_frac_fix_ko_1_env_g_c_j_s[loc]    = ba_fix_ko_1_env_g_c_j_s[loc, 1]    / sum(ba_fix_ko_1_env_g_c_j_s[loc]);
        ba_frac_fix_ko_1_env_h[loc]          = ba_fix_ko_1_env_h[loc, 1]          / sum(ba_fix_ko_1_env_h[loc]);
        ba_frac_fix_ko_1_env_h_c_a[loc]      = ba_fix_ko_1_env_h_c_a[loc, 1]      / sum(ba_fix_ko_1_env_h_c_a[loc]);
        ba_frac_fix_ko_1_env_l[loc]          = ba_fix_ko_1_env_l[loc, 1]          / sum(ba_fix_ko_1_env_l[loc]);
        ba_frac_fix_ko_1_env_r[loc]          = ba_fix_ko_1_env_r[loc, 1]          / sum(ba_fix_ko_1_env_r[loc]);
        ba_frac_fix_ko_1_env_s[loc]          = ba_fix_ko_1_env_s[loc, 1]          / sum(ba_fix_ko_1_env_s[loc]);

        ba_frac_fix_ko_2_env_b[loc]          = ba_fix_ko_2_env_b[loc, 1]          / sum(ba_fix_ko_2_env_b[loc]);
        ba_frac_fix_ko_2_env_b_c_b[loc]      = ba_fix_ko_2_env_b_c_b[loc, 1]      / sum(ba_fix_ko_2_env_b_c_b[loc]);
        ba_frac_fix_ko_2_env_c_a[loc]        = ba_fix_ko_2_env_c_a[loc, 1]        / sum(ba_fix_ko_2_env_c_a[loc]);
        ba_frac_fix_ko_2_env_c_b[loc]        = ba_fix_ko_2_env_c_b[loc, 1]        / sum(ba_fix_ko_2_env_c_b[loc]);
        ba_frac_fix_ko_2_env_c_j[loc]        = ba_fix_ko_2_env_c_j[loc, 1]        / sum(ba_fix_ko_2_env_c_j[loc]);
        ba_frac_fix_ko_2_env_g[loc]          = ba_fix_ko_2_env_g[loc, 1]          / sum(ba_fix_ko_2_env_g[loc]);
        ba_frac_fix_ko_2_env_g_c_j_s[loc]    = ba_fix_ko_2_env_g_c_j_s[loc, 1]    / sum(ba_fix_ko_2_env_g_c_j_s[loc]);
        ba_frac_fix_ko_2_env_h[loc]          = ba_fix_ko_2_env_h[loc, 1]          / sum(ba_fix_ko_2_env_h[loc]);
        ba_frac_fix_ko_2_env_h_c_a[loc]      = ba_fix_ko_2_env_h_c_a[loc, 1]      / sum(ba_fix_ko_2_env_h_c_a[loc]);
        ba_frac_fix_ko_2_env_l[loc]          = ba_fix_ko_2_env_l[loc, 1]          / sum(ba_fix_ko_2_env_l[loc]);
        ba_frac_fix_ko_2_env_r[loc]          = ba_fix_ko_2_env_r[loc, 1]          / sum(ba_fix_ko_2_env_r[loc]);
        ba_frac_fix_ko_2_env_s[loc]          = ba_fix_ko_2_env_s[loc, 1]          / sum(ba_fix_ko_2_env_s[loc]);
        
        real ba_frac_fix_loc_Fagus = ba_frac_fix[loc, 1];
        ba_frac_diff_fix_ko_1_env_b[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_b[loc];
        // ba_frac_diff_fix_ko_1_env_b_other_s[loc]  =  ba_frac_fix_other_s[loc] - ba_frac_fix_ko_1_env_b_other_s; //#
        ba_frac_diff_fix_ko_1_env_b_c_b[loc]      =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_b_c_b[loc];
        ba_frac_diff_fix_ko_1_env_c_a[loc]        =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_c_a[loc];
        ba_frac_diff_fix_ko_1_env_c_b[loc]        =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_c_b[loc];
        ba_frac_diff_fix_ko_1_env_c_j[loc]        =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_c_j[loc];
        ba_frac_diff_fix_ko_1_env_g[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_g[loc];
        ba_frac_diff_fix_ko_1_env_g_c_j_s[loc]    =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_g_c_j_s[loc];
        ba_frac_diff_fix_ko_1_env_h[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_h[loc];
        ba_frac_diff_fix_ko_1_env_h_c_a[loc]      =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_h_c_a[loc];
        ba_frac_diff_fix_ko_1_env_l[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_l[loc];
        ba_frac_diff_fix_ko_1_env_r[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_r[loc];
        ba_frac_diff_fix_ko_1_env_s[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_1_env_s[loc];
        
        ba_frac_diff_fix_ko_2_env_b[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_b[loc];
        ba_frac_diff_fix_ko_2_env_b_c_b[loc]      =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_b_c_b[loc];
        ba_frac_diff_fix_ko_2_env_c_a[loc]        =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_c_a[loc];
        ba_frac_diff_fix_ko_2_env_c_b[loc]        =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_c_b[loc];
        ba_frac_diff_fix_ko_2_env_c_j[loc]        =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_c_j[loc];
        ba_frac_diff_fix_ko_2_env_g[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_g[loc];
        ba_frac_diff_fix_ko_2_env_g_c_j_s[loc]    =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_g_c_j_s[loc];
        ba_frac_diff_fix_ko_2_env_h[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_h[loc];
        ba_frac_diff_fix_ko_2_env_h_c_a[loc]      =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_h_c_a[loc];
        ba_frac_diff_fix_ko_2_env_l[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_l[loc];
        ba_frac_diff_fix_ko_2_env_r[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_r[loc];
        ba_frac_diff_fix_ko_2_env_s[loc]          =  ba_frac_fix_loc_Fagus - ba_frac_fix_ko_2_env_s[loc];

        
        major_fix_ko_1_env_b[loc]          = ba_fix_ko_1_env_b[loc,1]        > ba_fix_ko_1_env_b[loc,2];
        major_fix_ko_1_env_b_c_b[loc]      = ba_fix_ko_1_env_b_c_b[loc,1]    > ba_fix_ko_1_env_b_c_b[loc,2];
        major_fix_ko_1_env_c_a[loc]        = ba_fix_ko_1_env_c_a[loc,1]      > ba_fix_ko_1_env_c_a[loc,2];
        major_fix_ko_1_env_c_b[loc]        = ba_fix_ko_1_env_c_b[loc,1]      > ba_fix_ko_1_env_c_b[loc,2];
        major_fix_ko_1_env_c_j[loc]        = ba_fix_ko_1_env_c_j[loc,1]      > ba_fix_ko_1_env_c_j[loc,2];
        major_fix_ko_1_env_g[loc]          = ba_fix_ko_1_env_g[loc,1]        > ba_fix_ko_1_env_g[loc,2];
        major_fix_ko_1_env_g_c_j_s[loc]    = ba_fix_ko_1_env_g_c_j_s[loc,1]  > ba_fix_ko_1_env_g_c_j_s[loc,2];
        major_fix_ko_1_env_h[loc]          = ba_fix_ko_1_env_h[loc,1]        > ba_fix_ko_1_env_h[loc,2];
        major_fix_ko_1_env_h_c_a[loc]      = ba_fix_ko_1_env_h_c_a[loc,1]    > ba_fix_ko_1_env_h_c_a[loc,2];
        major_fix_ko_1_env_l[loc]          = ba_fix_ko_1_env_l[loc,1]        > ba_fix_ko_1_env_l[loc,2];
        major_fix_ko_1_env_r[loc]          = ba_fix_ko_1_env_r[loc,1]        > ba_fix_ko_1_env_r[loc,2];
        major_fix_ko_1_env_s[loc]          = ba_fix_ko_1_env_s[loc,1]        > ba_fix_ko_1_env_s[loc,2];
        
        major_fix_ko_2_env_b[loc]          = ba_fix_ko_2_env_b[loc,1]        > ba_fix_ko_2_env_b[loc,2];
        major_fix_ko_2_env_b_c_b[loc]      = ba_fix_ko_2_env_b_c_b[loc,1]    > ba_fix_ko_2_env_b_c_b[loc,2];
        major_fix_ko_2_env_c_a[loc]        = ba_fix_ko_2_env_c_a[loc,1]      > ba_fix_ko_2_env_c_a[loc,2];
        major_fix_ko_2_env_c_b[loc]        = ba_fix_ko_2_env_c_b[loc,1]      > ba_fix_ko_2_env_c_b[loc,2];
        major_fix_ko_2_env_c_j[loc]        = ba_fix_ko_2_env_c_j[loc,1]      > ba_fix_ko_2_env_c_j[loc,2];
        major_fix_ko_2_env_g[loc]          = ba_fix_ko_2_env_g[loc,1]        > ba_fix_ko_2_env_g[loc,2];
        major_fix_ko_2_env_g_c_j_s[loc]    = ba_fix_ko_2_env_g_c_j_s[loc,1]  > ba_fix_ko_2_env_g_c_j_s[loc,2];
        major_fix_ko_2_env_h[loc]          = ba_fix_ko_2_env_h[loc,1]        > ba_fix_ko_2_env_h[loc,2];
        major_fix_ko_2_env_h_c_a[loc]      = ba_fix_ko_2_env_h_c_a[loc,1]    > ba_fix_ko_2_env_h_c_a[loc,2];
        major_fix_ko_2_env_l[loc]          = ba_fix_ko_2_env_l[loc,1]        > ba_fix_ko_2_env_l[loc,2];
        major_fix_ko_2_env_r[loc]          = ba_fix_ko_2_env_r[loc,1]        > ba_fix_ko_2_env_r[loc,2];
        major_fix_ko_2_env_s[loc]          = ba_fix_ko_2_env_s[loc,1]        > ba_fix_ko_2_env_s[loc,2];

      }
      
      
      //// Simulate marginal fix point with average population starts
      //
      // array[N_fix] vector[N_species] Fix_avg;
      // Fix_avg = iterateFix(avg_state_init, // !!!
      //                     B[loc], C_a[loc], C_b[loc], C_j[loc], G[loc], H[loc], L_loc[loc, ], R[loc], S[loc],
      //                     ba_a_avg, ba_a_upper,
      //                     N_species, i_j, i_a, i_b,
      //                     tolerance_fix, fixiter_max, fixiter_min, N_fix);
      //                                     
      //                                
      // iterations_fix_avg[loc] = Fix_avg[6, 1]; // the 6th element is the vector: [n_iter, n_iter]'
      // converged_fix_avg[loc] = iterations_fix_avg[loc] < fixiter_max; // (i starts at 0), when fixiter_max is reached the model ran 5001 times
      // 
      // if (converged_fix_avg[loc]) { // && convergent[loc]
      // 
      //   //// unpack Fix_avg
      //   // J_fix_avg[loc] = Fix[1];
      //   // A_fix_avg[loc] = Fix[2];
      //   // B_fix_avg[loc] = Fix[3];
      //   // ba_fix_avg[loc] = Fix[4];
      //   eps_ba_fix_avg[loc] = Fix[5];
      //   
      //   // Fix[6] is unpacked before
      //
      //
      // } // end if (converged_fix_avg[loc])
  
    } // end for(loc in 1:N_locs)

  } // end if(generateposteriorq)

}
