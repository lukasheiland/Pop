functions {
  
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
      State[i_j, t]  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      State[i_a, t]  =  (g .* J + (A - h .*A )) ./ (1 + c_a*BA_sum);
      State[i_b, t]  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
    
    }
    
    return State;
  }
  


  //// simulate_1(): Difference equations of one-step model
  vector simulate_1(vector J, vector A, vector B,
                    vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                    vector ba_a_avg, real ba_a_upper,
                    int N_spec) {
    
    vector[N_spec] BA = (A .* ba_a_avg) + B;
    real BA_sum = sum(BA);

    /// Model (1 iteration)
    /// First run: get limiting states
    vector[N_spec] J_1  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
    vector[N_spec] A_1  =  (g .* J_1 + (A - h .*A )) ./ (1 + c_a*BA_sum);
    vector[N_spec] B_1  =  (1+b) .* ((h .* A_1 * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
    
    real BA_sum_1 = sum(A_1 .* ba_a_avg + B_1);
        
    // Second run: use limiting states
    J_1  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J_1) + s*BA_sum_1); // use the states already generated before: J_1, BA_sum_1
    A_1  =  (g .* J_1 + (A - h .*A )) ./ (1 + c_a*BA_sum_1);
    B_1  =  (1+b).*((h .* A_1 * ba_a_upper) + B) ./ (1 + c_b*BA_sum_1);
    
    vector[N_spec] ba_1 = (A_1 .* ba_a_avg) + B_1;
    
    return ba_1;
  }
  
  
  
  //// unpack(): ODE integration and data assignment to one long vector is wrapped into a function here because
  // - the transformed parameters block does not allow declaring integers (necessary for ragged data structure indexing),
  // - the model block does not alllow for declaring variables within a loop.
  vector unpack(array[] vector state_init, array[] int time_max, array[] int times,    //* vector unpack(array[] vector state_init_log, array[] int time_max, array[] int times,
                // vector b_log, vector c_a_log, vector c_b_log, vector c_j_log, vector g_log, vector h_log, vector l_log, vector r_log, vector s_log,
                vector b_log, vector c_a_log, vector c_b_log, vector c_j_log, vector g_log, vector h_log, array[] vector L_loc, vector r_log, vector s_log,
                // vector b_log, vector c_a_log, vector c_b_log, matrix C_j_log, matrix G_log, vector h_log, array[] vector L_loc, matrix R_log, matrix S_log,
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
                        
                        simulate(state_init[loc], //* simulate(exp(state_init[loc]),
                                 time_max[loc],
                                 // exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), exp(l_log), exp(r_log), exp(s_log),
                                 exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log),
                                 // ... etc. ,
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
      
      // s_1[i_j]  =  ((r .* BA)/(1 + BA_sum) + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      J_1  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      A_1  =  (g .* J + (A - h .*A )) ./ (1 + c_a*BA_sum);
      B_1  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
      
      BA_1 = A_1 .* ba_a_avg + B_1; // New BA as additional state.

      eps_ba = fabs((BA_1 - BA) ./ BA_1);
      

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
    

  //// iterateFix(): Difference equations simulated up to the fix point (i.e. equilibrium) given a maximum tolerance over all states.
  // Expects a state vector[N_pops]
  // returns a state vector of the form [J1, J2, A1, A2, B1, B2, BA1, BA2, eps_BA1, eps_BA2, iterations]
  array[] vector iterateFix_contributions(vector state_0,
                                          vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                                          vector ba_a_avg, real ba_a_upper,
                                          int N_spec, array[] int i_j, array[] int i_a, array[] int i_b,
                                          real tolerance_fix, int fixiter_max, int fixiter_min, int N_fix) {
                       
        
    /// Summed up contributions, initialize with 0 to avoid NaNs
    vector[N_spec] sum_ko_1_b = rep_vector(0.0, N_spec);
    vector[N_spec] sum_ko_1_b_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_1_c_a = sum_ko_1_b;
    vector[N_spec] sum_ko_1_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_1_c_j = sum_ko_1_b;
    vector[N_spec] sum_ko_1_g = sum_ko_1_b;
    vector[N_spec] sum_ko_1_h = sum_ko_1_b;
    vector[N_spec] sum_ko_1_l = sum_ko_1_b;
    vector[N_spec] sum_ko_1_r = sum_ko_1_b;
    vector[N_spec] sum_ko_1_s = sum_ko_1_b;
    
    vector[N_spec] sum_ko_2_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_b_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_c_a = sum_ko_1_b;
    vector[N_spec] sum_ko_2_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_c_j = sum_ko_1_b;
    vector[N_spec] sum_ko_2_g = sum_ko_1_b;
    vector[N_spec] sum_ko_2_h = sum_ko_1_b;
    vector[N_spec] sum_ko_2_l = sum_ko_1_b;
    vector[N_spec] sum_ko_2_r = sum_ko_1_b;
    vector[N_spec] sum_ko_2_s = sum_ko_1_b;
    
    vector[N_spec] sum_ko_1_prop_b = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_b_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_c_a = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_c_j = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_g = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_h = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_l = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_r = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_s = sum_ko_1_b;
    
    vector[N_spec] sum_ko_2_prop_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_b_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_c_a = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_c_j = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_g = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_h = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_l = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_r = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_s = sum_ko_1_b;
    
    
    vector[N_spec] sum_switch_b = sum_ko_1_b;
    vector[N_spec] sum_switch_b_c_b = sum_ko_1_b;
    vector[N_spec] sum_switch_c_a = sum_ko_1_b;
    vector[N_spec] sum_switch_c_b = sum_ko_1_b;
    vector[N_spec] sum_switch_c_j = sum_ko_1_b;
    vector[N_spec] sum_switch_g = sum_ko_1_b;
    vector[N_spec] sum_switch_h = sum_ko_1_b;
    vector[N_spec] sum_switch_l = sum_ko_1_b;
    vector[N_spec] sum_switch_r = sum_ko_1_b;
    vector[N_spec] sum_switch_s = sum_ko_1_b;


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
      
      // s_1[i_j]  =  ((r .* BA)/(1 + BA_sum) + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      J_1  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      A_1  =  (g .* J + (A - h .*A )) ./ (1 + c_a*BA_sum);
      B_1  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
      
      BA_1 = A_1 .* ba_a_avg + B_1; // New BA as additional state.

      eps_ba = fabs((BA_1 - BA) ./ BA_1);
      
      //// This is here to only calculate the contributions for an earlier period
      if (i < fixiter_min) {
      
        vector[N_spec] ba_ko_none = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);

        vector[N_spec] ba_ko_1_b = simulate_1(J, A, B, [0, b[2]]', c_a, c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_b_c_b = simulate_1(J, A, B, [0, b[2]]', c_a, [0, c_b[2]]', c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_c_a  = simulate_1(J, A, B, b, [0, c_a[2]]', c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_c_b = simulate_1(J, A, B, b, c_a, [0, c_b[2]]', c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_c_j = simulate_1(J, A, B, b, c_a, c_b, [0, c_j[2]]', g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_g = simulate_1(J, A, B, b, c_a, c_b, c_j, [0, g[2]]', h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_h = simulate_1(J, A, B, b, c_a, c_b, c_j, g, [0, h[2]]', l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_l = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, [0, l[2]]', r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_r = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, [0, r[2]]', s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_s = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, r, [0, s[2]]', ba_a_avg, ba_a_upper, N_spec);
        
        vector[N_spec] ba_ko_2_b = simulate_1(J, A, B, [b[1], 0]', c_a, c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_b_c_b = simulate_1(J, A, B, [b[1], 0]', c_a, [c_b[1], 0]', c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_c_a  = simulate_1(J, A, B, b, [c_a[1], 0]', c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_c_b = simulate_1(J, A, B, b, c_a, [c_b[1], 0]', c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_c_j = simulate_1(J, A, B, b, c_a, c_b, [c_j[1], 0]', g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_g = simulate_1(J, A, B, b, c_a, c_b, c_j, [g[1], 0]', h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_h = simulate_1(J, A, B, b, c_a, c_b, c_j, g, [h[1], 0]', l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_l = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, [l[1], 0]', r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_r = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, [r[1], 0]', s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_s = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, r, [s[1], 0]', ba_a_avg, ba_a_upper, N_spec);
        

        vector[2] switch_b = [b[2], b[1]]'; // reversal of the parameter vector with range 2:1 does not seem to work
        vector[2] switch_c_a = [c_a[2], c_a[1]]';
        vector[2] switch_c_b = [c_b[2], c_b[1]]';
        vector[2] switch_c_j = [c_j[2], c_j[1]]';
        vector[2] switch_g = [g[2], g[1]]';
        vector[2] switch_h = [h[2], h[1]]';
        vector[2] switch_l = [l[2], l[1]]';
        vector[2] switch_r = [r[2], r[1]]';
        vector[2] switch_s = [s[2], s[1]]';
        
        vector[N_spec] ba_switch_b = simulate_1(J, A, B, switch_b, c_a, c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_b_c_b = simulate_1(J, A, B, switch_b, c_a, switch_c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_c_a  = simulate_1(J, A, B, b, switch_c_a, c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_c_b = simulate_1(J, A, B, b, c_a, switch_c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_c_j = simulate_1(J, A, B, b, c_a, c_b, switch_c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_g = simulate_1(J, A, B, b, c_a, c_b, c_j, switch_g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_h = simulate_1(J, A, B, b, c_a, c_b, c_j, g, switch_h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_l = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, switch_l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_r = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, switch_r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_switch_s = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, r, switch_s, ba_a_avg, ba_a_upper, N_spec);
        
        
        /// Summed up contributions
        sum_ko_1_b += ba_ko_none - ba_ko_1_b;
        sum_ko_1_b_c_b += ba_ko_none - ba_ko_1_b_c_b;
        sum_ko_1_c_a += ba_ko_none - ba_ko_1_c_a;
        sum_ko_1_c_b += ba_ko_none - ba_ko_1_c_b;
        sum_ko_1_c_j += ba_ko_none - ba_ko_1_c_j;
        sum_ko_1_g += ba_ko_none - ba_ko_1_g;
        sum_ko_1_h += ba_ko_none - ba_ko_1_h;
        sum_ko_1_l += ba_ko_none - ba_ko_1_l;
        sum_ko_1_r += ba_ko_none - ba_ko_1_r;
        sum_ko_1_s += ba_ko_none - ba_ko_1_s;
        
        sum_ko_2_b += ba_ko_none - ba_ko_2_b;
        sum_ko_2_b_c_b += ba_ko_none - ba_ko_2_b_c_b;
        sum_ko_2_c_a += ba_ko_none - ba_ko_2_c_a;
        sum_ko_2_c_b += ba_ko_none - ba_ko_2_c_b;
        sum_ko_2_c_j += ba_ko_none - ba_ko_2_c_j;
        sum_ko_2_g += ba_ko_none - ba_ko_2_g;
        sum_ko_2_h += ba_ko_none - ba_ko_2_h;
        sum_ko_2_l += ba_ko_none - ba_ko_2_l;
        sum_ko_2_r += ba_ko_none - ba_ko_2_r;
        sum_ko_2_s += ba_ko_none - ba_ko_2_s;
        

        /// Summed up contributions of the difference
        sum_switch_b += ba_ko_none - ba_switch_b;
        sum_switch_b_c_b += ba_ko_none - ba_switch_b_c_b;
        sum_switch_c_a += ba_ko_none - ba_switch_c_a;
        sum_switch_c_b += ba_ko_none - ba_switch_c_b;
        sum_switch_c_j += ba_ko_none - ba_switch_c_j;
        sum_switch_g += ba_ko_none - ba_switch_g;
        sum_switch_h += ba_ko_none - ba_switch_h;
        sum_switch_l += ba_ko_none - ba_switch_l;
        sum_switch_r += ba_ko_none - ba_switch_r;
        sum_switch_s += ba_ko_none - ba_switch_s;
        
        
        /// Summed up proportional contributions
	      real ba_sum_rate = sum(ba_ko_none) / BA_sum;
        
	      sum_ko_1_prop_b     += (ba_ko_none ./ ba_ko_1_b) / ba_sum_rate;
	      sum_ko_1_prop_b_c_b += (ba_ko_none ./ ba_ko_1_b_c_b) / ba_sum_rate;
	      sum_ko_1_prop_c_a   += (ba_ko_none ./ ba_ko_1_c_a) / ba_sum_rate;
	      sum_ko_1_prop_c_b   += (ba_ko_none ./ ba_ko_1_c_b) / ba_sum_rate;
	      sum_ko_1_prop_c_j   += (ba_ko_none ./ ba_ko_1_c_j) / ba_sum_rate;
	      sum_ko_1_prop_g     += (ba_ko_none ./ ba_ko_1_g) / ba_sum_rate;
	      sum_ko_1_prop_h     += (ba_ko_none ./ ba_ko_1_h) / ba_sum_rate;
	      sum_ko_1_prop_l     += (ba_ko_none ./ ba_ko_1_l) / ba_sum_rate;
	      sum_ko_1_prop_r     += (ba_ko_none ./ ba_ko_1_r) / ba_sum_rate;
	      sum_ko_1_prop_s     += (ba_ko_none ./ ba_ko_1_s) / ba_sum_rate;
        
	      sum_ko_2_prop_b     += (ba_ko_none ./ ba_ko_2_b) / ba_sum_rate;
	      sum_ko_2_prop_b_c_b += (ba_ko_none ./ ba_ko_2_b_c_b) / ba_sum_rate;
	      sum_ko_2_prop_c_a   += (ba_ko_none ./ ba_ko_2_c_a) / ba_sum_rate;
	      sum_ko_2_prop_c_b   += (ba_ko_none ./ ba_ko_2_c_b) / ba_sum_rate;
	      sum_ko_2_prop_c_j   += (ba_ko_none ./ ba_ko_2_c_j) / ba_sum_rate;
	      sum_ko_2_prop_g     += (ba_ko_none ./ ba_ko_2_g) / ba_sum_rate;
	      sum_ko_2_prop_h     += (ba_ko_none ./ ba_ko_2_h) / ba_sum_rate;
	      sum_ko_2_prop_l     += (ba_ko_none ./ ba_ko_2_l) / ba_sum_rate;
	      sum_ko_2_prop_r     += (ba_ko_none ./ ba_ko_2_r) / ba_sum_rate;
	      sum_ko_2_prop_s     += (ba_ko_none ./ ba_ko_2_s) / ba_sum_rate;
          
        
      }  // end if i < fixiter_min
      
      
      /// !
      J = J_1;
      A = A_1;
      B = B_1;
      
      i += 1;

    } // end while i < fixiter_max
    
    // array with 3 (states) + 1 (BA) + 1 (eps) + 1 (n_iter) + some variables (overall N_fix)
    array[N_fix] vector[N_spec] fix = {J_1, A_1, B_1, BA_1,
                                       eps_ba, rep_vector(i, N_spec), // int i gets cast to real
                                       
                                       //// when considering the whole period, use i, when considering only the first period use fixiter_min as a denominator here
                                       // indices 7–26
                                       sum_ko_1_b/fixiter_min, sum_ko_1_b_c_b/fixiter_min, sum_ko_1_c_a/fixiter_min, sum_ko_1_c_b/fixiter_min, sum_ko_1_c_j/fixiter_min, sum_ko_1_g/fixiter_min, sum_ko_1_h/fixiter_min, sum_ko_1_l/fixiter_min, sum_ko_1_r/fixiter_min, sum_ko_1_s/fixiter_min,
                                       sum_ko_2_b/fixiter_min, sum_ko_2_b_c_b/fixiter_min, sum_ko_2_c_a/fixiter_min, sum_ko_2_c_b/fixiter_min, sum_ko_2_c_j/fixiter_min, sum_ko_2_g/fixiter_min, sum_ko_2_h/fixiter_min, sum_ko_2_l/fixiter_min, sum_ko_2_r/fixiter_min, sum_ko_2_s/fixiter_min,
                                       
                                       // indices 27–46
                                       sum_ko_1_prop_b/fixiter_min, sum_ko_1_prop_b_c_b/fixiter_min, sum_ko_1_prop_c_a/fixiter_min, sum_ko_1_prop_c_b/fixiter_min, sum_ko_1_prop_c_j/fixiter_min, sum_ko_1_prop_g/fixiter_min, sum_ko_1_prop_h/fixiter_min, sum_ko_1_prop_l/fixiter_min, sum_ko_1_prop_r/fixiter_min, sum_ko_1_prop_s/fixiter_min,
                                       sum_ko_2_prop_b/fixiter_min, sum_ko_2_prop_b_c_b/fixiter_min, sum_ko_2_prop_c_a/fixiter_min, sum_ko_2_prop_c_b/fixiter_min, sum_ko_2_prop_c_j/fixiter_min, sum_ko_2_prop_g/fixiter_min, sum_ko_2_prop_h/fixiter_min, sum_ko_2_prop_l/fixiter_min, sum_ko_2_prop_r/fixiter_min, sum_ko_2_prop_s/fixiter_min,
                                       
                                       // indices 47–56
                                       sum_switch_b/fixiter_min, sum_switch_b_c_b/fixiter_min, sum_switch_c_a/fixiter_min, sum_switch_c_b/fixiter_min, sum_switch_c_j/fixiter_min, sum_switch_g/fixiter_min, sum_switch_h/fixiter_min, sum_switch_l/fixiter_min, sum_switch_r/fixiter_min, sum_switch_s/fixiter_min};
                                    
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
  int<lower=0> N_beta;
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
  vector[2] prior_l_log;
  vector[2] prior_s_log;

  // Specific priors
  array[2] vector[N_species] prior_g_log;
  array[2] vector[N_species] prior_h_log;
  // array[2] vector[N_species] prior_l_log;
  array[2] vector[N_species] prior_r_log;

}


///////////////////////////////////////////////////////////////////////////////



transformed data {
  
  // Times are all assumed to be shifted to start at 1!
  
  //// Data for generated quantities
  int N_fix = 6; // an array of vectors[N_species] { J, A, B, BA, eps, n_iter}
  int N_fix_contributions = 56; // an array of vectors[N_species] { J, A, B, BA, eps, n_iter, x * 2 * (9+1) diff_ko_parameter }
  
}


///////////////////////////////////////////////////////////////////////////////



parameters {
  //// Model parameters  
  vector[N_species] b_log;
  vector[N_species] c_a_log;
  vector[N_species] c_b_log;
  vector[N_species] c_j_log;
  vector<upper=0>[N_species] g_log;
  vector<upper=0>[N_species] h_log;
  vector[N_species] l_log;
  vector[N_species] r_log;
  vector[N_species] s_log;
  

  //// Dispersion
  vector<lower=0>[N_protocolTax] phi_obs_inv; // error in neg_binomial per tax and stage
  
  //// Initial state
  array[N_locs] vector<lower=0, upper=1>[N_pops] state_init_raw;

}


///////////////////////////////////////////////////////////////////////////////



transformed parameters {
      
  //// Local variables
  array[N_locs] vector<lower=0>[N_species] L_loc;
  array[N_locs] vector<lower=0>[N_pops] state_init;

  for(loc in 1:N_locs) {
    
    L_loc[loc, ] = exp(l_log + L_smooth_log[loc, ]); /// l * L_smooth == exp(l_log + L_smooth_log)
    
    state_init[loc] = state_init_raw[loc] .* upper_init;
    
  }
  
  vector<lower=0>[L_y] y_hat = unpack(state_init, time_max, times,
                                b_log, c_a_log, c_b_log, c_j_log, g_log, h_log, L_loc, r_log, s_log,
                                ba_a_avg, ba_a_upper,
                                n_obs, n_yhat,
                                N_species, N_pops, L_y, N_locs, // fixed numbers
                                i_j, i_a, i_b);
                                
  vector[L_y] y_hat_offset = y_hat .* offset_data;
  
  vector<lower=0>[N_protocolTax] phi_obs = inv(phi_obs_inv);
  vector[L_y] phi_obs_rep = phi_obs[rep_protocolTax2y];
  
}


///////////////////////////////////////////////////////////////////////////////



model {

  //—————————————————————————————————————————————————————————————————————//
  // Priors       ------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  //// "Hyperprior"
  phi_obs_inv ~ normal(0, 10);

  for(l in 1:N_locs) { 

    //// Prior for initial state
    state_init[l] ~ gamma(alpha_init[l], beta_init[l]); // state_init was just a linearly transformed. -> No Jacobian correction necessary.

  }
  
  
  //// Priors for Parameters  
  b_log ~ normal(prior_b_log[1], prior_b_log[2]);
  c_a_log ~ normal(prior_c_a_log[1], prior_c_a_log[2]);
  c_b_log ~ normal(prior_c_b_log[1], prior_c_b_log[2]);
  c_j_log ~ normal(prior_c_j_log[1], prior_c_j_log[2]);
  g_log ~ normal(prior_g_log[1], prior_g_log[2]); // species-specific! 
  h_log ~ normal(prior_h_log[1], prior_h_log[2]); // species-specific! 
  l_log ~ normal(prior_l_log[1], prior_l_log[2]);
  r_log ~ normal(prior_r_log[1], prior_r_log[2]); // species-specific! 
  s_log ~ normal(prior_s_log[1], prior_s_log[2]);

  
  //—————————————————————————————————————————————————————————————————————//
  // Model       -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  y ~ neg_binomial_2(y_hat_offset, phi_obs_rep);

}


///////////////////////////////////////////////////////////////////////////////



generated quantities {
  
  //—————————————————————————————————————————————————————————————————————//
  // Print MC state       ----------------------------------------------//
  //———————————————————————————————————————————————————————————————————// 
  
  // print("b_log: ", b_log);
  // print("c_a_log: ", c_a_log);
  // print("c_b_log: ", c_b_log);
  // print("c_j_log: ", c_j_log);
  // print("g_log: ", g_log);
  // print("h_log: ", h_log);
  // print("l_log: ", l_log);
  // print("r_log: ", r_log);
  // print("s_log: ", s_log);


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
  for (s in 1:N_species) avg_L_loc[s] = mean(L_loc[, s]);


  //—————————————————————————————————————————————————————————————————————//
  // Priors                  -------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  ///// Priors -----------------------------------
  real b_log_prior = normal_rng(prior_b_log[1], prior_b_log[2]);
  real c_a_log_prior = normal_rng(prior_c_a_log[1], prior_c_a_log[2]);
  real c_b_log_prior = normal_rng(prior_c_b_log[1], prior_c_b_log[2]);
  real c_j_log_prior = normal_rng(prior_c_j_log[1], prior_c_j_log[2]);
  vector<upper=0>[N_species] g_log_prior = -sqrt(square(to_vector(normal_rng(prior_g_log[1,], prior_g_log[2,]))));
  vector<upper=0>[N_species] h_log_prior = -sqrt(square(to_vector(normal_rng(prior_h_log[1,], prior_h_log[2,]))));
  real l_log_prior = normal_rng(prior_l_log[1], prior_l_log[2]);
  // vector[N_species] l_log_prior = to_vector(normal_rng(prior_l_log[1,], prior_l_log[2,]));
  vector[N_species] r_log_prior = to_vector(normal_rng(prior_r_log[1,], prior_r_log[2,]));  
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
  
  
  //// Declarations of posterior quantites (as global variables).
  // … are directly initiated with zeroes or 9, so that there are never NaNs in generated quantities.
  
  array[N_locs, N_fix_contributions] vector[N_species] Fix = rep_array(rep_vector(0, N_species), N_locs, N_fix_contributions); // N_locs arrays of vectors[N_specices] { J, A, B, BA, eps, n_iter, 2 * 9 * diff_ko_parameter }
  
  array[N_locs] vector[N_species] J_init = rep_array(rep_vector(0.0, N_species), N_locs);
  array[N_locs] vector[N_species] A_init = J_init;
  array[N_locs] vector[N_species] B_init = J_init;
  array[N_locs] vector[N_species] J_fix = J_init;
  array[N_locs] vector[N_species] A_fix = J_init;
  array[N_locs] vector[N_species] B_fix = J_init;
  array[N_locs] vector[N_species] ba_init = J_init;
  array[N_locs] vector[N_species] ba_fix = J_init;
  
  array[N_locs] vector[N_species] eps_ba_fix = J_init;
  array[N_locs] real iterations_fix = rep_array(0.0, N_locs);

  array[N_locs] vector[N_species] sum_ko_1_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_b_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_g_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_h_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_l_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_r_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_s_fix = J_init;
  
  array[N_locs] vector[N_species] sum_ko_2_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_b_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_g_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_h_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_l_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_r_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_s_fix = J_init;
  
  array[N_locs] vector[N_species] sum_ko_1_prop_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_b_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_g_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_h_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_l_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_r_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_s_fix = J_init;
  
  array[N_locs] vector[N_species] sum_ko_2_prop_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_b_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_g_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_h_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_l_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_r_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_s_fix = J_init;
  
  
  array[N_locs] vector[N_species] sum_switch_b_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_b_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_g_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_h_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_l_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_r_fix = J_init;
  array[N_locs] vector[N_species] sum_switch_s_fix = J_init;

  int fixiter_max = 5000;
  int fixiter_min = 250;

  
  array[N_locs] int converged_fix = rep_array(9, N_locs); // tolerance has been reached

  array[N_locs] int dominant_init = converged_fix;
  array[N_locs] int dominant_fix = converged_fix;
  array[N_locs] int major_init = converged_fix;
  array[N_locs] int major_fix = converged_fix;
  
  // array[N_locs] int major_fix_ko_b = converged_fix;
  // array[N_locs] int major_fix_ko_s = converged_fix;
  // array[N_locs] int major_fix_ko_2_b = converged_fix;
  // array[N_locs] int major_fix_ko_2_s = converged_fix;
  
  // array[N_locs] int major_fix_switch_b = converged_fix;
  array[N_locs] int major_fix_switch_b_c_b = converged_fix;
  // array[N_locs] int major_fix_switch_b_c_a_c_b_h = converged_fix;
  array[N_locs] int major_fix_switch_c_a = converged_fix;
  // array[N_locs] int major_fix_switch_c_b = converged_fix;
  array[N_locs] int major_fix_switch_c_j = converged_fix;
  array[N_locs] int major_fix_switch_g = converged_fix;
  // array[N_locs] int major_fix_switch_g_l_r_s = converged_fix;
  array[N_locs] int major_fix_switch_h = converged_fix;
  // array[N_locs] int major_fix_switch_l = converged_fix;
  array[N_locs] int major_fix_switch_l_r = converged_fix;
  array[N_locs] int major_fix_switch_s = converged_fix;
  // 
  // //// Declarations of counterfactual posterior quantities
  // array[N_locs, N_fix] vector[N_species] Fix_ko_b = Fix;
  // array[N_locs, N_fix] vector[N_species] Fix_ko_s = Fix;
  // array[N_locs, N_fix] vector[N_species] Fix_ko_2_b = Fix;
  // array[N_locs, N_fix] vector[N_species] Fix_ko_2_s = Fix;
  // 
  // array[N_locs] vector[N_species] ba_fix_ko_b = J_init;
  // array[N_locs] vector[N_species] ba_fix_ko_s = J_init;
  // array[N_locs] vector[N_species] ba_fix_ko_2_b = J_init;
  // array[N_locs] vector[N_species] ba_fix_ko_2_s = J_init;
   
  // array[N_locs, N_fix] vector[N_species] Fix_switch_b = Fix;
  array[N_locs, N_fix] vector[N_species] Fix_switch_b_c_b = Fix;
  // array[N_locs, N_fix] vector[N_species] Fix_switch_b_c_a_c_b_h = Fix;
  array[N_locs, N_fix] vector[N_species] Fix_switch_c_a = Fix;
  // array[N_locs, N_fix] vector[N_species] Fix_switch_c_b = Fix;
  array[N_locs, N_fix] vector[N_species] Fix_switch_c_j = Fix;
  array[N_locs, N_fix] vector[N_species] Fix_switch_g = Fix;
  // array[N_locs, N_fix] vector[N_species] Fix_switch_g_l_r_s = Fix;
  array[N_locs, N_fix] vector[N_species] Fix_switch_h = Fix;
  // array[N_locs, N_fix] vector[N_species] Fix_switch_l = Fix;
  array[N_locs, N_fix] vector[N_species] Fix_switch_l_r = Fix;
  array[N_locs, N_fix] vector[N_species] Fix_switch_s = Fix;
   
  // array[N_locs] vector[N_species] ba_fix_switch_b = J_init;
  // array[N_locs] vector[N_species] ba_fix_switch_b_c_a_c_b_h = J_init;
  array[N_locs] vector[N_species] ba_fix_switch_b_c_b = J_init;
  array[N_locs] vector[N_species] ba_fix_switch_c_a = J_init;
  // array[N_locs] vector[N_species] ba_fix_switch_c_b = J_init;
  array[N_locs] vector[N_species] ba_fix_switch_c_j = J_init;
  array[N_locs] vector[N_species] ba_fix_switch_g = J_init;
  array[N_locs] vector[N_species] ba_fix_switch_h = J_init;
  // array[N_locs] vector[N_species] ba_fix_switch_g_l_r_s = J_init;
  // array[N_locs] vector[N_species] ba_fix_switch_l = J_init;
  array[N_locs] vector[N_species] ba_fix_switch_l_r = J_init;
  array[N_locs] vector[N_species] ba_fix_switch_s = J_init;


  //———————————————————————————————————————————————————————————————————//
  // Generate Posterior quantities conditioned on setting  -----------//
  //—————————————————————————————————————————————————————————————————//
  
  if (generateposteriorq) {

  
    //// Fix point iteration -------------------------------------------
    for(loc in 1:N_locs) {
      
      J_init[loc] = state_init[loc, 1:N_species];
      A_init[loc] = state_init[loc, (N_species+1):(N_species+N_species)];
      B_init[loc] = state_init[loc, (N_pops-N_species+1):N_pops];
      
      ba_init[loc] = state_init[loc, (N_pops-N_species+1):N_pops] + // State B
                     ba_a_avg .* state_init[loc, (N_species+1):(N_species+N_species)]; // State A * ba

      
      //// Booleans at init
      dominant_init[loc] = (ba_init[loc, 1]/ba_init[loc, 2]) > 3; // ba_1 > 75%
      major_init[loc] = ba_init[loc, 1] > ba_init[loc, 2]; // ba_1 > 50%
      

      //// Simulate fix point, given parameters
      Fix[loc] = iterateFix_contributions(state_init[loc],
                                          exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log),
                                          ba_a_avg, ba_a_upper,
                                          N_species, i_j, i_a, i_b,
                                          tolerance_fix, fixiter_max, fixiter_min, N_fix);
                                     
      iterations_fix[loc] = Fix[loc, 6, 1]; // the 6th element is the vector: [n_iter, n_iter]'
      converged_fix[loc] = iterations_fix[loc] < fixiter_max; // (i starts at 0), when fixiter_max is reached the model ran 5001 times
      
      if (converged_fix[loc]) { // && convergent[loc]
      
        //// unpack Fix
        J_fix[loc] = Fix[loc, 1];
        A_fix[loc] = Fix[loc, 2];
        B_fix[loc] = Fix[loc, 3];
        ba_fix[loc] = Fix[loc, 4];
        eps_ba_fix[loc] = Fix[loc, 5];        
        // Fix[loc, 6] is unpacked before
        
        //// Booleans at fixpoint
        dominant_fix[loc] = (ba_fix[loc, 1]/ba_fix[loc, 2]) > 3; // ba_1 > 75%
        major_fix[loc] = ba_fix[loc, 1] > ba_fix[loc, 2]; // ba_1 > 50%
        
        
        //// Unpack the cumulative contributions into variables
        sum_ko_1_b_fix[loc] = Fix[loc, 7];
        sum_ko_1_b_c_b_fix[loc] = Fix[loc, 8];
        sum_ko_1_c_a_fix[loc] = Fix[loc, 9];
        sum_ko_1_c_b_fix[loc] = Fix[loc, 10];
        sum_ko_1_c_j_fix[loc] = Fix[loc, 11];
        sum_ko_1_g_fix[loc] = Fix[loc, 12];
        sum_ko_1_h_fix[loc] = Fix[loc, 13];
        sum_ko_1_l_fix[loc] = Fix[loc, 14];
        sum_ko_1_r_fix[loc] = Fix[loc, 15];
        sum_ko_1_s_fix[loc] = Fix[loc, 16];
        
        sum_ko_2_b_fix[loc] = Fix[loc, 17];
        sum_ko_2_b_c_b_fix[loc] = Fix[loc, 18];
        sum_ko_2_c_a_fix[loc] = Fix[loc, 19];
        sum_ko_2_c_b_fix[loc] = Fix[loc, 20];
        sum_ko_2_c_j_fix[loc] = Fix[loc, 21];
        sum_ko_2_g_fix[loc] = Fix[loc, 22];
        sum_ko_2_h_fix[loc] = Fix[loc, 23];
        sum_ko_2_l_fix[loc] = Fix[loc, 24];
        sum_ko_2_r_fix[loc] = Fix[loc, 25];
        sum_ko_2_s_fix[loc] = Fix[loc, 26];
        
        sum_ko_1_prop_b_fix[loc] = Fix[loc, 27];
        sum_ko_1_prop_b_c_b_fix[loc] = Fix[loc, 28];
        sum_ko_1_prop_c_a_fix[loc] = Fix[loc, 29];
        sum_ko_1_prop_c_b_fix[loc] = Fix[loc, 30];
        sum_ko_1_prop_c_j_fix[loc] = Fix[loc, 31];
        sum_ko_1_prop_g_fix[loc] = Fix[loc, 32];
        sum_ko_1_prop_h_fix[loc] = Fix[loc, 33];
        sum_ko_1_prop_l_fix[loc] = Fix[loc, 34];
        sum_ko_1_prop_r_fix[loc] = Fix[loc, 35];
        sum_ko_1_prop_s_fix[loc] = Fix[loc, 36];

        sum_ko_2_prop_b_fix[loc] = Fix[loc, 37];
        sum_ko_2_prop_b_c_b_fix[loc] = Fix[loc, 38];
        sum_ko_2_prop_c_a_fix[loc] = Fix[loc, 39];
        sum_ko_2_prop_c_b_fix[loc] = Fix[loc, 40];
        sum_ko_2_prop_c_j_fix[loc] = Fix[loc, 41];
        sum_ko_2_prop_g_fix[loc] = Fix[loc, 42];
        sum_ko_2_prop_h_fix[loc] = Fix[loc, 43];
        sum_ko_2_prop_l_fix[loc] = Fix[loc, 44];
        sum_ko_2_prop_r_fix[loc] = Fix[loc, 45];
        sum_ko_2_prop_s_fix[loc] = Fix[loc, 46];
        
        
        sum_switch_b_fix[loc] = Fix[loc, 47];
        sum_switch_b_c_b_fix[loc] = Fix[loc, 48];
        sum_switch_c_a_fix[loc] = Fix[loc, 49];
        sum_switch_c_b_fix[loc] = Fix[loc, 50];
        sum_switch_c_j_fix[loc] = Fix[loc, 51];
        sum_switch_g_fix[loc] = Fix[loc, 52];
        sum_switch_h_fix[loc] = Fix[loc, 53];
        sum_switch_l_fix[loc] = Fix[loc, 54];
        sum_switch_r_fix[loc] = Fix[loc, 55];
        sum_switch_s_fix[loc] = Fix[loc, 56];
        
        
        //// Counterfactual fix point iteration
        vector[N_pops] state_fix = append_row(append_row(J_fix[loc],  A_fix[loc]),  B_fix[loc]);
        
        ////// ... with knocked out parameters
        // vector[N_species] ko = [0.0, 0.0]';
        // vector[N_species] ko_2_b = [exp(b_log[1]), 0]';
        // vector[N_species] ko_2_s = [exp(s_log[1]), 0]';
        // 
        // Fix_ko_b[loc] = iterateFix(state_init[loc], ko, exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        // Fix_ko_s[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), ko, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        // Fix_ko_2_b[loc] = iterateFix(state_init[loc], ko_2_b, exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        // Fix_ko_2_s[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), ko_2_s, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        // 
        // ba_fix_ko_b[loc] = Fix_ko_b[loc, 4];
        // ba_fix_ko_s[loc] = Fix_ko_s[loc, 4];
        // ba_fix_ko_2_b[loc] = Fix_ko_2_b[loc, 4];
        // ba_fix_ko_2_s[loc] = Fix_ko_2_s[loc, 4];


        // ////// ... with switched parameters
        vector[N_species] switch_b = exp([b_log[2], b_log[1]]'); // b_log[2:1] does not work at runtime
        vector[N_species] switch_c_a = exp([c_a_log[2], c_a_log[1]]');
        vector[N_species] switch_c_b = exp([c_b_log[2], c_b_log[1]]');
        vector[N_species] switch_c_j = exp([c_j_log[2], c_j_log[1]]');
        vector[N_species] switch_g = exp([g_log[2], g_log[1]]');
        vector[N_species] switch_h = exp([h_log[2], h_log[1]]');
        vector[N_species] switch_l = [L_loc[loc, 2], L_loc[loc, 1]]';
        vector[N_species] switch_r = exp([r_log[2], r_log[1]]');
        vector[N_species] switch_s = exp([s_log[2], s_log[1]]');
        
        // Fix_switch_b[loc] = iterateFix(state_init[loc], switch_b, exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        //ba_fix_switch_b[loc] = Fix_switch_b[loc, 4];
        
        Fix_switch_b_c_b[loc] = iterateFix(state_init[loc], switch_b, exp(c_a_log), switch_c_b, exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        ba_fix_switch_b_c_b[loc] = Fix_switch_b_c_b[loc, 4];
        
        // Fix_switch_b_c_a_c_b_h[loc] = iterateFix(state_init[loc], switch_b, switch_c_a, switch_c_b, exp(c_j_log), exp(g_log), switch_h, L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        // ba_fix_switch_b_c_a_c_b_h[loc] = Fix_switch_b_c_a_c_b_h[loc, 4];
        
        Fix_switch_c_a[loc] = iterateFix(state_init[loc], exp(b_log), switch_c_a, exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        ba_fix_switch_c_a[loc] = Fix_switch_c_a[loc, 4];
        
        // Fix_switch_c_b[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), switch_c_b, exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        // ba_fix_switch_c_b[loc] = Fix_switch_c_b[loc, 4];
        
        Fix_switch_c_j[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), switch_c_j , exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        ba_fix_switch_c_j[loc] = Fix_switch_c_j[loc, 4];
        
        Fix_switch_g[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), switch_g, exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        ba_fix_switch_g[loc] = Fix_switch_g[loc, 4];
        
        // Fix_switch_g_l_r_s[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), switch_g, exp(h_log), switch_l, switch_r, switch_s, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        // ba_fix_switch_g_l_r_s[loc] = Fix_switch_g_l_r_s[loc, 4];

        Fix_switch_h[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), switch_h, L_loc[loc, ], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        ba_fix_switch_h[loc] = Fix_switch_h[loc, 4];
        
        // Fix_switch_l[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, 2:1], exp(r_log), exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        // ba_fix_switch_l[loc] = Fix_switch_l[loc, 4];
        
        Fix_switch_l_r[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), switch_l, switch_r, exp(s_log), ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
        ba_fix_switch_l_r[loc] = Fix_switch_l_r[loc, 4];
        
        Fix_switch_s[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), switch_s, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);       
        ba_fix_switch_s[loc] = Fix_switch_s[loc, 4];

        // //// Counterfactual Booleans at fixpoint
        // major_fix_ko_b[loc] = ba_fix_ko_b[loc, 1] > ba_fix_ko_b[loc, 2]; // ba_1 > 50%
        // major_fix_ko_s[loc] = ba_fix_ko_s[loc, 1] > ba_fix_ko_s[loc, 2]; // ba_1 > 50%
        // major_fix_ko_2_b[loc] = ba_fix_ko_2_b[loc, 1] > ba_fix_ko_2_b[loc, 2]; // ba_1 > 50%
        // major_fix_ko_2_s[loc] = ba_fix_ko_2_s[loc, 1] > ba_fix_ko_2_s[loc, 2]; // ba_1 > 50%
         
        // major_fix_switch_b[loc] = ba_fix_switch_b[loc, 1] > ba_fix_switch_b[loc, 2]; // ba_1 > 50%
        major_fix_switch_b_c_b[loc] = ba_fix_switch_b_c_b[loc, 1] > ba_fix_switch_b_c_b[loc, 2]; // ba_1 > 50%
        // major_fix_switch_b_c_a_c_b_h[loc] = ba_fix_switch_b_c_a_c_b_h[loc, 1] > ba_fix_switch_b_c_a_c_b_h[loc, 2]; // ba_1 > 50%
        major_fix_switch_c_a[loc] = ba_fix_switch_c_a[loc, 1] > ba_fix_switch_c_a[loc, 2]; // ba_1 > 50%
        // major_fix_switch_c_b[loc] = ba_fix_switch_c_b[loc, 1] > ba_fix_switch_c_b[loc, 2]; // ba_1 > 50%
        major_fix_switch_c_j[loc] = ba_fix_switch_c_j[loc, 1] > ba_fix_switch_c_j[loc, 2]; // ba_1 > 50%
        major_fix_switch_g[loc] = ba_fix_switch_g[loc, 1] > ba_fix_switch_g[loc, 2]; // ba_1 > 50%
        // major_fix_switch_g_l_r_s[loc] = ba_fix_switch_g_l_r_s[loc, 1] > ba_fix_switch_g_l_r_s[loc, 2]; // ba_1 > 50%
        major_fix_switch_h[loc] = ba_fix_switch_h[loc, 1] > ba_fix_switch_h[loc, 2]; // ba_1 > 50%
        // major_fix_switch_l[loc] = ba_fix_switch_l[loc, 1] > ba_fix_switch_l[loc, 2]; // ba_1 > 50%
        major_fix_switch_l_r[loc] = ba_fix_switch_l_r[loc, 1] > ba_fix_switch_l_r[loc, 2]; // ba_1 > 50%
        major_fix_switch_s[loc] = ba_fix_switch_s[loc, 1] > ba_fix_switch_s[loc, 2]; // ba_1 > 50%

      }
  
    }

  } // end if(generateposteriorq)
}
