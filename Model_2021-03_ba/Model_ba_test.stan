functions {
  
  //// Difference equations
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
      // Note: log1p(expm1(a) + exp(b)) == log(exp(a) + exp(b)); It is important to expm1() some state (here J), because the rates are positive anyway
      // State[i_j, t]  =  ((r .* BA)/(1 + BA_sum) + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      State[i_j, t]  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      State[i_a, t]  =  (g .* J + (A - h .*A )) ./ (1 + c_a*BA_sum);
      State[i_b, t]  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
    
    }
    
    return State;
  }
  
  ////# Return states without altering them for state debugging
  // matrix simulate_null(vector initialstate, vector state_2, vector state_3, int N_pops) {
  //   
  //   // State matrix with [species, times]. Times columns is sensible in order to have col-major access in matrices and for to_vector() later on
  //   matrix[N_pops,3] State;
  //   State[,1] = initialstate;
  //   State[,2] = state_2;
  //   State[,3] = state_3;
  //   
  //   return State;
  // }
  

  //// Difference equations
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
    vector[N_spec] B_1  =  (1+b).*((h .* A_1 * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
    
    real BA_sum_1 = sum(A_1 .* ba_a_avg + B_1);
		
	// Second run: use limiting states
    J_1  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J_1) + s*BA_sum_1); // use the states already generated before: J_1, BA_sum_1
    A_1  =  (g .* J_1 + (A - h .*A )) ./ (1 + c_a*BA_sum_1);
    B_1  =  (1+b).*((h .* A_1 * ba_a_upper) + B) ./ (1 + c_b*BA_sum_1);
    
    vector[N_spec] ba_1 = (A_1 .* ba_a_avg) + B_1;
    
    return ba_1;
  }
  
  
  //// ODE integration and data assignment to one long vector is wrapped into a function here because
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
    
  
  //// Difference equations simulated up to the fix point given a maximum tolerance over all states.
  // Expects a state vector[N_pops]
  // returns a state vector of the form [J1, J2, A1, A2, B1, B2, BA1, BA2, eps_BA1, eps_BA2, iterations]
  array[] vector iterateFix(vector state_0,
                      vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                      vector ba_a_avg, real ba_a_upper,
                      int N_spec, array[] int i_j, array[] int i_a, array[] int i_b,
                      real tolerance_fix, int fixiter_max, int fixiter_min, int N_fix) {
                       
        
    /// Summed up contributions, initialize with 0 to avoid NaNs
    vector[N_spec] sum_ko_1_b = rep_vector(0.0, N_spec);
    vector[N_spec] sum_ko_1_c_a = sum_ko_1_b;
    vector[N_spec] sum_ko_1_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_1_c_j = sum_ko_1_b;
    vector[N_spec] sum_ko_1_g = sum_ko_1_b;
    vector[N_spec] sum_ko_1_h = sum_ko_1_b;
    vector[N_spec] sum_ko_1_l = sum_ko_1_b;
    vector[N_spec] sum_ko_1_r = sum_ko_1_b;
    vector[N_spec] sum_ko_1_s = sum_ko_1_b;
    
    vector[N_spec] sum_ko_2_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_c_a = sum_ko_1_b;
    vector[N_spec] sum_ko_2_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_c_j = sum_ko_1_b;
    vector[N_spec] sum_ko_2_g = sum_ko_1_b;
    vector[N_spec] sum_ko_2_h = sum_ko_1_b;
    vector[N_spec] sum_ko_2_l = sum_ko_1_b;
    vector[N_spec] sum_ko_2_r = sum_ko_1_b;
    vector[N_spec] sum_ko_2_s = sum_ko_1_b;
    
    vector[N_spec] sum_ko_1_prop_b = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_c_a = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_c_j = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_g = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_h = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_l = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_r = sum_ko_1_b;
    vector[N_spec] sum_ko_1_prop_s = sum_ko_1_b;
    
    vector[N_spec] sum_ko_2_prop_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_c_a = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_c_b = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_c_j = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_g = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_h = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_l = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_r = sum_ko_1_b;
    vector[N_spec] sum_ko_2_prop_s = sum_ko_1_b;

    
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
    // int notconvergent = 1;
    // vector[N_spec] KO = rep_vector(0.0, N_spec);
    
    while ( i < fixiter_min || (max(eps_ba) > tolerance_fix && i < fixiter_max) ) { // if notconvergent were a good criterion: (notconvergent && max(eps_ba) > tolerance_fix)
            
      vector[N_spec] BA = A .* ba_a_avg + B;
      real BA_sum = sum(BA);
      
      // s_1[i_j]  =  ((r .* BA)/(1 + BA_sum) + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      J_1  =  (r .* BA + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
      A_1  =  (g .* J + (A - h .*A )) ./ (1 + c_a*BA_sum);
      B_1  =  (1+b).*((h .* A * ba_a_upper) + B) ./ (1 + c_b*BA_sum);
      
      BA_1 = A_1 .* ba_a_avg + B_1; // New BA as additional state.

      // notconvergent = (1 <= norm(jacobian(s_1[i_j[1]], s_1[i_j[2]], s_1[i_a[1]], s_1[i_a[2]], s_1[i_b[1]], s_1[i_b[2]], b[1], b[2], c_b[1], c_b[2], c_j[1], c_j[2], g[1], g[2], h[1], h[2], l[1], l[2], r[1], r[2], s[1], s[2], ba_a_avg, ba_a_upper)) );
      eps_ba = fabs((BA_1 - BA) ./ BA_1);
      
      if (i < fixiter_min) {
      
        vector[N_spec] ba_ko_none = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);

        vector[N_spec] ba_ko_1_b = simulate_1(J, A, B, [0, b[2]]', c_a, c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_c_a  = simulate_1(J, A, B, b, [0, c_a[2]]', c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_c_b = simulate_1(J, A, B, b, c_a, [0, c_b[2]]', c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_c_j = simulate_1(J, A, B, b, c_a, c_b, [0, c_j[2]]', g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_g = simulate_1(J, A, B, b, c_a, c_b, c_j, [0, g[2]]', h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_h = simulate_1(J, A, B, b, c_a, c_b, c_j, g, [0, h[2]]', l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_l = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, [0, l[2]]', r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_r = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, [0, r[2]]', s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_1_s = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, r, [0, s[2]]', ba_a_avg, ba_a_upper, N_spec);
        
        vector[N_spec] ba_ko_2_b = simulate_1(J, A, B, [b[1], 0]', c_a, c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_c_a  = simulate_1(J, A, B, b, [c_a[1], 0]', c_b, c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_c_b = simulate_1(J, A, B, b, c_a, [c_b[1], 0]', c_j, g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_c_j = simulate_1(J, A, B, b, c_a, c_b, [c_j[1], 0]', g, h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_g = simulate_1(J, A, B, b, c_a, c_b, c_j, [g[1], 0]', h, l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_h = simulate_1(J, A, B, b, c_a, c_b, c_j, g, [h[1], 0]', l, r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_l = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, [l[1], 0]', r, s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_r = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, [r[1], 0]', s, ba_a_avg, ba_a_upper, N_spec);
        vector[N_spec] ba_ko_2_s = simulate_1(J, A, B, b, c_a, c_b, c_j, g, h, l, r, [s[1], 0]', ba_a_avg, ba_a_upper, N_spec);
        
        
        /// Summed up contributions
        sum_ko_1_b += ba_ko_none - ba_ko_1_b;
        sum_ko_1_c_a += ba_ko_none - ba_ko_1_c_a;
        sum_ko_1_c_b += ba_ko_none - ba_ko_1_c_b;
        sum_ko_1_c_j += ba_ko_none - ba_ko_1_c_j;
        sum_ko_1_g += ba_ko_none - ba_ko_1_g;
        sum_ko_1_h += ba_ko_none - ba_ko_1_h;
        sum_ko_1_l += ba_ko_none - ba_ko_1_l;
        sum_ko_1_r += ba_ko_none - ba_ko_1_r;
        sum_ko_1_s += ba_ko_none - ba_ko_1_s;
        
        sum_ko_2_b += ba_ko_none - ba_ko_2_b;
        sum_ko_2_c_a += ba_ko_none - ba_ko_2_c_a;
        sum_ko_2_c_b += ba_ko_none - ba_ko_2_c_b;
        sum_ko_2_c_j += ba_ko_none - ba_ko_2_c_j;
        sum_ko_2_g += ba_ko_none - ba_ko_2_g;
        sum_ko_2_h += ba_ko_none - ba_ko_2_h;
        sum_ko_2_l += ba_ko_none - ba_ko_2_l;
        sum_ko_2_r += ba_ko_none - ba_ko_2_r;
        sum_ko_2_s += ba_ko_none - ba_ko_2_s;
        
        /// Summed up proportional contributions
	    real ba_rate = sum(ba_ko_none)/BA_sum;
  
	    sum_ko_1_prop_b   += (ba_ko_none ./ ba_ko_1_b) / ba_rate; // 1. proportional ko decline rate (how much does the basal area grow in the non-ko model compared to the k.o. model), in proportion to 2. the total basal area proportional increment 
	    sum_ko_1_prop_c_a += (ba_ko_none ./ ba_ko_1_c_a) / ba_rate; // ... will later be divided by number of iterations
	    sum_ko_1_prop_c_b += (ba_ko_none ./ ba_ko_1_c_b) / ba_rate;
	    sum_ko_1_prop_c_j += (ba_ko_none ./ ba_ko_1_c_j) / ba_rate;
	    sum_ko_1_prop_g   += (ba_ko_none ./ ba_ko_1_g) / ba_rate;
	    sum_ko_1_prop_h   += (ba_ko_none ./ ba_ko_1_h) / ba_rate;
	    sum_ko_1_prop_l   += (ba_ko_none ./ ba_ko_1_l) / ba_rate;
	    sum_ko_1_prop_r   += (ba_ko_none ./ ba_ko_1_r) / ba_rate;
	    sum_ko_1_prop_s   += (ba_ko_none ./ ba_ko_1_s) / ba_rate;
  
	    sum_ko_2_prop_b   += (ba_ko_none ./ ba_ko_1_b) / ba_rate; // proportion of the ko growth rate to the overall ba growth rate, i.e. 
	    sum_ko_2_prop_c_a += (ba_ko_none ./ ba_ko_2_c_a) / ba_rate;
	    sum_ko_2_prop_c_b += (ba_ko_none ./ ba_ko_2_c_b) / ba_rate;
	    sum_ko_2_prop_c_j += (ba_ko_none ./ ba_ko_2_c_j) / ba_rate;
	    sum_ko_2_prop_g   += (ba_ko_none ./ ba_ko_2_g) / ba_rate;
	    sum_ko_2_prop_h   += (ba_ko_none ./ ba_ko_2_h) / ba_rate;
	    sum_ko_2_prop_l   += (ba_ko_none ./ ba_ko_2_l) / ba_rate;
	    sum_ko_2_prop_r   += (ba_ko_none ./ ba_ko_2_r) / ba_rate;
	    sum_ko_2_prop_s   += (ba_ko_none ./ ba_ko_2_s) / ba_rate;
      
      } // end if (i < fixiter_min)
      
      
      /// !
      J = J_1;
      A = A_1;
      B = B_1;
      
      i += 1;

    } // end while i < fixiter_max
    
    // array with 3 (states) + 1 (BA) + 1 (eps) + 1 (n_iter) + 9 (parameters) variables
    array[N_fix] vector[N_spec] fix = {J_1, A_1, B_1, BA_1,
                                       eps_ba, rep_vector(i, N_spec), // int i gets cast to real
                                       //// when considering the whole period, use i as a denominator here
                                       sum_ko_1_b/fixiter_min, sum_ko_1_c_a/fixiter_min, sum_ko_1_c_b/fixiter_min, sum_ko_1_c_j/fixiter_min, sum_ko_1_g/fixiter_min, sum_ko_1_h/fixiter_min, sum_ko_1_l/fixiter_min, sum_ko_1_r/fixiter_min, sum_ko_1_s/fixiter_min,
                                       sum_ko_2_b/fixiter_min, sum_ko_2_c_a/fixiter_min, sum_ko_2_c_b/fixiter_min, sum_ko_2_c_j/fixiter_min, sum_ko_2_g/fixiter_min, sum_ko_2_h/fixiter_min, sum_ko_2_l/fixiter_min, sum_ko_2_r/fixiter_min, sum_ko_2_s/fixiter_min,
                                       //
                                       sum_ko_1_prop_b/fixiter_min, sum_ko_1_prop_c_a/fixiter_min, sum_ko_1_prop_c_b/fixiter_min, sum_ko_1_prop_c_j/fixiter_min, sum_ko_1_prop_g/fixiter_min, sum_ko_1_prop_h/fixiter_min, sum_ko_1_prop_l/fixiter_min, sum_ko_1_prop_r/fixiter_min, sum_ko_1_prop_s/fixiter_min,
                                       sum_ko_2_prop_b/fixiter_min, sum_ko_2_prop_c_a/fixiter_min, sum_ko_2_prop_c_b/fixiter_min, sum_ko_2_prop_c_j/fixiter_min, sum_ko_2_prop_g/fixiter_min, sum_ko_2_prop_h/fixiter_min, sum_ko_2_prop_l/fixiter_min, sum_ko_2_prop_r/fixiter_min, sum_ko_2_prop_s/fixiter_min};
                                    
    return fix;
  }
  
  
  
  //// Transforms the vertex form into normal polynomial.
  array[] vector transformToNormal(array[] vector V) {
    
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
//  array[] vector transformToVertex(matrix P) {
//    
//    // P is a matrix[N_beta, N_species]
//    array[2] vector[5] V;
//    
//    for (i in 1:2) {
//      
//      real c = P[1, i];
//      real b_1 = P[2, i];
//      real a_1 = P[3, i];
//      real b_2 = P[4, i];
//      real a_2 = P[5, i];
//      
//      // assumes vertex form f(x,y) == a_1*(x−p)^2 + a_2*(y−q)^2 + z
//      // output vector[z, p, a_1, q, a_2]
//      V[i, 2] = -b_1 / (2*a_1); // p
//      V[i, 4] = -b_2 / (2*a_2); // q
//      V[i, 1] = c - a_2*V[i, 4]^2 - a_1*V[i, 2]^2; // z
//      V[i, 3] = a_1;
//      V[i, 5] = a_2;
//    }
//    return V;
//  }
//  
//  // Implementation of gamma probability density with zero hurdle model
//  real gamma_0_lpdf(vector y, vector y_hat_rep, vector alpha_rep, real theta, int L_y) {
//    
//    real t;
//    vector[L_y] beta_rep = alpha_rep ./ y_hat_rep;
//    for (l in 1:L_y) {
//	
//	// From a process point of view:
//	// rbinom(100, 1, 0.2) * rgamma(100, 3, 2)
//    
//      if (y[l] == 0) {
//        // Likelihood of 0 coming from probability theta; synonymous to t += bernoulli_lpmf(1 | theta);
//        t += log(theta);
//      }
//      else {
//        t += log1m(theta) + // synonymous to bernoulli_lpmf(0 | theta)
//             gamma_lpdf(y[l] | alpha_rep[l], beta_rep[l]);
//      }
//    }
//    return t;
//  }


// Implementation of negbinomial probability density with zero inflation
real neg_binomial_0_lpmf(int y, real y_hat, real phi_obs, real theta) {

 real t; // target

 if (y == 0) {
   // Joint Likelihood of 0 coming from probability theta or negbinonial
 	t = log_sum_exp(bernoulli_lpmf(1 | theta),
                      bernoulli_lpmf(0 | theta) + neg_binomial_2_lpmf(y | y_hat, phi_obs));
 } else {
// Joint Likelihood of 0 coming from probability theta_rep or negbinonial
 	t = bernoulli_lpmf(0 | theta) +  // log1m(theta) synonymous to bernoulli_lpmf(0 | theta_rep)?
 		neg_binomial_2_lpmf(y | y_hat, phi_obs);
 }
 return t; // target wich will get summed up at each run
}
 
 
//// Implementation of zi negbinomial random number generator
//array[] real neg_binomial_0_rng(vector y_hat_rep, vector phi_obs_rep, vector theta_rep, int L_y) {
//  
//  //
//}


//// Implementation of poisson probability density with zero inflation
real poisson_0_lpmf(int y, real y_hat, real theta) {
   
  real t; // target
  
  if (y == 0) {
    // Joint Likelihood of 0 coming from probability theta or poisson
  	t = log_sum_exp(bernoulli_lpmf(1 | theta),
                    bernoulli_lpmf(0 | theta) + poisson_lpmf(y | y_hat));
  } else {
	// Joint Likelihood
  	t = bernoulli_lpmf(0 | theta) +
  		poisson_lpmf(y | y_hat);
  }
  return t; // target wich will get summed up at each run
 }


//// Implementation of zi negbinomial random number generator
//  array[] real poisson_0_rng(vector y_hat_rep, vector theta_rep, int L_y) {
//	//
//  }


}

data {
  
  //// On assumed data structure
  // The most comprehensive data set is L_y with grouping.
  // Everything is subset from the master subsets with these groupings (*_obs, *_y0) and thus consistently sorted.
  // Ratio: "resurveys" includes "pops" due of the return structure in ode_*(); "pops" includes "plots" because of the loc == population assumption (randowm effect from location/resurvey/pop to .../plots).
  // (With population interactions and ODE integration, n_species, n_stages, and n_times have to be constant within the process model level, i.e. here location.)
  // Factor "pops" is structured stages/species.
  
  //// N — number of observations/groups; L - lengths of ragged vectors
  int<lower=0> L_times; // locations/obsid
  int<lower=0> L_y; // locations/resurveys/pops/locs
  // int<lower=0> L_plots; // locations/plots
  // int<lower=0> L_init;
  // int<lower=0> L_a2b; // locations/obsid-1/n_species

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
  array[N_species] int<lower=1> i_j; // e.g. 1:4
  array[N_species] int<lower=1> i_a; // e.g. 5:9, etc.
  array[N_species] int<lower=1> i_b;

  //// rep - repeat indices within groups for broadcasting to the subsequent hierarchical levels
  //+ array[L_y] int<lower=1> rep_yhat2y; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  array[L_y] int<lower=1> rep_pops2y; // factor (1:6)
  array[L_y] int<lower=1> rep_protocolTax2y;
  array[L_y] int<lower=1> rep_obsidPop2y; // factor (1:12)
  // array[L_y] int<lower=1> rep_protocol2y; // factor (1:5)

  // int<lower=1> rep_locs2plots[L_plots]; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  // int<lower=1> rep_yhat2a2b[L_a2b];
  // int<lower=1> rep_species2a2b[L_a2b];
  // int<lower=1> rep_pops2init[L_init];
  
  //// sigma for regularizing phi
  //: vector<lower=0>[N_pops] sigma_phi;
  // vector<lower=0>[N_protocolTax] sigma_phi;
  
  //// actual data
  array[N_locs] int time_max;
  array[L_times] int times; // locations/observations
  // vector<lower=1>[L_a2b] timediff;
  
  matrix[N_locs, N_beta] X; // design matrix
  array[N_locs] vector[N_species] L_smooth_log;
  
  vector<lower=0>[L_y] offset_data;
  
  //// The response.
  array[L_y] int y;
  //  array[L_a2b] int a2b;

  
  //// Settings
  real<upper=0.5> tolerance_fix;
  // real dbh_lower_a; // 100
  // real dbh_lower_b; // 200
  real ba_a_upper;
  vector[N_species] ba_a_avg;
  int<lower=0,upper=1> generateposteriorq;
  
  ////$ Potential methods for altering the timestep
  // real<lower=0> parfactor;
  // real<lower=0> timestep;
  
  //// Priors. The 2 reflect the two parameters mu and sigma
  // environmentally-dependent priors are species-agnostic on purpose
  // assumes vertex form f(x,y) == a_1*(x−p)^2 + a_2*(y−q)^2 + z
  // and parameters vector[z, p, q, a_1, a_2]
  
  vector<lower=0>[N_pops] upper_init;
  
  /// Gamma version
  array[N_locs] vector<lower=0>[N_pops] alpha_init;
  array[N_locs] vector<lower=0>[N_pops] beta_init;
  
  ///* Lognormal version
  //* array[N_locs] vector[N_pops] Prior_state_init_log;

  vector[2] prior_b_log;
  vector[2] prior_c_a_log;
  vector[2] prior_c_b_log;
  vector[2] prior_c_j_log;

  array[2] vector[N_species] prior_g_log;
  array[2] vector[N_species] prior_h_log;
    
  // array[2] vector[N_species] prior_k_log;
  // vector[N_species] prior_k_log;
  // array[2] vector[N_species] prior_l_log;
  vector[2] prior_l_log;
  array[2] vector[N_species] prior_r_log;
  // vector[N_species] prior_r_log;
  
  vector[2] prior_s_log;
  
  // array[2] vector[N_beta] prior_Vertex_c_j;
  // array[2] vector[N_beta] prior_Vertex_g;
  // array[2] vector[N_beta] prior_Vertex_r;
  // array[2] vector[N_beta] prior_Vertex_s;
  
  // array[2] vector[N_species] prior_b_log;
  // array[2] vector[N_species] prior_c_a_log;
  // array[2] vector[N_species] prior_c_b_log;
  
  // array[2] vector[N_species] prior_l_log;
  
  //// nu from student-t for weak priors
  // real nu_student;
  
  ////# Data for state debugging
  //@ array[N_locs] vector<lower=0>[N_pops] state_init_data;
  // array[N_locs] vector<lower=0>[N_pops] state_2;
  // array[N_locs] vector<lower=0>[N_pops] state_3;
  

}

transformed data {
  
  // Times are all assumed to be shifted to start at 1!

  // priors
  // array[2] vector[N_beta] prior_Beta_c_j = transformToNormal(prior_Vertex_c_j);
  // array[2] vector[N_beta] prior_Beta_g = transformToNormal(prior_Vertex_g);
  // array[2] vector[N_beta] prior_Beta_r = transformToNormal(prior_Vertex_r);
  // array[2] vector[N_beta] prior_Beta_s = transformToNormal(prior_Vertex_s);
  
  //// Data for separate fitting of the initial state
  // vector[N_pops] y0 [N_locs, N_plots] = y[ , , 1, ];
  
  //// Data for generated quantities
  int N_fix = 42; // an array of vectors[N_species] { J, A, B, BA, eps, n_iter, 2 * 2 * 9 diff_ko_parameter }

  ////$ Timestep alternation
  // real factor_log = log(parfactor/timestep);
  
}


parameters {
  //// Level 1, (global, species): species-specific rates, indepndent of environment.
  
  // … independent of environment
  vector[N_species] b_log;
  vector[N_species] c_a_log;
  vector[N_species] c_b_log;
  vector[N_species] c_j_log;
  vector<upper=0>[N_species] g_log;
  vector<upper=0>[N_species] h_log;
  vector[N_species] s_log;
  
  // vector[N_species] k_log;
  vector[N_species] l_log;
  vector[N_species] r_log;

  
  // … dependent on environment. Matrix for convenient matrix multiplication
  // matrix[N_beta, N_species] Beta_c_j;
  // matrix[N_beta, N_species] Beta_g; // (J-, A+) transition rate from J to A
  // matrix[N_beta, N_species] Beta_r; // // (J+) flow into the system, dependent on env
  // matrix[N_beta, N_species] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A
  
  //// Random intercepts for input k
  // matrix[N_locs, N_species] K_loc_log_raw; // array[N_locs] vector[N_species] K_loc_log_raw; // here real array is used for compatibility with to_vector
  // vector<lower=0>[N_species] sigma_k_loc;
  
  
  //// Errors
  vector<lower=0>[N_protocolTax] phi_obs_inv; // error in neg_binomial per tax and stage
  //: vector<lower=0>[N_pops] phi_obs_inv_sqrt; // error in neg_binomial per tax and stage
  // vector<lower=0>[N_obsidPop] phi_obs_inv_sqrt; // error in neg_binomial per tax and stage
  
    // vector<lower=0>[N_protocol] zeta; // zero-offset_data parameter
	// real<lower=0> kappa_inv; // error in beta for h_log
    // vector<lower=0>[3] alpha_obs_inv; // observation error in gamma
    // vector<lower=0>[2] sigma_obs; // observation error
    // vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  
  //.. vector<lower=0,upper=1>[N_pops] theta; // zi probability
  
  // matrix[N_pops, timespan_max] u[N_locs];
  
  array[N_locs] vector<lower=0, upper=1>[N_pops] state_init_raw; //@ 
  //* array[N_locs] vector[N_pops] state_init_log_raw; // version with non-central
  // vector<lower=0>[N_pops] sigma_state_init;
  
  ////#
  // array[N_locs] vector<lower=0, upper=1>[N_pops] state_2_raw;
  // array[N_locs] vector<lower=0, upper=1>[N_pops] state_3_raw;
}


transformed parameters {
  
  ////# Fixed parameters for debugging
  // vector[N_species] b_log = rep_vector(-3, N_species);
  // vector[N_species] c_a_log = rep_vector(-5, N_species);
  // vector[N_species] c_b_log = rep_vector(-5, N_species);
  // vector[N_species] c_j_log = rep_vector(-10, N_species);
  // vector<upper=0>[N_species] g_log = rep_vector(-3, N_species);
  // vector<upper=0>[N_species] h_log = rep_vector(-2, N_species);
  // vector[N_species] s_log = rep_vector(-3, N_species);
  // vector[N_species] l_log = rep_vector(8, N_species);
  // vector[N_species] r_log = rep_vector(8, N_species);

  //// Random intercepts for input k (assumes both K_loc and L_loc)
  // array[N_locs] vector<lower=0>[N_species] K_loc;
  // array[N_locs] vector<lower=0>[N_species] L_loc;
  
  //// L version: Local input
  array[N_locs] vector<lower=0>[N_species] L_loc;
  
  //@ array[N_locs] vector<lower=0>[N_pops] state_init = state_init_data;
  array[N_locs] vector<lower=0>[N_pops] state_init;
  //* array[N_locs] vector[N_pops] state_init_log;
  // vector[L_y] offset_zeta;
  
  ////#
  //# array[N_locs] vector<lower=0>[N_pops] state_2;
  //# array[N_locs] vector<lower=0>[N_pops] state_3;
  
  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  // matrix[N_locs, N_species] C_j_log = X * Beta_c_j;
  // matrix[N_locs, N_species] G_log = X * Beta_g;
  // matrix[N_locs, N_species] R_log = X * Beta_r;
  // matrix[N_locs, N_species] S_log = X * Beta_s;

  vector<lower=0>[N_protocolTax] phi_obs = inv(phi_obs_inv);
  //: vector<lower=0>[N_pops] phi_obs = inv_square(phi_obs_inv_sqrt); // inv_square == square_inv;
  // vector<lower=0>[N_protocolTax] phi_obs = inv_square(phi_obs_inv_sqrt);
    // vector<lower=0>[3] alpha_obs = inv(alpha_obs_inv);
  
  
  for(loc in 1:N_locs) {
  
  	//// Random intercept version
  	//// K_loc[loc, ] = exp(k_log + sigma_k_loc .* K_loc_log_raw[loc, ]');
    //// L_loc[loc, ] = K_loc[loc, ] + exp(l_log + L_smooth_log[loc, ]);
    	// k + l * L_smooth
    	// k == exp(k_log + normal(k_loc, sigma)) // intercept, with non-centered loc-level random intercepts
    	// l * L_smooth == exp(l_log + L_smooth_log)
    	
    
    //// L version
    L_loc[loc, ] = exp(l_log + L_smooth_log[loc, ]); /// l * L_smooth == exp(l_log + L_smooth_log)
    
    
    state_init[loc] = state_init_raw[loc] .* upper_init; //@
    ///* Lognormal version, with ~ exp(Normal()) 
    //* state_init_log[loc] = Prior_state_init_log[loc] + state_init_log_raw[loc] .* [3, 3, 3, 3, 3, 3]'; //* [2, 2, 2, 2, 1, 1]';
    
    ////#
    // state_2[loc] = state_2_raw[loc] .* upper_init;
    // state_3[loc] = state_3_raw[loc] .* upper_init;
  }
  
  //  vector[L_y] zeta_rep = zeta[rep_protocol2y];
  //  for(j in 1:L_y) {
  //  	if (offset_data[j] == 0) {
  //  		offset_zeta[j] = zeta_rep[j];
  //  	} else {
  //  		offset_zeta[j] = offset_data[j];
  //  	}
  //  }
  
  
  vector<lower=0>[L_y] y_hat = unpack(state_init, time_max, times, //* unpack(state_init_log, time_max, times,
                                // b_log, c_a_log, c_b_log, c_j_log, g_log, h_log, l_log, r_log, s_log, // rates matrix[N_locs, N_species]; will have to be transformed
                                b_log, c_a_log, c_b_log, c_j_log, g_log, h_log, L_loc, r_log, s_log, // rates matrix[N_locs, N_species]; will have to be transformed
                                // b_log, c_a_log, c_b_log, C_j_log, G_log, h_log, L_loc, R_log, S_log, // rates matrix[N_locs, N_species]; will have to be transformed
                                ba_a_avg, ba_a_upper,
                                n_obs, n_yhat, // varying numbers per loc
                                N_species, N_pops, L_y, N_locs, // fixed numbers
                                i_j, i_a, i_b);
                                
  vector[L_y] y_hat_offset = y_hat .* offset_data; // offset_zeta
  
  vector[L_y] phi_obs_rep = phi_obs[rep_protocolTax2y];
  //: vector[L_y] phi_obs_rep = phi_obs[rep_pops2y];
  //.. vector[L_y] theta_rep = theta[rep_pops2y];
  
}


model {

  //—————————————————————————————————————————————————————————————————————//
  // Priors       -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  //// Hyperpriors

  phi_obs_inv ~ normal(0, 10);
  //: phi_obs_inv_sqrt ~ normal(rep_vector(0.0, N_pops), sigma_phi);
  //. phi_obs_inv ~ normal(rep_vector(0.0, N_pops), sigma_phi);
  // phi_obs_inv_sqrt ~ normal(rep_vector(0.0, N_obsidPop), sigma_phi);
  	// On prior choice for the overdispersion in negative binomial 2: https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#story-when-the-generic-prior-fails-the-case-of-the-negative-binomial
  
  //.. theta ~ beta(0.5, 1);
  
  // Random intercepts for input k
  // to_vector(K_loc_log_raw) ~ std_normal(); // Random intercept for l
  // sigma_k_loc ~ std_normal(); // Regularizing half-cauchy on sigma for random slope for l  ## cauchy(0, 2);
  
  // sigma_state_init ~ std_normal();

  // sigma_process ~ normal(0, 0.01);
  // sigma_obs ~ normal(0, [0.5, 0.1]); // for observations from predictions
  // alpha_obs_inv ~ normal(0, 0.1); // Observation error for gamma
  // zeta ~ normal(0, 1);
  

  //// Prior for initial state
  for(l in 1:N_locs) { 
    state_init[l] ~ gamma(alpha_init[l], beta_init[l]); // state_init is just a linear transform. -> No Jacobian correction necessary.
  }
  
  
  //// Priors for Parameters
  
  // prior_*[2, N_species]
  
  // b_log   ~ normal(prior_b_log[1,], prior_b_log[2,]);
  // c_a_log ~ normal(prior_c_a_log[1,], prior_c_a_log[2,]);
  // c_b_log ~ normal(prior_c_b_log[1,], prior_c_b_log[2,]);
  
  b_log ~ normal(prior_b_log[1], prior_b_log[2]); // b_log ~ student_t(nu_student, prior_b_log[1], prior_b_log[2]); // 
  c_a_log ~ normal(prior_c_a_log[1], prior_c_a_log[2]); // c_a_log ~ student_t(nu_student, prior_c_a_log[1], prior_c_a_log[2]);
  c_b_log ~ normal(prior_c_b_log[1], prior_c_b_log[2]); // c_b_log ~ student_t(nu_student, prior_c_b_log[1], prior_c_b_log[2]);
  c_j_log ~ normal(prior_c_j_log[1], prior_c_j_log[2]); // c_j_log ~ student_t(nu_student, prior_c_j_log[1], prior_c_j_log[2]);
  
  g_log ~ normal(prior_g_log[1,], prior_g_log[2,]);
  h_log ~ normal(prior_h_log[1,], prior_h_log[2,]);

  // k_log ~ normal(prior_k_log[1], prior_k_log[2]);
  l_log ~ normal(prior_l_log[1], prior_l_log[2]);
  r_log ~ normal(prior_r_log[1], prior_r_log[2]);
  
  s_log ~ normal(prior_s_log[1], prior_s_log[2]); // s_log ~ student_t(nu_student, prior_s_log[1], prior_s_log[2]);

  
  // same priors for both species
  // for (spec in 1:N_species) {
  //   // prior_Beta_*[2, N_beta]
  //   Beta_c_j[, spec] ~ normal(prior_Beta_c_j[1, ], prior_Beta_c_j[2, ]);
  //   Beta_g[, spec]  ~ normal(prior_Beta_g[1, ], prior_Beta_g[2, ]);
  //   Beta_r[, spec] ~ normal(prior_Beta_r[1, ], prior_Beta_r[2, ]);
  //   Beta_s[, spec] ~ normal(prior_Beta_s[1, ], prior_Beta_s[2, ]);
  // }
  
  //—————————————————————————————————————————————————————————————————————//
  // Model       -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    

  // a2b ~ poisson(y_hat[rep_yhat2a2b] .* exp(h_log)[rep_species2a2b] .* timediff);
  
  y ~ neg_binomial_2(y_hat_offset, phi_obs_rep);
  // y ~ poisson(y_hat_offset);

  
  // y ~ gamma_0(y_hat[rep_yhat2y], alpha_obs[rep_pops2y], theta, L_y);
  
  // for(l in 1:L_y) {
  //   //..
  //   y[l] ~ neg_binomial_0(y_hat_offset[l], phi_obs_rep[l], theta_rep[l]);
  //   //& y[l] ~ poisson_0(y_hat_offset[l], theta_rep[l]);
  // }
}


generated quantities {

  //—————————————————————————————————————————————————————————————————————//
  // Prediction  -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    

  array[L_y] real<lower=0> y_sim;
  
  //// y_hat_rep generated in transformed parameters
  
  y_sim = neg_binomial_2_rng(y_hat_offset, phi_obs_rep);
  //.. y_sim = neg_binomial_0_rng(y_hat_rep, phi_obs_rep, theta_rep, L_y);
  //& y_sim = poisson_0_rng(y_hat_offset, theta_rep, L_y);
  // y_sim = poisson_rng(y_hat_offset);


  //—————————————————————————————————————————————————————————————————————//
  // Averages over locations  ------------------------------------------//
  //———————————————————————————————————————————————————————————————————//
  vector[N_pops] avg_state_init;
  vector[N_species] avg_L_loc;
  
  for (p in 1:N_pops) avg_state_init[p] = mean(state_init[, p]);
  //* for (p in 1:N_pops) avg_state_init[p] = mean(exp(state_init_log[, p]));
  for (s in 1:N_species) avg_L_loc[s] = mean(L_loc[, s]);


  //—————————————————————————————————————————————————————————————————————//
  // Prior predictive checks -------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  ///// Priors -----------------------------------
  real b_log_prior = normal_rng(prior_b_log[1], prior_b_log[2]);
  real c_a_log_prior = normal_rng(prior_c_a_log[1], prior_c_a_log[2]);
  real c_b_log_prior = normal_rng(prior_c_b_log[1], prior_c_b_log[2]);
  real c_j_log_prior = normal_rng(prior_c_j_log[1], prior_c_j_log[2]);
  
  vector<upper=0>[N_species] g_log_prior = -sqrt(square(to_vector(normal_rng(prior_g_log[1,], prior_g_log[2,]))));
  vector<upper=0>[N_species] h_log_prior = -sqrt(square(to_vector(normal_rng(prior_h_log[1,], prior_h_log[2,]))));
  
  // vector[N_species] k_log_prior = to_vector(normal_rng(prior_k_log[1,], prior_k_log[2,]));
  // vector[N_species] l_log_prior = to_vector(normal_rng(prior_l_log[1,], prior_l_log[2,]));
  real l_log_prior = normal_rng(prior_l_log[1], prior_l_log[2]);
  vector[N_species] r_log_prior = to_vector(normal_rng(prior_r_log[1,], prior_r_log[2,]));
  // real r_log_prior = normal_rng(prior_r_log[1], prior_r_log[2]);

  
  
  real s_log_prior = normal_rng(prior_s_log[1], prior_s_log[2]); // real s_log_prior = student_t_rng(nu_student, prior_s_log[1], prior_s_log[2]); // 
  
  vector[N_species] vector_b_log_prior = to_vector(normal_rng(rep_vector(prior_b_log[1], N_species), rep_vector(prior_b_log[2], N_species)));
  vector[N_species] vector_c_a_log_prior = to_vector(normal_rng(rep_vector(prior_c_a_log[1], N_species), rep_vector(prior_c_a_log[2], N_species)));
  vector[N_species] vector_c_b_log_prior = to_vector(normal_rng(rep_vector(prior_c_b_log[1], N_species), rep_vector(prior_c_b_log[2], N_species)));
  vector[N_species] vector_c_j_log_prior = to_vector(normal_rng(rep_vector(prior_c_j_log[1], N_species), rep_vector(prior_c_j_log[2], N_species)));
  vector[N_species] vector_l_log_prior = to_vector(normal_rng(rep_vector(prior_l_log[1], N_species), rep_vector(prior_l_log[2], N_species)));
  // vector[N_species] vector_r_log_prior = to_vector(normal_rng(rep_vector(prior_r_log[1], N_species), rep_vector(prior_r_log[2], N_species)));
  vector[N_species] vector_s_log_prior = to_vector(normal_rng(rep_vector(prior_s_log[1], N_species), rep_vector(prior_s_log[2], N_species)));
  
  
  array[N_protocolTax] real<lower=0> phi_obs_prior = inv(sqrt(square(normal_rng(rep_vector(0.0, N_protocolTax), 5)))); //! generation of distribution probably wrong (not consistent with density transformations)???
  //: array[N_pops] real<lower=0> phi_obs_prior = inv_square(normal_rng(rep_vector(0.0, N_pops), sigma_phi)); //! generation of distribution probably wrong (not consistent with density transformations)???
  //. array[N_pops] real<lower=0> phi_obs_prior = inv(sqrt(square(normal_rng(rep_vector(0.0, N_pops), sigma_phi)))); //! generation of distribution probably wrong (not consistent with density transformations)???
  // array[N_obsidPop] real<lower=0> phi_obs_prior = inv(sqrt(square(normal_rng(rep_vector(0.0, N_obsidPop), sigma_phi)))); //! generation of distribution probably wrong (not consistent with density transformations)???
  
  //// Random intercept version for input k
  // array[N_locs, N_species] real K_loc_log_raw_prior;
  // array[N_locs] vector<lower=0>[N_species] L_loc_prior;  

  //// L version
  // array[N_locs] vector<lower=0>[N_species] L_loc_prior;
  
  // vector<lower=0>[N_species] sigma_k_loc_prior = sqrt(square(to_vector(normal_rng(rep_vector(0, N_species), rep_vector(1, N_species)))));
  // vector<lower=0>[N_protocol] zeta_prior = sqrt(square(to_vector(normal_rng(rep_vector(0, N_protocol), rep_vector(0.2, N_protocol)))));
  
  // for(loc in 1:N_locs) {
 	
    //// Random intercept version
    // K_loc_log_raw_prior[loc,] = normal_rng(rep_vector(0, N_species), rep_vector(1, N_species));
    // L_loc_prior[loc, ] = exp(k_log_prior + sigma_k_loc_prior .* to_vector(K_loc_log_raw_prior[loc, ])) + exp(l_log_prior + L_smooth_log[loc, ]);

    //// L version
    // L_loc_prior[loc, ] = exp(l_log_prior + L_smooth_log[loc, ]);

  // }
  
  
//  //// Prior simulation
//  vector<lower=0>[L_y] y_hat_prior;
//  vector<lower=0>[L_y] y_hat_prior_rep_offset;
//  // vector<lower=0>[L_y] offset_zeta_prior;
//  array[L_y] real<lower=0> y_prior_sim;
//                   
//  //  vector<lower=0>[L_y] zeta_prior_rep = zeta_prior[rep_protocol2y];
//  
//  //  for(j in 1:L_y) {
//  //  	if (offset_data[j] == 0) {
//  //  		offset_zeta_prior[j] = zeta_prior_rep[j];
//  //  	} else {
//  //  		offset_zeta_prior[j] = offset_data[j];
//  //  	}
//  //  }
// 
//  y_hat_prior = unpack(state_init, time_max, times, //* unpack(state_init_log, time_max, times,
//                       vector_b_log_prior, vector_c_a_log_prior, vector_c_b_log_prior, vector_c_j_log_prior, g_log_prior, h_log_prior, L_loc_prior, r_log_prior, vector_s_log_prior, // rates matrix[N_locs, N_species]; will have to be transformed
//                       ba_a_avg, ba_a_upper,
//                       n_obs, n_yhat, // varying numbers per loc
//                       N_species, N_pops, L_y, N_locs, // fixed numbers
//                       i_j, i_a, i_b);
// 
//  y_hat_prior_offset = y_hat_prior .* offset_data; // offset_zeta_prior
//
//
//  y_prior_sim = neg_binomial_2_rng(y_hat_prior_offset, phi_obs_prior[rep_pops2y]);
  
  
  //—————————————————————————————————————————————————————————————————————————//
  // pgq -------------------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————————//
  
  
  //// Rate tests -------------------------------------
  int greater_b = b_log[1] > b_log[2];
  int greater_c_a = c_a_log[1] > c_a_log[2];
  int greater_c_b = c_b_log[1] > c_b_log[2];
  int greater_c_j = c_j_log[1] > c_j_log[2];
  int greater_g = g_log[1] > g_log[2];
  int greater_h = h_log[1] > h_log[2];
  int greater_l = l_log[1] > l_log[2];
  // int greater_k = k_log[1] > k_log[2];
  int greater_r = r_log[1] > r_log[2];
  int greater_s = s_log[1] > s_log[2];
  
  
  //// Declarations of posterior quantites (as global variables).
  // … are directly initiated with zeroes or 9, so that there are never NaNs in generated quantities.
  
  array[N_locs, N_fix] vector[N_species] Fix = rep_array(rep_vector(0, N_species), N_locs, N_fix); // N_locs arrays of vectors[N_specices] { J, A, B, BA, eps, n_iter, 2 * 9 * diff_ko_parameter }
  
  // array[N_locs] vector[N_pops] state_fix = rep_array(rep_vector(0.0, N_pops), N_locs); // state_fix is a vector [J1, …, A1, …, B1, …, BA1, …]
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
  array[N_locs] vector[N_species] sum_ko_1_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_g_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_h_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_l_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_r_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_s_fix = J_init;
  
  array[N_locs] vector[N_species] sum_ko_2_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_g_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_h_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_l_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_r_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_s_fix = J_init;
  
  array[N_locs] vector[N_species] sum_ko_1_prop_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_g_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_h_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_l_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_r_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_1_prop_s_fix = J_init;
  
  array[N_locs] vector[N_species] sum_ko_2_prop_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_c_a_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_c_b_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_c_j_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_g_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_h_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_l_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_r_fix = J_init;
  array[N_locs] vector[N_species] sum_ko_2_prop_s_fix = J_init;

  
  int fixiter_max = 5000;
  int fixiter_min = 500;

  
  array[N_locs] int converged_fix = rep_array(9, N_locs); // tolerance has been reached

  array[N_locs] int dominant_init = converged_fix;
  array[N_locs] int dominant_fix = converged_fix;
  array[N_locs] int major_init = converged_fix;
  array[N_locs] int major_fix = converged_fix;
  array[N_locs] int major_fix_ko_s = converged_fix;
  array[N_locs] int major_fix_switch_s = converged_fix;

  
  //// Declarations of counterfactual posterior quantities
  array[N_locs, N_fix] vector[N_species] Fix_ko_s = Fix;
  array[N_locs] vector[N_species] ba_fix_ko_s = J_init;
  
  array[N_locs, N_fix] vector[N_species] Fix_switch_s = Fix;
  array[N_locs] vector[N_species] ba_fix_switch_s = J_init;
  
  
  //// Declarations of quantities for sensitivity checks (as global variables).
  // real log_prior = 0; // this is zero to prevent NaNs from being in the sum.
  // vector[L_y] log_lik = rep_vector(0, L_y);
  

  //// The conditional generation -------------------------------------
  if (generateposteriorq) {

    //———————————————————————————————————————————————————————————————————//
    // Posterior quantities --------------------------------------------//
    //—————————————————————————————————————————————————————————————————//
    
    
    //// Fix point iteration -------------------------------------------
    for(loc in 1:N_locs) {
      
      J_init[loc] = state_init[loc, 1:N_species];
      A_init[loc] = state_init[loc, (N_species+1):(N_species+N_species)];
      B_init[loc] = state_init[loc, (N_pops-N_species+1):N_pops];
      
      ba_init[loc] = state_init[loc, (N_pops-N_species+1):N_pops] + // State B
                     ba_a_avg .* state_init[loc, (N_species+1):(N_species+N_species)]; // State A * ba
                     
      //* ba_init[loc] = exp(state_init_log[loc, (N_pops-N_species+1):N_pops]) + // State B
      //*               ba_a_avg .* exp(state_init_log[loc, (N_species+1):(N_species+N_species)]); // State A * ba
      
      /// Booleans at init
      dominant_init[loc] = (ba_init[loc, 1]/ba_init[loc, 2]) > 3; // ba_1 > 75%
      major_init[loc] = ba_init[loc, 1] > ba_init[loc, 2]; // ba_1 > 50%
      

      //// Simulate fix point, given parameters
      Fix[loc] = iterateFix(state_init[loc], //* Fix[loc] = iterateFix(exp(state_init_log[loc]),
                            // exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), exp(l_log), exp(r_log), exp(s_log),
                            exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log),
                            // exp(b_log), exp(c_a_log), exp(c_b_log), exp(C_j_log[loc,]'), exp(G_log[loc,]'), exp(h_log), exp(L_loc[loc, ]), exp(R_log[loc,]'), exp(S_log[loc,]'),
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
        
        sum_ko_1_b_fix[loc] = Fix[loc, 7];
        sum_ko_1_c_a_fix[loc] = Fix[loc, 8];
        sum_ko_1_c_b_fix[loc] = Fix[loc, 9];
        sum_ko_1_c_j_fix[loc] = Fix[loc, 10];
        sum_ko_1_g_fix[loc] = Fix[loc, 11];
        sum_ko_1_h_fix[loc] = Fix[loc, 12];
        sum_ko_1_l_fix[loc] = Fix[loc, 13]; // k
        sum_ko_1_r_fix[loc] = Fix[loc, 14];
        sum_ko_1_s_fix[loc] = Fix[loc, 15];
        
        sum_ko_2_b_fix[loc] = Fix[loc, 16];
        sum_ko_2_c_a_fix[loc] = Fix[loc, 17];
        sum_ko_2_c_b_fix[loc] = Fix[loc, 18];
        sum_ko_2_c_j_fix[loc] = Fix[loc, 19];
        sum_ko_2_g_fix[loc] = Fix[loc, 20];
        sum_ko_2_h_fix[loc] = Fix[loc, 21];
        sum_ko_2_l_fix[loc] = Fix[loc, 22];
        sum_ko_2_r_fix[loc] = Fix[loc, 23];
        sum_ko_2_s_fix[loc] = Fix[loc, 24];
        
        sum_ko_1_prop_b_fix[loc] = Fix[loc, 25];
		sum_ko_1_prop_c_a_fix[loc] = Fix[loc, 26];
		sum_ko_1_prop_c_b_fix[loc] = Fix[loc, 27];
		sum_ko_1_prop_c_j_fix[loc] = Fix[loc, 28];
		sum_ko_1_prop_g_fix[loc] = Fix[loc, 29];
		sum_ko_1_prop_h_fix[loc] = Fix[loc, 30];
		sum_ko_1_prop_l_fix[loc] = Fix[loc, 31];
		sum_ko_1_prop_r_fix[loc] = Fix[loc, 32];
		sum_ko_1_prop_s_fix[loc] = Fix[loc, 33];
		
		sum_ko_2_prop_b_fix[loc] = Fix[loc, 34];
		sum_ko_2_prop_c_a_fix[loc] = Fix[loc, 35];
		sum_ko_2_prop_c_b_fix[loc] = Fix[loc, 36];
		sum_ko_2_prop_c_j_fix[loc] = Fix[loc, 37];
		sum_ko_2_prop_g_fix[loc] = Fix[loc, 38];
		sum_ko_2_prop_h_fix[loc] = Fix[loc, 39];
		sum_ko_2_prop_l_fix[loc] = Fix[loc, 40];
		sum_ko_2_prop_r_fix[loc] = Fix[loc, 41];
		sum_ko_2_prop_s_fix[loc] = Fix[loc, 42];
        

        //// Counterfactual fix point iteration ---------------------------
        vector[N_pops] state_fix = append_row(append_row(J_fix[loc],  A_fix[loc]),  B_fix[loc]);
        
        vector[N_species] ko_s_2 = [exp(s_log[1]), 0]';
		// Fix_ko_b[loc]
		// Fix_ko_c_a[loc]
		// Fix_ko_c_j[loc]
		// Fix_ko_k[loc]
		// Fix_ko_l[loc]
		// Fix_ko_r[loc]
		
		Fix_ko_s[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), ko_s_2, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
		//* Fix_ko_s[loc] = iterateFix(exp(state_init_log[loc]), exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), ko_s_2, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);
		
		ba_fix_ko_s[loc] = Fix_ko_s[loc, 4];
		
		
		vector[N_species] switch_s = exp(s_log[2:1]);		
		Fix_switch_s[loc] = iterateFix(state_init[loc], exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), switch_s, ba_a_avg, ba_a_upper, N_species, i_j, i_a, i_b, tolerance_fix, fixiter_max, fixiter_min, N_fix);		
		ba_fix_switch_s[loc] = Fix_switch_s[loc, 4];
		
		
		//// Booleans at fixpoint
        dominant_fix[loc] = (ba_fix[loc, 1]/ba_fix[loc, 2]) > 3; // ba_1 > 75%
        major_fix[loc] = ba_fix[loc, 1] > ba_fix[loc, 2]; // ba_1 > 50%
        major_fix_ko_s[loc] = ba_fix_ko_s[loc, 1] > ba_fix_ko_s[loc, 2]; // ba_1 > 50%
        major_fix_switch_s[loc] = ba_fix_switch_s[loc, 1] > ba_fix_switch_s[loc, 2]; // ba_1 > 50%
      
      }
  
    }
    

    //———————————————————————————————————————————————————————————————————//
    // Sensitivity analysis --------------------------------------------//
    //—————————————————————————————————————————————————————————————————//
  
    // for(loc in 1:N_locs) {
    //   log_prior += gamma_lpdf(alpha_init[loc], beta_init[loc]); //* log_prior += normal_lpdf(state_init_log[loc,] | Prior_state_init_log[loc,], [2, 2, 2, 2, 1, 1]);
    // }
    // 
    // // joint prior specification, sum of all logpriors
    // log_prior = log_prior +
    // 		  normal_lpdf(phi_obs_inv_sqrt | rep_vector(0.0, N_pops), sigma_phi) +
	//   		  //// Random intecept K version version
	//            // normal_lpdf(sigma_k_loc | 0, 1) +		  
	//   		  // normal_lpdf(to_vector(K_loc_log_raw) | 0, 1) +
	//   		  
    //            normal_lpdf(b_log | prior_b_log[1], prior_b_log[2]) + // student_t_lpdf(b_log | nu_student, prior_b_log[1], prior_b_log[2]) +
	//   		  normal_lpdf(c_a_log | prior_c_a_log[1], prior_c_a_log[2]) + // student_t_lpdf(c_a_log | nu_student, prior_c_a_log[1], prior_c_a_log[2]) +
	//   		  normal_lpdf(c_b_log | prior_c_b_log[1], prior_c_b_log[2]) + // student_t_lpdf(c_b_log | nu_student, prior_c_b_log[1], prior_c_b_log[2]) +
	//   		  normal_lpdf(c_j_log | prior_c_j_log[1], prior_c_j_log[2]) + // student_t_lpdf(c_j_log | nu_student, prior_c_j_log[1], prior_c_j_log[2]) +
	//   		  normal_lpdf(g_log | prior_g_log[1,], prior_g_log[2,]) +
	//   		  normal_lpdf(h_log | prior_h_log[1,], prior_h_log[2,]) +
	//   		  // normal_lpdf(k_log | prior_k_log[1], prior_k_log[2]) +
	//   		  normal_lpdf(l_log | prior_l_log[1], prior_l_log[2]) +
	//   		  normal_lpdf(r_log | prior_r_log[1], prior_r_log[2]) +
	//   		  normal_lpdf(s_log | prior_s_log[1], prior_s_log[2]); // student_t_lpdf(s_log | nu_student, prior_s_log[1], prior_s_log[2]);
	//       			  
    // // for(l in 1:L_y) {
    // //   log_lik[l] = neg_binomial_0_lpmf(y[l] | y_hat_offset[l], phi_obs_rep[l], theta_rep[l]);
    // // }
    // 
    // for(i in 1:L_y) {
    //   log_lik[i] = neg_binomial_2_lpmf(y[i] | y_hat_offset[i], phi_obs_rep[i]); // offset_zeta
    // }
  
  } // END: if(generateposteriorq)
}
