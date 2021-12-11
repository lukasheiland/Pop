functions {
  
  //// Difference equations
  matrix simulate(vector initialstate, int time_max,
                  vector b, vector c_a, vector c_b, vector c_j, vector g,  vector h, vector l, vector r, vector s, 
                  vector ba_a_avg, real ba_a_upper,
                  int N_spec, int N_pops,
                  int[] i_j, int[] i_a, int[] i_b) {
    
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
  
  //// ODE integration and data assignment to one long vector is wrapped into a function here because
  // - the transformed parameters block does not allow declaring integers (necessary for ragged data structure indexing),
  // - the model block does not alllow for declaring variables within a loop.
  vector unpack(vector[] state_init_log, int[] time_max, int[] times,
                vector b_log, vector c_a_log, vector c_b_log, vector c_j_log, vector g_log, vector h_log, vector[] L_loc, vector r_log, vector s_log,
                // vector b_log, vector c_a_log, vector c_b_log, matrix C_j_log, matrix G_log, vector h_log, vector[] L_loc, matrix R_log, matrix S_log,
                vector ba_a_avg, real ba_a_upper,
                int[] n_obs, int[] n_yhat,
                int N_species, int N_pops, int L_yhat, int N_locs,
                int[] i_j, int[] i_a, int[] i_b) {

    int pos_times = 1; // segmenting times[L_y]
    int pos_yhat = 1;
    vector[L_yhat] y_hat;

    for (loc in 1:N_locs) {
      int n_o = n_obs[loc];
      int n_y = n_yhat[loc];
      
      //// Returns a matrix State[p, t]
      // print("r", R_log);
      matrix[N_pops, time_max[loc]] States =
                        
                        simulate(exp(state_init_log[loc]),
                                 time_max[loc],
                                 exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log),
                                 ba_a_avg, ba_a_upper,
                                 N_species, N_pops,
                                 i_j, i_a, i_b);
      
  
      // Flattening the matrix into a vector for the location and append it to y_hat local vector yhat[m], and then into function-global vector y_hat[L_yhat].
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
      
      // s_1[i_j]  =  ((r .* BA)/(1 + BA_sum) + l + (J - g .* J)) ./ (1 + c_j*sum(J) + s*BA_sum);
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
//  vector[] transformToVertex(matrix P) {
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


//// Implementation of negbinomial probability density with zero inflation
//real neg_binomial_0_lpmf(int y, real y_hat, real phi_obs, real theta) {
//   
//  real t; // target
//  // From a process point of view this is just a negbinom model, with some values multiplied by zero.
//  // rbinom(100, 1, prob = 0.2) * rnbinom(100, size = 1100, mu = 10)
//  
//  if (y == 0) {
//    // Joint Likelihood of 0 coming from probability theta or negbinonial
//  	t = log_sum_exp(bernoulli_lpmf(1 | theta),
//                       bernoulli_lpmf(0 | theta) + neg_binomial_2_lpmf(y | y_hat, phi_obs));
//  } else {
//	// Joint Likelihood of 0 coming from probability theta_rep or negbinonial
//  	t = bernoulli_lpmf(0 | theta) +  // log1m(theta) synonymous to bernoulli_lpmf(0 | theta_rep)?
//  		neg_binomial_2_lpmf(y | y_hat, phi_obs);
//  }
//  return t; // target wich will get summed up at each run
// }
 
 
//// Implementation of zi negbinomial random number generator
//real[] neg_binomial_0_rng(vector y_hat_rep, vector phi_obs_rep, vector theta_rep, int L_y) {
//   
//   array[L_y] real n_io;
//   array[L_y] real io = bernoulli_rng(theta_rep); // variable assignment operator to force array real instead of int
//   array[L_y] real n = neg_binomial_2_rng(y_hat_rep, phi_obs_rep);
//   n_io = to_array_1d(to_vector(n) .* to_vector(io));
//
//   return n_io;
// }



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
  int<lower=0> L_yhat; // locations/surveys/pops
  int<lower=0> L_y; // locations/resurveys/pops/plots
  // int<lower=0> L_plots; // locations/plots
  // int<lower=0> L_init;
  // int<lower=0> L_a2b; // locations/obsid-1/n_species

  int<lower=0> N_locs; // overall number of locations. This is the major running variable to iterate over the model vectors.
  int<lower=0> N_species; // overall number of unique species across plot and locations (not nested!)
  int<lower=0> N_pops; // (species*stages) within loc; this is the length of initial values!
  int<lower=0> N_beta;
  int<lower=0> N_obsmethodTax;
  int<lower=0> N_protocol;


  //// n - number of levels within locs for more ragged models
  int<lower=0> n_obs[N_locs]; // n solution times within locs
  int<lower=0> n_yhat[N_locs]; // number of pops*obs within loc

  //// i — indices of stages
  array[N_species] int<lower=1> i_j; // e.g. 1:4
  array[N_species] int<lower=1> i_a; // e.g. 5:9, etc.
  array[N_species] int<lower=1> i_b;

  //// rep - repeat indices within groups for broadcasting to the subsequent hierarchical levels
  int<lower=1> rep_yhat2y[L_y]; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  int<lower=1> rep_obsmethodTax2y[L_y]; // factor (1, 2, 3), repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  int<lower=1> rep_protocol2y[L_y]; // factor (1:5)

  // int<lower=1> rep_locs2plots[L_plots]; // repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
  // int<lower=1> rep_yhat2a2b[L_a2b];
  // int<lower=1> rep_species2a2b[L_a2b];
  // int<lower=1> rep_pops2init[L_init];


  //// actual data
  int time_max[N_locs];
  int times[L_times]; // locations/observations
  // vector<lower=1>[L_a2b] timediff;
  
  matrix[N_locs, N_beta] X; // design matrix
  array[N_locs] vector[N_species] L_smooth_log;
  
  vector<lower=0>[L_y] offset;
  
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
  real<lower=0> parfactor;
  real<lower=0> timestep;
  
  //// Priors. The 2 reflect the two parameters mu and sigma
  // environmentally-dependent priors are species-agnostic on purpose
  // assumes vertex form f(x,y) == a_1*(x−p)^2 + a_2*(y−q)^2 + z
  // and parameters vector[z, p, q, a_1, a_2]
  
  array[N_locs] vector[N_pops] Prior_state_init_log;

  vector[2] prior_b_log;
  vector[2] prior_c_a_log;
  vector[2] prior_c_b_log;
  vector[2] prior_c_j_log;

  array[2] vector[N_species] prior_g_log;
  array[2] vector[N_species] prior_h_log;
  
  array[2] vector[N_species] prior_l_log;
  array[2] vector[N_species] prior_r_log;
  
  // vector[2] prior_l_log;
  // vector[2] prior_r_log;
  
  vector[2] prior_s_log;
  
  // array[2] vector[N_beta] prior_Vertex_c_j;
  // array[2] vector[N_beta] prior_Vertex_g;
  // array[2] vector[N_beta] prior_Vertex_r;
  // array[2] vector[N_beta] prior_Vertex_s;
  
  // array[2] vector[N_species] prior_b_log;
  // array[2] vector[N_species] prior_c_a_log;
  // array[2] vector[N_species] prior_c_b_log;
  
  // array[2] vector[N_species] prior_l_log;
  

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
  int N_genstates = N_pops + N_species;
  real factor_log = log(parfactor/timestep);
  
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
  
  vector[N_species] r_log;
  vector[N_species] l_log;
  
  // … dependent on environment. Matrix for convenient matrix multiplication
  // matrix[N_beta, N_species] Beta_c_j;
  // matrix[N_beta, N_species] Beta_g; // (J-, A+) transition rate from J to A
  // matrix[N_beta, N_species] Beta_r; // // (J+) flow into the system, dependent on env
  // matrix[N_beta, N_species] Beta_s; // (J-), here still unconstrained on log scale, shading affectedness of juveniles from A
  
  //// Special case l  
  matrix[N_locs, N_species] L_random_log; // array[N_locs] vector[N_species] L_random_log; // here real array is used for compatibility with to_vector
  vector<lower=0>[N_species] sigma_l;
  
  
  //// Errors
  vector<lower=0>[N_obsmethodTax] phi_obs_inv_sqrt; // error in neg_binomial per tax and stage
  
    // vector<lower=0>[N_protocol] zeta; // zero-offset parameter
	// real<lower=0> kappa_inv; // error in beta for h_log
    // vector<lower=0>[3] alpha_obs_inv; // observation error in gamma
    // vector<lower=0>[2] sigma_obs; // observation error
    // vector<lower=0>[3] sigma_process; // lognormal error for observations from predictions
  
  // vector<lower=0,upper=1>[3] theta_obs; // observation error in neg_binomial
  
  // matrix[N_pops, timespan_max] u[N_locs];
  
  ///
  array[N_locs] vector[N_pops] state_init_log_raw;
}


transformed parameters {

  array[N_locs] vector<lower=0>[N_species] L_loc;
  array[N_locs] vector[N_pops] state_init_log;
  // vector[L_y] offset_zeta;
  
  
  //// Level 1 (species) to 2 (locs). Environmental effects on population rates.
  // Location level parameters (unconstrained because on log scale!)
  // matrix[N_locs, N_species] C_j_log = X * Beta_c_j;
  // matrix[N_locs, N_species] G_log = X * Beta_g;
  // matrix[N_locs, N_species] R_log = X * Beta_r;
  // matrix[N_locs, N_species] S_log = X * Beta_s;

  vector<lower=0>[N_obsmethodTax] phi_obs = inv_square(phi_obs_inv_sqrt); // inv_square == square_inv
    // vector<lower=0>[3] alpha_obs = inv(alpha_obs_inv);
  
  for(loc in 1:N_locs) {
    
    L_loc[loc, ] = exp(l_log + L_smooth_log[loc, ] + // The smooth effect
                       sigma_l .* L_random_log[loc, ]'); // non-centered loc-level random intercept 
                   
    state_init_log[loc] = state_init_log_raw[loc] * 1;
  }
  
  //  vector[L_y] zeta_rep = zeta[rep_protocol2y];
  //  for(j in 1:L_y) {
  //  	if (offset[j] == 0) {
  //  		offset_zeta[j] = zeta_rep[j];
  //  	} else {
  //  		offset_zeta[j] = offset[j];
  //  	}
  //  }
  
  
  vector<lower=0>[L_yhat] y_hat = unpack(state_init_log, time_max, times,
                                b_log, c_a_log, c_b_log, c_j_log, g_log, h_log, L_loc, r_log, s_log, // rates matrix[N_locs, N_species]; will have to be transformed
                                // b_log, c_a_log, c_b_log, C_j_log, G_log, h_log, L_loc, R_log, S_log, // rates matrix[N_locs, N_species]; will have to be transformed
                                ba_a_avg, ba_a_upper,
                                n_obs, n_yhat, // varying numbers per loc
                                N_species, N_pops, L_yhat, N_locs, // fixed numbers
                                i_j, i_a, i_b);
                                
  vector[L_y] y_hat_rep = y_hat[rep_yhat2y];
  vector[L_y] y_hat_rep_offset = y_hat_rep .* offset; // offset_zeta
  vector[L_y] phi_obs_rep = phi_obs[rep_obsmethodTax2y];
  // vector[L_y] theta_obs_rep = theta_obs[rep_obsmethodTax2y];
  
}


model {

  //—————————————————————————————————————————————————————————————————————//
  // Priors       -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  //// Hyperpriors

  phi_obs_inv_sqrt ~ normal(rep_vector(0.0, 6), [0.2, 0.6, 0.05, 0.2, 0.6, 0.02]); // Observation error for neg_binomial
  	// Levels of obsmethodTax: j.F a.F ba.F j.o a.o ba.o
  	// On prior choice for the overdispersion in negative binomial 2: https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#story-when-the-generic-prior-fails-the-case-of-the-negative-binomial
  
  // ... for special offset L
  to_vector(L_random_log) ~ std_normal(); // Random intercept for l
  sigma_l ~ std_normal(); // Regularizing half-cauchy on sigma for random slope for l  ## cauchy(0, 2);


  // sigma_process ~ normal(0, 0.01);
  // sigma_obs ~ normal(0, [0.5, 0.1]); // for observations from predictions
  // alpha_obs_inv ~ normal(0, 0.1); // Observation error for gamma
  // zeta ~ normal(0, 1);

  
  
  //// Priors for Parameters
  
  // prior_*[2, N_species]
  
  // b_log   ~ normal(prior_b_log[1,], prior_b_log[2,]);
  // c_a_log ~ normal(prior_c_a_log[1,], prior_c_a_log[2,]);
  // c_b_log ~ normal(prior_c_b_log[1,], prior_c_b_log[2,]);
  
  b_log   ~ normal(prior_b_log[1], prior_b_log[2]);
  c_a_log ~ normal(prior_c_a_log[1], prior_c_a_log[2]);
  c_b_log ~ normal(prior_c_b_log[1], prior_c_b_log[2]);
  c_j_log ~ normal(prior_c_j_log[1], prior_c_j_log[2]); // strong believe that c is smaller than s in trees
  
  g_log ~ normal(prior_g_log[1,], prior_g_log[2,]);
  h_log ~ normal(prior_h_log[1,], prior_h_log[2,]);

  
  l_log ~ normal(prior_l_log[1], prior_l_log[2]);
  r_log ~ normal(prior_r_log[1], prior_r_log[2]); // wanna constrain this a bit, otherwise the model will just fill up new trees and kill them off with g
  
  s_log ~ normal(prior_s_log[1], prior_s_log[2]);

  
  // same priors for both species
  // for (spec in 1:N_species) {
  //   // prior_Beta_*[2, N_beta]
  //   Beta_c_j[, spec] ~ normal(prior_Beta_c_j[1, ], prior_Beta_c_j[2, ]);
  //   Beta_g[, spec]  ~ normal(prior_Beta_g[1, ], prior_Beta_g[2, ]);
  //   Beta_r[, spec] ~ normal(prior_Beta_r[1, ], prior_Beta_r[2, ]);
  //   Beta_s[, spec] ~ normal(prior_Beta_s[1, ], prior_Beta_s[2, ]);
  // }
  
   for(loc in 1:N_locs) {
     state_init_log_raw[loc] ~ std_normal();
   }


  //—————————————————————————————————————————————————————————————————————//
  // Model       -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    

  // a2b ~ poisson(y_hat[rep_yhat2a2b] .* exp(h_log)[rep_species2a2b] .* timediff);
  // y ~ neg_binomial_2(y_hat[rep_yhat2y], phi_obs[rep_obsmethodTax2y]);
  
  y ~ neg_binomial_2(y_hat_rep_offset, phi_obs_rep);

  
  // y ~ gamma_0(y_hat[rep_yhat2y], alpha_obs[rep_obsmethodTax2y], theta_obs, L_y);
  
  //  for(l in 1:L_y) {
  //	
  //	y[l] ~ neg_binomial_0(y_hat_rep[l], phi_obs_rep[l], theta_obs_rep[l]);
  //    // y[l] ~ neg_binomial_0(y_hat[rep_yhat2y][l], phi_obs[rep_obsmethodTax2y][l], theta_obs[rep_obsmethodTax2y][l]);
  //  }  

}


generated quantities {

  //—————————————————————————————————————————————————————————————————————//
  // Prediction  -------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
	
  // vector[N_locs] y_hat_temp;
  // y_hat_temp = y_hat;
  array[L_y] real<lower=0> y_sim;
  
  //// y_hat_rep in transformed parameters
  // y_sim = neg_binomial_0_rng(y_hat_rep, phi_obs_rep, theta_obs_rep, L_y);
  y_sim = neg_binomial_2_rng(y_hat_rep_offset, phi_obs_rep);



  //—————————————————————————————————————————————————————————————————————//
  // Prior predictive checks -------------------------------------------//
  //———————————————————————————————————————————————————————————————————//    
  
  ///// Priors -----------------------------------
  real b_log_prior = normal_rng(prior_b_log[1], prior_b_log[2]);
  real c_a_log_prior = normal_rng(prior_c_a_log[1], prior_c_a_log[2]);
  real c_b_log_prior = normal_rng(prior_c_b_log[1], prior_c_b_log[2]);
  real c_j_log_prior = normal_rng(prior_c_j_log[1], prior_c_j_log[2]); // strong believe that c is smaller than s in trees
  
  vector<upper=0>[N_species] g_log_prior = -sqrt(square(to_vector(normal_rng(prior_g_log[1,], prior_g_log[2,]))));
  vector<upper=0>[N_species] h_log_prior = -sqrt(square(to_vector(normal_rng(prior_h_log[1,], prior_h_log[2,]))));
  
  vector[N_species] l_log_prior = to_vector(normal_rng(prior_l_log[1,], prior_l_log[2,]));
  vector[N_species] r_log_prior = to_vector(normal_rng(prior_r_log[1,], prior_r_log[2,]));
  
  // real l_log_prior = normal_rng(prior_l_log[1], prior_l_log[2]);
  // real r_log_prior = normal_rng(prior_r_log[1], prior_r_log[2]);
  
  real s_log_prior = normal_rng(prior_s_log[1], prior_s_log[2]);
  
  vector[N_species] vector_b_log_prior = to_vector(normal_rng(rep_vector(prior_b_log[1], N_species), rep_vector(prior_b_log[2], N_species)));
  vector[N_species] vector_c_a_log_prior = to_vector(normal_rng(rep_vector(prior_c_a_log[1], N_species), rep_vector(prior_c_a_log[2], N_species)));
  vector[N_species] vector_c_b_log_prior = to_vector(normal_rng(rep_vector(prior_c_b_log[1], N_species), rep_vector(prior_c_b_log[2], N_species)));
  vector[N_species] vector_c_j_log_prior = to_vector(normal_rng(rep_vector(prior_c_j_log[1], N_species), rep_vector(prior_c_j_log[2], N_species)));
  // vector[N_species] vector_l_log_prior = to_vector(normal_rng(rep_vector(prior_l_log[1], N_species), rep_vector(prior_l_log[2], N_species)));
  // vector[N_species] vector_r_log_prior = to_vector(normal_rng(rep_vector(prior_r_log[1], N_species), rep_vector(prior_r_log[2], N_species)));
  vector[N_species] vector_s_log_prior = to_vector(normal_rng(rep_vector(prior_s_log[1], N_species), rep_vector(prior_s_log[2], N_species)));
  
  array[N_obsmethodTax] real<lower=0> phi_obs_prior = inv_square(normal_rng(rep_vector(0.0, 6), [1, 1, 0.1, 1, 1, 0.1])); // [0.2, 0.6, 0.05, 0.2, 0.6, 0.02]
  
  // special case L
  array[N_locs] vector<lower=0>[N_species] L_loc_prior;
  array[N_locs, N_species] real L_random_log_prior;  
  vector<lower=0>[N_species] sigma_l_prior = sqrt(square(to_vector(normal_rng(rep_vector(0, N_species), rep_vector(1, N_species)))));
  // vector<lower=0>[N_protocol] zeta_prior = sqrt(square(to_vector(normal_rng(rep_vector(0, N_protocol), rep_vector(0.2, N_protocol)))));
  
  for(loc in 1:N_locs) {
  
    L_random_log_prior[loc,] = normal_rng(rep_vector(0, N_species), rep_vector(1, N_species));
    L_loc_prior[loc, ] = exp(l_log_prior + L_smooth_log[loc, ] +
                             sigma_l_prior .* to_vector(L_random_log_prior[loc, ]));
                         
    
  }
  
  
  //// Prior simulation
  vector<lower=0>[L_yhat] y_hat_prior;
  vector<lower=0>[L_y] y_hat_prior_rep;
  vector<lower=0>[L_y] y_hat_prior_rep_offset;
  // vector<lower=0>[L_y] offset_zeta_prior;
  array[L_y] real<lower=0> y_prior_sim;
                   
  //  vector<lower=0>[L_y] zeta_prior_rep = zeta_prior[rep_protocol2y];
  
  //  for(j in 1:L_y) {
  //  	if (offset[j] == 0) {
  //  		offset_zeta_prior[j] = zeta_prior_rep[j];
  //  	} else {
  //  		offset_zeta_prior[j] = offset[j];
  //  	}
  //  }
 
  y_hat_prior = unpack(Prior_state_init_log, time_max, times,
                       vector_b_log_prior, vector_c_a_log_prior, vector_c_b_log_prior, vector_c_j_log_prior, g_log_prior, h_log_prior, L_loc_prior, r_log_prior, vector_s_log_prior, // rates matrix[N_locs, N_species]; will have to be transformed
                       ba_a_avg, ba_a_upper,
                       n_obs, n_yhat, // varying numbers per loc
                       N_species, N_pops, L_yhat, N_locs, // fixed numbers
                       i_j, i_a, i_b);
 
  y_hat_prior_rep = y_hat_prior[rep_yhat2y];
  y_hat_prior_rep_offset = y_hat_prior_rep .* offset; // offset_zeta_prior


  y_prior_sim = neg_binomial_2_rng(y_hat_prior_rep_offset, phi_obs_prior[rep_obsmethodTax2y]); // , [0, 0, 0]', L_y
  
  
  //—————————————————————————————————————————————————————————————————————————//
  // pgq -------------------------------------------------------------------//
  //———————————————————————————————————————————————————————————————————————//
  
  //// Declarations of posterior quantites (as global variables).
  // … are directly initiated with zeroes or 9, so that there are never NaNs in generated quantities.
  int fixiter_max = 5000; 
  array[N_locs] int converged = rep_array(9, N_locs); // tolerance has been reached
  array[N_locs] real iterations_fix = rep_array(0, N_locs);
  array[N_locs] vector[N_genstates+N_species+1] state_fix = rep_array(rep_vector(0, N_genstates+N_species+1), N_locs); // state_fix is a vector [J1, …, A1, …, B1, …, BA1, …, eps_ba1, …, iterations]
  array[N_locs] int dominant_fix = converged;
  array[N_locs] int major_fix = converged;
  

  //// Declarations of quantities for sensitivity checks (as global variables).
  real log_prior = 0; // this is zero to prevent NaNs from being in the sum.
  vector[L_y] log_lik = rep_vector(0, L_y);


  //// The conditional generation -------------------------------------
  if (generateposteriorq) {

    //———————————————————————————————————————————————————————————————————//
    // Posterior quantities --------------------------------------------//
    //—————————————————————————————————————————————————————————————————//
    
    //// Fix point iteration -------------------------------------------
    for(loc in 1:N_locs) {
      
      //// fix point, given parameters
      state_fix[loc] = iterateFix(exp(state_init_log[loc]),
                                  exp(b_log), exp(c_a_log), exp(c_b_log), exp(c_j_log), exp(g_log), exp(h_log), L_loc[loc, ], exp(r_log), exp(s_log),
                                  // exp(b_log), exp(c_a_log), exp(c_b_log), exp(C_j_log[loc,]'), exp(G_log[loc,]'), exp(h_log), exp(L_loc[loc, ]), exp(R_log[loc,]'), exp(S_log[loc,]'),
                                  ba_a_avg, ba_a_upper,
                                  N_species, N_pops,
                                  i_j, i_a, i_b,
                                  tolerance_fix, fixiter_max);
                                     
      iterations_fix[loc] = state_fix[loc, N_genstates+N_species+1];
      converged[loc] = iterations_fix[loc] < fixiter_max;
      
      if (converged[loc]) { // && convergent[loc]
        
        dominant_fix[loc] = state_fix[loc, N_pops+1]/state_fix[loc, N_genstates] > 3; // BA_1 > 75%
        major_fix[loc] = state_fix[loc, N_pops+1] > state_fix[loc, N_genstates]; // BA_1 > 50%
        
      }
      
     //    if (y_hat_temp[loc] > 1e+07) {
     //    	y_hat_temp[loc] = 1e+07;
     //    } 
  
    }
  
  
  
    //———————————————————————————————————————————————————————————————————//
    // Sensitivity analysis --------------------------------------------//
    //—————————————————————————————————————————————————————————————————//
  
    for(loc in 1:N_locs) {
      log_prior += normal_lpdf(state_init_log[loc,] | Prior_state_init_log[loc,], 1);
    }
    
    // [0.2, 0.6, 0.05, 0.2, 0.6, 0.02]
    log_prior = log_prior +
    		  normal_lpdf(phi_obs_inv_sqrt | rep_vector(0.0, 6), [1, 1, 0.1, 1, 1, 0.1]) +
	  		  normal_lpdf(sigma_l | 0, 1) +		  
	  		  normal_lpdf(to_vector(L_random_log) | 0, 1) +
	  		  normal_lpdf(b_log | prior_b_log[1], prior_b_log[2]) +
	  		  normal_lpdf(c_a_log | prior_c_a_log[1], prior_c_a_log[2]) +
	  		  normal_lpdf(c_b_log | prior_c_b_log[1], prior_c_b_log[2]) +
	  		  normal_lpdf(c_j_log | prior_c_j_log[1], prior_c_j_log[2]) +
	  		  normal_lpdf(g_log | prior_g_log[1,], prior_g_log[2,]) +
	  		  normal_lpdf(h_log | prior_h_log[1,], prior_h_log[2,]) +
	  		  normal_lpdf(l_log | prior_l_log[1], prior_l_log[2]) +
	  		  normal_lpdf(r_log | prior_r_log[1], prior_r_log[2]) +
	  		  normal_lpdf(s_log | prior_s_log[1], prior_s_log[2]); // joint prior specification, sum of all logpriors // normal_lpdf(zeta | rep_vector(0.0, 5), rep_vector(0.2, 5)) + log(2) +
	      			  
    // for(l in 1:L_y) {
    //   log_lik[l] = neg_binomial_0_lpmf(y[l] | y_hat_rep[l], phi_obs_rep[l], theta_obs_rep[l]);
    // }
    
    for(i in 1:L_y) {
      log_lik[i] = neg_binomial_2_lpmf(y[i] | y_hat_rep_offset[i], phi_obs_rep[i]); // offset_zeta
    }
  
  } // END: if(generateposteriorq)
}
