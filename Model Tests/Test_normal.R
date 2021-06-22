## At some point I experienced a bias towards smaller estimates, which I wanted to explore here.
## Applied fix of an indexing error which might have fixed the bias.

library(deSolve)
library(cmdstanr)

n <- 20
times <- 1:20


calcODE <- function(t,
                    y, # A vector of species states.
                    r){
  r <- r[[1]]
  dy <- y * r
  return(list(c(dy)))
}

r <- 0.3
sigma <- 1
y_0_hat <- runif(n, 1, 2)
Y_hat <- ode(y_0_hat, c(0, times), calcODE, list(r))[,-1] # drop time column
Y <- matrix(rnorm(Y_hat, Y_hat, sigma), ncol = ncol(Y_hat), nrow = nrow(Y_hat))
y_0 <- Y[1,]
Y <- Y[-1,]


## Even a small overestimation if sigma leads to a large bias towards smaller r.
## This is dependent on the size of r, where with smaller r < 1 the effect gets more pronounced.
## This is why a strong regularization prior on sigma is needed.

stancode <- (" functions {
  vector ds_dt(real time, vector state, real r, int n_pops) {
    vector[n_pops] ds = r * state;
    return ds;
  }

}

data {
  int<lower=0> N;
  int<lower=0> N_times;

  real time[N_times];
  real time_init;
  vector<lower=0>[N] y_init;

  vector<lower=0>[N] Y[N_times];
}

parameters {
  real<lower=0> r;
  real<lower=0> sigma;
  vector<lower=0>[N] state_init;
}

transformed parameters {
  vector<lower=0>[N] State[N_times];

  State = ode_rk45(ds_dt,
              state_init, time_init, time,
              r, // model parameters
              N);
}

model {
  // this isn't enough as regularizing prior:
  sigma ~ cauchy(0, 0.001); // sigma is constrained to be > 0, thus half-cauchy

  // Model
  y_init ~ normal(state_init, sigma); // Separately fitting initial state to feed to integrator.

  for (t in 1:N_times) {
    Y[t] ~ normal(State[t], sigma); // Lognormal model of states.
  }
  }

")


# f_normal <- write_stan_file(stancode)
model_normal <- cmdstan_model(f_normal)

fit_normal <- model_normal$optimize(data = list(N = n,
                                              N_times = length(times),
                                              time = times,
                                              time_init = 0,
                                              y_init = y_0,
                                              Y = Y),
                                  iter = 10000)

fit_normal$summary()


