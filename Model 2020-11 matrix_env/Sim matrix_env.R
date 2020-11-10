# Library -----------------------------------------------------------------
library(here)

library(tidyverse)
library(magrittr)
library(glue)

library(deSolve)
library(cmdstanr)
# install_cmdstan(cores = 3)


# Orientation -------------------------------------------------------------

modelname <- "matrix_env"

setwd(here())
modeldir <- list.files(pattern = glue("^(Model).*{modelname}$"))
modelpath <- file.path(modeldir, glue('Model {modelname}.stan'))


# Multi-species matrix model -------------------------------------------------------------------

calcModel <- function(t,
                    m, # A vector of species states.
                    par){
  f <- par[[1]] # Length n vector of growth rates.
  A <- par[[2]] # Full n x n matrix of interaction factors.

  dm <- f + A %*% m
  return(list(c(dm)))
}


# Simulation -------------------------------------------------------------------

## Return indices for the layout of an interaction matrix.
layoutMatrix <- function(n_species, includespeciescolumn = T){
  n_pops <- n_species * 2
  species <- 1:n_species
  isjuvenile <- rep(c(T, F), length.out = n_pops)

  i_j <- which(isjuvenile)
  i_a <- which(!isjuvenile)

  # Species specific parameter vectors with length n_species
  I_s <- cbind(spec = rep(species, each = n_species),
               row = rep(i_j, each = n_species),
               col = i_a)

  I_0 <- cbind(spec = rep(species, each = n_species-1),
               row = i_a,
               col = rev(rep(i_j, each = n_species-1))) # indives with no competition of J on A

  I_g <- cbind(spec = species, row = i_a, col = i_j)

  # These are also in the small matrices
  I_lg <- cbind(spec = species,
                row = i_j,
                col = i_j)
  I_l <- cbind(spec = species,
               row = i_a,
               col = i_a)

  # Map the juvenile and adult matrices to a big matrix
  Map_j_intra <- cbind(row_a = species, col_a = species,
                       row_A = i_j, col_A = i_j)
  Map_a_intra <- cbind(row_a = species, col_a = species,
                       row_A = i_a, col_A = i_a)

  Full <- expand.grid(col = 1:n_species, row = 1:n_species)
  # Withoutdiag <- Full[-seq(1, (n_species^2), by = n_species),]

  Map_j <- cbind(row_a = Full$row, col_a = Full$col,
                 row_A = rep(i_j, each = n_species), col_A = i_j)

  Map_a <- cbind(row_a = Full$row, col_a = Full$col,
                 row_A = rep(i_a, each = n_species), col_A = i_a)

  # Test <- matrix("", n_pops, n_pops)
  # a_juvenile <- matrix(paste0("j", 1:n_species^2), n_species, n_species, byrow = 2)
  # a_adult <- matrix(paste0("a", 1:n_species^2), n_species, n_species, byrow = 2)
  #
  # Test[Map_j[,3:4]] <- a_juvenile[Map_j[,1:2]]
  # Test[Map_a[,3:4]] <- a_adult[Map_a[,1:2]]
  #
  # Test[I_s[,2:3]] <- rep(paste0("s", 1:n_species), each = n_species)
  # Test[I_0[,2:3]] <- "0"
  # Test[I_g[,2:3]] <- paste0("g", 1:n_species)
  # Test[I_lg[,2:3]] <- paste0("lg", 1:n_species)
  # Test[I_l[,2:3]] <- paste0("l", 1:n_species)

  i <- c(includespeciescolumn, T, T)
  return(list(I_s = I_s[, i],
              I_0 = I_0[, i],
              I_g = I_g[, i],
              I_lg = I_lg[, i],
              I_l = I_l[, i],
              Map_j = Map_j,
              Map_a = Map_a,
              i_j = i_j,
              i_a = i_a,
              isjuvenile = isjuvenile
  )
  )
}

#### Get environmental variables
simulateEnv <- function(n_env, n_locs){
  n <- n_env * n_locs
  Env <- matrix(runif(n, -1, 1), nrow = n_locs, ncol = n_env)
  colnames(Env) <- paste("env", 1:n_env, sep = "_")
  return(Env)
}


## Beta_* are effect of environment on species' parameters matrix[1 + n_env, n_species]
simulateParametersInEnv <- function(Env, Beta_r, Beta_s, Beta_g, returndf = F) {
  n_species <- ncol(Beta_r)
  M <- model.matrix(~ poly(Env, 2))
  r_hat <- M %*% Beta_r
  s_hat <- M %*% Beta_s
  g_hat <- M %*% Beta_g
  pars <- list(r_hat = r_hat, s_hat = s_hat, g_hat = g_hat)

  if (returndf) {
    pars <- lapply(pars, as.data.frame)
    pars <- do.call(cbind, pars)
    pars <- pivot_longer(cbind(pars, Env),
                         cols = starts_with(c("r_hat", "s_hat", "g_hat")),
                         names_sep = "\\.",
                         names_to = c("parameter", "species"),
                         values_to = c("q")) %>%
      mutate(species = as.integer(as.factor(species)))
  }

  return(pars)
}


#### Set parameters for the population model
## three parameters dependent on environment: r, s, g.
getParameters <- function(r_hat, s_hat, g_hat, returndf = F) {
  n_species <- length(r_hat)
  n_pops <- n_species * 2
  I <- layoutMatrix(n_species, includespeciescolumn = T)

  ## 1. regneration rate r
  r <- rlnorm(r_hat, r_hat, 0.0001) # exp(log(m)), log link!

  f <- rlnorm(n_pops, r, 0.0001) # Vector of intrinsic growth rates.
  f[I$isjuvenile] <- r # Vector of intrinsic growth rates.
  f[!I$isjuvenile] <- - f[!I$isjuvenile]/10 # output rates

  A <- matrix(NA, n_pops, n_pops)

  ## Competition of J on A
  A[I$I_0[,2:3]] <- 0

  ## (-) Competition of A on J: s[n_species]
  ## 2. juveniles experience same shading s from all species
  s <- -rlnorm(s_hat, s_hat, 0.0001)
  A[I$I_s[,2:3]] <- rep(s, each = n_species) # -> s[n_species^2]

  ## (+) Transition rate from J to A: g[n_species]
  ## 3. g
  g <- rlnorm(g_hat, g_hat, 0.0001)
  A[I$I_g[,2:3]] <- g

  ## (-) Competition matrix of adults: a_A
  A[I$Map_a[,3:4]] <- -rlnorm(n_species^2, log(0.05), 0.05)

  ## (-) Competition matrix of juveniles: a_J
  A[I$Map_j[,3:4]] <- -rlnorm(n_species^2, log(0.02), 0.01)

  ## (-) Intra competition of adults (the diagonal of a_A): l[n_species]
  A[I$I_l[,2:3]] <- -rlnorm(n_species, log(0.1), 0.1)

  ## (-) Intra competition of juveniles - transition from juveniles: l - g = lg[n_species]
  A[I$I_lg[,2:3]] <- A[I$I_lg[,2:3]] - g

  par <- list(f = f, A = A) # Parameters list, including a matrix of alpha values.

  if(returndf) {
    par <- data.frame(
      r = t(r_hat),
      s = t(s_hat),
      g = t(g_hat),
      species = 1:n_species
    )
  }

  return(par)
}


#### Simulate vector of initial states
simulateInitialState <- function(n_pops) {
  runif(n_pops, log(10), 3) * c(2, 1) # Initial state matrix.
}

#### Integrate model to simulate states
simulateSeries <- function(state_init, times, par, sigma) {

  m0 <- state_init

  sim0 <- c(time = 0, pop = m0)
  Sim <- ode(m0, times, calcModel, par)

  sim0[2:(n_pops+1)] <- rlnorm(sim0, log(sim0), sigma)[2:(n_pops+1)]
  Sim[, 2:(n_pops+1)] <- matrix(rlnorm(Sim, log(Sim), sigma), nrow = nrow(Sim))[, 2:(n_pops+1)]

  Sim[is.nan(Sim) | Sim < 0] <- 0

  return(rbind(sim0, Sim))
}


#### Multiple simulations
simulateSeriesInEnv <- function(n_species,
                                n_locs,
                                n_seriesperloc = 3,
                                n_times = 4,
                                Env,
                                Beta_r,
                                Beta_s,
                                Beta_g,
                                sigma = 0.1,
                                format = c("long", "wide", "list")) {

  times <- seq(10, 40, length.out = n_times)
  n_pops <- n_species * 2

  parsinenv <- simulateParametersInEnv(Env, Beta_r, Beta_s, Beta_g)

  simulateSeriesAtLoc <- function(r_hat, s_hat, g_hat) {
    par <- getParameters(r_hat, s_hat, g_hat)
    m0 <- simulateInitialState(n_pops = n_pops)
    seriesatloc <- replicate(n_seriesperloc, simulateSeries(m0, times, par, sigma), simplify = F)
    names(seriesatloc) <- 1:n_seriesperloc
    return(seriesatloc)
  }

  # Sims is a: list[n_loc] <- list[n_seriesperloc] <- matrix[n_times, n_pops]
  Sims <- mapply(simulateSeriesAtLoc,
                 as.data.frame(t(parsinenv$r_hat)),
                 as.data.frame(t(parsinenv$s_hat)),
                 as.data.frame(t(parsinenv$g_hat)),
                 SIMPLIFY = F)

  if (match.arg(format) %in% c("wide", "long")) {
    Sims <- lapply(Sims, function(l) do.call(rbind, l)) # bind the series within a location into one table
    Sims <- lapply(Sims, function(m) cbind(m, seriesid = rep(1:n_seriesperloc, each = n_times+1))) # attach an id to series per location

    Sims <- do.call(rbind, Sims) # bind all lists of locations into one table
    Sims <- cbind(Sims, locid = rep(1:n_locs, each = n_seriesperloc * (n_times+1))) # attach an id to series per location

    Sims <- cbind(Sims, Env[Sims[,"locid"],])
  }

  if (match.arg(format) == "long") {
    Sims <- tidyr::pivot_longer(as.data.frame(Sims),
                                cols = all_of(paste0("pop", 1:n_pops)),
                                names_to = "pop",
                                values_to = "abundance") %>%
      mutate(pop = as.integer(as.factor(pop))) %>%
      mutate(species = rep(1:n_species, each = 2, length.out = nrow(.)),
             stage = rep(c("j", "a"), length.out = nrow(.))
      )
  }
  return(Sims)
}



# Simulate! ---------------------------------------------------------------
n_species <- 4
n_locs <- 500
n_seriesperloc = 3
n_times <- 4

## Beta_* are effect of environment on species' parameters matrix[n_env, n_species]
## These are on the log scale!
n_env <- 3
n_beta <- 1 + ncol(poly(Env, 2)) # + intercept
Beta_r <- matrix(rlnorm(n_beta*n_species, log(0.05), 0.01), n_beta, n_species)
Beta_r[1,] <- log(2)
Beta_s <- matrix(runif(n_beta*n_species, log(0.05), 0.01), n_beta, n_species)
Beta_s[1,] <- log(0.1)
Beta_g <- matrix(runif(n_beta*n_species, log(0.05), 0.01), n_beta, n_species)
Beta_g[1,] <- log(0.2)

sigma <- 0.1

## Simulate one time series
initialstate1 <- simulateInitialState(n_species*2)
par1 <- getParameters(Beta_r[1,], Beta_s[1,], Beta_g[1,])
Sim1 <- simulateSeries(initialstate1, times = 1:100, par1, sigma = 0.02)
matplot(Sim1[, 1], Sim1[, -1], type = "b", ylab="N") # log='y'

## Simulate multiple time series in environmental space
E <- simulateEnv(n_env, n_locs)
P <- simulateParametersInEnv(E, Beta_r, Beta_s, Beta_g, returndf = T)
S <- simulateSeriesInEnv(n_species, n_locs, n_seriesperloc, n_times,
                         E, Beta_r, Beta_s, Beta_g,
                         sigma,
                         format = "long")


ggplot(filter(P, parameter == "r_hat"),
       mapping = aes(x = env_1, y = q, color = species, group = species)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(facets = c("species"))

ggplot(S,
       mapping = aes(x = time, y = abundance, color = env_1, lty = stage, group = interaction(locid, stage, seriesid))) +
  geom_line() +
  facet_wrap(facets = c("species"))



# Prepare stan model data --------------------------------------------------------------

#### Get a list of lists of data, where each second-level list belongs to one time series
## The lists will be unpacked as arrays of doubles, matrices, vectors etc. in stan
getStanData <- function(locations, Env) {

  ## environmental data
  X <- model.matrix(~ poly(Env, 2))
  envdata <- list(
    X = X,
    N_beta = ncol(X)
  )

  ## global numbers (For now, n_obs, and n_species have to be the same over all series. Thus, only first element extraction)
  n_pops <- ncol(locations[[1]][[1]]) - 1
  ndata <- list(
    N_locs = length(locations),
    N_seriesperloc = length(locations[[1]]),
    N_obs = nrow(locations[[1]][[1]]) - 1, # observations, other than at t_0
    N_species = n_pops/2,
    N_pops = n_pops
  )

  ## location-wide data
  ldata <- list(
    time_init = t(sapply(locations,
                       function(l) sapply(l, function(s) s[1, "time"]))), # stan expects a two-dimensional array time_init[N_locs, N_seriesperloc]
    times = aperm(sapply(locations,
                   function(l) sapply(l, function(s) s[2:nrow(s), "time"]),
                   simplify = "array"),
                  c(3, 2, 1)), # stan expects a three-dimensional array real times[N_locs,N_seriesperloc,N_obs]; // arrays are row major!
    y_init = aperm(sapply(locations,
                    function(l) sapply(l, function(s) s[1, -1]),
                    simplify = "array"),
                   c(3, 2, 1)), # stan expects three-dimensional array: real y_init[N_locs,N_seriesperloc,N_species];
    Y = aperm(sapply(locations,
               function(l) sapply(l, function(s) s[-1, -1], simplify = "array"),
               simplify = "array"),
              c(4, 3, 1, 2)) # stan expects three-dimensional array of n-dimensional vectors: vector[N_species] Y[N_locs,N_seriesperloc,N_obs] ## provided indexing order is opposite in R and stan!
  )

  return(c(envdata, ndata, ldata, layoutMatrix(n_species = n_pops/2)))
}


#### generate acceptable start values
getInits <- function() {
  list(
    g = rlnorm(data$N_species, log(1), 0.01),
    s = - rlnorm(data$N_species, log(2), 0.01),

    A_j = - matrix(rlnorm(data$N_species^2, log(0.5), 0.01), nrow = data$N_species),
    A_a = - matrix(rlnorm(data$N_species^2, log(0.5), 0.01), nrow = data$N_species),

    r = rlnorm(data$N_species, log(20), 0.01),
    o = - rlnorm(data$N_species, log(0.2), 0.01),

    state_init =  data$y_init, # t(replicate(data$N_series, runif(data$N_pops, 10, 20), simplify = "array")), # stan expects vector<lower=0>[N_species] state_init[N_series]

    sigma = 2, # rep(0.1, data$N_species)
    theta = 1
  )
}


# getInitsTrue <- function() {
#   list(
#     g = A[I$I_g[,2:3]], # (+)
#     s = A[I$I_s[,2:3]], # (-)
#
#     A_j = A[I$Map_j[,3:4]],
#     A_a = A[I$Map_a[,3:4]],
#
#     r = f[I$isjuvenile],
#     o = f[!I$isjuvenile],
#
#     state_init =  data$y_init, # stan expects vector<lower=0>[N_species] state_init[N_series]
#
#     sigma = 0.1 # rep(0.1, data$N_species)
#   )
# }



# Simulate and fit! -------------------------------------------------------
Env <- simulateEnv(n_env, n_locs)

sims <- simulateSeriesInEnv(n_species, n_locs, n_seriesperloc, n_times,
                    Env, Beta_r, Beta_s, Beta_g,
                    sigma,
                    format = "list")
data <- getStanData(sims, Env)


#### Compile model
model <- cmdstan_model(modelpath)
# model$exe_file() ## path to compiled model


#### Optimize model.
fit_optim <- model$optimize(data = data,
                            iter = 10^9,
                            output_dir = "Fits/",
                            init = getInitsTrue)

Estimates <- fit_optim$summary()

Compare <- cbind(estimate = Estimates[c(92:134, 32:91),],
                 true = c(sigma, c(A), f, c(data$y_init)))
View(Compare)


#### Sample.
n_chains <- 1
fit <- model$sample(data = data,
                    output_dir = "Fits",
                    init = getInits,
                    iter_warmup = 200, iter_sampling = 500,
                    chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))

fit$summary(variables = 'state_init')
data$y_init

# shinystanfit <- rstan::read_stan_csv(fit$output_files())
# shinystan::launch_shinystan(shinystanfit)

