# Library -----------------------------------------------------------------
library(here)

library(tidyverse)
library(magrittr)
library(glue)

library(deSolve)
library(cmdstanr)
# install_cmdstan(cores = 3)


# Orientation -------------------------------------------------------------

modelname <- "matrix-lognormal"

setwd(here())
modeldir <- dir(pattern = glue("^(Model).*{modelname}$"))
modelpath <- file.path(modeldir, glue('Model_{modelname}.stan'))


# Multi-species matrix model -------------------------------------------------------------------

calcModel <- function(t,
                    m, # A vector of species states.
                    par){
  f <- par[[1]] # Length n vector of growth rates.
  A <- par[[2]] # Full n x n matrix of interaction factors.

  dm <- f + A %*% m # Here, limitation emerges only from the diagonal.
  return(list(c(dm)))
}


# Simulation -------------------------------------------------------------------

## Return indices for the layout of an interaction matrix.
layoutMatrix <- function(n_species, includespeciescolumn = T){
  n <- n_species * 2
  species <- 1:n_species
  isjuvenile <- rep(c(T, F), length.out = n)

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

  # Test <- matrix("", n, n)
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


#### Set parameters for the population model
n_species <- 3
n <- n_species * 2
I <- layoutMatrix(n_species, includespeciescolumn = T)

f <- rlnorm(n, log(10), 0.1) # Vector of intrinsic growth rates.
# f[!I$isjuvenile] <- 0
f[!I$isjuvenile] <- -f[!I$isjuvenile]/10 # death rates

A <- matrix(NA, n, n)

## Competition of J on A
A[I$I_0[,2:3]] <- 0

## (-) Competition of A on J: s[n_species^2]
## juveniles experience same shading from all species
A[I$I_s[,2:3]] <- rep(-rlnorm(n_species, log(0.5), 0.01), each = n_species)

## (+) Transition rate from J to A: g[n_species]
g <- rlnorm(n_species, log(0.2), 0.01)
A[I$I_g[,2:3]] <- g

## (-) Competition matrix of adults: a_A
A[I$Map_a[,3:4]] <- -rlnorm(n_species^2, log(0.1), 0.05)

## (-) Competition matrix of juveniles: a_J
A[I$Map_j[,3:4]] <- -rlnorm(n_species^2, log(0.05), 0.01)

## (-) Intra competition of adults (the diagonal of a_A): l[n_species]
A[I$I_l[,2:3]] <- -rlnorm(n_species, log(0.25), 0.1)

## (-) Intra competition of juveniles - transition from juveniles: l - g = lg[n_species]
A[I$I_lg[,2:3]] <- A[I$I_lg[,2:3]] - g

## error
sigma <- 0.1

par <- list(f, A) # Parameters list, including a matrix of alpha values.

#### Integrate model to simulate states
simulateSeries <- function(times = seq(0, 30, by = 2), sigma = 0) {

  m0 <- runif(n, log(10), 3) * c(2, 1) # Initial state matrix.

  Sim <- ode(m0, times, calcModel, par)
  Sim[, 2:(n+1)] <- matrix(rlnorm(Sim, log(Sim), sigma), nrow = nrow(Sim))[, 2:(n+1)]

  Sim[is.nan(Sim) | Sim < 0] <- 0

  return(Sim)
}

Sim <- simulateSeries(sigma = 0.02)
matplot(Sim[, 1], Sim[, -1], type = "b", ylab="N") # log='y'

#### Get multiple time series


## Multiple simulations
##
simulateMultipleSeries <- function(n_series = 10, n_times = 10, sigma = 0.1, format = c("long", "wide", "list")) {
  Sims <- replicate(n_series,
                    simulateSeries(0:(n_times-1), sigma = sigma),
                    simplify = F)

  if (match.arg(format) %in% c("wide", "long")) {
    Sims <- cbind(do.call(rbind, Sims), series = rep(1:n_series, each = n_times))
  }

  if (match.arg(format) == "long") {
    Sims <- tidyr::pivot_longer(as.data.frame(Sims),
                                cols = all_of(paste(1:n)),
                                names_to = "pop",
                                values_to = "abundance") %>%
      mutate(species = rep(1:n_species, each = 2, length.out = nrow(.)),
             stage = rep(c("j", "a"), length.out = nrow(.))
      )
  }
  return(Sims)
}


S <- simulateMultipleSeries(format = "long", n_times = 40, sigma = 0.05)

ggplot(S,
       mapping = aes(x = time, y = abundance, color = stage, group = interaction(stage, series))) +
  geom_line() +
  facet_wrap(facets = c("species"))



# Stan model --------------------------------------------------------------

#### Get a list of lists of data, where each second-level list belongs to one time series
## The lists will be unpacked as arrays of doubles, matrices, vectors etc. in stan
getStanData <- function(sims){

  ## replace NaNs with 0
  sims <- lapply(sims, function(s) { s[is.nan(s)] <- 0; return(s)})
  n_pops <- sapply(sims, function(s) ncol(s) - 1)[[1]]
  isjuvenile <- rep(c(T, F), length.out = n_pops)

  stanlist <- list(
    N_series = length(sims),

    ## For now, n_obs, and n_species have to be the same over all series. Thus, only first element extraction
    N_obs = sapply(sims, function(s) nrow(s)-1)[[1]], # observations, other than at t_0
    N_pops = n_pops,
    N_species = n_pops/2,

    time_init = sapply(sims, function(s) s[1, "time"]), # stan expects a one-dimensional array
    times = t(sapply(sims, function(s) s[2:nrow(s), "time"])), # stan expects a two-dimensional array real times[N_series,N_obs]; // arrays are row major!
    y_init = t(sapply(sims, function(s) s[1, -1])), # stan expects two-dimensional array: real y_init[N_series,N_species];
    Y = aperm(sapply(sims, function(s) s[-1, -1], simplify = "array"), c(3, 1, 2)) # stan expects two-dimensional array of vectors: vector[N_species] Y[N_series,N_obs] ## provided indexing order is opposite in R and stan!
  )
  return(c(stanlist, layoutMatrix(n_species = n_pops/2)))
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

    sigma = 2 # rep(0.1, data$N_species)
  )
}


getInitsTrue <- function() {
  list(
    g = A[I$I_g[,2:3]], # (+)
    s = A[I$I_s[,2:3]], # (-)

    A_j = A[I$Map_j[,3:4]],
    A_a = A[I$Map_a[,3:4]],

    r = f[I$isjuvenile],
    o = f[!I$isjuvenile],

    state_init =  data$y_init, # stan expects vector<lower=0>[N_species] state_init[N_series]

    sigma = 0.1 # rep(0.1, data$N_species)
  )
}


sims <- simulateMultipleSeries(format = "list", sigma = sigma)
data <- getStanData(sims)


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
                    output_dir = "Fits.nosync",
                    init = getInits,
                    iter_warmup = 200, iter_sampling = 500,
                    chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))

fit$summary(variables = 'state_init')
data$y_init

# shinystanfit <- rstan::read_stan_csv(fit$output_files())
# shinystan::launch_shinystan(shinystanfit)

