# Library -----------------------------------------------------------------
library(here)

library(tidyverse)
library(ggformula)
library(magrittr)
library(glue)

library(tictoc)

library(deSolve)
library(cmdstanr)
# install_cmdstan(cores = 3)


# Orientation -------------------------------------------------------------

modelname <- "m"

setwd(here())
modeldir <- dir(pattern = glue("^(Model).*{modelname}$"))
modelpath <- file.path(modeldir, glue('Model_{modelname}.stan'))
modelpath_continuous <- file.path(modeldir, glue('Model_{modelname}_continuous.stan'))



# Distribution ------------------------------------------------------------

rnbinom2 <- function(n, mu, phi) {
  rnbinom(n, mu = mu, size = phi)
}

dnbinom2 <- function(x, mu, phi) {
  dnbinom(x, mu = mu, size = phi)
}

dgamma2 <- function(x, mean, shape) {
  dgamma(x, shape = shape, rate = shape/mean)
}

rgamma2 <- function(n, mean, shape) {
  rgamma(n, shape = shape, rate = shape/mean)
}

######################################################################################
# Model formulations  ----------------------------------------------------------------
######################################################################################

## Continuous model -------------------------------------------------------------------

#### ODE continuous model version
## returns derivatives ("change rates")

calcModel <- function(t,
                      state, # A vector of species states.
                      pars){
  
  ## Just unpacking for readability.
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  g <- pars$g # Length n_species vector of transition rates.
  m_j <- pars$m_j # Length n_species vector of mortalities
  h <- pars$h # Length n_species vector of transition rates.
  c_j <- pars$c_j # Length n_species vector
  c_ab <- pars$c_ab # Length n_species vector
  
  n <- length(s)
  whichstate <- rep(1:3, each = n)
  
  J <- state[whichstate == 1]
  A <- state[whichstate == 2]
  B <- state[whichstate == 3]
  AB <- sum(A + B) 
  

  dJ <- r - (c_j*sum(J)  + s*AB + m_j)*J - g*J
  dA <- g*J - (c_ab * AB + h)*A
  dB <- h*A - (c_ab * AB)*B
  
  return(list(c(dJ, dA, dB)))
}



## Discrete model -------------------------------------------------------------------

#### Euler method, discretizing the above continuous model
## returns states for times
calcDiscreteModel <- function(times,
                              initialstate, # A vector of species states.
                              pars,
                              stepspertime = NULL # internal number of discrete time steps within a unit of time
                              ) { 
  
  
  ## Just unpacking for readability.
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  g <- pars$g # Length n_species vector of transition rates.
  h <- pars$h # Length n_species vector of transition rates.
  m_j <- pars$m_j # Length n_species vector of mortalities
  c_j <- pars$c_j # Length n_species vector
  c_ab <- pars$c_ab # Length n_species vector
  sigma <- pars$sigma_process
  
  ## Set the count variables
  n <- length(s)
  if(is.null(stepspertime)) {
    times_intern <- times
  } else {
    timediffmin <- min(times[2:(length(times))] - times[1:(length(times)-1)])
    times_intern <- seq(min(times), max(times), by = timediffmin/stepspertime)
    times_intern <- sort(unique(c(times, times_intern))) # Wow, what an Aufwand
  }
  n_times <- length(times_intern)
  dt <- times_intern[2:n_times] - times_intern[1:(n_times-1)]
  
  # Prepare a state matrix
  whichstate <- rep(1:3, each = n)
  State <- matrix(rep(initialstate, times = n_times), nrow = n_times, byrow = T)
  
  ## Here comes the model.
  for (i in 2:n_times) {
    J <- State[i-1,whichstate == 1]
    A <- State[i-1,whichstate == 2]
    B <- State[i-1,whichstate == 3]
    AB <- sum(A + B)
    
    e_j <- rnorm(s, 0, sigma[1])
    J_i <- J + (r - (c_j*sum(J)  + s*AB + m_j + g)*J + e_j) * dt[i-1]
    
    e_a <- rnorm(s, 0, sigma[2])
    A_i <- A + (g*J - (c_ab * AB + h)*A + e_a) * dt[i-1]
    
    e_b <- rnorm(s, 0, sigma[3])
    B_i <- B + (h*A - (c_ab * AB)*B + e_b) * dt[i-1]
    
    state_i <- c(J_i, A_i, B_i)
    State[i, ] <- replace(state_i, state_i<0, 0) # it is enough to replace the state here, in the next iteration only the state at t-1 is taken anyway
  }
  returnwhichtimes <- which(times_intern %in% times)
  return(State[returnwhichtimes,])
}



# Simulation -------------------------------------------------------------------

#### Get environmental variables
simulateEnv <- function(n_env, n_locs){
  n <- n_env * n_locs
  Env <- matrix(runif(n, -1, 1), nrow = n_locs, ncol = n_env)
  colnames(Env) <- paste("env", 1:n_env, sep = "_")
  return(Env)
}


generateParameters <- function(seed = 1,
                               n_species = 4, n_locs = 40, n_plotsperloc = 3, n_times = 10, n_env = 2,
                               sigma_process = c(10, 6, 5), shape_par = 10, phi_obs = c(500, 1000),
                               ...
                               )
  {
  
  set.seed(seed)
  
  ## Beta_* are effect of environment on species' parameters matrix[n_env, n_species]
  ## These are on the log scale!
  n_beta <- 1 + 2 * n_env # degree 2 polynomial + intercept
  
  Beta_r <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  r_log <- rnorm(n_species, 3.8, 0.1)
  Beta_r[1,] <- r_log
  Beta_r[3,] <- rnorm(n_species, -1, 0.2)
  
  
  Beta_s <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  s_log <- rnorm(n_species, -2.2, 0.1)
  Beta_s[1,] <- s_log
  Beta_s[3,] <- rnorm(n_species, -2, 0.2)
  
  Beta_g <- matrix(rnorm(n_beta*n_species, 1, 0.5), n_beta, n_species)
  g_log <- rnorm(n_species, -1.3, 0.1)
  Beta_g[1,] <- g_log
  Beta_g[3,] <- rnorm(n_species, -2, 0.2)
  
  Beta_m_j <- matrix(rnorm(n_beta*n_species, -4, 0.5), n_beta, n_species)
  m_j_log <- rnorm(n_species, -5, 0.1)
  Beta_m_j[1,] <- m_j_log
  Beta_m_j[3,] <- rnorm(n_species, -3, 0.2)
  
  ## env-independent parameters
  h_log <- rnorm(n_species, -1.7, 0.5)
  c_j_log <- rnorm(n_species, -4, 0.5)
  c_ab_log <- rnorm(n_species, -3, 0.5)
  

  pars <- c(as.list(environment()), list(...)) ## Carries along all parameter values within the named list.
  return(pars)
}

#### Transform parameters
### 1. Generate environmen parameters for the population model from expected values on the log scale
## three parameters dependent on environment: r, s, g.
## Beta_* are effect of environment on species' parameters matrix[1 + n_env, n_species]
### 2. other link functions

transformParameters <- function(pars, Env, returndf = F) {
  
  Env <- as.data.frame(Env)
  
  ## 1. Envrionmentally-dependent transformations
  polyformula <- as.formula(paste("~", paste("poly(", colnames(Env), ", 2)", collapse = "+")))
  M <- model.matrix(polyformula, data = as.data.frame(Env))
  R_log <- M %*% pars$Beta_r
  S_log <- M %*% pars$Beta_s
  G_log <- M %*% pars$Beta_g
  M_j_log <- M %*% pars$Beta_m_j
  
  ## 2. Environmentally-independent transformation
  h <- t(replicate(pars$n_locs, rgamma2(pars$h_log, exp(pars$h_log), pars$shape_par)))
  c_j <- t(replicate(pars$n_locs, rgamma2(pars$c_j_log, exp(pars$c_j_log), pars$shape_par)))
  c_ab <- t(replicate(pars$n_locs, rgamma2(pars$c_ab_log, exp(pars$c_ab_log), pars$shape_par)))
  
  
  transpars <- list(r = exp(R_log), s = exp(S_log), g = exp(G_log), m_j = exp(M_j_log),
                    h = h, c_j = c_j, c_ab = c_ab)
  logpars <- list(R_log = R_log, S_log = S_log, G_log = G_log, M_j_log = M_j_log)
  
  if (returndf) {
    transpars <- lapply(transpars, as.data.frame)
    transpars <- do.call(cbind, transpars)
    transpars <- pivot_longer(cbind(transpars, Env),
                         cols = starts_with(c("r", "s", "g", "h", "c_", "m_")),
                         names_sep = "\\.",
                         names_to = c("parameter", "species"),
                         values_to = c("q")) %>%
      mutate(species = as.integer(as.factor(species)))
  } else {
    transpars <- c(pars, transpars, logpars)
  }
  
  return(transpars)
}



#### Simulate vector of initial states
generateInitialState <- function(n_species, n_stages = 3, generatelog = T) {
  s <- rnbinom2(n_species*n_stages, 5, 3) * rep(n_stages:1, each = n_species) # Initial state matrix.
  if(generatelog) return(log1p(s)) else return(s)
}


#### Integrate model to simulate states
## simpars is a subset of pars
simulateOneSeries <- function(state_init, times, simpars, phi_obs = c(500, 5000),
                              internaltimes = rep(F, length(times)),
                              discrete = F, processerror = T, obserror = T, expinit = T, stepspertime = NULL) {

  ## ODE solving
  ## Different in ode interfaces in deSolve and stan! deSolve::ode() returns state at t0, stan only at the solution times.
  n_pops <- length(state_init)
  if(expinit) state_init <- exp(state_init)
  if(discrete) {
    if(!processerror) simpars$sigma_process <- c(0, 0, 0)
    Sim <- calcDiscreteModel(initialstate = state_init, times = times, pars = simpars, stepspertime = stepspertime)
    ## Make Sim look like from ode
    Sim <- set_colnames(cbind(times, Sim), c("time", paste(1:ncol(Sim))))
  } else {
    Sim <- ode(state_init, times, calcModel, simpars)
  }

  if (obserror) {
    ## phi for juveniles
    Sim[, 2:(n_pops/3+1)] <- matrix(rnbinom2(Sim, Sim, phi_obs[1]), nrow = nrow(Sim))[,  2:(n_pops/3+1)]
    Sim[, (n_pops/3+1):(n_pops+1)] <- matrix(rnbinom2(Sim, Sim, phi_obs[2]), nrow = nrow(Sim))[, (n_pops/3+1):(n_pops+1)]
  }
  Sim <- Sim[!internaltimes,]
  return(Sim)
}


#### Multiple simulations
simulateMultipleSeriesInEnv <- function(pars,
                                        Env,
                                        seed = 1,
                                        discrete = F,
                                        processerror = T, obserror = T, stepspertime = NULL,
                                        format = c("long", "wide", "locs", "init", "y0", "y", "list", "standatalist")) {
  
  set.seed(seed)
  n_cutoffstarttimes <- 1
  n_times_intern <- pars$n_times + n_cutoffstarttimes ## internally there are more times, to ditch the first n_times from output
  areinternaltimes <- c(rep(T, n_cutoffstarttimes), rep(F, pars$n_times))

  times <- round(seq(0, 3, length.out = n_times_intern))
  n_stages <- 3

  pars <- transformParameters(pars, Env)

  simulateOneSeriesInEnv <- function(r_loc, s_loc, g_loc, m_j_loc,
                                     h_loc, c_j_loc, c_ab_loc) {
    simpars <- within(pars, {r <- r_loc; s <- s_loc; g <- g_loc; m_j <- m_j_loc;
                             h <- h_loc; c_j <- c_j_loc; c_ab <- c_ab_loc;}) # list of dependent and independent parameters for one environment

    m0 <- generateInitialState(n_species = pars$n_species)
    seriesatloc <- replicate(pars$n_plotsperloc,
                             simulateOneSeries(m0, times, simpars, pars$phi_obs, areinternaltimes,
                                               discrete = discrete, processerror = processerror, obserror = obserror, expinit = T), simplify = F)
    names(seriesatloc) <- 1:pars$n_plotsperloc
    return(seriesatloc)
  }

  ## CASE "list" etc.
  ## Sims is a: list[n_loc] <- list[n_plotsperloc] <- matrix[n_times, n_species]
  Sims <- mapply(simulateOneSeriesInEnv,
                 as.data.frame(t(pars$r)),
                 as.data.frame(t(pars$s)),
                 as.data.frame(t(pars$g)),
                 as.data.frame(t(pars$m_j)),

                 as.data.frame(t(pars$h)),
                 as.data.frame(t(pars$c_j)),
                 as.data.frame(t(pars$c_ab)),
                 
                 SIMPLIFY = F)
  

  ## CASE: wide etc.
  if (match.arg(format) %in% c("long", "wide", "locs", "init", "y0", "y", "standatalist")) {
    Sims <- lapply(Sims, function(l) do.call(rbind, l)) # bind the series within a location into one table
    Sims <- lapply(Sims, function(m) cbind(m, plot = rep(1:pars$n_plotsperloc, each = pars$n_times))) # attach an id to series per location. Each n_times, because simulateOneSeries() does only return n_times (n_times_intern/2) observations

    Sims <- do.call(rbind, Sims) # bind all lists of locations into one table
    Sims <- cbind(Sims, loc = rep(1:pars$n_locs, each = pars$n_plotsperloc * (pars$n_times))) # attach an id to series per location

    Sims <- cbind(Sims, Env[Sims[,"loc"],])
  }
  
  ## CASE: long etc.
  ## Here come all different long formats used in different hierarchical levels of the stan fit.
  
  
  if (match.arg(format) %in% c("long", "locs", "init", "y0", "y", "standatalist")) {
    Sims <- tidyr::pivot_longer(as.data.frame(Sims),
                                cols = all_of(paste(1:(pars$n_species*n_stages))),
                                names_to = "pop",
                                values_to = "abundance") %>%
      mutate(pop = as.integer(pop)) %>%
      mutate(species = rep(1:pars$n_species, times = n_stages, length.out = nrow(.))) %>%
      mutate(stage = rep(c("j", "a", "b"), each = pars$n_species, length.out = nrow(.))) %>%
      mutate(obsmethod = rep(c(1, 2, 2), each = pars$n_species, length.out = nrow(.)))
    
    ## Sims <- complete(Sims, group = nesting("loc", "time", species", "stage"), fill = 0) # is assumed!
    Sims <- group_by(Sims, loc) %>%
      mutate(isy0 = time == min(time)) %>%
      ungroup() %>%
      mutate(stage = factor(stage, levels = c("j", "a", "b"), ordered = T)) %>%
      arrange(loc, time, pop, plot)
    
    
    
    #### Data sets structure
    ## The most comprehensive data set is N_y (resp. N_y0) with grouping locations/resurveys/pops/plots.
    ## NOTE THAT plots HAS TO BE THE LAST GROUP IN SORTING
    ## Everything is subset from the master subsets with these groupings (*_reobs, *_y0) and thus consistently sorted.
    ## Ratio: "resurveys" includes "pops" due of the return structure in ode_*(); "pops" includes "plots" because of the loc == population assumption (randowm effect from location/resurvey/pop to .../plots).
    ## Factor "pops" is structured stages/species.
    
    
    ## Format: [N_y0] —   locations/pops(/stage/species)/plots
    Sims_y0 <- filter(Sims, isy0) %>% select(-isy0) %>%
      arrange(loc, time, pop, plot)
    
    ## Format: [N_y] — locations/resurveys/pops(/stage/species)/plots
    Sims_reobs <- filter(Sims, !isy0) %>% select(-isy0) %>%
      arrange(loc, time, pop, plot)
    
    ## Format: [N_init] — locations/pops
    Sims_init <- Sims_y0 %>%
      group_by(loc, pop, stage, species) %>%
      summarize(n_plots = n_distinct(plot), .groups = "drop")
    
    ## Format: [N_yhat] — locations/resurveys/pops
    Sims_yhat <- Sims_reobs %>%
      group_by(loc, time, pop, stage, species) %>%
      summarize(n_plots = n_distinct(plot), .groups = "drop")
    
    ## Format: [N_times] — locations/resurveys
    Sims_times <- Sims_reobs %>% # reobs to  count all times
      group_by(loc) %>%
      summarize(time = unique(time), .groups = "drop")

    ## Format: [N_locs] — locations
    Sims_locs <- Sims_reobs %>% # reobs to  count all times
      group_by(loc) %>%
      ## assumes completion within locations!
      summarize(time_max = max(time),
                n_species = n_distinct(species),
                n_plots = n_distinct(plot),
                n_pops = n_distinct(pop),
                n_reobs = n_distinct(time),
                n_yhat = n_distinct(interaction(pop, time)), .groups = "drop")

  }
  
  
  ## CASE: locs
  if (match.arg(format) == "locs") Sims <- Sims_locs
  
  ## CASE: init
  if (match.arg(format) == "init") Sims <- Sims_init
  
  ## CASE: y0
  if (match.arg(format) == c("y0")) Sims <- Sims_y0

  ## CASE: y
  if (match.arg(format) == "y") Sims <- Sims_reobs
  
  if (match.arg(format) == "standatalist") {
    Sims_species <- Sims_init[c("loc", "species")] %>% unique()
    polyformula <- as.formula(paste("~", paste("poly(", colnames(as.data.frame(Env)), ", 2)", collapse = "+")))
    X <- model.matrix(polyformula, data = as.data.frame(Env))
    N_init <- nrow(Sims_init)
    N_times <- nrow(Sims_times)
    N_yhat <- nrow(Sims_yhat)
    N_locs <-  nrow(Sims_locs)
    time_init <-  Sims$time[match(Sims_locs$loc, Sims$loc)] # N_locs!
    time_init <-  Sims$time[match(Sims_locs$loc, Sims$loc)] # N_locs!
    
  
    vrep <- function(...) unlist( c( Vectorize(rep.int)(...) ) ) # :)))))
    
    Sims <- list(
      N_locs = N_locs,
      N_init = N_init,
      N_times = N_times,
      N_yhat = N_yhat,
      N_y0 = nrow(Sims_y0),
      N_y = nrow(Sims_reobs),
      N_species = nrow(Sims_species), 
      
      
      N_totalspecies = length(unique(Sims_init$species)),
      
      N_beta = ncol(X),
      
      n_species = Sims_locs$n_species, 
      n_pops = Sims_locs$n_pops,
      n_reobs = Sims_locs$n_reobs,
      n_yhat = Sims_locs$n_yhat,
      
      ## repeat predictions on level "locations/pops" n_plots times to "locations/pops/plots"
      rep_init2y0 = vrep(1:N_init, Sims_init$n_plots),
      ## repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
      rep_yhat2y = vrep(1:N_yhat, Sims_yhat$n_plots),
      ## repeat locations on level "locations" n_reobs times to "locations/pops/resurveys"
      rep_locs2times = vrep(1:N_locs, Sims_locs$n_reobs),
      
      obsmethod_y0 = Sims_y0$obsmethod,
      obsmethod_y = Sims_reobs$obsmethod,
      
      species = Sims_species$species,
      stage = Sims_init$stage,
      pops = Sims_init$pop,
      time_init = time_init,
      time_max_data = Sims_locs$time_max,
      times_data = Sims_times$time,
      
      X = X,
      
      y0 = Sims_y0$abundance,
      y = Sims_reobs$abundance

    )
  }

  attr(Sims, "pars") <- pars
  return(Sims)
}



### Demo simulation: one time series ---------------------------------------------------------------

#### prerequisites for a single time series
formals(generateParameters)
pars1 <- generateParameters(seed = 1)
initialstate1 <- generateInitialState(n_species = pars1$n_species)
pars1 <- within(pars1, {
  r = exp(pars1$r_log); s = exp(pars1$s_log); g = exp(pars1$g_log); m_j = exp(pars1$m_j);
  h = exp(pars1$h_log); c_j = exp(pars1$c_j_log); c_ab = exp(pars1$c_ab_log);
  sigma_process = c(10, 5, 4);
})
times1 <- seq(0, 1, length.out = 20)

#### 1. Continuous
Sim1 <- simulateOneSeries(initialstate1, times = times1, simpars = pars1, phi_obs = c(10, 100), # pars1$phi_obs                       rep(F, length(times1)), obserror = F)
                          internaltimes = rep(F, length(times1)), discrete = F, obserror = F)

matplot(Sim1[-1, 1], Sim1[-1, -1], type = "b", ylab="N",
        pch = rep(c("J", "A", "B"), each = pars1$n_species),
        col = 1:pars1$n_species, xlab = "time") # log='y'

#### 2. Discrete
# matplot(calcDiscreteModel(times1, exp(initialstate), pars1), type = "b")

Sim1_discr <- simulateOneSeries(initialstate1, times = times1, simpars = pars1,
                                internaltimes = rep(F, length(times1)),
                                discrete = T, processerror = T, obserror = F,
                                stepspertime = NULL) # alternative: NULL. number of steps increases accuracy compared to ODE

matplot(Sim1_discr[-1, 1], Sim1_discr[-1, -1], type = "b", ylab="N",
        pch = rep(c("J", "A", "B"), each = pars1$n_species),
        col = 1:pars1$n_species, xlab = "time") # log='y'



### Demo simulations: multiple time series ---------------------------------------------------------------
pars_demo <- generateParameters(n_times = 50, n_species = 4, n_locs = 10, phi_obs = c(10000, 1000000))
E_demo <- simulateEnv(pars_demo$n_env, pars_demo$n_locs)
P_env_demo <- transformParameters(pars_demo, E_demo, returndf = T)
S_demo <- simulateMultipleSeriesInEnv(pars_demo, E_demo, discrete = F, obserror = F)

#### Plot time series
S_demo %>%
  # filter(stage %in% c("a")) %>%
  ggplot(mapping = aes(x = time, y = abundance, color = env_1, lty = stage, group = interaction(loc, stage, plot))) +
  geom_line() +
  facet_wrap(facets = c("species"))


#### Plot parameters
P_env_demo %>%
  filter(parameter %in% c("g")) %>%
  gf_point(q ~ env_1 | species) %>%
  gf_smooth()

ggplot(filter(P_env_demo, parameter == "r"),
       mapping = aes(x = env_1, y = q, color = species, group = species)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(facets = c("species"))

P_env_demo %>%
  filter(parameter %in% c("s")) %>%
  gf_point(q ~ env_1 | species) %>%
  gf_smooth()



# Fit stan model  --------------------------------------------------------------

#### Simulate stan model data --------------------------------------------------------------
fitseed <- 1

pars <- generateParameters(n_times = 3, seed = fitseed, n_locs = 50)
Env <- simulateEnv(n_env = pars$n_env, n_locs = pars$n_locs)
data <- simulateMultipleSeriesInEnv(pars, discrete = F,
                                    Env, seed = fitseed, format = "standata")
Data_long <- simulateMultipleSeriesInEnv(pars, Env,
                                         discrete = F, # use the continuous version!
                                         stepspertime = 10,
                                         seed = fitseed, format = "long", obserror = F) %>%
  mutate(sortid = paste(loc, plot, species, stage, sep = "_")) %>%
  arrange(sortid, time) %>%
  group_by(loc, species, stage, time) %>%
  mutate(abundance_loc = mean(abundance)) %>%
  group_by(loc, plot, species, stage) %>%
  mutate(ist0 = time == time[1], ist1 = time == time[2]) %>%
  ungroup()

####### Plot time series
Data_long %>%
  # filter(stage %in% c("a")) %>%
  ggplot(mapping = aes(x = time, y = abundance, color = env_1, lty = stage, group = interaction(loc, stage, plot))) +
  geom_line() +
  facet_wrap(facets = c("species"))


##### Start values ------------------------------------------------------------

getTrueInits <- function() {
  
  truepars <- attr(data, "pars")
  
  inits <- list(
    ## from global environment
    Beta_s = truepars$Beta_s,
    Beta_g = truepars$Beta_g,
    Beta_r = truepars$Beta_r,
    Beta_m_j = truepars$Beta_m_j,

    state_init_log = log1p(aggregate(data$y0, by = list(data$rep_init2y0), mean)$x +
      rlnorm(data$N_y0)), # to void zero inits
    
    ## Version with random rates.
    h        = truepars$h,
    c_j      = truepars$c_j,
    c_ab     = truepars$c_ab,
    h_log    = truepars$h_log,
    c_j_log  = truepars$c_j,
    c_ab_log = truepars$c_ab,
    
    sigma_process   = truepars$sigma_process,
    phi_obs   = truepars$phi_obs,
    shape_par = truepars$shape_par
  )
  
  return(inits)
}



#### Compile model. --------------------------------------------------------------

model <- cmdstan_model(modelpath)
# model <- cmdstan_model(modelpath, force_recompile = T)
# model <- cmdstan_model(modelpath_continuous)

compiledpath <-  model$exe_file()


#### Variational Bayes sampling. --------------------------------------------------------------
fit_var <- model$variational(data = data,
                             output_dir = "Fits.nosync",
                             init = getTrueInits,
                             eta = 1,
                             iter = 10**4) # convergence after 1500 iterations

# fit_var$summary()
fit_var$summary(variables = c("shape_par", "phi_obs", "h_log", "c_j", "state_init", "u")) %>% # "o", "state_init"
  View()


# Extract draws -----------------------------------------------------------
library(rstan)
library(bayesplot)


plotFitVsTrue <- function(parname = c(sim = "R_log",
                                      stan = "r_log"),
                              data,
                              rstandraws) {
  
  simpar <- attr(data, "pars")[[parname[1]]]
  # simpar <- c(simpar)
  
  stanpar <- rstan::extract(rstandraws, pars = parname[2])[[1]] %>% aperm(c(2,3,1)) %>% c()
  
  # cbind(simpar, stanpar
  plot(rep(c(simpar), length.out = length(stanpar)), stanpar, cex = 0.1, col = alpha("black", alpha = 0.03))
  abline(a = 0, b = 1, col = "green")
}



drawpath <- fit_var$output_files()

## alternatively
# drawpath <- file.info(list.files("Fits.nosync", full.names = T)) %>%
#   arrange(desc(mtime)) %>%
#   slice(1) %>%
#   rownames()

## get draws with cmdstanr
# draws_var <- read_cmdstan_csv(recentfitpath)

## get draws with rstan
draws_var <- rstan::read_stan_csv(fit_var$output_files())


## Fitted parameters vs. true
plotFitVsTrue(c("c_j", "c_j"), data, draws_var)
plotFitVsTrue(c("c_ab", "c_ab"), data, draws_var)

plotFitVsTrue(c("Beta_m_j", "Beta_m_j"), data, draws_var)

plotFitVsTrue(c("Beta_g", "Beta_g"), data, draws_var)
plotFitVsTrue(c("G_log", "g_log"), data, draws_var)

plotFitVsTrue(c("Beta_r", "Beta_r"), data, draws_var)
plotFitVsTrue(c("R_log", "r_log"), data, draws_var)

plotFitVsTrue(c("Beta_s", "Beta_s"), data, draws_var) #???
plotFitVsTrue(c("S_log", "s_log"), data, draws_var)



# Plot parameter vs. env -------------------------------------------------------

### predict from Beta
predictB <- function(B, x = seq(-1, 1, by = 0.01), n_beta = 2, species = 1) {
  E <- as.data.frame(replicate(n_beta, x))
  polyformula <- as.formula(paste("~", paste("poly(", colnames(E), ", 2)", collapse = "+")))
  X <- model.matrix(polyformula, data = E)
  exp(X %*% B)[,1]
}

poly2 <- function(x, b) {b[1] + x*b[2] + x^2*b[3]}
getVertex <- function(b) {-b[2]/(2*b[3])}


x <- Env[,1]
truepars <- attr(data, "pars")
Beta_r_true <- truepars$Beta_r
Beta_r_mean <- matrix(fit_var$summary(variables = "Beta_r")$mean, data$N_beta, data$N_totalspecies)
plot(x, predictB(Beta_r_mean, x))
points(x, predictB(Beta_r_true, x), col = "red")


Beta_g_true <- truepars$Beta_g
Beta_g_mean <- matrix(fit_var$summary(variables = "Beta_g")$mean, data$N_beta, data$N_totalspecies)
plot(x, predictB(Beta_g_true, x), col = "red")
plot(x, predictB(Beta_g_mean, x))



#### NUTS sampling. --------------------------------------------------------------

n_chains <- 3
fit <- model$sample(data = data,
                    output_dir = "Fits.nosync",
                    init = getTrueInits,
                    iter_warmup = 200, iter_sampling = 200,
                    chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))


draws <- read_cmdstan_csv(fit$output_files())
fit$summary(variables = 'c_j')


# shinystanfit <- rstan::read_stan_csv(fit$output_files())
# shinystan::launch_shinystan(shinystanfit)

