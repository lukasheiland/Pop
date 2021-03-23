# Library -----------------------------------------------------------------
library(here)

library(tidyverse)
library(ggformula)
library(magrittr)
library(glue)

library(deSolve)
library(cmdstanr)
# install_cmdstan(cores = 3)


# Orientation -------------------------------------------------------------

modelname <- "alphamatrix-threestage"

setwd(here())
modeldir <- dir(pattern = glue("^(Model).*{modelname}$"))
modelpath <- file.path(modeldir, glue('Model_{modelname}.stan'))


# Distribution ------------------------------------------------------------

rnbinom2 <- function(n, mu, phi) {
  rnbinom(n, mu = mu, size = phi)
}

dnbinom2 <- function(x, mu, phi) {
  dnbinom(x, mu = mu, size = phi)
}


# Multi-species matrix model -------------------------------------------------------------------


## Explicit third state version
calcModel <- function(t,
                      state, # A vector of species states.
                      pars){
  
  
  ## Just unpacking for readability.
  Alpha_j <- pars$Alpha_j
  Alpha_ab <- pars$Alpha_ab
  
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  g <- pars$g # Length n_species vector of transition rates.
  h <- pars$h # Length n_species vector of transition rates.
  
  n <- length(s)
  whichstate <- rep(1:3, each = n)
  
  J <- state[whichstate == 1]
  A <- state[whichstate == 2]
  B <- state[whichstate == 3]
  
  
  ## Here comes the model.
  AB <- A+B
  dJ <- r - (Alpha_j %*% J + s*sum(AB))*J - g*J
  dA <- g*J - (Alpha_ab %*% AB)*A - h*A
  dB <- h*A - (Alpha_ab %*% AB)*B

  return(list(c(dJ, dA, dB)))
}



# Simulation -------------------------------------------------------------------

#### Get environmental variables
simulateEnv <- function(n_env, n_locs){
  n <- n_env * n_locs
  Env <- matrix(runif(n, -1, 1), nrow = n_locs, ncol = n_env)
  colnames(Env) <- paste("env", 1:n_env, sep = "_")
  return(Env)
}

#### Return an interaction matrix.
generateAlpha <- function(n_species, negmu, phi, factor_diag, seed = 1){
  n <- n_species^2
  Alpha <- matrix(rlnorm(n, log(negmu), phi), nrow = n_species)
  diag(Alpha) <- diag(Alpha) * factor_diag
  return(Alpha)
}

generateParameters <- function(seed = 1,
                               n_species = 4, n_locs = 40, n_plotsperloc = 3, n_times = 10, n_env = 2,
                               
                               m_alpha = 0.2, phi_alpha = 0.1,
                               factor_diag_Alpha = 1.2, # factor for diagonal compared to others
                               factor_AB_Alpha = 0.5, # factor for adult matrix compared to juveniles)
                               
                               phi_obs = 100,
                               sigma_par = 0.1,
                               
                               ...
                               )
  {
  
  set.seed(seed)
  
  ## Beta_* are effect of environment on species' parameters matrix[n_env, n_species]
  ## These are on the log scale!
  n_beta <- 1 + ncol(poly(matrix(rnorm(100*n_env), 100, n_env), 2)) # + intercept
  
  Beta_r <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  r_log <- rnorm(n_species, 5, 0.1)
  Beta_r[1,] <- r_log
  Beta_r[3,] <- rnorm(n_species, -1, 0.2)
  
  
  Beta_s <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  s_log <- rnorm(n_species, -3, 0.1)
  Beta_s[1,] <- s_log
  Beta_s[3,] <- rnorm(n_species, -2, 0.2)
  
  
  Beta_g <- matrix(rnorm(n_beta*n_species, 1, 0.5), n_beta, n_species)
  g_log <- rnorm(n_species, -1.7, 0.1)
  Beta_g[1,] <- g_log
  Beta_g[3,] <- rnorm(n_species, -2, 0.2)
  
  
  h_log <- rnorm(n_species, -2, 0.5)
  
  Alpha_j <-  generateAlpha(n_species, m_alpha, phi_alpha, factor_diag_Alpha, seed = seed)
  Alpha_ab <- factor_AB_Alpha * generateAlpha(n_species, m_alpha, phi_alpha, factor_diag_Alpha, seed = seed+1)

  pars <- c(as.list(environment()), list(...))
  return(pars)
}

#### Transform parameters
### 1. Generate environmen parameters for the population model from expected values on the log scale
## three parameters dependent on environment: r, s, g.
## Beta_* are effect of environment on species' parameters matrix[1 + n_env, n_species]
### 2. other link functions

transformParameters <- function(pars, Env, returndf = F) {
  
  ## 1. Envrionmentally-dependent transformations
  M <- model.matrix(~ poly(Env, 2))
  R_log <- M %*% pars$Beta_r
  S_log <- M %*% pars$Beta_s
  G_log <- M %*% pars$Beta_g
  
  ## 2. Environmentally-independent transformation
  h <- t(replicate(pars$n_locs, rlnorm(pars$h_log, pars$h_log, pars$sigma_par)))
  
  transpars <- list(r = exp(R_log), s = exp(S_log), g = exp(G_log), h = h)
  
  if (returndf) {
    transpars <- lapply(transpars, as.data.frame)
    transpars <- do.call(cbind, transpars)
    transpars <- pivot_longer(cbind(transpars, Env),
                         cols = starts_with(c("r", "s", "g", "h")),
                         names_sep = "\\.",
                         names_to = c("parameter", "species"),
                         values_to = c("q")) %>%
      mutate(species = as.integer(as.factor(species)))
  } else {
    transpars <- c(pars, transpars)
  }
  
  return(transpars)
}



#### Simulate vector of initial states
generateInitialState <- function(n_species, n_stages = 3) {
  rnbinom2(n_species*n_stages, 8, 0.8) * rep(n_stages:1, each = n_species) # Initial state matrix.
}


#### Integrate model to simulate states
## simpars is a subset of pars
simulateOneSeries <- function(state_init, times, simpars, phi_obs = 10,
                              internaltimes = rep(F, length(times)), obserror = T) {

  ## ODE solving
  ## Different in ode interfaces in deSolve and stan! deSolve::ode() returns state at t0, stan only at the solution times.
  n_pops <- length(state_init)
  
  Sim <- ode(state_init, times, calcModel, simpars)

  if (obserror) { Sim[, 2:(n_pops+1)] <- matrix(rnbinom2(Sim, Sim, phi_obs), nrow = nrow(Sim))[, 2:(n_pops+1)] }
  Sim <- Sim[!internaltimes,]
  return(Sim)
}


#### Multiple simulations
simulateMultipleSeriesInEnv <- function(pars,
                                        Env,
                                        seed = 1,
                                        obserror = T,
                                        format = c("long", "wide", "locs", "init", "y0", "y", "list", "standatalist")) {
  
  set.seed(seed)
  n_times_intern <- pars$n_times + 1 ## internally there are more times, to ditch the first n_times from output
  areinternaltimes <- c(rep(T, 1), rep(F, pars$n_times))

  times <- seq(0.1, 3, length.out = n_times_intern)
  n_stages <- 3

  pars <- transformParameters(pars, Env)

  simulateOneSeriesInEnv <- function(r_loc, s_loc, g_loc, h_loc) {
    simpars <- within(pars, {r <- r_loc; s <- s_loc; g <- g_loc; h <- h_loc}) # list of dependent and independent parameters for one environment

    m0 <- generateInitialState(n_species = pars$n_species)
    seriesatloc <- replicate(pars$n_plotsperloc, simulateOneSeries(m0, times, simpars, pars$phi_obs, areinternaltimes, obserror = obserror), simplify = F)
    names(seriesatloc) <- 1:pars$n_plotsperloc
    return(seriesatloc)
  }

  ## CASE "list" etc.
  ## Sims is a: list[n_loc] <- list[n_plotsperloc] <- matrix[n_times, n_species]
  Sims <- mapply(simulateOneSeriesInEnv,
                 as.data.frame(t(pars$r)),
                 as.data.frame(t(pars$s)),
                 as.data.frame(t(pars$g)),
                 as.data.frame(t(pars$h)),
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
      mutate(stage = rep(c("j", "a", "b"), each = pars$n_species, length.out = nrow(.)))
    
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

    ## Format: [N_locs] — locations
    Sims_locs <- Sims_reobs %>% # reobs to  count all times
      group_by(loc, across(starts_with("env"))) %>%
      ## assumes completion within locations!
      summarize(n_species = n_distinct(species),
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
    X <- model.matrix(~ poly(Env, 2))
    N_init <- nrow(Sims_init)
    N_yhat <- nrow(Sims_yhat)
    time_init <-  Sims_y0[c("loc", "time")] %>% unique() %$% time # N_locs!
    

    vrep <- function(...) c( Vectorize(rep.int)(...) ) # :)))))
    
    Sims <- list(
      N_locs = nrow(Sims_locs),
      N_init = N_init,
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
      
      species = Sims_species$species,
      time_init = time_init,
      times = Sims_reobs$time,
      
      X = X,
      
      y0 = Sims_y0$abundance,
      y = Sims_reobs$abundance

    )
  }

  attr(Sims, "pars") <- pars
  return(Sims)
}




### Demo simulation: one time series ---------------------------------------------------------------
formals(generateParameters)
pars1 <- generateParameters(seed = 5)
initialstate1 <- generateInitialState(n_species = pars1$n_species)
pars1 <- within(pars1, {r = exp(pars1$r_log); s = exp(pars1$s_log); g = exp(pars1$g_log); h = exp(pars1$h_log)})
times1 <- seq(0, 3, length.out = 40)
Sim1 <- simulateOneSeries(initialstate1, times = times1, simpars = pars1, phi_obs = 10, # pars1$phi_obs
                          rep(F, length(times1)), obserror = F)
matplot(Sim1[, 1], Sim1[, -1], type = "b", ylab="N",
        pch = as.character(rep(1:3, each = pars1$n_species)),
        col = 1:pars1$n_species) # log='y'


### Demo simulations: multiple time series ---------------------------------------------------------------
pars_demo <- generateParameters(n_times = 50, n_species = 4, n_locs = 10, phi_obs = 10000)
E_demo <- simulateEnv(pars_demo$n_env, pars_demo$n_locs)
P_env_demo <- transformParameters(pars_demo, E_demo, returndf = T)
S_demo <- simulateMultipleSeriesInEnv(pars_demo, E_demo, obserror = F)

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

P_env_demo %>%
  filter(parameter %in% c("s")) %>%
  gf_point(q ~ env_1 | species) %>%
  gf_smooth()

ggplot(filter(P_env_demo, parameter == "r"),
       mapping = aes(x = env_1, y = q, color = species, group = species)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(facets = c("species"))






# Fit stan model  --------------------------------------------------------------

#### Simulate stan model data --------------------------------------------------------------
fitseed <- 1

pars <- generateParameters(n_times = 10, seed = fitseed)
Env <- simulateEnv(n_env = pars$n_env, n_locs = pars$n_locs)
data <- simulateMultipleSeriesInEnv(pars, Env, seed = fitseed, format = "standata")
Data_long <- simulateMultipleSeriesInEnv(pars, Env, seed = fitseed, format = "long", obserror = T) %>%
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
    
    state_init = aggregate(data$y0, by = list(data$rep_init2y0), mean)$x +
      rlnorm(data$N_y0), # to void zero inits
    
    ## Version with random rates.
    h = truepars$h,
    h_log = truepars$h_log,
    
    
    Alpha_j = truepars$Alpha_j,
    Alpha_ab = truepars$Alpha_ab,
    
    theta = 0.5,
    
    phi_obs = truepars$phi_obs,
    sigma_par = truepars$sigma_par
  )
  
  return(inits)
}


#### Calculate initial values based on data
## NOT YET WORKING!
getInits <- function() {
  
  n_s <- data$N_species
  n_ts <- data$N_totalspecies
  n_i <- data$N_inits
  n_l <- data$N_locs
  
  ## Preparing the rate matrices
  Beta_s <- Beta_g <- Beta_r <- matrix(rnorm(data$N_beta*n_ts, 0, 0.00001), data$N_beta, n_ts); # zero effect of environment
  # Beta_s[1,] <- 0
  
  #### Picking a solution for parameters with the given states, and
  ## given interactions s == 0, and Alpha_* == 0:
  
  #### h ----------
  ## h = dB/A
  dB <- filter(Data_long, stage == "b") %>%
    group_by(species, loc) %>%
    summarize(dB = first(abundance_loc[ist1] - abundance_loc[ist0])) %>%
    pivot_wider(id_cols = "loc", names_from = "species", values_from = "dB") %>%
    select(-loc) %>%
    as.matrix()
  A0 <- filter(Data_long, stage == "a") %>%
    group_by(species, loc) %>%
    summarize(A0 = first(abundance_loc[ist0])) %>%
    pivot_wider(id_cols = "loc", names_from = "species", values_from = "A0") %>%
    select(-loc) %>%
    as.matrix()
  h <- dB/A0
  h <- h + rlnorm(h, 10**-16)
  h[h<0 | is.nan(h)] <- 0
  h_log = log(apply(h, 2, mean))
  
  #### g ----------
  ## g = (dA + h*A)/J
  J0 <- filter(Data_long, stage == "j") %>%
    group_by(species, loc) %>%
    summarize(J0 = first(abundance_loc[ist0])) %>%
    pivot_wider(id_cols = "loc", names_from = "species", values_from = "J0") %>%
    select(-loc) %>%
    as.matrix()
  dA <- filter(Data_long, stage == "a") %>%
    group_by(species, loc) %>%
    summarize(dA = first(abundance_loc[ist1] - abundance_loc[ist0])) %>%
    pivot_wider(id_cols = "loc", names_from = "species", values_from = "dA") %>%
    select(-loc) %>%
    as.matrix()
  g <- (dA + h*A0)/J0
  g[is.nan(g)] <- 0
  Beta_g[1,] <- apply(g, 2, mean)
  
  #### r ----------
  ## r = dJ + g*J
  dJ <-  filter(Data_long, stage == "j") %>%
    group_by(species, loc) %>%
    summarize(dJ = first(abundance_loc[ist1] - abundance_loc[ist0])) %>%
    pivot_wider(id_cols = "loc", names_from = "species", values_from = "dJ") %>%
    select(-loc) %>%
    as.matrix()
  r <- dJ + g*J0
  r[is.nan(r)] <- 0
  Beta_r[1,] <- apply(r, 2, mean)
  
  
  inits <- list(
    Beta_s = Beta_s,
    Beta_g = Beta_g,
    Beta_r = Beta_r,
    
    state_init = aggregate(data$y0, by = list(data$rep_init2y0), mean)$x +
      rlnorm(data$N_y0), # to void zero inits
    
    ## Version with random rates.
    # r = matrix(rlnorm(n_l*n_s, log(0.5), 0.1), nrow = n_l, ncol = n_s),
    # s = matrix(rlnorm(n_l*n_s, log(0.5), 0.1), nrow = n_l, ncol = n_s),
    # g = matrix(rlnorm(n_l*n_s, log(0.5), 0.1), nrow = n_l, ncol = n_s),
    h = h,
    h_log = h_log,
    
    
    Alpha_j = matrix(rlnorm(n_ts*n_ts, log(0), 0.00001), nrow = n_ts, ncol = n_ts),
    Alpha_ab = matrix(rlnorm(n_ts*n_ts, log(0), 0.00001), nrow = n_ts, ncol = n_ts),
    
    theta = 0.5,
    
    phi_obs = 0.001,
    sigma_par = 10
    )
  return(inits)
}

compareParameters <- function(parname = c(sim = "h",
                              stan = "h"),
                              data,
                              stanfit) {
  
  simpar <- attr(data, "pars")[[parname[1]]]
  simpar <- c(simpar)
  
  stanpar <- stanfit$summary(variables = parname[2])$mean
  
  cbind(simpar, stanpar)
}


#### Compile model. --------------------------------------------------------------
# model$check_syntax()

model <- cmdstan_model(modelpath)
# model$exe_file() ## path to compiled model


#### Variational Bayes sampling. --------------------------------------------------------------

fit_var <- model$variational(data = data,
                         output_dir = "Fits.nosync",
                         init = getTrueInits,
                         iter = 2*10**2)

# fit_var$init()

# fit_var$summary()

Estimates <- fit_var$summary(variables = c("sigma_par", "phi_obs", "state_init", "r_log", "h", "Alpha_j")) # "o", "state_init"

## true parameters vs. estimates
Comp_h <- compareParameters(c("r", "r_log"), data, fit_var)
plot(exp(Comp_h[,2]), Comp_h[,1])

Comp_R <- compareParameters(c("Beta_r", "Beta_r"), data, fit_var)
plot(Comp_R[,2], Comp_R[,1])


## parameter vs. env
r <- exp(fit_var$summary(variables = "r_log")$mean)
plot(Env[,1], r[seq(1, length(r), by = 4)])

Beta_r <- fit_var$summary(variables = "Beta_r")$mean
r_hat <- data$X %*% matrix(Beta_r, data$N_beta, data$N_totalspecies)
points(Env[,1], exp(r_hat[,1]), col = "red")

r_hat_true <- data$X %*% attr(data, "pars")$Beta_r
points(Env[,1], exp(r_hat_true[,1]), col = "blue")


#### NUTS sampling. --------------------------------------------------------------


n_chains <- 3
fit <- model$sample(data = data,
                    output_dir = "Fits.nosync",
                    init = getTrueInits,
                    iter_warmup = 200, iter_sampling = 500,
                    chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))

fit <- model$sample(data = data,
                    fixed_param = T,
                    output_dir = "Fits.nosync",
                    init = getTrueInits,
                    iter_warmup = 10, iter_sampling = 100,
                    chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))

draws <- read_cmdstan_csv(fit$output_files())

fit$init()
fit$summary(variables = 'state_init')
truepars <- attr(data, "pars")
log(truepars$r)
data$y_init

# shinystanfit <- rstan::read_stan_csv(fit$output_files())
# shinystan::launch_shinystan(shinystanfit)

