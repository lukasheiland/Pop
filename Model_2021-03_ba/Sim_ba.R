# Library -----------------------------------------------------------------
library(here)

library(tidyverse)
library(ggformula)
library(magrittr)
library(glue)

library(tictoc)

library(cmdstanr)
# install_cmdstan(cores = 3)
library(rstan)
library(bayesplot)


# Orientation -------------------------------------------------------------
setwd(here())


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
# Model formulation  ----------------------------------------------------------------
######################################################################################

## Discrete model -------------------------------------------------------------------

#### Returns log states for times -------------------
calcModel <- function(times,
                      initialstate_log, # A vector of species states.
                      pars # internal number of discrete time steps within a unit of time
) { 
  
  
  ## Just unpacking for readability. Order: 1. environmentally-dependent, 2. -independent
  
  g <- pars$g # Length n_species vector of transition rates.
  m_j <- pars$m_j # Length n_species vector of mortalities
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  
  b <- pars$b # Length n_species vector of basal area increment rates.
  c_j <- pars$c_j # Length n_species vector
  c_a <- pars$c_a # Length n_species vector
  c_b <- pars$c_b # Length n_species vector
  h <- pars$h # Length n_species vector of transition rates.
  m_a <- pars$m_a # Length n_species vector of mortalities
  
  dbh_lower_a <- pars$dbh_lower_a
  dbh_lower_b <- pars$dbh_lower_b
  
  sigma <- pars$sigma_process
  
  ## Set the count variables
  n <- length(r) # no species
  
  ## 
  times_intern <- 1:max(times)
  n_times <- length(times_intern)
  
  
  radius_a_upper <- dbh_lower_a/2
  radius_a_avg <- (dbh_lower_a + dbh_lower_b)/2/2 # [mm]
  ba_a_upper <-  pi * radius_a_upper^2 * 1e-6
  ba_a_avg <- pi * radius_a_avg^2 * 1e-6 # mm^2 to m^2
  
  # Prepare a state matrix
  whichstate <- rep(1:3, each = n)
  State_log <- matrix(rep(initialstate_log, times = n_times), nrow = n_times, byrow = T)
  
  ## Here comes the model.
  for (t in 2:n_times) {
    
    ## States at t-1: *State*_log, and *State*
    J_log <- State_log[t-1,whichstate == 1]
    A_log <- State_log[t-1,whichstate == 2]
    B_log <- State_log[t-1,whichstate == 3]
    
    J <- exp(J_log)
    A <- exp(A_log)
    B <- exp(B_log)
    
    ## The total basal area of big trees
    BA <- sum(A * ba_a_avg + B)
    
    # observation error with for stages/times/species
    u <- array(0, dim = c(3, n_times-1, length(s))) + rnorm(length(s)*3*(n_times-1), 0, sigma)
    
    ## Two kinds of processes acting ot State_log
    ## 1. All processes that add additively to State, are added within log(State)
    ## 2. All processes that act multiplicatively on the state are added in log space
    
    J_t_log <- log(J + r) - c_j*sum(J) - s*BA - m_j - g + u[1, t-1, ] # count of juveniles J
    
    A_t_log <- log(A + exp(J_log - g)) - c_a*BA - m_a - h + u[2, t-1, ] # count of small adults A
    A_ba <- exp(A_log - h)*ba_a_upper # Basal area of small adults A. Conversion by multiplication with basal area of State exit (based on upper dhh boundary of the class)
    
    B_t_log <- log(B + A_ba) + b - c_b*sum(B) + u[3, t-1, ] # basal area of big adults B
    ## b is the net basal area increment (including density-independent m) basically equivalent to a Ricker model, i.e. constant increment rate leading to exponential growth, negative density dependent limitation scaled with total BA.
    
    State_log[t, ] <- c(J_t_log, A_t_log, B_t_log)
  }
  
  whichtimes <- which(times_intern %in% times)
  return(State_log[whichtimes,])
}



######################################################################################
# Simulation functions ---------------------------------------------------------------
######################################################################################


#### Get environmental variables -------------------
simulateEnv <- function(n_env, n_locs){
  n <- n_env * n_locs
  Env <- matrix(runif(n, -1, 1), nrow = n_locs, ncol = n_env)
  colnames(Env) <- paste("env", 1:n_env, sep = "_")
  return(Env)
}


#### Transform parameters -------------------
## Accepts and returns a list of parameters. No parameter is lost.
##
### 1. Generate environment parameters for the population model from expected values on the log scale
## Beta_* are effect of environment on species' parameters matrix[1 + n_env, n_species]
### 2. others: either a: random gamma, b: just replicate to loc dimensions.
transformParameters <- function(pars, Env, ranef = F, returndf = F) {
  
  Env <- as.data.frame(Env)
  
  ## 1. Envrionmentally-dependent transformations
  polyformula <- as.formula(paste("~", paste("poly(", colnames(Env), ", 2)", collapse = "+")))
  M <- model.matrix(polyformula, data = as.data.frame(Env))
  
  ### 1. env-dependent parameters  ###
  g_log <- M %*% pars$Beta_g
  m_j_log <- M %*% pars$Beta_m_j
  r_log <- M %*% pars$Beta_r
  s_log <- M %*% pars$Beta_s
  
  transpars <- list(g_loc = exp(g_log), m_j_loc = exp(m_j_log), r_loc = exp(r_log), s_loc = exp(s_log)) # env-dependent
  
  ## 2. Environmentally-independent transformation to a matrix with nrow == nrow(Env)
  if(ranef) {
    ## 2a. With random effects
    b_loc <- t(replicate(pars$n_locs, rgamma2(pars$b, pars$b, pars$shape_par)))
    c_j_loc <- t(replicate(pars$n_locs, rgamma2(pars$c_j, pars$c_j, pars$shape_par)))
    c_a_loc <- t(replicate(pars$n_locs, rgamma2(pars$c_a, pars$c_a, pars$shape_par)))
    c_b_loc <- t(replicate(pars$n_locs, rgamma2(pars$c_b, pars$c_b, pars$shape_par)))
    h_loc <- t(replicate(pars$n_locs, rgamma2(pars$h, pars$h, pars$shape_par)))
    m_a_loc <- t(replicate(pars$n_locs, rgamma2(pars$m_a, pars$m_a, pars$shape_par)))
    
  } else {
    ## 2b. Just a replicate without random effects
    b_loc <- t(replicate(pars$n_locs, pars$b))
    c_j_loc <- t(replicate(pars$n_locs, pars$c_j))
    c_a_loc <- t(replicate(pars$n_locs, pars$c_a))
    c_b_loc <- t(replicate(pars$n_locs, pars$c_b))
    h_loc <- t(replicate(pars$n_locs, pars$h))
    m_a_loc <- t(replicate(pars$n_locs, pars$m_a))
  }
  
  ## add env independent parameters to transpars
  transpars <- c(transpars,
                 list(b_loc = b_loc, c_j_loc = c_j_loc, c_a_loc = c_a_loc, c_b_loc = c_b_loc, h_loc = h_loc, m_a_loc = m_a_loc))
  
  ## untransformed env-dependent
  logpars <- list(g_log = g_log, m_j_log = m_j_log, r_log = r_log, s_log = s_log)
  
  if (returndf) {
    transpars <- lapply(transpars, as.data.frame)
    transpars <- do.call(cbind, transpars)
    transpars <- pivot_longer(cbind(transpars, Env),
                              cols = starts_with(c("g", "m_", "r", "s",
                                                   "b", "c_", "h")),
                              names_sep = "\\.",
                              names_to = c("parameter", "species"),
                              values_to = c("q")) %>%
      mutate(species = as.integer(as.factor(species)))
  } else {
    transpars <- c(pars, transpars, logpars)
  }
  
  return(transpars)
}


#### Simulate vector of initial states -------------------
generateInitialState <- function(n_species, n_stages = 3, generatelog = T) {
  s <-  rep(2:0, each = n_species) + rnorm(n_stages, 0, 0.1) # Initial state matrix.
  if(generatelog) return(s) else return(exp(s))
}


#### Simulate one time series for a specific fixed set of parameters. -------------------
simulateOneSeries <- function(state_init, times, pars, processerror = T, obserror = T, log = T) {
  
  
  n_pops <- length(state_init)
  
  if(!processerror) pars$sigma_process <- c(0, 0, 0)
  
  Sim <- calcModel(state_init, times = times, pars = pars)
  if(!log) Sim <- exp(Sim)
  
  Sim <- set_colnames(cbind(times, Sim), c("time", paste(1:ncol(Sim)))) # Make matrix look like from ode integrator for compatibility
  
  if (obserror) {
    Sim[, 2:(n_pops/3+1)] <- matrix(rnorm(Sim, Sim, pars$sigma_obs[1]), nrow = nrow(Sim))[,  2:(n_pops/3+1)]
    Sim[, (n_pops/3+1):(n_pops+1)] <- matrix(rnorm(Sim, Sim, pars$sigma_obs[2]), nrow = nrow(Sim))[, (n_pops/3+1):(n_pops+1)]
  }
  
  return(Sim)
}


#### Multiple simulations -------------------
simulateMultipleSeriesInEnv <- function(pars,
                                        Env,
                                        times = 2:32, # internal model will still start at 1
                                        modelstructure = c("ba", "ba-rect", "ba-rag", "ba-rag-ranef"),
                                        format = c("long", "stan"),
                                        obserror = T, processerror = NULL, ranef = NULL) {
  
  modelstructure <- match.arg(modelstructure)
  format <- match.arg(format)
  
  n_times <- length(times)
  
  ## set default values for modelstructures when not explicitly overridden
  if (is.null(ranef)) {
    ranef <- if (modelstructure == "ba-rag-ranef") T else F
  }
  
  if (is.null(processerror)) {
    processerror <- if (modelstructure %in% "ba") T else F
  }
  
  plotlevelinit <- F ## whether initial values are different within one cluster (= loc)
  
  ## ---- Set parameters ----
  pars <- transformParameters(pars, Env, ranef = ranef)
  meta <- pars # for info like n_stages, dbh_lower_a
  
  #### Internal function to mapply over locations with different environments and therefore possibly different parameters
  simulateOneSeriesInEnv <- function(g_, m_j_, r_, s_,
                                     b_, c_j_, c_a_, c_b_, h_, m_a_) {
    
    ## replace global pars with local variants
    simpars <- within(pars, {g <- g_; m_j <- m_j_; r <- r_; s <- s_;
    b <- b_; c_j <- c_j_; c_a <- c_a_; c_b <- c_b_; h <- h_; m_a <- m_a_;}
    ) # list of dependent and independent parameters for one environment
    
    
    ## Switch whether init is the sane for all plots within cluster
    if (plotlevelinit) {
      seriesatloc <- replicate(pars$n_plotsperloc,
                               simulateOneSeries(generateInitialState(n_species = pars$n_species), times, simpars,
                                                 processerror = processerror, obserror = obserror), simplify = F)
    } else {
      # Post initial state replication
      m0 <- generateInitialState(n_species = pars$n_species)
      seriesatloc <- replicate(pars$n_plotsperloc,
                               simulateOneSeries(m0, times, simpars,
                                                 processerror = processerror, obserror = obserror), simplify = F)
    }
    
    
    names(seriesatloc) <- 1:pars$n_plotsperloc
    return(seriesatloc)
  }
  
  ## Sims is a: list[n_loc] <- list[n_plotsperloc] <- matrix[n_times, n_species]
  sims <- mapply(simulateOneSeriesInEnv,
                 ## env-dependent
                 as.data.frame(t(pars$g_loc)),
                 as.data.frame(t(pars$m_j_loc)),
                 as.data.frame(t(pars$r_loc)),
                 as.data.frame(t(pars$s_loc)),
                 
                 ## env-independent
                 as.data.frame(t(pars$b_loc)),
                 as.data.frame(t(pars$c_j_loc)),
                 as.data.frame(t(pars$c_a_loc)),
                 as.data.frame(t(pars$c_b_loc)),
                 as.data.frame(t(pars$h_loc)),
                 as.data.frame(t(pars$m_a_loc)),
                 
                 SIMPLIFY = F)
  
  
  
  ###### Fomat the data from multiple simulations -------------------
  ## This function is designed to work only within environment simulateMultipleSeriesInEnv
  formatSims <- function() {
    
    isrectangular <- grepl("rect$", modelstructure)
    
    #### DECISION structure for formatting
    ## 0. Do general long casting
    ##    if format == "long": return(Sims) data.frame!
    ##
    ## 1. else if: format == "stan":
    ##      1a. if: is rectangular: reuse original sims list for rectangular casting, then construct stansims list
    ##      1b. else: construct stansims list
    ##   return(stansims)
    
    ######### (0) ##########
    sims_loc <- lapply(sims, function(l) do.call(rbind, l)) %>% # bind the series within a location into one table
      lapply(function(m) cbind(m, plot = rep(1:pars$n_plotsperloc, each = n_times))) # attach an id to series per location. Each n_times, because simulateOneSeries() does only return n_times (n_times_intern/2) observations
    
    ## Long casting is needed for both format == "long" and "stan"
    Sims <- do.call(rbind, sims_loc) # bind all lists of locations into one table
    Sims <- cbind(Sims, loc = rep(1:pars$n_locs, each = pars$n_plotsperloc * n_times)) # attach an id to series per location
    Sims <- cbind(Sims, Env[Sims[,"loc"],])
    
    ## Here come all different long formats used in different hierarchical levels of the stan fit
    Sims <- tidyr::pivot_longer(as.data.frame(Sims),
                                cols = all_of(paste(1:(pars$n_species*meta$n_stages))),
                                names_to = "pop",
                                values_to = "abundance") %>%
      mutate(pop = as.integer(pop)) %>%
      mutate(species = rep(1:pars$n_species, times = meta$n_stages, length.out = nrow(.))) %>%
      mutate(stage = rep(c("j", "a", "b"), each = pars$n_species, length.out = nrow(.))) %>%
      mutate(obsmethod = rep(c(1, 2, 2), each = pars$n_species, length.out = nrow(.)))
    
    ## Sims <- complete(Sims, group = nesting("loc", "time", species", "stage"), fill = 0) # is assumed!
    Sims <- group_by(Sims, loc) %>%
      mutate(isy0 = time == min(time)) %>%
      ungroup() %>%
      mutate(stage = factor(stage, levels = c("j", "a", "b"), ordered = T)) %>%
      arrange(loc, time, pop, plot)
    
    
    ## Get log and exped version of the response
    Sims <- mutate(Sims, abundance_log = abundance,
                   abundance = exp(abundance)
    )
    
    
    if (format == "long") {
      
      return(Sims)
      
    }
    
    ######### (1) ##########
    else if (format == "stan"){
      
      #### Special long data set structures with different amounts of long casting along different ids
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
      
      
      ## Prepare design matrix
      polyformula <- as.formula(paste("~", paste("poly(", colnames(as.data.frame(Env)), ", 2)", collapse = "+")))
      X <- model.matrix(polyformula, data = as.data.frame(Env))
      
      ######### (1a) ##########
      if (isrectangular) {
        
        ## Go back to the original sim list for rectangular casting
        time_rect <- sims[[1]][[1]][,1] # get the time column  from the very first table
        sims_rect <- lapply(sims, function(loc) lapply(loc, function (plot) as.data.frame(t(plot)))) # make everything into lists!
        sims_rect <- array(unlist(sims_rect),
                           dim = c(length(sims_rect[[1]][[1]][[1]]), length(sims_rect[[1]][[1]]), length(sims_rect[[1]]), length(sims_rect))
        ) 
        sims_rect <- aperm(sims_rect)
        # str(S) ## # 'loc', 'plot', 'pop', 'time'
        sims_rect <- sims_rect[ , , , -1] # drop the first 'column', which is not a population but the times
        
        Times <- t(replicate(pars$n_locs, times)) - min(times) + 1
        
        stansims <- list(N_locs = pars$n_locs,
                         N_plots = pars$n_plotsperloc,
                         N_times = n_times,
                         N_species = pars$n_species,
                         N_pops = pars$n_stages * pars$n_species,
                         N_beta = pars$n_beta,
                         
                         rep_obsmethod2pops = rep(c(1, 2, 2), each = pars$n_species),
                         i_j = 1:pars$n_species,
                         i_a = (1:pars$n_species) + pars$n_species,
                         i_b = (1:pars$n_species) + 2*pars$n_species,
                         
                         dbh_lower_a = meta$dbh_lower_a,
                         dbh_lower_b = meta$dbh_lower_b,
                         
                         times = t(replicate(pars$n_locs, times)) - min(times) + 1,
                         time_max = Times[,ncol(Times)],
                         timespan_max = max(Times[,ncol(Times)]) - 1, # difference to the first time
                         
                         X = X,
                         
                         y_log = sims_rect,
                         y = exp(sims_rect)
        )
        
      }
      
      ######### (1b) ##########
      else if (modelstructure %in% c("ba-rag-ranef", "ba-rag", "ba")) {
        
        Sims_species <- Sims_init[c("loc", "species")] %>% unique()
        N_init <- nrow(Sims_init)
        N_times <- nrow(Sims_times)
        N_yhat <- nrow(Sims_yhat)
        N_locs <-  nrow(Sims_locs)
        time_init <-  Sims$time[match(Sims_locs$loc, Sims$loc)] # N_locs!
        
        vrep <- function(...) unlist( c( Vectorize(rep.int)(...) ) ) # :)))))
        
        stansims <- list(
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
          
          dbh_lower_a = meta$dbh_lower_a,
          dbh_lower_b = meta$dbh_lower_b,
          
          species = Sims_species$species,
          stage = Sims_init$stage,
          pops = Sims_init$pop,
          time_init = time_init,
          time_max_data = Sims_locs$time_max,
          times_data = Sims_times$time,
          
          X = X,
          
          y0 = Sims_y0$abundance,
          y = Sims_reobs$abundance,
          
          y0_log = Sims_y0$abundance_log,
          y_log = Sims_reobs$abundance_log
        )
        
      }
    }
    
    return(stansims)
  }
  
  Sims <- formatSims()
  
  attr(Sims, "pars") <- pars
  return(Sims)
}



#### Returns a list from within-function environment -------------------
generateParameters <- function(seed = 1,
                               
                               ## "meta data":
                               n_species = 4, n_stages = 3, n_locs = 70, n_plotsperloc = 4, n_env = 2,
                               dbh_lower_a = 100, dbh_lower_b = 200,
                               
                               ## errors:
                               sigma_process = c(0.05, 0.05, 0.03), shape_par = 10, sigma_obs = c(1, 0.2),
                               ...) {
  set.seed(seed)
  
  ## Beta_* are effect of environment on species' parameters matrix[n_env, n_species]
  ## These are on the log scale!
  n_beta <- 1 + 2 * n_env # degree 2 polynomial + intercept
  
  ### 1. env-dependent parameters (log scale) ### 
  Beta_g <- matrix(rnorm(n_beta*n_species, -0.5, 0.5), n_beta, n_species)
  g_log <- rnorm(n_species, -12, 0.2)
  Beta_g[1,] <- g_log
  Beta_g[3,] <- rnorm(n_species, -1, 0.2)
  
  Beta_m_j <- matrix(rnorm(n_beta*n_species, -3, 0.5), n_beta, n_species)
  m_j_log <- rnorm(n_species, -1.2, 0.1)
  Beta_m_j[1,] <- m_j_log
  Beta_m_j[3,] <- rnorm(n_species, -2.5, 0.2)
  
  Beta_r <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  r_log <- rnorm(n_species, 2.2, 0.2)
  Beta_r[1,] <- r_log
  Beta_r[3,] <- rnorm(n_species, -0.5, 0.2)
  
  Beta_s <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  s_log <- rnorm(n_species, -2.9, 0.1)
  Beta_s[1,] <- s_log
  Beta_s[3,] <- rnorm(n_species, -1, 0.2)
  
  
  ### 2. env-independent parameters ("state model scale", positive) ### 
  b <- exp(rnorm(n_species, -1, 0.01))
  c_j <- exp(rnorm(n_species, -10, 0.5))
  c_a <- exp(rnorm(n_species, -3.3, 0.2))
  c_b <- exp(rnorm(n_species, -3, 0.1))
  h <- exp(rnorm(n_species, -0.5, 0.3))
  m_a <- exp(rnorm(n_species, -1.5, 0.5))
  
  
  pars <- c(as.list(environment()), list(...)) ## Carries along all parameter values within the named list.
  return(pars)
}




######################################################################################
# Demo Simulations ------------------------------------------------------------------
######################################################################################

## One time series ---------------------------------------------------------------

## Direct plotting
# matplot(exp(calcModel(times1, exp(initialstate1), pars1)[times1[-c(1:3)],]), type = "b")

#### prerequisites for a single time series
formals(generateParameters)
pars1 <- generateParameters(seed = 1)
initialstate1 <- generateInitialState(n_species = pars1$n_species)
## transform parameters manually, use *_log because they are vectors[n_species] (instead of matrices as produced by transformParameters)
pars1 <- within(pars1,
                {g = exp(pars1$g_log); m_j = exp(m_j_log); r = exp(pars1$r_log); s = exp(pars1$s_log);}
)
times1 <- 1:30

Sim1 <- simulateOneSeries(initialstate1, times = times1, pars = pars1,
                          processerror = F, obserror = F, log = F)

matplot(Sim1[,-1], type = "b", ylab="N",
        pch = rep(c("J", "A", "B"), each = pars1$n_species),
        col = 1:pars1$n_species, xlab = "time") # log='y'



## Multiple time series in env ---------------------------------------------------------------
pars_demo <- generateParameters(n_species = 4, n_locs = 50, sigma_obs = c(0.2, 0.1))
E_demo <- simulateEnv(pars_demo$n_env, pars_demo$n_locs)
P_env_demo <- transformParameters(pars_demo, E_demo, ranef = F, returndf = T)
S_demo <- simulateMultipleSeriesInEnv(pars_demo, E_demo, ranef = F, processerror = T, obserror = F)

#### Plot time series
S_demo %>%
  # filter(stage %in% c("a")) %>%
  ggplot(mapping = aes(x = time, y = abundance, color = env_1, lty = stage, group = interaction(loc, stage, plot))) +
  geom_line() +
  facet_wrap(facets = c("species"))


#### Plot parameters
P_env_demo %>%
  filter(parameter %in% c("g_loc")) %>%
  gf_point(q ~ env_1 | species) %>%
  gf_smooth()

ggplot(filter(P_env_demo, parameter == "r_loc"),
       mapping = aes(x = env_1, y = q, color = species, group = species)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(facets = c("species"))

P_env_demo %>%
  filter(parameter %in% c("s_loc")) %>%
  gf_point(q ~ env_1 | species) %>%
  gf_smooth()



######################################################################################
# Fit simulated data with stan    ----------------------------------------------------
######################################################################################

## Orient and compile model. --------------------------------------------------------------

modeldir <- dir(pattern = glue("^(Model).*ba$"))

modelname <- c("ba",
               "ba-rect", # fully rectangular, without process error
               "ba-rag", # ragged, clusterwise-initialization
               "ba-rag-ranef" # like ba-rag but with random demographic pars
)[2]

modelpath <- file.path(modeldir, glue('Model_{modelname}.stan'))

model <- cmdstan_model(modelpath)
# model <- cmdstan_model(modelpath, force_recompile = T)
# compiledpath <-  model$exe_file()



## Simulate stan model data --------------------------------------------------------------
fitseed <- 4

pars <- generateParameters(seed = fitseed, n_locs = 70)

Env <- simulateEnv(n_env = pars$n_env, n_locs = pars$n_locs)
data <- simulateMultipleSeriesInEnv(pars, Env, times = c(3, 15, 25), modelstructure = modelname, format = "stan")
Data_long <- simulateMultipleSeriesInEnv(pars, Env, times = c(3, 15, 25),
                                         modelstructure = modelname, format = "long") %>%
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


#### Returns +- the true start values ----------------------
getTrueInits <- function() {
  
  isragged <- grepl("^ba-rag", modelname)
  
  truepars <<- attr(data, "pars")
  
  inits <- list(
    ## from global environment
    Beta_g = truepars$Beta_g,
    Beta_m_j = truepars$Beta_m_j,
    Beta_r = truepars$Beta_r,
    Beta_s = truepars$Beta_s,
    
    state_init_log = if (isragged) rnorm(data$y0_log, data$y0_log, 0.01) else
                                   data$y_log[,,1,] + rnorm(data$y_log[,,1,], 0, 0.01),
    
    ## Version with random rates.
    b     = truepars$b,
    c_j   = truepars$c_j,
    c_a   = truepars$c_a,
    c_b   = truepars$c_b,
    h     = truepars$h,
    m_a   = truepars$m_a,
    
    shape_par      = truepars$shape_par,
    sigma_process  = truepars$sigma_process,
    sigma_obs      = truepars$sigma_obs,
    
    u = replicate(pars$n_locs, matrix(rnorm(pars$n_species*3, 0, 0), nrow = pars$n_species, ncol = data$timespan_max)),
    
    b_loc     = truepars$b_loc,
    c_j_loc   = truepars$c_j_loc,
    c_a_loc   = truepars$c_a_loc,
    c_b_loc   = truepars$c_b_loc,
    h_loc     = truepars$h_loc,
    m_a_loc   = truepars$m_a_loc
  )
  
  return(inits)
}

#### Returns viable start values ---------------------------
getInits <- function() {
  
  isragged <- grepl("^ba-rag", modelname)
  
  truepars <<- attr(data, "pars")
  
  inits <- list(
    ## from global environment
    Beta_g = matrix(c(1, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species),
    Beta_m_j = matrix(c(-2, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species),
    Beta_r = matrix(c(4, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species),
    Beta_s = matrix(c(-3, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species),
    
    state_init_log = if (isragged) rnorm(data$y0_log, data$y0_log, 0.01) else
      data$y_log[,,1,] + rnorm(data$y_log[,,1,], 0, 0.01),

    
    ## Version with random rates.
    b     = rgamma(truepars$n_species, 10, 10/0.2),
    c_j   = rgamma(truepars$n_species, 10, 10/0.001),
    c_a   = rgamma(truepars$n_species, 10, 10/0.001),
    c_b   = rgamma(truepars$n_species, 10, 10/0.001),
    h     = rgamma(truepars$n_species, 10, 10/0.1), #!
    m_a   = rgamma(truepars$n_species, 10, 10/0.001),
    
    shape_par      = c(10, 10, 10),
    sigma_process  = c(0.1, 0.1, 0.1),
    sigma_obs      = c(10, 10),
    
    u = replicate(pars$n_locs, matrix(rnorm(pars$n_species*3, 0, 0), nrow = pars$n_species, ncol = data$timespan_max)),
    
    b_loc     = truepars$b_loc,
    c_j_loc   = truepars$c_j_loc,
    c_a_loc   = truepars$c_a_loc,
    c_b_loc   = truepars$c_b_loc,
    h_loc     = truepars$h_loc,
    m_a_loc   = truepars$m_a_loc
  )
  
  return(inits)
}



## Draw from model --------------------------------------------------------------

####  Do the fit ----------------------
drawSamples <- function(model, data, variational = F, n_chains = 1, initfunc = getInits) {
  if(variational) {
    fit <- model$variational(data = data,
                             output_dir = "Fits.nosync",
                             init = initfunc,
                             iter = 20**4) # convergence after 1500 iterations
    
  } else {
    
    fit <- model$sample(data = data,
                        output_dir = "Fits.nosync",
                        init = initfunc,
                        iter_warmup = 1500, iter_sampling = 500,
                        chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
  }
  return(fit)
}


fit <- drawSamples(model, data, variational = F, initfunc = getTrueInits)
fit$output()
fit$time()
fit$init()
fit$draws(variables = NULL, inc_warmup = F)
fit$return_codes()



## Extract draws -----------------------------------------------------------
drawpath <- fit$output_files() # dput(fit$output_files())

## alternatively
# drawpath <- file.info(list.files("Fits", full.names = T)) %>%
#   arrange(desc(mtime)) %>%
#   slice(1:6) %>%
#   rownames()


### Get draws
## rstan format for use in other packages
draws <- rstan::read_stan_csv(drawpath)
draws

## get draws with cmdstanr
# draws <- read_cmdstan_csv(drawpath, variables = NULL)



## Summary -----------------------------------------------------------------

fit$summary() %>%
  View()
# fit$summary(variables = c("shape_par", "sigma_obs", "h", "c_j", "state_init", "u")) %>% # "o", "state_init"

shinystan::launch_shinystan(draws)



# Inspection --------------------------------------------------------------
truepars <- attr(data, "pars")

#### Draws vs true ----------------------
plotDrawVsSim <- function(parname = "h",
                          simdata = data,
                          rstandraws = draws) {
  
  
  simpar <- attr(simdata, "pars")[[parname]]
  
  if(is.vector(simpar)) {
    stanpar <- rstan::extract(rstandraws, pars = parname)[[1]] %>% t() %>% c()
    
  } else {
    stanpar <- rstan::extract(rstandraws, pars = parname)[[1]] %>% aperm(c(2,3,1)) %>% c()
  }
  
  
  # cbind(simpar, stanpar
  plot(rep(c(simpar), length.out = length(stanpar)), stanpar,
       cex = 0.1, col = alpha("black", alpha = 0.03), main = parname, ylab = "Draw", xlab = "Sim")
  abline(a = 0, b = 1, col = "green")
}


## Fitted parameters vs. true
plotDrawVsSim("b")
plotDrawVsSim("c_j")
plotDrawVsSim("c_a")
plotDrawVsSim("c_b")
plotDrawVsSim("h")
plotDrawVsSim("m_a")

plotDrawVsSim("Beta_g")
plotDrawVsSim("Beta_m_j")
plotDrawVsSim("Beta_r")
plotDrawVsSim("Beta_s")



#### Parameter vs Env ----------------------
predictParameterInEnv <- function(B, x = seq(-1, 1, by = 0.01), n_beta = 2, species = 1) {
  E <- as.data.frame(replicate(n_beta, x))
  polyformula <- as.formula(paste("~", paste("poly(", colnames(E), ", 2)", collapse = "+")))
  X <- model.matrix(polyformula, data = E)
  X %*% B[,1]
}

# poly2 <- function(x, b) {b[1] + x*b[2] + x^2*b[3]}
# getVertex <- function(b) {-b[2]/(2*b[3])}

plotParameterInEnv <- function(betaname, x = Env[,1]) {
  truepars <- attr(data, "pars")
  Beta_true <- truepars[[betaname]]
  Beta_mean <- matrix(fit$summary(variables = betaname)$mean, data$N_beta, data$N_totalspecies)
  
  X <- cbind(true = predictParameterInEnv(Beta_true, x), meandraw = predictParameterInEnv(Beta_mean, x))
  matplot(x, X, col = c("blue", "black"), pch = c("T", "D"), type = "p")
}


plotParameterInEnv("Beta_g")
plotParameterInEnv("Beta_m_j")
plotParameterInEnv("Beta_r")
plotParameterInEnv("Beta_s")

