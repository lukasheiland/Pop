
# Stats ------------------------------------------------------------

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



# Simulation --------------------------------------------------------------


#### Get environmental variables -------------------
simulateEnv <- function(n_env, n_locs){
  n <- n_env * n_locs
  Env <- matrix(runif(n, -1, 1), nrow = n_locs, ncol = n_env)
  colnames(Env) <- paste("env", 1:n_env, sep = "_")
  return(Env)
}



#### Simulate vector of initial states -------------------
generateInitialState <- function(n_species, n_stages = 3, generatelog = T) {
  s <-  rep(2:0, each = n_species) + rnorm(n_stages, 0, 0.1) # Initial state matrix.
  if(generatelog) return(s) else return(exp(s))
}



#### Transform parameters -------------------
## Accepts and returns a list of parameters. No parameter is lost.
##
### 1. Generate environment parameters for the population model from expected values on the log scale
## Beta_* are effect of environment on species' parameters matrix[1 + n_env, n_species]
### 2. others: either a: random gamma, b: just replicate to loc dimensions.
transformParameters <- function(pars, Env,
                                envdependent = c(b = F, c_a = F, c_b = F, c_j = F, g = F, h = F, m_a  = F, m_j = F, r = F, s = F),
                                ranef = F, returndf = F) {
  
  
  parname_envdep <- names(envdependent)[envdependent]
  parname_envindep <- names(envdependent)[!envdependent]
  
  Env <- as.data.frame(Env)
  polyformula <- as.formula(paste("~", paste("poly(", colnames(Env), ", 2)", collapse = "+")))
  M <- model.matrix(polyformula, data = as.data.frame(Env))
  
  ### 1. env-dependent parameters  ###
  multiplyEnv <- function(parname) {
    Par <- M %*% pars[[glue("Beta_{parname}")]]
    return(Par)
  }
  
  pars_envdep <- sapply(parname_envdep, multiplyEnv, USE.NAMES = T, simplify = F) %>%
    setNames(glue("{str_to_title(parname_envdep)}_log"))  # ! name env-depedent matrices to upper
  
  pars_envdep <- lapply(pars_envdep, exp) %>%
    setNames(glue("{str_to_title(parname_envdep)}_loc")) %>%
    c(pars_envdep)
  
  ## 2. Environmentally-independent transformation to a matrix with nrow == nrow(Env)
  replicateWithRanef <- function(parname) {
    Par <- t(replicate(pars$n_locs, rgamma2(pars[[parname]], pars[[parname]], pars$shape_par)))
    return(Par)
  }
  
  replicateToLoc <- function(parname) {
    Par <- t(replicate(pars$n_locs, pars[[parname]]))
    return(Par)
  }
  
  if(ranef) {
    ## 2a. With random effects
    pars_envindep <- sapply(parname_envindep, replicateWithRanef, USE.NAMES = T, simplify = F) %>%
      setNames(glue("{str_to_title(parname_envindep)}_loc"))  # ! name env-indepedent matrices to upper
    
  } else {
    ## 2b. Just a replicate without random effects
    pars_envindep <- sapply(parname_envindep, replicateToLoc, USE.NAMES = T, simplify = F) %>%
      setNames(glue("{str_to_title(parname_envindep)}_loc"))
  }
  
  
  if (returndf) {
    
    transpars <- lapply(pars_envdep, as.data.frame)
    transpars <- do.call(cbind, transpars)
    transpars <- cbind(transpars, Env)
    
    transpars <- pivot_longer(transpars,
                              cols = starts_with(names(envdependent), ignore.case = T),
                              names_sep = "\\.",
                              names_to = c("parameter", "species"),
                              values_to = c("q")) %>%
      mutate(parameter = str_to_lower(parameter)) %>%
      mutate(species = as.integer(as.factor(species)))
    
  } else {
    
    transpars <- c(pars, pars_envdep, pars_envindep, list(envdependent = envdependent))
    
  }
  
  return(transpars)
}



#### Fomat the data from multiple simulations -------------------
## This function is designed to work only within environment simulateMultipleSeriesInEnv()
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
  
  if (format == "long") { } # Sims <- Sims
  
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
      sims_rect <- aperm(sims_rect, perm = c(4, 2, 3, 1))
      # str(S) ## # 'loc', 'time', 'plot', 'pop'
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
                       
                       y0_loc_log = apply(exp(sims_rect), c(1, 2, 4), mean)[, 1, ] %>% log(),
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
    
    ## if (format == "stan")
    Sims <- stansims
    
  } ## closing if (format == "stan")
  
  return(Sims)
}



# Fitting -----------------------------------------------------------------

#### Returns +- the true start values ----------------------
getTrueInits <- function() {
  
  responsescaleerror <- 0
  
  isragged <- grepl("^ba-rag", modelname)
  
  truepars <- attr(data, "pars")
  truepars <- truepars[sapply(truepars, is.numeric)]
  parnames <- names(truepars)
  names(truepars) <- ifelse(str_ends(parnames, "_loc"), str_to_lower(parnames), parnames)
  
  newpars <- list(
    state_init_log = if (isragged) rnorm(data$y0_log, data$y0_log, 0.01) else
      data$y_log[,1,1,] + rnorm(data$y_log[,1,1,], 0, responsescaleerror),
    
    u = replicate(truepars$n_locs, matrix(rnorm(truepars$n_species*3, 0, 0.001), nrow = pars$n_species*3, ncol = data$timespan_max))
  )
  
  inits <- c(newpars, truepars)
  truepars <<- truepars
  
  return(inits)
}

#### Returns viable start values ---------------------------
getInits <- function() {
  
  responsescaleerror <- 1
  
  isragged <- grepl("^ba-rag", modelname)
  
  truepars <<- attr(data, "pars")
  n_species <- truepars$n_species
  n_locs <- truepars$n_locs
  
  
  b_log <- rnorm(n_species, -1, 0.01)
  c_a_log <- rnorm(n_species, -3.3, 0.2)
  c_b_log <- rnorm(n_species, -3, 0.1)
  c_j_log <- rnorm(n_species, -10, 0.5)
  g_log <- rnorm(n_species, -12, 0.2)
  h_log <- rnorm(n_species, -0.5, 0.3)
  m_a_log <- rnorm(n_species, -1.5, 0.5)
  m_j_log <- rnorm(n_species, -1.2, 0.1)
  r_log <- rnorm(n_species, 2.2, 0.2)
  s_log <- rnorm(n_species, -2.9, 0.1)
  
  beta_null <- function() { matrix(c(-1, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.001) }
  
  
  inits <- list(
    
    b_log = b_log,
    c_a_log = c_a_log,
    c_b_log = c_b_log,
    c_j_log = c_j_log,
    g_log = g_log,
    h_log = h_log,
    m_a_log = m_a_log,
    m_j_log = m_j_log,
    r_log = r_log,
    s_log = s_log,
    
    ## from global environment
    Beta_b = matrix(c(1, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species),
    Beta_c_a = beta_null(),
    Beta_c_b = beta_null(),
    Beta_c_j = beta_null(),
    Beta_g = matrix(c(1, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_h = matrix(c(1, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_m_a = matrix(c(-4, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_m_j = matrix(c(-2, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_r = matrix(c(3, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_s = matrix(c(-3, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    
    b = exp(b_log),
    c_a = exp(c_a_log),
    c_b = exp(c_b_log),
    c_j = exp(c_j_log),
    g = exp(g_log),
    h = exp(h_log),
    m_a = exp(m_a_log),
    m_j = exp(m_j_log),
    r = exp(r_log),
    s = exp(s_log),
    
    b_loc = matrix(rep(exp(b_log), n_locs), nrow = n_locs, byrow = T),
    c_a_loc =  matrix(rep(exp(c_a_log),n_locs), nrow = n_locs, byrow = T),
    c_b_loc =  matrix(rep(exp(c_b_log),n_locs), nrow = n_locs, byrow = T),
    c_j_loc =  matrix(rep(exp(c_j_log),n_locs), nrow = n_locs, byrow = T),
    g =  matrix(rep(exp(g_log),n_locs), nrow = n_locs, byrow = T),
    h =  matrix(rep(exp(h_log),n_locs), nrow = n_locs, byrow = T),
    m_a_loc =  matrix(rep(exp(m_a_log),n_locs), nrow = n_locs, byrow = T),
    m_j_loc =  matrix(rep(exp(m_j_log),n_locs), nrow = n_locs, byrow = T),
    r_loc =  matrix(rep(exp(r_log),n_locs), nrow = n_locs, byrow = T),
    s_loc =  matrix(rep(exp(s_log),n_locs), nrow = n_locs, byrow = T),
    
    shape_par      = c(10, 10, 10),
    sigma_process  = c(0.01),
    sigma_obs      = c(1.1, 1.1),
    
    state_init_log = if (isragged) rnorm(data$y0_log, data$y0_log, 0.01) else
      data$y_log[,1,1,] + rnorm(data$y_log[,1,1,], 0, responsescaleerror),
    
    u = replicate(truepars$n_locs, matrix(rnorm(truepars$n_species*3, 0, 0.001), nrow = pars$n_species*3, ncol = data$timespan_max))
    )
  
  return(inits)
}
