# Library -----------------------------------------------------------------
library(tidyverse)
library(ggformula)
library(magrittr)
library(glue)
library(tictoc)

library(rstan)
rstan_options(javascript = FALSE)
library(cmdstanr)
# install_cmdstan(cores = 3)

library(bayesplot)


# ——————————————————————————————————————————————————————————————————————————————————— #
# Model formulation  ----------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #


#### Returns states for times -------------------
calcModel <- function(times,
                      initialstate, # A vector of species states.
                      pars # internal number of discrete time steps within a unit of time
) {
  
  ## Just unpacking for readability.
  
  b <- pars$b # Length n_species vector of basal area increment rates.
  c_a <- pars$c_a # Length n_species vector ## only used in model versions, where there are distinct parameters for competition on b an a
  c_b <- pars$c_b # Length n_species vector
  c_j <- pars$c_j # Length n_species vector
  g <- pars$g # Length n_species vector of transition rates.
  h <- pars$h # Length n_species vector of transition rates.
  l <- pars$l # Length n_species vector of input rates.
  m_a <- pars$m_a # Length n_species vector of mortalities
  m_j <- pars$m_j # Length n_species vector of mortalities
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  
  ba_a_upper <- pars$ba_a_upper
  ba_a_avg <- pars$ba_a_avg
  
  sigma <- pars$sigma_process
  
  ## Set the count variables
  n <- length(r) # no. of species
  
  ## 
  times_intern <- 1:max(times)
  n_times <- length(times_intern)
  
  # Prepare a state matrix
  whichstate <- rep(1:3, each = n)
  State <- matrix(rep(initialstate, times = n_times), nrow = n_times, byrow = T)
  
  ## Here comes the model.
  for (t in 2:n_times) {
    
    ## States at t-1: *State*, and *State*
    J <- State[t-1,whichstate == 1]
    A <- State[t-1,whichstate == 2]
    B <- State[t-1,whichstate == 3]
    
    ## The total basal area of big trees
    BA <- A * ba_a_avg + B
    BA_sum <- sum(BA)
    
    # observation error for stages/times/species
    u <- matrix(0, nrow = length(s), ncol = 3) + rnorm(length(s)*3, 0, sigma)
    
    ## Two kinds of processes acting ot State. multiplicative in exp(...) and additive at the end
    
    J_trans <- r*BA + l + (J - g*J) # use rlnorm + u[ ,1]
    J_t <- J_trans * 1/(1 + c_j*sum(J) + s*BA_sum + m_j) # count of juveniles J
    # with all variables > 0, and 0 < g, m_j < 1
    
    A_trans <- g*J + (A - h*A) # + u[ ,2]
    A_t <-  A_trans * 1/(1 + c_a*BA_sum + m_a) # count of small adults
    # with all variables > 0, and 0 < h, m_a < 1
    
    A_ba <- A * h * ba_a_upper # Basal area of small adults A. Conversion by multiplication with basal area of State exit (based on upper dhh boundary of the class)
    B_trans <- A_ba + B # + u[ ,3]
    B_t <- (1+b)*B_trans * 1/(1 + c_b*BA_sum)  # basal area of big adults B
    ## b is the net basal area increment (including density-independent m) basically equivalent to a Ricker model, i.e. constant increment rate leading to exponential growth, negative density dependent limitation scaled with total BA_sum.
    
    State[t, ] <- c(J_t, A_t, B_t)
  }
  
  whichtimes <- which(times_intern %in% times)
  
  return(State[whichtimes,])
  
}


# ——————————————————————————————————————————————————————————————————————————————————— #
# Simulation functions ---------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #


#### Simulate one time series for a specific fixed set of parameters. -------------------
simulateOneSeries <- function(state_init, times, pars,
                              processerror = T, obserror = T, exping = F, nominaltimes = NULL) {
  
  
  n_pops <- length(state_init)
  
  if(!processerror) pars$sigma_process <- 0
  
  Sim <- calcModel(state_init, times = times, pars = pars)
  
  if(exping) Sim <- exp(Sim)
  
  if(is.null(nominaltimes)) nominaltimes <- times ## this is for setting the time column in the returned matrix
  
  Sim <- set_colnames(cbind(nominaltimes, Sim), c("time", paste(1:ncol(Sim)))) # Make matrix look like from ode integrator for compatibility
  
  if (obserror) {
    
    robs <- function(n, mu, errorgroup) {
      
      switch(pars$obsprocess,
             normal = rnorm(n, mu, pars$sigma_obs[errorgroup]),
             negbinomial = rnbinom2(n, mu, pars$phi_obs[errorgroup]),
             gamma = rgamma2(n, mu, pars$alpha_obs[errorgroup])
      )
    }
    Sim[, 2:(n_pops/3+1)] <- matrix(robs(Sim, Sim, errorgroup = 1), nrow = nrow(Sim))[,  2:(n_pops/3+1)]
    Sim[, (n_pops/3+1):(n_pops+1)] <- matrix(robs(Sim, Sim, errorgroup = 2), nrow = nrow(Sim))[, (n_pops/3+1):(n_pops+1)]
    
  }
  
  return(Sim)
}


#### Multiple simulations -------------------
simulateMultipleSeriesInEnv <- function(pars,
                                        Env,
                                        times = 2:22, # internal model will still start at 1
                                        envdependent = c(b = F, c_a = F, c_b = F, c_j = F, g = F, h = F, l = F, m_a  = F, m_j = F, r = F, s = F),
                                        logstate = F,
                                        modelstructure = c("ba", "ba-rect", "ba-rag", "ba-rag-ranef"),
                                        format = c("long", "stan"), priorfactor = 10,
                                        obserror = T, processerror = NULL, ranef = NULL, independentstart = F) {
  
  modelstructure <- match.arg(modelstructure)
  format <- match.arg(format)
  
  n_times <- length(times)
  
  ## Times are just used internally to let simulations start independently, the data still includes the times vector
  Times <- rep(times, pars$n_locs) %>% matrix(nrow = n_times) %>% as.data.frame()
  nominaltimes <- NULL
  
  if(independentstart) {
    nominaltimes <- times
    Times <- lapply(Times, function(t) t + sample(0:max(times*5), 1, replace = T)) # prob = (max(times*3):0)^6)
  }
  
  ## set default values for modelstructures when not explicitly overridden
  if (is.null(ranef)) {
    ranef <- if (modelstructure == "ba-rag-ranef") T else F
  }
  
  if (is.null(processerror)) {
    processerror <- if (modelstructure %in% "ba") T else F
  }
  
  plotlevelinit <- F ## whether initial values are different within one cluster (= loc)
  
  ## ---- Set parameters ----
  pars <- transformParameters(pars, Env, envdependent = envdependent, ranef = ranef)
  meta <- pars # for info like n_stages, dbh_lower_a
  
  #### Internal function to mapply over locations with different environments and therefore possibly different parameters
  simulateOneSeriesInEnv <- function(b_, c_a_, c_b_, c_j_, g_, h_, l_, m_a_, m_j_, r_, s_, times_) {
    
    
    ## list of env-dependent and -independent parameters for one environment, replacing the original global parameters in "pars"
    simpars <- within(pars, {
      b <- b_;
      c_a <- c_a_; c_b <- c_b_; c_j <- c_j_;
      g <- g_; h <- h_;
      l <- l_;
      m_a <- m_a_; m_j <- m_j_;
      r <- r_;
      s <- s_;})
    
    
    ## Switch whether init is the sane for all plots within cluster
    if (plotlevelinit) {
      
      seriesatloc <- replicate(pars$n_plotsperloc,
                               simulateOneSeries(generateInitialState(n_species = pars$n_species, logstate = logstate),
                                                 times_, simpars,
                                                 processerror = processerror, obserror = obserror, nominaltimes = nominaltimes),
                               simplify = F)
    } else {
      
      ## cluster level init
      # Replication is post initial state generation
      m0 <- generateInitialState(n_species = pars$n_species, logstate = logstate)
      seriesatloc <- replicate(pars$n_plotsperloc,
                               simulateOneSeries(m0, times_, simpars,
                                                 processerror = processerror, obserror = obserror, nominaltimes = nominaltimes),
                               simplify = F)
    }
    
    
    names(seriesatloc) <- 1:pars$n_plotsperloc
    return(seriesatloc)
  }
  
  ## Sims is a: list[n_loc] <- list[n_plotsperloc] <- matrix[n_times, n_species]
  sims <- mapply(simulateOneSeriesInEnv,
                 as.data.frame(t(pars$B_loc)),
                 as.data.frame(t(pars$C_a_loc)),
                 as.data.frame(t(pars$C_b_loc)),
                 as.data.frame(t(pars$C_j_loc)),
                 as.data.frame(t(pars$G_loc)),
                 as.data.frame(t(pars$H_loc)),
                 as.data.frame(t(pars$L_loc)),
                 as.data.frame(t(pars$M_a_loc)),
                 as.data.frame(t(pars$M_j_loc)),
                 as.data.frame(t(pars$R_loc)),
                 as.data.frame(t(pars$S_loc)),
                 Times,
                 
                 SIMPLIFY = F)
  
  ## function formatSims()
  ## Only to be called within environment of simulateMultipleSeriesInEnv()
  environment(formatSims) <- environment()
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
                               sigma_process = 0.05, shape_par = 10,
                               alpha_obs = c(100, 500), sigma_obs = c(1, 0.2), phi_obs = c(10, 100),# observation error, both for a possible normal, negbinomial, and response are provided
                               obsprocess = c("gamma", "normal", "negbinomial"),
                               
                               knockoutparname = c("m_a", "m_j"),
                               
                               ## to make log parameters quickly settable
                               b_log = rnorm(n_species, c(-1, -2), 0.5),
                               c_a_log = rnorm(n_species, -2, 1),
                               c_b_log = rnorm(n_species, -2.5, 1),
                               c_j_log = rnorm(n_species, -3, 0.5),
                               g_logit = rnorm(n_species, -0.5, 0.5),
                               h_logit = rnorm(n_species, 0.2, 0.5),
                               l_log = rnorm(n_species, 3, 0.8),
                               m_a_log = rnorm(n_species, -1.5, 0.5),
                               m_j_log = rnorm(n_species, -1, 0.5),
                               r_log = rnorm(n_species, 2, 0.5),
                               s_log = rnorm(n_species, -1.6, 0.5),
                               
                               ## Other parameters can be passed ...
                               
                               ...) {
  set.seed(seed)
  
  parnames <- c("b", "c_a", "c_b", "c_j", "g", "h", "l", "m_a", "m_j", "r", "s")
  obsprocess <- match.arg(obsprocess)
  
  alpha_obs_inv <- 1/alpha_obs
  phi_obs_inv <- 1/phi_obs
  
  ## Beta_* are effect of environment on species' parameters matrix[n_env, n_species]
  ## These are on the log scale!
  n_beta <- 1 + 2 * n_env # degree 2 polynomial + intercept
  
  ### 1. log-scale parameters and env-dependence
  Vertex_b <- rbind(z = b_log,
                    p = rnorm(n_species, 0, 0.2),
                    a_1 = -rlnorm(n_species, log(0.1), 0.01),
                    q = rnorm(n_species, 0, 0.2),
                    a_2 = -rlnorm(n_species, log(0.2), 0.01))
  Beta_b <- transformVertexToNormal(Vertex_b)
  
  
  # Usually thought of as env-indepedent, but generated with zero effects here for consistency
  Vertex_c_a <- rbind(z = c_a_log,
                      p = rep(0, n_species),
                      a_1 = rep(0, n_species),
                      q = rep(0, n_species),
                      a_2 = rep(0, n_species))
  Beta_c_a <- transformVertexToNormal(Vertex_c_a)
  
  Vertex_c_b <- Vertex_c_a
  Vertex_c_b[1,] <- c_b_log
  Beta_c_b <- transformVertexToNormal(Vertex_c_b)
  
  Vertex_c_j <- rbind(z = c_j_log,
                      p = rnorm(n_species, 0.1, 0.2),
                      a_1 = rlnorm(n_species, log(0.1), 0.01),
                      q = rnorm(n_species, -0.1, 0.2),
                      a_2 = rlnorm(n_species, log(0.1), 0.01))
  Beta_c_j <- transformVertexToNormal(Vertex_c_j)
  
  
  Vertex_g <- rbind(z = g_logit,
                    p = rnorm(n_species, -0.8, 0.1),
                    a_1 = -rlnorm(n_species, log(0.8), 0.01),
                    q = rnorm(n_species, 0, 1),
                    a_2 = -rlnorm(n_species, log(0.8), 0.01))
  Beta_g <- transformVertexToNormal(Vertex_g)
  
  
  Vertex_h <- rbind(z = h_logit,
                    p = rnorm(n_species, -0.8, 0.1),
                    a_1 = -rlnorm(n_species, log(0.001), 0.01),
                    q = rnorm(n_species, 0, 1),
                    a_2 = -rlnorm(n_species, log(0.001), 0.01))
  Beta_h <- transformVertexToNormal(Vertex_g)
  
  
  Vertex_l <- rbind(z = l_log,
                    p = rnorm(n_species, 1, 1),
                    a_1 = -rlnorm(n_species, log(0.1), 0.1),
                    q = rnorm(n_species, -0.5, 0.5),
                    a_2 = -rlnorm(n_species, log(0.2), 0.2))
  Beta_l <- transformVertexToNormal(Vertex_l)
  
  Vertex_m_a <- rbind(z = m_a_log,
                      p = rnorm(n_species, 0, 1),
                      a_1 = rlnorm(n_species, log(0.1), 0.1),
                      q = rnorm(n_species, -0.2, 0.5),
                      a_2 = rlnorm(n_species, log(0.2), 0.2))
  Beta_m_a <- transformVertexToNormal(Vertex_m_a)
  
  Vertex_m_j <- rbind(z = m_j_log,
                      p = rnorm(n_species, 0, 1),
                      a_1 = rlnorm(n_species, log(2), 0.1),
                      q = rnorm(n_species, 0.2, 0.5),
                      a_2 = rlnorm(n_species, log(2), 0.2))
  Beta_m_a <- transformVertexToNormal(Vertex_m_j)
  
  Vertex_r <- rbind(z = r_log,
                    p = rnorm(n_species, 0.1, 1),
                    a_1 = -rlnorm(n_species, log(1.2), 0.2),
                    q = rnorm(n_species, -0.3, 0.5),
                    a_2 = -rlnorm(n_species, log(1.2), 0.2))
  Beta_r <- transformVertexToNormal(Vertex_r)
  
  Vertex_s <- rbind(z = s_log,
                    p = rnorm(n_species, 0, 1),
                    a_1 = rlnorm(n_species, log(1), 0.2),
                    q = rnorm(n_species, 0, 0.5),
                    a_2 = rlnorm(n_species, log(1), 0.2))
  Beta_s <- transformVertexToNormal(Vertex_s)
  
  ## old method for Beta
  # s_log <- rnorm(n_species, -1.6, 0.5)
  # Beta_s <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  # Beta_s[1,] <- s_log
  # Beta_s[3,] <- rnorm(n_species, -1, 0.2)
  
  
  ### 2. "state model scale"-parameters, positive; for env-independent models
  b <- exp(b_log)
  
  c_a <- exp(c_a_log)
  c_b <- exp(c_b_log)
  c_j <- exp(c_j_log)
  
  g <- plogis(g_logit)
  h <- plogis(h_logit)
  
  l <- exp(l_log)
  m_a <- exp(m_a_log)
  m_j <- exp(m_j_log)
  r <- exp(r_log)
  s <- exp(s_log)
  
  ## 
  radius_a_upper <- dbh_lower_b/2
  radius_a_avg <- (dbh_lower_a + dbh_lower_b)/2/2 # [mm]
  ba_a_upper <-  pi * radius_a_upper^2 * 1e-6
  ba_a_avg <- rep(pi * radius_a_avg^2 * 1e-6, times = n_species) # mm^2 to m^2
  
  pars <- c(as.list(environment()), list(...)) ## Carries along all parameter values within the named list.
  pars <- knockOutParameters(pars, knockoutparname)
  
  
  return(pars)
}




# ——————————————————————————————————————————————————————————————————————————————————— #
# Simulation helpers  ----------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

#### Get environmental variables -------------------
simulateEnv <- function(n_env, n_locs){
  n <- n_env * n_locs
  Env <- matrix(runif(n, -1, 1), nrow = n_locs, ncol = n_env)
  colnames(Env) <- paste("env", 1:n_env, sep = "_")
  return(Env)
}


#### Simulate vector of initial states -------------------
generateInitialState <- function(n_species, n_stages = 3, logstate = F) {
  s <-  rep(c(2, 0, 1), each = n_species) + rnorm(n_stages, 0, 1) # Initial state matrix.
  if(logstate) return(s) else return(exp(s))
}


#### Knockout parameters -------------------
knockOutParameters <- function(pars, parname) {
  parname_log <- c(paste0(str_to_title(parname), "_log"), paste0(parname, "_log"), paste0("Beta_", parname))
  parname <- c(parname, paste0(parname, "_log"), paste0(str_to_title(parname), "_loc"))
  
  ## knockout, but by multiplication to retain data structure (matrices, etc.)
  pars[parname] <-  lapply(pars[parname], function(p) p*0L)
  pars[parname_log] <-  lapply(pars[parname_log], function(p) p*NA)
  
  ## purge empty list items that got created
  pars <- pars[lapply(pars,length)>0]
  
  return(pars)
}


#### Transform parameters -------------------
## Accepts and returns a list of parameters. No parameter is lost.
## returns only envdependent parameters for returndf!
##
### 1. Generate environment parameters for the population model from expected values on the log scale
## Beta_* are effect of environment on species' parameters matrix[1 + n_env, n_species]
### 2. others: either a: random gamma, b: just replicate to loc dimensions.
transformParameters <- function(pars, Env,
                                envdependent = c(b = F, c_a = F, c_b = F, c_j = F, g = F, h = F, l = F, m_a  = F, m_j = F, r = F, s = F),
                                ranef = F, returndf = F, l01 = T) {
  
  
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
  
  transformLocPar <- function(L, name) {
    if (name %in% c("g", "h")) {
      return(plogis(L))
    } else return(exp(L))
  }
  
  pars_envdep <- mapply(transformLocPar, pars_envdep, parname_envdep, SIMPLIFY = F) %>%
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
  
  ## fake a smooth with l in [0,1]
  pars_envdep$L_smooth <- scale01(pars_envdep$L_loc)
  
  ## handle knockoutpars even after random effects
  pars_envdep <- knockOutParameters(pars_envdep, pars$knockoutparname)
  pars_envindep <- knockOutParameters(pars_envindep, pars$knockoutparname)
  
  ## handle data structure.
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
    
    transpars <- c(pars, pars_envdep, pars_envindep, list(envdependent = envdependent, factor_L = attr(pars_envdep$L_smooth, "factor")))
    
  }

  return(transpars)
}


#### Transform Vertex matrix to Normal Polynomial matrix -------------------
transformVertexToNormal <- function(Vertex) {
  
  ## assumes vertex form f(x,y) == a_1*(x−p)^2 + a_2*(y−q)^2 + z
  # and input matrix[N_beta, N_species] parameters [z, p, a_1, q, a_2]
  z <- Vertex[1,]
  p <- Vertex[2,]
  a_1 = Vertex[3,]
  q = Vertex[4,]
  a_2 = Vertex[5,]
  
  # and output polynomial parameters [c, b_1, a_1, b_2, a_2]
  # as in f(x, y) = a_1*x^2 + b_1*x + a_2*y^2 + b_2*y + c
  Vertex[1,] = a_2*q^2 + a_1*p^2 + z # replace z with c
  Vertex[2,] = -2*a_1*p # replace p and q with b1 and b2
  Vertex[4,] = -2*a_2*q
  
  return(Vertex)
}


#### Generate some sensible prior distribution from true parameterss --------------
generatePriors <- function(truepars, factor = 10) {
  
  ## Environmental effects vertices
  bindVertex <- function(Vertex, factor) {
    mu_Vertex <- rowMeans(Vertex)
    sigma_Vertex <- abs(mu_Vertex + log(factor))
    rbind(mu_Vertex, sigma_Vertex)
  }
  
  vertexpriors <- sapply(paste0("Vertex_", truepars$parnames), function(p) bindVertex(truepars[[p]], factor = factor),
                         USE.NAMES = T, simplify = F) %>%
    set_names(paste0("prior_", names(.)))
  
  bindGlobal <- function(global, factor) {
    rbind(mu = global, sigma = abs(global + log(factor)))
  }
  
  globalpriors <- sapply(paste0(truepars$parnames, ifelse(truepars$parnames %in% c("g", "h"), "_logit", "_log")),
                         function(p) bindGlobal(truepars[[p]], factor = factor),
                         USE.NAMES = T, simplify = F) %>%
    set_names(paste0("prior_", names(.)))

  ## Special case l
  globalpriors$prior_l_log[2,] <- 10
  
  ## purge empty list items that got created
  priors <- c(vertexpriors, globalpriors)
  priors <- priors[!sapply(priors, anyNA)]
  
  return(priors)
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
  
  
  if (format == "long") { } # Sims <- Sims
  
  ######### (1) ##########
  else if (format == "stan"){
    
    #### Special long data set structures with different amounts of long casting along different ids
    ## The most comprehensive data set is L_y (resp. L_y0) with grouping locations/resurveys/pops/plots.
    ## NOTE THAT plots HAS TO BE THE LAST GROUP IN SORTING
    ## Everything is subset from the master subsets with these groupings (*_reobs, *_y0) and thus consistently sorted.
    ## Ratio: "resurveys" includes "pops" due of the return structure in ode_*(); "pops" includes "plots" because of the loc == population assumption (randowm effect from location/resurvey/pop to .../plots).
    ## Factor "pops" is structured stages/species.
    
    ## scale times within group to start with 1
    Sims %<>%
      group_by(loc) %>%
      mutate(time = time - min(time) + 1) %>% 
      ungroup() %>%
      arrange(loc, time, pop, stage, species, plot)
    
    ## Format: [L_y0] —   locations/pops(/stage/species)/plots
    Sims_y0 <- filter(Sims, isy0) %>%
      select(-isy0) %>%
      arrange(loc, time, pop, stage, species, plot)
    
    ## Format: [L_y_reobs] — locations/resurveys/pops(/stage/species)/plots
    Sims_reobs <- filter(Sims, !isy0) %>%
      select(-isy0) %>%
      arrange(loc, time, pop, stage, species, plot)
    
    ## Format: [L_y] — locations/obs/pops(/stage/species)/plots
    Sims <- select(Sims, -isy0)
    
    ## Format: [L_init] — locations/pops
    Sims_init <- Sims_y0 %>%
      group_by(loc, pop, stage, species) %>%
      summarize(n_plots = n_distinct(plot), .groups = "drop")
    
    ## Format: [L_yhat_reobs] — locations/resurveys/pops
    Sims_yhat_reobs <- Sims_reobs %>%
      group_by(loc, time, pop, stage, species) %>%
      summarize(n_plots = n_distinct(plot), .groups = "drop")
    
    ## Format: [L_yhat] — locations/resurveys/pops
    Sims_yhat <- Sims %>%
      group_by(loc, time, pop, stage, species) %>%
      summarize(n_plots = n_distinct(plot), .groups = "drop")
    
    ## Format: [L_times_reobs] — locations/resurveys
    Sims_times_reobs <- Sims_reobs %>% # reobs to  count all times
      group_by(loc) %>%
      summarize(time = unique(time), .groups = "drop")
    
    ## Format: [L_times] — locations/resurveys
    Sims_times <- Sims %>%
      group_by(loc) %>%
      summarize(time = unique(time), .groups = "drop") ## DANGER!
    
    ## Format: [N_locs] — locations
    Sims_locs_reobs <- Sims_reobs %>% # reobs to  count all times
      group_by(loc) %>%
      ## assumes completion within locations!
      summarize(time_max = max(time),
                timespan = diff(range(time)) + 1,
                n_species = n_distinct(species),
                n_plots = n_distinct(plot),
                n_pops = n_distinct(pop),
                n_reobs = n_distinct(time),
                n_yhat = n_distinct(interaction(pop, time)), .groups = "drop")
    
    ## Format: [N_locs] — locations
    Sims_locs <- Sims %>%
      group_by(loc) %>%
      ## assumes completion within locations!
      summarize(time_max = max(time),
                timespan = diff(range(time)) + 1,
                n_species = n_distinct(species),
                n_plots = n_distinct(plot),
                n_pops = n_distinct(pop),
                n_obs = n_distinct(time),
                n_yhat = n_distinct(interaction(pop, time)), .groups = "drop")
    
    
    ## Prepare design matrix
    polyformula <- as.formula(paste("~", paste("poly(", colnames(as.data.frame(Env)), ", 2)", collapse = "+")))
    X <- model.matrix(polyformula, data = as.data.frame(Env))
    
    ## prepare rect sims, which might also be used in others
    ## Go back to the original sim list for rectangular casting
    time_rect <- sims[[1]][[1]][,1] # get the time column  from the very first table
    sims_rect <- lapply(sims, function(loc) lapply(loc, function (plot) as.data.frame(t(plot)))) # make everything into lists!
    sims_rect <- array(unlist(sims_rect),
                       dim = c(length(sims_rect[[1]][[1]][[1]]), length(sims_rect[[1]][[1]]), length(sims_rect[[1]]), length(sims_rect))
    ) 
    sims_rect <- aperm(sims_rect, perm = c(4, 2, 3, 1))
    # str(S) ## # 'loc', 'time', 'plot', 'pop'
    sims_rect <- sims_rect[ , , , -1] # drop the first 'column', which is not a population but the times
    
    
    ######### (1a) ##########
    if (isrectangular) {
      

      
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
                       
                       ba_a_upper = meta$ba_a_upper,
                       ba_a_avg = meta$ba_a_avg,
                       
                       times = t(replicate(pars$n_locs, times)) - min(times) + 1,
                       time_max = Times[,ncol(Times)],
                       timespan_max = max(Times[,ncol(Times)]) - 1, # difference to the first time
                       
                       X = X,
                       
                       y0_loc = apply(sims_rect, c(1, 2, 4), mean)[, 1, ],
                       y0_loc_log = apply(exp(sims_rect), c(1, 2, 4), mean)[, 1, ] %>% log(), # only for cases with log response
                       
                       y = sims_rect
      )
      
    }
    
    ######### (1b) ##########
    else if (modelstructure %in% c("ba-rag-ranef", "ba-rag", "ba")) {
      
      ## commented out are former versions with separate fitting of initials
      
      Sims_species <- Sims_init[c("loc", "species")] %>% unique()
      L_init <- nrow(Sims_init)
      L_times <- nrow(Sims_times)
      L_yhat <- nrow(Sims_yhat)
      
      N_locs <-  nrow(Sims_locs)
      time_init <-  Sims$time[match(Sims_locs$loc, Sims$loc)] # N_locs!
      
      vrep <- function(...) unlist( c( Vectorize(rep.int)(...) ) ) # :)))))
      
      stansims <- list(
        N_locs = N_locs,
        N_plots = pars$n_plotsperloc,
        
        L_init = L_init, # !!!
        L_times = L_times,
        L_yhat = L_yhat,
        L_y0 = nrow(Sims_y0),
        L_y = nrow(Sims),
        L_species = nrow(Sims_species), 
        
        N_species = length(unique(Sims$species)),
        N_pops = length(unique(Sims$pop)),
        N_beta = ncol(X),
        
        n_species = Sims_locs$n_species,
        n_pops = Sims_locs$n_pops,
        n_obs = Sims_locs$n_obs,
        n_reobs = Sims_locs_reobs$n_reobs,
        n_yhat = Sims_locs$n_yhat,
        
        i_j = 1:pars$n_species,
        i_a = (1:pars$n_species) + pars$n_species,
        i_b = (1:pars$n_species) + 2*pars$n_species,
        
        ## repeat predictions on level "locations/pops" n_plots times to "locations/pops/plots"
        rep_init2y0 = vrep(1:L_init, Sims_init$n_plots),
        ## repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
        rep_yhat2y = vrep(1:L_yhat, Sims_yhat$n_plots),
        ## repeat locations on level "locations" n_reobs times to "locations/pops/resurveys"
        rep_locs2times = vrep(1:N_locs, Sims_locs$n_obs),
        rep_locs2times_reobs = vrep(1:N_locs, Sims_locs_reobs$n_reobs),
        
        rep_obsmethod2y0 = Sims_y0$obsmethod,
        rep_obsmethod2y = Sims$obsmethod,
        rep_obsmethod2yreobs = Sims_reobs$obsmethod,
        
        ba_a_upper = meta$ba_a_upper,
        ba_a_avg = meta$ba_a_avg,
        
        species = Sims_species$species,
        stage = Sims_init$stage,
        pops = Sims_init$pop,
        
        time_init = time_init,
        
        times = Sims_times$time,
        time_max = Sims_locs$time_max,
        timespan_max = max(Sims_locs$timespan),
        
        X = X,
        
        y0 = Sims_y0$abundance,
        y = Sims$abundance,
        
        ## rectangular inits for parameter recovery
        state_init = apply(sims_rect, c(1, 2, 4), mean)[, 1, ],
        state_init_log = apply(exp(sims_rect), c(1, 2, 4), mean)[, 1, ] %>% log() # only for cases with log response

      )
      
    }
    
    ### Some "settings" for the model
    ## if (format == "stan")
    stansims$tolerance_fix <- 0.001
    stansims$L_loc <- pars$L_loc # get L_log, from transformed parameters to represent predicted dispersal prob
    stansims$L_loc_log <- pars$L_log # get L_log, from transformed parameters to represent predicted dispersal prob
    stansims$L_smooth <- pars$L_smooth
    
    ## Get normal priors from parameters with sigma = mu + log(factor)
    stansims <- c(stansims, generatePriors(pars, factor = priorfactor))
    
    stansims$Long <- Sims
    Sims <- stansims
    
  } ## closing if (format == "stan")
  
  return(Sims)
}



####  Do the fit ----------------------
drawSamples <- function(model, data,
                        method = c("variational", "mcmc", "sim"), n_chains = 3, initfunc = 0,
                        iter_warmup = 400, iter_sampling = 600,
                        dirpath = "Fits.nosync", ...) {
  
  if (!dir.exists(dirpath)) {
    dir.create(dirpath)
  }
  
  if(match.arg(method) == "variational") {
    fit <- model$variational(data = data,
                             output_dir = dirpath,
                             init = initfunc,
                             eta = 0.001,
                             iter = 20**4,
                             ...) # convergence after 1500 iterations
    
  } else if (match.arg(method) == "mcmc") {
    
    fit <- model$sample(data = data,
                        output_dir = dirpath,
                        # output_basename = ,
                        init = initfunc,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = 0.99,
                        max_treedepth = 16,
                        chains = n_chains, parallel_chains = getOption("mc.cores", n_chains),
                        ...)
    
  } else if (match.arg(method) == "sim") {
    
    fit <- model$sample(data = data,
                        output_dir = NULL,
                        init = initfunc,
                        iter_sampling = 500,
                        fixed_param = TRUE,
                        ...)
  }
  
  return(fit)
}



# ——————————————————————————————————————————————————————————————————————————————————— #
# General helpers  ----------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

repArray <- function(A, n) {
  dims <- if(is.array(A)) length(dim(A)) else 1 # else arm handles vectors
  B <- aperm(replicate(n, A), perm = c(dims+1, 1:dims))
  return(B)
}


# ——————————————————————————————————————————————————————————————————————————————————— #
# Stats helpers  ----------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

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


scale01 <- function(x) {
  x <- x - min(x)
  factor_x <- 1/max(x)
  x <- x * factor_x
  attr(x, "factor") <- factor_x
  return(x)
}




# ——————————————————————————————————————————————————————————————————————————————————— #
# Alternaitve model formulation  ------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

calcModel_ricker <- function(times,
                             initialstate_log, # A vector of species states.
                             pars # internal number of discrete time steps within a unit of time
) { 
  
  
  ## Just unpacking for readability.
  
  b <- pars$b # Length n_species vector of basal area increment rates.
  c_a <- pars$c_a # Length n_species vector
  c_b <- pars$c_b # Length n_species vector
  c_j <- pars$c_j # Length n_species vector
  g <- pars$g # Length n_species vector of transition rates.
  h <- pars$h # Length n_species vector of transition rates.
  l <- pars$l 
  m_a <- pars$m_a # Length n_species vector of mortalities
  m_j <- pars$m_j # Length n_species vector of mortalities
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  
  ba_a_upper <- pars$ba_a_upper
  ba_a_avg <- pars$ba_a_avg
  
  sigma <- pars$sigma_process
  
  ## Set the count variables
  n <- length(r) # no species
  
  ## 
  times_intern <- 1:max(times)
  n_times <- length(times_intern)
  
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
    BA <- A * ba_a_avg + B
    BA_sum <- sum(BA)
    
    # observation error for stages/times/species
    u <- matrix(0, nrow = length(s), ncol = 3) + rnorm(length(s)*3, 0, sigma)
    
    ## Two kinds of processes acting ot State_log
    ## 1. All processes that add additively to State, are added within log(State)
    ## 2. All processes that act multiplicatively on the state are added in log space
    
    J_trans <- A + r*BA + l + (J - J*g)
    J_t_log <- log(J_trans) + (1 - c_j*sum(J_t_add) - s*BA_sum - m_j)  + u[ ,1] # count of juveniles J
    
    A_trans <- A + J*g + (A - A*h)
    A_t_log <- log(A_trans) + (1 - c_b*BA_sum - m_a) + u[ ,2] # count of small adults A
    
    A_ba <-  A * h * ba_a_upper # Basal area of small adults A. Conversion by multiplication with basal area of State exit (based on upper dhh boundary of the class)
    B_trans <- B + A_ba
    B_t_log <- log(B_trans) + b - c_b*BA_sum + u[ ,3] # basal area of big adults B
    ## b is the net basal area increment (including density-independent m) basically equivalent to a Ricker model, i.e. constant increment rate leading to exponential growth, negative density dependent limitation scaled with total BA_sum.
    
    State_log[t, ] <- c(J_t_log, A_t_log, B_t_log)
  }
  
  whichtimes <- which(times_intern %in% times)
  return(State_log[whichtimes,])
}


