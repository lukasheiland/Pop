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
rstan_options(javascript = FALSE)

library(bayesplot)


# Orientation -------------------------------------------------------------
setwd(here())
modeldir <- dir(pattern = glue("^(Model).*ba$"))

source(glue("{modeldir}/Sim_ba_helpers.R"))


######################################################################################
# Model formulation  ----------------------------------------------------------------
######################################################################################

#### Returns log states for times -------------------
calcModel <- function(times,
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
  m_a <- pars$m_a # Length n_species vector of mortalities
  m_j <- pars$m_j # Length n_species vector of mortalities
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  
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
    
    # observation error for stages/times/species
    u <- matrix(0, nrow = length(s), ncol = 3) + rnorm(length(s)*3, 0, sigma)
    
    ## Two kinds of processes acting ot State_log
    ## 1. All processes that add additively to State, are added within log(State)
    ## 2. All processes that act multiplicatively on the state are added in log space
    
    J_t_log <- log(J + r) - c_j*sum(J) - s*BA - m_j - g + u[ ,1] # count of juveniles J
    
    A_t_log <- log(A + exp(J_log - g)) - c_a*BA - m_a - h + u[ ,2] # count of small adults A
    A_ba <- exp(A_log - h)*ba_a_upper # Basal area of small adults A. Conversion by multiplication with basal area of State exit (based on upper dhh boundary of the class)
    
    B_t_log <- log(B + A_ba) + b - c_b*sum(B) + u[ ,3] # basal area of big adults B
    ## b is the net basal area increment (including density-independent m) basically equivalent to a Ricker model, i.e. constant increment rate leading to exponential growth, negative density dependent limitation scaled with total BA.
    
    State_log[t, ] <- c(J_t_log, A_t_log, B_t_log)
  }
  
  whichtimes <- which(times_intern %in% times)
  return(State_log[whichtimes,])
}



######################################################################################
# Simulation functions ---------------------------------------------------------------
######################################################################################


#### Simulate one time series for a specific fixed set of parameters. -------------------
simulateOneSeries <- function(state_init, times, pars, processerror = T, obserror = T, log = T) {
  
  
  n_pops <- length(state_init)
  
  if(!processerror) pars$sigma_process <- 0
  
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
                                        times = 2:22, # internal model will still start at 1
                                        envdependent = c(b = F, c_a = F, c_b = F, c_j = F, g = F, h = F, m_a  = F, m_j = F, r = F, s = F),
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
  pars <- transformParameters(pars, Env, envdependent = envdependent, ranef = ranef)
  meta <- pars # for info like n_stages, dbh_lower_a
  
  #### Internal function to mapply over locations with different environments and therefore possibly different parameters
  simulateOneSeriesInEnv <- function(b_, c_a_, c_b_, c_j_, g_, h_, m_a_, m_j_, r_, s_) {
    
    
    ## list of env-dependent and -independent parameters for one environment, replacing the original global parameters in "pars"
    simpars <- within(pars, {
      b <- b_;
      c_a <- c_a_; c_b <- c_b_; c_j <- c_j_;
      g <- g_; h <- h_;
      m_a <- m_a_; m_j <- m_j_;
      r <- r_;
      s <- s_;})
    
    
    ## Switch whether init is the sane for all plots within cluster
    if (plotlevelinit) {
      
      seriesatloc <- replicate(pars$n_plotsperloc,
                               simulateOneSeries(generateInitialState(n_species = pars$n_species), times, simpars,
                                                 processerror = processerror, obserror = obserror), simplify = F)
    } else {
      
      ## cluster level init
      # Replication is post initial state generation
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
                 as.data.frame(t(pars$B_loc)),
                 as.data.frame(t(pars$C_a_loc)),
                 as.data.frame(t(pars$C_b_loc)),
                 as.data.frame(t(pars$C_j_loc)),
                 as.data.frame(t(pars$G_loc)),
                 as.data.frame(t(pars$H_loc)),
                 as.data.frame(t(pars$M_a_loc)),
                 as.data.frame(t(pars$M_j_loc)),
                 as.data.frame(t(pars$R_loc)),
                 as.data.frame(t(pars$S_loc)),
                 
                 SIMPLIFY = F)
  
  ## function formatSims() in Sim_ba_helpers.R
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
                               sigma_process = 0.05, shape_par = 10, sigma_obs = c(1, 0.2),
                               ...) {
  set.seed(seed)
  
  parnames <- c("b", "c_a", "c_b", "c_j", "g", "h", "m_a", "m_j", "r", "s")
  
  ## Beta_* are effect of environment on species' parameters matrix[n_env, n_species]
  ## These are on the log scale!
  n_beta <- 1 + 2 * n_env # degree 2 polynomial + intercept
  
  ### 1. log-scale parameters and env-dependence
  
  b_log <- rnorm(n_species, -1, 0.01)
  Beta_b <- matrix(rnorm(n_beta*n_species, -0.5, 0.5), n_beta, n_species)
  Beta_b[1,] <- b_log
  Beta_b[3,] <- rnorm(n_species, -1, 0.2)
  
  # Usually thought of as env-indepedent, but generated with zero effects here for consistency
  Beta_c_a <- Beta_c_b <- Beta_c_j <- matrix(0, n_beta, n_species)
  c_a_log <- rnorm(n_species, -3.3, 0.2)
  Beta_c_a[1,] <- c_a_log
  c_b_log <- rnorm(n_species, -3, 0.1)
  Beta_c_b[1,] <- c_b_log
  c_j_log <- rnorm(n_species, -10, 0.5)
  Beta_c_j[1,] <- c_j_log
  
  g_log <- rnorm(n_species, -5, 0.2)
  Beta_g <- matrix(rnorm(n_beta*n_species, -0.5, 0.5), n_beta, n_species)
  Beta_g[1,] <- g_log
  Beta_g[3,] <- rnorm(n_species, -1, 0.2)
  
  h_log <- rnorm(n_species, -0.5, 0.3)
  Beta_h <- matrix(rnorm(n_beta*n_species, -0.5, 0.5), n_beta, n_species)
  Beta_h[1,] <- h_log
  Beta_h[3,] <- rnorm(n_species, -1, 0.2)
  
  m_a_log <- rnorm(n_species, -1.5, 0.5)
  Beta_m_a <- matrix(rnorm(n_beta*n_species, -2, 0.5), n_beta, n_species)
  Beta_m_a[1,] <- m_a_log
  Beta_m_a[3,] <- rnorm(n_species, -2, 0.2)
  
  m_j_log <- rnorm(n_species, -1.2, 0.1)
  Beta_m_j <- matrix(rnorm(n_beta*n_species, -3, 0.5), n_beta, n_species)
  Beta_m_j[1,] <- m_j_log
  Beta_m_j[3,] <- rnorm(n_species, -2.5, 0.2)
  
  r_log <- rnorm(n_species, 2.2, 0.2)
  Beta_r <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  Beta_r[1,] <- r_log
  Beta_r[3,] <- rnorm(n_species, -0.5, 0.2)
  
  s_log <- rnorm(n_species, -2.9, 0.1)
  Beta_s <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  Beta_s[1,] <- s_log
  Beta_s[3,] <- rnorm(n_species, -1, 0.2)
  
  
  ### 2. "state model scale"-parameters, positive; for env-independent models
  b <- exp(b_log)
  
  c_a <- exp(c_a_log)
  c_b <- exp(c_b_log)
  c_j <- exp(c_j_log)
  
  g <- exp(g_log)
  h <- exp(h_log)
  m_a <- exp(m_a_log)
  m_j <- exp(m_j_log)
  r <- exp(r_log)
  s <- exp(s_log)
  
  
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
times1 <- 1:15

Sim1 <- simulateOneSeries(initialstate1, times = times1, pars = pars1,
                          processerror = F, obserror = F, log = F)

matplot(Sim1[,-1], type = "b", ylab="N",
        pch = rep(c("J", "A", "B"), each = pars1$n_species),
        col = 1:pars1$n_species, xlab = "time") # log='y'



## Multiple time series in env ---------------------------------------------------------------
envdep_demo <- c(b = F, c_a = F, c_b = F, c_j = F, g = T, h = F, m_a = F, m_j = F, r = T, s = T)
pars_demo <- generateParameters(n_species = 4, n_locs = 50, sigma_obs = c(0.2, 0.1))
E_demo <- simulateEnv(pars_demo$n_env, pars_demo$n_locs)
P_env_demo <- transformParameters(pars_demo, E_demo, envdep_demo, ranef = F, returndf = T)
S_demo <- simulateMultipleSeriesInEnv(pars_demo, E_demo, times = 2:22,
                                      envdependent = envdep_demo, ranef = F, processerror = T, obserror = T)

#### Plot time series
S_demo %>%
  # filter(stage %in% c("a")) %>%
  ggplot(mapping = aes(x = time, y = abundance, color = env_1, lty = stage, group = interaction(loc, stage, plot))) +
  geom_line() +
  facet_wrap(facets = c("species"))


#### Plot parameters
P_env_demo %>%
  filter(parameter %in% c("s_loc")) %>%
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
fitseed <- 1

pars <- generateParameters(seed = fitseed, n_locs = 70)

Env <- simulateEnv(n_env = pars$n_env, n_locs = pars$n_locs)

data <- simulateMultipleSeriesInEnv(pars, Env, times = c(2, 5, 8),
                                    envdependent = c(b = F, c_a = F, c_b = F, c_j = F, g = F, h = F, m_a = F, m_j = F, r = F, s = F),
                                    modelstructure = modelname, format = "stan",
                                    obserror = T, processerror = F)

Data_long <- simulateMultipleSeriesInEnv(pars, Env, times = c(2, 5, 8),
                                         envdependent = c(b = F, c_a = F, c_b = F, c_j = F, g = F, h = F, m_a = F, m_j = F,  r = F, s = F),
                                         modelstructure = modelname, format = "long",
                                         obserror = T, processerror = F) %>%
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
  
  truepars <- attr(data, "pars")
  truepars <- truepars[sapply(truepars, is.numeric)]
  parnames <- names(truepars)
  names(truepars) <- ifelse(str_ends(parnames, "_loc"), str_to_lower(parnames), parnames)
  
  newpars <- list(
    state_init_log = if (isragged) rnorm(data$y0_log, data$y0_log, 0.01) else
                                   data$y_log[,1,,] + rnorm(data$y_log[,1,,], 0, 0.01),
    
    u = replicate(pars$n_locs, matrix(rnorm(pars$n_species*3, 0, 0.001), nrow = pars$n_species, ncol = data$timespan_max))
  )
  
  inits <- c(newpars, truepars)
  truepars <<- truepars
  
  return(inits)
}

#### Returns viable start values ---------------------------
getInits <- function() {
  
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
    
    state_init_log = if (isragged) rnorm(data$y0_log, data$y0_log, 0.01) else
      data$y_log[,,1,] + rnorm(data$y_log[,,1,], 0, 0.01),

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
    
    u = replicate(pars$n_locs, matrix(rnorm(pars$n_species*3, 0, 0.001), nrow = pars$n_species, ncol = data$timespan_max))
  )
  
  return(inits)
}


## Draw from model --------------------------------------------------------------

####  Do the fit ----------------------
drawSamples <- function(model, data, method = c("variational", "mcmc", "sim"), n_chains = 6, initfunc = 0) {
  
  if(match.arg(method) == "variational") {
    fit <- model$variational(data = data,
                             output_dir = "Fits.nosync",
                             init = initfunc,
                             iter = 20**4) # convergence after 1500 iterations
    
    } else if (match.arg(method) == "mcmc") {
    
    fit <- model$sample(data = data,
                        output_dir = "Fits.nosync",
                        init = initfunc,
                        iter_warmup = 1500, iter_sampling = 500,
                        chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
    
    } else if (match.arg(method) == "sim") {
    
      fit <- model$sample(data = data,
                          output_dir = NULL,
                          init = initfunc,
                          iter_sampling = 200,
                          fixed_param = TRUE)
  }
  return(fit)
}

# model <- cmdstan_model(modelpath)

## Model fit
fit <- drawSamples(model, data, method = "variational", initfunc = getInits)

## Other diagnostics
# fit$output()
# fit$time()
# fit$init()
# fit$draws(variables = NULL, inc_warmup = F)
# fit$return_codes()

## Simulate with fixed parameters for debugging
# fit <- drawSamples(model, data, method = "sim", initfunc = getTrueInits)


## Extract draws -----------------------------------------------------------

drawpath <- fit$output_files() # dput(fit$output_files())
    ## alternatively
    # drawpath <- file.info(list.files("Fits", full.names = T)) %>%
    #   arrange(desc(mtime)) %>%
    #   slice(1:6) %>%
    #   rownames()


# drawnames <- unique(str_extract(s$variable, "[a-zA-Z_]*"))
# targetvarname <- intersect(drawnames, names(truepars)) # this drops boilerplate variables like u
# targetvarname <- c(targetvarname, "state_init_log") # "y_hat_log" "sim_log", "y_hat_log"
excludevarname <- c("u") # state_init_log

### Get draws
## rstan format for use in other packages
stanfit <- rstan::read_stan_csv(drawpath)
draws <- extract(stanfit, pars = c(excludevarname), include = F)
    ## get draws with cmdstanr
    # draws <- read_cmdstan_csv(drawpath, variables = NULL)



## Summary -----------------------------------------------------------------

s <- summary(stanfit) # , pars = excludevarname, include = F
s


# Inspection --------------------------------------------------------------
truepars <- attr(data, "pars")


#### Predictions vs. true --------------------

plot(c(draws$sim_log[1,,,,,drop =T]) ~ c(data$y_log))
plot(c(draws$y_hat_log[1,,,,drop =T]) ~ c(apply(data$y_log, c(1, 2, 4), mean))) # average the plots away
abline(0, 1)


#### Draws vs true ----------------------
plotDrawVsSim <- function(parname = "h",
                          simdata = data,
                          rstandraws = stanfit) {
  
  simpar <- attr(simdata, "pars")[[parname]]
  
  if(is.vector(simpar)) {
    
    Stanpar <- rstan::extract(rstandraws, pars = parname)[[1]] %>%
      as.data.frame() %>%
      pivot_longer(cols = everything(), names_to = "species", values_to = "draw") %>%
      mutate(species = as.integer(as.factor(species))) %>%
      bind_cols(true = rep(simpar, length.out = nrow(.)))
    
    Stanpar %>%
      # sample_n(1000) %>%
      gf_point(draw ~ true, alpha = 0.3, size = 0.1) %>%
      gf_abline(slope = 1, intercept = 0, gformula = NULL)
    
  } else if (str_starts(parname, "Beta")) {
    
    print("Not yet implemented.")
    
  } else {
    
    # Stanpar <- rstan::extract(rstandraws, pars = parname)[[1]] %>%
    #   as.data.frame() %>%
    #   pivot_longer(cols = everything(), names_to = "species", values_to = "draw") %>%
    #   mutate(species = as.integer(as.factor(species))) %>%
    #   bind_cols(true = rep(simpar, length.out = nrow(.)))
    # 
    # Stanpar %>%
    #   # sample_n(1000) %>%
    #   gf_point(draw ~ true | species, alpha = 0.3, size = 0.1) %>%
    #   gf_abline(slope = 1, intercept = 0, gformula = NULL)
    
  }
}


## Fitted parameters vs. true
plotDrawVsSim("b")
plotDrawVsSim("c_j")
plotDrawVsSim("c_a")
plotDrawVsSim("c_b")
plotDrawVsSim("h")
plotDrawVsSim("m_a")
plotDrawVsSim("g")
plotDrawVsSim("m_j")
plotDrawVsSim("r")
plotDrawVsSim("s")

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

