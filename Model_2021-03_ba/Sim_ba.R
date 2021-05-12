# Orientation -------------------------------------------------------------
library(here)

setwd(here())
modeldir <- dir(pattern = glue("^(Model).*ba$"))

source(file.path(modeldir, "Sim_ba_helpers.R"))

# startsource —————————————————————————————————————————————————————————————
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



######################################################################################
# Model formulation  ----------------------------------------------------------------
######################################################################################


#### Returns states for times -------------------
calcModel <- function(times,
                      initialstate, # A vector of species states.
                      pars # internal number of discrete time steps within a unit of time
                      ) {
  
  ## Just unpacking for readability.
  
  b <- pars$b # Length n_species vector of basal area increment rates.
  c_a <- pars$c_a # Length n_species vector
  c_b <- pars$c_b # Length n_species vector
  c_j <- pars$c_j # Length n_species vector
  g <- pars$g # Length n_species vector of transition rates.
  h <- pars$h # Length n_species vector of transition rates.
  l <- pars$l # Length n_species vector of input rates.
  m_a <- pars$m_a # Length n_species vector of mortalities
  m_j <- pars$m_j # Length n_species vector of mortalities
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  
  dbh_lower_a <- pars$dbh_lower_a
  dbh_lower_b <- pars$dbh_lower_b
  
  sigma <- pars$sigma_process
  
  ## Set the count variables
  n <- length(r) # no. of species
  
  ## 
  times_intern <- 1:max(times)
  n_times <- length(times_intern)
  
  radius_a_upper <- dbh_lower_b/2
  radius_a_avg <- (dbh_lower_a + dbh_lower_b)/2/2 # [mm]
  ba_a_upper <-  pi * radius_a_upper^2 * 1e-6
  ba_a_avg <- pi * radius_a_avg^2 * 1e-6 # mm^2 to m^2
  
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
    A_t <-  A_trans * 1/(1 + c_b*BA_sum + m_a) # count of small adults
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


## Wrapper for choosing a model and hard-resetting parameters
calcModelWrapper <- function(times, initialstate, pars, ...) { 
  
  calcModel(times, initialstate,
            within(pars, { m_a <- 0;  m_j <- 0; l <- 0;}), # 
            ...)
  }


######################################################################################
# Simulation functions ---------------------------------------------------------------
######################################################################################


#### Simulate one time series for a specific fixed set of parameters. -------------------
simulateOneSeries <- function(state_init, times, pars,
                              processerror = T, obserror = T, exping = F, nominaltimes = NULL) {
  
  
  n_pops <- length(state_init)
  
  if(!processerror) pars$sigma_process <- 0
  
  Sim <- calcModelWrapper(state_init, times = times, pars = pars)
  
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
                                        format = c("long", "stan"),
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
                               sigma_process = 0.05, shape_par = 10,
                               alpha_obs = c(100, 500), sigma_obs = c(1, 0.2), phi_obs = c(10, 100),# observation error, both for a possible normal, negbinomial, and response are provided
                               obsprocess = c("gamma", "normal", "negbinomial"),
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
  
  b_log <- rnorm(n_species, c(-1, -2), 0.5)
  Beta_b <- matrix(rnorm(n_beta*n_species, -0.5, 0.5), n_beta, n_species)
  Beta_b[1,] <- b_log
  Beta_b[3,] <- rnorm(n_species, -1, 0.2)
  
  # Usually thought of as env-indepedent, but generated with zero effects here for consistency
  Beta_c_a <- Beta_c_b <- Beta_c_j <- matrix(0, n_beta, n_species)
  c_a_log <- rnorm(n_species, -2, 1)
  Beta_c_a[1,] <- c_a_log
  c_b_log <- rnorm(n_species, -2.5, 1)
  Beta_c_b[1,] <- c_b_log
  c_j_log <- rnorm(n_species, -3, 2)
  Beta_c_j[1,] <- c_j_log
  
  g_logit <- rnorm(n_species, c(-0.2, -1), 1)
  Beta_g <- matrix(rnorm(n_beta*n_species, -0.3, 0.2), n_beta, n_species)
  Beta_g[1,] <- g_logit
  Beta_g[3,] <- rnorm(n_species, -0.4, 0.2)
  
  h_logit <- rnorm(n_species, 0.2, 0.5)
  Beta_h <- matrix(rnorm(n_beta*n_species, -0.3, 0.2), n_beta, n_species)
  Beta_h[1,] <- h_logit
  Beta_h[3,] <- rnorm(n_species, -0.4, 0.2)
  
  l_log <- rnorm(n_species, 3, 0.8)
  Beta_l <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  Beta_l[1,] <- l_log
  Beta_l[3,] <- rnorm(n_species, -0.5, 0.2)
  
  m_a_log <- rnorm(n_species, -1.5, 0.5)
  Beta_m_a <- matrix(rnorm(n_beta*n_species, -2, 0.5), n_beta, n_species)
  Beta_m_a[1,] <- m_a_log
  Beta_m_a[3,] <- rnorm(n_species, -2, 0.2)
  
  m_j_log <- rnorm(n_species, -1.2, 0.1)
  Beta_m_j <- matrix(rnorm(n_beta*n_species, -3, 0.5), n_beta, n_species)
  Beta_m_j[1,] <- m_j_log
  Beta_m_j[3,] <- rnorm(n_species, -2.5, 0.2)
  
  r_log <- rnorm(n_species, c(3, 1), 0.5)
  Beta_r <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  Beta_r[1,] <- r_log
  Beta_r[3,] <- rnorm(n_species, -0.5, 0.2)
  
  s_log <- rnorm(n_species, -1.6, 0.5)
  Beta_s <- matrix(rnorm(n_beta*n_species, 0, 0.5), n_beta, n_species)
  Beta_s[1,] <- s_log
  Beta_s[3,] <- rnorm(n_species, -1, 0.2)
  
  
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
  
  
  pars <- c(as.list(environment()), list(...)) ## Carries along all parameter values within the named list.
  return(pars)
}


# endsource —————————————————————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————————————————————————————————— #
# Demo Simulations ------------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## One time series ---------------------------------------------------------------

## Direct plotting
# matplot(exp(calcModel(times1, exp(initialstate1), pars1)[times1[-c(1:3)],]), type = "b")

#### prerequisites for a single time series
formals(generateParameters)
pars1 <- generateParameters(seed = 2, n_species = 2, obsprocess = "gamma")
initialstate1 <- generateInitialState(n_species = pars1$n_species)
times1 <- 2:40

Sim1 <- simulateOneSeries(initialstate1, times = times1, pars = pars1,
                          processerror = F, obserror = F)

matplot(Sim1[,-1], type = "b", ylab="N",
        pch = rep(c("J", "A", "B"), each = pars1$n_species),
        col = 1:pars1$n_species, xlab = "time") # log='y'



## Multiple time series in env ---------------------------------------------------------------
envdep_demo <- c(b = F, c_a = F, c_b = F, c_j = F, g = T, h = F, l = T, m_a = F, m_j = F, r = F, s = T)
pars_demo <- generateParameters(n_species = 3, n_locs = 50)
E_demo <- simulateEnv(pars_demo$n_env, pars_demo$n_locs)
P_env_demo <- transformParameters(pars_demo, E_demo, envdep_demo, ranef = F, returndf = T)
S_demo <- simulateMultipleSeriesInEnv(pars_demo, E_demo, times = 2:100,
                                      envdependent = envdep_demo, ranef = F, processerror = F, obserror = F, independentstart = F)

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

ggplot(filter(P_env_demo, parameter == "g_loc"),
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
parseed <- 2

pars <- generateParameters(seed = parseed, n_locs = 100, n_species = 2, n_plotsperloc = 4)

Env <- simulateEnv(n_env = pars$n_env, n_locs = pars$n_locs)


envdependent_ba_rect <- c(b = F, c_a = F, c_b = F, c_j = F, g = F, h = F, l = T, m_a = F, m_j = F, r = F, s = F)
envdependent_ba_rag <- c(b = F, c_a = F, c_b = F, c_j = F, g = T, h = F, l = T, m_a = F, m_j = F, r = T, s = T)
envdependent_ba <- c(b = F, c_a = F, c_b = F, c_j = F, g = T, h = F, l = T, m_a = F, m_j = F, r = T, s = T)
envdependent_ba_rag_ranef <- envdependent_ba_rag

data <- simulateMultipleSeriesInEnv(pars, Env, times = c(2, 8, 15),
                                    envdependent = get(paste0("envdependent_", sub("-", "_", modelname))),
                                    logstate = F,
                                    modelstructure = modelname, # !!! this determines the data layout
                                    format = "stan",
                                    obserror = T, processerror = F, independentstart = T)


Data_long <- simulateMultipleSeriesInEnv(pars, Env, times = c(2, 8, 15),
                                         envdependent = c(b = F, c_a = F, c_b = F, c_j = F, g = F, h = F, l = T, m_a = F, m_j = F, r = F, s = F),
                                         modelstructure = modelname, format = "long",
                                         obserror = F, processerror = F, independentstart = T) %>%
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





## Draw from model --------------------------------------------------------------

####  Do the fit ----------------------
drawSamples <- function(model, data, method = c("variational", "mcmc", "sim"), n_chains = 6, initfunc = getTrueInits) {
  
  if(match.arg(method) == "variational") {
    fit <- model$variational(data = data,
                             output_dir = "Fits.nosync",
                             init = initfunc,
                             eta = 0.001,
                             iter = 20**4) # convergence after 1500 iterations
    
    } else if (match.arg(method) == "mcmc") {
    
      fit <- model$sample(data = data,
                          output_dir = "Fits.nosync",
                          init = initfunc,
                          iter_warmup = 500, iter_sampling = 800,
                          adapt_delta = 0.99,
                          max_treedepth = 16,
                          chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
    
    } else if (match.arg(method) == "sim") {
    
      fit <- model$sample(data = data,
                          output_dir = NULL,
                          init = initfunc,
                          iter_sampling = 500,
                          fixed_param = TRUE)
  }
  return(fit)
}

# model <- cmdstan_model(modelpath)

## Model fit
fit <- drawSamples(model, data, method = "variational", initfunc = 0)

recoverysetup <- list(drawfile = basename(fit$output_files()),
                      metadata = fit$metadata(),
                      truepars = pars,
                      data = data,
                      model = calcModel)

fitbasename <- str_split(recoverysetup$drawfile[1], "-")[[1]]
fitbasename <- paste(fitbasename[1:(length(fitbasename)-2)], collapse = "-")
saveRDS(recoverysetup, file.path("Fits.nosync", glue("{fitbasename}.rds")))
  
## Other diagnostics
# fit$output()
# fit$time()
# fit$init()
# fit$draws(variables = NULL, inc_warmup = F)
# fit$return_codes()


## Extract draws -----------------------------------------------------------

drawpath <- fit$output_files() # dput(fit$output_files())
    ## alternatively
    # drawpath <- file.info(list.files("Fits", full.names = T)) %>%
    #   arrange(desc(mtime)) %>%
    #   slice(1:6) %>%
    #   rownames()


# drawnames <- unique(str_extract(fit$metadata()$stan_variables, "[a-zA-Z_]*"))
# targetvarname <- intersect(drawnames, names(truepars)) # this drops boilerplate variables like u
# targetvarname <- c(targetvarname, "state_init_log") # "y_hat_log" "sim_log", "y_hat_log"
excludevarname <- c("u") # state_init_log

### Get draws
## rstan format for use in other packages
stanfit <- rstan::read_stan_csv(drawpath)
draws <- rstan::extract(stanfit, pars = c(excludevarname), include = F)
    ## get draws with cmdstanr
    # draws <- read_cmdstan_csv(drawpath, variables = NULL)



## Summary -----------------------------------------------------------------

## Explore stanfit object
s <- summary(stanfit) # , pars = excludevarname, include = F
s

traceplot(stanfit, pars = c("r_log", "g_logit", "h_logit", "c_j_log", "alpha_obs"))
traceplot(stanfit, pars = c("b_log", "c_b_log", "Beta_l"))

varname <- fit$metadata()$stan_variables
pairs(stanfit, pars = varname[1:floor(0.2*length(varname))])
pairs(stanfit, pars = varname[(length(varname)-1):length(varname)])


## Shinystan
# library(shinystan)
# launch_shinystan(stanfit, rstudio = F)


