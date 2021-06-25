# Library -------------------------------------------------------------
library(here)
library(glue)

library(tidyverse)
library(ggformula)
library(magrittr)
library(tictoc)

library(rstan)
rstan_options(javascript = FALSE)
library(cmdstanr)
# install_cmdstan(cores = 3)
library(bayesplot)


# Orientation -------------------------------------------------------------
setwd(here())
modeldir <- dir(pattern = glue("^(Model).*ba$"))

# Source -------------------------------------------------------------
source(file.path(modeldir, "Sim_ba_functions.R"))
source(file.path(modeldir, "Sim_ba_recovery_functions.R"))
source(file.path(modeldir, "Fit_ba_functions.R"))
source(file.path(modeldir, "Fit_ba_test_functions.R"))


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
times1 <- 2:20

Sim1 <- simulateOneSeries(initialstate1, times = times1, pars = pars1,
                          processerror = F, obserror = F)

matplot(Sim1[,-1], type = "b", ylab="N",
        pch = rep(c("J", "A", "B"), each = pars1$n_species),
        col = 1:pars1$n_species, xlab = "time") # log='y'



## Multiple time series in env ---------------------------------------------------------------
envdep_demo <- c(b = F, c_a = F, c_b = F, c_j = T, g = T, h = F, l = T, m_a = F, m_j = F, r = T, s = T)
pars_demo <- generateParameters(n_species = 2, n_locs = 50)
E_demo <- simulateEnv(pars_demo$n_env, pars_demo$n_locs)
P_env_demo <- transformParameters(pars_demo, E_demo, envdep_demo, ranef = F, returndf = T)
S_demo <- simulateMultipleSeriesInEnv(pars_demo, E_demo, times = 2:20,
                                      envdependent = envdep_demo, ranef = F, processerror = F, obserror = F, independentstart = F)

#### Plot time series
S_demo %>%
  # filter(stage %in% c("a")) %>%
  ggplot(mapping = aes(x = time, y = abundance, color = env_1, lty = stage, group = interaction(loc, stage, plot))) +
  geom_line() +
  facet_wrap(facets = c("species"))

#### Plot parameters
P_env_demo %>%
  filter(parameter %in% c("c_j_loc")) %>%
  gf_point(q ~ env_1 | species) %>%
  gf_smooth()

P_env_demo %>%
  filter(parameter %in% c("g_loc")) %>%
  gf_point(q ~ env_1 | species) %>%
  gf_smooth()

P_env_demo %>%
  filter(parameter %in% c("r_loc", "l_loc")) %>%
  gf_point(q ~ env_1 | parameter*species) %>%
  gf_smooth()

ggplot(filter(P_env_demo, parameter == "s_loc"),
       mapping = aes(x = env_1, y = q, color = species, group = species)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(facets = c("species"))



# ————————————————————————————————————————————————————————————————————————————————— #
# Fit simulated data with stan    --------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## Orient and compile model. -------------------------------------------------------

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



## Simulate stan model data --------------------------------------------------------
parseed <- 1

pars <- generateParameters(seed = parseed, n_locs = 100, n_species = 2, n_plotsperloc = 4, knockoutparname = c("m_a", "m_j"))
Env <- simulateEnv(n_env = pars$n_env, n_locs = pars$n_locs)


envdependent_ba_rect <- c(b = F, c_a = F, c_b = F, c_j = T, g = T, h = F, l = T, m_a = F, m_j = F, r = T, s = T)
envdependent_ba_rag <- c(b = F, c_a = F, c_b = F, c_j = F, g = T, h = F, l = T, m_a = F, m_j = F, r = T, s = T)
envdependent_ba <- c(b = F, c_a = F, c_b = F, c_j = F, g = T, h = F, l = T, m_a = F, m_j = F, r = T, s = T)
envdependent_ba_rag_ranef <- envdependent_ba_rag

times <- c(2, 3, 4)
data <- simulateMultipleSeriesInEnv(pars, Env, times = times,
                                    envdependent = get(paste0("envdependent_", sub("-", "_", modelname))),
                                    logstate = F,
                                    modelstructure = modelname, # !!! this determines the data layout
                                    format = "stan", priorfactor = 10,
                                    obserror = T, processerror = F, independentstart = F)

## just a data set with a long time span to check fix point (equilibrium) recovery
pseudofixpointdata <- simulateMultipleSeriesInEnv(pars, Env, times = c(1, 500),
                                                   envdependent = get(paste0("envdependent_", sub("-", "_", modelname))),
                                                   logstate = F,
                                                   modelstructure = modelname, # !!! this determines the data layout
                                                   format = "stan",
                                                   obserror = F, processerror = F, independentstart = F)


Data_long <- data$Long

####### Plot time series
Data_long %>%
  # filter(stage %in% c("a")) %>%
  ggplot(mapping = aes(x = time, y = abundance, color = env_1, lty = stage, group = interaction(loc, stage, plot))) +
  geom_line() +
  facet_wrap(facets = c("species"))




## Draw from model --------------------------------------------------------------
## Model fit
fit <- drawSamples(model, data, method = "mcmc", initfunc = 0)


## Other diagnostics
# fit$output()
# fit$time()
# fit$init()
# fit$draws(variables = NULL, inc_warmup = F)
# fit$return_codes()


## Extract draws and save setup object -----------------------------------------------------------

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


fitbasename <- str_split(drawpath[1], "-")[[1]]
fitbasename <- paste(fitbasename[1:(length(fitbasename)-2)], collapse = "-")
fitbasename <- basename(fitbasename)

setup <- list(drawfile = fitbasename,
              metadata = fit$metadata(),
              truepars = pars,
              times = times, # also in data of course
              data = data,
              pseudofixpointdata = pseudofixpointdata,
              envdependent = get(paste0("envdependent_", sub("-", "_", modelname)))
              )

saveRDS(setup, file.path(modeldir, "Sim_ba_recovery", "Fits.nosync", glue("{fitbasename}.rds")))


## Summary -----------------------------------------------------------------

## Explore stanfit object
s <- summary(stanfit) # , pars = excludevarname, include = F
s

# traceplot(stanfit, pars = c("r_log", "g_logit", "h_logit", "c_j_log", "alpha_obs"))
traceplot(stanfit, pars = c("b_log", "c_b_log", "Beta_g"))

varname <- fit$metadata()$stan_variables
pairs(stanfit, pars = varname[1:floor(0.2*length(varname))])
pairs(stanfit, pars = varname[(length(varname)-1):length(varname)])


## Shinystan
# library(shinystan)
# launch_shinystan(stanfit, rstudio = F)


## Quick recovery check -----------------------------------------------------------------

## Recovery check expected format:
## lists (data frames of compatible-length variables)

#### 0. Predicted vs. true
plotPredictedVsTrue(draws, setup$data)

### 1. Estimates vs. priors
n_draws <- length(draws$lp__)
priors <- drawPriors(data, n = n_draws)

plotEstimates(priors = priors$prior_b_log,
              posteriors = unpackDrawsArray(draws$b_log),
              expectedvalues = as.list(setup$truepars$b_log))

plotEstimates(priors = priors$prior_Vertex_c_j,
              posteriors = unpackDrawsArray(draws$vertex_c_j),
              expectedvalues = as.list(c(setup$truepars$Beta_c_j)))


### 1a. Estimates line plots vs. true

## here are the drawn parameters: # setup$metadata$stan_variables
truepars <- attr(setup$data, "pars") # includes transformed

plotEstimateLine(priors = priors$prior_h_logit,
                 posteriors = unpackDrawsArray(draws$h_logit),
                 expectedvalues = setup$truepars$h_logit)

plotEstimateLine(priors = cbind(a = priors$prior_c_b_log, b = priors$prior_b_log),
                 posteriors = cbind(a = unpackDrawsArray(draws$c_b_log), b = unpackDrawsArray(draws$b_log)),
                 expectedvalues = c(setup$truepars$c_b_log, setup$truepars$b_log))

plotEstimateLine(priors = unpackDrawsArray(draws$G_logit),
                 posteriors = unpackDrawsArray(draws$G_logit),
                 expectedvalues = c(truepars$G_log))


## 2. Posterior pairs
# varnames <- setup$metadata$stan_variables
# corpars_a <- varnames[c(1, 2, 3, 6, 7)]
# corpars_b <- varnames[c(2, 4, 5, 7, 8)]

# ?rstan:::pairs.stanfit
pairs(stanfit, pars = c("h_logit", "c_b_log"))

## 3. pairs for states in data to inspect how states are correlated
plotStatePairs(data)

## 4. Beta in Env
plotBetaInEnv(Betadraws = draws$Beta_g, Beta_true = truepars$Beta_g)
plotBetaInEnv(Betadraws = draws$Beta_g, Beta_true = truepars$Beta_g, env = 2)
plotBetaInEnv(Betadraws = draws$Beta_c_j, Beta_true = truepars$Beta_c_j)
plotBetaInEnv(Betadraws = draws$Beta_s, Beta_true = truepars$Beta_s, env = 2) # Check priors with marginal plot!

## 5. Calculate RMSE for all parameters
rmse <- calculateRMSE(standraws = draws, standata = data)
boxplotRMSE(rmse) # boxplot(rmse)



## Quick test power checks -------------------------------------------------


## 0. inspect iterations
hist(draws$iterations_fix)
table(draws$iterations_fix)


## 1. Compare fix point
pseudofixpoint <- setup$pseudofixpointdata$y[, 2, 1,]
drawfixpoint <- apply(draws$state_fix[,, 1:6], c(2, 3), mean)
plot(drawfixpoint ~ pseudofixpoint, col = rowMeans(draws$iterations_fix))

## 3. Contribution of parameter given dominance for env bins


