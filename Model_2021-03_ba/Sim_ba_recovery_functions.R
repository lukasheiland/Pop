# Library -----------------------------------------------------------------
# library(GGally) ##?
library(cmdstanr)


# ——————————————————————————————————————————————————————————————————————————————————— #
# Fitting functions  ----------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

setupRecovery <- function(pars,
                          times,
                          independentstart,
                          priorfactor,
                          recoveryname = "",
                          model,
                          Env,
                          envdependent = c(b = F, c_a = F, c_b = F, c_j = T, g = T, h = F, l = T, m_a = F, m_j = F, r = T, s = T),
                          ...) # passed on to draw samples and cmdstanr methods
{
  
  modelname <- "ba-rect"
  modeldir <- dir(pattern = glue("^(Model).*ba$"))
  
  data <- simulateMultipleSeriesInEnv(pars = pars, Env = Env, times = times,
                                      envdependent = envdependent,
                                      logstate = F,
                                      modelstructure = modelname, # !!! this determines the data layout
                                      format = "stan", priorfactor = priorfactor,
                                      obserror = T, processerror = F, independentstart = independentstart)
  
  pseudofixpointdata <- simulateMultipleSeriesInEnv(pars = pars, Env = Env, times = c(1, 500),
                                                    envdependent = envdependent,
                                                    logstate = F,
                                                    modelstructure = modelname, # !!! this determines the data layout
                                                    format = "stan",
                                                    obserror = F, processerror = F, independentstart = F)
  
  
  require(cmdstanr)
  
  ## try multiple times for when initial values lead to errors  
  fit <- NULL
  attempt <- 1
  while( is.null(fit) && attempt <= 5 ) {
    attempt <- attempt + 1
    try({
      fit <- drawSamples(model, data, method = "mcmc", initfunc = 0,
                         iter_warmup = 400, iter_sampling = 600,
                         dirpath = file.path(modeldir, "Sim_ba_recovery", "Fits.nosync"), ...)
    })
  }
  
  setup <- list(drawfile = basename(fit$output_files()),
                recoveryname = recoveryname,
                metadata = fit$metadata(),
                truepars = pars,
                times = times, # also in data of course
                data = data,
                pseudofixpointdata = pseudofixpointdata,
                envdependent = envdependent)
  
  ## draws will get saved under fitbasename….csv
  fitbasename <- basename(fit$output_files()[1])
  fitbasename_generic <- str_replace(fitbasename, pattern = "-\\d-", "-x-")
  id_generic <- tools::file_path_sans_ext(fitbasename_generic)
  
  saveRDS(setup, file.path(modeldir, "Sim_ba_recovery", "Fits.nosync", paste(id_generic, recoveryname, "setup.rds", sep = "_")))
  
  return(setup)
  
}


#### Returns +- the true start values ----------------------
getTrueInits <- function() {
  
  isragged <- grepl("^ba-rag", modelname) # || modelname == "ba"
  
  state_init <- if (isragged) {
      data$y0[which(!duplicated(data$rep_init2y0))]
    } else if(modelname %in% c("ba", "ba_test")) {
      data$state_init
    } else { data$y[,1,1,] }
  
  truepars <- attr(data, "pars")
  truepars <- truepars[sapply(truepars, function(x) { is.numeric(x) && !anyNA(x) })]
  parnames <- names(truepars)
  names(truepars) <- ifelse(str_ends(parnames, "_loc"), str_to_lower(parnames), parnames)
  
  newpars <- list(
    
    state_init = state_init,
    state_init_log = log(state_init),
    sigma_l = rep(1e-16, truepars$n_species),
    
    u = replicate(truepars$n_locs, matrix(rnorm(truepars$n_species*3, 0, 0.001), nrow = truepars$n_species*3, ncol = data$timespan_max))
  )
  
  inits <- c(newpars, truepars)
  truepars <<- truepars
  
  return(inits)
}

#### Returns viable start values ---------------------------
getInits <- function() {
  
  responsescaleerror <- 0.1
  
  isragged <- grepl("^ba-rag", modelname)
  
  truepars <<- attr(data, "pars")
  n_species <- truepars$n_species
  n_locs <- truepars$n_locs
  
  state_init <- if (isragged) data$y0[which(!duplicated(data$rep_init2y0))] + rnorm(data$N_init, 0, responsescaleerror) else data$y[,1,1,] + rnorm(data$y[,1,1,], 0, responsescaleerror)
  
  b_log <- rnorm(n_species, -1, 0.01)
  c_a_log <- rnorm(n_species, -3.3, 0.2)
  c_b_log <- rnorm(n_species, -3, 0.1)
  c_j_log <- rnorm(n_species, -5, 0.5)
  g_logit <- rnorm(n_species, -1, 0.2)
  h_logit <- rnorm(n_species, -0.5, 0.3)
  m_a_log <- rnorm(n_species, -1.5, 0.5)
  m_j_log <- rnorm(n_species, -1.2, 0.1)
  r_log <- rnorm(n_species, 1.5, 0.2)
  l_log <- rnorm(n_species, -1, 0.2)
  s_log <- rnorm(n_species, -2.9, 0.1)
  
  beta_null <- function() { matrix(c(-1, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.001) }
  
  
  inits <- list(
    
    b_log = b_log,
    c_a_log = c_a_log,
    c_b_log = c_b_log,
    c_j_log = c_j_log,
    
    g_logit = g_logit,
    h_logit = h_logit,
    
    l_log = l_log,
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
    Beta_l = matrix(c(1, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_m_a = matrix(c(-4, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_m_j = matrix(c(-2, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_r = matrix(c(3, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    Beta_s = matrix(c(-3, rep(0, truepars$n_beta-1)), ncol = truepars$n_species, nrow = truepars$n_beta) + rnorm(truepars$n_beta*truepars$n_species, 0, 0.3),
    
    b = exp(b_log),
    c_a = exp(c_a_log),
    c_b = exp(c_b_log),
    c_j = exp(c_j_log),
    
    g = plogis(g_logit),
    h = plogis(h_logit),
    
    l = exp(l_log),
    m_a = exp(m_a_log),
    m_j = exp(m_j_log),
    r = exp(r_log),
    s = exp(s_log),
    
    b_loc = matrix(rep(exp(b_log), n_locs), nrow = n_locs, byrow = T),
    c_a_loc =  matrix(rep(exp(c_a_log),n_locs), nrow = n_locs, byrow = T),
    c_b_loc =  matrix(rep(exp(c_b_log),n_locs), nrow = n_locs, byrow = T),
    c_j_loc =  matrix(rep(exp(c_j_log),n_locs), nrow = n_locs, byrow = T),
    g_loc =  matrix(rep(plogis(g_logit),n_locs), nrow = n_locs, byrow = T),
    h_loc =  matrix(rep(plogis(h_logit),n_locs), nrow = n_locs, byrow = T),
    l_loc =  matrix(rep(exp(l_log),n_locs), nrow = n_locs, byrow = T),
    m_a_loc =  matrix(rep(exp(m_a_log),n_locs), nrow = n_locs, byrow = T),
    m_j_loc =  matrix(rep(exp(m_j_log),n_locs), nrow = n_locs, byrow = T),
    r_loc =  matrix(rep(exp(r_log),n_locs), nrow = n_locs, byrow = T),
    s_loc =  matrix(rep(exp(s_log),n_locs), nrow = n_locs, byrow = T),
    
    shape_par      = c(10, 10, 10),
    sigma_process  = c(0.01),
    
    sigma_obs      = c(1.1, 1.1),
    alpha_obs      = c(10, 20),
    alpha_obs_inv   = c(0.1, 0.2),
    phi_obs      = c(10, 20),
    
    state_init = state_init,
    state_init_log = log(state_init),
    
    u = replicate(truepars$n_locs, matrix(rnorm(truepars$n_species*3, 0, 0.001), nrow = pars$n_species*3, ncol = data$timespan_max))
  )
  
  return(inits)
}


# ——————————————————————————————————————————————————————————————————————————————————— #
# Plotting functions  ----------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

# Plot predicted (generated quantities) vs. true values -----------------------------------------
plotPredictedVsTrue <- function(draws, data, modelname = "ba-rect") {
  
  truepars <- attr(data, "pars")
  
  if(modelname == "ba-rect") {
    
    par(mfrow = c(2, 2))
        
    # predictions per life stages
    plot(c(draws$y_hat_rep[100,,3,,data$i_j]) ~ c(data$y[,3,,data$i_j]))
    abline(0, 1)
    plot(c(draws$y_hat_rep[100,,3,,data$i_a]) ~ c(data$y[,3,,data$i_a]))
    abline(0, 1)
    plot(c(draws$y_hat_rep[1,,3,,data$i_b]) ~ c(data$y[,3,,data$i_b]))
    abline(0, 1)
    
    # plot(c(draws$y_sim[1,,3,,data$i_j]) ~ c(data$y[,3,,data$i_j]))
    # abline(0, 1)
    # plot(c(draws$y_sim[1,,3,,data$i_a]) ~ c(data$y[,3,,data$i_a]))
    # abline(0, 1)
    # plot(c(draws$y_sim[1,,3,,data$i_b]) ~ c(data$y[,3,,data$i_b]))
    # abline(0, 1)
        
    par(mfrow = c(1, 1))
        
  }
  
  if(modelname %in% c("ba", "ba-rag", "ba-rag-ranef", "ba_test")) {
    R <- data.frame(y_hat_rep = draws$y_hat_rep[1,], y = data$Long$abundance, # y_sim = draws$y_sim[1,],
                    pops = data$pops[rep(data$rep_init2y0, each = data$n_obs[1])])
    
    R %>% ggformula::gf_point(y_hat_rep ~ y) %>%
      gf_abline(gformula = NULL, slope = 1, intercept = 0)
    
    # plot(draws$state_init[1, data$rep_init2y0] ~ data$y0)
    
  }
  
}


## Arrange estimates on a line with expected values ---------------------------------------
plotEstimateLine <- function(priors,
                             posteriors,
                             expectedvalues) {
  
  Prior <- pivot_longer(as.data.frame(priors), cols = everything(), names_to = "i", values_to = "prior") %>%
    mutate(i = as.integer(fct_inorder(i)))
  Posterior <-  pivot_longer(as.data.frame(posteriors), cols = everything(), names_to = "i", values_to = "posterior") %>%
    select(-i)
  Data <- bind_cols(Prior, Posterior, expected = unlist(expectedvalues)[Prior$i]) %>%
    pivot_longer(cols = c("prior", "posterior"), names_to = "when")
  
  Data %>%
    sample_n(10000) %>%
    gf_point(value ~ expected, alpha = 0.2, size = 0.5, color = ~ when) %>%
    gf_abline(slope = 1, intercept = 0, gformula = NULL, alpha = 0.1)
}



#### Parameter vs Env ----------------------
predictBeta <- function(Beta, x = seq(-1, 1, by = 0.01), env = 1, species = 1) {
  E <- as.data.frame(x)
  polyformula <- as.formula(paste("~", paste("poly(", colnames(E), ", 2)", collapse = "+")))
  X <- model.matrix(polyformula, data = E)
  X %*% Beta[c(1, 2:3 + {if(env == 2) env else 0} ),species]
}


plotBetaInEnv <- function(Betadraws, Beta_true = NULL, x = seq(-1, 1, by = 0.01), env = 1, species = 1) {
  
  n_draws <- dim(Betadraws)[1]
  
  plot(x = x, y = predictBeta(apply(Betadraws, c(2, 3), mean), x = x, env = env, species = species), type = "n")
  
  for(i in 1:n_draws) {
    lines(x = x, y = predictBeta(Betadraws[i, ,], x = x, env = env, species = species), col = alpha("red", 0.3), lwd = 0.1)
  }
  
  if (!is.null(Beta_true)) lines(x, predictBeta(Beta_true, x = x, env = env, species = species), col = "blue", lwd = 1)
  
}


# plotParameterInEnv("Beta_g")

#### Draws vs true ----------------------
plotStatePairs <- function(standata) {
  
  ba_a_avg <- pi * ((standata$dbh_lower_a + standata$dbh_lower_b)/2/2)^2 * 1e-6
  
  D <- standata$Long %>%
    mutate(id = interaction(time, plot, loc)) %>%
    dplyr::select(stage, species, abundance, id, time) %>%
    pivot_wider(id_cols = c("id", "time"), names_from = c("stage", "species"), values_from = c("abundance"))
  
  plotlevel <- as.integer(as.factor(D$time))
  
  D %>% 
    mutate(ba_sum = (a_1*ba_a_avg + b_1) + (a_2*ba_a_avg + b_2)) %>%
    dplyr::select(-c("id", "time")) %>%
    pairs(col = plotlevel, pch = plotlevel)
  
}


# ——————————————————————————————————————————————————————————————————————————————————— #
# Recovery summary statistics  -------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

calculateRMSE <- function(standraws, standata) {
  
  truepars <- attr(standata, "pars")
  true <- c(truepars, list(y_hat_rep = data$y))
  
  calc <- function(d, p) {
    sparedims <- 2:length(dim(d))
    n_draws <- dim(d)[1]
    
    mse <- apply((d - repArray(p, n_draws))^2, sparedims, mean)
    return(sqrt(mse))
  }
  
  names <- intersect(names(standraws), names(true))
  
  rmse <- mapply(calc, standraws[names], true[names], SIMPLIFY = F, USE.NAMES = T)
  return(rmse)
}


boxplotRMSE <- function(rmse, n_species = 2) {
  rmse$y_hat_rep <- matrix(c(rmse$y_hat_rep), ncol = n_species) # y_hat structure does not play a role, recast into a 2xX matrix like the others
  R <- do.call(rbind, rmse)
  
  names_rmse <- rep(names(rmse), sapply(rmse, function(x) length(x)/n_species))
  R <- data.frame(R, par = names_rmse) %>%
    pivot_longer(cols = starts_with("X"), names_to = "species")
  
  boxplot(value ~ par, data = R)
}

