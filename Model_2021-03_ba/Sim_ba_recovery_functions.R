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
                          model,
                          Env,
                          envdependent = c(b = F, c_a = F, c_b = F, c_j = T, g = T, h = F, l = T, m_a = F, m_j = F, r = T, s = T)
                          ) {
  
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
  
  fit <- drawSamples(model, data, method = "mcmc", initfunc = 0,
                     dirpath = file.path(modeldir, "Sim_ba_recovery", "Fits.nosync"))
  
  ## draws will get saved unter fitbasename….csv
  fitbasename <- str_split(recoverysetup$drawfile[1], "-")[[1]]
  fitbasename <- paste(fitbasename[1:(length(fitbasename)-2)], collapse = "-")
  
  setup <- list(drawfile = basename(fit$output_files()),
                metadata = fit$metadata(),
                truepars = pars,
                times = times, # also in data of course
                data = data,
                pseudofixpointdata = pseudofixpointdata,
                envdependent = get(paste0("envdependent_", sub("-", "_", modelname)))
                )
  
  saveRDS(setup, file.path(modeldir, "Sim_ba_recovery", "Fits.nosync", fitbasename, "setup.rds")
  
  return(setup)
  
}


#### Returns +- the true start values ----------------------
getTrueInits <- function() {
  
  isragged <- grepl("^ba-rag", modelname) || modelname == "ba"
  
  state_init <- if (isragged) data$y0[which(!duplicated(data$rep_init2y0))] else data$y[,1,1,]
  
  truepars <- attr(data, "pars")
  truepars <- truepars[sapply(truepars, is.numeric)]
  parnames <- names(truepars)
  names(truepars) <- ifelse(str_ends(parnames, "_loc"), str_to_lower(parnames), parnames)
  
  newpars <- list(
    
    state_init = state_init,
    state_init_log = log(state_init),
    
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
  
  if(modelname %in% c("ba", "ba-rag", "ba-rag-ranef")) {
    R <- data.frame(y_hat_rep = draws$y_hat_rep[1,], y = data$y,
                    pops = data$pops[rep(data$rep_init2y0, each = data$n_reobs[1])])
    
    R %>% ggformula::gf_point(y_hat_rep ~ y) %>%
      gf_abline(gformula = NULL, slope = 1, intercept = 0)
    
    # plot(draws$state_init[1, data$rep_init2y0] ~ data$y0)
    
  }
  
}


plotEstimateVsTrue <- function(prior,
                               posterior,
                               true) {
  
  simpar <- attr(simdata, "pars")[[parname]]
  
  if(is.vector(simpar)) {
    
    Stanpar <- rstan::extract(stanfit, pars = parname)[[1]] %>%
      as.data.frame() %>%
      pivot_longer(cols = everything(), names_to = "species", values_to = "draw") %>%
      mutate(species = as.integer(as.factor(species))) %>%
      bind_cols(true = rep(simpar, length.out = nrow(.)))
    
    Stanpar %>%
      # sample_n(1000) %>%
      gf_point(draw ~ true, alpha = 0.3, size = 0.5, title = parname) %>%
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
  
  plot(density(prior), col = "lightblue")
  lines(density(posterior))
  abline(v = 1, col = "lightblue")
  
  
}




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


# plotParameterInEnv("Beta_g")

#### Draws vs true ----------------------
plotStatePairs <- function(standata) {
  
  ba_a_avg <- pi * ((standata$dbh_lower_a + standata$dbh_lower_b)/2/2)^2 * 1e-6
  
  D <- standata$Long %>%
    mutate(id = interaction(time, plot, loc)) %>%
    select(stage, species, abundance, id, time) %>%
    pivot_wider(id_cols = c("id", "time"), names_from = c("stage", "species"), values_from = c("abundance"))
  
  plotlevel <- as.integer(as.factor(D$time))
  
  D %>% 
    mutate(ba_sum = (a_1*ba_a_avg + b_1) + (a_2*ba_a_avg + b_2)) %>%
    select(-c("id", "time")) %>%
    pairs(col = plotlevel, pch = plotlevel)
  
}





