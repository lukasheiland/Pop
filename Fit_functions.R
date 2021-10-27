# ——————————————————————————————————————————————————————————————————————————————————#
# Functions for fitting  ----------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————#

## formatStanData --------------------------------
# Stages  <- tar_read("Stages_scaled")
# Stages_transitions  <- tar_read("Stages_transitions")
# taxon_s  <- tar_read("taxon_s")
# threshold_dbh <- tar_read("threshold_dbh")

formatStanData <- function(Stages, Stages_transitions, taxon_s, threshold_dbh) {
  
  #### Essential worker
  ## reps each element of the first vector by the corresponding number in the second vector (c and unlist make sure that a vector is returned)
  vrep <- function(...) unlist( c( Vectorize(rep.int)(...) ) ) # :)))))
  
  ## Prepare ba_a_upper, the upper basal area of any one tree in A.
  radius_a_upper <- threshold_dbh/2
  ba_a_upper <-  pi * radius_a_upper^2 * 1e-6
  
  #### Prepare ba_a_avg
  ## vector[N_species] for ba_a_avg
  ## For multiplication with count A, ba_a_avg is the average basal area of an individual tree in class A and has unit $\mathit{areaunit} \cdot 1^-1$, so that A * ba_a_avg has unit $\mathit{areaunit} \cdot ha^-1$.
  ba_a_avg <- Stages %>%
    filter(stage == "A") %>%
    group_by(tax) %>%
    summarize(ba_a_avg = mean(ba_ha/count_ha, na.rm = T)) %>% ## [area/ha * ha/1 == ba/1]
    pull(ba_a_avg, name = tax)
  
  
  #### Prepare special long data set structures with different amounts of long casting along different ids
  ## The most comprehensive data set is L_y with grouping locations/resurveys/pops/plots.
  ## NOTE THAT plots HAS TO BE THE LAST GROUP IN SORTING
  ## Everything is subset from the master subsets with these groupings (*_reobs, *_y0) and thus consistently sorted.
  ## Factor "pops" is structured stages/species.
  
  Stages_transitions %<>%
    ungroup() %>%
    dplyr::select("plotid", "tax", "obsid",
                  "timediff_plot",
                  "count_A2B_plot", "count_J2A_plot",
                  "count_A_sum_before", "count_J_integr_plot", "count_J_sum_before", "count_A_integr_plot",
                  "h_plot", "g_plot") %>%
    # mutate_at(c("count_A2B_plot", "count_J2A_plot", "count_J_integr_plot", "count_A_integr_plot"), round) %>%
    # mutate(joinid = interaction(plotid, tax)) %>%
    # dplyr::select(-plotid, -tax) %>%
    
    group_by(plotid, obsid, tax) %>%
    slice(1)
  
  G <- Stages_transitions %>%
    # filter(!is.na(g_plot))
    filter(!is.na(count_J2A_plot) & isTRUE(count_J_integr_plot > 0)) ## also drops NAs

  H <- Stages_transitions %>%
    # filter(!is.na(h_plot))
    filter(!is.na(count_A2B_plot) & isTRUE(count_A_integr_plot > 0)) ## also drops NAs
  
  Stages %<>%
    
    ## Join Stages_transitions
    # bind_cols(Stages_transitions[match(interaction(.$plotid, .$tax), Stages_transitions$joinid), ]) %>%
    
    ## Synonyms for consistency wiht model lingo
    group_by(clusterid) %>%
    mutate(plot = match(plotid, unique(plotid))) %>% ## for numbers designating which corner it is, extract the substring at the end of plotid
    ungroup() %>%
    mutate(loc = as.integer(as.factor(clusterid))) %>%
    
    ## factor ordering!!!
    mutate(stage = factor(stage, levels = c("J", "A", "B", "BA")),
           tax = factor(tax, levels = levels(taxon_s))) %>%
    
    ## Stages are measured in different terms: ba or count
    mutate(y = case_when(
      stage == "J" ~ count_ha_r,
      stage == "A" ~ count_ha_r,
      stage == "B" ~ ba_ha_r,
      stage == "BA" ~ ba_ha_r
    )) %>%
    
    ## Different levels of observation error were assumed for J (area count sampling), the stage A (counts from fixed angle sampling), and stage B (basal area from fixed angle sampling).
    mutate(obsmethod = fct_recode(stage, "j" = "J", "a" = "A", "ba" = "B", "ba" = "BA")) %>% 
    
    ## pop is just an id for the initial states vector
    mutate(pop = interaction(tax, stage)) %>%
    
    arrange(loc, time, pop, stage, tax, plot) %>%
    
    ## scale times within group to start with 1
    group_by(loc) %>%
    mutate(t = as.integer(lubridate::year(time))) %>%
    mutate(t = t - min(t) + 1) %>%
    mutate(isy0 = (t == 1)) %>%
    ungroup()
  
  
  ## Format: [L_y] — locations/obs/pops(/stage/species)/plots
  ## for fitting just use stages J, A, B
  S <- Stages %>% 
    filter(stage %in% c("J", "A", "B")) %>%
    mutate_at(c("stage", "pop"), droplevels) %>%
    arrange(loc, time, pop, stage, tax, plot) ## safety first
  
  
  ## Format: [L_init] — locations/pops
  S_init <- filter(S, isy0) %>%
    st_drop_geometry() %>%
    dplyr::select(-isy0) %>%
    group_by(loc, pop, stage, tax) %>%
    summarize(n_plots = n_distinct(plot), .groups = "drop")
  
  S <- S %>%
    st_drop_geometry()
  
  ## Format: [L_yhat] — locations/obs/pops
  S_yhat <- S %>%
    group_by(loc, time, pop, stage, tax) %>%
    summarize(n_plots = n_distinct(plot), .groups = "drop")
  
  ## Format: [L_times] — locations/resurveys
  S_times <- S %>%
    group_by(loc, obsid) %>% # first, group by obsid to extract one time per obsid within loc (in case there are multiple days somehow, which is not the case atm!)
    summarize(t = first(t)) %>%
    group_by(loc) %>%
    summarize(t = unique(t), .groups = "drop") ## This way its certain, that there is one time per obsid
  
  ## Format: [N_locs] — locations
  S_locs <- S %>%
    group_by(loc) %>%
    ## assumes completion within locations!
    summarize(time_max = max(t),
              timespan = diff(range(t)) + 1,
              n_tax = n_distinct(tax),
              n_plots = n_distinct(plot),
              n_pops = n_distinct(pop),
              n_obs = n_distinct(t),
              n_yhat = n_distinct(interaction(pop, t)),
              alt_loc_s = first(alt_loc_s),
              phCaCl_esdacc_s = first(phCaCl_esdacc_s),
              waterLevel_loc_s = first(waterLevel_loc_s),
              clusterid = first(clusterid),
              .groups = "drop")
  
  
  ## Format: [L_init] — locations/pops
  # S_a2b <- S %>%
  #   # st_drop_geometry() %>%
  #   ## !!! sic! Use A for matching in S_yhat later.
  #   filter(stage == "A") %>% 
  #   
  #   ## !!! sic! Assign to preceding state for matching in S_yhat later.
  #   mutate(a2b = case_when(obsid == "DE_BWI_1987" ~ count_A2B_2002_plot,
  #                          obsid == "DE_BWI_2002" ~ count_A2B_2012_plot)) %>%
  #   mutate(timediff = case_when(obsid == "DE_BWI_1987" ~ timediff_2002,
  #                               obsid == "DE_BWI_2002" ~ timediff_2012)) %>%
  #   drop_na(a2b, timediff) %>%
  #   arrange(loc, time, pop, stage, tax, plot) %>% ## safety first
  #   mutate(rep_yhat2a2b = match(interaction(loc, time, pop), with(S_yhat, interaction(loc, time, pop))),
  #          rep_species2a2b = as.integer(factor(tax)))
  
  
  #### Prepare design matrix
  Env <- S_locs[c("phCaCl_esdacc_s", "waterLevel_loc_s")]
  if (anyNA(Env)) {
    X <- matrix(0, nrow = nrow(Env), ncol = ncol(Env))
    message("With Env having NAs, an empty design matrix was produced.")
  } else {
    polyformula <- as.formula(paste("~", paste("poly(", colnames(as.data.frame(Env)), ", 2)", collapse = "+")))
    X <- model.matrix(polyformula, data = as.data.frame(Env))
  }
  
  
  #### Prepare ldd smooth. Predicted with log-link 
  L_smooth_log <- S %>%
    group_by(loc) %>%
    summarize_at(paste("s", taxon_s, sep = "_"), first) %>%
    column_to_rownames(var = "loc") %>%
    as.matrix()
  
  #### Some data for multiple reuse in the list
  N_species <-  length(unique(S$tax))
  L_yhat <- nrow(S_yhat)
  
  
  #### Prepare priors
  prior_state_init_log <- filter(S, isy0) %>%
    filter(stage %in% c("J", "A", "B")) %>%
    group_by(pop) %>%
    summarize(y = log(mean(y, na.rm = T))) %>% ## log(mean()) is imperfect and itroduces a small bias, but median won't work because of the many zeroes
    pull(y)
  
  
  #### The stan-formatted list
  data <- list(
    
    L_times = nrow(S_times),
    L_yhat = L_yhat,
    L_init = nrow(S_init),
    L_y = nrow(S),
    # L_a2b = nrow(S_a2b),
    
    N_locs = nrow(S_locs),
    N_species = N_species,
    N_pops = length(unique(S$pop)),
    N_beta = ncol(X),
    
    n_obs = S_locs$n_obs,
    n_yhat = S_locs$n_yhat,
    
    i_j = 1:N_species,
    i_a = (1:N_species) + N_species,
    i_b = (1:N_species) + 2*N_species,
    
    rep_yhat2y = vrep(1:L_yhat, S_yhat$n_plots), ## repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
    rep_obsmethod2y = as.integer(S$obsmethod),
    # rep_yhat2a2b = S_a2b$rep_yhat2a2b,
    # rep_species2a2b = S_a2b$rep_species2a2b,
    rep_pops2init =  as.integer(S_init$pop),
    rep_pops2y =  as.integer(S$pop),
    
    time_max = S_locs$time_max,
    times = S_times$t,
    # timediff = S_a2b$timediff,
    
    X = X,
    L_smooth_log = L_smooth_log, ## array[N_locs] vector<lower=0>[N_species] L_smooth_log;
    L_smooth = exp(L_smooth_log),
    
    y = S$y,
    # a2b = S_a2b$a2b,
    
    ## Settings corner
    tolerance_fix = 0.001,
    ba_a_upper = ba_a_upper,
    ba_a_avg = ba_a_avg,
    
    #### priors
    prior_state_init_log = prior_state_init_log,
    
    ## the transitions.
    ## no values between 0 and 1, so rounding is fine (H$count_A_integr_plot > 0 & H$count_A_integr_plot < 1) %>% any()
    L_g = nrow(G),
    L_h = nrow(H),
    y_j2a = round(G$count_J2A_plot),
    y_a2b = round(H$count_A2B_plot),
    y_j = round(G$count_J_integr_plot),
    y_a = round(H$count_A_integr_plot),
    # y_g = G$g_plot,
    # y_h = H$h_plot,
    species_g = as.integer(factor(G$tax, levels = levels(taxon_s))),
    species_h = as.integer(factor(H$tax, levels = levels(taxon_s)))
  )
  
  return(data)
}


## fitTransition --------------------------------
# data_stan  <- tar_read("data_stan")
# which  <- "g"
# model_transitions  <- tar_read("model_transitions")
# fitpath <- "Fits.nosync/"

fitTransition <- function(data_stan, which, model_transitions, fitpath = "Fits.nosync/") { # priors!
  
  isg <- (which == "g")
  
  d <- with(data_stan, list(
    ## the transitions
    L = if(isg) L_g else L_h,
    y_trans = if(isg) y_j2a else y_a2b,
    y_base = if(isg) y_j else y_a,
    rep_species = if(isg) species_g else species_h,
    N_species = N_species
  ))
    
  n_chains <- 3
  fit_transition <- model_transitions$sample(data = d,
                                             output_dir = fitpath,
                                             # iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                                             chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
  
  
  # bayesplot::mcmc_trace(fit_transition$draws())
  # bayesplot::mcmc_pairs(fit_transition$draws())
  # bayesplot::mcmc_areas(fit_transition$draws(variables = c("prop_logit")), area_method = "scaled height")
  
  message("Summary of the the fit for parameter ", which, ":")
  print(fit_transition$summary())
  
  return(fit_transition)
}


## formatPriors --------------------------------
# data_stan  <- tar_read("data_stan")
# fit_g  <- tar_read("fit_g")
# fit_h  <- tar_read("fit_h")
# fits_seedlings  <- tar_read("fits_seedlings")

formatPriors <- function(data_stan, weakpriors, fit_g, fit_h, fits_seedlings, widthfactor = 1) {
  
  ## Matrices[draws, species]
  Draws_g <- fit_g$draws(variables = "prop_logit", format = "draws_matrix") %>% as.data.frame()
  Draws_h <- fit_h$draws(variables = "prop_logit", format = "draws_matrix") %>% as.data.frame()
  Draws_seedlings <- posterior_samples(fits_seedlings, fixed = F, pars = c("Intercept$", "b_ba_ha$"))
  
  pars_g <- lapply(Draws_g, function(d) MASS::fitdistr(d, "normal")$estimate)
  pars_h <- lapply(Draws_h, function(d) MASS::fitdistr(d, "normal")$estimate)
  pars_seedlings <- lapply(Draws_seedlings, function(d) MASS::fitdistr(d, "normal")$estimate)
  
  if(widthfactor != 1) {
    pars_g <- lapply(pars_g, function(p) c(p["mean"], widthfactor*p["sd"]))
    pars_h <- lapply(pars_h, function(p) c(p["mean"], widthfactor*p["sd"]))
    pars_seedlings <- lapply(pars_seedlings, function(p) c(p["mean"], widthfactor*p["sd"]))
  }
  
  ## The model assumes array[2] vector[N_species] prior_*; which means that the vectors stretch over rows!
  
  priors <- list(
    prior_g_logit = bind_cols(pars_g), ## Matrix[N_species, (mu, sigma)]
    prior_h_logit = bind_cols(pars_h),
    prior_l_log = bind_cols(dplyr::select(as.data.frame(pars_seedlings), ends_with("Intercept"))),
    prior_r_log = bind_cols(dplyr::select(as.data.frame(pars_seedlings), ends_with("b_ba_ha")))
  )
  
  return(c(data_stan, weakpriors, priors))
}



#### Returns viable start values ---------------------------
# getInits <- function() {
#   
#   responsescaleerror <- 0.1
#   data <- tar_read(data_stan) ## This is fine in this case, as getInits is only called inside targets that would be invalidated after change of data_stan
#   n_species <- data$N_species
#   n_beta <- data$N_beta
#   n_locs <- data$N_locs
#   
#   state_init <- matrix(rlnorm(data$L_yhat, log(0.01), 0.01), nrow = n_locs, ncol = data$N_pops) # array[N_locs] vector<lower=0>[N_pops] state_init;
#   
#   b_log <- rnorm(n_species, -1, 0.01)
#   c_a_log <- rnorm(n_species, -3.3, 0.2)
#   c_b_log <- rnorm(n_species, -3, 0.1)
#   c_j_log <- rnorm(n_species, -5, 0.5)
#   g_logit <- rnorm(n_species, -1, 0.2)
#   h_logit <- rnorm(n_species, -0.5, 0.3)
#   m_a_log <- rnorm(n_species, -1.5, 0.5)
#   m_j_log <- rnorm(n_species, -1.2, 0.1)
#   r_log <- rnorm(n_species, 1.5, 0.2)
#   l_log <- rnorm(n_species, -1, 0.2)
#   s_log <- rnorm(n_species, -2.9, 0.1)
#   
#   beta_null <- function() { matrix(c(-1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.001) }
#   
#   
#   inits <- list(
#     
#     # b_log = b_log,
#     # c_a_log = c_a_log,
#     # c_b_log = c_b_log,
#     # c_j_log = c_j_log,
#     # 
#     # g_logit = g_logit,
#     # h_logit = h_logit,
#     # 
#     # l_log = l_log,
#     # m_a_log = m_a_log,
#     # m_j_log = m_j_log,
#     # r_log = r_log,
#     # s_log = s_log,
#     
#     ## from global environment
#     Beta_b = matrix(c(1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species),
#     Beta_c_a = beta_null(),
#     Beta_c_b = beta_null(),
#     Beta_c_j = beta_null(),
#     Beta_g = matrix(c(1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
#     Beta_h = matrix(c(1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
#     Beta_l = matrix(c(1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
#     Beta_m_a = matrix(c(-4, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
#     Beta_m_j = matrix(c(-2, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
#     Beta_r = matrix(c(3, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
#     Beta_s = matrix(c(-3, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
#     
#     b = exp(b_log),
#     c_a = exp(c_a_log),
#     c_b = exp(c_b_log),
#     c_j = exp(c_j_log),
#     
#     g = plogis(g_logit),
#     h = plogis(h_logit),
#     
#     l = exp(l_log),
#     m_a = exp(m_a_log),
#     m_j = exp(m_j_log),
#     r = exp(r_log),
#     s = exp(s_log),
#     
#     # b_loc = matrix(rep(exp(b_log), n_locs), nrow = n_locs, byrow = T),
#     # c_a_loc =  matrix(rep(exp(c_a_log),n_locs), nrow = n_locs, byrow = T),
#     # c_b_loc =  matrix(rep(exp(c_b_log),n_locs), nrow = n_locs, byrow = T),
#     # c_j_loc =  matrix(rep(exp(c_j_log),n_locs), nrow = n_locs, byrow = T),
#     # g_loc =  matrix(rep(plogis(g_logit),n_locs), nrow = n_locs, byrow = T),
#     # h_loc =  matrix(rep(plogis(h_logit),n_locs), nrow = n_locs, byrow = T),
#     # l_loc =  matrix(rep(exp(l_log),n_locs), nrow = n_locs, byrow = T),
#     # m_a_loc =  matrix(rep(exp(m_a_log),n_locs), nrow = n_locs, byrow = T),
#     # m_j_loc =  matrix(rep(exp(m_j_log),n_locs), nrow = n_locs, byrow = T),
#     # r_loc =  matrix(rep(exp(r_log),n_locs), nrow = n_locs, byrow = T),
#     # s_loc =  matrix(rep(exp(s_log),n_locs), nrow = n_locs, byrow = T),
#     
#     shape_par      = c(10, 10, 10),
#     sigma_process  = c(0.01),
#     
#     sigma_obs      = c(1.1, 1.1, 1.1),
#     alpha_obs      = c(1, 1, 1),
#     alpha_obs_inv   = c(1, 1, 1),
#     phi_obs      = c(10, 20, 20),
#     
#     theta = 0.5
#     
#     # state_init = state_init,
#     # state_init_log = log(state_init)
#     
#     # u = replicate(n_locs, matrix(rnorm(n_species*3, 0, 0.001), nrow = pars$n_species*3, ncol = data$timespan_max))
#   )
#   
#   return(inits)
# }




## drawTest --------------------------------
# tar_make("data_stan")
# data_stan <- tar_read("data_stan")
# tar_make("testmodel")
# model <- testmodel <- tar_read("testmodel")

drawTest <- function(model, data_stan, initfunc = 0.5,
                     method = c("mcmc", "variational"), n_chains = 4, iter_warmup = 1000, iter_sampling = 500, # openclid = c(0, 0),
                     fitpath = "Fits.nosync/") {
  
  require(cmdstanr)
  
  if (!dir.exists(fitpath)) {
    dir.create(fitpath)
  }
  
  
  if(match.arg(method) == "variational") {
    
    fit <- model$variational(data = data_stan,
                             output_dir = fitpath,
                             init = initfunc,
                             eta = 0.001,
                             iter = 20**4)
    
  } else if (match.arg(method) == "mcmc") {
    
    ## https://mc-stan.org/cmdstanr/articles/opencl.html
    ## system("clinfo -l")
    
    fit <- model$sample(data = data_stan,
                        output_dir = fitpath,
                        # output_basename = ,
                        init = initfunc,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        # opencl_ids = openclid,
                        # adapt_delta = 0.99,
                        # max_treedepth = 16,
                        chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
  }
  
  return(fit)
}


## draw --------------------------------
# data_stan_prior <- tar_read("data_stan_prior")

draw <- function(model, data_stan, initfunc) {
  
  return(fit)
}



## summarizeFit --------------------------------
# fit <- tar_read("fit")
# fit <- tar_read("fit_test")

summarizeFit <- function(fit) {
  
  summary <- fit$summary()
  
  summarypath <- fit$output_files()[1] %>%
    stringr::str_replace("-1-", "-x-") %>%
    stringr::str_replace(".csv", "_summary.csv")
  
  write.csv(summary, summarypath)
  
  return(summary)
}


## readStanfit --------------------------------
# fit  <- tar_read("fit")
# fit  <- tar_read("fit_test")
## stanfit <- readStanfit(fit, purge = T)

readStanfit <- function(fit, purge = FALSE) {
  
  outputfile <- fit$output_files()
  stanfit <- rstan::read_stan_csv(outputfile)
  
  ## This is to purge draws with NaN values from the model for plotting. NaNs can arise in generated quantities ...
  if (purge) {
    draws <- stanfit@sim$samples
    completerow <- lapply(draws, function(d) complete.cases(d)) %>% as.data.frame() %>% apply(1, all)
    stanfit@sim$samples <- lapply(draws, function(D) {attr(D, "sampler_params") <-  attr(D, "sampler_params")[completerow,]; D[completerow,] })
    
    n_draws_complete <- sum(completerow)
    n_draws <- stanfit@sim$n_save
    diff_complete <- stanfit@sim$n_save - n_draws_complete
    stanfit@sim$n_save <- rep(n_draws_complete, stanfit@sim$chains)
    stanfit@sim$iter <- stanfit@sim$iter - diff_complete[1]
    stanfit@sim$permutation <- lapply(stanfit@sim$permutation, function(p) p[!(p %in% n_draws[1]:n_draws_complete)])
    
    message("There were ", diff_complete[1], " draws with NaNs, that were purged from the fit.")
  }
  
  return(stanfit)
}


## extractDraws --------------------------------
# stanfit  <- tar_read("stanfit")
# stanfit  <- tar_read("stanfit_test")
# helpers_exclude  <- tar_read("helpers_exclude")


extractDraws <- function(stanfit, exclude = helpers_exclude) {
  
  draws <- rstan::extract(stanfit, pars = exclude, include = F)
  
  return(draws)
}


## plotStanfit --------------------------------
# stanfit  <- tar_read("stanfit_plotting")
# stanfit  <- tar_read("stanfit_test_plotting")
# exclude <- tar_read("exclude")

## plotStanfit(stanfit, exclude)

plotStanfit <- function(stanfit, exclude) {
  
  plotRidges <- function(startswith, fit = stanfit) {
    bayesplot::mcmc_areas_ridges(fit, pars = vars(starts_with(startswith, ignore.case = F)))
  }
  
  usedmcmc <- "sample" == attr(stanfit, "stan_args")[[1]]$method
  basename <- attr(stanfit, "model_name") %>%
    str_replace("-[1-9]-", "-x-")
  excludeplus <- c("lp__", "phi_obs")
  parname <- setdiff(stanfit@model_pars, c(exclude, excludeplus))
  parnamestart <- na.omit(unique(str_extract(parname, "^[a-z]_[ljab]"))) # Everything that starts with a small letter, and has the right index after that to be a meaningful parameter.
  parname_sansprior <- parname[!grepl("prior$", parname)]

  traceplot <- rstan::traceplot(stanfit, pars = parname_sansprior, include = T)
  areasplot <- bayesplot::mcmc_areas(stanfit, area_method = "scaled height", pars = vars(!matches(c(exclude, "log_", "lp_", "prior"))))
  ridgeplots <- parallel::mclapply(parnamestart, plotRidges, mc.cores = getOption("mc.cores", 7L))
  ridgeplotgrid <- cowplot::plot_grid(plotlist = ridgeplots)
  
  # parallelplot_c <- bayesplot::mcmc_parcoord(stanfit, pars = vars(starts_with(c("c_", "s_"))))
  # parallelplot_others <- bayesplot::mcmc_parcoord(stanfit, pars = vars(!matches(c(exclude, "c_", "log_", "phi_", "lp_", "s_", "_prior"))))
  
  plots <- list(traceplot = traceplot,
                ridgeplotgrid = ridgeplotgrid,
                areasplot = areasplot) # parallelplot_c = parallelplot_c, parallelplot_others = parallelplot_others,
  
  mapply(function(p, n) ggsave(paste0("Fits.nosync/", basename, "_", n, ".pdf"), p), plots, names(plots))

  if(usedmcmc) {
    
    png(paste0("Fits.nosync/", basename, "_", "pairsplot", ".png"), width = 2600, height = 2600)
    pairs(stanfit, pars = c(parname_sansprior, "phi_obs_inv", "phi_obs", "lp__"), include = T)
    dev.off()
    
  }

  return(plots)
}


## plotDensCheck --------------------------------
# cmdstanfit  <- tar_read("priorsim_test")
# cmdstanfit  <- tar_read("fit_test")
# data_stan_priors <- tar_read("data_stan_priors")
# draws <- tar_read("draws_test") ## this is here as an option for plotting draw objects if the fit has NaNs in generated quantities

plotDensCheck <- function(cmdstanfit, data_stan_priors, draws = NULL, check = c("prior", "posterior")) {
  
  data <- data_stan_priors$y
  pop <- data_stan_priors$rep_pops2y
  
  if(match.arg(check) == "prior") {
    
    if (is.null(draws)) {
      Sim <- cmdstanfit$draws(variables = "y_prior_sim", format = "draws_matrix")
    } else {
      Sim <- draws$y_prior_sim
    }
    
  } else if (match.arg(check) == "posterior") {
    
    if (is.null(draws)) {
      Sim <- cmdstanfit$draws(variables = "y_hat_rep", format = "draws_matrix")
      Fixpoint <- cmdstanfit$draws(variables = "state_fix", format = "draws_matrix")
      fixpointconverged <- cmdstanfit$draws(variables = "converged", format = "draws_matrix")
      
    } else {
      Sim <- draws$y_hat_rep
      
      ## untested:
      Fixpoint <- draws$state_fix
      fixpointconverged <- draws$converged
    }
  }
  
  completerows <- complete.cases(Sim)
  Sim <- Sim[completerows,]
  attr(Sim, "dimnames")$draw <- attr(Sim, "dimnames")$draw[completerows]
  densplot <- bayesplot::ppc_dens_overlay_grouped(log(data), log(Sim), group = pop)
  
  if (match.arg(check) == "posterior") {
    Fixpoint <- Fixpoint[completerows,]
    attr(Fixpoint, "dimnames")$draw <- attr(Sim, "dimnames")$draw[completerows]
    # fixpointconverged <- fixpointconverged[completerows,]
    # attr(fixpointconverged, "dimnames")$draw <- attr(fixpointconverged, "dimnames")$draw[completerows]
    popstatesinfixpoint <- rep(c(rep(T, data_stan_priors$N_pops + data_stan_priors$N_species), rep(F, data_stan_priors$N_species + 1)), data_stan_priors$N_locs)
    Fixpoint <- Fixpoint[, popstatesinfixpoint]
    attr(Fixpoint, "dimnames")$variable <- rep(c(paste("log pop", 1:data_stan_priors$N_pops), paste("log ba", 1:data_stan_priors$N_species)), data_stan_priors$N_locs)
    
    fixdensplot <- bayesplot::mcmc_areas_ridges(log(Fixpoint))
    densplot <- cowplot::plot_grid(densplot, fixdensplot, labels = c("States", "Equilibria")) # , ncol = 1, axis = "b", align = "h"
  }

  
  basename <- cmdstanfit$metadata()$model_name %>%
    str_replace("-[1-9]-", "-x-")
  name <- paste0("dens_", check)
  ggsave(paste0("Fits.nosync/", basename, "_", name, ".pdf"), densplot)
  
  return(densplot)
}


## scaleResiduals --------------------------------
# cmdstanfit  <- tar_read("fit_test")
# cmdstanfit  <- tar_read("fit")
# data_stan_priors  <- tar_read("data_stan_priors")

scaleResiduals <- function(cmdstanfit, data_stan_priors) {
  
  Sim <- cmdstanfit$draws(variables = "y_sim", format = "draws_matrix") %>% t()# matrix of observations simulated from the fitted model - row index for observations and colum index for simulations
  Sim[is.na(Sim)] <- 0
  y <- data_stan_priors$y
  y_hat <- cmdstanfit$draws(variables = "y_hat_rep", format = "draws_matrix") %>% apply(2, median, na.rm = T)
  y_hat[is.na(y_hat)] <- 0
  
  residuals <- DHARMa::createDHARMa(simulatedResponse = Sim, observedResponse = y, fittedPredictedResponse = y_hat, integerResponse = T)

  basename <- cmdstanfit$metadata()$model_name %>%
    str_replace("-[1-9]-", "-x-")
  
  png(paste0("Fits.nosync/", basename, "_", "DHARMa", ".png"), width = 2000, height = 1200)
  plotResiduals(residuals, quantreg = T, smoothScatter = F)
  dev.off()
  
  return(residuals)
}


# ——————————————————————————————————————————————————————————————————————————————————#
# Seedlings: Functions for fitting a regeneration model to SK seedling data --------
# ——————————————————————————————————————————————————————————————————————————————————#

## wrangleSeedlings --------------------------------
# Data_seedlings  <- tar_read("Data_seedlings")
# taxon_select <- tar_read("taxon_select")
# threshold_dbh <- tar_read("threshold_dbh")

wrangleSeedlings <- function(Data_seedlings, taxon_select = taxon_select, threshold_dbh = threshold_dbh) {
  
  if (taxon_select != "Fagus.sylvatica") stop("Prior for seedling regeneration rate r is only implemented for Fagus.sylvatica!")

  Data_seedlings <- Data_seedlings %>%
    mutate(tax = str_replace_all(taxon, " ", replacement = ".")) %>%
    mutate(tax = forcats::fct_other(tax, !!!as.list(taxon_select), other_level = "other")) %>%
    mutate(taxid = if_else(tax == "Fagus.sylvatica", "kew.83891", "other")) %>%
    
    # mutate(taxid = if_else(tax == "other", "other", as.character(taxid))) %>% ## if not only Fagus would be used
    droplevels()
    
  D_select <- Data_seedlings %>%
    group_by(plotid, year) %>%
    mutate(drop = any(regeneration == "Artificial", na.rm = T)) %>%
    ungroup() %>%
    filter(!drop)
    ## Any observations on plots were removed that had unnatural regeneration.
  
  D_filtered <- D_select %>%
    filter((sizeclass == "big" & dbh >= 100) |
           (sizeclass == "small" & height < 0.2) |
            is.na(sizeclass)) %>%
    filter(status != "Dead tree" | is.na(status)) %>%
    filter(!is.na(tax))
    ## The species-specific basal area of a plot was confined to living trees above the dbh >= 100mm.
    ## Only small trees with size class [10cm, 20cm), were included in the seedling counts.
  
  
  D_count <- D_filtered %>%
    group_by(plotid, year, tax, taxid) %>%
    dplyr::summarize(count_ha = sum(count_ha, na.rm = T),
                     ba_ha = sum(ba_ha, na.rm = T)) %>%
    mutate(count_ha = as.integer(round(count_ha)))
    
  
  # library(glmmTMB)
  # m <- glmmTMB::glmmTMB(count_ha ~ ba_ha * tax + 0, data = D_count, family = nbinom2)
  # summary(m)
  # res <- DHARMa::simulateResiduals(m)
  # plot(res)
  
  return(D_count)
}


## fitSeedlings --------------------------------
# Seedlings  <- tar_read("Seedlings")

fitSeedlings <- function(Seedlings) {
  
  fit_seedlings <- brms::brm(count_ha ~ ba_ha, data = Seedlings[Seedlings$tax == "Fagus.sylvatica",], family = negbinomial,
                             cores = getOption("mc.cores", 4))
  
  fit_seedlings_other <- brms::brm(count_ha ~ ba_ha, data = Seedlings[Seedlings$tax == "other",], family = negbinomial,
                             cores = getOption("mc.cores", 4))

  message("Summary of the the fit for Fagus seedlings:")
  print(summary(fit_seedlings))
  
  message("Summary of the the fit for other seedlings:")
  print(summary(fit_seedlings_other))
  
  return(list(Fagus.sylvatica = fit_seedlings, other = fit_seedlings_other))
  
}

