# ——————————————————————————————————————————————————————————————————————————————————#
# Functions for fitting  ----------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————#

## constructPriors --------------------------------
# draws_direct  <- tar_read("draws_direct")

constructPriors <- function(draws_direct) {
  
  x <- log(rgamma(150,5))
  df <- approxfun(density(x))
  plot(density(x))
  xnew <- seq(-1, 3, by = 0.1)
  points(xnew,df(xnew),col=2)
  
  return(fit)
}



## formatStanData --------------------------------
# Stages  <- tar_read("Stages_scaled")
# taxon_select  <- tar_read("taxon_select")
# threshold_dbh <- tar_read("threshold_dbh")

formatStanData <- function(Stages, Rates, taxon_select, threshold_dbh) { # priors!
  
  taxon_selectother <- c(taxon_select, "other")
  
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
  Stages %<>%
    
    ## Synonyms for consistency wiht model lingo
    group_by(clusterid) %>%
    mutate(plot = match(plotid, unique(plotid))) %>% ## for numbers designating which corner it is, extract the substring at the end of plotid
    ungroup() %>%
    mutate(loc = as.integer(as.factor(clusterid))) %>%
    
    ## factor ordering!!!
    mutate(stage = factor(stage, levels = c("J", "A", "B", "BA")),
           tax = factor(tax, levels = taxon_selectother)
                          ) %>%
    
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
    select(-isy0) %>%
    group_by(loc, pop, stage, tax) %>%
    summarize(n_plots = n_distinct(plot), .groups = "drop")
  
  S <- select(S, -isy0) %>%
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
  
  #### Prepare design matrix
  Env <- S_locs[c("phCaCl_esdacc_s", "waterLevel_loc_s")]
  polyformula <- as.formula(paste("~", paste("poly(", colnames(as.data.frame(Env)), ", 2)", collapse = "+")))
  X <- model.matrix(polyformula, data = as.data.frame(Env))
  
  #### Prepare ldd smooth
  L_smooth <- S %>%
    group_by(loc) %>%
    summarize_at(paste("s", taxon_selectother, sep = "_"), first) %>%
    column_to_rownames(var = "loc") %>%
    as.matrix()
  
  #### Some data for multiple reuse in the list
  N_species <-  length(unique(S$tax))
  L_yhat <- nrow(S_yhat)
  
  
  #### The stan-formatted list
  data <- list(
    
    L_times = nrow(S_times),
    L_yhat = L_yhat,
    L_y = nrow(S),
    
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
    
    time_max = S_locs$time_max,
    times = S_times$t,

    X = X,
    L_smooth = L_smooth, ## array[N_locs] vector<lower=0>[N_species] L_smooth;
    
    y = S$y,
    
    ## Settings corner
    tolerance_fix = 0.001,
    ba_a_upper = ba_a_upper,
    ba_a_avg = ba_a_avg
    
    # Long = S
  )
  
  return(data)
}


#### Returns viable start values ---------------------------
getInits <- function() {
  
  responsescaleerror <- 0.1
  data <- tar_read(data_stan) ## This is fine in this case, as getInits is only called inside targets that would be invalidated after change of data_stan
  n_species <- data$N_species
  n_beta <- data$N_beta
  n_locs <- data$N_locs
  
  state_init <- matrix(rlnorm(data$L_yhat, log(0.01), 0.01), nrow = n_locs, ncol = data$N_pops) # array[N_locs] vector<lower=0>[N_pops] state_init;
  
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
  
  beta_null <- function() { matrix(c(-1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.001) }
  
  
  inits <- list(
    
    # b_log = b_log,
    # c_a_log = c_a_log,
    # c_b_log = c_b_log,
    # c_j_log = c_j_log,
    # 
    # g_logit = g_logit,
    # h_logit = h_logit,
    # 
    # l_log = l_log,
    # m_a_log = m_a_log,
    # m_j_log = m_j_log,
    # r_log = r_log,
    # s_log = s_log,
    
    ## from global environment
    Beta_b = matrix(c(1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species),
    Beta_c_a = beta_null(),
    Beta_c_b = beta_null(),
    Beta_c_j = beta_null(),
    Beta_g = matrix(c(1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
    Beta_h = matrix(c(1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
    Beta_l = matrix(c(1, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
    Beta_m_a = matrix(c(-4, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
    Beta_m_j = matrix(c(-2, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
    Beta_r = matrix(c(3, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
    Beta_s = matrix(c(-3, rep(0, n_beta-1)), ncol = n_species, nrow = n_beta) + rnorm(n_beta*n_species, 0, 0.3),
    
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
    
    # b_loc = matrix(rep(exp(b_log), n_locs), nrow = n_locs, byrow = T),
    # c_a_loc =  matrix(rep(exp(c_a_log),n_locs), nrow = n_locs, byrow = T),
    # c_b_loc =  matrix(rep(exp(c_b_log),n_locs), nrow = n_locs, byrow = T),
    # c_j_loc =  matrix(rep(exp(c_j_log),n_locs), nrow = n_locs, byrow = T),
    # g_loc =  matrix(rep(plogis(g_logit),n_locs), nrow = n_locs, byrow = T),
    # h_loc =  matrix(rep(plogis(h_logit),n_locs), nrow = n_locs, byrow = T),
    # l_loc =  matrix(rep(exp(l_log),n_locs), nrow = n_locs, byrow = T),
    # m_a_loc =  matrix(rep(exp(m_a_log),n_locs), nrow = n_locs, byrow = T),
    # m_j_loc =  matrix(rep(exp(m_j_log),n_locs), nrow = n_locs, byrow = T),
    # r_loc =  matrix(rep(exp(r_log),n_locs), nrow = n_locs, byrow = T),
    # s_loc =  matrix(rep(exp(s_log),n_locs), nrow = n_locs, byrow = T),
    
    shape_par      = c(10, 10, 10),
    sigma_process  = c(0.01),
    
    sigma_obs      = c(1.1, 1.1, 1.1),
    alpha_obs      = c(1, 1, 1),
    alpha_obs_inv   = c(1, 1, 1),
    phi_obs      = c(10, 20, 20),
    
    theta = 0.5
    
    # state_init = state_init,
    # state_init_log = log(state_init)
    
    # u = replicate(n_locs, matrix(rnorm(n_species*3, 0, 0.001), nrow = pars$n_species*3, ncol = data$timespan_max))
  )
  
  return(inits)
}




## drawTest --------------------------------
# tar_make("data_stan")
# data_stan <- tar_read("data_stan")
# tar_make("testmodel")
# model <- testmodel <- tar_read("testmodel")


drawTest <- function(model, data_stan, initfunc = 0.5,
                     method = c("mcmc", "variational"), n_chains = 3, iter_warmup = 200, iter_sampling = 300,
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

    fit <- model$sample(data = data_stan,
                        output_dir = fitpath,
                        # output_basename = ,
                        init = initfunc,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = 0.99,
                        max_treedepth = 16,
                        chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
  }
  
  return(fit)
}


## draw --------------------------------
# data_stan <- tar_read("data_stan")

draw <- function(model, data_stan) {
  
  return(fit)
}


## summarizeDraws --------------------------------
# fit  <- tar_read("fit")
# fit  <- tar_read("testfit")

summarizeDraws <- function(fit) {
  
  summary <- fit$summary()
  
  summarypath <- fit$output_files()[1] %>%
    stringr::str_replace("-1-", "-x-") %>%
    stringr::str_replace(".csv", "_summary.csv")
  
  write.csv(summary, summarypath)
  
  return(summary)
}


## extractDraws --------------------------------
# fit  <- tar_read("fit")
# testfit  <- tar_read("testfit")

extractDraws <- function(fit) {
  
  outputfile <- fit$output_files()
  draws <- rstan::read_stan_csv(outputfile) %>%
    rstan::extract()
  
  return(draws)
}


