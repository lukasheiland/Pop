# ——————————————————————————————————————————————————————————————————————————————————#
# Functions for fitting  ----------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————#

## formatStanData --------------------------------
# Stages  <- tar_read("Stages_loc")
# Stages_transitions  <- tar_read("Stages_transitions")
# taxon_s  <- tar_read("taxon_s")
# predictor_select  <- tar_read("predictor_select")
# threshold_dbh <- tar_read("threshold_dbh")
# loclevel <- tar_read(loc)
# timestep <- 1

formatStanData <- function(Stages, Stages_transitions, taxon_s, threshold_dbh, predictor_select,
                           loc = c("plot", "nested", "cluster"), timestep = 1,
                           smoothmediantax = c("none", "other", "Fagus.sylvatica")) {
  
  predictor_select_s <- paste0(predictor_select, "_s")
  loclevel <- match.arg(loc)
  
  ## Prepare ba_a_upper, the upper basal area of any one tree in A.
  radius_a_upper <- threshold_dbh/2
  ba_a_upper <-  pi * radius_a_upper^2 * 1e-6
  
  ## reparameterizeModeGamma
  ## http://doingbayesiandataanalysis.blogspot.com/2012/01/parameterizing-gamma-distribution-by.html
  ## shape == alpha; inverse scale == rate == beta
  ## v <- reparameterizeModeGamma(200, 100); curve(dgamma(x, v["alpha"] , v["beta"]), 0, 300)
  reparameterizeModeGamma <- function(mode, sd) {
    rate <- ( mode + sqrt( mode^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
    shape <- 1 + mode * rate
    return(c(alpha = shape, beta = rate))
  }
  
  #### Prepare ba_a_avg
  ## vector[N_species] for ba_a_avg
  ## For multiplication with count A, ba_a_avg is the average basal area of an individual tree in class A and has unit $\mathit{areaunit} \cdot 1^-1$, so that A * ba_a_avg has unit $\mathit{areaunit} \cdot ha^-1$.
  ba_a_avg <- Stages %>%
    st_drop_geometry() %>%
    filter(stage == "A") %>%
    group_by(tax) %>%
    summarize(ba_a_avg = mean(ba_obs/count_obs, na.rm = T)) %>% ## [(ba/ha)/(1/ha) == ba/ha * ha == ba/1]
    pull(ba_a_avg, name = tax)
  
  
  #### Prepare special long data set structures with different amounts of long casting along different ids
  ## The most comprehensive data set is L_y with grouping locations/resurveys/pops/plots.
  ## NOTE THAT plots HAS TO BE THE LAST GROUP IN SORTING
  ## Everything is subset from the master subsets with these groupings (*_reobs, *_y0) and thus consistently sorted.
  ## Factor "pops" is structured stages/species.
  
  Stages_transitions %<>%
    ungroup() %>%
    dplyr::select("plotid", "tax", "obsid",
                  # "timediff_plot", "count_A_sum_before", "count_J_sum_before",
                  # "h_plot", "g_plot"
                  "count_A2B_plot", "count_J2A_plot", ## [1/ha]
                  "count_A2B_plot_obs", "count_J2A_plot_obs", ## true observed ...
                  "area_A2B", "area_J2A", ## on this area
                  "count_J_integr_plot", "count_A_integr_plot") %>%
    
    group_by(plotid, obsid, tax) %>%
    slice(1)
  
  G <- Stages_transitions %>%
    # filter(!is.na(g_plot))
    filter(!is.na(count_J2A_plot_obs) & !is.na(count_J_integr_plot)) %>% ## also drops NAs
    
    ### Version with y_base from data:
    filter(!is.na(count_J2A_plot_obs) & isTRUE(count_J_integr_plot > 0)) ## also drops NAs
  
  H <- Stages_transitions %>%
    # filter(!is.na(h_plot))
    filter(!is.na(count_A2B_plot_obs) & !is.na(count_A_integr_plot)) %>% ## also drops NAs
    
    ### Version with y_base from data:
    filter(!is.na(count_A2B_plot_obs) & isTRUE(count_A_integr_plot > 0)) ## also drops NAs
  
  
  if(loclevel == "nested") {
    
    Stages %<>%
      ## Synonyms for consistency with model lingo
      group_by(clusterid) %>%
      mutate(plot = match(plotid, unique(plotid))) %>% ## for numbers designating which corner it is, extract the substring at the end of plotid
      ungroup()
  }
  
  Stages %<>%
    
    ## Join Stages_transitions (for direct use in the model on the loc level)
    # bind_cols(Stages_transitions[match(interaction(.$plotid, .$tax), Stages_transitions$joinid), ]) %>%
    
    ## Synonyms for consistency with model lingo
    mutate(loc = as.integer(as.factor(loc))) %>%
    
    ## factor ordering!!!
    mutate(stage = factor(stage, levels = c("J", "A", "B", "BA")),
           tax = factor(tax, levels = levels(taxon_s))) %>%
    
    arrange(loc, time, stage, tax, plotid) ## loc might be == plotid
  
  ## Generate ids
  Stages %<>% 
    ## Stages are measured in different terms: ba or count
    mutate(y = case_when(
      stage == "J" ~ count_obs, # count_ha_r,
      stage == "A" ~ count_obs, # count_ha_r,
      stage == "B" ~ count_obs, # ba_ha_r
      stage == "BA" ~ as.double(ba_ha_r)
    )) %>%
    
    mutate(y_prior = case_when(
      stage == "J" ~ count_ha, # count_ha_r,
      stage == "A" ~ count_ha, # count_ha_r,
      stage == "B" ~ as.double(ba_ha), # ba_ha_r
      stage == "BA" ~ as.double(ba_ha)
    )) %>%
    
    ## pop is just an id for the initial states vector
    mutate(pop = interaction(tax, stage)) %>%
    
    ## Different levels of observation error were assumed ...
    ## - different obsmethod for J (area count sampling), the stage A (counts from fixed angle sampling), and stage B (basal area from fixed angle sampling).
    mutate(obsmethod = fct_recode(stage, "j" = "J", "a" = "A", "ba" = "B", "ba" = "BA")) %>% # obsmethodTax == pop
    ## - different methods for A, and B in 1987 vs. 2002/2012 and different in J over all iterations (different areas)
    mutate(protocol = methodid) %>%
    mutate(protocolTax = interaction(substr(tax, 1, 1), stage, methodid)) %>% ##! stage has to be included here!
    ## - per species and per survey
    # mutate(obsidPop = interaction(obsid, pop)) %>%
    ## - per species and per method (BWI1 vs. BWI2/3)
    mutate(obsidPop = interaction(if_else(obsid == "DE_BWI_1987", "init", "later"), pop)) %>%
    
    ## scale times within group to start with 1
    group_by(loc) %>%
    mutate(t = as.integer(lubridate::year(time))) %>%
    mutate(t = t - min(t) + 1) %>%
    mutate(isy0 = (t == 1)) %>%
    
    mutate(t1 = t, t_min = min(t)) %>%
    mutate(t = round((t - t_min)/timestep) + t_min) %>%
    mutate(res_t = (t - t_min)*timestep + t_min - t1) %>%
    
    mutate(t_which = match(t, sort(unique(t)))) %>%
    ungroup() %>%
    
    mutate(yhat2y = interaction(loc, t_which, pop, drop = T))
  
  if (loclevel == "nested") {
    ## Generate variables for "nested" structure, e.g. yhat2y matching
    Stages %<>% 
      group_by(yhat2y) %>%
      mutate(y_hat_prior = mean(y_prior, na.rm = T)) %>% ## For state debugging, per ha
      ungroup() %>%
      arrange(loc, time, pop, stage, tax, plot)
    
  } else {
    
    Stages %<>% 
      mutate(y_hat_prior = y_prior) %>%
      arrange(loc, time, pop, stage, tax)
  }
  
  
  #### Here: different subsets of the data in "Stages", with different formats to be used at different points in the model
  
  ## Format: [L_y] — locations/obs/pops(/stage/species)/(plots)
  ## This is the same as Stages, but without BA. Only stages J, A, B for fitting
  S <- Stages %>%
    st_drop_geometry %>%
    filter(stage %in% c("J", "A", "B")) %>%
    mutate_at(c("stage", "pop", "obsidPop", "protocolTax"), droplevels) ## !!!
  
  ## Format: [L_noninit] — locations/resurveys/pops
  # S_noninit <- bind_cols(S, yhat2y = 1:nrow(S)) %>%
  #   filter(!isy0)
  
  ## Format: [L_init] — locations/pops
  S_init <- filter(S, isy0) %>%
    dplyr::select(-isy0) %>%
    group_by(loc, pop, stage, tax) %>%
    summarize(n_plots = n_distinct(plotid), .groups = "drop") ## will only be used in "nested", if at all
  
  ## Format: [L_yhat] — locations/obs/pops  ## will only be used in "nested"
  S_yhat <- S %>%
    group_by(loc, time, pop, stage, tax) %>%
    summarize(n_plots = n_distinct(plotid),
              y_hat_prior = first(y_hat_prior), ## == y_prior for anything but loclevel == "nested"
              yhat2y = first(yhat2y),
              .groups = "drop")
  
  ## Format: [L_times] — locations/resurveys
  S_times <- S %>%
    group_by(loc, obsid) %>% # first, group by obsid to extract one time per obsid within loc (in case there are multiple days somehow, which is not the case atm!)
    summarize(t = first(t), t1 = first(t1), res_t = first(res_t), .groups = "drop") %>%
    group_by(loc) %>%
    mutate(t = unique(t), .groups = "drop") ## This way its certain, that there is one time per obsid
  
  if (timestep != 1) {
    Timetable <- table(S_times$res_t)
    message("Table of residual times in years throug timestep reduction: ")
    write.csv("Publish/Residual_times.csv")
    print(Timetable)
  }
  
  ## Format: [N_locs] — locations
  S_locs <- S %>%
    group_by(loc) %>%
    ## assumes completion within locations!
    summarize(time_max = max(t),
              timespan = diff(range(t)) + 1,
              n_tax = n_distinct(tax),
              n_plots = n_distinct(plotid),
              n_pops = n_distinct(pop),
              n_obs = n_distinct(t),
              n_yhat = n_distinct(interaction(pop, t)),
              
              across(all_of(predictor_select_s), first),
              
              clusterid = first(clusterid),
              .groups = "drop")
  
  
  ## Format: [L_init] — locations/pops
  # S_a2b <- S %>%
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
  Env <- S_locs[predictor_select_s]
  if (anyNA(Env)) {
    X <- matrix(0, nrow = nrow(Env), ncol = ncol(Env))
    message("With Env having NAs, an empty design matrix was produced, assuming that it will not be used.")
  } else {
    polyformula <- as.formula(paste("~", paste("poly(", colnames(as.data.frame(Env)), ", 2, raw = T)", collapse = "+")))
    X <- model.matrix(polyformula, data = as.data.frame(Env))
  }
  
  
  #### Prepare ldd smooth. Predicted with log-link 
  smt <- match.arg(smoothmediantax)
  L_smooth_log <- S %>%
    group_by(loc) %>%
    summarize_at(paste("s", taxon_s, sep = "_"), first) %>%
    ungroup() %>%
    { if (!smt == "none") mutate_at(., paste("s", smt, sep = "_"), .funs = median)  else . } %>%
    tibble::column_to_rownames(var = "loc") %>%
    as.matrix()
  
  #### Some data for multiple reuse in the list
  N_species <-  length(unique(S$tax))
  L_yhat <- nrow(S_yhat)
  
  #### Prepare some summary statistics that inform priors
  # Phi_empirical <- S %>%
  #   filter(stage %in% c("J", "A", "B")) %>%
  #   # group_by(pop) %>%
  #   group_by(protocolTax) %>%
  #   summarize(var = var(y_prior, na.rm = T), mean = mean(y_prior, na.rm = T), phi = mean^2/(var - mean)) %>%
  #   mutate(phi_inv = 1/phi) %>%
  #   mutate(phi_inv_sqrt = 1/sqrt(phi)) %>%
  #   mutate(sigma_phi = signif(phi_inv_sqrt, digits = 1)) %>%
  #   # arrange(pop)
  #   arrange(protocolTax)
  # 
  # message("Empirical estimates of phi per species and protocol:")
  # print(Phi_empirical)
  
  upper_init <- S %>%
    group_by(pop) %>%
    summarize(upper = max(y_prior, na.rm = T) * 1.5) %>%
    pull(upper, name = pop)
  
  Y_init <-  filter(S, isy0) %>%
    group_by(pop) %>%
    ## Together with the offset, the minimum observation is always == 1! This way we construct a prior for the zeroes, that has the most density around zero, but an expected value at 1, assuming that 5% of the zero observations are actually wrong.
    mutate(min_pop = case_when(stage == "J" ~ min(y_prior[y_prior != 0], na.rm = T) * 0.02,
                               stage == "A" ~ min(y_prior[y_prior != 0], na.rm = T) * 0.01,
                               stage == "B" ~ min(y_prior[y_prior != 0], na.rm = T) * 0.005)) %>%
    
    group_by(loc, pop) %>%
    ## The summaries here are only effectual for loclevel == "nested", because otherwise the grouping group_by(loc, pop) is identical to the original id structure 
    summarize(lower = min(y_prior, na.rm = T) * 0.01, ## the lower boundary, 1% of the minimum value, is not used in the model
              y_prior = if (loclevel == "nested") first(y_hat_prior) else y_prior,
              min_pop = first(min_pop),
              
              y_prior_0 = if_else(y_prior == 0, min_pop, y_prior),
              
              ## Setting alpha = 1 for 0, so that the most density is towards zero
              ## parameters of the prior gamma density, alpha and beta, are computed based on the data.
              ## Computation happens for two different parameterizations: by mean and by model.
              # alpha = if_else(y_prior == 0, 1, 10),
              alpha_mean = case_when(stage == "J" ~ 1 + as.integer(count_obs > 0) * 9, # 4 + 2 * count_obs, ## this will assign 1 to count_obs == 0
                                     stage == "A" ~ 1 + as.integer(count_obs > 0) * 9, # 4 + 10 * count_obs,
                                     stage == "B" ~ 1 + as.integer(count_obs > 0) * 9  # 9 + 10 * count_obs
                                     ),
              beta_mean = alpha_mean/y_prior_0, ## this is equivalent to beta
              
              sd_mode = case_when(stage == "J" ~ 50 + 30 * count_obs, ## median 3000, mean 9000, min 200
                                  stage == "A" ~ 5 + 3 * count_obs, ## min prob around 80, long tail
                                  stage == "B" ~ 0.5 + 0.3 * count_obs), ## median 20, mean 22, min 4, short tail
              alpha_mode = reparameterizeModeGamma(y_prior_0, sd_mode)["alpha"],
              beta_mode = reparameterizeModeGamma(y_prior_0, sd_mode)["beta"],
              ## The priors concern the N/ha
              ## example, minimum value for A:  v <- reparameterizeModeGamma(130, 20); curve(dgamma(x, v["alpha"] , v["beta"]), 0, 200)
              ## example, mean value for J:  v <- reparameterizeModeGamma(200, 200); curve(dgamma(x, v["alpha"] , v["beta"]), 0, 600)
              
              alpha = if_else(y_prior == 0, 1, alpha_mode), ## for case zero, alpha_mode is like in alpha_mean case
              beta = if_else(y_prior == 0, 1/min_pop, beta_mode),
              
              ## alphaByE = alpha/y_prior_0, ## this is equivalent to beta in the mean parameterization
              .groups = "drop")
  
  Lower_init <- Y_init %>%
    dplyr::select(pop, loc, lower) %>%
    pivot_wider(names_from = pop, values_from = lower, id_cols = c("loc")) %>%
    arrange(loc) %>%
    tibble::column_to_rownames(var = "loc")
  
  State_init <- Y_init %>%
    dplyr::select(pop, loc, y_prior) %>%
    pivot_wider(names_from = pop, values_from = y_prior, id_cols = c("loc")) %>%
    arrange(loc) %>%
    tibble::column_to_rownames(var = "loc")
  
  Alpha_init <- Y_init %>%
    dplyr::select(pop, loc, alpha) %>%
    pivot_wider(names_from = pop, values_from = alpha, id_cols = c("loc")) %>%
    arrange(loc) %>%
    tibble::column_to_rownames(var = "loc")
  
  Beta_init <- Y_init %>%
    dplyr::select(pop, loc, beta) %>%
    pivot_wider(names_from = pop, values_from = beta, id_cols = c("loc")) %>%
    arrange(loc) %>%
    tibble::column_to_rownames(var = "loc")
  
  ## Zero
  # curve(dgamma(x, Alpha_init[1,1], Beta_init[1,1]), 0, 10)
  ## High J observation
  # curve(dgamma(x, Alpha_init[10,1], Beta_init[10,1]), 0, 8000)
  
  ## State debugging
  # State_1 <- filter(S, isy0) %>%
  #   filter(stage %in% c("J", "A", "B")) %>%
  #   group_by(pop, loc) %>%
  #   summarize(s = mean(y_prior, na.rm = T)) %>%
  #   ungroup() %>%
  #   pivot_wider(names_from = pop, values_from = s) %>%
  #   arrange(loc) %>%
  #   tibble::column_to_rownames(var = "loc") %>%
  #   add(1e-12)
  # 
  # State_2 <- filter(S, t_which == 2) %>%
  #   filter(stage %in% c("J", "A", "B")) %>%
  #   group_by(pop, loc) %>%
  #   summarize(s = mean(y_prior, na.rm = T)) %>%
  #   ungroup() %>%
  #   pivot_wider(names_from = pop, values_from = s) %>%
  #   arrange(loc) %>%
  #   tibble::column_to_rownames(var = "loc") %>%
  #   add(1e-12)
  # 
  # State_3 <- filter(S, t_which == 3) %>%
  #   filter(stage %in% c("J", "A", "B")) %>%
  #   group_by(pop, loc) %>%
  #   summarize(s = mean(y_prior, na.rm = T)) %>%
  #   ungroup() %>%
  #   pivot_wider(names_from = pop, values_from = s) %>%
  #   bind_rows(data.frame(loc = setdiff(as.integer(rownames(State_1)), .$loc),
  #                        Fagus.sylvatica.J = 0, other.J = 0, Fagus.sylvatica.A = 0, other.A = 0, Fagus.sylvatica.B = 0, other.B = 0)) %>%
  #   arrange(loc) %>%
  #   tibble::column_to_rownames(var = "loc") %>%
  #   add(1e-12)
  
  #### Prior method of determining the log mean as an informative prior for initial state
  # Prior_state_init_log <- filter(S, isy0) %>%
  #   filter(stage %in% c("J", "A", "B")) %>%
  #   group_by(pop, loc) %>%
  #   summarize(sd = sd(log(y_prior)), y_prior = log(mean(y_prior, na.rm = T))) %>% ## log(mean) is imperfect and introduces a small bias, median won't work because of the many zeroes
  #   ungroup() %>%
  #   dplyr::mutate(y_prior = na_if(y_prior, -Inf)) %>%
  #   group_by(pop) %>%
  #   dplyr::mutate(y_prior_min = min(y_prior, na.rm = T)) %>%
  #   dplyr::mutate(y_prior = dplyr::coalesce(y_prior, y_prior_min)) %>%
  #   dplyr::select(-y_prior_min, -sd) %>%
  #   ungroup() %>%
  #   pivot_wider(names_from = pop, values_from = y_prior) %>%
  #   arrange(loc) %>%
  #   dplyr::select(-loc) %>%
  #   as.matrix()
  
  
  #### The stan-formatted list
  data <- list(
    
    L_times = nrow(S_times),
    L_yhat = L_yhat,
    L_init = nrow(S_init),
    L_y = nrow(S),
    # L_noninit = nrow(S_noninit),
    # L_a2b = nrow(S_a2b),
    
    N_locs = nrow(S_locs),
    N_species = N_species,
    N_pops = length(unique(S$pop)),
    N_beta = ncol(X),
    N_env = length(predictor_select_s),
    N_obsmethod = length(unique(S$obsmethod)),
    N_protocol = length(unique(S$protocol)), ## different sampling area levels for J even in 2002 and 2012
    N_protocolTax = length(unique(S$protocolTax)), ## different sampling area levels for J even in 2002 and 2012
    N_obsidPop = length(unique(S$obsidPop)),
    N_protocolPop = length(unique(S$methodid)), ## different sampling area levels for J even in 2002 and 2012
    
    n_obs = S_locs$n_obs,
    n_yhat = S_locs$n_yhat,
    
    i_j = 1:N_species,
    i_a = (1:N_species) + N_species,
    i_b = (1:N_species) + 2*N_species,
    
    rep_yhat2y = match(S$yhat2y, S_yhat$yhat2y),
    rep_obsmethod2y = as.integer(S$obsmethod),
    rep_protocol2y = as.integer(S$protocol),
    rep_protocolTax2y = as.integer(S$protocolTax),
    rep_obsidPop2y = as.integer(S$obsidPop), ## S$obsidPop %>% levels()
    rep_pops2y =  as.integer(S$pop),
    rep_pops2init =  as.integer(S_init$pop), ## S_init$pop %>% levels()
    ## unused reps
    # rep_noninit2y = S_noninit$yhat2y,
    # rep_yhat2a2b = S_a2b$rep_yhat2a2b,
    # rep_species2a2b = S_a2b$rep_species2a2b,

    # sigma_phi = Phi_empirical$sigma_phi,
    
    time_max = S_locs$time_max,
    times = S_times$t,
    # timediff = S_a2b$timediff,
    
    X = X,
    L_smooth_log = L_smooth_log, ## array[N_locs] vector<lower=0>[N_species] L_smooth_log;
    L_smooth = exp(L_smooth_log),
    
    y = as.integer(S$y),
    offset_data = S$offset,
    # a2b = S_a2b$a2b,
    
    ## Settings corner
    tolerance_fix = 0.001,
    ba_a_upper = ba_a_upper,
    ba_a_avg = ba_a_avg,
    generateposteriorq = 0,
    
    
    ##$ Altering the timestep
    # timestep = timestep,
    # parfactor = parfactor,
    
    
    #### priors
    # Prior_state_init_log = Prior_state_init_log,
    # Prior_state_init = Prior_state_init,
    
    state_init_data = State_init + 1e-12,
    upper_init = upper_init,
    lower_init = Lower_init,
    alpha_init = Alpha_init,
    beta_init = Beta_init,
    
    ## State debugging
    # state_init = State_1,
    # state_2 = State_2,
    # state_3 = State_3,
    y_hat_prior = S_yhat$y_hat_prior,
    
    ## the transitions.
    L_g = nrow(G),
    L_h = nrow(H),
    y_j2a = G$count_J2A_plot_obs, # integer
    y_a2b = H$count_A2B_plot_obs, # integer
    area_log_j2a = log(G$area_J2A),
    area_log_a2b = log(H$area_A2B),
    
    ## Version with y_base as latent population fitted to counts:
    # y_j = as.integer(round(G$count_J_integr_plot)), # [1/ha]
    # y_a = as.integer(round(H$count_A_integr_plot)), # # [1/ha]
    
    ## Version with y_base as data:
    y_j = G$count_J_integr_plot, # [1/ha]
    y_a = H$count_A_integr_plot, # # [1/ha]
    
    species_g = as.integer(factor(G$tax, levels = levels(taxon_s))),
    species_h = as.integer(factor(H$tax, levels = levels(taxon_s)))
  )
  
  if (!all(data$n_obs * data$N_pops == data$n_yhat)) message("Unexpected lengths of y_hat per locations. Assuming completion of all possible taxa/stages within plot/times went wrong.")
  
  attr(data, "Long") <- S
  attr(data, "Long_BA") <- Stages
  # attr(data, "Phi_empirical") <- Phi_empirical
  message("Phi levels are:", levels(S$protocolTax))
  
  return(data)
}



## fitTransition --------------------------------
# data_stan  <- tar_read("data_stan")
# which  <- "h" # "g"
# model_transitions  <- tar_read("model_transitions")
# fitpath <- tar_read("dir_fit")
fitTransition <- function(data_stan, which, model_transitions, prior_rate = c(g = -6.0, h = -3.0), fitpath = dir_fit) { # priors!
  
  isg <- (which == "g")
  
  d <- with(data_stan, list(
    ## the transitions
    L = if(isg) L_g else L_h,
    y_trans = if(isg) y_j2a else y_a2b,
    y_base = if(isg) y_j else y_a,
    rep_species = if(isg) species_g else species_h,
    N_species = N_species,
    area_log = if(isg) area_log_j2a else area_log_a2b,
    
    prior_rate = if(isg) prior_rate["g"] else prior_rate["h"]
    
  ))
  
  n_chains <- 4
  fit_transition <- model_transitions$sample(data = d,
                                             output_dir = fitpath,
                                             init = 0.1,
                                             # init = lapply(1:n_chains, function(x) list(phi_inv = c(1, 1), rate_loc = c(-1, -1))),
                                             # iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                                             adapt_delta = 0.95, ## difficult geometry
                                             chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
  
  #### Pairs
  ggsave(paste0(fitpath, "/Pairs_transitions_", which, ".png"),
         bayesplot::mcmc_pairs(fit_transition$draws(variables = c("rate_log"))), # "phi", "rate_global", "rate_contrast", "sigma_raw"
         device = "png", width = 12, height = 12)
  
  # bayesplot::mcmc_trace(fit_transition$draws())
  # bayesplot::mcmc_areas(fit_transition$draws(variables = c("rate_log")), area_method = "scaled height")
  
  #### Summary
  s <- fit_transition$summary(variables = c("rate_log")) # "theta", "phi"
  write.csv(s, file.path(fitpath, paste0("summary_tansitions_", which, ".csv")))  
  message("Summary of the the fit for parameter ", which, ":")
  print(s)
  
  #### Residuals
  ## these do not include zero inflation
  y_sim <- fit_transition$draws(variables = "y_sim", format = "draws_matrix") %>% t()# matrix of observations simulated from the fitted model - row index for observations and colum index for simulations
  residuals <- DHARMa::createDHARMa(simulatedResponse = y_sim, observedResponse = d$y_trans, integerResponse = T) # fittedPredictedResponse = y_hat
  # testZeroInflation(residuals)
  # testDispersion(residuals)
  png(paste0(fitpath, "/Transitions_DHARMa_", which, ".png"), width = 1600, height = 1000)
  plot(residuals, quantreg = T, smoothScatter = F)
  dev.off()
  
  return(fit_transition)
}


## formatPriors --------------------------------
# data_stan  <- tar_read("data_stan")
# fit_g  <- tar_read("fit_g")
# fit_h  <- tar_read("fit_h")
# fit_Seedlings  <- tar_read("fit_Seedlings")
# af(formatPriors)
formatPriors <- function(data_stan, weakpriors,
                         fit_g = NULL, fit_h = NULL, fit_Seedlings = NULL,
                         widthfactor_trans = 1, widthfactor_reg = 1) {
  
  anyinferredprior <- NULL
  pars_g <- NULL
  pars_h <- NULL
  
  ## CASE: inferred priors are provided
  if( !is.null(fit_Seedlings)) {
    seedlingpar <- c("r_log") # , "k_log", "l_log"
    draws_seedlings <- sapply(seedlingpar,
                              function(v) fit_Seedlings$draws(variables = v, format = "draws_matrix") %>% as.data.frame(),
                              simplify = F, USE.NAMES = T)
    pars_r <- lapply(draws_seedlings[["r_log"]], function(d) MASS::fitdistr(d, "normal")$estimate)
    # pars_l <- lapply(draws_seedlings[["l_log"]], function(d) MASS::fitdistr(d, "normal")$estimate)
    # pars_k <- lapply(draws_seedlings[["k_log"]], function(d) MASS::fitdistr(d, "normal")$estimate)
    ## pars_seedlings <- lapply(Draws_seedlings, function(d) MASS::fitdistr(d, "normal")$estimate)
    
    anyinferredprior <- T
  }
  
  if( !is.null(fit_g) & !is.null(fit_h) ) {
    
    ## Matrices[draws, species]
    Draws_g <- fit_g$draws(variables = "rate_log", format = "draws_matrix") %>% as.data.frame()
    Draws_h <- fit_h$draws(variables = "rate_log", format = "draws_matrix") %>% as.data.frame()
    
    
    
    # Draws_seedlings <- posterior_samples(fit_Seedlings, fixed = F, pars = c("b_ba_ha$", "b_s_"))
    
    pars_g <- lapply(Draws_g, function(d) MASS::fitdistr(d, "normal")$estimate)
    pars_h <- lapply(Draws_h, function(d) MASS::fitdistr(d, "normal")$estimate)
    
    anyinferredprior <- T
  }
  
  if(isTRUE(anyinferredprior)) {
    
    if(widthfactor_trans != 1) {
      pars_g <- lapply(pars_g, function(p) c(p["mean"], widthfactor_trans*p["sd"]))
      pars_h <- lapply(pars_h, function(p) c(p["mean"], widthfactor_trans*p["sd"]))
    }
    
    if(widthfactor_reg != 1) {
      pars_r <- lapply(pars_r, function(p) c(p["mean"], widthfactor_reg*p["sd"]))
      # pars_l <- lapply(pars_l, function(p) c(p["mean"], widthfactor_reg*p["sd"]))
      # pars_k <- lapply(pars_k, function(p) c(p["mean"], widthfactor_reg*p["sd"]))
      ## pars_seedlings <- lapply(pars_seedlings, function(p) c(p["mean"], widthfactor_reg*p["sd"]))
    }
    
    ## The model assumes array[2] vector[N_species] prior_*; which means that the vectors stretch over rows!
    
    priors <- list(
      prior_g_log = bind_cols(pars_g), ## Matrix[N_species, (mu, sigma)]
      prior_h_log = bind_cols(pars_h),
      prior_r_log = bind_cols(pars_r)
      # prior_l_log = bind_cols(pars_l),
      # prior_k_log = bind_cols(pars_k)
    )
    
    data_stan_priors <- c(data_stan, priors)
  
  } else {
    
    data_stan_priors <- data_stan
  
  }
  
  ## ANY CASE
  data_stan_priors <- utils::modifyList(data_stan_priors, weakpriors)
  
  attr(data_stan_priors, "Long") <- attr(data_stan, "Long")
  attr(data_stan_priors, "Long_BA") <- attr(data_stan, "Long_BA")
  
  return(data_stan_priors)
}


## selectOffset --------------------------------
# data_stan_priors  <- tar_read("data_stan_priors")
# offsetname  <- tar_read("offsetname")
selectOffset <- function(offsetname, data_stan_priors) {
  
  L <- attr(data_stan_priors, "Long")
  
  Offset_0 <- L %>%
    filter(count_ha == 0) %>%
    group_by(obsid, tax, stage) %>%
    summarize_at(c("offset", "offset_avg", "offset_q1", "offset_q3"), .funs = function(x) mean(x, na.rm = T))
  
  write.csv(Offset_0, paste0("Publish.nosync/", Sys.Date(), "_", offsetname, "_", "Offset_0_averages.csv"))
  
  data_stan_priors_offset <- data_stan_priors
  data_stan_priors_offset$offset_data <- L[,offsetname[1], drop = T]
  
  attr(data_stan_priors_offset, "Long") <- L
  attr(data_stan_priors_offset, "offsetname") <- offsetname
  
  return(data_stan_priors_offset)
}


## fitModel --------------------------------
# data_stan <- tar_read("data_stan_priors_offset")
# model <- testmodel <- tar_read("model_test"); model <- tar_read("model_env")
# fitpath  <- tar_read("dir_fit")
fitModel <- function(model, data_stan, gpq = FALSE,
                     method = c("mcmc", "chkptstanr","variational", "sim", "diagnose"), n_chains = 4, iter_warmup = 1000, iter_sampling = 500, # openclid = c(0, 0),
                     fitseed = tar_seed(),
                     fitpath = dir_fit, ...) {
  
  if(!is.null(fitseed)) {
    ## Setting seed is necessary to reproduce the same targets when the targets pipeline has changed in another place, e.g. after tidying up for publication
    ## this function can be called within a target with fitseed = tar_seed() to retrieve the seed that is generated based on the name
    set.seed(fitseed)
  }
  
  require(cmdstanr)
  
  data_stan$generateposteriorq <- as.integer(gpq)
  
  inits <- 1e-1
  # inits <- list(state_init_raw = apply(data_stan$state_init_data, 2, function(x) scales::rescale(x, to = c(1e-12, 0.7))),
  #               b_log = data_stan$prior_b_log[1], c_a_log = data_stan$prior_c_a_log[1], c_b_log = data_stan$prior_c_b_log[1], c_j_log = data_stan$prior_c_j_log[1],
  #               g_log = unlist(data_stan$prior_g_log[1,], use.names = F), h_log = unlist(data_stan$prior_h_log[1,], use.names = F), l_log = data_stan$prior_l_log[1], r_log = unlist(data_stan$prior_r_log[1,], use.names = F), s_log = data_stan$prior_s_log[1], 
  #               phi_obs_inv_sqrt = 10
  #               )
  
  if (!dir.exists(fitpath)) {
    dir.create(fitpath)
  }
  
  if(match.arg(method) == "variational") {
    
    # inits <- replicate(1, inits, simplify = F)
    
    fit <- model$variational(data = data_stan,
                             output_dir = fitpath,
                             init = inits,
                             eta = 0.001,
                             iter = 20**4)
    
  } else if (match.arg(method) == "mcmc") {
    
    ## https://mc-stan.org/cmdstanr/articles/opencl.html
    ## system("clinfo -l")
    
    # inits <- replicate(n_chains, inits, simplify = F)
    
    fit <- model$sample(data = data_stan,
                        output_dir = fitpath,
                        # output_basename = ,
                        init = inits,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        # opencl_ids = openclid,
                        # adapt_delta = 0.99,
                        max_treedepth = 12,
                        chains = n_chains, parallel_chains = getOption("mc.cores", n_chains),
                        ...)
    
  } else if (match.arg(method) == "chkptstanr") {
    
    # inits <- replicate(n_chains, inits, simplify = F)
    message("chkptstanr returns a draws object.")
    
    chkptpath <- create_folder(paste0(Sys.Date(), "_checkpoints"), path = fitpath) ## creates nested directory structure
    fit <- chkpt_stan(model_code = model$code(),
                      data = data_stan,
                      path = chkptpath,
                      init = inits,
                      iter_per_chkpt = 100, iter_typical = 200,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                      chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
    fit <- combine_chkpt_draws(object = fit)
    
  } else if (match.arg(method) == "sim") {
    
    # inits <- replicate(1, inits, simplify = F)
    
    fit <- model$sample(data = data_stan,
                        fixed_param = TRUE,
                        output_dir = fitpath,
                        # output_basename = ,
                        init = inits, iter_sampling = iter_sampling)
  }
  
  else if (match.arg(method) == "diagnose") {
    
    fit <- model$diagnose(data = data_stan,
                          output_dir = fitpath,
                          # output_basename = ,
                          init = inits)
  }
  
  basename <- fit$output_files()[1] %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    str_replace("-[1-9]-", "-x-")
  
  ### Write out fit data
  
  ## Write a txt file with the seed
  if(!is.null(fitseed)) {
    as.character(fitseed) %>% writeLines(file(file.path(fitpath, paste0(basename, "_seed", ".txt"))))
  }
  
  ## Write an empty file to indicate the used offset
  if ( !is.null(attr(data_stan, "offsetname")) ) {
    file.create(file.path(fitpath, paste0(basename, "_", attr(data_stan, "offsetname"), ".txt")), showWarnings = TRUE)
  }
  
  ## Write code
  model$code() %>% writeLines(file(file.path(fitpath, paste0(basename, "_code", ".stan"))))
  
  ## Write median for L_p
  l_other_unique <- unique(data_stan$L_smooth[,"s_other"])
  if( isTRUE(length(l_other_unique) == 1) ) {
    l_other_unique %>% as.character %>% writeLines(file(file.path(fitpath, paste0(basename, "_median_l_other", ".txt"))))
  }
  
  ## Write data
  saveRDS(data_stan, file.path(fitpath, paste0(basename, "_data", ".rds")))
  
  attr(fit, "basename") <- basename
  
  return(fit)
}


## getBaseName --------------------------------
# cmdstanfit <- tar_read("fit_test")
getBaseName <- function(cmdstanfit) {
  
  basename <- cmdstanfit$output_files()[1] %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    str_replace("-[1-9]-", "-x-")
  
  if(is.na(basename)) basename <- "Model_failed"
  
  return(basename)
}

