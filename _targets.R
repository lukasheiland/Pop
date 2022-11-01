# ————————————————————————————————————————————————————————————————————————————————— #
# Targets setup ------------------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

### Installation
# remotes::install_github("ropensci/targets")
# remotes::install_github("ropensci/tarchetypes")
# remotes::install_github("ropensci/stantargets")
# cmdstanr::install_cmdstan()

### Library
library(targets)
library(tarchetypes)
# library(stantargets)

### Future
library(future)
library(future.callr)
plan(callr) ## "It is crucial that future::plan() is called in the target script file itself"

### Source the functions
source("Wrangling_functions.R")
source("Map_functions.R")
source("Fit_seedlings_functions.R")
source("Fit_functions.R")
source("Posterior_functions.R")
source("Theme/Theme.R")


### Options
options(tidyverse.quiet = TRUE)
tar_option_set(error = "abridge")

package <- c("dplyr", "ggplot2", "tidyr", "magrittr", "glue", "forcats", "vctrs", "tibble", "stringr", "knitr", # "multidplyr" ## extended tidyverse
             "lubridate", "DescTools",
             "sf", "raster", "rasterVis", ## for correct loading of environmental data
             "MASS", "mgcv", "glmnet", "itsadug", "stargazer",
             "cmdstanr", "rstan", "brms", "posterior", "bayesplot", "tidybayes", "parallel", "DHARMa", "priorsense", # "chkptstanr",
             "cowplot", "hrbrthemes", "showtext", "ggallin", "ggridges", "elementalist",  "ggspatial", "GGally", "scales", "gganimate", "metR",
             "future.apply")
tar_option_set(packages = package)

### Operation system
onserver <- Sys.info()["sysname"] != "Darwin"
if (!onserver) tar_config_set(store = "_targets.nosync")


# ————————————————————————————————————————————————————————————————————————————————— #
# Inner pipelines -----------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## Settings pipeline ------------------------------------------------------
targets_settings <- list(
  
  ## Whether to fit the model with a population ("loc", i.e. location) structure ...
  ##  "plot" — where the populations correspond to a plot in the data (n_locs == n_plots); environment joined by plot; counts/offsets correspond to a plot
  ##  "nested" — where the populations correspond to a cluster, but get fitted to data on the plot level (n_locs == n_clusters); environment joined by cluster; counts/offsets correspond to a plot
  ##  "cluster" — where the populations correspond to a cluster and get fitted to the sum within a cluster (n_locs == n_clusters); environment joined by cluster; counts/offsets correspond to the sum of the plots within a cluster
  tar_target(loc, c("plot", "nested", "cluster")[1]),
  
  ## No. of locations to subset (currently only for loc == "plot")
  tar_target(n_locations, 1000),
  
  ## Threshold to discriminate A and B [mm]
  ## 160 is the 10%tile, 185 is the 15%tile, 207 is the 20%tile, 228 is the 25%tile
  ## lower in the data is 100, so that: 100mm > A > threshold_dbh > B
  tar_target(threshold_dbh, 180), ## [mm]
  
  ## Upper sampling radius
  ## 	- All trees above a sampling radius of 14m were dropped, which is about the 98%tile (14.08m). The radius of 14m corresponds to the threshold radius of trees with dbh = 56cm
  ##    - dbh_threshold = radius_threshold/c with c == 25
  ##    - Alternatives: 99% radius == 1571 cm, 95% radius == 1188,  96% radius == 1242, 97% 1310.91
  tar_target(radius_max, 14000), ## [mm]
  
  ## Vector of taxa to select. All others will be lumped into "other".
  tar_target(taxon_select, c("Fagus.sylvatica")),
  
  ## Regeneration classes to include
  ## here we select all 20cm height <= trees < 7cm dbh: regglass_select <- c("h[20,50)" = 1, "h[50,130)" = 2, "hd[130,Inf)[0,5)" = 3, "d[5,6)" = 4, "d[6,7)" = 5)
  ## These are all size classes that are consistent across the three surveys.
  tar_target(regclass_select, c("h[20,50)" = 1, "h[50,130)" = 2, "hd[130,Inf)[0,5)" = 3, "d[5,6)" = 4, "d[6,7)" = 5, "hd[130,Inf)[0,7)" = 9)), ## [mm]; the last level (9) is a special category of counts in the second survey DE_BWI_2 over all classes
  
  ## Weakly informative priors.
  ## If provided, they are prioritzed over the fitted priors
  tar_target(weakpriors,
             ## Priors are organized like the parameter data structure but with an additional dimension in the case of a vector row of sds.
             list(
               prior_b_log = c(-3, 2),
               prior_c_a_log = c(-8, 3),
               prior_c_b_log = c(-7, 3),
               prior_c_j_log = c(-9, 3),
               # prior_g_log = cbind(Fagus = c(-6, 0.1), others = c(-6, 0.1)),
               # prior_h_log = cbind(Fagus = c(-4, 1), others = c(-4, 1)),
               # prior_l_log = cbind(Fagus = c(4, 1), others = c(5, 1)),
               prior_l_log = c(4, 1),
               # prior_r_log = cbind(Fagus = c(4, 1), others = c(4, 1)),
               prior_s_log = c(-6, 2)
             )
  ),
  
  tar_target(weakpriors_env,
             list(
               prior_b_log = c(-3.5, 1),
               prior_c_a_log = c(-7, 2),
               prior_c_b_log = c(-7, 2),
               prior_c_j_log = c(-12, 3),
               prior_g_log = c(-6, 2),
               prior_h_log = c(-3, 2),
               prior_l_log = c(4, 2),
               prior_r_log = c(4, 2),
               prior_s_log = c(-5, 2)
             )
  ),
  
  tar_target(twocolors, c("#119973", "#FFCC11")), ## Other greens include: Spanish Green '#208E50', Feldgrau '#3A7867' # Jade: #0FA564 # See green: #0E8C55
  tar_target(themefunction, theme_fagus)
  
)



## Paths pipeline ------------------------------------------------------
targets_paths <- list(
  
  ## Directories
  tar_target(dir_fit, file.path("Fits.nosync/")),
  tar_target(dir_publish, file.path("Publish.nosync/")),
  
  ## Data files
  tar_target(file_DE_big,
             "Inventory.nosync/DE BWI/Data/DE_BWI_big_abund.rds",
             format = "file"),
  tar_target(file_DE_big_status,
             "Inventory.nosync/DE BWI/Data/DE_BWI_big_2-3_status.rds",
             format = "file"),
  tar_target(file_DE_small,
             "Inventory.nosync/DE BWI/Data/DE_BWI_small_abund.rds",
             format = "file"),
  tar_target(file_DE_env,
             'Inventory.nosync/DE BWI/Data/DE_BWI_Env_sf.rds',
             format = "file"),
  tar_target(file_DE_geo,
             'Inventory.nosync/DE BWI/Data/DE_BWI_geo.rds',
             format = "file"),
  tar_target(file_SK,
             "Inventory.nosync/SK NIML/Data/SK_NIML_complete.rds",
             format = "file"),
  tar_target(file_SK_fullgrid,
             "Inventory.nosync/SK NIML/Data/SK_NIML_complete_fullgrid.rds",
             format = "file"),
  # tar_target(file_Taxa,
  #            'Inventory.nosync/Taxa/Taxa.csv',
  #            format = "file")
  
  ## References
  # dput(.packages()) ## does not work directly because of incomplete metadata in some
  tar_target(bibliography,
             knitr::write_bib(c("scales", "ggspatial", "DHARMa",
                                "tidybayes", "bayesplot", "posterior", "brms", "rstan", "cmdstanr", 
                                "MASS", "mgcv", "nlme", "rasterVis", "lattice", "raster", "sp", 
                                "sf", "DescTools", "lubridate", "stringr", "tibble", "vctrs"), file.path(dir_publish, "Packages.bib"))
             )
  
)



## Parameter spec pipeline ------------------------------------------------------
targets_parname <- list(
  
  tar_target(par_exclude,
             c("y_hat", "L_loc_log", "K_loc_log_raw", "L_loc", "K_loc", "state_init", "state_init_raw", "state_init_log", "phi_obs_inv", "phi_obs_inv_sqrt", "m")),
  
  tar_target(helper_exclude,
             c("Fix", "Fix_ko_s",
               paste0("Fix_switch_", c(names(parname_plotorder), "b_c_b", "b_c_a_c_b_h", "l_r", "g_l_r_s")),
               "B_log_raw", "C_a_log_raw", "C_b_log_raw", "C_j_log_raw", "G_log_raw", "H_log_raw", "L_log_raw", "R_log_raw", "S_log_raw",
               "vector_b_log_prior", "vector_c_a_log_prior", "vector_c_b_log_prior", "vector_c_j_log_prior", "vector_s_log_prior",
               "phi_obs_rep", "phi_obs_rep_prior", "theta_rep",
               "avg_state_init", "avg_L_loc",
               "fixiter_max", "fixiter_min", "eps_ba_fix",
               "c_j_log_10",
               "log_prior", "log_lik", "lp__", "state_init_log_raw", "state_2", "state_2_raw", "state_3", "state_3_raw")),
  
  ## do not exclude: sort(apply(expand.grid(paste0(setdiff(parname_plotorder, "l_log"), "_spread_env"), 1:2, "_100"), 1, paste0, collapse = ""))

  tar_target(rep_exclude,
             c("phi_obs_rep", "phi_obs_rep_prior", "theta_rep",
               "y_hat_rep", "y_hat_prior_rep", "y_hat_offset", "y_hat_rep_offset", "y_hat_prior_rep_offset")),
  
  tar_target(simname_prior,
             c("y_hat_prior", "y_hat_prior_rep", "y_hat_prior_offset", "y_hat_prior_rep_offset", "y_prior_sim")),
  
  tar_target(statename,
             c("ba_init", "ba_fix", "J_init", "J_fix", "A_init", "A_fix", "B_init", "B_fix",
               "ba_tot_init", "ba_tot_fix", "ba_frac_init", "ba_frac_fix",
               "ba_fix_ko_b", "ba_fix_ko_s", "ba_fix_ko_2_b", "ba_fix_ko_2_s",
               "major_init", "major_fix", "dominant_init", "dominant_fix",
               "major_fix_ko_b", "major_fix_ko_s", "major_fix_ko_2_b", "major_fix_ko_2_s",
               paste0("major_fix_switch_", c(names(parname_plotorder), "b_c_b", "b_c_a_c_b_h", "l_r", "g_l_r_s")),
               "ba_fix_ko_b", "ba_fix_ko_s", "ba_fix_ko_2_b", "ba_fix_ko_2_s",
               paste0("ba_fix_switch_", c(names(parname_plotorder), "b_c_b", "b_c_a_c_b_h", "l_r", "g_l_r_s")))
             ),
  
  tar_target(basalareaname, statename[str_starts(statename, "ba_")]),
  
  tar_target(majorname, statename[str_starts(statename, "major_")]),

  tar_target(simname_posterior,
             c(statename,
               paste0("Fix_switch_", c(names(parname_plotorder), "b_c_b", "b_c_a_c_b_h", "l_r", "g_l_r_s")),
               
               "converged_fix", "iterations_fix", "fixiter_max", "fixiter_min", "eps_ba_fix",
               "sum_ko_1_b_fix", "sum_ko_1_b_c_b_fix", "sum_ko_1_c_a_fix", "sum_ko_1_c_b_fix", "sum_ko_1_c_j_fix", "sum_ko_1_g_fix", "sum_ko_1_h_fix", "sum_ko_1_l_fix", "sum_ko_1_r_fix", "sum_ko_1_s_fix",
               "sum_ko_2_b_fix", "sum_ko_2_b_c_b_fix", "sum_ko_2_c_a_fix", "sum_ko_2_c_b_fix", "sum_ko_2_c_j_fix", "sum_ko_2_g_fix", "sum_ko_2_h_fix", "sum_ko_2_l_fix", "sum_ko_2_r_fix", "sum_ko_2_s_fix",
               
               "sum_switch_b_fix", "sum_switch_b_c_b_fix", "sum_switch_c_a_fix", "sum_switch_c_b_fix", "sum_switch_c_j_fix", "sum_switch_g_fix", "sum_switch_h_fix", "sum_switch_l_fix", "sum_switch_r_fix", "sum_switch_s_fix",
               
               "sum_ko_1_prop_b_fix", "sum_ko_1_prop_b_c_b_fix", "sum_ko_1_prop_c_a_fix", "sum_ko_1_prop_c_b_fix", "sum_ko_1_prop_c_j_fix", "sum_ko_1_prop_g_fix", "sum_ko_1_prop_h_fix", "sum_ko_1_prop_l_fix", "sum_ko_1_prop_r_fix", "sum_ko_1_prop_s_fix",
               "sum_ko_2_prop_b_fix", "sum_ko_2_prop_b_c_b_fix", "sum_ko_2_prop_c_a_fix", "sum_ko_2_prop_c_b_fix", "sum_ko_2_prop_c_j_fix", "sum_ko_2_prop_g_fix", "sum_ko_2_prop_h_fix", "sum_ko_2_prop_l_fix", "sum_ko_2_prop_r_fix", "sum_ko_2_prop_s_fix",
               
               "greater_b", "greater_c_a", "greater_c_b", "greater_c_j", "greater_g", "greater_h", "greater_l", "greater_k", "greater_r", "greater_s")
             ),
  
  tar_target(parname,
             c("phi_obs", # "sigma_k_loc", # "k_log", # "theta", 
               "b_log", "c_a_log", "c_b_log", "c_j_log", "g_log", "h_log", "l_log", "r_log", "s_log")),
  
  tar_target(parname_plotorder,
             c(l = "l_log", r = "r_log", c_j = "c_j_log", s = "s_log", g = "g_log", c_a = "c_a_log", h = "h_log", b = "b_log", c_b = "c_b_log" )),
  
  tar_target(parname_loc,
             c("state_init", "L_loc")),
  
  tar_target(parname_sigma,
             c(l = "sigma_l", r = "sigma_r", c_j = "sigma_c_j", s = "sigma_s", g = "sigma_g", c_a = "sigma_c_a", h = "sigma_h", b = "sigma_b", c_b = "sigma_c_b" )),
  
  tar_target(parname_alpha,
             c(l = "alpha_l", r = "alpha_r", c_j = "alpha_c_j", s = "alpha_s", g = "alpha_g", c_a = "alpha_c_a", h = "alpha_h", b = "alpha_b", c_b = "alpha_c_b" )),
  
  tar_target(parname_sim,
             setdiff(parname, c("theta", "phi_obs", "sigma_k_loc"))),
  
  tar_target(exclude,
             c(par_exclude, parname_loc_env, parname_sigma, parname_alpha, helper_exclude, rep_exclude, simname_prior, simname_posterior))
)



## Wrangling pipeline ------------------------------------------------------
targets_wrangling <- list(

  ## Read data files
  # tar_load(starts_with("Data"))
  list(
    tar_target(Data_big, readRDS(file_DE_big)),
    tar_target(Data_big_status, readRDS(file_DE_big_status)),
    tar_target(Data_small, readRDS(file_DE_small)),
    tar_target(Data_env, readRDS(file_DE_env)),
    tar_target(Data_geo, readRDS(file_DE_geo)),
    
    tar_target(Data_seedlings, readRDS(file_SK)),
    tar_target(Data_seedlings_fullgrid, readRDS(file_SK_fullgrid))
    
    # tar_target(Taxa, read.csv(file = file_Taxa, colClasses = c('factor')) %>% filter(!duplicated(tax.id)))
    ## tax.id is not unique in Taxa! Unique is however needed for left_join by tax.id (not by inventory specific ids)!
  ),
  
  
  ## Env data
  # tar_load(starts_with("Env"))
  list(
    tar_target(predictor_select,
               c("phCaCl_esdacc", "waterLevel_loc")), # "alt_loc"
    tar_target(Env_clean,
               cleanEnv(Data_env, predictor_select)),
    ## Summarize (mutate) predictors per cluster, but keep plot-level disturbance data etc.
    tar_target(Env_cluster,
               summarizeEnvByCluster(Env_clean, predictor_select))
  ),
  

  ## Stage abundance data
  list(
    tar_target(Data_big_area,
               prepareBigData(Data_big, Data_big_status,
                              taxon_select = taxon_select, threshold_dbh = threshold_dbh, radius_max = radius_max,
                              tablepath = dir_publish)),
    
    tar_target(Data_small_area,
               prepareSmallData(Data_small, taxon_select = taxon_select, regclass_select = regclass_select)),
    
    tar_target(Stages_transitions,
               countTransitions(Data_big, Data_big_status, Env_cluster, Stages_select,
                                taxon_select = taxon_select, threshold_dbh = threshold_dbh, radius_max = radius_max,
                                loc = c("plot", "nested", "cluster"))),
    
    tar_target(Stages,
               joinStages(Data_big_area, Data_small_area, taxon_select = taxon_select, threshold_dbh = threshold_dbh)),
    tar_target(Stages_env,
               joinEnv(Stages, Env_clean)), # Env_cluster ## always joins by plotid, even if Env is aggregated by cluster
    
    list(
      tar_target(taxon_s,
                 factor(c(sort(taxon_select), "other"))),
      tar_target(BA_s,
                 constructConstantGrid(taxon_s, Stages_env, Data_geo),
                 pattern = map(taxon_s),
                 iteration = "list"),
      tar_target(fits_s,
                 fitS(BA_s, path = dir_publish),
                 pattern = map(BA_s),
                 iteration = "list"),
      ## fits_s each have an attribute "taxon"
      tar_target(Stages_s,
                 predictS(fits_s, Stages_env),
                 iteration = "list"),
        ## Former workaround for Linux installlations, where the geo libraries did not work. Assumes upload of "Data/Stages_s.rds"
        # tar_target(file_Stages_s, if (!! onserver) "Data/Stages_s.rds" else saveStages_s(Stages_s), format = "file"),
        # tar_target(Data_Stages_s, readRDS(file_Stages_s)), ## this would be loaded in downstream targets instead of Stages_s
      
      tar_target(surfaces_s,
                 predictSurfaces(fits_s),
                 iteration = "list"),
      tar_target(surfaceplots_s,
                 plotSurfaces(surfaces_s, themefun = themefunction, path = dir_publish),
                 iteration = "list")
    ),
    
    ## The seedlings pipeline
    list(
      tar_target(Seedlings,
                 wrangleSeedlings(Data_seedlings, taxon_select = taxon_select, threshold_dbh = threshold_dbh)),
      tar_target(seedlings_fullgrid,
                 wrangleSeedlings_s(Data_seedlings_fullgrid, taxon_select = taxon_select, threshold_dbh = threshold_dbh),
                 iteration = "list"),
      tar_target(fits_Seedlings_s, ## fits_s each have an attribute "taxon"
                 fitS(seedlings_fullgrid, path = NULL), ##! setting path to NULL prevents overwriting of the summary of the gam
                 pattern = map(seedlings_fullgrid),
                 iteration = "list"),
      tar_target(surfaces_Seedlings_s,
                 predictSeedlingsSurfaces(fits_Seedlings_s),
                 iteration = "list"),
      tar_target(surfaceplots_Seedlings_s,
                 plotSurfaces(surfaces_Seedlings_s, path = dir_publish),
                 iteration = "list"),
      tar_target(Seedlings_s,
                 predictS(fits_Seedlings_s, Seedlings)),
      tar_target(file_model_Seedlings,
                 "Model_2021-03_ba/Model_seedlings.stan",
                 format = "file"),
      tar_target(model_Seedlings,
                 cmdstan_model(file_model_Seedlings, stanc_options = list("O1"))),
      tar_target(fit_Seedlings,
                 fitSeedlings(model_Seedlings, Seedlings_s, fitpath = dir_fit))
    ),
    
    
    tar_target(Stages_select,
               selectLocs(Stages_s, predictor_select,
                          selectspec = T, selectpred = F, stratpred = F, n_locations = n_locations, loc = "plot", tablepath = dir_publish)), # Subsetting after smooth, so that smooth can be informed by all plots.
    
    tar_target(Stages_scaled,
               scaleData(Stages_select, predictor_select)), # After selection, so that scaling includes selected plots .
    
    tar_target(Stages_loc,
               setLocLevel(Stages_scaled, Env_cluster, loc = c("plot", "nested", "cluster"))),

    
    ## Publishing
    tar_target(Summary_taxa,
               summarizeTaxa(Data_big, Data_seedlings, Stages_select, Seedlings_s, tablepath = dir_publish)),
    tar_target(Summary_NFIs,
               summarizeNFIs(Data_big, Data_seedlings, Stages_select, Seedlings_s, tablepath = dir_publish)),
    tar_target(map_select,
               mapLocations(Stages_select, path = dir_publish, themefun = themefunction),
               packages = c(package, "eurostat", "elevatr", "terrainr", "rayshader", "ggspatial", "elementalist"))
    
  )
)



## Fitting pipelines ------------------------------------------------
#### general ----------
targets_fit_general <- list(
  tar_target(data_stan,
             formatStanData(Stages_loc, Stages_transitions, taxon_s, threshold_dbh, predictor_select,
                            loc = c("plot", "nested", "cluster"))),
  
  tar_target(data_stan_transitions,
             formatStanData(Stages_loc, Stages_transitions, taxon_s, threshold_dbh, predictor_select,
                            loc = "plot")),
  
  tar_target(file_model_transitions,
             "Model_2021-03_ba/Model_transitions.stan",
             format = "file"),
  
  tar_target(model_transitions,
             cmdstan_model(file_model_transitions, stanc_options = list("O1"))),
  
  tar_target(fit_g,
             fitTransition(data_stan_transitions, which = "g", model_transitions, fitpath = dir_fit)),
  
  tar_target(fit_h,
             fitTransition(data_stan_transitions, which = "h", model_transitions, fitpath = dir_fit)),
  
  tar_target(data_stan_priors,
             formatPriors(data_stan, weakpriors, fit_g, fit_h, fit_Seedlings,
                          widthfactor_trans = 3, widthfactor_reg = 4)),
  
  tar_target(offsetname,
             c("offset", "offset_avg", "offset_q1", "offset_q3")),
  tar_target(offsetname_select,
             offsetname[1]),
  tar_target(data_stan_priors_offset,
             selectOffset(offsetname_select, data_stan_priors)),
  tar_target(data_stan_priors_offsets,
             selectOffset(offsetname, data_stan_priors),
             pattern = map(offsetname), iteration = "list"),
  tar_target(data_stan_priors_offsets_1,
             data_stan_priors_offsets[1], iteration = "list")
)

#### fit_test ----------
targets_fit_test <- list(

  tar_target(file_model_test,
             "Model_2021-03_ba/Model_ba_test.stan",
             format = "file"),
  tar_target(model_test,
             cmdstan_model(file_model_test, stanc_options = list("O1")) #, cpp_options = list(stan_opencl = TRUE)
             ),
  ## Prior predictive tests that rely on currently out-commented generated quantities
  # tar_target(priorsim_test,
  #            drawTest(model = model_test, data_stan = data_stan_priors_offset, method = "sim", gpq = FALSE, fitpath = dir_fit)),
  # tar_target(plots_predictions_priorsim_test,
  #            plotPredictions(cmdstanfit = priorsim_test, data_stan_priors_offset, check = "prior")),

  tar_target(fit_test_sansgq,
             fitModel(model = model_test, data_stan = data_stan_priors_offset, gpq = FALSE,
                      method = "mcmc", n_chains = 4, iter_warmup = 800, iter_sampling = 500, fitpath = dir_fit)),
  tar_target(fit_test,
             fitModel(model = model_test, data_stan = data_stan_priors_offset, gpq = TRUE,
                      method = "mcmc", n_chains = 4, iter_warmup = 800, iter_sampling = 500, fitpath = dir_fit)),
  tar_target(basename_fit_test,
             getBaseName(fit_test))
)

#### fit ----------
targets_fit <- list(
  tar_target(file_model,
             "Model_2021-03_ba/Model_ba.stan",
             format = "file"),
  tar_target(model,
             cmdstan_model(file_model, stanc_options = list("O1"))),
  tar_target(fit,
             fitModel(model = model, data_stan = data_stan_priors_offsets_1, gpq = TRUE, ## data_stan_priors_offsets ## for testing all offsets
                      method = "mcmc", n_chains = 4, iter_warmup = 1000, iter_sampling = 1000, fitpath = dir_fit),
             pattern = map(data_stan_priors_offsets_1), iteration = "list"), ## data_stan_priors_offsets ## for testing all offsets
  tar_target(basename_fit,
             getBaseName(fit),
             pattern = map(fit), iteration = "list")
)

#### fit_env ----------
targets_fit_env <- list(
  
  ## Prepare data
  tar_target(Stages_select_env,
             selectLocs(Stages_s, predictor_select,
                        selectspec = F, selectpred = T, stratpred = T, n_locations = n_locations, loc = "plot", tablepath = dir_publish)), # Selection based on whether environmental variables are present
  
  tar_target(Stages_scaled_env,
             scaleData(Stages_select_env, predictor_select)), # After selection, so that scaling includes selected plots 
  
  tar_target(Stages_loc_env,
             setLocLevel(Stages_scaled_env, Env_cluster, loc = c("plot", "nested", "cluster"))),
  
  tar_target(Stages_transitions_env,
             countTransitions(Data_big, Data_big_status, Env_cluster, Stages_select_env,
                              taxon_select = taxon_select, threshold_dbh = threshold_dbh, radius_max = radius_max,
                              loc = c("plot", "nested", "cluster"))),
  
  tar_target(data_stan_env,
             formatStanData(Stages_loc_env, Stages_transitions_env, taxon_s, threshold_dbh,  predictor_select,
                            loc = c("plot", "nested", "cluster"))),
  
  tar_target(data_stan_transitions_env,
             formatStanData(Stages_loc, Stages_transitions_env, taxon_s, threshold_dbh, predictor_select,
                            loc = "plot")),
  
  tar_target(fit_g_env,
             fitTransition(data_stan_transitions_env, which = "g", model_transitions, fitpath = dir_fit)),
  
  tar_target(fit_h_env,
             fitTransition(data_stan_transitions_env, which = "h", model_transitions, fitpath = dir_fit)),
  
  tar_target(data_stan_priors_env,
             formatPriors(data_stan_env, weakpriors_env)),

  tar_target(data_stan_priors_offset_env,
             selectOffset(offsetname_select, data_stan_priors_env)),

  ## Fit
  tar_target(file_model_env,
             "Model_2021-03_ba/Model_ba_env.stan",
             format = "file"),
  tar_target(file_model_env_vertex,
             "Model_2021-03_ba/Model_ba_env_vertex.stan",
             format = "file"),
  tar_target(model_env,
             cmdstan_model(file_model_env_vertex, stanc_options = list("O1"))),
  tar_target(fit_env,
             fitModel(model = model_env, data_stan = data_stan_priors_offset_env, gpq = TRUE,
                      method = "mcmc", n_chains = 6, iter_warmup = 800, iter_sampling = 600, fitpath = dir_fit,
                      adapt_delta = 0.95)
             ),
  tar_target(basename_fit_env,
             getBaseName(fit_env))
)



## Posterior pipeline ------------------------------------------------------
#### posterior_test -----------
targets_posterior_test <- list(
  
  ## Extract
  tar_target(stanfit_test,
             extractStanfit(cmdstanfit = fit_test)),
  tar_target(draws_test,
             extractDraws(stanfit = stanfit_test, exclude = helper_exclude)),
  # tar_target(stanfit_test_plotting,
  #            extractStanfit(fit_test_sansgq, purge = TRUE)),
  
  
  ## Summarize
  tar_target(summary_test,
             summarizeFit(cmdstanfit = fit_test, exclude = c(helper_exclude, rep_exclude, par_exclude, simname_prior, parname_loc),
                          publishpar = parname_plotorder, path = dir_publish)),
  tar_target(summary_states_test,
             summarizeStates(States = States_test, data_stan = data_stan, basename = basename_fit_test, path = dir_publish)),
  tar_target(Freq_converged_test,
             summarizeFreqConverged(cmdstanfit = fit_test, data_stan_priors, path = dir_publish)),
  
  
  ## Generate
  tar_target(residuals_test,
             generateResiduals(cmdstanfit = fit_test, data_stan_priors, path = dir_publish)),
  # tar_target(Trajectories_test,
  #            generateTrajectories(cmdstanfit = fit_test, data_stan_priors, parname, locparname = parname_loc,
  #                                 time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 50, average = "none")),
  tar_target(Trajectories_avg_test,
             generateTrajectories(cmdstanfit = fit_test, data_stan_priors, parname, locparname = parname_loc,
                                  time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1, average = "locsperdraws_all")),
  # tar_target(Trajectories_quantiles_test,
  #            generateTrajectories(cmdstanfit = fit_test, data_stan_priors, parname, locparname = parname_loc,
  #                                 time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1, average = "locsperdraws_avgL_qInit")),
            ### average options:
            ## c("none",            # — no averaging. Paraneters in simulations vary per loc and draw
            ## "locsperdraws_all",  # — average all loc-wise parameters per draw
            ## "drawsperlocs_all",  # — average all parameters per loc, so that there is only one trajectory per cluster
            ## "locsperdraws_avgL", # — average only "L_loc" per draw, so that there are initial values that vary with cluster
            ## "locsperdraws_avgL_qInit", # — average only "L_loc" per draw, so that there are initial values that vary with cluster
  
  ## Formatted posterior data stuctures
  tar_target(States_test,
             formatStates(cmdstanfit = fit_test, statename = statename, data_stan_priors = data_stan_priors)),
  
  ## Plot
  tar_target(plots_test,
             plotStanfit(stanfit = stanfit_test, exclude = exclude, path = dir_publish, basename = basename_fit_test, color = twocolors, themefun = themefunction)),
  tar_target(plots_parameters_test,
             plotParameters(stanfit = stanfit_test, parname = parname_plotorder, exclude = exclude, path = dir_publish, basename = basename_fit_test, color = twocolors, themefun = themefunction)),
  ## Prior predictive tests that rely on currently out-commented generated quantities
  # tar_target(plots_predictions_prior_test,
  #            plotPredictions(cmdstanfit = fit_test, data_stan_priors, check = "prior")),
  tar_target(plots_predictions_posterior_test,
             plotPredictions(cmdstanfit = fit_test, data_stan_priors_offset, check = "posterior", path = dir_publish)),
  tar_target(plots_trace_test,
             plotTrace(cmdstanfit = fit_test, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plots_pairs_test,
             plotPairs(cmdstanfit = fit_test, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plots_conditional_test,
             plotConditional(cmdstanfit = fit_test, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plot_contributions_test,
             plotContributions(cmdstanfit = fit_test, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, color = twocolors, themefun = themefunction)),
  # tar_target(plot_contributions_prop_test,
  #            plotContributions(cmdstanfit = fit_test, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, contribution = "sum_ko_prop", color = twocolors, themefun = themefunction)),
  tar_target(plots_states_test,
             plotStates(States_test, allstatevars = basalareaname,
                        path = dir_publish, basename = basename_fit_test, color = twocolors, themefun = themefunction)),
  tar_target(plot_trajectories_avg_test,
             plotTrajectories(Trajectories_avg_test, thicker = T, path = dir_publish, basename = basename_fit_test, color = twocolors, themefun = themefunction)),
  tar_target(animation_trajectories_avg_test,
             animateTrajectories(plot_trajectories_avg_test, path = dir_publish, basename = basename_fit_test)),
  
  # tar_target(plot_powerscale_test,
  #            plotSensitivity(cmdstanfit = fit_test, include = parname, path = dir_publish)),
  
  ## Test
  tar_target(sensitivity_test,
             testSensitivity(fit_test, include = parname, path = dir_publish))
)

#### posterior -----------
targets_posterior <- list(
  ## Extract
  tar_target(stanfit,
             extractStanfit(cmdstanfit = fit),
             pattern = map(fit), iteration = "list"),
  tar_target(draws,
             extractDraws(stanfit = stanfit, exclude = helper_exclude),
             pattern = map(stanfit), iteration = "list"),
  
  ## Summarize
  tar_target(summary,
             summarizeFit(cmdstanfit = fit, exclude = c(helper_exclude, rep_exclude, par_exclude, simname_prior, parname_loc),
                          publishpar = parname_plotorder, path = dir_publish),
             pattern = map(fit), iteration = "list"),
  tar_target(summary_states,
             summarizeStates(States = States, data_stan = data_stan, basename = basename_fit, path = dir_publish),
             pattern = map(States, basename_fit), iteration = "list"),
  tar_target(Freq_converged,
             summarizeFreqConverged(cmdstanfit = fit, data_stan_priors, path = dir_publish),
             pattern = map(fit), iteration = "list"),
  
  ## Generate
  tar_target(residuals,
             generateResiduals(cmdstanfit = fit, data_stan_priors, path = dir_publish),
             pattern = map(fit), iteration = "list"),
  tar_target(Trajectories_avg,
             generateTrajectories(cmdstanfit = fit, data_stan_priors, parname, locparname = parname_loc,
                                  time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1, average = "locsperdraws_all"),
             pattern = map(fit), iteration = "list"),
  
  ## Formatted posterior data stuctures
  tar_target(States,
             formatStates(cmdstanfit = fit, statename = statename, data_stan_priors = data_stan_priors),
             pattern = map(fit), iteration = "list"),
  
  ## Plot
  tar_target(plots,
             plotStanfit(stanfit = stanfit, exclude = exclude, path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(stanfit, basename_fit), iteration = "list"),
  tar_target(plots_parameters,
             plotParameters(stanfit = stanfit, parname = parname_plotorder, exclude = exclude, path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(stanfit, basename_fit), iteration = "list"),
  tar_target(plots_trace,
             plotTrace(cmdstanfit = fit, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction),
             pattern = map(fit, basename_fit), iteration = "list"),
  tar_target(plots_pairs,
             # plotPairs(cmdstanfit = fit, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction),
             plotConditional(cmdstanfit = fit, parname = parname_plotorder, conditional = F, path = dir_publish, themefun = themefunction),
             pattern = map(fit, basename_fit), iteration = "list"),
  tar_target(plots_conditional,
             plotConditional(cmdstanfit = fit, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction),
             pattern = map(fit), iteration = "list"),
  tar_target(plot_contributions,
             plotContributions(cmdstanfit = fit, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, plotlog = T, color = twocolors, themefun = themefunction),
             pattern = map(fit), iteration = "list"),
  tar_target(plot_contributions_supp,
             plotContributions(cmdstanfit = fit, parname = c(parname_plotorder[1:7], b_c_b = "b_c_b_log"), path = dir_publish, plotlog = T, color = twocolors, themefun = themefunction),
             pattern = map(fit), iteration = "list"),
  tar_target(plot_contributions_prop,
             plotContributions(cmdstanfit = fit, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, contribution = "sum_ko_prop", plotlog = F, color = twocolors, themefun = themefunction),
             pattern = map(fit), iteration = "list"),
  tar_target(plot_contributions_switch,
             plotContributions(cmdstanfit = fit, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, contribution = "sum_switch", plotlog = T, color = twocolors, themefun = themefunction),
             pattern = map(fit), iteration = "list"),
  tar_target(plots_states,
             plotStates(States, allstatevars = basalareaname,
                        path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(States, basename_fit), iteration = "list"),
  tar_target(plot_predominant,
             plotPredominant(States, majorname,
                             path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(States, basename_fit), iteration = "list"),
  tar_target(plot_trajectories_avg,
             plotTrajectories(Trajectories_avg, thicker = T, path = dir_publish, basename = basename_fit, plotpdf = TRUE, color = twocolors, themefun = themefunction),
             pattern = map(Trajectories_avg, basename_fit), iteration = "list"),
  tar_target(animation_trajectories_avg,
             animateTrajectories(plot_trajectories_avg[[1]], path = dir_publish, basename = basename_fit[[1]]))
)

#### posterior_env -----------
targets_posterior_env <- list(
  
  ## Parameter names (in addition to the general ones above)
  tar_target(parname_env, c(setdiff(parname_loc_env, "state_init"),
                            "ba_init", "ba_fix", "ba_frac_init", "ba_frac_fix", "major_fix", "major_init")),
  tar_target(contribname_env, paste("sum_ko",
                                    rep(1:2, each = length(parname_plotorder)),
                                    names(parname_plotorder), "fix", sep = "_")),
  
  tar_target(parname_environmental, selectParnameEnvironmental(parname_env, Environmental_env)), ## selects the variables, that are actually in the fit
  tar_target(contribname_environmental, selectParnameEnvironmental(contribname_env, Environmental_env)), ## selects the variables, that are actually in the fit
  
  tar_target(parname_loc_env,
             c(parname_loc, "B_log", "C_a_log", "C_b_log", "C_j_log", "G_log", "H_log", "L_log", "R_log", "S_log")),
  
  tar_target(parname_env_vertex_center,
             sort(apply(expand.grid(paste0(setdiff(parname_plotorder, "l_log"), "_center_env"),
                                    1:2),
                        1, paste0, collapse = "")
                  )
             ),
  
  tar_target(parname_env_vertex_spread,
             sort(apply(expand.grid(paste0(setdiff(parname_plotorder, "l_log"), "_spread_env"),
                                      1:2),
                          1, paste0, collapse = "")
                  )
             ),
  
  tar_target(parname_env_vertex,
             sort(c(parname_env_vertex_center, parname_env_vertex_spread))),
  
  tar_target(parname_Beta,
             paste0("Beta_", setdiff(parname_plotorder, "l_log"))),
  
  tar_target(vertexname_env_vertex,
             sort(c(setdiff(parname_plotorder, "l_log"), parname_env_vertex_center, parname_env_vertex_spread))),
  
  
  
  ## Extract
  tar_target(stanfit_env,
             extractStanfit(cmdstanfit = fit_env)),
  tar_target(draws_env,
             extractDraws(stanfit = stanfit_env, exclude = helper_exclude)),
  
  ## Summarize
  tar_target(summary_env,
             summarizeFit(cmdstanfit = fit_env, exclude = c(helper_exclude, rep_exclude, par_exclude, simname_prior, parname_loc_env),
                          publishpar = c(parname_plotorder, parname_env_vertex), path = dir_publish)),
  tar_target(summary_states_env,
           summarizeStates(States = States_env, data_stan = data_stan_env, basename = basename_fit_env, path = dir_publish)),
  tar_target(Freq_converged_env,
           summarizeFreqConverged(cmdstanfit = fit_env, data_stan_priors_env, path = dir_publish)),

  ## Generate
  tar_target(residuals_env,
             generateResiduals(cmdstanfit = fit_env, data_stan_priors_env, path = dir_publish)),
  tar_target(Trajectories_avg_env,
             generateTrajectories(cmdstanfit = fit_env, data_stan_priors_env, parname, locparname = parname_loc_env,
                                  time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1, average = "drawsperlocs_all")),

  
  ## Formatted posterior data stuctures
  tar_target(States_env,
             formatStates(cmdstanfit = fit_env, statename = statename, data_stan_priors = data_stan_priors_env)),
  tar_target(Environmental_env,
             formatEnvironmental(cmdstanfit = fit_env, parname = c(parname_env, contribname_env),
                                 data_stan = data_stan_priors_offset_env, envname = predictor_select, locmeans = F)),
  tar_target(Envgrid_env, formatEnvgrid(data_stan_priors_offset_env, envname = predictor_select, res = 500)),
  tar_target(Marginal_env,
             formatMarginal(Environmental_env)),

  
  ## Post-hoc inference
  tar_target(parname_environmental_gaussian, setdiff(parname_environmental, c("major_init", "major_fix", "ba_init", "ba_fix", "ba_frac_init", "ba_frac_fix"))),
  tar_target(parname_environmental_ba, c("ba_init", "ba_fix", "ba_frac_init", "ba_frac_fix")), # contribname_environmental
  tar_target(parname_environmental_binomial, c("major_init", "major_fix")),
  
  tar_target(Surfaces_poly_env, predictPoly(cmdstanfit = fit_env, parname_Beta = parname_Beta, Envgrid = Envgrid_env)),
  tar_target(plot_poly_env, plotPoly(Surfaces_poly_env, Environmental = Environmental_env,
                                     basename = basename_fit_env, path = dir_publish, color = twocolors, themefun = theme_fagus)),
  
  tar_target(fit_environmental_env_gaussian,
             fitEnvironmental_glm(Environmental_env, parname = parname_environmental_gaussian, envname = predictor_select, taxon = taxon_s, fam = "gaussian", path = dir_publish),
             pattern = cross(parname_environmental_gaussian, taxon_s),
             iteration = "list"),
  tar_target(fit_environmental_env_ba,
             fitEnvironmental_gam(Environmental_env, parname = parname_environmental_ba, envname = predictor_select, taxon = taxon_s, fam = "gaussian", path = dir_publish),
             pattern = cross(parname_environmental_ba, taxon_s),
             iteration = "list"),
  tar_target(fit_environmental_env_binomial,
             fitEnvironmental_gam(Environmental_env, parname = parname_environmental_binomial, envname = predictor_select, taxon = 0, fam = "binomial", path = dir_publish),
             pattern = map(parname_environmental_binomial),
             iteration = "list"),
  tar_target(fit_environmental_env,
             c(fit_environmental_env_ba, fit_environmental_env_binomial), # excluded: fit_environmental_env_gaussian, i.e. all the pure parameters
             iteration = "list"),  
  
  
  tar_target(surface_environmental_env,
             predictEnvironmental(fit_environmental_env, envname = predictor_select,
                                  path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
             pattern = map(fit_environmental_env),
             iteration = "list"),
  
  # tar_target(surface_environmental_env_gaussian,
  #            predictEnvironmental(fit_environmental_env_gaussian, envname = predictor_select,
  #                                 path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
  #            pattern = map(fit_environmental_env_gaussian),
  #            iteration = "list"),
  # tar_target(surface_environmental_env_ba,
  #            predictEnvironmental(fit_environmental_env_ba, envname = predictor_select,
  #                                 path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
  #            pattern = map(fit_environmental_env_ba),
  #            iteration = "list"),
  # tar_target(surface_environmental_env_binomial,
  #            predictEnvironmental(fit_environmental_env_binomial, envname = predictor_select,
  #                                 path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
  #            pattern = map(fit_environmental_env_binomial),
  #            iteration = "list"),

  
  ## Plot
  tar_target(plot_environmental_env,
             plotEnvironmental(surfaces = surface_environmental_env, binaryname = "major_fix",
                               basename = basename_fit_env, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plot_environmental_env_gaussian,
             plotEnvironmental(surfaces = surface_environmental_env_gaussian, binaryname = "major_fix",
                               basename = basename_fit_env, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plot_environmental_env_ba,
             plotEnvironmental(surfaces = surface_environmental_env_ba, binaryname = "major_fix",
                               basename = basename_fit_env, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plot_environmental_env_binomial,
             plotEnvironmental(surfaces = surface_environmental_env_binomial, binaryname = "major_fix",
                               basename = basename_fit_env, path = dir_publish, color = twocolors, themefun = themefunction)),
  
  
  tar_target(plot_binary_par_env,
             plotBinary(Environmental = Environmental_env,
                        parname = str_to_sentence(parname_plotorder),
                        fit_bin = fit_environmental_env_binomial[[sapply(fit_environmental_env_binomial, function(x) attr(x, "par") == "major_fix") %>% which()]],
                        binarythreshold = 0.9,
                        path = dir_publish, basename = basename_fit_env,  color = twocolors, themefun = themefunction)),
  
  tar_target(plot_binary_contrib_env,
             plotBinary(Environmental = Environmental_env,
                        parname = contribname_env,
                        fit_bin = fit_environmental_env_binomial[[sapply(fit_environmental_env_binomial, function(x) attr(x, "par") == "major_fix") %>% which()]],
                        binarythreshold = 0.9,
                        path = dir_publish, basename = basename_fit_env,  color = twocolors, themefun = themefunction)),

  tar_target(plots_trace_env,
             plotTrace(cmdstanfit = fit_env, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  
  tar_target(plots_pairs_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plots_pairs_center_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_env_vertex_center)),
  tar_target(plots_pairs_spread_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_env_vertex_spread)),
  tar_target(plots_pairs_phi_env,
             plotPairs(cmdstanfit = fit_env, parname = "phi_obs_inv_sqrt")),
  
  tar_target(plots_parameters_env,
             plotParameters(stanfit = stanfit_env, parname = parname_plotorder, exclude = exclude, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction)),
  tar_target(plot_posterior_center_env,
             plotPosterior(cmdstanfit = fit_env, varname = parname_env_vertex_center)),
  tar_target(plot_posterior_spread_env,
             plotPosterior(cmdstanfit = fit_env, varname = parname_env_vertex_spread)),
  tar_target(plot_posterior_phi_env,
             plotPosterior(cmdstanfit = fit_env, varname = "phi_obs_inv_sqrt")),
  
  tar_target(plot_marginal_env,
             plotMarginal(Marginal = Marginal_env, parname = str_to_sentence(parname_plotorder),
                          path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = theme_fagus)),
  
  tar_target(plot_contributions_env,
             plotContributions(cmdstanfit = fit_env, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, plotlog = T,
                               color = twocolors, themefun = themefunction)),
  
  tar_target(plots_states_env,
             plotStates(States_env, allstatevars = c("ba_init", "ba_fix"),
                        path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction)),
  
  tar_target(plot_trajectories_avg_env,
             plotTrajectories(Trajectories_avg_env, thicker = T, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction))
  
  # tar_target(plots_env,
  #            plotStanfit(stanfit = stanfit_env, exclude = exclude, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction)),
  # tar_target(plots_predictions_posterior_env,
  #            plotPredictions(cmdstanfit = fit_env, data_stan_priors_offset_env, check = "posterior", path = dir_publish)),
  # tar_target(plots_conditional_env,
  #            plotConditional(cmdstanfit = fit_env, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  # tar_target(plot_contributions_prop_env,
  #            plotContributions(cmdstanfit = fit_env, parname = parname_plotorder, path = dir_publish, contribution = "sum_ko_prop", color = twocolors, themefun = themefunction)),
  # tar_target(plot_contributions_switch_env,
  #            plotContributions(cmdstanfit = fit_env, parname = parname_plotorder, plotlog = T, path = dir_publish, contribution = "sum_switch", color = twocolors, themefun = themefunction)),
  # tar_target(animation_trajectories_avg_env,
  #            animateTrajectories(plot_trajectories_avg_env, path = dir_publish, basename = basename_fit_env))
)



## Standalone generated quantities pipeline -----------------------------------
# targets_sgq <- list(
#   
#   tar_target(file_gq,
#              paste0(tools::file_path_sans_ext(file_model),"_gq.stan"), format = "file"),
#   tar_target(model_gq,
#              cmdstan_model(file_gq)),
#   tar_target(fit_gq,
#              model_gq$generate_quantities(fitted_params = fit$output_files(),
#                                           data = data_stan_priors,
#                                           output_dir = dir_fit,
#                                           parallel_chains = getOption("mc.cores", 4)))
#   
#   # tar_target(rstanfit_gq_test,
#   #            extractStanfit(fit_gq_test)),
#   # tar_target(draws_gq_test,
#   #            extractDraws(rstanfit_gq_test, exclude = helper_exclude)),
#   
# )



# ————————————————————————————————————————————————————————————————————————————————— #
# Outer pipeline -----------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

list(
  targets_settings,
  targets_paths,
  targets_wrangling,
  targets_parname,
  targets_fit_general,
  targets_fit_test,
  targets_fit,
  targets_fit_env,
  targets_posterior_test,
  targets_posterior,
  targets_posterior_env
)

