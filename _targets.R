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

### Options
options(tidyverse.quiet = TRUE)
package <- c("dplyr", "ggplot2", "tidyr", "magrittr", "glue", "forcats", "vctrs", "tibble", "stringr", # "multidplyr" ## extended tidyverse
             "lubridate", "DescTools", # "zoo",
             "sf", "raster", "rasterVis", ## for correct loading of environmental data
             "eurostat", "elevatr", "rayshader", ## for mapping
             "mgcv", "MASS",
             "cmdstanr", "rstan", "brms", "posterior", "bayesplot", "parallel", "DHARMa", "priorsense",
             "cowplot",
             "future.apply")
tar_option_set(packages = package)

### Operation system
onserver <- Sys.info()["sysname"] != "Darwin"



# ————————————————————————————————————————————————————————————————————————————————— #
# Inner pipelines -----------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## Settings pipeline ------------------------------------------------------
targets_settings <- list(
  
  ## Threshold to discriminate A and B [mm]
  # quantile(B$dbh, seq(0, 1, by = 1e-1), na.rm = T): 160 is the 10%tile, 206 is the 20%tile
  ## lower in the data is 100, so that: 100mm > A > 200mm > B
  tar_target(threshold_dbh, 200), ## [mm]
  
  ## Upper sampling radius
  ## 	- All trees above a sampling radius of 14m were dropped, which is about the 98%tile (14.08m). The radius of 14m corresponds to the threshold radius of trees with dbh = 56cm
  ##    - dbh_threshold = radius_threshold/c with c == 25
  ##    - Alternatives: 99% radius == 1571 cm, 95% radius == 1188,  96% radius == 1242, 97% 1310.91
  tar_target(radius_max, 14000), ## [mm]
  
  ## Vector of taxa to select. All others will be lumped into "other".
  tar_target(taxon_select, c("Fagus.sylvatica")),
  
  ## Regeneration classes to include
  ## here we select all 20cm height <= trees < 7mm dbh: regglass_select <- c("h[20,50)" = 1, "h[50,130)" = 2, "hd[130,Inf)[0,5)" = 3, "d[5,6)" = 4, "d[6,7)" = 5)
  ## These are all size classes that are consistent across the three surveys.
  tar_target(regclass_select, c("h[20,50)" = 1, "h[50,130)" = 2, "hd[130,Inf)[0,5)" = 3, "d[5,6)" = 4, "d[6,7)" = 5)), ## [mm]
  
  ## Weakly informative priors.
  ## If provided, they are prioritzed over the fitted priors
  tar_target(weakpriors,
             ## Priors are organized like the parameter data structure but with an additional dimension in the case of a vector row of sds.
             list(
               prior_b_log = c(-3, 3),
               prior_c_a_log = c(-5, 3),
               prior_c_b_log = c(-5, 3),
               prior_c_j_log = c(-10, 3),
               # prior_g_log = cbind(Fagus = c(-1, 2), others = c(0, 2)),
               # prior_h_log = cbind(Fagus = c(-2, 3), others = c(-2, 3)),
               # prior_k_log = cbind(Fagus = c(4, 3), others = c(4, 3)),
               prior_l_log = c(9, 2),
               # prior_r_log = cbind(Fagus = c(4, 3), others = c(4, 3)),
               prior_s_log = c(-5, 3)
             )
          ),
  
  tar_target(twocolors, c("#208E50", "#FFC800"))
)



## Paths pipeline ------------------------------------------------------
targets_paths <- list(
  
  ## Directories
  tar_target(dir_fit, file.path("Fits.nosync/")),
  tar_target(dir_publish, file.path("Publish/")),
  
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
             format = "file")
  # tar_target(file_Taxa,
  #            'Inventory.nosync/Taxa/Taxa.csv',
  #            format = "file")
  
)



## Parameter spec pipeline ------------------------------------------------------
targets_parname <- list(
  
  tar_target(pars_exclude,
             c("y_hat", "L_loc_log", "K_loc_log_raw", "L_loc", "K_loc", "state_init_log", "phi_obs_inv", "phi_obs_inv_sqrt")),
  tar_target(helpers_exclude,
             c("Fix",
               "vector_b_log_prior", "vector_c_a_log_prior", "vector_c_b_log_prior", "vector_c_j_log_prior", "vector_s_log_prior",
               "phi_obs_rep", "phi_obs_rep_prior",
               "log_prior", "log_lik", "lp__", "state_init_log_raw")),
  tar_target(rep_exclude,
             c("phi_obs_rep", "phi_obs_rep_prior",
               "y_hat_rep", "y_hat_prior_rep", "y_hat_rep_offset", "y_hat_prior_rep_offset")),
  tar_target(simnames_prior,
             c("y_hat_prior", "y_hat_prior_rep", "y_hat_prior_rep_offset", "y_prior_sim")),
  tar_target(simnames_posterior,
             c("y_hat_rep", "y_hat_rep_offset", "y_sim", "y_hat_prior",
               "dominant_init", "major_init", "ba_init",
               "Fix", "dominant_fix", "major_fix", "ba_fix", "J_fix", "A_fix", "B_fix", 
               "converged_fix", "iterations_fix", "fixiter_max", "eps_ba_fix",
               "sum_ko_b_fix", "sum_ko_c_a_fix", "sum_ko_c_b_fix", "sum_ko_c_j_fix", "sum_ko_g_fix", "sum_ko_h_fix", "sum_ko_l_fix", "sum_ko_r_fix", "sum_ko_s_fix",
               "sum_ko_prop_b_fix", "sum_ko_prop_c_a_fix", "sum_ko_prop_c_b_fix", "sum_ko_prop_c_j_fix", "sum_ko_prop_g_fix", "sum_ko_prop_h_fix", "sum_ko_prop_l_fix", "sum_ko_prop_r_fix", "sum_ko_prop_s_fix",
               "greater_b", "greater_c_a", "greater_c_b", "greater_c_j", "greater_g", "greater_h", "greater_l", "greater_k", "greater_r", "greater_s")),
  tar_target(exclude,
             c(pars_exclude, helpers_exclude, rep_exclude, simnames_prior, simnames_posterior)),
  tar_target(parname,
             c("phi_obs", # "sigma_k_loc", # "k_log",
               "b_log", "c_a_log", "c_b_log", "c_j_log", "g_log", "h_log", "l_log", "r_log", "s_log")),
  tar_target(parname_loc,
             c("state_init_log", "L_loc")),
  tar_target(parname_sim,
             setdiff(parname, c("phi_obs", "sigma_k_loc")))
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
               c("alt_loc", "waterLevel_loc", "phCaCl_esdacc")),
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
                                taxon_select = taxon_select, threshold_dbh = threshold_dbh, radius_max = radius_max)),
    
    tar_target(Stages,
               joinStages(Data_big_area, Data_small_area, taxon_select = taxon_select, threshold_dbh = threshold_dbh)),
    tar_target(Stages_env,
               joinEnv(Stages, Env_cluster)),
    
    list(
      tar_target(taxon_s,
                 factor(c(sort(taxon_select), "other"))),
      tar_target(BA_s,
                 constructConstantGrid(taxon_s, Stages_env, Data_geo),
                 pattern = map(taxon_s),
                 iteration = "list"),
      tar_target(fits_s,
                 fitS(BA_s),
                 pattern = map(BA_s),
                 iteration = "list"),
      ## fits_s each have an attribute "taxon"
      tar_target(Stages_s,
                 predictS(fits_s, Stages_env),
                 iteration = "list"),
      
      tar_target(file_Stages_s, if (!! onserver) "Data/Stages_s.rds" else saveStages_s(Stages_s), format = "file"),
      tar_target(Data_Stages_s, readRDS(file_Stages_s)),
      
      tar_target(surfaces_s,
                 predictSurfaces(fits_s),
                 iteration = "list"),
      tar_target(surfaceplots_s,
                 plotSurfaces(surfaces_s, path = dir_publish),
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
                 fitS(seedlings_fullgrid),
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
      
      tar_target(file_Seedlings_s,  if(!! onserver) "Data/Seedlings_s.rds" else saveSeedlings_s(Seedlings_s), format = "file"),
      tar_target(Data_Seedlings_s, readRDS(file_Seedlings_s)),
      
      tar_target(fits_Seedlings,
                 fitSeedlings(Data_Seedlings_s, fitpath = dir_fit),
                 iteration = "list")
    ),
    
    
    tar_target(Stages_select,
               selectClusters(Data_Stages_s, predictor_select, selectpred = F)), # Data_Stages_s, After smooth, so that smooth can be informed by all plots.
    
    tar_target(Stages_select_pred,
               selectClusters(Data_Stages_s, predictor_select, selectpred = T)), ## Selection based on whether environmental variables are there
    
    tar_target(Stages_scaled,
               scaleData(Stages_select, predictor_select)), # After selection, so that scaling includes selected plots .
    
    tar_target(Stages_scaled_pred,
               scaleData(Stages_select_pred, predictor_select)), # After selection, so that scaling includes selected plots 
    
    ## Publishing
    tar_target(Summary_taxa,
               summarizeTaxa(Data_big, Data_seedlings, Stages_select, Seedlings_s, tablepath = dir_publish)),
    tar_target(Map_select,
               mapClusters(Stages_select, path = dir_publish))
    
  )
)



## Fitting pipeline ------------------------------------------------------
targets_fits <- list(
  tar_target(data_stan,
             formatStanData(Stages_scaled, Stages_transitions, taxon_s, threshold_dbh, timestep = 1, parfactor = 1)),
  
  tar_target(file_model_transitions,
             "Model_2021-03_ba/Model_transitions.stan",
             format = "file"),
  
  tar_target(model_transitions,
             cmdstan_model(file_model_transitions)),
  
  tar_target(fit_g,
             fitTransition(data_stan, which = "g", model_transitions, fitpath = dir_fit)),
  
  tar_target(fit_h,
             fitTransition(data_stan, which = "h", model_transitions, fitpath = dir_fit)),
  
  tar_target(data_stan_priors,
             formatPriors(data_stan, weakpriors, fit_g, fit_h, fits_Seedlings,
                          widthfactor_trans = 1, widthfactor_reg = 2)),
  
  tar_target(offsetname,
             c("offset", "offset_avg", "offset_q1", "offset_q3")[1]),
  
  tar_target(data_stan_priors_offset,
             selectOffset(offsetname, data_stan_priors)),
  
  tar_target(file_model_test,
             "Model_2021-03_ba/Model_ba_test.stan",
             format = "file"),
  tar_target(file_model,
             "Model_2021-03_ba/Model_ba.stan",
             format = "file"),
  
  tar_target(model_test,
             cmdstan_model(file_model_test) #, cpp_options = list(stan_opencl = TRUE)
             ),
  tar_target(model,
             cmdstan_model(file_model)),
  
  ## Prior predictive tests that rely on currently out-commented generated quantities
  # tar_target(priorsim_test,
  #            drawTest(model = model_test, data_stan = data_stan_priors, method = "sim", initfunc = 0.1, gpq = FALSE, fitpath = dir_fit)),
  # tar_target(plots_predictions_priorsim_test,
  #            plotPredictions(cmdstanfit = priorsim_test, data_stan_priors, check = "prior")),

  tar_target(fit_test_sansgq,
             fitModel(model = model_test, data_stan = data_stan_priors_offset, initfunc = 0.1, gpq = FALSE,
                      method = "mcmc", n_chains = 4, iter_warmup = 800, iter_sampling = 500, fitpath = dir_fit)),
  tar_target(fit_test,
             fitModel(model = model_test, data_stan = data_stan_priors_offset, initfunc = 0.1, gpq = TRUE,
                      method = "mcmc", n_chains = 4, iter_warmup = 800, iter_sampling = 500, fitpath = dir_fit)),
  tar_target(basename_fit_test,
             getBaseName(fit_test))
)



## Posterior pipeline ------------------------------------------------------
targets_posterior <- list(
  
  ## Extract
  tar_target(stanfit_test,
             extractStanfit(cmdstanfit = fit_test)),
  tar_target(draws_test,
             extractDraws(stanfit = stanfit_test, exclude = helpers_exclude)),
  # tar_target(stanfit_test_plotting,
  #            extractStanfit(fit_test_sansgq, purge = TRUE)),
  
  
  ## Summarize
  tar_target(summary_test,
             summarizeFit(cmdstanfit = fit_test, exclude = c(helpers_exclude, rep_exclude), path = dir_publish)),
  tar_target(Freq_converged_test,
             summarizeFreqConverged(cmdstanfit = fit_test, data_stan_priors, path = dir_publish)),
  
  
  ## Generate
  tar_target(residuals_test,
             generateResiduals(cmdstanfit = fit_test, data_stan_priors, path = dir_publish)),
  tar_target(Trajectories_test,
             generateTrajectories(cmdstanfit = fit_test, data_stan_priors, parname, locparname = parname_loc,
                                  time = c(seq(1, 491, by = 10), seq(500, 5000, by = 100)), thinstep = 50, usemean = F)),
  tar_target(Trajectories_mean_test,
             generateTrajectories(cmdstanfit = fit_test, data_stan_priors, parname, locparname = parname_loc,
                                  time = c(seq(1, 491, by = 10), seq(500, 5000, by = 100)), thinstep = 25, usemean = T)),
  
  ## Formatted posterior data stuctures
  tar_target(Twostates_test,
             formatTwoStates(cmdstanfit = fit_test, data_stan_priors)),
  
  
  ## Plot
  tar_target(plots_test,
             plotStanfit(stanfit = stanfit_test, exclude = exclude, path = dir_publish, basename = basename_fit_test)),
  ## Prior predictive tests that rely on currently out-commented generated quantities
  # tar_target(plots_predictions_prior_test,
  #            plotPredictions(cmdstanfit = fit_test, data_stan_priors, check = "prior")),
  tar_target(plots_predictions_posterior_test,
             plotPredictions(cmdstanfit = fit_test, data_stan_priors, check = "posterior", path = dir_publish)),
  tar_target(plots_conditional_test,
             plotConditional(cmdstanfit = fit_test, parname = parname_sim, path = dir_publish)),
  tar_target(plot_contributions_test,
             plotContributions(cmdstanfit = fit_test, parname = parname_sim, path = dir_publish, color = twocolors)),
  tar_target(plot_contributions_prop_test,
             plotContributions(cmdstanfit = fit_test, parname = parname_sim, path = dir_publish, plotprop = T, color = twocolors)),
  tar_target(plots_twostates_test,
             plotTwoStates(Twostates_test, path = dir_publish, basename = basename_fit_test)),
  tar_target(plot_trajectories_test,
             plotTrajectories(Trajectories_test, path = dir_publish, basename = basename_fit_test)),
  tar_target(plot_trajectories_mean_test,
             plotTrajectories(Trajectories_mean_test, thicker = T, path = dir_publish, basename = basename_fit_test)),
  
  # tar_target(plot_powerscale_test,
  #            plotSensitivity(cmdstanfit = fit_test, include = parname, path = dir_publish)),
  
  ## Test
  tar_target(sensitivity_test,
             testSensitivity(fit_test, include = parname, path = dir_publish))
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
#   #            extractDraws(rstanfit_gq_test, exclude = helpers_exclude)),
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
  targets_fits,
  targets_posterior
  
  )

