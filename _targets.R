# Targets setup -----------------------------------------------------------
### Installation
# remotes::install_github("ropensci/targets")
# remotes::install_github("ropensci/tarchetypes")
# remotes::install_github("ropensci/stantargets")
# cmdstanr::install_cmdstan()


### Library
library(targets)
library(tarchetypes)
# library(stantargets)

### Source the functions
source("Wrangling_functions.R")
source("Seedlings_functions.R")
source("Fit_functions.R")

### Options
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("dplyr", "ggplot2", "tidyr", "magrittr", "glue", "forcats", "vctrs", "tibble", "stringr", # "multidplyr" ## extended tidyverse
                            "lubridate", # "zoo",
                            "sf", "raster", ## for correct loading of environmental data
                            "mgcv", "MASS",
                            "cmdstanr", "rstan", "brms", "bayesplot", "cowplot", "parallel", "DHARMa", "priorsense"))
addPackage <- function(name) { c(targets::tar_option_get("packages"), as.character(name)) }

### Future
library(future)
future::plan(future::multisession, workers = 6)


# Pipeline ----------------------------------------------------------------

## Inner pipelines

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
  
  ## Weakly informative priors.
  tar_target(weakpriors,
             ## Priors are organized like the parameter data structure but with an additional dimension in the case of a vector row of sds.
             list(
               prior_b_log = c(-2, 2),
               prior_c_a_log = c(-5, 2),
               prior_c_b_log = c(-5, 2),
               prior_c_j_log = c(-6, 2),
               ## prior_g_log,
               ## prior_h_log,
               prior_l_log = cbind(Fagus = c(0.5, 2), others = c(0.5, 2)),
               # prior_r_log = cbind(Fagus = c(0.5, 2), others = c(0.5, 2)),
               prior_s_log = c(-2, 2)
             )
  )
)


targets_parname <- list(
  
  tar_target(pars_exclude,
             c("y_hat", "L_loc_log", "L_loc", "state_init_log", "phi_obs_inv", "phi_obs_inv_sqrt")),
  tar_target(helpers_exclude,
             c("vector_b_log_prior", "vector_c_a_log_prior", "vector_c_b_log_prior", "vector_c_j_log_prior", "vector_s_log_prior",
               "area_zeta", "area_zeta_prior", "phi_obs_rep", "phi_obs_rep_prior", "zeta_prior_rep", "zeta_rep",
               "log_prior", "log_lik")),
  tar_target(rep_exclude,
             c("phi_obs_rep", "phi_obs_rep_prior", "zeta_prior_rep", "zeta_rep",
               "y_hat_rep", "y_hat_prior_rep", "y_hat_rep_offset", "y_hat_prior_rep_offset")),
  tar_target(simnames_prior,
             c("y_hat_prior", "y_hat_prior_rep", "y_hat_prior_rep_offset", "y_prior_sim")),
  tar_target(simnames_posterior,
             c("y_hat_rep", "y_hat_rep_offset", "y_sim", "y_hat_prior",
               "converged", "iterations_fix", "state_fix", "dominant_fix", "major_fix", "iterations_fix", "state_fix", "dominant_fix", "major_fix", "fixiter_max")),
  tar_target(exclude,
             c(pars_exclude, helpers_exclude, rep_exclude, simnames_prior, simnames_posterior)),
  tar_target(parname,
             c("phi_obs", "sigma_l", "zeta",
               "b_log", "c_a_log", "c_b_log", "c_j_log", "g_log", "h_log", "l_log", "r_log", "s_log"))
)


## The master pipeline
list(
  
  targets_settings,
  
  ## State data files
  # tar_load(starts_with("file"))
  list(
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
               format = "file")
    # tar_target(file_Taxa,
    #            'Inventory.nosync/Taxa/Taxa.csv',
    #            format = "file")
  ),
  
  ## Read data files
  # tar_load(starts_with("Data"))
  list(
    tar_target(Data_big, readRDS(file_DE_big)),
    tar_target(Data_big_status, readRDS(file_DE_big_status)),
    tar_target(Data_small, readRDS(file_DE_small)),
    tar_target(Data_env, readRDS(file_DE_env)),
    tar_target(Data_geo, readRDS(file_DE_geo)),
    
    tar_target(Data_seedlings, readRDS(file_SK))
    
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
  # tar_load(starts_with("S"))
  list(
    tar_target(Data_big_area,
               prepareBigData(Data_big, Data_big_status,
                              taxon_select = taxon_select, threshold_dbh = threshold_dbh, radius_max = radius_max)),
    
    tar_target(Data_small_area,
               prepareSmallData(Data_small, taxon_select = taxon_select)),
    
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
      
      tar_target(file_Stages_s,
                 "Data/Stages_s.rds", ## on other machine: saveStages_s(Stages_s),
                 format = "file"),
      tar_target(Data_Stages_s,
                 readRDS(file_Stages_s)),
      ## explicit side effect for later use on other machines
      
      tar_target(surfaces_s,
                 predictSurfaces(fits_s),
                 iteration = "list"),
      tar_target(surfaceplots_s,
                 plotSurfaces(surfaces_s),
                 iteration = "list")
    ),
    
    ## The seedlings pipeline
    list(
      tar_target(Seedlings,
                 wrangleSeedlings(Data_seedlings, taxon_select = taxon_select, threshold_dbh = threshold_dbh)),
      
      ## for inclusion of splines:
      # tar_target(Seedlings_BA_s,
      #            constructConstantGrid_SK(taxon_s, Seedlings),
      #            pattern = map(taxon_s),
      #            iteration = "list"),
      # tar_target(fits_Seedlings_s, ## fits_s each have an attribute "taxon"
      #            fitS(Seedlings_BA_s),
      #            pattern = map(Seedlings_BA_s),
      #            iteration = "list"),
      # tar_target(Seedlings_s,
      #            predictS(fits_Seedlings_s, Seedlings),
      #            iteration = "list"),
      # tar_target(file_Seedlings_s,
      #            saveStages_s(Seedlings_s),
      #            format = "file"),
      # tar_target(Data_Seedlings_s, ## explicit side effect for later use on other machines
      #            readRDS(file_Seedlings_s)),
      # tar_target(surfaces_Seedlings_s,
      #            predictSurfaces(fits_Seedlings_s),
      #            iteration = "list"),
      # tar_target(surfaceplots_Seedlings_s,
      #            plotSurfaces(surfaces_Seedlings_s),
      #            iteration = "list"),
      
      tar_target(fits_Seedlings,
                 fitSeedlings(Seedlings))
    ),
    
    
    tar_target(Stages_select,
               selectClusters(Data_Stages_s, predictor_select, selectpred = F)), # Data_Stages_s, After smooth, so that smooth can be informed by all plots.
    
    tar_target(Stages_select_pred,
               selectClusters(Data_Stages_s, predictor_select, selectpred = T)), ## Selection based on whether environmental variables are there
    
    tar_target(Stages_scaled,
               scaleData(Stages_select, predictor_select)), # After selection, so that scaling includes selected plots .
    
    tar_target(Stages_scaled_pred,
               scaleData(Stages_select_pred, predictor_select)) # After selection, so that scaling includes selected plots .
  ),
  
  
  ## Model fit
  list(
    
    tar_target(data_stan,
               formatStanData(Stages_scaled, Stages_transitions, taxon_s, threshold_dbh)),
    
    tar_target(file_model_transitions,
               "Model_2021-03_ba/Model_transitions.stan",
               format = "file"),
    
    tar_target(model_transitions,
               cmdstan_model(file_model_transitions)),
    
    tar_target(fit_g,
               fitTransition(data_stan, which = "g", model_transitions)),
    
    tar_target(fit_h,
               fitTransition(data_stan, which = "h", model_transitions)),
    
    tar_target(data_stan_priors,
               formatPriors(data_stan, weakpriors, fit_g, fit_h, fits_Seedlings, widthfactor = 2)), # priors
    
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
    
    tar_target(priorsim_test,
               drawTest(model = model_test, data_stan_priors, method = "sim", initfunc = 0.8)),
    tar_target(plot_denscheck_priorsim_test,
               plotDensCheck(cmdstanfit = priorsim_test, data_stan_priors, check = "prior")),
    
    tar_target(fit_test,
               drawTest(model = model_test, data_stan_priors, method = "mcmc",
                        n_chains = 6, iter_warmup = 800, iter_sampling = 500, initfunc = 0.8)),
    tar_target(fit,
               draw(model = model, data_stan_priors, method = "mcmc",
                    n_chains = 6, iter_warmup = 800, iter_sampling = 500, initfunc = 0.8)),
    
    tar_target(stanfit_test,
               readStanfit(fit_test)),
    tar_target(stanfit,
               readStanfit(fit)),
    
    targets_parname,
    
    tar_target(summary_test,
               summarizeFit(fit_test, exclude = c(helpers_exclude, rep_exclude))),
    tar_target(summary,
               summarizeFit(fit, exclude = c(helpers_exclude, rep_exclude))),
    
    tar_target(draws_test,
               extractDraws(stanfit_test, exclude = helpers_exclude)),
    tar_target(draws,
               extractDraws(stanfit, exclude = helpers_exclude)),
    
    ## Posterior plots
    tar_target(stanfit_test_plotting,
               readStanfit(fit_test, purge = TRUE)),
    tar_target(stanfit_plotting,
               readStanfit(fit, purge = TRUE)),
    
    tar_target(plots_test,
               plotStanfit(stanfit_test_plotting, exclude = exclude)),
    tar_target(plots,
               plotStanfit(stanfit_plotting, exclude = exclude)),
    
    ## Posterior predictive tests
    tar_target(residuals_test,
               scaleResiduals(cmdstanfit = fit_test, data_stan_priors)),
    tar_target(plot_denscheck_prior_test,
               plotDensCheck(cmdstanfit = fit_test, data_stan_priors, check = "prior")),
    tar_target(plot_denscheck_posterior_test,
               plotDensCheck(cmdstanfit = fit_gq_test, data_stan_priors, check = "posterior")),
    
    ## Sensitivity analysis
    tar_target(sensitivity_test, testSensitivity(fit_gq_test, include = parname)),
    tar_target(plot_powerscale_test, plotSensitivity(fit_gq_test, include = parname))
    
  ),
  
  
  ## Generated quantities
  list(
    
    tar_target(model_gq_test,
               cmdstan_model(paste0(tools::file_path_sans_ext(file_model_test),"_gq.stan"))),
    tar_target(fit_gq_test,
               model_gq_test$generate_quantities(fitted_params = stanfit_test,
                                                 data = data_stan_priors,
                                                 output_dir = "Fits.nosync/",
                                                 parallel_chains = getOption("mc.cores", 6))),
    # tar_target(rstanfit_gq_test,
    #            readStanfit(fit_gq_test)),
    # tar_target(draws_gq_test,
    #            extractDraws(rstanfit_gq_test, exclude = helpers_exclude)),
    
    tar_target(model_gq,
               cmdstan_model(paste0(tools::file_path_sans_ext(file_model),"_gq.stan"))),
    tar_target(fit_gq,
               model_gq$generate_quantities(fitted_params = stanfit,
                                            data = data_stan_priors,
                                            output_dir = "Fits.nosync/",
                                            parallel_chains = getOption("mc.cores", 6)))
    
  )
)

