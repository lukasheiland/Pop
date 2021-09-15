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
# source("Fit_direct_functions.R")
source("Fit_functions.R")
# source("Model_2021-03_ba/Fit_ba_functions.R")

### Options
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("dplyr", "ggplot2", "tidyr", "magrittr", "glue", "forcats", "vctrs", "tibble", "stringr", # "multidplyr" ## extended tidyverse
                            "lubridate", # "zoo",
                            "sf", "raster", ## for correct loading of environmental data
                            "MASS",
                            "cmdstanr", "rstan", "bayesplot"))
addPackage <- function(name) { c(targets::tar_option_get("packages"), as.character(name)) }

### Future
library(future)
future::plan(future::multisession, workers = 6)


# Pipeline ----------------------------------------------------------------

list(
  
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
    tar_target(Data_geo, readRDS(file_DE_geo))
    
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
    ## Threshold to discriminate A and B [mm]
    # quantile(B$dbh, seq(0, 1, by = 1e-1), na.rm = T): 160 is the 10%tile, 206 is the 20%tile
    ## lower in the data is 100, so that: 100mm > A > 200mm > B
    tar_target(threshold_dbh, 200),
    tar_target(taxon_select, c("Fagus.sylvatica")),
    
    tar_target(Stages_transitions,
               countTransitions(Data_big, Data_big_status, Env_cluster, Stages_select, taxon_select = taxon_select, threshold_dbh = threshold_dbh)),
    
    tar_target(Stages,
               joinStages(Data_big, Data_small, taxon_select = taxon_select, threshold_dbh = threshold_dbh)),
    tar_target(Stages_env,
               joinEnv(Stages, Env_cluster)),
    
    list(
      tar_target(taxon_s,
                 c(taxon_select, "other")),
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
                 "Data/Stages_s.rds",
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
    
    tar_target(Stages_select,
               selectClusters(Data_Stages_s, predictor_select)), # Stages_s, After smooth, so that smooth can be informed by all plots.
               ## there is some random sampling here. Note: a target's name determines its random number generator seed. 
    tar_target(Stages_scaled,
               scaleData(Stages_select, predictor_select)) # After selection, so that scaling includes selected plots .
  ),
  
  
  ## Direct calibration
  # list(
  #   tar_target(Stages_direct,
  #              formatDirect(Stages_scaled)),
  #   tar_target(fit_direct,
  #              drawDirect(Stages_direct, predictor_select),
  #              packages = addPackage("brms")),
  #   tar_target(draws_direct,
  #              extractDrawsDirect(fit_direct)),
  #   tar_target(priors,
  #              constructPriors(draws_direct))
  # ),
  
  
  ## Model fit
  list(
    
    tar_target(data_stan,
               formatStanData(Stages_scaled, Stages_transitions, taxon_select, threshold_dbh)),
    
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
               formatPriors(data_stan, fit_g, fit_h)), # priors
    
    tar_target(file_model_test,
               "Model_2021-03_ba/Model_ba_test.stan",
               format = "file"),
    tar_target(file_model,
               "Model_2021-03_ba/Model_ba.stan",
               format = "file"),
    
    tar_target(model_test,
               cmdstan_model(file_model_test)),
    tar_target(model,
               cmdstan_model(file_model)),
    
    tar_target(fit_test,
               drawTest(model = model_test, data_stan_priors, method = "mcmc", initfunc = 0.5)),
    tar_target(fit,
               draw(model = model, data_stan_priors, initfunc = 0.5)),
    
    
    tar_target(summary_test,
               summarizeFit(fit_test)),
    tar_target(summary,
               summarizeFit(fit)),
    
    tar_target(stanfit_test,
               readStanfit(fit_test)),
    tar_target(stanfit,
               readStanfit(fit)),
    
    tar_target(pars_exclude,
               c("y_hat", "L_loc_log", "L_loc", "state_init_log",
                 "y_hat_rep", "converged", "iterations_fix", "state_fix", "dominant_fix", "major_fix", "fixiter_max")),
    
    tar_target(draws_test,
               extractDraws(stanfit_test)),
    tar_target(draws,
               extractDraws(stanfit)),
    
    tar_target(plots_test,
               plotStanfit(stanfit_test, exclude = pars_exclude)),
    tar_target(plots,
               plotStanfit(stanfit, exclude = pars_exclude))

  ),
  
  
  ## Generated quantities
  list(
    tar_target(gqmodel,
               cmdstan_model(paste0(tools::file_path_sans_ext(modelpath_stan),"_gq.stan"))),
    tar_target(gqfit,
               draw(gqmodel_stan, draws)),
    tar_target(gq,
               extractDraws(gqfit))
    # tar_stan_gq(gq,
    #             stan_files = gqmodelpath_stan,
    #             data = stages,
    #             fitted_params = draws)
  )
)

