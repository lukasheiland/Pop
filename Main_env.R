# Library -----------------------------------------------------------------
library(targets)
library(visNetwork)
library(future)
library(future.callr)


# Sourcing ----------------------------------------------------
source("_targets.R")
sapply(package, require, character.only = TRUE) ## package is a vector of all packages required in targets

## This is, how the system is determined in _targets.R
# onserver <- Sys.info()["sysname"] != "Darwin"

## assumed to have run and produced output files:
# source("Inventory.nosync/Main.R", chdir = T)


# Make targets pipeline -----------------------------------------------------------------
# tar_glimpse()
# M <- tar_manifest(fields = c("name", "command"))

## Wrangling pipeline
tar_make_future(c("data_stan_priors_offset_env"), workers = if(onserver) 12 else 3, reporter = "verbose_positives")

## Fitting parallelized internally
tar_make(c("fit_env"))

## Posterior targets
tar_make_future(c("summary_env",
                  "summary_states_env",
                  "summary_marginal_env",
                  
                  "residuals_env",
                  # "plots_trace_env",

                  
                  "plots_parameters_env",
                  "plot_binary_par_env",
                  "plot_poly_env",
                  # "plot_marginal_env",
                  # "plot_posterior_center_env",
                  # "plot_posterior_spread_env",
                  # "plot_posterior_phi_env",
                  
                  
                  "plots_pairs_J_Fagus_env",
                  "plots_pairs_J_others_env",
                  "plots_pairs_A_env",
                  "plots_pairs_B_Fagus_env",
                  "plots_pairs_B_others_env"
                  # "plots_pairs_phi_env",
                  # "plots_pairs_env", ## includes parname_plotorder
                  
                  # "plots_conditional_env",
                  # "plots_states_env",
                  # "plot_predominant_env"
                  ),
                workers = if(onserver) 12 else 3, reporter = "timestamp_positives")

## Environmental targets
tar_make_future(c("plot_environmental_env", ## currently includes *_ba and *_binomial (both init and fix)
                  "plot_triptych_env",
                  "plot_diff_env",

                  # "plot_diff_supp_env",

                  "plot_diff_lim_env",
                  "plot_diff_L_env",
                  "plot_environmental_lim_env",
                  "plot_environmental_L_env"),
                workers = if(onserver) 12 else 3, reporter = "timestamp_positives")


## Trajectories: simulation parallelized internally
tar_make(c("plot_trajectories_avg_env"))


## Publishing targets
tar_make_future(c("Summary_NFI_env",
                  "Summary_range",
                  "Summary_taxa_env",
                  "plot_range"),
                  workers = if(onserver) 12 else 3, reporter = "verbose_positives")



# tar_make(c("map_select", "surfaceplots_s")) ## Make sure to have the development version of rayshader installed! For some reason, future does not work with rasterVis and ggplot addition (some problem with correct overloading of `+`)



# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T,
                          names = ends_with("_env"),
                          exclude = contains(c("_test", "dir_", "threshold_", "taxon_", "predictor_",
                                               "pars_", "exclude", "parname", "simnames", "basename"))
                          )

network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T,
                        blockShifting = T, parentCentralization = T)


# Load results ----------------------------------------------------------------
tar_load(c("summary_env", "fit_env"))
