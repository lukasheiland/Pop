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

## Posterior
tar_make_future(c("summary_env",
                  "summary_states_env",
                  "summary_marginal_env",
                  "residuals_env",
                  # "plots_trace_env",
                  
                  # "plot_environmental_env_gaussian",
                  # "plot_environmental_env_ba",
                  # "plot_environmental_env_binomial",
                  "plot_environmental_env", ## currently includes *_ba and *_binomial (both init and fix)
                  
                  "plot_poly_env",
                  
                  "plots_pairs_env",
                  "plots_pairs_env1_env",
                  "plots_pairs_env2_env",
                  # "plots_pairs_phi_env",
                  
                  "plot_marginal_env",
                  "plot_binary_par_env",
                  "plot_binary_contrib_env",
                  "plots_parameters_env",
                  "plot_posterior_center_env",
                  "plot_posterior_spread_env",
                  "plot_posterior_phi_env",
                  
                  "plot_contributions_env",
                  "plots_states_env",
                  "plots_conditional_env"),
                workers = if(onserver) 32 else 3, reporter = "verbose_positives")

## Simulations parallelized internally
tar_make(c("plot_trajectories_avg_env"))

## Publishing
tar_make_future(c("Summary_taxa",
                  "Summary_NFIs",
                  "Summary_range",
                  "plot_range"),
                workers = if(onserver) 12 else 3, reporter = "verbose_positives")



# tar_make(c("map_select", "surfaceplots_s")) ## For some reason, future does not work with rasterVis and ggplot addition (some problem with correct overloading of `+`)



# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T, exclude = contains(c("_test", "dir_", "threshold_", "taxon_", "predictor_", "pars_", "exclude", "parname", "simnames", "basename")))
network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T, blockShifting = T, parentCentralization = T)


# Load results ----------------------------------------------------------------
tar_load(c("summary", "fit"))
