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
tar_make_future(c("data_stan_priors_offset_env",
                  # "surfaceplots_s"
                  ), workers = if(onserver) 12 else 3, reporter = "timestamp")

## Fitting, parallelized internally
tar_make(c("fit_env", "summary_env"))

## Posterior
tar_make_future(c("summary_env",
                  "summary_states_env",
                  "residuals_env",
                  "plot_environmental_env",
                  "plots_trace_env",
                  "plots_pairs_env",
                  "plots_parameters_env",
                  # "plot_contributions_env",
                  # "plot_contributions_log_env",
                  # "plots_states_env",
                  # "plot_trajectories_avg_env",
                  "plots_conditional_env"),
                workers = if(onserver) 24 else 3, reporter = "timestamp")

## Publishing
## Publishing
tar_make_future(c("Summary_taxa",
                  "Summary_NFIs",
                  # "map_select",
                  # "surfaceplots_s"
                  ),
workers = if(onserver) 12 else 3, reporter = "timestamp")


# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T, exclude = contains(c("_test", "dir_", "threshold_", "taxon_", "predictor_", "pars_", "exclude", "parname", "simnames", "basename")))
network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T, blockShifting = T, parentCentralization = T)


# Load results ----------------------------------------------------------------
tar_load(c("summary", "fit"))
