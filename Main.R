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
tar_make_future(c("data_stan_priors_offsets"), workers = if(onserver) 12 else 3, reporter = "timestamp")

## Fitting, parallelized internally
tar_make(c("fit", "summary"))

## Posterior
tar_make_future(c("summary",
                  "summary_states",
                  "residuals",
                  "plot_contributions",
                  "plot_contributions_log",
                  "plots_parameters",
                  "plots_states",
                  "plot_trajectories_avg",
                  "animation_trajectories_avg",
                  "plots",
                  "plots_conditional"),
                workers = if(onserver) 24 else 3, reporter = "timestamp")

## Publishing
tar_make_future(c("Summary_taxa",
                  "Summary_NFIs",
                  "map_select"),
                workers = if(onserver) 12 else 3, reporter = "timestamp")


# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T, exclude = contains(c("_test", "dir_", "threshold_", "taxon_", "predictor_", "pars_", "exclude", "parname", "simnames", "basename")))
network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T, blockShifting = T, parentCentralization = T)


# Load results ----------------------------------------------------------------
tar_load(c("summary", "fit"))
