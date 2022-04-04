# Library -----------------------------------------------------------------
library(targets)
library(visNetwork)
library(future)
library(future.callr)
library(aria)

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
# M$name
# tar_watch(seconds = 5, outdated = FALSE, targets_only = TRUE)

## Wrangling pipeline
tar_make_future(c("data_stan_priors_offset"), workers = if(onserver) 12 else 3, reporter = "timestamp")

## Fitting, parallelized internally
tar_make(c("fit_test", "summary_test"))

## Posterior
# tar_make(c("Trajectories")) ## parallelized internally
tar_make_future(c("summary_test",
                  "summary_states_test",
                  "residuals_test",
                  "plots_test",
                  "plots_parameters_test",
                  "plots_conditional_test",
                  "plot_contributions_test",
                  "plot_contributions_prop_test",
                  "plots_states_test",
                  "plot_trajectories_avg_test"),
                workers = if(onserver) 12 else 3, reporter = "timestamp")


# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T, exclude = contains(c("dir_", "threshold_", "taxon_", "predictor_", "pars_", "exclude", "parname", "simnames", "basename")))
network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T, blockShifting = T, parentCentralization = T)


# Commit target meta to current branch ---------------------------------------------
onbranch <- system("git branch --show-current", intern = T)
oncorrectbranch <- xor(onserver & onbranch == "server",
                     !onserver & onbranch != "server")
if (oncorrectbranch) system("git commit -m 'targets meta' _targets/meta/meta") else message("Wrong branch for comitting meta!")


# Load results ----------------------------------------------------------------
tar_load(c("summary_test", "fit_test"))

