# Library -----------------------------------------------------------------
library(here)
library(magrittr)
library(glue)
library(dplyr)
library(stringr)
library(targets)
library(visNetwork)
library(future)

library(sf)
library(MASS)
library(cmdstanr)
library(rstan)
library(bayesplot)

# Orientation -------------------------------------------------------------
setwd(here())


# Pre-targets sourcing ----------------------------------------------------
## assumed to have run and produced output files:
# source("Inventory.nosync/Main.R")


# Make targets pipeline -----------------------------------------------------------------
# tar_glimpse()
# M <- tar_manifest(fields = c("name", "command"))
# M$name
# tar_watch(seconds = 5, outdated = FALSE, targets_only = TRUE)
tar_make(c("summary_test", "plots_test", "stanfit_test"))
## alternatives
# plan(multisession)
# future(tar_make(names = "predict_splines")) # just as a future

# tar_make_future(names = "Stages_s") # parallel


# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T, exclude = contains(c("file_", "threshold_", "taxon_", "predictor_", "pars_")))
network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T, blockShifting = T, parentCentralization = T)

## More
# tar_visnetwork(targets_only = F, exclude = starts_with("file"))


# Inspect results ----------------------------------------------------------------
tar_load(c("summary_test", "stanfit_test"))
shinystan::launch_shinystan(stanfit_test)


