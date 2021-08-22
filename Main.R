# Library -----------------------------------------------------------------
library(here)
library(magrittr)
library(glue)
library(dplyr)
library(targets)
library(visNetwork)
library(future)

library(sf)
library(cmdstanr)
library(rstan)

# Orientation -------------------------------------------------------------
setwd(here())

# modelname <- ""
# modeldir <- dir(pattern = glue("^(Model).*{modelname}$"))
# modelpath <- file.path(modeldir, glue('Model_{modelname}.stan'))
# 
# datadir <- dir(pattern = "Data")


# Pre-targets sourcing ----------------------------------------------------
## assumed to have run and produced output files:
# source("Inventory.nosync/Main.R")

# Make targets pipeline -----------------------------------------------------------------
tar_glimpse()
M <- tar_manifest(fields = c("name", "command"))
# M$name

tar_make(c("data_stan"))
  ## alternatives
  # plan(multisession)
  # future(tar_make(names = "predict_splines")) # just as a future
  
  # tar_make_future(names = "Stages_splines") # parallel

# tar_load("Stages_splines")
# Stages %>% View()
# tar_read("Stages_env") %>% View()

# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T, exclude = contains(c("file_", "threshold_", "taxon_", "predictor_")))
network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T, blockShifting = T, parentCentralization = T)

tar_visnetwork(targets_only = F, exclude = starts_with("file"))
