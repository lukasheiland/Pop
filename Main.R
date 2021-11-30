# Library -----------------------------------------------------------------
library(targets)
source("_targets.R")
sapply(package, require, character.only = TRUE) ## package is a vector of all packages required in targets
library(visNetwork)
library(future)
library(future.callr)


# Pre-targets sourcing ----------------------------------------------------
## assumed to have run and produced output files:
# source("Inventory.nosync/Main.R")


# Make targets pipeline -----------------------------------------------------------------
# tar_glimpse()
# M <- tar_manifest(fields = c("name", "command"))
# M$name
# tar_watch(seconds = 5, outdated = FALSE, targets_only = TRUE)

tar_make("fit_test_pq")

tar_make_future(c("summary_test",
           "residuals_test",
           "plots_test",
           "plots_denscheck_prior_test",
           "plots_denscheck_posterior_test",
           "sensitivity_test",
           "plot_powerscale_test"),
           workers = 7, reporter = "timestamp")


## alternatives
# plan(multisession)
# future(tar_make(names = "predict_splines")) # just as a future

# tar_make_future(names = "Stages_s") # parallel


# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T, exclude = contains(c("file_", "threshold_", "taxon_", "predictor_", "pars_", "exclude", "parname", "simnames")))
network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T, blockShifting = T, parentCentralization = T)

## More
# tar_visnetwork(targets_only = F, exclude = starts_with("file"))


# Commit target meta to current branch ---------------------------------------------
system("git branch --show-current")
system("git commit -m 'targets meta' _targets/meta/meta")

# Inspect results ----------------------------------------------------------------
tar_load(c("summary_test", "stanfit_test"))
shinystan::launch_shinystan(stanfit_test)


