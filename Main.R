# Library -----------------------------------------------------------------
library(targets)
source("_targets.R")
sapply(package, require, character.only = TRUE) ## package is a vector of all packages required in targets
library(visNetwork)
library(future)
library(future.callr)

# Orientation -----------------------------------------------------------------
onserver <- Sys.info()["sysname"] != "Darwin"

# Pre-targets sourcing ----------------------------------------------------
## assumed to have run and produced output files:
# source("Inventory.nosync/Main.R")

# Make targets pipeline -----------------------------------------------------------------
# tar_glimpse()
# M <- tar_manifest(fields = c("name", "command"))
# M$name
# tar_watch(seconds = 5, outdated = FALSE, targets_only = TRUE)

tar_make(c("file_Seedlings_s", "file_Stages_s"))

tar_make_future(c("data_stan_priors_offset"),
                workers = if(onserver) 12 else 3, reporter = "timestamp")

tar_make(c("fit_test", "summary_test")) ## parallelized internally
# tar_make(c("Trajectories")) ## parallelized internally

tar_make_future(c("summary_test",
                  "residuals_test",
                  "plots_test",
                  "plots_denscheck_posterior_test",
                  "plots_conditional_test",
                  "plots_contributions_test",
                  "plots_predictions_posterior_test",
                  "plots_twostates_test",
                  "plot_trajectories_mean_test"),
                workers = if(onserver) 12 else 3, reporter = "timestamp")


# Inspect pipeline ----------------------------------------------------------------
network <- tar_visnetwork(targets_only = T, exclude = contains(c("dir_", "file_", "threshold_", "taxon_", "predictor_", "pars_", "exclude", "parname", "simnames", "basename")))
network %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 100, nodeSpacing = 120, edgeMinimization = T, blockShifting = T, parentCentralization = T)


# Commit target meta to current branch ---------------------------------------------
onbranch <- system("git branch --show-current", intern = T)
oncorrectbranch <- xor(onserver & onbranch == "server",
                     !onserver & onbranch != "server")
if (oncorrectbranch) system("git commit -m 'targets meta' _targets/meta/meta") else message("Wrong branch for comitting meta!")


# Load results ----------------------------------------------------------------
tar_load(c("summary_test", "fit_test"))

