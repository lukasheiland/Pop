# Library -----------------------------------------------------------------
library(here)

library(tidyverse)
library(ggformula)
library(magrittr)
library(glue)
library(sf)

# install_cmdstan(cores = 3)
library(cmdstanr)


# Orientation -------------------------------------------------------------

modelname <- "c-threestage"

setwd(here())
modeldir <- dir(pattern = glue("^(Model).*{modelname}$"))
modelpath <- file.path(modeldir, glue('Model_{modelname}.stan'))

datadir <- dir(pattern = "Data")



# Load data ------------------------------------------------------------

## This is long and tidy. No NAs except for "time", completed times in "t"
Countdata <- readRDS(glue("{datadir}//Treestagecounts_sub.rds")) %>%
  
  ### PRELIMINARY
  group_by(clusterid) %>%
  filter(n_distinct(time) > 2)


# Make standata -----------------------------------------------------------


#### Return a differently long data sets structure
## The most comprehensive data set is N_y (resp. N_y0) with grouping locations/resurveys/pops/plots.
## NOTE THAT plots HAS TO BE THE LAST GROUP IN SORTING
## Everything is subset from the master subsets with these groupings (*_reobs, *_y0) and thus consistently sorted.
## Ratio: "resurveys" includes "pops" due of the return structure in ode_*(); "pops" includes "plots" because of the loc == population assumption (randowm effect from location/resurvey/pop to .../plots).
## Factor "pops" is structured stages/species.

## README!
## Not yet tested with new model structure
## Probably has to be updated. Look into sim_threestage.

pivotData <- function(D,
                      envnames = c("phCaCl_esdacc", "alt_loc"),
                      roundresponse = T,
                      format = c("standatalist", "long", "wide", "locs", "init", "y0", "y", "list")) {
  
  ## Correct completion within loc is assumed!
  
  if (roundresponse) {
    roundvars <- "abundance"
    D <- mutate_at(D, all_of(roundvars), round)
  }
  
  D <- group_by(D, loc) %>%
    mutate(isy0 = time == min(time)) %>%
    ungroup() %>%
    arrange(loc, time, pop, plot) # !!!
  
  ## Format: [N_y0] —   locations/pops(/stage/species)/plots
  D_y0 <- filter(D, isy0) %>% select(-isy0) %>%
    arrange(loc, time, pop, plot)
  
  ## Format: [N_y] — locations/resurveys/pops(/stage/species)/plots
  D_reobs <- filter(D, !isy0) %>% select(-isy0) %>%
    arrange(loc, time, pop, plot)
  
  ## Format: [N_init] — locations/pops
  D_init <- D_y0 %>%
    group_by(loc, pop, stage, species) %>%
    summarize(n_plots = n_distinct(plot), .groups = "drop")
  
  ## Format: [N_yhat] — locations/resurveys/pops
  D_yhat <- D_reobs %>%
    group_by(loc, time, pop, stage, species) %>%
    summarize(n_plots = n_distinct(plot), .groups = "drop")
  
  ## Format: [N_times] — locations/resurveys
  D_times <- D_reobs %>% # reobs to  count all times
    group_by(loc) %>%
    summarize(time = unique(time), .groups = "drop")
  
  ## Format: [N_locs] — locations
  D_locs <- D_reobs %>% # reobs to  count all times
    group_by(loc) %>%
    ## assumes completion within locations!
    summarize(time_max = max(time),
              n_species = n_distinct(species),
              n_plots = n_distinct(plot),
              n_pops = n_distinct(pop),
              n_reobs = n_distinct(time),
              n_yhat = n_distinct(interaction(pop, time)), .groups = "drop")
  
  ## Format: [N_locs]. Will be returned as an attribute of Data
  Env <- D[match(D_locs$loc, D$loc), envnames] %>%
    scale()

  ## CASE: locs
  if (match.arg(format) == "locs") D <- D_locs
  
  ## CASE: init
  if (match.arg(format) == "init") D <- D_init
  
  ## CASE: y0
  if (match.arg(format) == c("y0")) D <- D_y0
  
  ## CASE: y
  if (match.arg(format) == "y") D <- D_reobs
  
  if (match.arg(format) == "standatalist") {
    D_species <- D_init[c("loc", "species")] %>% unique()
    polyformula <- as.formula(paste("~", paste("poly(", colnames(as.data.frame(Env)), ", 2)", collapse = "+")))
    X <- model.matrix(polyformula, data = as.data.frame(Env))
    N_init <- nrow(D_init)
    N_yhat <- nrow(D_yhat)
    time_init <-  D$time[match(D_locs$loc, D$loc)] # N_locs!
    
    
    vrep <- function(...) unlist( Vectorize(rep.int)(...) ) # :)))))
    
    D <- list(
      N_locs = N_locs,
      N_init = N_init,
      N_times = N_times,
      N_yhat = N_yhat,
      N_y0 = nrow(Sims_y0),
      N_y = nrow(Sims_reobs),
      N_species = nrow(Sims_species), 
      
      
      N_totalspecies = length(unique(Sims_init$species)),
      
      N_beta = ncol(X),
      
      n_species = Sims_locs$n_species, 
      n_pops = Sims_locs$n_pops,
      n_reobs = Sims_locs$n_reobs,
      n_yhat = Sims_locs$n_yhat,
      
      ## repeat predictions on level "locations/pops" n_plots times to "locations/pops/plots"
      rep_init2y0 = vrep(1:N_init, Sims_init$n_plots),
      ## repeat predictions on level "locations/resurveys/pops" n_plots times to "locations/pops/resurveys/plots"
      rep_yhat2y = vrep(1:N_yhat, Sims_yhat$n_plots),
      ## repeat locations on level "locations" n_reobs times to "locations/pops/resurveys"
      rep_locs2times = vrep(1:N_locs, Sims_locs$n_reobs),
      
      species = Sims_species$species,
      stage = Sims_init$stage,
      time_init = time_init,
      time_max = Sims_locs$time_max,
      times = Sims_times$time,
      
      X = X,
      
      y0 = Sims_y0$abundance,
      y = Sims_reobs$abundance
      
    )
  }
  
  attr(D, "Env") <- Env

  return(D)
}


## rename data consistent with simulated data. ---------------
C <- Countdata %>%
  rename(loc = "clusterid", plot = "plotid", species = "agg", abundance = "count_ha") %>%
  mutate(time = t) %>%
  mutate(pop = interaction(species, stage))

## sample locs (clusters) ------------------------------------
set.seed(5)
loc_unique <- unique(C$loc)
loc_select <- sample(loc_unique, size = 50)

C_sub <- C %>%
  filter(loc %in% loc_select)
  
## pivot data into list format ------------------------------------
data <- pivotData(C_sub)


# Fit model -----------------------------------------------------------

## Model prep ---------------------------------
model <- cmdstan_model(modelpath) # model <- cmdstan_model(modelpath, force_recompile = T)


## Variational fit ----------------------------

## no convergence!
fit_var <- model$variational(data = data,
                             output_dir = "Fits",
                             init = 0,
                             # eta = 0.1,
                             iter = 5*10**3)



