# Library -----------------------------------------------------------------
library(tidyverse)
library(ggformula)
library(magrittr)
library(glue)

library(rstan)
rstan_options(javascript = FALSE)
library(cmdstanr)
# install_cmdstan(cores = 3)

library(bayesplot)


# ——————————————————————————————————————————————————————————————————————————————————— #
# Wrangling functions  ----------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

## Make a matrix of draws into a plottable list ----------------------------
unpackDrawsArray <- function(M) {
  nrow <- dim(M)[2]
  ncol <- dim(M)[3]
  posteriors <- as.data.frame(t(apply(M, 1, c)))
  return(posteriors)
}


## tidyData ----------------------------------------------------------------
tidyData <- function(draws, setup) {
  
  Data <- setup$data$Long
  haslocdimensions <- sapply(draws, function(l) is.matrix(l) && ncol(l) == data$N_locs)
  draws_loc <- draws[haslocdimensions] # generated quantities on the loc scale
  draws_loc <- lapply(draws_loc, function(l) as.data.frame(l) %>% pivot_longer(everything(), names_to = "loc", names_prefix = "V"))
  Data_plot <- bind_rows(draws_loc, .id = "test")
  Data <- cbind(Data_plot, Data[match(Data_plot$loc, Data$loc), c("env_1", "env_2")])
  
  return(Data)
}

## Draw the priors specified in data, for use in marginal plots ----------------------
drawPriors <- function(data, n = 1000) {
  priorpars <- data[str_starts(names(data), "prior_")]
  isVertex <- str_starts(names(priorpars), "prior_Vertex")
  priorpars[isVertex] <- lapply(priorpars[isVertex], function(p) cbind(p, p)) # the cbind is for compatibility with the posteriors
  priorpars <- lapply(priorpars, function(P) apply(P, 2, function(col) rnorm(n, col[1], col[2])))
  priorpars <- lapply(priorpars, as.data.frame)
  return(priorpars)
}


# ——————————————————————————————————————————————————————————————————————————————————— #
# Plotting functions  ----------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————— #

## Dominance ---------------------------------------------------------------
plotDominance <- function(L) {
  
  L_sub <- L %>% filter(test %in% c("dominant_fix", "dominant_fix_s_0"))
  L_sub %>%
    ggplot(aes(x = env_2, y=value, colour = test)) +
    geom_point() +
    geom_smooth(formula = y ~ x)
  # geom_smooth(
  #   method="glm",
  #   method.args=list(family="binomial"),
  #   formula = cbind(dom, sub) ~ x
  # )
  
}


## plot time series in Data ----------------------
plotTimeSeries <- function(longdata) {
  
  longdata %>%
    mutate(sortid = paste(loc, plot, species, stage, sep = "_")) %>%
    arrange(sortid, time) %>%
    group_by(loc, species, stage, time) %>%
    mutate(abundance_loc = mean(abundance)) %>%
    group_by(loc, plot, species, stage) %>%
    mutate(ist0 = time == time[1], ist1 = time == time[2]) %>%
    ungroup() %>% 
    
    ggplot(mapping = aes(x = time, y = abundance, color = env_1, lty = stage, group = interaction(loc, stage, plot))) +
    geom_line() +
    facet_wrap(facets = c("species"))
  
}


#### Posterior vs prior vs. true ----------------------
plotEstimate <- function(prior,
                         posterior,
                         expected = NULL) {
  
  plot(density(posterior), xlim = range(prior))
  lines(density(prior), col = "lightblue")
  if(!is.null(expected)) abline(v = expected, col = "lightblue")
}


plotEstimates <- function(priors,
                          posteriors,
                          truevalues = NULL) {
  
  l <- length(priors)
  par(mfrow = c(ceiling(l/2), 2))
  for (i in 1:l) {
    plotEstimate(priors[[i]], posteriors[[i]], truevalues[[i]]) 
  }
  par(mfrow = c(1, 1))
}


