
# Library -----------------------------------------------------------------
library(here)
library(tidyverse)
library(ggformula)
library(magrittr)
library(glue)

library(rstan)
rstan_options(javascript = FALSE)
library(cmdstanr)
# install_cmdstan(cores = 3)

library(bayesplot)

# Orientation -------------------------------------------------------------
setwd(here())
modeldir <- dir(pattern = glue("^(Model).*ba$"))

source(file.path(modeldir, "Sim_ba_helpers.R"))


# Load data ---------------------------------------------------------------
setupbasename <- "Model_ba-rect-202106092040"
setup <- readRDS(file.path("Fits.nosync", paste0(setupbasename, ".rds")))
stanfit <- rstan::read_stan_csv(file.path("Fits.nosync", setup$drawfile))
draws <- rstan::extract(stanfit)



# One series for viz ------------------------------------------------------
initialstate <- generateInitialState(n_species = pars$n_species)
times <- seq(1:100)

Sim <- simulateOneSeries(initialstate, times = times, pars = setup$truepars,
                          processerror = F, obserror = F)

matplot(Sim[,-1], type = "b", ylab="N",
        pch = rep(c("J", "A", "B"), each = setup$truepars$n_species),
        col = 1:setup$truepars$n_species, xlab = "time") # log='y'



# tidyData ----------------------------------------------------------------

tidyData <- function(draws, setup) {
  
  Data <- setup$data$Long
  haslocdimensions <- sapply(draws, function(l) is.matrix(l) && ncol(l) == data$N_locs)
  draws_loc <- draws[haslocdimensions] # generated quantities on the loc scale
  draws_loc <- lapply(draws_loc, function(l) as.data.frame(l) %>% pivot_longer(everything(), names_to = "loc", names_prefix = "V"))
  Data_plot <- bind_rows(draws_loc, .id = "test")
  Data <- cbind(Data_plot, Data[match(Data_plot$loc, Data$loc), c("env_1", "env_2")])
  
  return(Data)
  
}

L <- tidyData(draws, setup)



# Test fix point recovery -------------------------------------------------

## 0. inspect iterations
table(draws$iterations_fix)

## 1. Compare fix point
pseudofixpoint <- setup$pseudofixpointdata$y[, 2, 1,]
drawfixpoint <- apply(draws$state_fix[,, 1:6], c(2, 3), mean)
plot(drawfixpoint ~ pseudofixpoint, col = rowMeans(draws$iterations_fix))

## 3. rho_3 and rho_fix
hist(draws$rho_3)
hist(draws$rho_fix, add = T, col = 3)
diff_3 <- draws$y_hat[,,3,6] - draws$y_hat[,,2,6]
plot(draws$rho_3 ~ diff_3)
plot(draws$rho_fix ~ draws$state_fix[,,10])


# Dominance ---------------------------------------------------------------
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

