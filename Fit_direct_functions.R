
# ——————————————————————————————————————————————————————————————————————————————————#
# Functions for fitting some rates directly ----------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————#


## formatDirect --------------------------------
# Stages  <- tar_read("Stages_scaled")

formatDirect <- function(Stages) {
  
  return(Stages)
}


## drawDirect --------------------------------
# Stages_direct  <- tar_read("Stages_direct")

drawDirect <- function(Stages_direct, predictor_select) {
  
  return(fit)
}


## extractDrawsDirect --------------------------------
# fit_direct  <- tar_read("fit_direct")

extractDrawsDirect <- function(fit_direct) {
  
  return(draws)
}


## constructPriors --------------------------------
# draws_direct  <- tar_read("draws_direct")

constructPriors <- function(draws_direct) {
  
  x <- log(rgamma(150,5))
  df <- approxfun(density(x))
  plot(density(x))
  xnew <- seq(-1, 3, by = 0.1)
  points(xnew,df(xnew),col=2)
  
  return(fit)
}

