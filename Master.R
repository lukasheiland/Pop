# Library -----------------------------------------------------------------
library(here)

library(magrittr)
library(glue)


# Orientation -------------------------------------------------------------
setwd(here())

modelname <- ""
modeldir <- dir(pattern = glue("^(Model).*{modelname}$"))
modelpath <- file.path(modeldir, glue('Model_{modelname}.stan'))

datadir <- dir(pattern = "Data")
