# Library -----------------------------------------------------------------
library(here)

library(magrittr)
library(glue)


# Orientation -------------------------------------------------------------

modelname <- ""

setwd(here())
modeldir <- list.files(pattern = glue("^(Model).*{modelname}$"))
modelpath <- file.path(modeldir, glue('Model {modelname}.stan'))