# Library -----------------------------------------------------------------
library(here)

library(magrittr)
library(glue)


# Orientation -------------------------------------------------------------
setwd(here())

# modelname <- ""
# modeldir <- dir(pattern = glue("^(Model).*{modelname}$"))
# modelpath <- file.path(modeldir, glue('Model_{modelname}.stan'))
# 
# datadir <- dir(pattern = "Data")


# Make targets pipeline -----------------------------------------------------------------
tar_make()

# Inspect pipeline ----------------------------------------------------------------
tar_manifest(fields = "command")
tar_glimpse()
tar_visnetwork()