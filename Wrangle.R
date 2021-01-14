# Library -----------------------------------------------------------------
library(here)

library(tidyverse)
library(magrittr)
library(glue)


# Orientation -------------------------------------------------------------
setwd(here())

# this references the forest inventory data directory, which is  another project
invdir <- '../Inventories/' 

datadir <- list.files(pattern = "Data")


# Load NFI data -----------------------------------------------------------
Small <- readRDS(glue('{invdir}/DE BWI/Data/DE_BWI_small_abund.rds')) # this is long and tidy and zero-completed.
Big <- readRDS(glue('{invdir}/DE BWI/Data/DE_BWI_big_abund.rds'))


# Load environmental data -------------------------------------------------
Env <- readRDS(glue("{invdir}/DE BWI/Data/DE_BWI_Env_sf.rds")
Env$lat <-  st_coordinates(Env)[,2]

predictor_select <- c("cwbYear_aclim",
                      # "cwbGdd0_aclim",
                      "gdd0_aclim" #,
                      # "tPeriodic2010_mh",
                      # "precPeriodic2010_mh",
                      # "clay_esdact",
                      # "sand_esdact",
                      # "octop_esdacoc",
                      # "cn_esdacc",
                      # "p_esdacc",
                      # "phCaCl_esdacc",
                      # "tMinColdMonth_wc2",
                      # "wwpi_cop"
                      # "lat"
)

idname <- c("envjoinid", "plotid", "plotobsid", "obsid", "clusterobsid", "clusterid", "time")

Env %<>% dplyr::select(all_of(c(idname, predictor_select)))

## plotid is not necessarily the id for unique environments (plotid == envjoinid as it is DE_BWI).
## As is, the predictors are available on the plotobsid level in SE_NFI. But here Env is only needed for plotids.
## The coordinates do not change per observation of the same plot (extracted variables are dependent on coordinates), but the local variables aspect_loc and slope_loc do.

Env %<>% dplyr::filter(!duplicated(plotid)) #!!! It's ok in this case. But only if environment is joined by plotid (see below).

## Save before NA handling.
saveRDS(Env, glue('{datadir}/Env.rds'))

## Count NAs
Env %>% is.na() %>% colSums() # plotobsid/time are purposefully NA in DE_BWI as the environment does only change with plotid, i.e. is constant over times (time).

## Explore NA locations
require(leaflet)
basemap <- leaflet(Env) %>% fitBounds(3, 62, 10, 40) %>% addTiles(urlTemplate = 'https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png')
nacoord <- st_coordinates(Env[is.na(Env$gdd0_aclim), "gdd0_aclim"])
addCircles(basemap, lng = nacoord[,1], lat = nacoord[,2], color = 'red') # NAs in most datasets are at the coasts, on land close to water.
## _esdacc NA are mostly close to water.
## _wc2 NA for coastline, etc.
## _cgiar NA for coastline
## _esdact NA close to coastlines and to Swiss border
## _mh NA for coastline

## IMPUTATION.
## Problem: Water-logged places often are NA in rasters.
## Solution: those should be easily predictable by features and coordinates.
imputeSpatial <-  function(Sf, varname_lhs, varname_rhs, ...){
  imputationformula <- paste(paste(varname_lhs, collapse = '+'), '~ X + Y +', paste(varname_rhs, collapse = '+')) %>%
    as.formula()
  Sf %<>% bind_cols(as.data.frame(st_coordinates(.)))
  geom <<- st_geometry(Sf)
  Sf %<>% st_drop_geometry() %>%
    missRanger(imputationformula, ...) %>%
    dplyr::select(-X, -Y) %>%
    st_set_geometry(geom)
  return(Sf)
}

varname_lhs <- predictor_select # The variables to be imputed. All variables have waterfront NAs.
varname_rhs <- predictor_select # The features to be predicted from. X, Y coordinates will be added by imputeSpatial().

Env_imputed <- imputeSpatial(Env, varname_lhs,
                             varname_rhs,
                             returnOOB = T,
                             seed = 1,
                             num.trees = 25)

## So many NAs
Env %>% is.na() %>% colSums() # plotobsid/time are purposefully NA in DE_BWI (Env changes only on plot level, not in time.)
## of
varname_lhs
## have been imputed with
varname_rhs
##. Has taken 4 iterations with num.trees = 25.
## With final average out of bag prediction error (oob):
# dput(attr(Env_imputed, "oob"))
# c(cwbYear_aclim = 0.000703498097850094, cwbGdd0_aclim = 0.000839791244127095, 
#   gdd0_aclim = 0.00070294745930576, tMinColdMonth_wc2 = 0.0006023741105646, 
#   wwpi_cop = 0.012081488887649, clay_esdact = 0.00190083851873737, 
#   sand_esdact = 0.00160083954955943, p_esdacc = 0.00581990457037878, 
#   cn_esdacc = 0.0019057081428244, phCaCl_esdacc = 0.00253252059989526)

## Save version with imputation.
saveRDS(Env_imputed, glue('{datadir}/Env_imputed.rds'))

## IMPUTED DATA IS NOT UTILIZED!
## !!!
# Env <- Env_imputed
rm(Env_imputed)


# Join environmental data -------------------------------------------------



# Subset NFI data ------------------------------------------------------------------

## Subset to target plots ----------------------------------------------------------

### filter DE
### filter species
### filter plots based on management


## Subset to consistent size classes over time ----------------------------------------------------------

### filter small size classes: DE consistent
### filter big size: DE consistent Kluppschwelle

# Bin counts into target size classes ------------------------------------------------------------------


