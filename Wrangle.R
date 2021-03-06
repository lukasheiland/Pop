# Library -----------------------------------------------------------------
library(here)

library(tidyverse)
library(magrittr)
library(glue)

require(sf)


# Orientation -------------------------------------------------------------
setwd(here())

# this references the forest inventory data directory, which is  another project
invdir <- '../Inventories/' 

datadir <- list.files(pattern = "Data")



# Load environmental data -------------------------------------------------
Env <- readRDS(glue("{invdir}/DE BWI/Data/DE_BWI_Env_sf.rds"))


predictor_select <- c("alt_loc", "phCaCl_esdacc", "wwpi_cop")
idname <- c("envjoinid", "plotid", "plotobsid", "obsid", "clusterobsid", "clusterid", "time", "treeid")


Env %<>% dplyr::select(any_of(c(idname, predictor_select)))

## plotid is not necessarily the id for unique environments (plotid == envjoinid as it is DE_BWI).
## As is, the predictors are available on the plotobsid level in SE_NFI. But here Env is only needed for plotids.
## The coordinates do not change per observation of the same plot (extracted variables are dependent on coordinates), but the local variables aspect_loc and slope_loc do.

Env %<>% dplyr::filter(!duplicated(plotid)) #!!! It's ok in this case. But only if environment is joined by plotid (see below).



# Env NA handling -------------------------------------------------------------
## Save before NA handling.
saveRDS(Env, glue('{datadir}/Env.rds'))

## Count NAs
Env %>% is.na() %>% colSums() # plotobsid/time are purposefully NA in DE_BWI as the environment does only change with plotid, i.e. is constant over times (time).

## Explore NA locations
require(leaflet)

basemap <- leaflet(Env) %>% fitBounds(3, 62, 10, 40) %>% addTiles(urlTemplate = 'https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png')
nacoord <- st_coordinates(Env[is.na(Env$phCaCl_esdacc),])
addCircles(basemap, lng = nacoord[,1], lat = nacoord[,2], color = 'red') # NAs in most datasets are at the coasts, on land close to water.


# Env NA Imputation -----------------------------------------------------------

## Problem: Water-logged places often are NA in rasters.
## Solution: those should be easily predictable by features and coordinates.

# imputeSpatial <-  function(Sf, varname_lhs, varname_rhs, ...){
#   imputationformula <- paste(paste(varname_lhs, collapse = '+'), '~ X + Y +', paste(varname_rhs, collapse = '+')) %>%
#     as.formula()
#   Sf %<>% bind_cols(as.data.frame(st_coordinates(.)))
#   geom <<- st_geometry(Sf)
#   Sf %<>% st_drop_geometry() %>%
#     missRanger(imputationformula, ...) %>%
#     dplyr::select(-X, -Y) %>%
#     st_set_geometry(geom)
#   return(Sf)
# }

varname_lhs <- predictor_select # The variables to be imputed. All variables have waterfront NAs.
varname_rhs <- predictor_select # The features to be predicted from. X, Y coordinates will be added by imputeSpatial().

# Env_imputed <- imputeSpatial(Env, varname_lhs,
#                              varname_rhs,
#                              returnOOB = T,
#                              seed = 1,
#                              num.trees = 25)

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
# saveRDS(Env_imputed, glue('{datadir}/Env_imputed.rds'))

## IMPUTED DATA IS NOT UTILIZED!
## !!!
# Env <- Env_imputed
# rm(Env_imputed)


# Load NFI data -----------------------------------------------------------

####  these are long and tidy and zero-completed.
Small <- readRDS(glue('{invdir}/DE BWI/Data/DE_BWI_small_abund.rds'))
Big <- readRDS(glue('{invdir}/DE BWI/Data/DE_BWI_big_abund.rds'))



# Make size stages -------------------------------------------------------
## First subsetting here to enable rbinding the big data.

Small %<>%
  ### subset to constant regclasses, NOTE: There are more here.
  ## levels(Small$regclass)
  ## here we select all 20cm height <= trees < 5mm dbh
  filter(as.integer(regclass) %in% 1:3) # c("h[20,50)", "h[50,130)", "hd[130,Inf)[0,5)")

Small %<>% 
  group_by(plotid, obsid, taxid, # real groupings
           plotobsid, tax, inventory, methodid, time # redundant but desirable factors
  ) %>%
  summarize(count_ha = sum(count/countarea)) %>%
  ungroup() %>%
  mutate(stage = factor("j"), size = "small") %>%
  droplevels()

Big %<>%
  ### make stages
  # quantile(Big$dbh, seq(0, 1, by = 1e-1), na.rm = T): 160 is the 10%tile, 206 is the 20%tile
  ## 100mm-200mm == a
  mutate(stage = as.factor(if_else(dbh < 200, "a", "b"))) %>%
  group_by(plotid, obsid, taxid, stage, # real groupings
           plotobsid, tax, inventory, methodid, time # redundant but desirable factors
  ) %>%
  summarize(count_ha = sum(count)) %>% # in big trees, count is already per hectare
  mutate(size = "big") %>%
  
  ## taxa also have to be complete within stage, but only within obsid!!!
  complete(nesting(obsid, plotid, taxid), stage, fill = list(count_ha = 0)) %>%
  droplevels()

C <- bind_rows(Small, Big)
saveRDS(C, glue("{datadir}/Treestagecounts.rds"))
rm(Small, Big)
## SHORTCUT:
# C <- readRDS(glue("{datadir}/Treestagecounts.rds"))


# Manage species ----------------------------------------------------------

## which levels are in Small and Big?
## which taxa are even constant over time?
# table(C$size, C$obsid, C$tax)

## Join genera.
Taxa <- read.csv(file = glue('{invdir}/Taxa/Taxa.csv'), colClasses = c('factor')) # colClasses gets recycled, everything is a factor
Taxa_uid <- Taxa[!duplicated(Taxa$tax.id),] # tax.id is not unique in Taxa! Unique is however needed for left_join by tax.id (not by inventory specific ids)!

## Lump taxa.
C %<>%
  left_join(Taxa_uid[c('tax.id', 'tax.genus', 'tax.species')], # Taxa_uid because it does only contain unique rows per tax.id! Otherwise left_joint would join more rows.
           by = c(taxid = 'tax.id')) %>%
  mutate(agg = fct_collapse(tax.genus, Fagus = c("Fagus"), Quercus = ("Quercus"), Pinus = c("Pinus"),
                            other_level = "other"))

## Aggregate again.
C %<>%
  group_by(plotid, obsid, agg, stage, # real groupings
           plotobsid, inventory, methodid, time # redundant but desirable factors
  ) %>%
  summarize(count_ha = sum(count_ha)) %>%
  ungroup()


# Join environmental data (with cluster info) and average by cluster (tract) ------------------------------

C %<>%
  left_join(Env, by = "plotid") %>%
  group_by(clusterid) %>%
  mutate_at(predictor_select, mean, na.rm = T) %>% # produces NaNs for all-NA clusters
  na.omit() # also drops NaNs


# Subset data based on community ------------------------------------------------------------------
C %<>%
  group_by(obsid, agg, stage, # real groupings
           plotobsid, inventory, methodid, time # redundant but desirable factors
  ) %>%
  summarize(count_ha = sum(count_ha)) %>%
  ungroup()



## Subset to target plots ----------------------------------------------------------

## filter plots based on management



