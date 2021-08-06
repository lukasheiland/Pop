
# ——————————————————————————————————————————————————————————————————————————————————#
# Stage abundances functions -------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————#

## joinStages --------------------------------
# B  <- tar_read("Data_big")
# S  <- tar_read("Data_small")
# taxon_select <- tar_read("taxon_select")

joinStages <- function(B, S,
                       taxon_select,
                       threshold_dbh,
                       id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time")
) {
  
  id_select_S <- intersect(names(S), id_select)
  id_select_B <- intersect(names(B), id_select) %>% setdiff("treeid") ## make sure to always exclude treeid for good grouping
    ## B doesn't include clusterid!
  
  selectTaxa <- function(Abundance, taxon_select) {
    Abundance %>%
      mutate(tax = forcats::fct_other(tax, !!!as.list(taxon_select), other_level = "other")) %>%
      mutate(taxid = if_else(tax == "other", "other", as.character(taxid))) %>%
      droplevels()
  }
  
  S %<>%
    ### Lump taxa first, to facilitate summarize() and complete() later
    selectTaxa(taxon_select) %>%
    
    ### subset to time-constant regclasses, NOTE: There are more here.
    ## levels(S$regclass)
    ## table(S$regclass,S$obsid)
    ## here we select all 20cm height <= trees < 7mm dbh:
    ## c("h[20,50)", "h[50,130)", "hd[130,Inf)[0,5)", "d[5,6)", "d[6,7)")
    dplyr::filter(as.integer(regclass) %in% 1:5) %>%  
    ## these classes are all present in all three obsids, but consider different amounts of plots: # %>% dplyr::group_by(regclass, obsid) %>% dplyr::summarize(count_ha = sum(count/countarea))
    
    ## drop size classes and summarize counts
    dplyr::group_by_at(id_select_S) %>%
    dplyr::summarize(count_ha = sum(count/countarea)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(stage = factor("J")) %>%
    droplevels()
  ## S is alredy already completed!
  
  gc()
  
  ## interesting stuff like age, height, will get dropped here
  
  B %<>%
    ### Lump taxa first, to facilitate summarize() and complete() later
    selectTaxa(taxon_select) %>%
    
    ### Calculate basal area per hectare
    mutate(ba = pi * (dbh/2)^2 * 1e-6) %>% # mm^2 to m^2)
    mutate(count_ha = count) %>% # count is already per hectare, countarea == 1
    mutate(ba_ha = ba * count_ha) %>%
    
    ### make stages
    # quantile(B$dbh, seq(0, 1, by = 1e-1), na.rm = T): 160 is the 10%tile, 206 is the 20%tile
    ## 100mm…200mm == a
    mutate(stage = as.factor(if_else(dbh < threshold_dbh, "A", "B"))) %>%
    
    ## summarize per all ids, ... "plotobsid", "tax", "stage", being the fastest-varying
    group_by_at(c(id_select_B, "stage")) %>%
    dplyr::summarize(count_ha = sum(count_ha, na.rm = T), ba_ha = sum(ba_ha, na.rm = T)) %>% # in big trees, count is already per hectare
      ## this will yield NaN for stages which 
    ungroup() %>% #!!!
    
    ## taxa also have to be completed to include all stages, but only within obsid!!!
    ## (to avoid that 1. absence observations are fabricated where no observations where made, and 2. include observations of species that were not looked for in a survey)
    ## Note: method is a level that is maskedly reflecting stage, thus will prevent completion
    mutate(accountfortaxamethods = obsid) %>%
    tidyr::complete(nesting(accountfortaxamethods, taxid, tax), stage, nesting(obsid, plotid, plotobsid, methodid, time), fill = list(count_ha = 0, ba_ha = 0)) %>%
    filter(accountfortaxamethods == obsid) %>%
    select(-accountfortaxamethods) %>%
    droplevels() %>%
    
    ## stage creation will yield NA stages for when dbh is NA (for completed species)
    ## but any(is.na(B$dbh) & (B$count != 0)) # so tha, after completiont:
    tidyr::drop_na(stage)
  
  BA <- B %>%
    group_by_at(id_select_B) %>% ## without stage -> summarizes across stages within plot and tax
    dplyr::summarize(count_ha = sum(count_ha, na.rm = T), ba_ha = sum(ba_ha, na.rm = T)) %>% # in big trees, count is already per hectare
    mutate(stage = "BA")
  
  Stages <- dplyr::bind_rows(S, B, BA) %>%
    mutate(ba_ha = replace(ba_ha, is.nan(ba_ha), NA), count_ha = replace(count_ha, is.nan(count_ha), NA))
  
  return(Stages)
}



## addAbundanceSmooth --------------------------------
# BA_splines <- tar_read("BA_splines")
# BA <- BA_splines[[1]]
fitSplines <- function(BA) {
  
  
  X <- st_coordinates(BA)
  
  ## Thin plate spline regression with the coordinates as prediictors
  require(fields)
  fit <- fields::fastTps(X, BA$ba_ha, aRange = 20, lon.lat = T) # aRange = 30
  attr(fit, "taxon") <- attr(BA, "taxon")
  
  return(fit)
}


# Stages_env  <- tar_read("Stages_env")
# fits  <- tar_read("fits_splines")
predictSplines <- function(fits, Stages_env) {
  
  pred <- function(fit) {
    
    require("fields")

    Predicted <- predict(fit, xnew = X_unique) %>% # ?predict.fastTps
      as.data.frame() ## one column
    
    Predicted <- Predicted[X_df$id, , drop = F] ## through order, this has the right nrows
    
    colnames(Predicted) <- paste0("s_", attr(fit, "taxon"))
    return(Predicted)
  }
  
  X_df <- st_coordinates(Stages_env) %>%
    as.data.frame() %>%
    mutate(id = vctrs::vec_group_id(.))
  
  X_unique <- filter(X_df, !duplicated(id)) %>%
    select(1:2) %>%
    as.matrix()
  
  predictions <- lapply(fits, pred) ## list of data_frames
  
  S <- bind_cols(Stages_env, predictions)

  return(S) # plot(S["s_Fagus.sylvatica...24"])
}



## selectClusters --------------------------------
# Stages <- tar_read("Stages_splines")
# predictor_select <- tar_read("predictor_select")
selectClusters <- function(Stages, predictor_select,
                           id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time")
                           ) {

  ## for reference
  disturbance_select = c("standage_DE_BWI_1",
                         "allNatRegen_DE_BWI_2", "allNatRegen_DE_BWI_3", "anyUnnatRegen_DE_BWI_3", "anyUnnatRegen_DE_BWI_2",
                         "anyHarvested_DE_BWI_2", "anyHarvested_DE_BWI_3",
                         "anyForestryDamage_DE_BWI_2", "anyForestryDamage_DE_BWI_3")
  
  # Stages[predictor_select] %>% is.na() %>% colSums()
  
  ## Filter based environment
  Stages_select <- Stages %>%
    
    filter(!is.na(waterLevel_loc)) %>%
    filter(!is.na(alt_loc)) %>%
    filter(!is.na(phCaCl_esdacc)) %>%
    
    ## Plots were excluded that had any record of unnatural regeneration in DE_BWI_2 or 3
    mutate(anyUnnatRegen_DE_BWI_2 = tidyr::replace_na(anyUnnatRegen_DE_BWI_2, FALSE),
           anyUnnatRegen_DE_BWI_3 = tidyr::replace_na(anyUnnatRegen_DE_BWI_3, FALSE)) %>%
    filter( !(anyUnnatRegen_DE_BWI_2 | anyUnnatRegen_DE_BWI_2)) %>%
    
    ## Anything harvested
    filter( !(anyHarvested_DE_BWI_2 | anyHarvested_DE_BWI_3)) %>%
    
    ## Any damage through forestry
    filter( !(anyForestryDamage_DE_BWI_2 | anyForestryDamage_DE_BWI_3))
    
    ## unique(Stages_select$clusterid) %>% length()
    
  
  ## Filter clusters based on tree observations
  Stages_select %<>%
    
    group_by(clusterid) %>%
    
    ## subset to clusters with at least two surveys
    mutate(n_surveys = n_distinct(obsid)) %>%
    filter(n_surveys == 3) %>% # table(Stages_select$n_surveys) ## 2: 53626, 3: 119998
    dplyr::select(-n_surveys) %>%
    
    ## subset to clusters with at least two plots
    mutate(n_plots = n_distinct(plotid)) %>%
    filter(n_plots > 2) %>%
    dplyr::select(-n_plots) %>%
    
    ## get clusters with at least some Fagus
    mutate(anyFagus = any(count_ha > 0 && tax == "Fagus.sylvatica")) %>%
    
    ungroup()
  
  ## Subsample
  Stages_Fagus <- Stages_select %>%
    filter(anyFagus)
  
  n_Fagus <- length(unique(Stages_Fagus$clusterid))
  
  Stages_other <- Stages_select %>%
    filter(!anyFagus) %>%
    filter(clusterid %in% base::sample(unique(.$clusterid), n_Fagus))
  
  Stages_select <- bind_rows(Stages_other, Stages_Fagus)
  
  Stages_select %<>%
    select(-any_of(setdiff(disturbance_select, "standage_DE_BWI_1")))
  
  if(anyNA(Stages_select$time)) warning("selectClusters(): There are missing values in variable `time`.")
  
  ## Filter clusters based on succession
  # Stages_select %<>%
  #  group_by(clusterid) %>%
  
  ## Check what's left
  # unique(Stages_select$clusterid) %>% length()
  # unique(Stages_select$plotid) %>% length()
  # plot(Stages_select["alt_loc"])
  # plot(Stages_select["phCaCl_esdacc"])
  # plot(Stages_select["waterLevel_loc"])
  # hist(log(Stages_select$count_ha))
  # hist(Stages_select$waterLevel_loc)
  # hist(Stages_select$phCaCl_esdacc)
  # hist(Stages_select$alt_loc)
  
  
  ## Check whether there are plots with Fagus majority
  # Check <- Stages_select %>%
  #   filter(stage == "BA") %>%
  #   group_by(clusterid, tax) %>%
  #   summarize(count_ha = sum(count_ha, na.rm = T)) %>%
  #   group_by(clusterid) %>%
  #   summarize(prop_Fagus = count_ha[tax == "Fagus.sylvatica"] /sum(count_ha))
  # 
  # hist(Check$prop_Fagus)
  # table(Check$prop_Fagus > 0.5)

  ## Are there similar counts_ha per survey?
  # Check_consistency <- Stages_select %>%
  #   group_by(plotid, obsid, stage) %>%
  #   mutate(count_ha = sum(count_ha, na.rm = T), ba_ha = sum(ba_ha, na.rm = T))
  # 
  # boxplot(log(count_ha) ~ obsid * stage, data = Check_consistency) ## looks good, some succession maybe?
  # boxplot(count_ha ~ obsid * stage, data = Check_consistency, outline = F) ## looks good, even some succession maybe?
  
  attr(Stages_select, "predictor_select") <- predictor_select ## mainly for invalidating the target, when they have changed
  
  return(Stages_select)
}


# ——————————————————————————————————————————————————————————————————————————————————#
# Environment functions -------------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————#

# Manual data sanity checks -----------------------------------------------
# library(sf)
# library(leaflet)
# any(st_is_empty(st_geometry(E))) # Fine!
# basemap <- leaflet(E) %>%
#   fitBounds(3, 62, 10, 40) %>%
#   addTiles(urlTemplate = 'https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png')

## alt_loc [m]
# E[which(is.na(E$alt_loc)),]$plotid ## 1 plot without alt: "DE_BWI_25882.1"
# E[E$alt_loc < 0,]$alt_loc # 6 plots below alt = 0: "DE_BWI_23848.1" "DE_BWI_24084.4" "DE_BWI_24084.3" "DE_BWI_24145.4" "DE_BWI_24145.1"
# coords <- na.omit(st_coordinates(E[E$alt_loc < 0,]))
# addCircles(basemap, lng = coords[,1], lat = coords[,2], color = 'red')
## -> plausible

## aspect_loc [%]
# table(E$aspect_loc)

## phCaCl_esdacc
# hist(E$phCaCl_esdacc_esdacc)

## NAs
# E %>% is.na() %>% colSums() # plotobsid/time are purposefully NA in DE_BWI as the environment does only change with plotid, i.e. is constant over times (time).
# nacoord <- st_coordinates(E[is.na(E$phCaCl_esdacc),])
# addCircles(basemap, lng = nacoord[,1], lat = nacoord[,2], color = 'red') # NAs in most datasets are at the coasts, on land close to water.

## plotid is not necessarily the id for unique environments (plotid == envjoinid as it is DE_BWI).
## The coordinates do not change per re-observation of the same plot (extracted variables are dependent on coordinates), but the local variables aspect_loc and slope_loc do.



# Get all Ecken from the regular grid -------------------------------------
# Stages_env <- tar_read("Stages_env")
# Data_geo <- tar_read("Data_geo")

constructConstantGrid <- function(taxon, Stages_env, Data_geo) {
  
  ## All the clusters on the consistently sampled 4x4 km grid.
  Coords <- Data_geo %>%
    mutate(clusterid = paste0("DE_BWI_", Tnr)) %>%
    filter(Netz == 16) %>%  # 16 = 4km x 4 km grid, subset of the net (structure of four), used in all of the Länder in 2012 (see BMEL_BWI_Methodenband_Web_BWI3.pdf), readRDS("Inventory.nosync/DE BWI/Data/DE_BWI_meta.rds")$x_netz
    select(clusterid) # plot(Coords)
  
  ## Basal area of both stages combined ("BA") were averaged by the same location
  S <- dplyr::filter(Stages_env, tax == taxon &
                                 stage == "BA" &
                                 obsid == "DE_BWI_2012") %>% # only use the latest obsid because there was consistent sampling
    st_drop_geometry()
    
  expandMean <- function(x) {
    expanded <- rep(0, 4) # has only zeroes
    expanded[1:length(x)] <- x
    return(mean(x))
  }
  
  BA_fullgrid <-
    merge(Coords, S, by = "clusterid", all.x = T, all.y = F) %>% # only merge at locations of the constant grid
    select(clusterid, ba_ha) %>%
    st_centroid() %>%
    st_transform(crs = st_crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')) %>% # transform to lon/lat after centroid (important for spatial methods)
    mutate(ba_ha = tidyr::replace_na(ba_ha, 0)) %>%
  
    group_by(clusterid) %>%
    summarize(ba_ha = expandMean(ba_ha)) %>% # within clusterid, expand plots to 4 and fills with 0 before averaging
    
    mutate(coord = paste(geometry)) %>%
    group_by(coord) %>%
    summarize(ba_ha = mean(ba_ha))
    
    ## plot(S_fullgrid["ba_ha"])
  
  attr(BA_fullgrid, "taxon") <- taxon

  ## Return average ba_a per cluster with NA plots == 0 and NA clusters == 0, but only of the clusters with proptery 16
  ## Also subset to abundance from 2012
  ## Check whether the grid is really complete

  return(BA_fullgrid)
}


## Subset and clean the selected predictors --------------------------------
# E <- tar_read("Data_env")
# predictor_select <- tar_read("predictor_select")
cleanEnv <- function(E, predictor_select,
                     id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time"),
                     disturbance_select = c("standage_DE_BWI_1",
                                            "allNatRegen_DE_BWI_2", "allNatRegen_DE_BWI_3", "anyUnnatRegen_DE_BWI_3", "anyUnnatRegen_DE_BWI_2", 
                                            "anyHarvested_DE_BWI_2", "anyHarvested_DE_BWI_3",
                                            "anyForestryDamage_DE_BWI_2", "anyForestryDamage_DE_BWI_3")
                     ) {
  
  
  E %<>%
    dplyr::select(any_of(c(id_select, predictor_select, disturbance_select)))
  
  if("aspect_loc" %in% predictor_select) {
    
    E %<>% dplyr::mutate(aspect_loc = replace(aspect_loc, aspect_loc == 360, 0)) %>%
      dplyr::mutate(southexposition_loc = if_else(!(aspect_loc == 0 | is.na(aspect_loc)), # If there is a slope at all
                                                  -cos(pi * aspect_loc/180), # curve(-cos(pi * x/180), 0, 359)
                                                  0) # else 0.
      )
  }
  
  if("wwpi_cop" %in% predictor_select) {
    
    E %<>%
      mutate(wwpi_cop = na_if(wwpi_cop, 255))
  }
  
  return(E)
}



## summarizeEnvByCluster --------------------------------
# E <- tar_read("Env_clean")
# predictor_select <- tar_read("predictor_select")
summarizeEnvByCluster <- function(E,
                                  predictor_select = NULL,
                                  id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time"),
                                  disturbance_select = c("standage_DE_BWI_1",
                                                         "allNatRegen_DE_BWI_2", "allNatRegen_DE_BWI_3", "anyUnnatRegen_DE_BWI_3", "anyUnnatRegen_DE_BWI_2", 
                                                         "anyHarvested_DE_BWI_2", "anyHarvested_DE_BWI_3",
                                                         "anyForestryDamage_DE_BWI_2", "anyForestryDamage_DE_BWI_3")
                                  ) {
  
  if (is.null(predictor_select)) predictor_select <- setdiff(names(E), c(id_select, disturbance_select, "geometry"))
  
  E %<>%
    select(any_of(c(id_select, disturbance_select, predictor_select))) %>%
    dplyr::group_by(clusterid) %>%
    dplyr::mutate_at(predictor_select, function(x) mean(x, na.rm = T)) %>% # mean(NA, na.rm = T) produces NaN for only NAs
    dplyr::mutate_at(predictor_select, function(x) replace(x, is.nan(x), NA)) %>% # that's why ... ## na_if doesn't work with NaNs
    ungroup()
  
  return(E)
}


## join Env to Stages data --------------------------------
# Stages <- tar_read("Stages")
# E <- tar_read("Env_cluster")
joinEnv <- function(Stages, E) {
  
  dropcols <- setdiff(intersect(names(Stages), names(E)), "plotid") # drop all common cols, but plotid
  
  E <- dplyr::select(E, -all_of(dropcols))
  Stages %<>%
    left_join(E, by = "plotid") %>% # geometries are joined here, but crs is dropped
    st_sf(crs = st_crs(E))
  
  return(Stages)
}


## scale Predictors in Stages data --------------------------------
# Stages <- tar_read("Stages_select")
# predictor_select <- tar_read("predictor_select")

scalePredictors <- function(Stages,
                            predictor_select
                            ) {
  
  scalecolumns <- c(predictor_select) # , names(Stages)[grepl("^s_", names(Stages))]
  M_stages <- scale(st_drop_geometry(Stages[scalecolumns]))
  colnames(M_stages) <- paste0(scalecolumns, "_s")
  Stages <- cbind(Stages, M_stages)
  attr(Stages, "scaled:center") <- attr(M_stages, "scaled:center")
  attr(Stages, "scaled:scale") <- attr(M_stages, "scaled:scale")
  
  return(Stages)
}
