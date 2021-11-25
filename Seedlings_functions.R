# ——————————————————————————————————————————————————————————————————————————————————#
# Seedlings: Functions for fitting a regeneration model to SK seedling data --------
# ——————————————————————————————————————————————————————————————————————————————————#

## wrangleSeedlings --------------------------------
# Data_seedlings  <- tar_read("Data_seedlings")
# taxon_select <- tar_read("taxon_select")
# threshold_dbh <- tar_read("threshold_dbh")

wrangleSeedlings <- function(Data_seedlings, taxon_select = taxon_select, threshold_dbh = threshold_dbh) {
  
  if (taxon_select != "Fagus.sylvatica") stop("Prior for seedling regeneration rate r is only implemented for Fagus.sylvatica!")
  
  Data_seedlings <- Data_seedlings %>%
    mutate(tax = str_replace_all(taxon, " ", replacement = ".")) %>%
    mutate(tax = forcats::fct_other(tax, !!!as.list(taxon_select), other_level = "other")) %>%
    mutate(taxid = if_else(tax == "Fagus.sylvatica", "kew.83891", "other")) %>%
    
    # mutate(taxid = if_else(tax == "other", "other", as.character(taxid))) %>% ## if not only Fagus would be used
    droplevels()
  
  D_geo <- Data_seedlings %>%
    dplyr::select(plotid, WGS_E, WGS_N) %>%
    na.omit() %>%
    unique() %>%
    st_as_sf(crs = "+proj=longlat +datum=WGS84 +no_defs", coords = c("WGS_E", "WGS_N"))
  
  D_select <- Data_seedlings %>%
    group_by(plotid, year) %>%
    mutate(drop = any(regeneration == "Artificial", na.rm = T)) %>%
    ungroup() %>%
    filter(!drop)
  ## Any observations on plots were removed that had unnatural regeneration.
  
  D_filtered <- D_select %>%
    filter((sizeclass == "big" & dbh >= 100) |
             (sizeclass == "small" & height < 0.2) |
             is.na(sizeclass)) %>%
    filter(status != "Dead tree" | is.na(status)) %>%
    filter(!is.na(tax))
  ## The species-specific basal area of a plot was confined to living trees above the dbh >= 100mm.
  ## Only small trees with size class [10cm, 20cm), were included in the seedling counts.
  
  
  D_count <- D_filtered %>%
    group_by(plotid, year) %>%
    dplyr::mutate(ba_ha_sum = replace_na(sum(ba_ha, na.rm = T), 0), ## replaces NaNs!
                  ba_ha_sum_p1_inv = 1/(1+ba_ha_sum)) %>%
    group_by(plotid, year, tax, taxid) %>%
    dplyr::summarize(count_ha = sum(count_ha, na.rm = T),
                     ba_ha = sum(ba_ha, na.rm = T),
                     ba_ha_sum = first(ba_ha_sum),
                     ba_ha_sum_p1_inv = first(ba_ha_sum_p1_inv)) %>%
    mutate(count_ha = as.integer(round(count_ha)))
  
  D_count <- left_join(D_count, D_geo, by = "plotid")
  ## There are very few empty geometries
  isempty <- sf::st_is_empty(D_count$geometry)
  D_count <- D_count[!isempty,]

  # library(glmmTMB)
  # m <- glmmTMB::glmmTMB(count_ha ~ ba_ha * tax + 0, data = D_count, family = nbinom2)
  # summary(m)
  # res <- DHARMa::simulateResiduals(m)
  # plot(res)
  
  return(D_count)
}


# constructConstantGrid_SK -------------------------------------
# Seedlings <- tar_read("Seedlings")
constructConstantGrid_SK <- function(Seedlings) {
  
  Seedlings <- st_as_sf(Seedlings)
  d <- st_distance(D_geo)
  
  
  st_make_grid(Seedlings)
  
  table(Data_seedlings$plotid)
  
  Data_seedlings <- st_as_sf(Data_seedlings, coords = c("WGS_E", "WGS_N"))
  
  
  taxon <- as.character(taxon)
  
  ## All the clusters on the consistently sampled 4x4 km grid.
  Coords <- Data_geo %>%
    mutate(clusterid = paste0("DE_BWI_", Tnr)) %>%
    filter(Netz == 16) %>%  # 16 = 4km x 4 km grid, subset of the net (structure of four), used in all of the Länder in 2012 (see BMEL_BWI_Methodenband_Web_BWI3.pdf), readRDS("Inventory.nosync/DE BWI/Data/DE_BWI_meta.rds")$x_netz
    dplyr::select(clusterid) # plot(Coords)
  
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
    dplyr::select(clusterid, ba_ha) %>%
    st_centroid() %>%
    ## transformed for consistency with other data
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



## fitSeedlings --------------------------------
# Seedlings  <- tar_read("Seedlings")

fitSeedlings <- function(Seedlings) {
  
  ## count_ha = r*BA / (1+BA_sum)
  ## log(count_ha) = log(r * BA) + log(1/1+BA_sum)
  
  fit_seedlings <- brms::brm(count_ha ~ ba_ha + 1 + offset(ba_ha_sum_p1_inv) + (1 | plotid),
                             family = negbinomial,
                             data = Seedlings[Seedlings$tax == "Fagus.sylvatica",],
                             cores = getOption("mc.cores", 4))
  
  fit_seedlings_other <- brms::brm(count_ha ~ ba_ha + 1 + offset(ba_ha_sum_p1_inv) + (1 | plotid),
                                   family = negbinomial,
                                   data = Seedlings[Seedlings$tax == "other",],
                                   cores = getOption("mc.cores", 4))
  
  message("Summary of the the fit for Fagus seedlings:")
  print(summary(fit_seedlings))
  
  message("Summary of the the fit for other seedlings:")
  print(summary(fit_seedlings_other))
  
  return(list(Fagus.sylvatica = fit_seedlings, other = fit_seedlings_other))
  
}

