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
  
  
  D_geo <- Data_seedlings %>%
    dplyr::select(plotid, WGS_E, WGS_N) %>%
    na.omit() %>%
    unique()
  
  D_count <- left_join(D_count, D_geo, by = "plotid") %>%
    filter(!(is.na(WGS_E)|is.na(WGS_N))) %>%
    st_as_sf(crs = "+proj=longlat +datum=WGS84 +no_defs", coords = c("WGS_E", "WGS_N"))   ## There are very few empty geometries, see in wrangleSeedlings_s()

  # library(glmmTMB)
  # m <- glmmTMB::glmmTMB(count_ha ~ ba_ha * tax + 0, data = D_count, family = nbinom2)
  # summary(m)
  # res <- DHARMa::simulateResiduals(m)
  # plot(res)
  
  return(D_count)
}


## wrangleSeedlings --------------------------------
# Data_seedlings_fullgrid  <- tar_read("Data_seedlings_fullgrid")
# taxon_select <- tar_read("taxon_select")
# threshold_dbh <- tar_read("threshold_dbh")

wrangleSeedlings_s <- function(Data_seedlings_fullgrid, taxon_select = taxon_select, threshold_dbh = threshold_dbh) {
  
  if (taxon_select != "Fagus.sylvatica") stop("Prior for seedling regeneration rate r is only implemented for Fagus.sylvatica!")
  taxon_s <- c(taxon_select, 'other')
  
  D <- Data_seedlings_fullgrid %>%
    mutate(tax = str_replace_all(taxon, " ", replacement = ".")) %>%
    mutate(tax = forcats::fct_other(tax, !!!as.list(taxon_select), other_level = "other")) %>%
    mutate(taxid = if_else(tax == "Fagus.sylvatica", "kew.83891", "other")) %>%
    # mutate(taxid = if_else(tax == "other", "other", as.character(taxid))) %>% ## if not only Fagus would be used
    droplevels() %>%
    filter((sizeclass == "big" & dbh >= 100) | is.na(sizeclass)) %>% # drop sizeclass == 'small' ## table(Data_seedlings_fullgrid$sizeclass)
    filter(status != "Dead tree" | is.na(status)) %>% ## %>%
  
    ##  filter(!is.na(tax))?
    ## Ba per hectare was first summed up per plot and Fagus/others
    group_by(plotid, year, tax, taxid) %>%
    dplyr::summarize(ba_ha = sum(ba_ha, na.rm = T)) %>%
    
    ## … and subsequently averaged per plot across years.
    group_by(plotid, tax, taxid) %>%
    dplyr::summarize(ba_ha = mean(ba_ha, na.rm = T)) %>%
    ungroup() %>%
    
    complete(plotid, nesting(tax, taxid), fill = list(ba_ha = 0)) %>% ## colSums(is.na()) ## there are only NAs for 1 tax
    drop_na()
    
  D_geo <- Data_seedlings_fullgrid %>%
    dplyr::select(plotid, WGS_E, WGS_N) %>%
    na.omit() %>%
    unique()
  
  D_s <- left_join(D, D_geo, by = "plotid") %>%
    filter(!(is.na(WGS_E)|is.na(WGS_N))) %>%
    st_as_sf(crs = "+proj=longlat +datum=WGS84 +no_defs", coords = c("WGS_E", "WGS_N"))
  
  ## table(is.na(D_s$WGS_E), is.na(D_s$WGS_N)) ## There are very few empty geometries
  ## isempty <- sf::st_is_empty(D_s$geometry) # table(isempty)
  ## D_s <- D_s[!isempty,]
  # plot(D_s["ba_ha"])
  
  D_s %<>%
    mutate(coord = paste(geometry)) %>%
    group_by(coord, taxid, tax) %>%
    summarize(ba_ha = mean(ba_ha))
  
  subsetTaxon <- function(t, D = D_s) {
    O <- filter(D, tax == t)
    attr(O, "taxon") <- t
    return(O)
  }
  
  data <- lapply(taxon_s, subsetTaxon)
  
  return(data)
}


## predictSeedlingsSurfaces --------------------------------
# fits  <- tar_read("fits_Seedlings_s")

predictSeedlingsSurfaces <- function(fits) {
  
  SK <- raster::getData("GADM", country = "SK", level = 0, path = "Data/")
  R <- raster::raster(raster::extent(bbox(SK)), resolution = c(0.01, 0.01))
  crs(R) <- crs(SK)
  Coords <- rasterToPoints(R, spatial = TRUE) %>%
    as.data.frame()
  R$X <- Coords$x
  R$Y <- Coords$y
  R <- raster::mask(R, SK) # plot(R)
  
  # P <- raster::predict(R, fit); plot(P, col = viridis::viridis(255))
  surfaces <- lapply(fits, function(f) raster::predict(R, f, type = "response"))
  names(surfaces) <- sapply(fits, function(f) attr(f, "taxon"))
  
  return(surfaces)
}

saveSeedlings_s <- function(Seedlings_s) {
  path <- "Data/Seedlings_s.rds"
  saveRDS(Seedlings_s, file = path)
  return(path)
}


## fitSeedlings --------------------------------
# Seedlings_s  <- tar_read("Seedlings_s")

fitSeedlings <- function(Seedlings_s, fitpath = "Fits.nosync") {
  
  ## count_ha = r*BA / (1+BA_sum)
  ## log(count_ha) = log(r * BA) + log(1/1+BA_sum)
  
  fit_seedlings <- brms::brm(count_ha ~ ba_ha + s_Fagus.sylvatica + 0, # + offset(log(ba_ha_sum_p1_inv)), # + (1 | plotid),
                             family = negbinomial,
                             data = Seedlings_s[Seedlings_s$tax == "Fagus.sylvatica",],
                             cores = getOption("mc.cores", 4))
  ggsave(file.path(fitpath, "Pairs_Seedlings_Fagus.sylvatica.png"), pairs(fit_seedlings))
  message("Summary of the the fit for Fagus seedlings:")
  print(summary(fit_seedlings))
  
  fit_seedlings_other <- brms::brm(count_ha ~ ba_ha + s_other + 0, # + offset(log(ba_ha_sum_p1_inv)), # + (1 | plotid),
                                   family = negbinomial,
                                   data = Seedlings_s[Seedlings_s$tax == "other",],
                                   cores = getOption("mc.cores", 4))
  
  ggsave(file.path(fitpath, "Pairs_Seedlings_other.png"), pairs(fit_seedlings_other))
  message("Summary of the the fit for other seedlings:")
  print(summary(fit_seedlings_other))
  
  return(list(Fagus.sylvatica = fit_seedlings, other = fit_seedlings_other))
  
}

