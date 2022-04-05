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
    dplyr::select(-taxon) %>%
    
    # mutate(taxid = if_else(tax == "other", "other", as.character(taxid))) %>% ## if not only Fagus would be used
    droplevels()


  D_filtered <- Data_seedlings %>%
    group_by(plotid, year) %>%
    mutate(drop = any(regeneration == "Artificial", na.rm = T)) %>%
    ungroup() %>%
    filter(!drop) %>%
    
    ## get only plots with at least some Fagus. This was discarded because in SK Fagus plots are much more often Fagus mono stands.
    # group_by(plotid) %>%
    # mutate(anyFagus = any(count_ha > 0 & tax == "Fagus.sylvatica", na.rm = T)) %>%
    # filter(anyFagus) %>%
    
    ## Any observations on plots were removed that had unnatural regeneration within a year and if they did not have any Fagus in any sizeclass and year
    
    ## implicitly I wanna filter to keep:
    ## EITHER (sizeclass == "small" & height < 0.2) -> this has been ensured in SK complete file before
    ## OR (sizeclass == "big" & dbh >= 100) -> to maintain completeness I replace count with 0
    mutate(dropSmallBig = (sizeclass == "big" & dbh < 100)) %>%
    ## and filter(status != "Dead tree")
    mutate(dropDead = (sizeclass == "big" & status == "Dead tree")) %>%
    
    mutate(drop = dropDead | dropSmallBig ) %>%
    droplevels()
  
  ## To maintain completenes within plot/year, instead of filtering
  D_filtered[which(D_filtered$drop), c("count", "count_ha", "ba", "ba_ha")] <- 0

  ## The species-specific basal area of a plot was confined to living trees above the dbh >= 100mm.
  ## Only small trees with size class [10cm, 20cm), were included in the seedling counts.
  
  ## According to the SK manual, the area of small tree counts was extended depending on tree density!
  
  D_count <- D_filtered %>%
    
    ## Here are several options for setting the offset for zero observations
    ## In the hurdle model, the fake offset area_0 will only get used as a factor level for determining the prob of getting a zero, not as an actual offset.
    
    group_by(sizeclass) %>%
    mutate(median_area = quantile(area, prob = 0.5, na.rm = T, type = 1)) %>%
    mutate(min_area = min(area, na.rm = T)) %>%
    mutate(max_area = max(area, na.rm = T)) %>% ## Zero observations are replaced with the max possible sampling area, because the circle of small tree was determined by small tree density: the greater the density the smaller the plot
    ungroup() %>%
    
    group_by(sizeclass, year, plotid) %>%
    mutate(n_area = n_distinct(area[sizeclass == "small"], na.rm = T),
           area_0 = case_when(n_area == 1 ~ first(area[!is.na(area)]), ## Here, even when there is a 0 the area of the other observed species is taken as a replacement.
                              n_area == 0 ~ first(max_area) ## This is here as an option to set an offset for 0-observation plots
                              # n_area > 1 ~ replace_na(weighted.mean(area, w = count_ha, na.rm = T), first(min_area)) ## This case is not present in the cleaned data. But this is an option to set the 0 to the average area, when there are multiple areas per plot.
                              )
           ) %>%
    
    ## In the subset dataset there are only plots with small tree observations with either 0 or 1 area
    ## D_count[D_count$sizeclass == "small",] %>% pull(n_area) %>% table()

    dplyr::mutate(count_ha_sum = replace_na(sum(count_ha, na.rm = T), 0) ) %>% ## replaces NaNs!
    dplyr::mutate(ba_ha_sum = replace_na(sum(ba_ha, na.rm = T), 0) ) %>%
    ungroup() %>%
    
    group_by(plotid, year, tax, taxid, sizeclass) %>%
    dplyr::summarize(count_obs = sum(count, na.rm = T),
                     area_0 = first(area_0),
                     area_obs = first(area[!is.na(area)]), ## with unique areas equivalent to: # area_obs = weighted.mean(area, w = count_ha, na.rm = T), # weighted.mean(area, w = count_ha, na.rm = T),
                     count_ha = sum(count_ha, na.rm = T),
                     ba_ha = sum(ba_ha, na.rm = T),
                     count_ha_sum = first(count_ha_sum),
                     ba_ha_sum = first(ba_ha_sum)
                     ) %>%
    ungroup() %>%
    
    mutate(offset = if_else(is.na(area_obs), area_0, area_obs)) %>%
    mutate(offset = if_else(is.na(area_obs), area_0, area_obs)) %>%
    
    ## for big trees, set the offset to NA, might be superflouous now, but here for variant reasons.
    mutate(offset =  if_else(sizeclass == "small", offset, as.numeric(NA)),
           area_obs = if_else(sizeclass == "small", area_obs, as.numeric(NA))) %>%
    
    mutate(count_ha_r = as.integer(round(count_ha))) %>%
    mutate(offset = round(offset, digits = 8)) ## For some reason the areas are sometimes slightly off, so that there are more unique offsets than there really are (compare table(offset) and n_distinct(offset))
  
  
  D_geo <- Data_seedlings %>%
    dplyr::select(plotid, WGS_E, WGS_N) %>%
    na.omit() %>%
    unique()
  
  D_count <- left_join(D_count, D_geo, by = "plotid") %>%
    filter(!(is.na(WGS_E)|is.na(WGS_N))) %>%
    st_as_sf(crs = "+proj=longlat +datum=WGS84 +no_defs", coords = c("WGS_E", "WGS_N"))   ## There are very few empty geometries, see in wrangleSeedlings_s()

  return(D_count)
}

## wrangleSeedlings_s --------------------------------
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
    filter( (sizeclass == "big" & dbh >= 100) | ## drop sizeclass == 'small' ## table(Data_seedlings_fullgrid$sizeclass)
            (sizeclass == "big" & is.na(dbh)) |  ## include the completed zero-observations
            is.na(sizeclass) ## include the cbound non-forest plots
            ) %>% 
    filter(status != "Dead tree" | is.na(status)) %>% ## %>%
  
    ##  filter(!is.na(tax))?
    ## Ba per hectare was first summed up per plot and Fagus/others
    group_by(plotid, year, tax, taxid) %>%
    dplyr::summarize(ba_ha = sum(ba_ha, na.rm = T)) %>%
    
    ## … and subsequently averaged per plot across years.
    group_by(plotid, tax, taxid) %>%
    dplyr::summarize(ba_ha = mean(ba_ha, na.rm = T)) %>%
    ungroup() %>%
    
    complete(plotid, nesting(tax, taxid), fill = list(ba_ha = 0)) %>% ## colSums(is.na()) ## there are only NAs in tax, one for each plot
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
  if(!dir.exists("Data")) dir.create("Data")
  saveRDS(Seedlings_s, file = path)
  
  message("Remember to upload: scp ~/Documents/Studium/Projects/Pop/Data/Seedlings_s.rds ssh 'lukasheiland@rhsbio7.uni-regensburg.de:/home/lukasheiland/Projects/Pop/Data'")
  
  return(path)
}


## fitSeedlings --------------------------------
# Seedlings_s  <- tar_read("Seedlings_s")
# fitpath  <- tar_read("dir_fit")
fitSeedlings <- function(Seedlings_s, fitpath) {

  #### fitModel() ---------------------
  #
  fitModel <- function(tax) {
    
    S <- Seedlings_s[Seedlings_s$tax == tax,] %>%
      st_drop_geometry() %>%
      dplyr::select(plotid, year, sizeclass, count_obs, count_ha, offset, ba_ha, count_ha_sum, ba_ha_sum, s_other, s_Fagus.sylvatica) %>%
      pivot_wider(id_cols = c("plotid", "year"),
                  names_from = "sizeclass",
                  values_from = c("count_obs", "count_ha", "ba_ha", "count_ha_sum", "ba_ha_sum", "offset", "s_other", "s_Fagus.sylvatica"))

    model_seedlings <- cmdstanr::cmdstan_model(stan_file = "Model_2021-03_ba/Model_seedlings.stan")

    data_seedlings <- list(
      N = nrow(S),
      y = as.integer(S$count_obs_small),
      
      ## In the hurdle model, the fake offset area_0 will only get used as a factor level for determining the prob of getting a zero. 
      offset_data = S$offset_small,
      # offset_scaled = c(scale(S$offset_small)),
      rep_offset = as.integer(as.factor(S$offset_small)),
      N_offset = n_distinct(S$offset_small),
      
      ba_sum = S$ba_ha_sum_big,
      l_smooth = exp(S[, paste0("s_", tax, "_big"), drop = T]), ## predictions are on the log scale!
      ba = S$ba_ha_big
    )
    
    getInits <- function(chain_id) return(list(l_log = 0.1, r_log = 0.1, theta_logit = rep(-1.0, data_seedlings$N_offset), m_logit = -2, phi_inv_sqrt = rep(10, data_seedlings$N_offset)))
    fit_Seedlings <- model_seedlings$sample(data = data_seedlings, parallel_chains = getOption("mc.cores", 4), output_dir = fitpath, init = getInits)
    
    attr(fit_Seedlings, "data") <- data_seedlings
    attr(fit_Seedlings, "taxon") <- tax
    
    var <- c("r_log", "k_log", "theta_logit", "m_logit", "phi")
    
    ggsave(file.path(fitpath, paste0("Pairs_Seedlings_", tax, ".png")),
           mcmc_pairs(fit_Seedlings$draws(variables = setdiff(var, "phi"))))
    
    s <- fit_Seedlings$summary(variables = var)
    write.csv(s, file.path(fitpath, paste0("Summary_Seedlings_", tax, ".csv")))

    message("Summary of the the fit for ", tax, ":")
    print(s)
    
    return(fit_Seedlings)
  }
  
  
  #### generateResiudals() ---------------------
  #
  generateResiduals <- function(taxon) {
    
    fit_Seedlings <- fits_Seedlings[[taxon]]
    data_seedlings <- attr(fit_Seedlings, "data")
    
    y_sim <- fit_Seedlings$draws(variables = "y_sim", format = "draws_matrix") %>% t()# matrix of observations simulated from the fitted model - row index for observations and colum index for simulations
    # y_hat <- fit_Seedlings$draws(variables = "y_hat", format = "draws_matrix") %>% apply(2, median, na.rm = T) ## This does not work for zero-altered models
    y <- data_seedlings$y
    
    # plot(y_hat/data_seedlings$offset ~ data_seedlings$ba)
    
    is0 <- !data_seedlings$y
    residuals <- DHARMa::createDHARMa(simulatedResponse = y_sim, observedResponse = y, integerResponse = T) # fittedPredictedResponse = y_hat
    
    # testZeroInflation(residuals)
    
    saveRDS(residuals, paste0(fitpath, "/Seedlings_", taxon, "_DHARMa", ".rds"))

    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa", ".png"), width = 1600, height = 1000)
    plot(residuals, quantreg = T, smoothScatter = F)
    dev.off()
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa_ba.sum", ".png"), width = 1600, height = 1000)
    plotResiduals(residuals, form = data_seedlings$ba_sum, quantreg = T, smoothScatter = F)
    dev.off()
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa_0", ".png"), width = 1600, height = 1000)
    plotResiduals(residuals, form = is0, quantreg = T, smoothScatter = F)
    dev.off()
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa_ba", ".png"), width = 1600, height = 1000)
    plotResiduals(residuals, form = data_seedlings$ba, quantreg = T, smoothScatter = F)
    dev.off()
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa_offset", ".png"), width = 1600, height = 1000)
    plotResiduals(residuals, form = data_seedlings$offset, quantreg = T, smoothScatter = F)
    dev.off()
    
    return(NULL)
  }
  
  fits_Seedlings <- sapply(c("Fagus.sylvatica", "other"), fitModel, USE.NAMES = T, simplify = F)
  sapply(names(fits_Seedlings), generateResiduals)

  return(fits_Seedlings)
  
}

