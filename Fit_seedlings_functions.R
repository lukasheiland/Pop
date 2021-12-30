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
    
    ## get only plots with at least some Fagus
    group_by(plotid) %>%
    mutate(anyFagus = any(count_ha > 0 & tax == "Fagus.sylvatica", na.rm = T)) %>%
    filter(anyFagus) %>%
    
    ## Any observations on plots were removed that had unnatural regeneration within a year and if they did not have any Fagus in any sizeclass and year
    
    ## implicitly I wanna filter to keep:
    ## EITHER (sizeclass == "small" & height < 0.2) -> this has been ensured in SK complete file before
    ## OR (sizeclass == "big" & dbh >= 100) -> to maintain completeness I replace count with 0
    mutate(dropSmallBig = (sizeclass == "big" & dbh < 100)) %>%
    ## and filter(status != "Dead tree")
    mutate(dropDead = (status == "Dead tree")) %>%
  
    ## !!!
    ## filter(!(area > 0.07)) %>% # !!! these are incredible huge areas only for small trees
    ## !!!
    mutate(dropHugeArea = area > 0.07) %>%
    
    mutate(drop = dropDead | dropSmallBig | dropHugeArea) %>%
    droplevels()
  
  ## Too maintain completeness within plot/year, instead of filtering
  D_filtered[which(D_filtered$drop), c("count", "count_ha", "ba", "ba_ha")] <- 0
  D_filtered[which(D_filtered$dropHugeArea), c("area")] <- NA
  
  ## The species-specific basal area of a plot was confined to living trees above the dbh >= 100mm.
  ## Only small trees with size class [10cm, 20cm), were included in the seedling counts.
  
  ## According to the SK manual, the area of small tree counts is dependent on tree density!
  
  
  D_count <- D_filtered %>%
    group_by(year, sizeclass) %>%
    mutate(mode_area = as.numeric(names(sort(table(area), decreasing = T)[1]))) %>% ## Get the most frequent value of area for 0 observations. This was chosen to accomodate the fact that there are mostly the same areas, only a few big area seedlings plots.
    mutate(median_area = median(area, na.rm = T)) %>%
    mutate(max_area = max(area, na.rm = T)) %>%
    ungroup() %>%
    
    group_by(sizeclass, year, plotid) %>% # in 'small' this comes down to the same as group_by(sizeclass, tax, year) %>%
    mutate(n_area = n_distinct(area, na.rm = T), # there are no 0s (due to full completion), mostly 1s, and a few 2s
           area_0 = case_when(n_area == 1 ~ first(area[!is.na(area)]),
                              n_area == 0 ~ first(median_area),
                              n_area > 1 ~ replace_na(weighted.mean(area, w = count_ha, na.rm = T), first(median_area)))
           ) %>%
    # fill(area_0, .direction = "downup") %>%
    
    dplyr::mutate(count_ha_sum = replace_na(sum(count_ha, na.rm = T), 0) ) %>% ## replaces NaNs!
    dplyr::mutate(ba_ha_sum = replace_na(sum(ba_ha, na.rm = T), 0) ) %>%
    ungroup() %>%
    
    group_by(plotid, year, tax, taxid, sizeclass) %>%
    dplyr::summarize(count_obs = sum(count, na.rm = T),
                     area_0 = first(area_0),
                     area_obs = weighted.mean(area, w = count_ha, na.rm = T), # weighted.mean(area, w = count_ha, na.rm = T),
                     count_ha = sum(count_ha, na.rm = T),
                     ba_ha = sum(ba_ha, na.rm = T),
                     count_ha_sum = first(count_ha_sum),
                     ba_ha_sum = first(ba_ha_sum)
                     ) %>%
    ungroup() %>%
    
    mutate(offset = if_else(is.na(area_obs), area_0, area_obs)) %>%
    
    ## for bit trees, set the offset to NA
    mutate(offset =  if_else(sizeclass == "small", offset, as.numeric(NA)),
           area_obs = if_else(sizeclass == "small", area_obs, as.numeric(NA))) %>%
    
    mutate(count_ha_r = as.integer(round(count_ha)))
  
  
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
      pivot_wider(id_cols = c("plotid", "year"), names_from = "sizeclass", values_from = c("count_obs", "count_ha", "ba_ha", "count_ha_sum", "ba_ha_sum", "offset", "s_other", "s_Fagus.sylvatica"))
       # filter(ba_ha_sum_big > 0 & ba_ha_big > 0)
     
    code_seedlings <- "
    data {
      int<lower=1> N;
      
      int<lower=1> N_offset;
      vector<lower=0>[N] offset;
      vector[N] offset_scaled;
      array[N] int<lower=1> rep_offset;
      
      // int<lower=1> N_locs;
      // array[N] int<lower=1> rep_loc;
      array[N] int y;
      
      vector[N] l_smooth;
      vector[N] ba;
      vector[N] ba_sum;
    }
    
    parameters {
      real k_log;
      real r_log;
      // real l_log;
      // real d_log;
      // real s_log;
      
      real phi_inv_sqrt;
      real slope_phi;
    }
    
    transformed parameters {
      // vector<lower=0>[N_offset] phi = inv_square(phi_inv_sqrt);
      real phi = inv_square(phi_inv_sqrt);
      // vector<lower=0>[N] phi_hat = exp(log(phi) + slope_phi * offset_scaled);
      
      vector<lower=0>[N] y_hat_ha = (exp(k_log) + exp(r_log) * ba); // .* inv(1 + exp(s_log) * ba_sum) // exp(l_log) * l_smooth)
      
      
                                  // exp(k_log + r_log * log(ba) + s_log * log(ba_sum) + l_log * log(l_smooth) + offset);
      
                                 // exp(k_log - (s_log + ba_sum_log)) +
                                 // exp(d_log + offset_log) +
                                 // inv(1.0 + exp( + )) .*
                                 // exp(l_log + l_smooth_log - (s_log + ba_sum_log)) + 
                                 // exp(r_log + ba_log - (s_log + ba_sum_log));
                                 
                                 // r * BA *  // -  - ba_sum_log1p
                                 // exp(log(a) + log(b)) + exp(log(b) + log(c)) == a*b + b*c
                                 // prior inference model: seedlings ~ k/(s * ba_sum) + l * smooth + r * smooth; as we want to infer the seedling inputdue to ba and r
      
      vector<lower=0>[N] y_hat =  y_hat_ha .* offset; //  * 100;
    }
    
    model {
      // priors
      target += normal_lpdf(phi_inv_sqrt | 0, 5);
      // target += normal_lpdf(slope_phi | 0, 1);
      
      target += normal_lpdf(k_log | 0, 2);
      target += normal_lpdf(r_log | 0, 2);
      // target += normal_lpdf(l_log | 0, 2);
      // target += normal_lpdf(s_log | 0, 2);
      
      // likelihood
      target += neg_binomial_2_lpmf(y | y_hat, phi); //phi_hat
    }
    
    generated quantities {
      array[N] int y_sim = neg_binomial_2_rng(y_hat, phi); //phi_hat
    }
    "
    model_seedlings <- cmdstanr::cmdstan_model(stan_file = write_stan_file(code_seedlings))
    
    
    
    data_seedlings <- list(
      N = nrow(S),
      # N_locs = length(unique(S$plotid)),
      # rep_loc = as.integer(as.factor(S$plotid)), ## These are all the same anyway!
      y_r = round(S$count_ha_small),
      y = as.integer(S$count_obs_small),
      
      offset = S$offset_small,
      offset_log = log(S$offset_small),
      offset_scaled = c(scale(S$offset_small)),
      rep_offset = as.integer(as.factor(S$offset_small)),
      N_offset = n_distinct(S$offset_small),
      
      ba_sum_log = log(S$ba_ha_sum_big),
      ba_sum_log1p = log1p(S$ba_ha_sum_big),
      ba_sum = S$ba_ha_sum_big,
      l_smooth_log = S[, paste0("s_", tax, "_big"), drop = T], ## predictions are on the log scale!
      l_smooth = exp(S[, paste0("s_", tax, "_big"), drop = T]), ## predictions are on the log scale!
      ba_log = log(S$ba_ha_big),
      ba = S$ba_ha_big
    )
    
    fit_Seedlings <- model_seedlings$sample(data = data_seedlings, parallel_chains = getOption("mc.cores", 4), output_dir = fitpath)
    
    attr(fit_Seedlings, "data") <- data_seedlings
    attr(fit_Seedlings, "taxon") <- tax
    
    var <- c("r_log", "k_log") # "l_log", "s_log"
    
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
    y_hat <- fit_Seedlings$draws(variables = "y_hat", format = "draws_matrix") %>% apply(2, median, na.rm = T)
    y <- data_seedlings$y
    
    residuals <- DHARMa::createDHARMa(simulatedResponse = y_sim, observedResponse = y, fittedPredictedResponse = y_hat, integerResponse = T)
    
    saveRDS(residuals, paste0(fitpath, "/Seedlings_", taxon, "_DHARMa", ".rds"))
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa", ".png"), width = 1600, height = 1000)
    plot(residuals, quantreg = T, smoothScatter = F)
    dev.off()
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa_smooth", ".png"), width = 1600, height = 1000)
    plotResiduals(residuals, form = data_seedlings$l_smooth, quantreg = T, smoothScatter = F)
    dev.off()
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa_ba.sum", ".png"), width = 1600, height = 1000)
    plotResiduals(residuals, form = data_seedlings$ba_sum, quantreg = T, smoothScatter = F)
    dev.off()
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa_ba", ".png"), width = 1600, height = 1000)
    plotResiduals(residuals, form = data_seedlings$ba, quantreg = T, smoothScatter = F)
    dev.off()
    
    png(paste0(fitpath, "/Seedlings_", taxon, "_DHARMa_offset", ".png"), width = 1600, height = 1000)
    plotResiduals(residuals, form = data_seedlings$offset, quantreg = T, smoothScatter = F)
    dev.off()
    
    return(NULL)
  }
  
  
  # fit_Seedlings <- brms::brm(count_ha ~ ba_ha + s_Fagus.sylvatica + 1 + (1 | plotid), # + offset(log(ba_ha_sum_p1_inv)), # + (1 | plotid),
  #                            family = negbinomial,
  #                            prior = set_prior("cauchy(0,10)", class = "sd", group = "plotid"),
  #                            data = Seedlings_s[Seedlings_s$tax == "Fagus.sylvatica",],
  #                            cores = getOption("mc.cores", 4))
  
  fits_Seedlings <- sapply(c("Fagus.sylvatica", "other"), fitModel, USE.NAMES = T, simplify = F)
  sapply(names(fits_Seedlings), generateResiduals)

  return(fits_Seedlings)
  
}

