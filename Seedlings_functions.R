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
  
  ## predictions are on the log scale!
  # Seedlings_s %<>% 
  #   mutate(s_Fagus.sylvatica = exp(s_Fagus.sylvatica), s_other = exp(s_other))

  fitModel <- function(tax) {
    
   S  <- Seedlings_s[Seedlings_s$tax == tax,] %>% st_drop_geometry()
    
    data_seedlings <- list(
      N = nrow(S),
      # N_locs = length(unique(S$plotid)),
      # rep_loc = as.integer(as.factor(S$plotid)), ## These are all the same anyway!
      y = round(as.integer(S$count_ha)),
      l_smooth_log = S[, paste0("s_", tax), drop = T], ## predictions are on the log scale!
      ba_log = log(S$ba_ha) ## predictions are on the log scale!
    )
    
    fit_seedlings <- model_seedlings$sample(data = data_seedlings, parallel_chains = getOption("mc.cores", 4))
    
    var <- c("k_log", "l_log", "r_log", "phi")
    
    ggsave(file.path(fitpath, paste0("Pairs_Seedlings_", tax, ".png")),
           mcmc_pairs(fit_seedlings$draws(variables = var)))
    message("Summary of the the fit for ", tax, ":")
    print(fit_seedlings$summary(variables = var))
    
    return(fit_seedlings)
  }
  

  code_seedlings <- "
    data {
      int<lower=1> N;
      // int<lower=1> N_locs;
      // array[N] int<lower=1> rep_loc;
      array[N] int y;
      vector[N] l_smooth_log;
      vector[N] ba_log;
    }
    
    parameters {
      real k_log;
      // real<lower=0> sigma_k_loc;
      // vector[N_locs] k_loc_log_raw;
      real l_log;
      real r_log;
      real<lower=0> phi_inv_sqrt;
    }
    
    transformed parameters {
      real<lower=0> phi = inv_square(phi_inv_sqrt);
      vector<lower=0>[N] y_hat = exp(k_log) + // + k_loc_log_raw[rep_loc] * sigma_k_loc
                                 exp(l_log + l_smooth_log) +
                                 exp(r_log + ba_log); // exp(log(a) + log(b)) + exp(log(b) + log(c)) == a*b + b*c
    }
    
    model {
      // priors
      target += normal_lpdf(phi_inv_sqrt | 0, 1);
      // target += std_normal_lpdf(sigma_k_loc);
      // target += std_normal_lpdf(k_loc_log_raw);
      target += normal_lpdf(k_log | 3, 2);
      target += normal_lpdf(l_log | 0, 2);
      target += normal_lpdf(r_log | 1, 2);
      
      // likelihood
      target += neg_binomial_2_lpmf(y | y_hat, phi);
    }
  "
  model_seedlings <- cmdstanr::cmdstan_model(stan_file = write_stan_file(code_seedlings))
  
  # fit_seedlings <- brms::brm(count_ha ~ ba_ha + s_Fagus.sylvatica + 1 + (1 | plotid), # + offset(log(ba_ha_sum_p1_inv)), # + (1 | plotid),
  #                            family = negbinomial,
  #                            prior = set_prior("cauchy(0,10)", class = "sd", group = "plotid"),
  #                            data = Seedlings_s[Seedlings_s$tax == "Fagus.sylvatica",],
  #                            cores = getOption("mc.cores", 4))
  
  
  fits_seedlings <- lapply(c("Fagus.sylvatica", "other"), fitModel)
  
  return(fits_seedlings)
  
}

