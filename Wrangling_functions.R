# ——————————————————————————————————————————————————————————————————————————————————#
# Tree data functions -------------------------------------------------------
# ——————————————————————————————————————————————————————————————————————————————————#

## prepareBigData --------------------------------
# B  <- tar_read("Data_big")
# B_status  <- tar_read("Data_big_status")
# taxon_select <- tar_read("taxon_select")
# threshold_dbh <- tar_read("threshold_dbh")
# radius_max <- tar_read("radius_max")
# tablepath  <- tar_read("dir_publish")
prepareBigData <- function(B, B_status,
                           taxon_select,
                           threshold_dbh, radius_max,
                           id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time"),
                           tablepath
                           ) {
  
  ## Areas and radii
  radius_max_B_cm <- radius_max/10 ## conversion from mm to cm
  radius_max_A_cm <- 25 * threshold_dbh/10 # [cm]
  dbh_max_B <- (radius_max_B_cm/25) * 10
  
  area_max_B <- pi * radius_max_B_cm^2 * 1e-8 # [ha]
  area_max_A <- pi * radius_max_A_cm^2 * 1e-8 # [ha]
  
  ## The average area is the integral on 0 to radius_max of the function f(r) = pi * r^2, divided by radius_max
  area_0_B <- (pi * radius_max_B_cm^3 * 1e-8) / (3 * radius_max_B_cm) ## [cm^2 * 1e-8 == ha]
  area_0_A <- (pi * radius_max_A_cm^3 * 1e-8) / (3 * radius_max_A_cm) ## [cm^2 * 1e-8 == ha]
  
  M <- matrix(c(radius_max_A_cm, radius_max_B_cm,
                area_max_A, area_max_B,
                area_0_A, area_0_B,
                threshold_dbh, dbh_max_B
                ),
              nrow = 4, byrow = T,
              dimnames = list(c("radius_max[cm]", "area_max [ha]", "area_0 [ha]", "dbh_max [mm]"), c("A", "B")))
  
  write.csv(M, file.path(tablepath, "Areas_truncation.csv"))
  print(M)
  
  id_select_B <- intersect(names(B), id_select) %>% setdiff("treeid") ## make sure to always exclude treeid for good grouping
  ## B doesn't include clusterid!
  
  
  selectTaxa <- function(Abundance, taxon_select) {
    Abundance %>%
      mutate(tax = forcats::fct_other(tax, !!!as.list(taxon_select), other_level = "other")) %>%
      mutate(taxid = if_else(tax == "other", "other", as.character(taxid))) %>%
      droplevels()
  }
  
  floor2 <- function(x) {
    if_else(x < 1, ceiling(x), floor(x))
  }
  
  ### Prepare joining status data for radii
  S <- B_status %>%
    rename(count = count_ha) %>% # just temporarily for consistency with Data_big
    # mutate(excluded = dead | harvested) %>%
    mutate(treejoinid = interaction(obsid, treeid))
  
  ## For formulas:
  # see: http://wiki.awf.forst.uni-goettingen.de/wiki/index.php/Bitterlich_sampling#Estimation_of_number_of_stems
  # with Zählfaktor k = 4 and c = 50/sqrt(k): c == 50/sqrt(4) == 25
  # expansion factor n_ha = 10000 * 1000^2/(pi * c^2 * dbh^2) [factor 1000^2 for dbh in mm instead of m]
  # r_virtualPlot == c * dbh

  B %<>%
    ### Lump taxa first, to facilitate summarize() and complete() later
    selectTaxa(taxon_select) %>%
    
    # join distances
    mutate(treejoinid = interaction(obsid, treeid)) %>%
    bind_cols(distance = S$distance[match(.$treejoinid, S$treejoinid)]) %>%
    filter( !(plotid %in% c("DE_BWI_55329.2", "DE_BWI_55329.3"))) %>%
    # 24:978433 trees have no associated radius, these are from two plots within one cluster.
    # B_temp$distance[B_temp$count > 0] %>% is.na() %>% table()
    # B_temp %>% filter(is.na(distance)) %>% filter(count > 0) %>%  View()
    # These will get dropped in selection later anyway.
    # B_temp %>% filter(is.na(distance)) %>% filter(count > 0) %>% pull(plotid) %>% unique() %>% dput()

    ## Explore quantile for potential maximum distances
    # B_temp %>% pull(distance) %>% quantile(0.98, na.rm = T)
    
    ## Drop everything above the maximum sampling distance and multiply trees above radius_max with their new area
    # filter( !distance > radius_max_B_cm ) %>% ## this might drop plots, therefore:
    mutate(count == if_else(distance > radius_max_B_cm, 0, count)) %>% ## this way, the trees get dropped in all calculations and averages!
    
    ## Rename
    mutate(count_ha = count) %>% # count is already per hectare, countarea == 1
    mutate(countarea = 1) %>%
    
    ### Calculate the area per tree
    # dplyr::mutate(countarea_tree = tidyr::replace_na(pi * (25*dbh)^2 * 1e-10, 0)) %>% ## [ha == 1e10 mm2] ## dbh has NAs for count 0
    dplyr::mutate(countarea_tree = pi * 25^2 * dbh^2 * 1e-10) %>% ## [ha == 1e10 mm2] ## dbh has NAs for count 0
      ## equivalent to pi * r^2

    ## Assign a new area to the trees sampled on r_max (truncated area) and correct their count_ha
    mutate(dbhisgreatermax = dbh > dbh_max_B) %>%
    mutate(count_ha = if_else(dbhisgreatermax, count_ha * (countarea_tree/area_max_B), count_ha)) %>% ## count_ha_new = count_ha_old * area_old/area_new
    mutate(countarea_tree = if_else(dbhisgreatermax, area_max_B, countarea_tree)) %>% ## [ha == 1e10 mm2] ## dbh has NAs for count 0
    
    ## calculate the count per tree after the correction of the count/ha
      ## this is equivalent to doing it before, because both the count_ha and countarea_tree were corrected
    dplyr::mutate(count_tree = count_ha * countarea_tree) %>%
      ## The distribution looks fine. Corresponds to the crazy distribution of zf == Zählfaktoren.
    dplyr::mutate(f = floor2(count_tree)/count_tree) %>% ## There are a lot of trees that have count 0.999999826790235 for some reason.
    dplyr::mutate(count_tree = f * count_tree, countarea_tree = f * countarea_tree) %>%
    
    ### Calculate basal area per hectare
    mutate(ba = pi * (dbh/2)^2 * 1e-6) %>% # mm^2 to m^2)
    mutate(ba_ha = ba * count_ha) %>%

    ### make stages
    # quantile(B$dbh, seq(0, 1, by = 1e-1), na.rm = T): 160 is the 10%tile, 206 is the 20%tile
    ## 100mm…200mm == a
    mutate(stage = as.factor(if_else(dbh < threshold_dbh, "A", "B"))) %>%

    ## summarize per all ids, ... "plotobsid", "tax", "stage", being the fastest-varying
    ## interesting stuff like age, height, will get dropped here
    group_by_at(c(id_select_B, "stage")) %>%
    dplyr::summarize(area_obs = weighted.mean(countarea_tree, w = count_ha, na.rm = TRUE), ## Counts with zeroes are automatically removed from the area here.
                     count_ha = sum(count_ha, na.rm = T),
                     ba_ha = sum(ba_ha, na.rm = T),
                     ba_obs = ba_ha * area_obs, ## ba_obs = sum(ba_ha, na.rm = T) * area_obs ##
                     count_obs = sum(count_tree, na.rm = TRUE),
                     offset_ba_ha = count_obs/ba_ha ## [ha/m2] the offset includes both area and conversion to basal area.
                        ##Equivalent to: offset_ba_ha == (count_obs/ba_obs) * area [ha/m2]
                     ) %>%
    ## this will yield NaN for stages which have not been observed
    ungroup() %>% #!!!
    
    ## recovery check:
    ## plot(I(count_obs/area_obs) ~ count_ha, data = B_temp[B$stage == "A",])
    ## plot(I(count_obs/area_obs) ~ count_ha, data = B[B$stage == "B",])
    ## plot(I(ba_obs/area_obs) ~ ba_ha, data = B[B$stage == "A",]); abline(0,1)
    ## plot(count_obs ~ I(ba_ha * offset_ba_ha), data = B[B$stage == "B",]); abline(0,1)
    
    ## taxa also have to be completed to include all stages, but only within obsid!!!
    ## (to avoid that 1. absence observations are fabricated where no observations where made, and 2. include observations of species that were not looked for in a survey)
    ## Note: method is a level that is maskedly reflecting stage, thus will prevent completion
    mutate(accountfortaxamethods = obsid) %>%
    tidyr::complete(nesting(accountfortaxamethods, taxid, tax), stage, nesting(obsid, plotid, plotobsid, methodid, time),
                    fill = list(count_ha = 0, count_obs = 0, ba_ha = 0, ba_obs = 0)) %>%
    filter(accountfortaxamethods == obsid) %>%
    dplyr::select(-accountfortaxamethods) %>%
    
    ## Offset after completion with zeroes.
    group_by(methodid, obsid, tax, taxid, stage) %>%
    mutate(count_ba_survey = mean(count_ha/ba_ha, na.rm = T)) %>%
    ## Average area per tax, stage, methodid
    dplyr::mutate(area_0_avg = weighted.mean(area_obs, w = count_ha, na.rm = TRUE),
                  area_0_q1 = DescTools::Quantile(area_obs, weights = count_ha, probs = 0.25, na.rm = T, names = F),
                  area_0_q3 = DescTools::Quantile(area_obs, weights = count_ha, probs = 0.75, na.rm = T, names = F)) %>%
    ungroup() %>%
    
    ## Distribute the 0-areas
    ##    1. Take the mean areas from the respective other species on the plot (for two species this will just be the first not NA)
    ##    2. When neither species is present, take the integral
    mutate(iszero = count_ha == 0) %>%
    group_by(plotobsid, stage) %>%
    mutate(n_area = n_distinct(area_obs, na.rm = T), area_obs_avg_plot = mean(area_obs, na.rm = T)) %>% # for two species this will just be first(area_obs[!is.na(area_obs)])
    ungroup() %>%
    mutate(area_0 = case_when(iszero & stage == "A" & n_area > 0 ~ area_obs_avg_plot,
                              iszero & stage == "B" & n_area > 0 ~ area_obs_avg_plot,
                              iszero & stage == "A" & n_area == 0 ~ area_0_A,
                              iszero & stage == "B" & n_area == 0 ~ area_0_B)) %>%
    
    mutate(offset = case_when(iszero & stage == "A" ~ area_0,
                              iszero & stage == "B" ~ area_0 * count_ba_survey,
                              !iszero & (stage == "A") ~ area_obs,
                              !iszero & (stage == "B") ~ offset_ba_ha)) %>%
    
    mutate(offset_avg = case_when(!iszero ~ offset,
                                  iszero & stage == "A" ~ area_0_avg,
                                  iszero & stage == "B" ~ area_0_avg * count_ba_survey)) %>%
                                  
    mutate(offset_q1 = case_when(!iszero ~ offset,
                                 iszero & stage == "A" ~ area_0_q1,
                                 iszero & stage == "B" ~ area_0_q1 * count_ba_survey)) %>%
    
    mutate(offset_q3 = case_when(!iszero ~ offset,
                                 iszero & stage == "A" ~ area_0_q3,
                                 iszero & stage == "B" ~ area_0_q3 * count_ba_survey)) %>%
    
    droplevels() %>%
    ## stage creation will yield NA stages for when dbh is NA (for completed species)
    ## but any(is.na(B$dbh) & (B$count != 0)) # so that, after completion:
    tidyr::drop_na(stage)
  
  return(B)
}


## prepareSmallData --------------------------------
# J  <- tar_read("Data_small")
# taxon_select <- tar_read("taxon_select")
# regclass_select <- tar_read("regclass_select")
prepareSmallData <- function(J,
                           taxon_select,
                           regclass_select,
                           id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time")
                           ) {
  
  id_select_S <- intersect(names(J), id_select)
  regclass_smaller <- c("h[20,50)", "h[50,130)") ## The smaller size classes, that have been sampled on more or less consistent areas
  
  selectTaxa <- function(Abundance, taxon_select) {
    Abundance %>%
      mutate(tax = forcats::fct_other(tax, !!!as.list(taxon_select), other_level = "other")) %>%
      mutate(taxid = if_else(tax == "other", "other", as.character(taxid))) %>%
      droplevels()
  }
  
  
  ## J is already completed!
  J %<>%
    ### Lump taxa first, to facilitate summarize() and complete() later
    selectTaxa(taxon_select) %>%
    
    ### subset to time-constant regclasses, NOTE: There are more here.
    ## levels(J$regclass)
    ## table(J$regclass,J$obsid)
    dplyr::filter(as.integer(regclass) %in% regclass_select)  %>% ## these classes are all present in all three obsids, but consider different amounts of plots: # %>% dplyr::group_by(regclass, obsid) %>% dplyr::summarize(count_ha = sum(count/countarea))
    dplyr::mutate(issmallerregclass = as.character(regclass) %in% regclass_smaller) %>% ## boolean for selection of smaller plots in DE_BWI1987

    dplyr::group_by(methodid) %>%
    dplyr::mutate(area_0_methodid = mean(unique(countarea)), na.rm = T) %>%
    dplyr::ungroup() %>%
    
    ## drop size classes and summarize counts
    dplyr::group_by_at(c(id_select_S, "area_0_methodid")) %>%
    dplyr::summarize(count_ha = sum(count/countarea, na.rm = TRUE),
                     area_obs = weighted.mean(countarea, count/countarea, na.rm = TRUE), ## this produces NaNs when count is 0
                     count_obs = sum(count, na.rm = TRUE), ## The average weighted by the count of trees within this size class.
                     anysmallerregclass = any(issmallerregclass, na.rm = T)
                     ) %>% 
    dplyr::ungroup() %>%
    
    ## Distribute the 0-areas
    ##    1. Take the mean areas from the respective other species on the plot
    ##    2. When neither species is present, take the integral
    mutate(iszero = count_ha == 0) %>%
    group_by(plotobsid) %>%
    mutate(n_area = n_distinct(area_obs, na.rm = T), area_0_plot = mean(area_obs, na.rm = T)) %>% # for two species this will just be first(area_obs[!is.na(area_obs)])
    ungroup() %>%
    mutate(area_0 = case_when(iszero & n_area > 0 ~ area_0_plot,
                              iszero & n_area == 0 ~ area_0_methodid)) %>% ## !any(J$n_area == 2 & J$count_ha == 0)
    mutate(offset = if_else(iszero, area_0, area_obs)) %>%
    mutate(offset_avg = offset, offset_q1 = offset, offset_q3 = offset) %>% ## for compatibility with different offsets in big trees

    dplyr::mutate(stage = factor("J")) %>%
    droplevels()
    
    ## recovery test of count/ha
    ## plot(count_obs/area_obs ~ count_ha, data = J)

  return(J)
}



## joinStages --------------------------------
# B  <- tar_read("Data_big_area")
# J  <- tar_read("Data_small_area")
# taxon_select <- tar_read("taxon_select")
# threshold_dbh <- tar_read("threshold_dbh")
# dir_publish  <- tar_read("dir_publish")
joinStages <- function(B, J,
                       taxon_select,
                       threshold_dbh,
                       id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time") #, tablepath = dir_publish
                       ) {

  id_select_B <- intersect(names(B), id_select) %>% setdiff("treeid") ## make sure to always exclude treeid for good grouping
    ## B doesn't include clusterid!
  
  BA <- B %>%
    group_by_at(id_select_B) %>% ## grouping without stage -> summarizes across stages within plot and tax
    dplyr::summarize(area_obs = weighted.mean(area_obs, w = count_ha, na.rm = T),
                     count_ha = sum(count_ha, na.rm = T),
                     count_obs = sum(count_obs, na.rm = T),
                     ba_ha = sum(ba_ha, na.rm = T),
                     ba_obs = sum(ba_obs, na.rm = T)
                     ) %>%
    mutate(stage = "BA")
    ## recovery check
    ## plot(I(count_obs/area_obs) ~ count_ha, data = BA)
    ## plot(I(offset_ba_ha * ba_ha) ~ count_obs, data = BA); abline(0,1)
  
  Stages <- dplyr::bind_rows(J, B, BA)
  Stages[is.na(Stages)] <- NA # replace NaNs (is.na returns TRUE for NaNs!)

  ## Print the averages
  # Print <- dplyr::select(Stages, stage, methodid, tax, area_avg_method) %>%
  #   group_by(stage, tax, methodid) %>%
  #   summarize(area_avg = first(area_avg_method)) %>% # == unique() %>%
  #   ungroup() %>%
  #   filter(!is.na(stage)) %>%
  #   arrange(stage, methodid)
  # 
  # write.csv(Print, file.path(tablepath, "Areas_average_method.csv"))
  # cat("The calculated average areas in ha for the methods are …")
  # print(as.data.frame(Print))
  
  return(Stages)
}



## fitS --------------------------------
# BA_s <- tar_read("BA_s")
# BA_s <- tar_read("seedlings_s")
# BA <- BA_s[[1]]
# path  <- tar_read("dir_publish")
fitS <- function(BA, path = NULL) {
  
  BA_coordinates <- cbind(BA, st_coordinates(BA))
  tax <- attr(BA, "taxon")
  attr(BA_coordinates, "taxon") <- tax ## cbind drops attr somehow
  
  ## Get measures for
  # box <- st_bbox(BA) %>%
  #   st_as_sfc() %>%
  #   st_as_sf() %>%
  #   st_cast(to = "POINT")
  # 
  # dist_lon <- st_distance(box[3:4,])[1,2] %>% # distance between the two points at the top of the box in meters ()
  #   as.numeric()
  # dist_lat <- st_distance(box[c(1, 4),])[1,2] %>% # distance between the two points at the left of the box
  #   as.numeric()
  # 
  # 
  # Grid <- st_bbox(BA) %>%
  #   st_make_grid(n = c(100, 200), what = "centers") %>%
  #   st_coordinates()
  # 
  # k_Y <- Grid[,"Y"] %>% unique()
  # k_X <- Grid[,"X"] %>% unique()
  # 
  # n_k_lon <- round(dist_lon / 8000)
  # n_k_lat <- round(dist_lat / 8000)
  
  ### mgcv
  fit <- gam(ba_ha ~ s(Y, X, bs = "sos", k = 600), family = nb, data = BA_coordinates) ## The first argument is taken to be latitude (in degrees) and the second longitude (in degrees).
  
  if(!is.null(path)) {
    s <- summary(fit)
    textext <- itsadug::gamtabs(s, caption = "Summary of the thin plate spline fit for the background basal area ...", label = paste0("tab:gam_", tax))
    cat(textext, file = file.path(path, paste0(tax, "_summary_gam.tex")), fill = T) %>% invisible()
  }
  
  attr(fit, "taxon") <- tax

  return(fit)
}


## predictS --------------------------------
# Stages_env  <- tar_read("Stages_env")
# fits  <- tar_read("fits_s"); fit <- fits[[1]]
predictS <- function(fits, Stages_env) {
  
  pred <- function(fit) {
    
    predicted <- predict(fit, newdata = X_unique, type = "link")

    Predicted <- as.data.frame(predicted) ## one column
    Predicted <- Predicted[X_df$id, , drop = F] ## through order, this has the right nrows
    
    colnames(Predicted) <- paste0("s_", attr(fit, "taxon"))
    return(Predicted)
  }
  
  X_df <- st_coordinates(Stages_env) %>%
    as.data.frame() %>%
    mutate(id = vctrs::vec_group_id(.))
  
  X_unique <- filter(X_df, !duplicated(id)) %>%
    dplyr::select(1:2)
  
  predictions <- lapply(fits, pred) ## list[N_taxa] of data_frames
  
  S <- bind_cols(Stages_env, predictions)

  return(S) # plot(S["s_Fagus.sylvatica...24"])
}


## predictSurfaces --------------------------------
# fits  <- tar_read("fits_s")
predictSurfaces <- function(fits) {

  Ger <- raster::getData("GADM", country = "DE", level = 0, path = "Data/")
  R <- raster::raster(raster::extent(bbox(Ger)), resolution = c(0.01, 0.01))
  crs(R) <- crs(Ger)
  Coords <- rasterToPoints(R, spatial = TRUE) %>%
    as.data.frame()
  R$X <- Coords$x
  R$Y <- Coords$y
  R <- raster::mask(R, Ger) # plot(R)
  
  # P <- raster::predict(R, fit); plot(P, col = viridis::viridis(255))
  surfaces <- lapply(fits, function(f) raster::predict(R, f, type = "response"))
  names(surfaces) <- sapply(fits, function(f) attr(f, "taxon"))
  
  return(surfaces)
}

## plotSurfaces ------------------------------------
# surfaces <- tar_read(surfaces_Seedlings_s)
# path  <- tar_read("dir_publish")
plotSurfaces <- function(surfaces, path, themefun) {
  
  plotRaster <- function(r) {
    require(rasterVis)
    
    gplot(r) +
      geom_tile(aes(fill = value)) +
      labs(fill = "Basal area [m^2 ha^-1]") +
      scale_fill_viridis_c(na.value = "transparent") +
      coord_quickmap() + ## for coord_map: projection = "lambert", lat1 = 48, lat2 = 54
      theme_map(themefun) +
      
      
      ## Fix at some time
      # annotation_scale(width_hint = 0.3,
      #                 style = "bar",
      #                 location = "br",
      #                 bar_cols = c("black", "white"),
      #                 line_col = "black",
      #                 text_col = "black",
      #                 text_cex = 1.1,
      #                 pad_x = unit(0.3, "cm")) +
      annotation_north_arrow(which_north = "true",
                            location = "tl",
                            pad_x = unit(0.3, "cm"),
                            style = north_arrow_minimal(line_col = "black", text_col = "black", fill = c("black"), text_size = 17))
  }

  surfaceplots <- lapply(surfaces, function(s) plotRaster(s))
  
  mapply(function(p, n) ggsave(paste0(path, "/", "Surface_", n, ".png"), p, device = "png", width = 8, height = 9), surfaceplots, names(surfaceplots))
  
  return(surfaceplots)
}


saveStages_s <- function(Stages_s) {
  
  path <- "Data/Stages_s.rds"
  if(!dir.exists("Data")) dir.create("Data")
  saveRDS(Stages_s, file = path)
  return(path)
}



## selectLocs --------------------------------
# Stages_s <- tar_read("Stages_s")
# predictor_select <- tar_read("predictor_select")
# loclevel <- tar_read("loc")
selectLocs <- function(Stages_s, predictor_select, selectpred = F,
                       id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time"),
                       loc = c("plot", "nested", "cluster"),
                       n_locations = 1000
                       ) {

  loclevel <- match.arg(loc)
  
  message(paste(Stages_s %>% pull(clusterid) %>% unique() %>% length(), "clusters, and",
                Stages_s %>% pull(plotid) %>% unique() %>% length(), "plots before selectLocs()."))

  ## for reference
  disturbance_select = c("standage_DE_BWI_1",
                         "allNatRegen_DE_BWI_2", "allNatRegen_DE_BWI_3", "anyUnnatRegen_DE_BWI_3", "anyUnnatRegen_DE_BWI_2",
                         "anyHarvested_DE_BWI_2", "anyHarvested_DE_BWI_3",
                         "anyForestryDamage_DE_BWI_2", "anyForestryDamage_DE_BWI_3")
  
  # Stages_s[predictor_select] %>% is.na() %>% colSums()
  ## NOTE: waterLevel_loc, phCaCl_esdacc are data on the plot level
  
  ## Filter based on environment
  if (selectpred) {
    Stages_s  %<>%
      filter_at(predictor_select, function(x) !is.na(x))
  }
  
  ## Filter based on management etc.
  Stages_select <- Stages_s %>%
    
    ## The any variables have been assembled per plotid!
    ## Plots were excluded that had any record of unnatural regeneration in DE_BWI_2 or 3
    mutate(anyUnnatRegen_DE_BWI_2 = tidyr::replace_na(anyUnnatRegen_DE_BWI_2, FALSE),
           anyUnnatRegen_DE_BWI_3 = tidyr::replace_na(anyUnnatRegen_DE_BWI_3, FALSE)) %>%
    filter( !(anyUnnatRegen_DE_BWI_2 | anyUnnatRegen_DE_BWI_2)) %>%
    
    ## Anything harvested
    filter( !(anyHarvested_DE_BWI_2 | anyHarvested_DE_BWI_3)) %>%
    
    ## Any damage through forestry
    filter( !(anyForestryDamage_DE_BWI_2 | anyForestryDamage_DE_BWI_3))
    
    # unique(Stages_select$clusterid) %>% length() ## 5236
    
  
  ## Filter plots based on tree observations
  Stages_select %<>%
    
    group_by(plotid, obsid) %>%
    
    ## subset to plots without any clear cut, i.e. only zero observations or NA after there had been observations
    mutate(isclear = sum(count_ha, na.rm = T) == 0) %>% ## sum will replace vectors with exclusively NAs with 0
    
    group_by(plotid) %>%
    
    ## subset to plots that have any observation in a smaller regclass (height 20--130cm)
    mutate(anysmallerregclass1987 = any(anysmallerregclass & obsid == "DE_BWI_1987")) %>%
    
    mutate(isclearcut_2002 = any( (!isclear[obsid == "DE_BWI_1987"]) & isclear[obsid == "DE_BWI_2002"]),
           isclearcut_2012 = any( (!isclear[obsid == "DE_BWI_2002"]) & isclear[obsid == "DE_BWI_2012"]),
           isclearcut = isclearcut_2002 | isclearcut_2012) %>%
    filter((!isclearcut) & anysmallerregclass1987)
      ## 5012 plots remaining; # Stages_select %>% filter(isclearcut) %>% pull(plotid) %>% unique()
    
  
  if (loclevel %in% c("nested", "cluster")) {
    
    ## Filter clusters based on arbitrary thresholds that ensure sufficient samples
    Stages_select %<>% 
      group_by(clusterid) %>%
      
      ## subset to clusters with at least three surveys
      mutate(n_surveys = n_distinct(obsid)) %>%
      filter(n_surveys >= 3) %>% # table(Stages_select$n_surveys) ## 2: 53626, 3: 119998
      dplyr::select(-n_surveys) %>%
      
      ## subset to clusters with at least two plots
      mutate(n_plots = n_distinct(plotid)) %>%
      filter(n_plots > 2) %>%
      dplyr::select(-n_plots) %>%
      
      ## get clusters with at least some of both in small trees
      mutate(anyFagus = any(count_ha > 0 & tax == "Fagus.sylvatica")) %>%
      mutate(anyOther = any(count_ha > 0 & tax == "other")) %>%
      mutate(anySmallFagus = any(count_ha > 0 & tax == "Fagus.sylvatica" & stage == "J")) %>%
      mutate(anySmallOther = any(count_ha > 0 & tax == "other" & stage == "J")) %>%
      mutate(anyBigFagus = any(count_ha > 0 & tax == "Fagus.sylvatica" & stage %in% c("A", "B"))) %>%
      mutate(anyBigOther = any(count_ha > 0 & tax == "other" & stage %in% c("A", "B"))) %>%
      ungroup()
    
    ## Confined to clusters with any observation of the taxa in defined sizeclasses
    Stages_select %<>%
      filter(anyFagus & anyOther) # %>% pull(clusterid) %>% unique() %>% length() ## 635
  
  } else { ## case loclevel == "plot"
    
    Stages_select %<>% 
      group_by(plotid) %>%
      
      ## subset to clusters with at least three surveys
      mutate(n_surveys = n_distinct(obsid)) %>%
      filter(n_surveys >= 3) %>% # table(Stages_select$n_surveys) ## 2: 53626, 3: 119998
      dplyr::select(-n_surveys) %>%
      
      ## get plots with at least some of both in small trees
      mutate(anyFagus = any(count_ha > 0 & tax == "Fagus.sylvatica")) %>%
      mutate(anyOther = any(count_ha > 0 & tax == "other")) %>%
      mutate(anySmallFagus = any(count_ha > 0 & tax == "Fagus.sylvatica" & stage == "J")) %>%
      mutate(anySmallOther = any(count_ha > 0 & tax == "other" & stage == "J")) %>%
      mutate(anyBigFagus = any(count_ha > 0 & tax == "Fagus.sylvatica" & stage %in% c("A", "B"))) %>%
      mutate(anyBigOther = any(count_ha > 0 & tax == "other" & stage %in% c("A", "B"))) %>%
      ungroup()
      
    
    ## Confined to plots with any observation of the taxa in defined sizeclasses
    Stages_select %<>%
      filter(anyFagus & anyOther) # %>% pull(clusterid) %>% unique() %>% length() ## 3468
    
    ## Select only one random plot per cluster
    Stages_select %<>%
      group_by(clusterid) %>%
      mutate(plotid_select = sample(unique(plotid), 1, replace = F)) %>%
      filter(plotid == plotid_select)
    
    ## Subset n_plots
    plotid_subset <- Stages_select$plotid %>% unique() %>% sample(n_locations, replace = FALSE)
    Stages_select %<>%
      filter(plotid %in% plotid_subset)
      
  }
  
  
  # ## Selecting an equal no. of plots with and without Fagus
  # Stages_Fagus <- Stages_select %>%
  #   filter(anyFagus)
  # 
  # n_Fagus <- length(unique(Stages_Fagus$clusterid)) # 95
  # 
  # Stages_other <- Stages_select %>%
  #   filter(!anyFagus) %>%
  #   filter(clusterid %in% base::sample(unique(.$clusterid), n_Fagus))
  # 
  # Stages_select <- bind_rows(Stages_other, Stages_Fagus)
  # # unique(Stages_select$clusterid) %>% length() ## 190
  
  
  Stages_select %<>%
    dplyr::select(-any_of(setdiff(disturbance_select, "standage_DE_BWI_1")))
  
  message(paste(Stages_select %>% pull(clusterid) %>% unique() %>% length(), "clusters, and",
                Stages_select %>% pull(plotid) %>% unique() %>% length(), "plots after selectLocs()."))
  
  if(anyNA(Stages_select$time)) {
    
    n_na <- Stages_select %>%
      filter(is.na(time)) %>%
      pull(clusterid) %>%
      n_distinct()
    
    Stages_select %<>%
      filter(!is.na(time))

    warning("selectLocs(): There were ", n_na, " clusters with missing variable `time`. These clusters have been dropped.")
  }
  
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



## countTransitions --------------------------------
# Data_big  <- tar_read("Data_big")
# Data_big_status <- tar_read("Data_big_status")
# Env_cluster <- tar_read(Env_cluster)
# Stages_select <- tar_read("Stages_select")
# taxon_select <- tar_read("taxon_select")
# threshold_dbh <- tar_read("threshold_dbh")
# radius_max <- tar_read("radius_max")

countTransitions <- function(Data_big, Data_big_status, Env_cluster, Stages_select, taxon_select, threshold_dbh, radius_max) {
  
  ### local functions
  selectTaxa <- function(Abundance, taxon_select) {
    Abundance %>%
      mutate(tax = forcats::fct_other(tax, !!!as.list(taxon_select), other_level = "other")) %>%
      mutate(taxid = if_else(tax == "other", "other", as.character(taxid))) %>%
      droplevels()
  }
  
  diffTime <- function(time_first, time_second) {
    timediff <- lubridate::year(time_second) - lubridate::year(time_first)
    return(timediff)
  }
  
  findPrecedingTime <- function(time) {
    timeseq <- sort(unique(time))
    i <- match(time, timeseq)-1
    i[i == 0] <- NA
    time <- timeseq[i]
    return(time)
  }
  
  matchSucceedingTime <- function(time) {
    timeseq <- sort(unique(time))
    i <- match(time, timeseq)+1
    i[i > length(timeseq)] <- NA
    return(i)
  }
  
  integrateCounts <- function(count_first, count_second, timediff) {
    
    ## handle case: no timediff because the prior date does not exist
    if(is.na(timediff)) {
      
      return(as.numeric(NA))
      
    } else {
      time <- 0:(timediff-1)
      slope <- (count_second - count_first) / timediff
      count <- time * slope + count_first
      return(sum(count))
    }
  }
  
  # followedAbyB <- function(stage) {
  #   i <- zoo::rollapply(as.character(stage), width = 2, identical, c("A", "B"))
  #   return(c(i, F)) # padding in the end for vector length
  # }
  
  countNew <- function(idB_first, idB_second, count_trans) {
    idB_first <- c(NA, idB_first) ## include NAs so they can be matched and NOT counted
    count <- sum(!(idB_second %in% idB_first), na.rm = T)
    return(count * count_trans)
  }
  
  ## Make explicit that I drop all values where there is a 0 in the base population. 0 in the transitioning pop get 0 results
  divide <- function(trans, base) {
    prop <- case_when(
      trans == 0 ~ 0,
      is.na(trans) ~ as.numeric(NA), ## also includes NaNs
      base == 0 | is.na(base) ~ as.numeric(NA),
      TRUE ~ trans/base
    )
    return(prop)
  }

  #### Prepare joining data sets
  S <- Data_big_status %>%
    dplyr::filter(obsid %in% c("DE_BWI_2002", "DE_BWI_2012")) %>% # !!!
    dplyr::rename(count = count_ha) %>% # just temporarily for consistency with Data_big
    dplyr::mutate(excluded = dead | harvested) %>%
    dplyr::mutate(treejoinid = interaction(obsid, treeid))
  
  
  E <- Env_cluster[c('plotid', 'clusterid')] %>%
    st_drop_geometry()
  
  clusterid_select <- unique(Stages_select$clusterid)
  
  ## Areas and radii
  radius_max_B_cm <- radius_max/10 ## conversion from mm to cm
  dbh_max_B <- radius_max_B_cm/25 * 10 ## [mm]
  area_max_B <- pi * radius_max_B_cm^2 * 1e-8 # [ha]
  
  ##### Stages with plot level A2B

  ## the count factor at stage transition
  countfactor_A2B <- 10^10/(pi * 25^2 * threshold_dbh^2)
  kluppschwelle <- 100 # [mm]
  countfactor_J2A <- 10^10/(pi * 25^2 * kluppschwelle^2)
  
  ## For count_ha formula:
  # see: http://wiki.awf.forst.uni-goettingen.de/wiki/index.php/Bitterlich_sampling#Estimation_of_number_of_stems
  # with Zählfaktor k = 4 and c = 50/sqrt(k): c == 50/sqrt(4) == 25
  # expansion factor n_ha = 10000 * 1000^2/(pi * c^2 * dbh^2) [factor 1000^2 for dbh in mm instead of m]
  
  ## the threshold distance (plot radius) at stage transition
  ## i.e. the radius within which trees in A are counted
  distance_threshold <- 25 * threshold_dbh/10 # [mm] to [cm conversion] with 1/10
  distance_kluppschwelle <- 25 * kluppschwelle/10
  
  area_A2B <- pi * distance_threshold^2 * 1e-8 ## [cm2 to ha]
  area_J2A <- pi * distance_kluppschwelle^2 * 1e-8 ## [cm2 to ha]
  
  ## For multidply parallelization
  # cluster <- multidplyr::new_cluster(3)
  # cluster_copy(cluster, c("countNew"))
  
  J <- Stages_select %>%
    st_drop_geometry() %>%
    mutate(obsid_1987 = obsid == "DE_BWI_1987", obsid_2002 = obsid == "DE_BWI_2002", obsid_2012 = obsid == "DE_BWI_2012") %>%
    mutate(year = lubridate::year(time)) %>%
    
    group_by(plotid) %>%
    mutate(yearbefore_plot = findPrecedingTime(year), timediff_plot = year - yearbefore_plot) %>% # table(Stages_transitions$timediff_plot, Stages_transitions$year)
    
    ## Find a state interpolation between two surveys by integration
    group_by(plotid, taxid) %>%
    filter(stage == "J") %>%
    mutate(count_J_sum_before = case_when(obsid_2002 ~ sum(count_ha[first(yearbefore_plot[obsid_2002]) == year & stage == "A"], na.rm = T),
                                          obsid_2012 ~ sum(count_ha[first(yearbefore_plot[obsid_2012]) == year & stage == "A"], na.rm = T))) %>%
    group_by(plotid, taxid, obsid) %>%
    mutate(count_J_integr_plot = integrateCounts(first(count_J_sum_before), sum(count_ha, na.rm = T), first(timediff_plot))) %>%
    ungroup()
    
    
  
  Stages_transitions <- Data_big %>%
    filter(!is.na(dbh)) %>% # drop non-sample/completed trees
    
    ## add radii of all trees via matching
    mutate(treejoinid = interaction(obsid, treeid)) %>%
    bind_cols(distance = S$distance[match(.$treejoinid, S$treejoinid)]) %>%
    ## all trees in 2002, and 2012 have a radius
    ## Stages_A2B$distance[Stages_A2B$obsid != "DE_BWI_1987"] %>% anyNA
    
    ## add dead/excluded trees
    # if necessary, plots with harvest will be excluded later.
    bind_rows(dplyr::select(filter(S, excluded), -clusterid)) %>%
    
    ## Drop everything above the maximum sampling distance, S also has distance!
    filter( !distance > radius_max_B_cm ) %>%
    
    ## add clusters via matching
    bind_cols(clusterid = E$clusterid[match(.$plotid, E$plotid)]) %>%
    
    ## select clusters
    filter(clusterid %in% clusterid_select) %>%
    
    ### All trees on a plot have an id!
    # group_by(plotid) %>%
    # mutate(allid = all(!is.na(treeid)))
    
    ### Lump taxa first, to facilitate summarize() and complete() later
    selectTaxa(taxon_select) %>%
    
    ### Handle counts. (Correct area truncation. See comments in prepareBigData(...) for details.)
    rename(count_ha = count) %>% # count is already per hectare, countarea == 1
    mutate(countarea_tree = pi * 25^2 * dbh^2 * 1e-10) %>% ## [ha == 1e10 mm2] ## Calculate the area per tree, equivalent to pi * r^2
    mutate(count_ha = if_else(dbh > dbh_max_B, count_ha * (countarea_tree/area_max_B), count_ha)) %>% ## ## Assign a new area to the trees sampled on r_max (truncated area) and correct their count_ha: count_ha_new = count_ha_old * area_old/area_new
    
    ### make stages, including dead etc. trees in B!
    mutate(stage = as.factor(if_else(dbh < threshold_dbh, "A", "B"))) %>%
    
    ### time differences.
    ## times are equal among clusterids, but different plots have different no of observations
    mutate(year = lubridate::year(time)) %>%
    
    group_by(plotid) %>%
    mutate(yearbefore_plot = findPrecedingTime(year), timediff_plot = year - yearbefore_plot) %>% # table(Stages_transitions$timediff_plot, Stages_transitions$year)
    
    group_by(clusterid) %>%
    ## some clusters have been more times than singular plots on them: table(Stages_transitions$timediff_cluster, Stages_transitions$timediff_plot)
    mutate(yearbefore_cluster = findPrecedingTime(year), timediff_cluster = year - yearbefore_cluster) %>% # 
    ungroup()
  
  
  Stages_transitions %<>%
    
    ## Count new trees in B, that are within the sampling radius of A
    # (preparing some booleans to improve speed? of the groupwise counting)
    mutate(obsid_1987 = obsid == "DE_BWI_1987", obsid_2002 = obsid == "DE_BWI_2002", obsid_2012 = obsid == "DE_BWI_2012") %>%
    group_by(plotid, taxid) %>%
    # multidplyr::partition(cluster) %>%
    
    ## Introduce new treeids that have NAs, depending on whether they are excluded based on them being outside the virtual fixed plot
    mutate(treeid_int = as.integer(factor(treeid))) %>% ## integer %in% seems to be much faster, no NAs in treeid
    mutate(treeid_B = replace(treeid_int, stage != "B", NA)) %>% ## integer %in% seems to be much faster, no NAs in treeid
    mutate(treeid_A = replace(treeid_int, stage != "A", NA)) %>%
    
    mutate(treeid_A2B = replace(treeid_B, distance > distance_threshold, NA)) %>% ## introducing a lot of NAs: table(is.na(Stages_A2B$treeid_A2B))
    mutate(treeid_J2A = replace(treeid_A, distance > distance_kluppschwelle, NA)) %>% ## introducing a lot of NAs: table(is.na(Stages_A2B$treeid_A2B))
      
    ## These treeids will only be counted as newly appearing, when they hadnt been in the prior survey AND NOT NA!
    ## countNew(): treeid_A2B/J2A (including NAs for trees outside the radius for A, minimal), will be matched in treeid_int of the former survey
    ## Some plots have the observations 1987 and 2002.
    
    mutate(count_A2B_plot = case_when(obsid_2002 ~ countNew(treeid_A2B[obsid_1987], treeid_A2B[obsid_2002], countfactor_A2B),
                                      obsid_2012 ~ countNew(treeid_A2B[obsid_2002 | obsid_1987], treeid_A2B[obsid_2012], countfactor_A2B))) %>%

    mutate(count_J2A_plot = case_when(obsid_2002 ~ countNew(treeid_J2A[obsid_1987], treeid_J2A[obsid_2002], countfactor_J2A),
                                      obsid_2012 ~ countNew(treeid_J2A[obsid_2002 | obsid_1987], treeid_J2A[obsid_2012], countfactor_J2A))) %>%
    
    ## Find the correct number of trees that have been actually observed and corresponding area offsets
    mutate(area_A2B = area_A2B) %>%
    mutate(area_J2A = area_J2A) %>%
    
    mutate(count_A2B_plot_obs = as.integer(count_A2B_plot * area_A2B)) %>% ## [1/ha * ha]
    mutate(count_J2A_plot_obs = as.integer(count_J2A_plot * area_J2A)) %>% ## [1/ha * ha]
    
    ## Find a state interpolation between two surveys by integration. Unit [1/ha]
    mutate(count_A_sum_before = case_when(obsid_2002 ~ sum(count_ha[first(yearbefore_plot[obsid_2002]) == year & stage == "A"], na.rm = T),
                                          obsid_2012 ~ sum(count_ha[first(yearbefore_plot[obsid_2012]) == year & stage == "A"], na.rm = T))) %>%
    
    left_join(J[c("count_J_integr_plot", "count_J_sum_before", "plotid", "taxid", "obsid")], by = c("plotid", "taxid", "obsid")) %>%
    
    group_by(plotid, taxid, obsid) %>%
    mutate(count_A_integr_plot = integrateCounts(first(count_A_sum_before), sum(count_ha[stage == "A"], na.rm = T), first(timediff_plot))) %>%
    mutate(h_plot = divide(count_A2B_plot, count_A_integr_plot)) %>%
      ## NAs come from NAs in count_A_integr_plot, where if(is.na(timediff)) NA, i.e. no observation before.
      ## NaNs come from 0 in count_A_integr_plot
    
    mutate(g_plot = divide(count_J2A_plot, count_J_integr_plot)) %>%
    
    ungroup()
    

  
  
  ### calculate h
  
  # Changes_A <- filter(Data_big) %>% ## use only living trees in A
  #   filter(dbh < threshold_dbh & !is.na(dbh)) %>%
  #   selectTaxa(taxon_select) %>%
  #   group_by(plotid, taxid) %>%
  #   summarize(count_ha_2002_plot = sum(count[obsid == "DE_BWI_2002"], na.rm = T))
  
  # Stages_A2B %<>% 
  # bind_rows(Changes_A) %>%
  # group_by(plotid, taxid) %>%
  # mutate(h_plot_2012 = first(count_A2B_2012_plot[!is.na(count_A2B_2012_plot)]) / first(count_ha_2002_plot[!is.na(count_ha_2002_plot)])/ timediff_2012) %>%
  # filter(stage != "A")
  
  return(Stages_transitions)
}


## summarizeTaxa --------------------------------
# Stages_select  <- tar_read("Stages_select")
# Seedlings_s <- tar_read("Seedlings_s")
# Data_seedlings <- tar_read(Data_seedlings)
# Data_big <- tar_read("Data_big")

summarizeTaxa <- function(Data_big, Data_seedlings, Stages_select, Seedlings_s, tablepath) {
  
  plot_select_DE <- Stages_select$plotid %>% unique()
  plot_select_SK <- Seedlings_s$plotid %>% unique()
  rm(Stages_select, Seedlings_s)
  
  ## Taxon frequencies were calculated by first averaging per plot over multiple surveys, and then calculating a total average over all plots.
  
  Summary_DE <- Data_big %>%
    filter(plotid %in% plot_select_DE) %>%
    
    mutate(tax = str_replace(tax, "\\.", " ")) %>%
    mutate(count_ha = count, ## count is already per hectare, countarea == 1
           ba = pi * (dbh/2)^2 * 1e-6, # mm^2 to m^2)
           ba_ha = ba * count_ha) %>% 

    group_by(plotid, tax, obsid) %>%
    summarize(ba_ha =  sum(ba_ha, na.rm = T)) %>%
    
    group_by(plotid, tax) %>%
    summarize(ba_ha_plot =  mean(ba_ha, na.rm = T)) %>%
    
    group_by(plotid) %>%
    mutate(ba_ha_total_plot = sum(ba_ha_plot, na.rm = T), frac_ba_ha_plot = ba_ha_plot/ba_ha_total_plot) %>%
    
    ungroup() %>%
    mutate(ba_ha_total_avg = mean(ba_ha_total_plot, na.rm = T)) %>%
    
    group_by(tax) %>%
    summarize(ba_ha_avg = mean(ba_ha_plot, na.rm = T),
              frac_ba_ha_avg = mean(frac_ba_ha_plot, na.rm = T),
              ba_ha_total_avg = first(ba_ha_total_avg)) %>%
    filter(ba_ha_avg != 0) %>%
    
    ## tidy up for publishing
    arrange(desc(frac_ba_ha_avg)) %>%
    dplyr::select(-ba_ha_total_avg) %>%
    mutate(ba_ha_avg = round(ba_ha_avg, digits = 3)) %>%
    mutate(frac_ba_ha_avg = scales::percent(frac_ba_ha_avg, accuracy = 0.001)) %>%
    # mutate(across(2:3, ~ round(.x, digits = 2))) %>%
    setNames(c("Taxon", "Mean basal area m^2 ha^-1", "Mean percentage of the total basal area"))
  
  Summary_SK <- Data_seedlings %>%
    filter(sizeclass == "big") %>%
    filter(plotid %in% plot_select_SK) %>%

    group_by(plotid, taxon, year) %>%
    summarize(ba_ha =  sum(ba_ha, na.rm = T)) %>%
    
    group_by(plotid, taxon) %>%
    summarize(ba_ha_plot =  mean(ba_ha, na.rm = T)) %>%
    
    group_by(plotid) %>%
    mutate(ba_ha_total_plot = sum(ba_ha_plot, na.rm = T), frac_ba_ha_plot = ba_ha_plot/ba_ha_total_plot) %>%
    
    ungroup() %>%
    mutate(ba_ha_total_avg = mean(ba_ha_total_plot, na.rm = T)) %>%
    
    group_by(taxon) %>%
    summarize(ba_ha_avg = mean(ba_ha_plot, na.rm = T),
              frac_ba_ha_avg = mean(frac_ba_ha_plot, na.rm = T),
              ba_ha_total_avg = first(ba_ha_total_avg)) %>%
    filter(ba_ha_avg != 0) %>%
    
    ## tidy up for publishing
    arrange(desc(frac_ba_ha_avg)) %>%
    dplyr::select(-ba_ha_total_avg) %>%
    mutate(ba_ha_avg = round(ba_ha_avg, digits = 3)) %>%
    mutate(frac_ba_ha_avg = scales::percent(frac_ba_ha_avg, accuracy = 0.001)) %>%
    # mutate(across(2:3, ~ round(.x, digits = 2))) %>%
    setNames(c("Taxon", "Mean basal area m^2 ha^-1", "Mean percentage of the total basal area"))
  
  write.csv(Summary_DE, file.path(tablepath, "Taxa_freq_DE.csv"))
  write.csv(Summary_SK, file.path(tablepath, "Taxa_freq_SK.csv"))
  
  return(list(DE = Summary_DE, SK = Summary_SK))
}


## summarizeNFIs --------------------------------
# Stages_select  <- tar_read("Stages_select")
# Seedlings_s <- tar_read("Seedlings_s")
# Data_big <- tar_read("Data_big")
# Data_seedlings <- tar_read("Data_seedlings")


summarizeNFIs <- function(Data_big, Data_seedlings, Stages_select, Seedlings_s, tablepath) {
  
  DE_before <- Data_big
  SK_before <- Data_seedlings
  DE <- Stages_select
  SK <- Seedlings_s
  
  n_surveys_DE <- n_distinct(DE$obsid)
  n_surveys_SK <- 1
  
  # 
  # char_years_DE <- str_extract_all(unique(DE$obsid), "[0-9]{4}") %>% unlist() %>%
  char_years_DE <- DE_before$time %>% year() %>% table() %>% names() %>%
    paste(collapse = ".")
  
  char_years_SK <- str_extract_all(unique(SK$year), "[0-9]{4}") %>%
    unlist() %>%
    setdiff("2017") %>% ## 2017 is only in there because there was a resurvey of the relevé
    sort() %>%
    paste(collapse = "—")
  
  n_plots_before_DE <- n_distinct(DE_before$plotid)
  n_plots_before_SK <- n_distinct(SK_before$plotid)
  n_plots_DE <- n_distinct(DE$plotid)
  n_plots_SK <- n_distinct(SK$plotid)
  
  n_clusters_before_DE <- n_distinct(str_extract(DE_before$plotid, "(?<=_)([0-9]*?)(?=\\.)"))
  n_clusters_before_SK <- "."
  n_clusters_DE <- DE$clusterid %>% n_distinct()
  n_clusters_SK <- "."
  
  cas <- function(...) { c(sapply(list(...), as.character)) }
  
  column_title <- c("No. of surveys", "... years", "No. of plots", "... after selection", "No. of clusters", "... after selection")
  column_DE <- cas(n_surveys_DE, char_years_DE, n_plots_before_DE, n_plots_DE, n_clusters_before_DE, n_clusters_DE)
  column_SK <- cas(n_surveys_SK, char_years_SK, n_plots_before_SK, n_plots_SK, n_clusters_before_SK, n_clusters_SK)
  
  Table <- data.frame("." = column_title, "German.NFI" = column_DE, "Slovakian.NFI" = column_SK)
  print(Table)
  
  write.csv(Table, file.path(tablepath, "Summary_NFIs.csv"))

  return(Table)
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
    dplyr::select(any_of(c(id_select, predictor_select, disturbance_select))) %>%
    dplyr::mutate_at(predictor_select, function(x) replace(x, is.nan(x), NA)) ## na_if doesn't work with NaNs
  
  
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
    dplyr::select(any_of(c(id_select, disturbance_select, predictor_select))) %>%
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
# Stages_select <- tar_read("Stages_select")
# predictor_select <- tar_read("predictor_select")

scaleData <- function(Stages,
                      predictor_select
                      ) {
  
  #### Round counts and ba
  if(any(Stages$count_ha < 0.5 & Stages$count_ha != 0, na.rm = T)) warning("scaleData(): There were values < 0.5 in count_ha (in addition to zeroes). They have still been rounded in variable count_ha_r.")
  Stages$count_ha_r <- as.integer(round(Stages$count_ha))
  
  if(any(Stages$ba_ha < 0.5 & Stages$ba_ha != 0, na.rm = T)) warning("scaleData(): There were values < 0.5 in ba_ha (in addition to zeroes). They have still been rounded in variable ba_ha_r.")
  Stages$ba_ha_r <- as.integer(round(Stages$ba_ha))
  
  #### Scale predictors
  scalecolumns <- c(predictor_select) # , names(Stages)[grepl("^s_", names(Stages))]
  M_stages <- scale(st_drop_geometry(Stages[scalecolumns]))
  colnames(M_stages) <- paste0(scalecolumns, "_s")
  Stages <- cbind(Stages, M_stages)
  attr(Stages, "scaled:center") <- attr(M_stages, "scaled:center")
  attr(Stages, "scaled:scale") <- attr(M_stages, "scaled:scale")
  
  return(Stages)
}

## setLocLevel --------------------------------
# Stages_scaled <- tar_read("Stages_scaled")
# Env_cluster <- tar_read("Env_cluster")
# loclevel <- tar_read("loc")
setLocLevel <- function(Stages_scaled, Env_cluster, loc = c("plot", "nested", "cluster"),
                        id_select = c("clusterid", "clusterobsid", "methodid", "obsid", "plotid", "plotobsid", "tax", "taxid", "time")
                        ) {
  
  loclevel <- match.arg(loc)
  # locid <- if(loclevel == "plot") "plotid" else if(loclevel %in% c("cluster", "nested")) "clusterid"
  
  if(loclevel %in% c("cluster", "nested")) {
    
    Stages_loc <- st_drop_geometry(Stages_scaled)
    
    dropcols <- setdiff(intersect(names(Stages_loc), names(Env_cluster)), "clusterid") ## drop all common cols but clusterid
    E <- dplyr::select(Env_cluster, -all_of(dropcols))
    
    if(loclevel == "cluster") {
      
      id_select_cluster <- setdiff(id_select, c("plotid", "plotobsid"))
                                   
      Stages_loc %<>%
        mutate(clusterobsid = interaction(clusterid, obsid)) %>% ## In case it is dropped somewhere before
        group_by_at(id_select_cluster) %>%
        
        ### Summary per cluster. Env gets dropped here, will be joined again by plotid
        summarize(n_plots = n_distinct(plotid),
                  
                  ## Splines by first
                  across(paste("s", taxon_s, sep = "_"), first),
                  
                  ## Abundances and areas are summed up
                  across(c(count_obs, count_ha, count_ha_r,
                           ba_obs, ba_ha, ba_ha_r,
                           offset, offset_avg, offset_q1, offset_q3),
                         function(x) sum(x, na.rm = T)),
                  
                  .groups = "drop") %>%
        ungroup()

    }
    
    Stages_loc %<>%
      left_join(E, by = "clusterid") %>% # geometries are joined here, but crs is dropped
      st_sf(crs = st_crs(E)) %>%
      mutate(loc = clusterid)
  
  } else { # loclevel == "plot"
    
    Stages_loc <- Stages_scaled %>%
      mutate(loc = plotid)
  
  }
  
  attr(Stages_loc, "scaled:center") <- attr(Stages_scaled, "scaled:center")
  attr(Stages_loc, "scaled:scale") <- attr(Stages_scaled, "scaled:scale")
  
  return(Stages_loc)
}
