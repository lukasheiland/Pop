# ————————————————————————————————————————————————————————————————————————————————— #
# Functions for comparing environmental ranges   ---------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #


## selectRange --------------------------------
# Stages  <- tar_read("Stages_select_env")
# file_env <- tar_read("Stages_select")
selectRange <- function(Stages, taxon = "Fagus.sylvatica") {
  
  ## Filtering by taxon
  # Stages %<>%
  #   filter(as.character(tax) == taxon) %>%
  #   group_by(plotid, tax, clusterid, phCaCl_esdacc, waterLevel_loc) %>%
  #   summarize(keep = any(anyFagus), .groups = "drop") %>%
  #   filter(keep)
  
  Stages %<>%
    distinct(plotid, .keep_all = TRUE)

  return(Stages)
}


## extractRange --------------------------------
# Raster_EAFTS  <- tar_read("Raster_EAFTS_range")
# rasters_env <- tar_read("rasters_env_range")
extractRange <- function(Raster_EAFTS, rasters_env) {
  
  # plot(rasters_env[[2]])
  
  R <- Raster_EAFTS > 0.01 # we extracted environmental variables used the occurrence prob of >1%
  ## plot(R)
  D <- as.data.frame(R, xy = TRUE)
  
  # crs_common <- crs(R)
  n_rasters <- length(rasters_env)
  rasters_reprojected <- mclapply(rasters_env, raster::projectRaster, to = R, mc.cores = n_rasters)
  rm(rasters_env)
  rasters_df <- lapply(rasters_reprojected, as.data.frame, xy = FALSE) ## xy are now the same as in R
  rasters_df <- mapply(function(r, n) setNames(r, n), rasters_df, names(rasters_df), SIMPLIFY = F)
  
  Range <- bind_cols(D, rasters_df) %>%
    rename(greater1percent = layer) %>%
    dplyr::filter(greater1percent) %>%  ## the condition that occurrence prob >1%, also dropping NAs
    select_at(names(rasters_df))

  return(Range)
}


## summarizeRange --------------------------------
# Ranges  <- tar_read("Ranges_range")
# predictor <-  names(tar_read(predictor_range))
summarizeRange <- function(Ranges, predictor, path = tar_read(dir_publish)) {
  
  S <- Ranges %>%
    dplyr::select(all_of(c("origin", predictor))) %>% 
    group_by(origin) %>%
    summarize_at(predictor, function(x) { paste0(formatNumber(mean(x, na.rm = T), signif.digits = 5), " ± ",
                                                 formatNumber(sd(x, na.rm = T), signif.digits = 5)) })
  
  write.csv(S, file.path(path, "Summary_range_Fagus.sylvatica.csv"))
  
  return(S)
}


## plotRange --------------------------------
# Ranges  <- tar_read("Ranges_range")
# predictor <-  names(tar_read(predictor_range))
plotRange <- function(Ranges, predictor,
                      path = tar_read(dir_publish), color = tar_read(twocolors), themefun = tar_read(themefunction)) {
  
  ranges <- Ranges %>%
    dplyr::select(all_of(c("origin", predictor))) %>%
    pivot_longer(cols = predictor, names_to = "variable", values_to = "value") %>%
    mutate(origin = if_else(origin == "DE", "German NFI", origin)) %>%
    mutate(variable = case_when(variable == "phCaCl_esdacc" ~ "soil pH",
                                variable == "cwbYear_aclim" ~ "climatic water balance [mm/a]",
                                TRUE ~ variable)) %>% 
    split(., .$variable)
    
  
  plotR <- function(R) {
    plot <- ggplot(R, mapping = aes(x = value, y = origin, group = origin)) +
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
      facet_wrap(~variable) + ## just to have a nice title, although there is only one variable
      themefun() +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  
  plots <- lapply(ranges, plotR)
  plotgrid <- cowplot::plot_grid(plotlist = plots, ncol = 1, rel_widths = 2, labels = "AUTO")
  ggsave(file.path(path, "Plot_range_Fagus.sylvatica.pdf"), plotgrid, device = "pdf", height = 6, width = 6)
  
  return(plots)
}

