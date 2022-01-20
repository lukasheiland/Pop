## mapClusters --------------------------------
# Stages_select <- tar_read(Stages_select)

mapClusters <- function(Stages_select, pointcolor = "white", path = dir_publish) {
  
  ## ETRS89 Lambert Conformal Conic CRS
  projection <- crs("EPSG:3034")
  
  Clusters <- Stages_select %>%
    filter(!duplicated(clusterid)) %>%
    st_transform(crs = projection)

  adminid <- c("DE")
  Area <- eurostat::get_eurostat_geospatial(resolution = "01", # hightest possible 1:1million
                                             nuts_level = 0) %>%  # Administrative level. There are quite a few small scale levels, all with different ids!
    dplyr::filter(id %in% !!adminid) %>%
    st_transform(crs = projection)
  
  Elevation <- elevatr::get_elev_raster(Area, z = 8) %>% # Zoom level: determines resolution, crs is determined by Area. 8 is minimum for not getting artefacts after projection
    raster::crop(as_Spatial(Area)) %>%
    raster::mask(as_Spatial(Area))

  ext <- raster::extent(Elevation)
  
  E <- Elevation %>%
    raster_to_matrix()
  
  ### Rayshader
  ## See tutorial at https://www.tylermw.com/adding-open-street-map-data-to-rayshader-maps-in-r/
  
  # spheretexture <- create_texture(lightcolor = "#FFFBEF", shadowcolor = "#B0B0B0", leftcolor = "#ECECE0", rightcolor = "#ECECEA", centercolor = "white")
  spheretexture <- "imhof1"
  elevtexture <- terrain.colors(256)
  
  
  map <- E %>%
    height_shade(texture = elevtexture, range = c(-500, 1400)) %>% # 
    add_overlay(sphere_shade(E, sunangle = 315, texture = spheretexture,
                             zscale = 4, colorintensity = 5), alphalayer = 0.5) %>%
    add_shadow(texture_shade(E, detail = 8/10, contrast = 4, brightness = 100), 0.1) %>%
    
    add_overlay(generate_point_overlay(Clusters, extent = ext, heightmap = E,
                                       pch = 1, color = pointcolor, size = 30), rescale_original = T) %>%
    
    # add_overlay(generate_line_overlay(Area, extent = ext, heightmap = E)) %>%
    
    ## add this after point overlay, otherwise it will make the background black
    add_shadow(lamb_shade(E, zscale = 10), 0)
  
    

  png(paste0(path, "/", "Map_clusters", ".png"), width = ncol(map), height = nrow(map))
  plot_map(map) ## Handed on to raster::plotRGB
  dev.off()
  
  return(map)
}

