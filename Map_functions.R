
## theme_map --------------------------------
# th <- tar_read(themefunction)
theme_map <- function(th = theme_fagus){
  th() +
    theme(# panel.grid = element_line(colour = "transparent"), ## allow for a grid
      panel.border = element_blank(),
      panel.background = element_blank(),
      # axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank()
    )
}


## mapClusters --------------------------------
# Stages_select <- tar_read(Stages_select)
# themefun <- tar_read(themefunction)
# pack <- c("dplyr", "ggplot2", "magrittr", "sf", "raster", "eurostat", "elevatr", "terrainr", "rayshader", "ggspatial", "elementalist")
# lapply(pack, require, character.only = TRUE)

mapClusters <- function(Stages_select, color = NULL, themefun, path = dir_publish) {
  
  if (is.null(color)) color <- c("#00A4DB", "#3EFC38", "#FFEE00", "#FFFFFF") # lighter versions of the twocolors, for shading
  
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
  # spheretexture <- "imhof3"
  spheretexture <- create_texture("#f5dfca","#63372c","#dfa283","#195f67","#c2d1cf", cornercolors = c("#ffc500", "#387642", "#d27441","#296176")) ## manual imhof
  # elevtexture <- terrain.colors(256)
  elevtexture <- grDevices::colorRampPalette(color)
  
  relief <- E %>%
    height_shade(texture = elevtexture(256), range = c(-700, 3500)) %>% # some elevations, like coalmines go rather low. The range is set up in a way that the green is mostly in the middle
    add_overlay(sphere_shade(E, texture = spheretexture, sunangle = 0, zscale = 1, colorintensity = 5),
                alphalayer = 0.75) %>%
    add_shadow(ray_shade(E, lambert = FALSE, sunangle = 315, sunaltitude = 55, zscale = 10, multicore = T), max_darken = 0.6) %>% # shadowmap E, will be rescaled to match hillshade
    add_shadow(lamb_shade(E, zscale = 10, sunangle = 315, sunaltitude = 55), max_darken = 0.5)
  
  ## add these before lamb_shade, otherwise it will make the background black
    # add_overlay(generate_point_overlay(Clusters, extent = ext, heightmap = E,
    #                                    pch = 1, color = "white", size = 30), rescale_original = T) %>%
    # add_overlay(generate_line_overlay(Area, extent = ext, heightmap = E)) %>%
  
  # plot_map(relief)

  ## Convert rayshader array to RasterStack with some tricks
  B <- brick(relief) %>% stack()
  extent(B) <- extent(Elevation)
  crs(B) <- crs(Elevation)
  names(B) <- c("red", "green", "blue", "alpha")
  # res(B) <- res(B) * 10
  
  # raster::plotRGB(B, r = 1, g = 2, b = 3, scale = 1.0)
  
  C <- st_coordinates(Clusters) %>% as.data.frame()
  
  map <- ggplot() +  ## annotations need aes, does not actually do anything
    ## Base raster layer
    # coord_sf(crs = 3034) + ## implicitly added through geom_sf
    geom_spatial_rgb(data = B, mapping = aes(x = x, y = y, r = red, g = green, b = blue), scale = 1.0) +

    ## For making a legend. (Include dataframe with all elevations in ggplot data.)
    # scale_color_continuous(elevtexture) +
    # guides(color = guide_colourbar(barwidth = 0.45, barheight = 10, ticks.linewidth = 0.8, draw.ulim = TRUE, draw.llim = TRUE, frame.linewidth = 1, nbin = 100, title = "Elevation")) +
    # theme(legend.position = c(0.1, 0.8)) +
    
    ## Points
    geom_sf(data = Clusters, col = "black", fill = "white", pch = 21, size = 1.4, stroke = 1.25) +
    
    annotation_scale(width_hint = 0.3,
                     style = "bar",
                     location = "br",
                     bar_cols = c("black", "white"),
                     line_col = "black",
                     text_col = "black",
                     text_cex = 1.1,
                     pad_x = unit(0.3, "cm")) +
    
    annotation_north_arrow(which_north = "true",
                           location = "tl",
                           pad_x = unit(0.6, "cm"),
                           style = north_arrow_minimal(line_col = "black", text_col = "black", fill = c("black"), text_size = 17)) +
    theme_map(themefun)
    
  
  ggsave(paste0(path, "/", "Map_clusters", ".pdf"), map, device = "pdf", height = 9, width = 8)
  ggsave(paste0(path, "/", "Map_clusters", ".png"), map, device = "png", height = 9, width = 8)
  
  return(map)
}

