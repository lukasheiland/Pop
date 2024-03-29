
# Prerequisites -----------------------------------------------------------
# remotes::install_version("Rttf2pt1", version = "1.3.8")
# library(showtext)
# font_add(family = "TGL17", regular = "Theme/tgl17.ttf")
# font_add(family = "SF Compact Rounded", regular = "Theme/SF-Compact-Rounded-Regular.otf") # bold = "Theme/SF-Compact-Rounded-Medium.otf", r
# font_add(family = "VAG Rounded", regular = "Theme/VAG Rounded Regular.ttf")

# Themes ------------------------------------------------------------------

theme_fagus <- function(...) {
  theme_linedraw(...) +   # hrbrthemes::theme_ipsum(...) +
    theme(strip.background = element_rect_round(radius = unit(4, "pt")),
          panel.background = element_rect_round(radius = unit(4, "pt")),
          panel.border = element_rect_round(radius = unit(4, "pt"), size = rel(1.5)),
          panel.spacing = unit(2, "lines"),
          
          panel.grid = element_line(colour = "#444444"), ## only slightly less weight than full black
          # panel.grid.major.y.right = element_blank(), ## secondary y-axis is not addressable in ggplot as of now, specified via scale_*(sec.axis=*), on the right
          # panel.grid.minor.y.right = element_blank(),
          axis.ticks = element_line(colour = "black", size = rel(0.8), lineend = "round"),
          axis.ticks.y.right = element_line(colour = "grey40"),
          
          axis.title.x = element_text(margin = margin(t = 3)), #add margin to x-axis title
          axis.title.y = element_text(margin = margin(r = 2.5)),
          axis.title.y.right = element_text(colour = "grey40"),
          
          axis.text.x = element_text(margin = margin(t = 6)),
          axis.text.y = element_text(margin = margin(r = 4.5)),
          axis.text.y.right = element_text(colour = "grey40")
          ) + 
    
    theme(text = element_text(size = 14, family = "Helvetica")) # theme(text = element_text(size = 14, family = "SF Compact Rounded")) # theme(text = element_text(size = 14, family = "TGL17"))
  }

theme_empty <- function() {
  theme_fagus() +
  theme(panel.grid = element_line(colour = "transparent"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
  }
