
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
          
          panel.grid = element_line(colour = "#222222"), ## only slightly less weight than full black
          axis.ticks = element_line(colour = "black", size = rel(0.8), lineend = "round"),
          axis.title.x = element_text(margin = margin(t = 3)), #add margin to x-axis title
          axis.title.y = element_text(margin = margin(r = 2.5)),
          axis.text.x = element_text(margin = margin(t = 6)),
          axis.text.y = element_text(margin = margin(r = 4.5))
          ) + 
    
    theme(text = element_text(size = 14, family = "Helvetica")) # theme(text = element_text(size = 14, family = "SF Compact Rounded")) # theme(text = element_text(size = 14, family = "TGL17"))
  }

theme_empty <- function(){
  theme(panel.grid = element_line(colour = "transparent"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  }
