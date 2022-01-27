
# Prerequisites -----------------------------------------------------------
# remotes::install_version("Rttf2pt1", version = "1.3.8")
library(extrafont)
extrafont::font_import(paths = "Theme/", prompt = F)
# fonts()

# Themes ------------------------------------------------------------------

theme_fagus <- function(...) {
  # hrbrthemes::theme_ipsum(...) +
  theme_linedraw(...) +
    theme(text = element_text(size = 15, family = "TGL 0-17"))
  }


theme_empty <- function(){
  theme(panel.grid = element_line(colour = "transparent"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  }
