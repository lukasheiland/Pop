# Themes ------------------------------------------------------------------

theme_fagus <- function(...) {
  hrbrthemes::theme_ipsum(...)
  }


theme_empty <- function(){
  theme(panel.grid = element_line(colour = "transparent"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  }
