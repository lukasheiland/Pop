# Targets setup -----------------------------------------------------------
### Library
library(targets)
library(tarchetypes)

### Source the functions
source("Model_2021-03_ba/Fit_ba_functions.R")

### Options
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("biglm", "dplyr", "ggplot2", "readr", "tidyr"))


# Pipeline ----------------------------------------------------------------

list(
  tar_target(
    inv_data_file,
    "Data Inv/",
    format = "file"
  ),
  
  tar_target(
    inv_data,
    read_csv(inv_data_file, col_types = cols())
  ),
  
  tar_target(
    data,
    inv_data %>%
      filter(!is.na(abundance))
  ),
  
  tar_target(hist, create_plot(data)),
  
  tar_target(fit, biglm(Ozone ~ Wind + Temp, data)),
  
  tar_render(report, "index.Rmd")
)



