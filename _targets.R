# ————————————————————————————————————————————————————————————————————————————————— #
# Targets setup ------------------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

### Installation
# remotes::install_github("ropensci/targets")
# remotes::install_github("ropensci/tarchetypes")
# remotes::install_github("ropensci/stantargets")
# cmdstanr::install_cmdstan()

### Library
library(targets)
library(tarchetypes)
# library(stantargets)

### Future
library(future)
library(future.callr)
plan(callr) ## "It is crucial that future::plan() is called in the target script file itself"

### Source the functions
source("Wrangling_functions.R")
source("Map_functions.R")
source("Range_functions.R")
source("Fit_seedlings_functions.R")
source("Fit_functions.R")
source("Posterior_functions.R")
source("Theme/Theme.R")
source("Helper_functions.R")


### Options
options(tidyverse.quiet = TRUE)
tar_option_set(error = "abridge")

package <- c("dplyr", "ggplot2", "tidyr", "magrittr", "glue", "forcats", "vctrs", "tibble", "stringr", "knitr", # "multidplyr" ## extended tidyverse
             "lubridate", "DescTools",
             "sf", "raster", "rasterVis", ## for correct loading of environmental data
             "MASS", "mgcv", "glmnet", "itsadug", "interp",
             "cmdstanr", "rstan", "brms", "posterior", "bayesplot", "tidybayes", "parallel", "DHARMa", "priorsense", # "chkptstanr",
             "stargazer", "cowplot", "gridExtra", "hrbrthemes", "showtext", "ggallin", "ggridges", "elementalist",  "ggspatial", "GGally", "ggforce", # "emojifont",
             "scales", "gganimate", "metR", "colorspace", "gt", "gtExtras", "svglite",
             "future.apply")
tar_option_set(packages = package)

### Operation system
onserver <- Sys.info()["sysname"] != "Darwin"
if (!onserver) tar_config_set(store = "_targets.nosync")


# ————————————————————————————————————————————————————————————————————————————————— #
# Inner pipelines -----------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## Settings pipeline ------------------------------------------------------
targets_settings <- list(
  
  ## Whether to fit the model with a population ("loc", i.e. location) structure ...
  ##  "plot" — where the populations correspond to a plot in the data (n_locs == n_plots); environment joined by plot; counts/offsets correspond to a plot
  ##  "nested" — where the populations correspond to a cluster, but get fitted to data on the plot level (n_locs == n_clusters); environment joined by cluster; counts/offsets correspond to a plot
  ##  "cluster" — where the populations correspond to a cluster and get fitted to the sum within a cluster (n_locs == n_clusters); environment joined by cluster; counts/offsets correspond to the sum of the plots within a cluster
  tar_target(loc, c("plot", "nested", "cluster")[1]),
  
  ## No. of locations to subset (currently only for loc == "plot")
  tar_target(n_locations, 1000),
  tar_target(n_locations_env, 1200),
  
  ## Threshold to discriminate A and B [mm]
  ## 160 is the 10%tile, 185 is the 15%tile, 207 is the 20%tile, 228 is the 25%tile of pure measured trees, i.e. without area standardization
  ## lower in the data is 100, so that: 100mm > A > threshold_dbh > B
  tar_target(threshold_dbh, 180), ## [mm]
  
  ## Upper sampling radius
  ## 	- All trees above a sampling radius of 15m were dropped, which is about the 99%tile of measured trees (15.71m).
  ## - The radius of 15m corresponds to the threshold radius of trees with dbh = 60cm
  ##    - dbh_threshold = radius_threshold/c with c == 25
  ##    - Alternatives: 99% radius == 1571 cm, 95% radius == 1188,  96% radius == 1242, 97% 1310.91
  tar_target(radius_max, 15000), ## [mm]
  
  ## Vector of taxa to select. All others will be lumped into "other".
  tar_target(taxon_select, c("Fagus.sylvatica")),
  
  ## Regeneration classes to include
  ## here we select all 20cm height <= trees < 7cm dbh: regglass_select <- c("h[20,50)" = 1, "h[50,130)" = 2, "hd[130,Inf)[0,5)" = 3, "d[5,6)" = 4, "d[6,7)" = 5)
  ## These are all size classes that are consistent across the three surveys.
  tar_target(regclass_select, c("h[20,50)" = 1, "h[50,130)" = 2, "hd[130,Inf)[0,5)" = 3, "d[5,6)" = 4, "d[6,7)" = 5, "hd[130,Inf)[0,7)" = 9)), ## [mm]; the last level (9) is a special category of counts in the second survey DE_BWI_2 over all classes
  
  ## Weakly informative priors.
  ## If provided, they are prioritzed over the fitted priors
  tar_target(weakpriors,
             ## Priors are organized like the parameter data structure but with an additional dimension in the case of a vector row of sds.
             list(
               prior_b_log = c(-3.2, 1),
               prior_c_a_log = c(-7, 2),
               prior_c_b_log = c(-7, 2),
               prior_c_j_log = c(-10, 2),
               prior_g_log = c(-5, 2),
               prior_h_log = c(-4, 2),
               prior_l_log = c(4, 2),
               # prior_r_log = cbind(Fagus = c(4, 1), others = c(4, 1)),
               prior_s_log = c(-6, 2)
             )
  ),
  
  tar_target(weakpriors_env,
             list(prior_b_log = c(-3, 1),
                  prior_c_a_log = c(-7, 2),
                  prior_c_b_log = c(-7, 2),
                  prior_c_j_log = c(-11, 2),
                  prior_g_log = c(-4, 1),
                  prior_h_log = c(-3, 1),
                  prior_l_log = c(5, 2),
                  prior_r_log = c(4, 1),
                  prior_s_log = c(-6, 2)
                  )
  ),
  
  tar_target(twocolors, c("#007A7F", "#FFCC11")), ## Formerly "#119973", Other greens include: Spanish Green '#208E50', Feldgrau '#3A7867' # Jade: #0FA564 # See green: #0E8C55
  tar_target(plotsettings,
             list(jitter_points = position_jitter(seed = 1, width = 0.1, height = 0.1),
                  divergingcolorscale = function(...) scale_color_continuous_divergingx(palette = "BrBG", ...), # function(n) rev(pals::ocean.curl(n))
                  divergingfillscale = function(...) scale_fill_continuous_divergingx(palette = "BrBG", ...),
                  lims_colorscale = 0:1,
                  axislabs = labs(x = "soil pH", y = "soil water level"),
                  aspect = theme(aspect.ratio = 1.1),
                  removeylabs = theme(axis.title.y = element_blank(), axis.text.y = element_blank()),
                  removeleftylabs = theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank()),
                  hlines = lapply(-1:1, function(y) geom_hline(aes(yintercept = y), linetype = 2, col = "grey40", size = 0.4)),
                  height_plot = 5.3)),
  tar_target(themefunction, theme_fagus)
  
)



## Paths pipeline ------------------------------------------------------
targets_paths <- list(
  
  ## Directories
  tar_target(dir_fit, file.path("Fits.nosync/")),
  tar_target(dir_publish, file.path("Publish.nosync/")),
  
  ## Data files
  tar_target(file_DE_big,
             "Inventory.nosync/DE BWI/Data/DE_BWI_big_abund.rds",
             format = "file"),
  tar_target(file_DE_big_status,
             "Inventory.nosync/DE BWI/Data/DE_BWI_big_2-3_status.rds",
             format = "file"),
  tar_target(file_DE_small,
             "Inventory.nosync/DE BWI/Data/DE_BWI_small_abund.rds",
             format = "file"),
  tar_target(file_DE_env,
             'Inventory.nosync/DE BWI/Data/DE_BWI_Env_sf.rds',
             format = "file"),
  tar_target(file_DE_geo,
             'Inventory.nosync/DE BWI/Data/DE_BWI_geo.rds',
             format = "file"),
  tar_target(file_SK,
             "Inventory.nosync/SK NIML/Data/SK_NIML_complete.rds",
             format = "file"),
  tar_target(file_SK_fullgrid,
             "Inventory.nosync/SK NIML/Data/SK_NIML_complete_fullgrid.rds",
             format = "file"),
  # tar_target(file_Taxa,
  #            'Inventory.nosync/Taxa/Taxa.csv',
  #            format = "file")
  
  ## References
  # dput(.packages()) ## does not work directly because of incomplete metadata in some
  tar_target(bibliography,
             knitr::write_bib(c("scales", "ggspatial", "DHARMa",
                                "tidybayes", "bayesplot", "posterior", "brms", "rstan", "cmdstanr", 
                                "MASS", "mgcv", "nlme", "rasterVis", "lattice", "raster", "sp", 
                                "sf", "DescTools", "lubridate", "stringr", "tibble", "vctrs"), file.path(dir_publish, "Packages.bib"))
             )
  
)



## Parameter spec pipeline ------------------------------------------------------
targets_parname <- list(

  #### Generated quantities (including those for env)
  tar_target(simname_prior,
             c("y_hat_prior", "y_hat_prior_rep", "y_hat_prior_offset", "y_hat_prior_rep_offset", "y_prior_sim")),
  
  tar_target(statename_environmentalko_fracdiff_env,
             paste0("ba_frac_diff_fix_ko_", 1:2,"_env_", rep(c(names(parname_plotorder), "g_c_j_s", "h_c_a", "b_c_b"), each = 2))), # "b_other_s"
  
  tar_target(statename_environmentalko_env,
             c(paste0("ba_fix_ko_", 1:2,"_env_", rep(c(names(parname_plotorder), "g_c_j_s", "h_c_a", "b_c_b"), each = 2)),
               paste0("major_fix_ko_", 1:2,"_env_", rep(c(names(parname_plotorder), "g_c_j_s", "h_c_a", "b_c_b"), each = 2)),
               paste0("ba_frac_fix_ko_", 1:2,"_env_", rep(c(names(parname_plotorder), "g_c_j_s", "h_c_a", "b_c_b"), each = 2)),
               statename_environmentalko_fracdiff_env)),
  
  tar_target(statename,
             c("ba_init", "ba_fix", "J_init", "J_fix", "A_init", "A_fix", "B_init", "B_fix", "ba_sum_fix_loc",
               "ba_tot_init", "ba_tot_fix", "ba_frac_init", "ba_frac_fix",
               "ba_fix_ko_b", "ba_fix_ko_s", "ba_fix_ko_2_b", "ba_fix_ko_2_s", "ba_fix_ko_b_l_r",
               "major_init", "major_fix", "dominant_init", "dominant_fix",
               "major_fix_ko_b", "major_fix_ko_s", "major_fix_ko_2_b", "major_fix_ko_2_s",
               paste0("major_fix_switch_", c(names(parname_plotorder), "l_r", "h_c_a", "g_c_j_l_r_s", "g_c_j_s", "g_l_r_s", "g_s", "b_c_b", "b_c_a_c_b_h")),
               "ba_fix_other_s", "ba_frac_fix_other_s", "major_fix_other_s",
               "ba_fix_ko_b", "ba_fix_ko_s", "ba_fix_ko_2_b", "ba_fix_ko_2_s",
               statename_environmentalko_fracdiff_env,
               paste0("ba_fix_switch_", c(names(parname_plotorder), "l_r", "h_c_a", "g_c_j_l_r_s", "g_c_j_s", "g_l_r_s", "g_s", "b_c_b", "b_c_a_c_b_h")))),
  tar_target(basalareaname, statename[str_starts(statename, "ba_")]),
  tar_target(majorname, statename[str_starts(statename, "major_")]),
  tar_target(majorname_main, c("major_init", "major_fix", paste0("major_fix_switch_", c(names(parname_plotorder), "l_r")))),
  tar_target(majorname_supp, c("major_init", "major_fix", setdiff(majorname, majorname_main))),
  tar_target(statename_select, setdiff(statename, c(statename_environmentalko_fracdiff_env))),
  
  # tar_target(contribname_init, paste0("sum_ko_", 1:2, "_", rep(c(names(parname_plotorder), "b_c_b"), each = 2), "_fix")),
  # tar_target(contribname_prop, paste0("sum_ko_", 1:2, "_prop_", rep(c(names(parname_plotorder), "b_c_b"), each = 2), "_fix")),
  # tar_target(contribname_switch, paste0("sum_switch_", c(names(parname_plotorder), "b_c_b"), "_fix")),
  # tar_target(contribname_init_env, paste0("sum_ko_", 1:2, "_", rep(c(names(parname_plotorder), "b_c_b"), each = 2), "_fix")),
  # tar_target(contribname_avg_env, paste(contribname_init_env, "avg", sep = "_")),
  # tar_target(contribname_counterfactual_env, paste0("sum_ko_", 1:2, "_", rep("sKoEnvB", each = 2), "_fix")), ## not used
  # tar_target(contribname_env, c(contribname_init_env, contribname_avg_env, contribname_counterfactual_env)),
  # tar_target(contribname, c(contribname_init, contribname_prop, contribname_switch, contribname_env)),

  tar_target(simname_posterior,
             c(statename, # contribname,
               "y_sim",
               paste0("greater_", names(parname_plotorder)),
               "converged_fix", "converged_fix_avg", "iterations_fix", "iterations_fix_avg", "eps_ba_fix", "eps_ba_fix_avg",
               "fixiter_max", "fixiter_min")),
  
  
  #### Parameters
  tar_target(parname,
             c("phi_obs",
               "b_log", "c_a_log", "c_b_log", "c_j_log", "g_log", "h_log", "l_log", "r_log", "s_log")),
  
  tar_target(parname_hyper,
             c("center_hyper", paste0(setdiff(parname_plotorder, "l_log"), "_spread_hyper"))),
  
  tar_target(parname_lim,
             c("b_lim_init_log", "g_lim_init_log", "h_lim_init_log",
               "b_lim_fix_log", "g_lim_fix_log", "h_lim_fix_log")),
  
  tar_target(parname_sim,
             setdiff(parname, c("theta", "phi_obs", "sigma_k_loc"))),

  tar_target(parname_plotorder,
             c(l = "l_log", r = "r_log", c_j = "c_j_log", s = "s_log", g = "g_log", c_a = "c_a_log", h = "h_log", b = "b_log", c_b = "c_b_log" )),
  
  tar_target(parname_lim_plotorder,
             c(l = "l_log", r = "r_log", c_j = "c_j_log", s = "s_log", g = "g_log", g_lim_init = "g_lim_init_log", g_lim_fix = "g_lim_fix_log",
               c_a = "c_a_log", h = "h_log", h_lim_init = "h_lim_init_log", h_lim_fix = "h_lim_fix_log",
               b = "b_log", b_lim_init = "b_lim_init_log", b_lim_init = "b_lim_fix_log", c_b = "c_b_log" )),
  
  tar_target(parname_lim_plotorder_plot,
             c(g = "g_log", g_lim_init = "g_lim_init_log", g_lim_fix = "g_lim_fix_log",
               h = "h_log", h_lim_init = "h_lim_init_log", h_lim_fix = "h_lim_fix_log",
               b = "b_log", b_lim_init = "b_lim_init_log", b_lim_init = "b_lim_fix_log")),
  
  tar_target(parname_loc,
             c("state_init", "L_loc")),
  
  tar_target(parname_sigma,
             c(l = "sigma_l", r = "sigma_r", c_j = "sigma_c_j", s = "sigma_s", g = "sigma_g", c_a = "sigma_c_a", h = "sigma_h", b = "sigma_b", c_b = "sigma_c_b" )),
  tar_target(parname_alpha,
             c(l = "alpha_l", r = "alpha_r", c_j = "alpha_c_j", s = "alpha_s", g = "alpha_g", c_a = "alpha_c_a", h = "alpha_h", b = "alpha_b", c_b = "alpha_c_b" )),
  
  
  #### Additional parameters and generated quantities, only used in the env fit
  tar_target(parname_center_env,
             sort(apply(expand.grid(paste0(setdiff(parname_plotorder, "l_log"), "_center_env"), 1:2), 1, paste0, collapse = ""))),
  tar_target(parname_spread_env,
             sort(apply(expand.grid(paste0(setdiff(parname_plotorder, "l_log"), "_spread_env"), 1:2), 1, paste0, collapse = ""))),
  tar_target(parname_vertex_env,
             sort(c(parname_center_env, parname_spread_env))),
  
  tar_target(parname_J_env,
             sort(c(parname_plotorder[1:4],
                    sapply(parname_plotorder[1:4], function(x) parname_vertex_env[str_detect(parname_vertex_env, x)]) %>% unlist()
             ))),
  tar_target(parname_A_env,
             sort(c(parname_plotorder[5:6],
                    sapply(parname_plotorder[5:6], function(x) parname_vertex_env[str_detect(parname_vertex_env, x)]) %>% unlist()
             ))),
  tar_target(parname_B_env,
             sort(c(parname_plotorder[7:9],
                    sapply(parname_plotorder[7:9], function(x) parname_vertex_env[str_detect(parname_vertex_env, x)]) %>% unlist()
             ))),
  
  tar_target(parname_Beta_env,
             paste0("Beta_", setdiff(parname_plotorder, "l_log"))),
  
  tar_target(parname_sim_environmental_env, ## JAB model parameters that are environmental
             c("B_log", "C_a_log", "C_b_log", "C_j_log", "G_log", "H_log", "L_loc", "L_log", "R_log", "S_log")),
  
  tar_target(parname_lim_environmental_env, ## JAB model parameters that are environmental
             c("B_lim_init_log", "G_lim_init_log", "H_lim_init_log", "B_lim_fix_log", "G_lim_fix_log", "H_lim_fix_log")),
  
  tar_target(parname_trajectories_environmental_env, ## environmental JAB model parameters + environmental initial values
             c(parname_sim_environmental_env, "state_init")),
  
  tar_target(parname_environmental_env,  ## anything that is environmentally-relevant, except for statename_environmentalko_.*_env
             c(parname_sim_environmental_env,
               parname_lim_environmental_env,
               # contribname_env,
               "ba_init", "ba_fix", "ba_frac_init", "ba_frac_fix", "ba_frac_fix_other_s",
               "major_fix", "major_init",
               "ba_fix_ko_b_l_r")),
  
  ## these are mainly for distinguishing different plots
  tar_target(parname_environmental_ba_env, parname_environmental_env[str_starts(parname_environmental_env, "ba")]),
  tar_target(parname_environmental_ba_select_env, setdiff(parname_environmental_ba_env, c("ba_init", "ba_frac_fix", "ba_frac_init", "ba_frac_fix_other_s"))), ## remove fractions that will be fitted with a binommial response
  tar_target(parname_environmental_binomial_env, c("major_init", "major_fix", "ba_frac_fix", "ba_frac_fix_other_s")),
  tar_target(parname_environmental_gaussian_env, setdiff(parname_sim_environmental_env, c("L_loc", "L_log"))), ## note: these have only tax %in% 1:2
  tar_target(parname_environmental_diff_env, setdiff(statename_environmentalko_fracdiff_env, "ba_frac_diff_fix_ko_2_env_b_other_s")), ## note: 1. these have only tax == 0
  tar_target(parname_environmental_lim_env, c(parname_lim_environmental_env, "G_log", "H_log", "B_log")), ## note: 1. these have only tax == 0
  tar_target(parname_environmental_par_env, c(parname_environmental_lim_env, "L_loc")),
  
  tar_target(parname_loc_env, ## anything that is local
             c(parname_environmental_env, statename_environmentalko_env,
               "state_init", str_to_sentence(names(parname_plotorder)))),
  
  
  
  #### Variables to exclude from any summary
  tar_target(par_exclude,
             c("y_hat", "L_loc_log", "K_loc_log_raw", "L_loc", "K_loc", "state_init", "state_init_raw", "state_init_log", "m")),
  
  tar_target(helper_exclude,
             c("Fix", "Fix_avg", "Fix_ko_s", "Fix_ko_1_b_l_r", "Fix_ko_2_b_l_r",
               paste0("Fix_ko_", 1:2,"_env_", rep(c(names(parname_plotorder), "b_c_b"), each = 2)),
               paste0("Fix_switch_", c(names(parname_plotorder), "b_c_b", "b_c_a_c_b_h", "g_c_j_l_r_s", "g_c_j_s", "g_l_r_s", "g_s", "h_c_a", "l_r")),
               "B_log_raw", "C_a_log_raw", "C_b_log_raw", "C_j_log_raw", "G_log_raw", "H_log_raw", "L_log_raw", "R_log_raw", "S_log_raw",
               str_to_sentence(names(parname_plotorder)),
               "vector_b_log_prior", "vector_c_a_log_prior", "vector_c_b_log_prior", "vector_c_j_log_prior", "vector_s_log_prior",
               "phi_obs_rep", "phi_obs_rep_prior", "theta_rep",
               "avg_state_init",
               paste0("avg_", names(parname_plotorder)), "avg_L_loc",
               "fixiter_max", "fixiter_min", "eps_ba_fix",
               "c_j_log_10",
               "log_prior", "log_lik", "lp__", "state_init_log_raw", "state_2", "state_2_raw", "state_3", "state_3_raw")),
  
  ## do not exclude: sort(apply(expand.grid(paste0(setdiff(parname_plotorder, "l_log"), "_spread_env"), 1:2, "_100"), 1, paste0, collapse = ""))
  
  tar_target(rep_exclude,
             c("phi_obs_rep", "phi_obs_rep_prior", "theta_rep",
               "y_hat_rep", "y_hat_prior_rep", "y_hat_offset", "y_hat_rep_offset", "y_hat_prior_rep_offset")),
  
  tar_target(exclude,
             c(par_exclude, helper_exclude, rep_exclude, parname_sigma, parname_alpha, simname_prior, simname_posterior, parname_loc_env))
)



## Wrangling pipeline ------------------------------------------------------
targets_wrangling <- list(

  ## Read data files
  # tar_load(starts_with("Data"))
  list(
    tar_target(Data_big, readRDS(file_DE_big)),
    tar_target(Data_big_status, readRDS(file_DE_big_status)),
    tar_target(Data_small, readRDS(file_DE_small)),
    tar_target(Data_env, readRDS(file_DE_env)),
    tar_target(Data_geo, readRDS(file_DE_geo)),
    
    tar_target(Data_seedlings, readRDS(file_SK)),
    tar_target(Data_seedlings_fullgrid, readRDS(file_SK_fullgrid))
    
    # tar_target(Taxa, read.csv(file = file_Taxa, colClasses = c('factor')) %>% filter(!duplicated(tax.id)))
    ## tax.id is not unique in Taxa! Unique is however needed for left_join by tax.id (not by inventory specific ids)!
  ),
  
  
  ## Env data
  # tar_load(starts_with("Env"))
  list(
    tar_target(predictor_select,
               c("phCaCl_esdacc", "waterLevel_loc")), # "alt_loc"
    
    tar_target(Waterlevel,
               data.frame(de = c("trocken", "mäßig trocken", "mäßig frisch", "frisch", "(mäßig feucht)", "feucht", "nass"), # "sehr nass" might not be in the data
                          en = c("dry", "slightly dry", "slightly damp", "damp", "(slightly moist)", "moist", "wet"), # "very wet" might not be in the data
                          value = c(-6, -4, -2, 0, 2, 4, 6))), ## highly acidic soils (pH<3.5); neutral = 7 ## 8 might not be in the data
    
    tar_target(Env_clean,
               cleanEnv(Data_env, c(predictor_select, "alt_loc"))),
    ## Summarize (mutate) predictors per cluster, but keep plot-level disturbance data etc.
    tar_target(Env_cluster,
               summarizeEnvByCluster(Env_clean, predictor_select))
  ),
  

  ## Stage abundance data
  list(
    tar_target(Data_big_area,
               prepareBigData(Data_big, Data_big_status,
                              taxon_select = taxon_select, threshold_dbh = threshold_dbh, radius_max = radius_max,
                              tablepath = dir_publish)),
    
    tar_target(Data_small_area,
               prepareSmallData(Data_small, taxon_select = taxon_select, regclass_select = regclass_select)),
    
    tar_target(Stages_transitions,
               countTransitions(Data_big, Data_big_status, Env_cluster, Stages_select,
                                taxon_select = taxon_select, threshold_dbh = threshold_dbh, radius_max = radius_max,
                                loc = c("plot", "nested", "cluster"))),
    
    tar_target(Stages,
               joinStages(Data_big_area, Data_small_area, taxon_select = taxon_select, threshold_dbh = threshold_dbh)),
    tar_target(Stages_env,
               joinEnv(Stages, Env_clean)), # Env_cluster ## always joins by plotid, even if Env is aggregated by cluster
    
    list(
      tar_target(taxon_s,
                 factor(c(sort(taxon_select), "other"))),
      tar_target(BA_s,
                 constructConstantGrid(taxon_s, Stages_env, Data_geo),
                 pattern = map(taxon_s),
                 iteration = "list"),
      tar_target(fits_s,
                 fitS(BA_s, path = dir_publish),
                 pattern = map(BA_s),
                 iteration = "list"),
      ## fits_s each have an attribute "taxon"
      tar_target(Stages_s,
                 predictS(fits_s, Stages_env),
                 iteration = "list"),
        ## Former workaround for Linux installlations, where the geo libraries did not work. Assumes upload of "Data/Stages_s.rds"
        # tar_target(file_Stages_s, if (!! onserver) "Data/Stages_s.rds" else saveStages_s(Stages_s), format = "file"),
        # tar_target(Data_Stages_s, readRDS(file_Stages_s)), ## this would be loaded in downstream targets instead of Stages_s
      
      tar_target(surfaces_s,
                 predictSurfaces(fits_s),
                 iteration = "list"),
      tar_target(surfaceplots_s,
                 plotSurfaces(surfaces_s, themefun = themefunction, path = dir_publish),
                 iteration = "list")
    ),
    
    ## The seedlings pipeline
    list(
      tar_target(Seedlings,
                 wrangleSeedlings(Data_seedlings, taxon_select = taxon_select, threshold_dbh = threshold_dbh)),
      tar_target(seedlings_fullgrid,
                 wrangleSeedlings_s(Data_seedlings_fullgrid, taxon_select = taxon_select, threshold_dbh = threshold_dbh),
                 iteration = "list"),
      tar_target(fits_Seedlings_s, ## fits_s each have an attribute "taxon"
                 fitS(seedlings_fullgrid, path = NULL), ##! setting path to NULL prevents overwriting of the summary of the gam
                 pattern = map(seedlings_fullgrid),
                 iteration = "list"),
      tar_target(surfaces_Seedlings_s,
                 predictSeedlingsSurfaces(fits_Seedlings_s),
                 iteration = "list"),
      tar_target(surfaceplots_Seedlings_s,
                 plotSurfaces(surfaces_Seedlings_s, path = dir_publish),
                 iteration = "list"),
      tar_target(Seedlings_s,
                 predictS(fits_Seedlings_s, Seedlings)),
      tar_target(file_model_Seedlings,
                 "Model_2023-02_bb/Model_seedlings.stan",
                 format = "file"),
      tar_target(model_Seedlings,
                 cmdstan_model(file_model_Seedlings, stanc_options = list("O1"))),
      tar_target(fit_Seedlings,
                 fitSeedlings(model_Seedlings, Seedlings_s, fitpath = dir_fit))
    ),
    
    
    tar_target(Stages_select,
               selectLocs(Stages_s, predictor_select,
                          selectspec = T, selectpred = F, stratpred = F, selectalt = c(100, 600), n_locations = n_locations, loc = "plot", tablepath = dir_publish)), # Subsetting after smooth, so that smooth can be informed by all plots.
    
    tar_target(Stages_scaled,
               scaleData(Stages_select, predictor_select)), # After selection, so that scaling includes selected plots .
    
    tar_target(Stages_loc,
               setLocLevel(Stages_scaled, Env_cluster, loc = c("plot", "nested", "cluster"))),

    
    ## Publishing
    tar_target(Summary_taxa,
               summarizeTaxa(Data_big, Data_seedlings, Stages_select, Seedlings_s, tablepath = dir_publish)),
    tar_target(Summary_NFIs,
               summarizeNFIs(Data_big, Data_seedlings, Stages_select, Seedlings_s, tablepath = dir_publish)),
    tar_target(map_select,
               mapLocations(Stages_select, path = dir_publish, themefun = themefunction),
               packages = c(package, "eurostat", "elevatr", "terrainr", "rayshader", "ggspatial", "elementalist"))
    
  )
)



## Range pipeline ------------------------------------------------------
targets_range <- list(
  
  ## Read data files
  list(
    tar_target(file_EAFTS_range,
               "Range.nosync/EAFTS/Fagus-sylvatica_rpp.tif",
               format = "file"),
    tar_target(file_env_range,
               "Inventory.nosync/DE BWI/Data/DE_BWI_big_2-3_status.rds",
               format = "file")
  ),
  list(
    tar_target(predictor_range, c("cwbYear_aclim" = "Range.nosync/Climate aclim/Climatologies/CWB_YEAR.tif",
                                  "phCaCl_esdacc" = "Range.nosync/Soil esdacc ESDAC topsoil chemical properties/pH_CaCl/pH_CaCl.tif")),
    
    tar_target(Raster_EAFTS_range, raster(file_EAFTS_range)),
    tar_target(rasters_env_range, lapply(predictor_range, raster))
  ),
  
  
  ## Compare ranges
  list(
    tar_target(Stages_select_range,
               selectRange(Stages_select_env, taxon = "Fagus.sylvatica")),
               
    tar_target(Stages_env_range,
              joinEnv(Stages_select_range,
                      st_drop_geometry(Data_env[ ,c("plotid", "clusterid", names(predictor_range))]))
              ),
    
    tar_target(EAFTS_range,
               extractRange(Raster_EAFTS_range, rasters_env = rasters_env_range)),
    
    tar_target(Ranges_range,
               bind_rows("DE" = st_drop_geometry(Stages_env_range), "EAFTS" = EAFTS_range, .id = "origin")),
    
    tar_target(Summary_range,
               summarizeRange(Ranges_range, predictor = names(predictor_range), path = dir_publish)),
    
    tar_target(plot_range,
               plotRange(Ranges_range, predictor = names(predictor_range), path = dir_publish, color = twocolors, themefun = themefunction))
  )
)



## Fitting pipelines ------------------------------------------------
#### general ----------
targets_fit_general <- list(
  tar_target(data_stan,
             formatStanData(Stages_loc, Stages_transitions, taxon_s, threshold_dbh, predictor_select,
                            loc = c("plot", "nested", "cluster"))),
  
  tar_target(data_stan_transitions,
             formatStanData(Stages_loc, Stages_transitions, taxon_s, threshold_dbh, predictor_select,
                            loc = "plot")),
  
  tar_target(file_model_transitions,
             "Model_2023-02_bb/Model_transitions.stan",
             format = "file"),
  
  tar_target(model_transitions,
             cmdstan_model(file_model_transitions, stanc_options = list("O1"))),
  
  tar_target(fit_g,
             fitTransition(data_stan_transitions, which = "g", model_transitions, fitpath = dir_fit)),
  
  tar_target(fit_h,
             fitTransition(data_stan_transitions, which = "h", model_transitions, fitpath = dir_fit)),
  
  tar_target(data_stan_priors,
             formatPriors(data_stan, weakpriors,
                          fit_g = NULL, fit_h = NULL, fit_Seedlings = fit_Seedlings,
                          widthfactor_trans = 1, widthfactor_reg = 4)),
  
  tar_target(offsetname,
             c("offset", "offset_avg", "offset_q1", "offset_q3")),
  tar_target(offsetname_select,
             offsetname[1]),
  tar_target(data_stan_priors_offset,
             selectOffset(offsetname_select, data_stan_priors)),
  tar_target(data_stan_priors_offsets,
             selectOffset(offsetname, data_stan_priors),
             pattern = map(offsetname), iteration = "list"),
  tar_target(data_stan_priors_offsets_1,
             data_stan_priors_offsets[1], iteration = "list")
)

#### fit_test ----------
targets_fit_test <- list(

  tar_target(file_model_test,
             "Model_2021-03_ba/Model_ba_test.stan",
             format = "file"),
  tar_target(model_test,
             cmdstan_model(file_model_test, stanc_options = list("O1")) #, cpp_options = list(stan_opencl = TRUE)
             ),
  ## Prior predictive tests that rely on currently out-commented generated quantities
  # tar_target(priorsim_test,
  #            drawTest(model = model_test, data_stan = data_stan_priors_offset, method = "sim", gpq = FALSE, fitpath = dir_fit)),
  # tar_target(plots_predictions_priorsim_test,
  #            plotPredictions(cmdstanfit = priorsim_test, data_stan_priors_offset, check = "prior")),

  tar_target(fit_test_sansgq,
             fitModel(model = model_test, data_stan = data_stan_priors_offset, gpq = FALSE,
                      method = "mcmc", n_chains = 4, iter_warmup = 800, iter_sampling = 500, fitpath = dir_fit)),
  tar_target(fit_test,
             fitModel(model = model_test, data_stan = data_stan_priors_offset, gpq = TRUE,
                      method = "mcmc", n_chains = 4, iter_warmup = 800, iter_sampling = 500, fitpath = dir_fit)),
  tar_target(basename_fit_test,
             getBaseName(fit_test))
)

#### fit ----------
targets_fit <- list(
  tar_target(file_model,
             "Model_2023-02_bb/Model_bb.stan",
             format = "file"),
  tar_target(model,
             cmdstan_model(file_model, stanc_options = list("O1"))),
  tar_target(fit,
             fitModel(model = model, data_stan = data_stan_priors_offsets, gpq = TRUE, ## use data_stan_priors_offsets for mapping over all offset variants; data_stan_priors_offsets_1 for only using the one default option
                      method = "mcmc", n_chains = 4, iter_warmup = 1000, iter_sampling = 1000, fitpath = dir_fit),
             pattern = map(data_stan_priors_offsets), iteration = "list"), ## data_stan_priors_offsets ## for testing all offsets
  tar_target(basename_fit,
             getBaseName(fit),
             pattern = map(fit), iteration = "list")
)

#### fit_env ----------
targets_fit_env <- list(
  
  ## Prepare data
  tar_target(Stages_select_env,
             selectLocs(Stages_s, predictor_select,
                        selectspec = F, selectpred = T, stratpred = T, selectalt = c(100, 600), selectseed = 0, n_locations = n_locations_env,
                        loc = "plot", tablepath = dir_publish)), # Selection based on whether environmental variables are present
  
  tar_target(Stages_scaled_env,
             scaleData(Stages_select_env, predictor_select)), # After selection, so that scaling includes selected plots 
  
  tar_target(Stages_loc_env,
             setLocLevel(Stages_scaled_env, Env_cluster, loc = c("plot", "nested", "cluster"))),
  
  tar_target(Stages_transitions_env,
             countTransitions(Data_big, Data_big_status, Env_cluster, Stages_select_env,
                              taxon_select = taxon_select, threshold_dbh = threshold_dbh, radius_max = radius_max,
                              loc = c("plot", "nested", "cluster"))),
  
  tar_target(data_stan_env,
             formatStanData(Stages_loc_env, Stages_transitions_env, taxon_s, threshold_dbh,  predictor_select,
                            loc = c("plot", "nested", "cluster"), smoothmediantax = "none")),
  
  tar_target(data_stan_transitions_env,
             formatStanData(Stages_loc, Stages_transitions_env, taxon_s, threshold_dbh, predictor_select,
                            loc = "plot")),
  
  tar_target(fit_g_env,
             fitTransition(data_stan_transitions_env, which = "g", model_transitions, fitpath = dir_fit)),
  
  tar_target(fit_h_env,
             fitTransition(data_stan_transitions_env, which = "h", model_transitions, fitpath = dir_fit)),
  
  tar_target(data_stan_priors_env,
             formatPriors(data_stan_env, weakpriors_env)),

  tar_target(data_stan_priors_offset_env,
             selectOffset(offsetname_select, data_stan_priors_env)),

  ## Fit
  tar_target(file_model_env_vertex,
             "Model_2023-02_bb/Model_bb_env_vertex.stan",
             format = "file"),
  tar_target(model_env,
             cmdstan_model(file_model_env_vertex, stanc_options = list("O1"))),
  tar_target(fit_env,
             fitModel(model = model_env, data_stan = data_stan_priors_offset_env, gpq = TRUE,
                      method = "mcmc", n_chains = 4, iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95, fitpath = dir_fit)
             ),
  tar_target(basename_fit_env,
             getBaseName(fit_env))
)



## Posterior pipeline ------------------------------------------------------
#### posterior_test -----------
targets_posterior_test <- list(
  
  ## Extract
  tar_target(stanfit_test,
             extractStanfit(cmdstanfit = fit_test)),
  tar_target(draws_test,
             extractDraws(stanfit = fit_test, exclude = helper_exclude)),
  # tar_target(stanfit_test_plotting,
  #            extractStanfit(fit_test_sansgq, purge = TRUE)),
  
  
  ## Summarize
  tar_target(summary_test,
             summarizeFit(cmdstanfit = fit_test, exclude = c(helper_exclude, rep_exclude, par_exclude, simname_prior, simname_posterior, parname_loc),
                          publishpar = parname_plotorder, path = dir_publish)),
  tar_target(summary_states_test,
             summarizeStates(States = States_test, data_stan = data_stan, basename = basename_fit_test, path = dir_publish)),
  tar_target(Freq_converged_test,
             summarizeFreqConverged(cmdstanfit = fit_test, data_stan_priors, path = dir_publish)),
  
  
  ## Generate
  tar_target(residuals_test,
             generateResiduals(cmdstanfit = fit_test, data_stan_priors, path = dir_publish)),
  # tar_target(Trajectories_test,
  #            generateTrajectories(cmdstanfit = fit_test, data_stan_priors, parname, locparname = parname_loc,
  #                                 time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 50, average = "none")),
  tar_target(Trajectories_avg_test,
             generateTrajectories(cmdstanfit = fit_test, data_stan_priors, parname, locparname = parname_loc,
                                  time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1, average = "locsperdraws_all")),
  # tar_target(Trajectories_quantiles_test,
  #            generateTrajectories(cmdstanfit = fit_test, data_stan_priors, parname, locparname = parname_loc,
  #                                 time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1, average = "locsperdraws_avgL_qInit")),
            ### average options:
            ## c("none",            # — no averaging. Paraneters in simulations vary per loc and draw
            ## "locsperdraws_all",  # — average all loc-wise parameters per draw
            ## "drawsperlocs_all",  # — average all parameters per loc, so that there is only one trajectory per cluster
            ## "locsperdraws_avgL", # — average only "L_loc" per draw, so that there are initial values that vary with cluster
            ## "locsperdraws_avgL_qInit", # — average only "L_loc" per draw, so that there are initial values that vary with cluster
  
  ## Formatted posterior data stuctures
  tar_target(States_test,
             formatStates(cmdstanfit = fit_test, statename = statename,
                          data_stan_priors = data_stan_priors)),
  
  ## Plot
  tar_target(plots_test,
             plotStanfit(stanfit = stanfit_test, exclude = exclude, path = dir_publish, basename = basename_fit_test, color = twocolors, themefun = themefunction)),
  tar_target(plots_parameters_test,
             plotParameters(draws = stanfit_test, parname = parname_plotorder, exclude = exclude, path = dir_publish, basename = basename_fit_test, color = twocolors, themefun = themefunction)),
  ## Prior predictive tests that rely on currently out-commented generated quantities
  # tar_target(plots_predictions_prior_test,
  #            plotPredictions(cmdstanfit = fit_test, data_stan_priors, check = "prior")),
  tar_target(plots_predictions_posterior_test,
             plotPredictions(cmdstanfit = fit_test, data_stan_priors_offset, check = "posterior", path = dir_publish)),
  tar_target(plots_trace_test,
             plotTrace(cmdstanfit = fit_test, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plots_pairs_test,
             plotPairs(cmdstanfit = fit_test, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plots_conditional_test,
             plotConditional(cmdstanfit = fit_test, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plot_contributions_test,
             plotContributions(cmdstanfit = fit_test, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, color = twocolors, themefun = themefunction)),
  # tar_target(plot_contributions_prop_test,
  #            plotContributions(cmdstanfit = fit_test, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, contribution = "sum_ko_prop", color = twocolors, themefun = themefunction)),
  tar_target(plots_states_test,
             plotStates(States_test, allstatevars = basalareaname,
                        path = dir_publish, basename = basename_fit_test, color = twocolors, themefun = themefunction)),
  tar_target(plot_trajectories_avg_test,
             plotTrajectories(Trajectories_avg_test, thicker = T, path = dir_publish, basename = basename_fit_test, color = twocolors, themefun = themefunction)),
  tar_target(animation_trajectories_avg_test,
             animateTrajectories(plot_trajectories_avg_test, path = dir_publish, basename = basename_fit_test)),
  
  # tar_target(plot_powerscale_test,
  #            plotSensitivity(cmdstanfit = fit_test, include = parname, path = dir_publish)),
  
  ## Test
  tar_target(sensitivity_test,
             testSensitivity(fit_test, include = parname, path = dir_publish))
)

#### posterior -----------
targets_posterior <- list(
  ## Extract
  tar_target(stanfit,
             extractStanfit(cmdstanfit = fit),
             pattern = map(fit), iteration = "list"),
  tar_target(draws,
             extractDraws(fit, exclude = helper_exclude),
             pattern = map(stanfit), iteration = "list"),
  
  ## Summarize
  tar_target(summary,
             summarizeFit(cmdstanfit = fit, exclude = exclude,
                          publishpar = parname_lim_plotorder, path = dir_publish),
             pattern = map(fit), iteration = "list"),
  tar_target(summary_states,
             summarizeStates(States = States, data_stan = data_stan, basename = basename_fit, path = dir_publish),
             pattern = map(States, basename_fit), iteration = "list"),
  tar_target(Freq_converged,
             summarizeFreqConverged(cmdstanfit = fit, data_stan_priors, path = dir_publish),
             pattern = map(fit), iteration = "list"),
  
  ## Generate
  tar_target(residuals,
             generateResiduals(cmdstanfit = fit, data_stan_priors, path = dir_publish),
             pattern = map(fit), iteration = "list"),
  tar_target(Trajectories_avg,
             generateTrajectories(cmdstanfit = fit, data_stan_priors, parname, locparname = parname_loc,
                                  time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1, average = "locsperdraws_all"),
             pattern = map(fit), iteration = "list"),
  
  ## Formatted posterior data stuctures
  tar_target(States,
             formatStates(cmdstanfit = fit, statename = statename, data_stan_priors = data_stan_priors),
             pattern = map(fit), iteration = "list"),
  
  ## Plot
  tar_target(plots,
             plotStanfit(stanfit = stanfit, exclude = exclude, path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(stanfit, basename_fit), iteration = "list"),
  tar_target(plots_parameters,
             plotParameters(draws = stanfit, parname = parname_plotorder, exclude = exclude, supp = FALSE, path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(stanfit, basename_fit), iteration = "list"),
  tar_target(plots_parameters_lim,
             plotParameters(draws = stanfit, parname = parname_lim_plotorder_plot, exclude = exclude, supp = TRUE, path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(stanfit, basename_fit), iteration = "list"),
  tar_target(plots_trace,
             plotTrace(cmdstanfit = fit, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction),
             pattern = map(fit, basename_fit), iteration = "list"),
  tar_target(plots_pairs,
             # plotPairs(cmdstanfit = fit, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction),
             plotConditional(cmdstanfit = fit, parname = parname_plotorder, conditional = F, path = dir_publish, themefun = themefunction),
             pattern = map(fit, basename_fit), iteration = "list"),
  tar_target(plots_conditional,
             plotConditional(cmdstanfit = fit, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction),
             pattern = map(fit), iteration = "list"),
  
  # tar_target(plot_contributions,
  #            plotContributions(cmdstanfit = fit, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, plotlog = T, color = twocolors, themefun = themefunction),
  #            pattern = map(fit), iteration = "list"),
  # tar_target(plot_contributions_supp,
  #            plotContributions(cmdstanfit = fit, parname = c(parname_plotorder[1:7], b_c_b = "b_c_b_log"), path = dir_publish, plotlog = T, color = twocolors, themefun = themefunction),
  #            pattern = map(fit), iteration = "list"),
  # tar_target(plot_contributions_prop,
  #            plotContributions(cmdstanfit = fit, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, contribution = "sum_ko_prop", plotlog = F, color = twocolors, themefun = themefunction),
  #            pattern = map(fit), iteration = "list"),
  # tar_target(plot_contributions_switch,
  #            plotContributions(cmdstanfit = fit, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, contribution = "sum_switch", plotlog = T, color = twocolors, themefun = themefunction),
  #            pattern = map(fit), iteration = "list"),
  
  tar_target(plots_states,
             plotStates(States,
                        mainstatevars = c("ba_init", "ba_fix", "ba_fix_switch_s"),
                        suppstatevars = c("ba_fix", "ba_fix_switch_l_r", "ba_fix_switch_g_s", "ba_fix_switch_c_j", "ba_fix_switch_c_a", "ba_fix_switch_b_c_b"),
                        allstatevars = basalareaname,
                        path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(States, basename_fit), iteration = "list"),
  tar_target(plot_predominant,
             plotPredominant(States, majorname_main,
                             path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(States, basename_fit), iteration = "list"),
  tar_target(plot_predominant_supp,
             plotPredominant(States, majorname_supp, supp = T,
                             path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction),
             pattern = map(States, basename_fit), iteration = "list"),
  tar_target(plot_trajectories_avg,
             plotTrajectories(Trajectories_avg, thicker = T, path = dir_publish, basename = basename_fit, plotpdf = TRUE, color = twocolors, themefun = themefunction),
             pattern = map(Trajectories_avg, basename_fit), iteration = "list"),
  tar_target(animation_trajectories_avg,
             animateTrajectories(plot_trajectories_avg[[1]], path = dir_publish, basename = basename_fit[[1]]))
)

#### posterior_env -----------
targets_posterior_env <- list(

  ## Extract
  
  # tar_target(stanfit_env,
  #            extractStanfit(cmdstanfit = fit_env)),
  
  tar_target(draws_env,
             extractDraws(fit_env, exclude = exclude)),
  
  ## Summarize
  tar_target(Summary_NFIs_env,
             summarizeNFIs(Data_big, Data_seedlings, Stages_select_env, Seedlings_s, tablepath = dir_publish)),
  tar_target(Summary_taxa_env,
             summarizeTaxa(Data_big, Data_seedlings, Stages_select_env, Seedlings_s, tablepath = dir_publish)),
  tar_target(summary_env,
             summarizeFit(cmdstanfit = fit_env, exclude = exclude,
                          publishpar = c(parname_lim_plotorder, parname_vertex_env, parname_hyper), path = dir_publish)),
  tar_target(summary_states_env,
             summarizeStates(States = States_env, data_stan = data_stan_env, basename = basename_fit_env, path = dir_publish)),
  tar_target(summary_marginal_env,
             summarizeMarginal(Marginal = Marginal_env, basename = basename_fit_env, path = dir_publish)),
  tar_target(Freq_converged_env,
             summarizeFreqConverged(cmdstanfit = fit_env, data_stan_priors_env, path = dir_publish)),

  ## Generate 
  tar_target(residuals_env,
             generateResiduals(cmdstanfit = fit_env, data_stan_priors_env, includeinit = F, path = dir_publish)),
  tar_target(Trajectories_avg_env,
             generateTrajectories(cmdstanfit = fit_env, data_stan_priors_env, parname, locparname = parname_trajectories_environmental_env,
                                  time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1, average = "drawsperlocs_all")),

  
  ## Formatted posterior data stuctures
  tar_target(States_env,
             formatStates(cmdstanfit = fit_env, statename = statename_select, data_stan_priors = data_stan_priors_env)),
  
  tar_target(Environmental_env,
             formatEnvironmental(cmdstanfit = fit_env, parname = parname_environmental_env, ## formatEnvironmental() always includes major_fix in the output variables, regardless of parname!
                                 data_stan = data_stan_priors_offset_env, envname = predictor_select, locmeans = F, jitter = 0.1)),
  tar_target(Diff_environmentalko_env,
             formatEnvironmental(cmdstanfit = fit_env, parname = parname_environmental_diff_env,
                                 data_stan = data_stan_priors_offset_env, envname = predictor_select, locmeans = F)),
  tar_target(Cred_env,
             formatCred(Diff_environmentalko_env, envname = predictor_select, credlevel = 0.9, basename = basename_fit_env, path = dir_publish)),
  
  tar_target(Envgrid_env,
             formatEnvgrid(data_stan_priors_offset_env, envname = predictor_select, res = 500)),
  tar_target(Marginal_env,
             formatMarginal(Environmental_env)),

  
  ## Post-hoc inference
  tar_target(Surfaces_poly_env,
             predictPoly(cmdstanfit = fit_env, parname_Beta = parname_Beta_env, Envgrid = Envgrid_env)),
  
  
  tar_target(alltaxa_enironmental_env, c("Fagus.sylvatica" = 1, "others" = 2, "both" = 0)),
  tar_target(comb_taxa_enironmental_binomial_env, c(0, 0, 1, 1, 2, 0)),
  tar_target(comb_par_enironmental_binomial_env, c("major_init", "major_fix", "ba_frac_init", "ba_frac_fix", "ba_frac_fix", "ba_frac_fix_other_s")),

  tar_target(fit_environmental_gaussian_env,
             fitEnvironmental_glm(Environmental_env, parname = parname_environmental_gaussian_env, envname = predictor_select, taxon = taxon_s, fam = "gaussian", path = dir_publish),
             pattern = cross(parname_environmental_gaussian_env, taxon_s),
             iteration = "list"),
  tar_target(fit_environmental_ba_env,
             fitEnvironmental_gam(Environmental_env, parname = parname_environmental_ba_select_env, envname = predictor_select, taxon = taxon_s, fam = "gaussian", path = dir_publish),
             pattern = cross(parname_environmental_ba_select_env, taxon_s),
             iteration = "list"),
  tar_target(fit_environmental_diff_env,
             fitEnvironmental_gam(Diff_environmentalko_env, parname = parname_environmental_diff_env, envname = predictor_select, taxon = 0, fam = "gaussian", path = dir_publish),
             pattern = map(parname_environmental_diff_env),
             iteration = "list"),
  tar_target(fit_environmental_binomial_env,
             fitEnvironmental_gam(Environmental_env, parname = comb_par_enironmental_binomial_env, envname = predictor_select, taxon = comb_taxa_enironmental_binomial_env, fam = "binomial", path = dir_publish),
             pattern = map(comb_par_enironmental_binomial_env, comb_taxa_enironmental_binomial_env), ## more efficient manual version, avoids superfical loading of Environmental_env
             # pattern = cross(parname_environmental_binomial_env, alltaxa_enironmental_env), ## although this is crossing alltaxa 0:2, only the cases in the data (either 0, e.g. for major_fix, or 1 and 2 for frac) will be used
             iteration = "list"),
  tar_target(fit_environmental_par_env,
             fitEnvironmental_gam(Environmental_env, parname = parname_environmental_par_env, envname = predictor_select, taxon = taxon_s, fam = "gaussian", path = dir_publish),
             pattern = cross(parname_environmental_par_env, taxon_s),
             iteration = "list"),
  tar_target(fit_environmental_env,
             c(fit_environmental_binomial_env, fit_environmental_ba_env), # excluded: fit_environmental_gaussian_env, i.e. all the pure parameters
             iteration = "list"),
  
  tar_target(surface_environmental_env,
             predictEnvironmental(fit_environmental_env, envname = predictor_select, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
             pattern = map(fit_environmental_env),
             iteration = "list"),

  tar_target(Surface_binary_env,
             surface_environmental_env[[ matchAttr(surface_environmental_env, "parname", "ba_frac_fix") ]]), ## major_fix
  
  tar_target(surface_par_env,
             predictEnvironmental(fit_environmental_par_env, envname = predictor_select, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
             pattern = map(fit_environmental_par_env),
             iteration = "list"),
  
  tar_target(surface_diff_env,
             predictEnvironmental(fit_environmental_diff_env, envname = predictor_select, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
             pattern = map(fit_environmental_diff_env),
             iteration = "list"),
  
  tar_target(surface_diff_interp_env,
             interpolateSurface(Diff_environmentalko_env, parname = parname_environmental_diff_env, aggr = "mean", envname = predictor_select, avg = c("env", "loc"))),
  
  tar_target(surface_diff_interp_sd_env,
             interpolateSurface(Diff_environmentalko_env, parname = parname_environmental_diff_env, aggr = "sd", envname = predictor_select, avg = c("env", "loc"))),

  # tar_target(surface_environmental_env_gaussian,
  #            predictEnvironmental(fit_environmental_gaussian_env, envname = predictor_select,
  #                                 path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
  #            pattern = map(fit_environmental_gaussian_env),
  #            iteration = "list"),
  # tar_target(surface_environmental_env_ba,
  #            predictEnvironmental(fit_environmental_ba_env, envname = predictor_select,
  #                                 path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
  #            pattern = map(fit_environmental_ba_env),
  #            iteration = "list"),
  # tar_target(surface_environmental_env_binomial,
  #            predictEnvironmental(fit_environmental_binomial_env, envname = predictor_select,
  #                                 path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction),
  #            pattern = map(fit_environmental_binomial_env),
  #            iteration = "list"),

  
  ## Plot
  tar_target(plot_poly_env,
             plotPoly(Surfaces_poly_env,
                      Environmental = NULL, # Environmental_env,
                      Binary = NULL, # = Surface_binary_env
                      Waterlevel = Waterlevel,
                      basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = theme_fagus)),

  tar_target(plot_triptych_env,
             plotTriptych(Environmental_env,
                          Surface_init = surface_environmental_env[[ matchAttr(surface_environmental_env, "parname", "ba_frac_init") ]],
                          Surface_fix = surface_environmental_env[[ matchAttr(surface_environmental_env, "parname", "ba_frac_fix") ]],
                          Binary = Surface_binary_env,
                          Summary = summary_env,
                          Scaling = Stages_loc_env,
                          Waterlevel = Waterlevel,
                          basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  
  tar_target(plot_environmental_env,
             plotEnvironmental(surfaces = surface_environmental_env, binaryname = "ba_frac_fix", Waterlevel = Waterlevel, commonscale = F, removevar = c("major_init", "major_fix", "ba_frac_fix[2]"), ## binaryname = "major_fix"
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  
  tar_target(plot_environmental_lim_env,
             plotEnvironmental(surfaces = surface_par_env[c(13, 3, 9, 14, 4, 10, 15, 5, 11, 16, 6, 12, 17, 1, 7, 18, 2, 8)],
                               binaryname = NULL, Waterlevel = Waterlevel, commonscale = F, gridcols = 3,
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  tar_target(plot_environmental_L_env,
             plotEnvironmental(surfaces = surface_par_env[19:20],
                               binaryname = NULL, Waterlevel = Waterlevel, commonscale = F, gridcols = 3,
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  
  tar_target(plot_diff_env,
             plotEnvironmental(surfaces = c(surface_diff_env[c(r = 3:4, c_J = 5:6, s = 7:8, b = 15:16, c_b = 17:18)], list(Binary = Surface_binary_env)), binaryname = "ba_frac_fix", Waterlevel = Waterlevel, Cred = Cred_env, commonscale = T, removevar = c("ba_frac_init", "ba_frac_fix[2]", "major_init", "major_fix"),
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  
  tar_target(plot_diff_supp_env,
             plotEnvironmental(surfaces = c(surface_diff_env[3:18], list(Binary = Surface_binary_env)), binaryname = "ba_frac_fix", Waterlevel = Waterlevel, Cred = Cred_env, commonscale = T, removevar = c("ba_frac_init", "ba_frac_fix[2]", "major_init", "major_fix"),
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  
  tar_target(plot_diff_lim_env,
             plotEnvironmental(surfaces = c(surface_diff_env[19:24], list(Binary = Surface_binary_env)), binaryname = "ba_frac_fix", Waterlevel = Waterlevel, Cred = Cred_env, commonscale = T, removevar = c("ba_frac_init", "ba_frac_fix[2]", "major_init", "major_fix"),
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  tar_target(plot_diff_L_env,
             plotEnvironmental(surfaces = c(surface_diff_env[1:2], list(Binary = Surface_binary_env)), binaryname = "ba_frac_fix", Waterlevel = Waterlevel, Cred = Cred_env, commonscale = T, removevar = c("ba_frac_init", "ba_frac_fix", "major_init"),
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  tar_target(plot_diff_interp_env,
             plotEnvironmental(surfaces = surface_diff_interp_env, binaryname = "ba_frac_fix", Waterlevel = Waterlevel, commonscale = T, removevar = c("ba_frac_init", "ba_frac_fix[2]", "major_init", "major_fix"),
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),
  tar_target(plot_diff_interp_sd_env,
             plotEnvironmental(surfaces = surface_diff_interp_sd_env, binaryname = "ba_frac_fix", Waterlevel = NULL, commonscale = T, removevar = c("ba_frac_init", "ba_frac_fix[2]", "major_init", "major_fix"), ## interpolated surfaces do not have both scaled and unscaled predictors, therefore Waterlevel matching does not work 
                               basename = basename_fit_env, path = dir_publish, color = twocolors, ps = plotsettings, themefun = themefunction)),

  tar_target(plot_predominant_env,
             plotPredominant(States_env, majorname = c("major_init", "major_fix"),
                             path = dir_publish, basename = basename_fit, color = twocolors, themefun = themefunction)),
  
  tar_target(plot_binary_par_env,
             plotBinary(Environmental = Environmental_env,
                        parname = str_to_sentence(parname_plotorder),
                        fit_bin = fit_environmental_binomial_env[[ matchAttr(surface_environmental_env, "par", "major_fix") ]],
                        binarythreshold = 0.5, facetscale = "free", plotlog = F,
                        path = dir_publish, basename = basename_fit_env,  color = twocolors, themefun = themefunction)),
  
  tar_target(plot_binary_par_classic_env, ## binary division of the posterior based on majority
             plotBinary(Environmental = Environmental_env,
                        parname = str_to_sentence(parname_plotorder),
                        fit_bin = NULL,
                        binarythreshold = 0.5, facetscale = "free", plotlog = F,
                        path = dir_publish, basename = basename_fit_env,  color = twocolors, themefun = themefunction)),
  
  # tar_target(plot_binary_contrib_env,
  #            plotBinary(Environmental = Environmental_env,
  #                       parname = c(contribname_init_env), ## + contribname_counterfactual_env
  #                       fit_bin = fit_environmental_binomial_env[[ matchAttr(surface_environmental_env, "par", "major_fix") ]],
  #                       binarythreshold = 0.5, facetscale = "free_x", plotlog = T, ## !!!
  #                       path = dir_publish, basename = basename_fit_env,  color = twocolors, themefun = themefunction)),
  
  tar_target(plots_trace_env,
             plotTrace(cmdstanfit = fit_env, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  
  tar_target(plots_pairs_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  tar_target(plots_pairs_J_Fagus_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_J_env, formatparname = T, selecttax = 1)),
  tar_target(plots_pairs_J_others_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_J_env, formatparname = T, selecttax = 2)),
  tar_target(plots_pairs_A_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_A_env, formatparname = T)),
  tar_target(plots_pairs_J_env,
             plotPairs(cmdstanfit = fit_env,
                       parname =  c("r_log", "l_log", "s_log",
                                     "c_j_log", "c_j_center_env1", "c_j_center_env2", "c_j_log_spread_env1", "c_j_log_spread_env2",
                                     "g_log", "g_log_center_env1", "g_log_center_env2", "g_log_spread_env1", "g_log_spread_env2"),
                       formatparname = T)),
  tar_target(plots_pairs_B_Fagus_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_B_env, formatparname = T, selecttax = 1)),
  tar_target(plots_pairs_B_others_env,
             plotPairs(cmdstanfit = fit_env, parname = parname_B_env, formatparname = T, selecttax = 2)),
  
  tar_target(plots_parameters_env,
             plotParameters(draws = draws_env, parname = parname_plotorder, exclude = exclude, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction)),
  tar_target(plot_posterior_center_env,
             plotPosterior(cmdstanfit = fit_env, varname = parname_center_env)),
  tar_target(plot_posterior_spread_env,
             plotPosterior(cmdstanfit = fit_env, varname = parname_spread_env)),
  tar_target(plot_posterior_phi_env,
             plotPosterior(cmdstanfit = fit_env, varname = "phi_obs_inv_sqrt")),
  
  tar_target(plot_marginal_env,
             plotMarginal(Marginal = Marginal_env, parname = str_to_sentence(parname_plotorder),
                          path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = theme_fagus)),
  
  # tar_target(plot_contributions_env,
  #            plotContributions(cmdstanfit = fit_env, parname = c(parname_plotorder, b_c_b = "b_c_b_log"), path = dir_publish, plotlog = T,
  #                              color = twocolors, themefun = themefunction)),
  
  tar_target(plots_states_env,
             plotStates(States_env,
                        mainstatevars = c("ba_init", "ba_fix", "ba_fix_switch_s"),
                        suppstatevars = c("ba_fix", "ba_fix_switch_s", "ba_fix_ko_b_l_r"),
                        allstatevars = c("ba_init", "ba_fix", "ba_fix_ko_b_l_r"),
                        path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction)),
  
  tar_target(plot_trajectories_avg_env,
             plotTrajectories(Trajectories_avg_env, thicker = T, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction))
  
  # tar_target(plots_env,
  #            plotStanfit(stanfit = stanfit_env, exclude = exclude, path = dir_publish, basename = basename_fit_env, color = twocolors, themefun = themefunction)),
  # tar_target(plots_predictions_posterior_env,
  #            plotPredictions(cmdstanfit = fit_env, data_stan_priors_offset_env, check = "posterior", path = dir_publish)),
  # tar_target(plots_conditional_env,
  #            plotConditional(cmdstanfit = fit_env, parname = parname_plotorder, path = dir_publish, color = twocolors, themefun = themefunction)),
  # tar_target(plot_contributions_prop_env,
  #            plotContributions(cmdstanfit = fit_env, parname = parname_plotorder, path = dir_publish, contribution = "sum_ko_prop", color = twocolors, themefun = themefunction)),
  # tar_target(plot_contributions_switch_env,
  #            plotContributions(cmdstanfit = fit_env, parname = parname_plotorder, plotlog = T, path = dir_publish, contribution = "sum_switch", color = twocolors, themefun = themefunction)),
  # tar_target(animation_trajectories_avg_env,
  #            animateTrajectories(plot_trajectories_avg_env, path = dir_publish, basename = basename_fit_env))
)



## Standalone generated quantities pipeline -----------------------------------
# targets_sgq <- list(
#   
#   tar_target(file_gq,
#              paste0(tools::file_path_sans_ext(file_model),"_gq.stan"), format = "file"),
#   tar_target(model_gq,
#              cmdstan_model(file_gq)),
#   tar_target(fit_gq,
#              model_gq$generate_quantities(fitted_params = fit$output_files(),
#                                           data = data_stan_priors,
#                                           output_dir = dir_fit,
#                                           parallel_chains = getOption("mc.cores", 4)))
#   
#   # tar_target(rstanfit_gq_test,
#   #            extractStanfit(fit_gq_test)),
#   # tar_target(draws_gq_test,
#   #            extractDraws(rstanfit_gq_test, exclude = helper_exclude)),
#   
# )



# ————————————————————————————————————————————————————————————————————————————————— #
# Outer pipeline -----------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

list(
  targets_settings,
  targets_paths,
  targets_parname,
  targets_wrangling,
  targets_range,
  targets_fit_general,
  targets_fit_test,
  targets_fit,
  targets_fit_env,
  targets_posterior_test,
  targets_posterior,
  targets_posterior_env
)

