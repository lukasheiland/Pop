# ————————————————————————————————————————————————————————————————————————————————— #
# Posterior verbs         -----------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## extract — read draws. Returns a draws, or a stanfit object.
## format — formats the extracted draws for later use.
## summarize — tabulating the posterior, Returns a data.framoid
## generate — generate quantities from the posterior.
## test — statistical tests on the posteriors.
## plot - plots. Side effect: file. Returns a plot object or a list of plot objects.
## animate - animates a plot object. Side effect: file. Returns a gganimate object


# ————————————————————————————————————————————————————————————————————————————————— #
# Extract posterior         -----------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## extractStanfit --------------------------------
# cmdstanfit  <- tar_read("fit_test")
## stanfit <- extractStanfit(fit_test, purge = T)
extractStanfit <- function(cmdstanfit, purge = FALSE) {
  
  outputfile <- cmdstanfit$output_files()
  stanfit <- rstan::read_stan_csv(outputfile)
  
  ## This is to purge draws with NaN values from the model for plotting. NaNs can arise in generated quantities ...
  if (purge) {
    draws <- stanfit@sim$samples
    completerow <- lapply(draws, function(d) complete.cases(d)) %>% as.data.frame() %>% apply(1, all)
    stanfit@sim$samples <- lapply(draws, function(D) {attr(D, "sampler_params") <-  attr(D, "sampler_params")[completerow,]; D[completerow,] })
    
    n_draws_complete <- sum(completerow)
    n_draws <- stanfit@sim$n_save
    diff_complete <- stanfit@sim$n_save - n_draws_complete
    stanfit@sim$n_save <- rep(n_draws_complete, stanfit@sim$chains)
    stanfit@sim$iter <- stanfit@sim$iter - diff_complete[1]
    stanfit@sim$permutation <- lapply(stanfit@sim$permutation, function(p) p[!(p %in% n_draws[1]:n_draws_complete)])
    
    message("There were ", diff_complete[1], " draws with NaNs, that were purged from the fit.")
  }
  
  return(stanfit)
}


## extractDraws --------------------------------
# fit <- tar_read(fit_env)
# fit  <- tar_read("stanfit_test")
# exclude  <- tar_read(exclude)
extractDraws <- function(fit, exclude = tar_read(exclude)) {
  
  if(is(fit, "CmdStanMCMC")) {
    
    varname <- fit$metadata()$stan_variables
    include <- setdiff(varname, exclude)
    draws <- fit$draws(variables = varname)
    
  } else if(is(fit, "stanfit")) {
    
    draws <- rstan::extract(fit, pars = exclude, include = F)
  }
  
  return(draws)
}


## formatLoc --------------------------------
## helper for all variables that have structure array[N_locs] int/real
# name <- tar_read("parname_loc")
# cmdstanfit_ <- tar_read("fit_test")
# data_stan_priors_ <- tar_read("data_stan_priors")
formatLoc <- function(name, locmeans = FALSE, cmdstanfit_ = cmdstanfit, data_stan_priors_ = data_stan_priors) {
  n_locs <- data_stan_priors_$N_locs
  
  Draws <- cmdstanfit_$draws(variables = name, format = "draws_matrix")
  n_draws <- if (locmeans) 1 else nrow(Draws)
  n_i <- cmdstanfit_$metadata()$stan_variable_dims[[name]][2] ## get the length of the second dimension of the data structure called by 'name'
  if(is.null(n_i)) n_i <- cmdstanfit_$metadata()$stan_variable_sizes[[name]][2] ## fallback for old cmdstanr? versions. Indexing non-existing list element returns NULL.
  if(is.na(n_i)) n_i <- 1 ## case: has only one dimension. Indexing a vector out of bounds returns NA here.
  
  ## rep party
  loc <- rep(rep(1:n_locs, n_i), each = n_draws)
  i <- rep(rep(1:n_i, each = n_locs), each = n_draws)
  draw <- rep(1:n_draws, n_i * n_locs)
  
  Draws %<>%
    as.matrix() %>%
    { if (locmeans) colMeans(.) else c(.) } %>%
    set_names(NULL) %>%
    data.frame(value = ., loc = loc, i = i, draw = draw, var = name) %>%
    { if (locmeans) dplyr::select(., -i, -draw) else . }
  
  return(Draws)
}


## formatEnvironmental --------------------------------
# cmdstanfit <- tar_read("fit_env")
# parname <- tar_read("parname_environmental_env")
# envname <- tar_read("predictor_select")
# data_stan <- tar_read("data_stan_priors_offset_env")
formatEnvironmental <- function(cmdstanfit, parname = parname_env, data_stan = data_stan_priors_offset,
                                envname = predictor_select, locmeans = F) {
  
  varname_draws <- cmdstanfit$metadata()$stan_variables
  parname <- unique(c(parname, "major_fix")) ## explicitly include major fix, because it is assumed later
  parname <- intersect(parname, varname_draws)
  
  onlyreals <- all(str_detect(parname, "(ba_frac_.*fix_ko_.+)|(major_.+)"))
  anyenvdiff <- any(str_detect(parname, "ba_frac_.*fix_ko_.+"))
  
  draws_env <- cmdstanfit$draws(parname) %>%
    posterior::as_draws()
  
  Draws_env_bin <- tidybayes::gather_draws(draws_env, `major_.+`[loc], regex = T) %>% 
    { if(anyenvdiff) bind_rows(., tidybayes::gather_draws(draws_env, `ba_frac_.*fix_ko_.+`[loc], regex = T)) else . } %>% 
    bind_cols(tax = 0) %>%
    suppressWarnings() ## package tidyr warns about using deprecated gather_()
  
  if(onlyreals) {
    Draws_env <- Draws_env_bin
  } else {
    Draws_env <- tidybayes::gather_draws(draws_env, `.*`[loc,tax], regex = T) %>%
      suppressWarnings() %>% ## package tidyr warns about using deprecated gather_()
      bind_rows(Draws_env_bin)
  }

  
  if("L_loc" %in% parname) {
    
    Draws_env_l <- dplyr::filter(Draws_env, .variable == "L_loc") %>%
      mutate(.value = log(.value)) %>%
      mutate(.variable = "L_log")
    
    Draws_env %<>% bind_rows(Draws_env_l)
    
  }
  
  if (isTRUE(locmeans)) {
    Draws_env %<>%
      group_by(tax, loc, .variable) %>%
      summarize(.value = mean(.value, na.rm = T),
                sd = sd(.value, na.rm = T),
                p10 = quantile(.value, probs = 0.1, na.rm = T),
                p25 = quantile(.value, probs = 0.25, na.rm = T),
                median = median(.value, na.rm = T),
                p75 = quantile(.value, probs = 0.75, na.rm = T),
                p90 = quantile(.value, probs = 0.9, na.rm = T)) %>%
      ungroup()
  }
  
  envname_all <- c(envname, paste0(envname, "_s"))
  
  ## Add environmental values by loc
  Env <- attr(data_stan, "Long") %>% ## use "Long_BA" for an sf with point coordinates
    group_by(loc) %>%
    summarize_at(envname_all, function(x) first(x[!is.na(x)])) %>%
    ungroup()
  
  n_row <- nrow(Draws_env)
  Draws_env <- bind_cols(Draws_env, Env[match(Draws_env$loc, Env$loc), envname_all]) %>%
    filter(!(.variable %in% c("major_fix") & .value == 9)) %>%
    ungroup()
  
  n_dropped <- n_row - nrow(Draws_env)
  if (n_dropped > 0) warning("There were ", n_dropped ," draws*variables dropped, because not in all iterations the simulation has converged to the the fix point (equilibrium).")
  
  return(Draws_env)
}


## formatMarginal --------------------------------
# Environmental <- tar_read("Environmental_env")
formatMarginal <- function(Environmental) {
  
  D <- Environmental %>%
    group_by(.chain, .iteration, .draw, .variable, tax) %>%
    summarize(.value = mean(.value, na.rm = T), median = median(.value, na.rm = T), .groups = "drop")
  
  return(D)
}



## formatStates --------------------------------
# cmdstanfit <- tar_read("fit_test")
# statename <- tar_read("statename")
# data_stan_priors <- tar_read("data_stan_priors")
formatStates <- function(cmdstanfit, statename, data_stan_priors) {
  
  varname_draws <- cmdstanfit$metadata()$stan_variables
  varname <- intersect(statename, varname_draws)
  majorname <- varname[str_starts(varname, "major")]
  
  States <- lapply(varname, formatLoc, cmdstanfit_ = cmdstanfit, data_stan_priors_ = data_stan_priors)
  States <- lapply(States, function(S) if( length(unique(S$i)) == 1 ) bind_rows(S, within(S, {i <- 2})) else S )
  States %<>%
    bind_rows() %>%
    mutate(tax = factor(c("Fagus", "other")[i]))
  
  States$value[States$value == 9 & States$var %in% majorname] <- NA
  # States$value[States$value == 0 & !(States$var %in% majorname)] <- NA
  
  Quantiles <- filter(States, var == "ba_init") %>%
    group_by(draw, tax) %>%
    mutate(avg_ba_init = mean(value, na.rm = T),
           median_ba_init = quantile(value, prob = 0.5, type = 1, na.rm = T),
           p10_ba_init = quantile(value, prob = 0.1, type = 1, na.rm = T),
           p90_ba_init = quantile(value, prob = 0.9, type = 1, na.rm = T)) %>%
    summarize(loc_median_draw = first(loc[value == median_ba_init]), ## just in case that there might be more than 1, which is currently not the case
              loc_p10_draw = first(loc[value == p10_ba_init]),
              loc_p90_draw = first(loc[value == p90_ba_init])
    ) %>%
    ## get the loc that is most frequently the median for both taxa
    group_by(tax) %>%
    mutate(loc_median = first(sort(table(loc_median_draw), decreasing = T)),
           loc_p10 = first(sort(table(loc_p10_draw), decreasing = T)),
           loc_p90 = first(sort(table(loc_p90_draw), decreasing = T)))
  
  
  # implement if used: match by draw and tax
  States %<>%
    left_join(Quantiles, by = c("draw", "tax")) %>%
    mutate(is_loc_median = loc == loc_median,
           is_loc_p10 = loc == loc_p10,
           is_loc_p90 = loc == loc_p90,
           # is_loc_avg_draw = loc == loc_avg_draw,
           is_loc_median_draw = loc == loc_median_draw,
           is_loc_p10_draw = loc == loc_p10_draw,
           is_loc_p90_draw = loc == loc_p90_draw)
  
  return(States)
}


## predictPoly --------------------------------
# cmdstanfit <- tar_read("fit_env")
# parname_Beta <- tar_read("parname_Beta")
# Envgrid <- tar_read("Envgrid_env")
# data_stan_priors <- tar_read("data_stan_priors")
predictPoly <- function(cmdstanfit, parname_Beta, Envgrid) {
  
  varname_draws <- cmdstanfit$metadata()$stan_variables
  varname <- intersect(parname_Beta, varname_draws)
  
  draws <- cmdstanfit$draws(variables = parname_Beta) %>%
    as_draws_rvars()
  
  envname <- names(Envgrid)
  formula <- paste0("~ 1 + ", paste0("poly(", envname, ", 2, raw = T)" , collapse = " + "))
  Matrix <- model.matrix(as.formula(formula), data = Envgrid)
  dimnames(Matrix) <- NULL ## character dimnames were fatal for %**%

  predictGrid <- function(M, Beta) {
    
    Beta1 <- as_draws_matrix(Beta[, 1, drop = T])
    Beta2 <- as_draws_matrix(Beta[, 2, drop = T])
    
    ## %**% matrix multiplication from rvars does not seem to work with the called function tensorA::mul.tensor()
    mult <- function(A, b) {
      return(A %*% b)
    }
    
    p1 <- apply(Beta1, MARGIN = 1, mult, A = M)
    p1_mean <- apply(p1, MARGIN = 1, mean)
    p1_sd <- apply(p1, MARGIN = 1, sd)
    
    p2 <- apply(Beta2, MARGIN = 1, mult, A = M)
    p2_mean <- apply(p2, MARGIN = 1, mean)
    p2_sd <- apply(p2, MARGIN = 1, sd)
    
    ## uncertainties could also be carried on, to ggplot
    P1 <- data.frame(Envgrid, z = p1_mean, sd_z = p1_sd, taxon = "Fagus")
    P2 <- data.frame(Envgrid, z = p2_mean, sd_z = p2_sd, taxon = "others")
    
    P <- bind_rows(P1, P2)

    return(P)
  }

  predictions <- mclapply(draws, predictGrid, M = Matrix, mc.cores = length(draws))
  names(predictions) <- str_remove(names(draws), "Beta_")
  Predictions <- bind_rows(predictions, .id = "parname")
  
  attr(Predictions, "name_x") <- attr(Envgrid, "name_x")
  attr(Predictions, "name_y") <- attr(Envgrid, "name_y")

  return(Predictions)
}


## interpolateSurface --------------------------------
# Environmental <- tar_read(Diff_envko_env)
# envname <-  tar_read(predictor_select)
# parname <- tar_read(parname_environmental_diff_env)
# parname <-  NULL
interpolateSurface <- function(Environmental, parname = NULL, aggr = c("mean", "sd", "p10", "p25", "median", "p75", "p90"),
                               envname = tar_read(predictor_select), avg = c("env", "loc")) {
  
  aggr <- first(match.arg(aggr))
  
  if(!is.null(parname)) {
    parname <- selectParnameEnvironmental(parname, Environmental)
    Environmental %<>%  dplyr::filter(.variable %in% parname)
  } else {
    parname <- unique(Environmental$.variable)
  }
  
  # koname <- parname[str_detect(parname, "_ko_")]
  # koname_ba <- koname[str_detect(koname, "ba_fix_")]
  # koname_frac <- koname[str_detect(koname, "ba_frac_fix_")]
  
  E <- distinct(Environmental, loc, .keep_all = T) ## to make the interaction() faster, will be joined to Environmental
  if(match.arg(avg) == "env") grp <- interaction(E[envname], drop = T) else grp <- E$loc ## interaction works with lists!
  E <- bind_cols(E, grp = grp) %>%
    dplyr::select(grp, loc)
  
  Environmental %<>%
    left_join(E, by = "loc") %>%
    group_by_at(c("grp", ".variable", "tax", envname)) %>%
    summarize(mean = mean(.value, na.rm = T),
                     sd = sd(.value, na.rm = T),
                     p10 = quantile(.value, probs = 0.1, na.rm = T),
                     p25 = quantile(.value, probs = 0.25, na.rm = T),
                     median = median(.value, na.rm = T),
                     p75 = quantile(.value, probs = 0.75, na.rm = T),
                     p90 = quantile(.value, probs = 0.9, na.rm = T),
                     .groups = "drop") %>%
    dplyr::mutate(.value = !!ensym(aggr), .aggr = aggr) ## just rename based on arg aggr
  
  interpolateVar <- function(v, t) {
    E <- filter(Environmental, .variable == v & tax == t)
    
    ## use akima?
    I <- interp::interp(x = E[[envname[1]]], y = E[[envname[2]]], z = E$.value, nx = 500, ny = 500, method = "linear") ## duplicate = "mean" would also work
    I <- data.frame(expand.grid(I[1:2]), c(I$z)) %>%
      setNames(c(envname[1:2], 'z'))
    
    ## Check
    # ggplot(data = I, aes_string(x = envname[1], y = envname[2], z = "z")) + geom_contour_filled() + scale_y_reverse()
    # print(paste(t, v))

    
    attr(I, "name_x") <- envname[1]
    attr(I, "name_y") <- envname[2]
    attr(I, "taxon") <- t
    attr(I, "parname") <- v
    
    
    return(I)
  }
  
  Comb <- distinct(Environmental[c(".variable", "tax")], .variable, tax)
  surfaces <- mapply(interpolateVar, Comb$.variable, Comb$tax, SIMPLIFY = F) # mc.cores = ceiling(nrow(Comb)/3)
  
  return(surfaces)
}


## formatNumber --------------------------------
# x <- 34364343.24324
# signif.digits <- 4
formatNumber <- function(x, signif.digits = 4) {
  formatC(signif(x, digits = signif.digits), digits = signif.digits,format="fg", flag="#")
}


## formatEnvgrid --------------------------------
# envname <- tar_read("predictor_select")
# data_stan <- tar_read("data_stan_priors_offset_env")
formatEnvgrid <- function(data_stan, envname, res = 500) {
  
  if(length(envname) != 2) stop("There are more or fewer than 2 environmental gradients.")
  
  L <- attr(data_stan, "Long")
  
  envname_s <- paste0(envname, "_s")
  name_x <- envname_s[1]
  name_y <- envname_s[2]
  
  O <- cbind(x = L[[name_x]], y = L[[name_y]])
  range_x <- range(L[[name_x]], na.rm = T)
  range_y <- range(L[[name_y]], na.rm = T)
  
  ch <- chull(O)
  Hullpoint <- O[c(ch, ch[1]),] ## close the polygon
  poly <- st_sfc(st_polygon(list(Hullpoint)))
  
  x <- seq(range_x[1], range_x[2], length.out = res)
  y <- seq(range_y[1], range_y[2], length.out = res)
  P <- expand.grid(x, y) %>% setNames(c(name_x, name_y))
  points <- st_sfc(lapply(1:nrow(P), function(i) st_point(c(as.matrix(P[i,])))))
  
  iscovered <- st_covered_by(points, poly, sparse = F)
  P_covered <- P[iscovered,]
  
  attr(P_covered, "name_x") <- name_x
  attr(P_covered, "name_y") <- name_y
  
  return(P_covered)
}


# ————————————————————————————————————————————————————————————————————————————————— #
# Summarize posterior         -----------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## summarizeFit --------------------------------
# cmdstanfit <- tar_read("fit_env")
# publishpar <- tar_read(parname_plotorder)
# publishpar <- c(tar_read(parname_plotorder), tar_read(parname_env_vertex))
# exclude <- tar_read(exclude)
# path <- tar_read("dir_publish")
summarizeFit <- function(cmdstanfit, exclude = NULL, publishpar, path) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  allpar <- cmdstanfit$metadata()$stan_variables
  includepar <- setdiff(allpar, exclude)
  phipar <- includepar[str_starts(includepar, "phi")]
  sigmapar <- includepar[str_starts(includepar, "sigma_")]
  publishpar_prior <- c(publishpar, paste0(publishpar, "_prior"))
  publishpar <- intersect(c(sigmapar, publishpar_prior), allpar)
  
  summary <- cmdstanfit$summary(includepar)
  write.csv(summary, paste0(path, "/", basename_cmdstanfit, "_summary.csv"))
  
  summary_publish <- cmdstanfit$summary(publishpar) %>%
    mutate(p = if_else(str_detect(variable, "_prior"), "prior", "posterior")) %>%
    mutate(tax = if_else(str_detect(variable, "\\[2\\]"), "other", "Fagus")) %>%
    mutate(var = str_extract(variable, ".*_log(_center_env[1-2]|_spread_env[1-2])?")) %>% # ".*_log" and optionally one of the terms in braces
    mutate(value = paste0(formatNumber(mean), " ± ", formatNumber(sd))) %>%
    dplyr::select(var, p, tax, value, ess_bulk) %>%
    pivot_wider(values_from = c("value", "ess_bulk"), names_from = c("p", "tax"), id_cols = "var") %>%
    dplyr::select(-starts_with("ess_bulk_prior"))
  # write.csv(summary_publish, paste0(path, "/", basename_cmdstanfit, "_summary_parameters.csv"))
  
  summary_phipar <- cmdstanfit$summary(phipar) %>%
    mutate(value = paste0(formatNumber(mean), " ± ", formatNumber(sd))) %>%
    dplyr::select(var = variable, value, ess_bulk)
  # write.csv(summary_phipar, paste0(path, "/", basename_cmdstanfit, "_summary_phi.csv"))
  
  summary_publish  <-  dplyr::bind_rows(summary_publish, summary_phipar) %>%
    as.data.frame() %>%
    apply(2, as.character) ## Flattens list columns, that seem to emerge invisibly
  
  write.csv(summary_publish, paste0(path, "/", basename_cmdstanfit, "_summary_parameters.csv"))
  
  ## Number of years until equilibrium
  Iter <- cmdstanfit$draws("iterations_fix") %>%
    as_draws_matrix()
  Iter <- data.frame(min = min(Iter), median = median(Iter), max = max(Iter))
  write.csv(Iter, paste0(path, "/", basename_cmdstanfit, "_summary_nyears.csv"))
  
  
  ## Console output
  
  drawconverged <- cmdstanfit$draws(variables = "converged_fix", format = "draws_matrix") %>%
    apply(1, all)
  
  if (all(drawconverged)) {
    message("All clusters in all draws have converged to the fix point, i.e. the population trajectories ran into an equilibrium.")
  } else {
    message("For ", sum(!drawconverged), " draws, not all of the clusters have converged to the fix point, i.e. not all population trajectories ran into an equilibrium.")
  }
  
  head(summary, 20) %>%
    as.data.frame() %>%
    print()
  
  summary_publish %>%
    as.data.frame() %>%
    print()
  
  return(summary)
}


## summarizeMarginal --------------------------------
# Marginal <- tar_read("Marginal_env")
# basename <- tar_read("basename_fit_env")
# path <- tar_read("dir_publish")
summarizeMarginal <- function(Marginal, basename, path) {
  
  M <- Marginal %>%
    group_by(tax, .variable) %>%
    summarize(m = mean(.value), sd = sd(.value)) %>%
    mutate(value = paste0(formatNumber(m), " ± ", formatNumber(sd))) %>%
    pivot_wider(values_from = "value", names_from = "tax", id_cols = ".variable")
  
  write.csv(M, paste0(path, "/", basename, "_summary_marginal.csv"))
  
  return(M)
}


## summarizeStates --------------------------------
# States <- tar_read("States")[[1]]
# data_stan <- tar_read("data_stan")
# path <- tar_read("dir_publish")
summarizeStates <- function(States, data_stan, basename, path) {
  
  D <- attr(data_stan, "Long_BA") %>%
    st_drop_geometry() %>% ## avoiding potential problems from sf methods
    mutate(obsid = as.numeric(factor(obsid, levels = c("DE_BWI_1987", "DE_BWI_2002", "DE_BWI_2012")))) %>% # mutate(obsid = str_extract(obsid, "[0-9]{4}")) %>%
    group_by(stage, tax, obsid) %>%
    summarize(mean = mean(y_prior, na.rm = T), sd = sd(y_prior, na.rm = T)) %>% ## y is count_ha for J and A, ba_ha for B and BA
    mutate(value = paste0(formatNumber(mean, signif.digits = 5), " ± ", formatNumber(sd, signif.digits = 5))) %>%
    pivot_wider(names_from = "tax", id_cols = c("stage", "obsid")) %>%
    mutate(stage = paste0(stage, "_data_", obsid)) %>%
    dplyr::select(var = stage, Fagus = Fagus.sylvatica, other)
  
  S <- States %>%
    mutate(value = if_else(tax == 'other' & (str_starts(var, "major")),
                           1 - value,
                           value)
           ) %>%
    group_by(var, tax) %>%
    summarize(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T)) %>%
    mutate(value = paste0(formatNumber(mean, signif.digits = 5), " ± ", formatNumber(sd, signif.digits = 5))) %>%
    pivot_wider(names_from = "tax", id_cols = "var") %>%
    bind_rows(D) %>%
    bind_rows(c(var = "ba_a_avg", setNames(formatNumber(data_stan$ba_a_avg, signif.digits = 5), c("Fagus", "other")))) %>%
    arrange(var)
  
  write.csv(S, paste0(path, "/", basename, "_summary_states.csv"))
  print(S)
  
  return(S)
}


## summarizeFreqConverged --------------------------------
# cmdstanfit <- tar_read("fit_test")
# data_stan_priors <- tar_read("data_stan_priors")
# path <- tar_read("dir_publish")
summarizeFreqConverged <- function(cmdstanfit, data_stan_priors, path) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  n_locs <- data_stan_priors$N_locs
  Freq_converged <- formatLoc("converged_fix", locmeans = T, cmdstanfit_ = cmdstanfit, data_stan_priors_ = data_stan_priors)
  
  write.csv(Freq_converged, file.path(path, paste0(basename_cmdstanfit, "Freq_fixpoint_converged.csv")))
  
  message(sum(Freq_converged$value != 1), " of ", n_locs, " locations have not always converged to fixpoint.")
  return(Freq_converged)
}



# ————————————————————————————————————————————————————————————————————————————————— #
# Generate posterior quantities         --------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## generateResiduals --------------------------------
# cmdstanfit  <- tar_read("fit_test")
# data_stan_priors  <- tar_read("data_stan_priors")
# path  <- tar_read("dir_fit")
# yhatvar <- "y_hat_offset"
# includeinit <- FALSE
generateResiduals <- function(cmdstanfit, data_stan_priors, yhatvar = c("y_hat_offset", "y_hat_rep_offset"), includeinit = FALSE, path) {
  
  yhatvar <- match.arg(yhatvar)
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  Sim <- cmdstanfit$draws(variables = "y_sim", format = "draws_matrix") %>% t()# matrix of observations simulated from the fitted model - row index for observations and column index for simulations
  Sim[is.na(Sim)] <- 0
  y <- data_stan_priors$y
  y_hat <- cmdstanfit$draws(variables = yhatvar, format = "draws_matrix") %>% apply(2, median, na.rm = T)
  
  Longdata <- attr(data_stan_priors, "Long")
  grp <- with(Longdata, interaction(as.integer(as.factor(obsid)), stage, substr(tax, 1, 1)))
  offset <- data_stan_priors$offset_data
  isnotinit <- !Longdata$isy0
  
  if(includeinit) {
    residuals <- DHARMa::createDHARMa(simulatedResponse = Sim, observedResponse = y, fittedPredictedResponse = y_hat, integerResponse = T)
  } else {
    residuals <- DHARMa::createDHARMa(simulatedResponse = Sim[isnotinit,], observedResponse = y[isnotinit], fittedPredictedResponse = y_hat[isnotinit], integerResponse = T)
    grp <- droplevels(grp[isnotinit])
  }
  
  png(paste0(path, "/", basename_cmdstanfit, "_", "DHARMa", ".png"), width = 1600, height = 1000)
  plot(residuals, quantreg = T, smoothScatter = F)
  dev.off()
  
  # residuals_recalc <- recalculateResiduals(residuals, group = grp)
  # plot(residuals_recalc, quantreg = T, smoothScatter = F)
  
  png(paste0(path, "/", basename_cmdstanfit, "_", "DHARMa_grouped", ".png"), width = 2200, height = 800)
  plot(residuals, form = grp, quantreg = T, smoothScatter = F)
  dev.off()
  
  # png(paste0(path, "/", basename_cmdstanfit, "_", "DHARMa_offset", ".png"), width = 2200, height = 800)
  # plot(residuals, form = offset, quantreg = T, smoothScatter = F)
  # dev.off()
  
  return(residuals)
}


## generatePredictiveChecks --------------------------------
## This is more of an opportunity for interactive checks
# cmdstanfit  <- tar_read("fit_test")
# data_stan_priors_offset  <- tar_read("data_stan_priors_offset")
generatePredictiveChecks <- function(cmdstanfit, data_stan_priors_offset, path) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  Draws <- cmdstanfit$draws(variables = c("y_sim", "y_hat", "y_hat_offset")) %>%
    as_draws_rvars()
  
  d <- data_stan_priors_offset
  L <- attr(data_stan_priors_offset, "Long") %>%
    
    bind_cols(y_sim = E(Draws$y_sim),
              y_hat = E(Draws$y_hat),
              y_hat_offset = E(Draws$y_hat_offset)) # %>%
  
  #   ## depending on what was generated in the model
  #   mutate(y_hat_rep_offset = y_hat_offset,
  #          y_hat_rep = y_hat)
  
  # grp <- with(Longdata, interaction(as.integer(as.factor(obsid)), stage, substr(tax, 1, 1)))
  
  # plot(y_sim ~ y_hat_rep_offset, data = L)
  #   ## -> does simulation work?
  
  # plot(y_hat_rep_offset ~ I(y_hat_rep * offset), data = L)
  #   ## -> does offset work?
  
  # plot(sqrt(L$y_hat_prior), sqrt(d$y_hat_prior[d$rep_yhat2y]))
  # plot(sqrt(y_hat_rep) ~ sqrt(y_hat_prior), data = L) ## the former is the extracted from the model, the latter is the correct data
  #   ## -> does rep_yhat2y and indexing work?
  
  # plot(sqrt(y_sim) ~ sqrt(y_hat_prior * offset), data = L)
  # plot(sqrt(y_hat_prior_rep) ~ sqrt(y_hat_prior * offset), data = L)
  #   ## -> the model works fine
  
  # plot(sqrt(y_sim) ~ sqrt(y), data = L)
  # plot(log1p(y_hat_rep_offset) ~ log1p(y), data = L)
  #   ## -> Averaging could be problematic. These should be exactly the same.
  
  # lattice::xyplot(log1p(y_sim) ~ log1p(y_hat_prior) | t_which, data = L)
  # lattice::xyplot(log1p(y_hat_rep) ~ log1p(y_prior) | pop + t_which, data = L)
  
  # plot(y_hat_rep_offset ~ y_hat_rep, data = L)
  
  return(L)
}


## generateTrajectories --------------------------------
# cmdstanfit <- tar_read("fit_env")
# parname <- tar_read("parname")
# data_stan_priors <- tar_read("data_stan_priors_env")
# locparname <- tar_read("parname_trajectories_env")

generateTrajectories <- function(cmdstanfit, data_stan_priors, parname, locparname = c("state_init", "L_loc"),
                                 time = c(1:25, seq(30, 300, by = 10), seq(400, 5000, by = 100)), thinstep = 1,
                                 average = c("none", "locsperdraws_all", "drawsperlocs_all", "locsperdraws_avgL", "locsperdraws_avgL_qInit")) {
  
  varname_draws <- cmdstanfit$metadata()$stan_variables
  
  ## Global parameters
  parname <- intersect(parname, varname_draws)
  parname_sans_log <- gsub("_log$", "", parname)
  
  ## Local variables
  locparname <- intersect(locparname, varname_draws)
  locparname_avg <- paste0("avg_", locparname) %>%
    intersect(varname_draws)
  # locparname_unavg <- setdiff(locparname, paste0(gsub("avg_", "", locparname_avg))) ## pars that were never averaged. currently not used
  locparname_all <- c(locparname_avg, locparname)
  
  
  Draws <- cmdstanfit$draws(variables = c(parname, locparname_all), format = "draws_list") %>%
    thin_draws(thin = thinstep)
  
  
  #### The goal of the following wrangling:
  ## to makes pars a nested list[draws[parameters[n_locs,n_species]] (list might be a data frame),
  ## The list will be structured here in two steps:
  ##     1. for global variables "p_log", which are exped to "p" here,
  ##     2. for local variables "P_log"/"P_loc" (which are exped later, see below).
  ## Later, when the simulation is applied over the nested parameter list, the structure will be unpacked, and parameters selected (local or global parameter?).
  ## In particular, the averaging at the loc-level over draws ("drawsperloc") will be done in the nested process.
  ## Note that exping (and renaming) must always happen after averaging! This is why the local variables get exped later.
  
  ## 1. Distinction between when to average for global variables
  if(match.arg(average) == "drawsperlocs_all") {
    
    ## averaging the draws of non-local parameters ("within 1 loc")!
    
    M <- cmdstanfit$summary(variables = parname) %>%
      dplyr::select(variable, mean)
    parmeans <- sapply(parname, function(n) M$mean[str_starts(M$variable, n)], USE.NAMES = T, simplify = F)
    pars <- sapply(parname, function(n) { if(str_ends(n, "_log")) exp(parmeans[[n]]) else parmeans[[n]] }, USE.NAMES = T, simplify = F)
    names(pars) <- parname_sans_log
    pars <- list(pars) ## nested list[1 draw[parameter[vector_2species]]]
    
  } else if (match.arg(average) %in% c("none", "locsperdraws_all", "locsperdraws_avgL", "locsperdraws_avgL_qInit")) {
    
    ## no averaging of the non-local parameters!
    
    D <- subset_draws(Draws, variable = c(parname)) %>%
      as_draws_matrix() %>%
      as.data.frame() %>%
      dplyr::mutate(across(contains("_log["), exp))
    pars <- split(D, 1:nrow(D))
    pars <- lapply(pars, function(p) matrix(c(as.matrix(p)), byrow = F, nrow = data_stan_priors$N_species, dimnames = list(NULL, parname_sans_log)))
    pars <- lapply(pars, as.data.frame) ## nested list[draws[parameter[vector_2species]]]
    
  }
  
  
  ## 2. Distinction between when to average for local variables
  if (match.arg(average) %in% c("none", "drawsperlocs_all", "locsperdraws_avgL", "locsperdraws_avgL_qInit")) {
    
    draws_loc <- subset_draws(Draws, variable = locparname) %>% # locparname important, because will get iterated over locs
      as_draws_rvars()
    
    n_locs <- data_stan_priors$N_locs
    
  }
  
  if (match.arg(average) %in% c("locsperdraws_all", "locsperdraws_avgL", "locsperdraws_avgL_qInit")) {
    
    warning("Averaging over locs has not been implemented for variables other than ('*_avg'), i.e. variables that were already averaged in the model.")
    
    draws_loc_avg <- subset_draws(Draws, variable = locparname_avg) %>%
      as_draws_rvars() %>%
      lapply(as.matrix) %>%
      lapply(t) ## make this a list of 1-rowed matrices
    names(draws_loc_avg) <- gsub("^avg_", "", locparname_avg)
    
    n_locs_avg <- 1
    
  }
  
  ## Generate quantiles in any case!
  draws_loc_q <- subset_draws(Draws, variable = locparname) %>% # c("state_init", "L_loc")
    posterior::as_draws()
  
  Quantiles_init <- draws_loc_q %>%
    tidybayes::gather_draws(state_init[loc,pop]) %>%
    group_by(pop, .draw, .iteration, .chain) %>%
    summarize(pop = first(pop),
              p10 = quantile(.value, prob = 0.1, type = 1, na.rm = T),
              median = quantile(.value, prob = 0.5, type = 1, na.rm = T),
              p90 = quantile(.value, prob = 0.9, type = 1, na.rm = T),
              loc_p10 = first(loc[.value == p10]),
              loc_median = first(loc[.value == median]),
              loc_p90 = first(loc[.value == p90]),
              loc_p10_b = replace(loc_p10, pop < 5, NA),
              loc_median_b = replace(loc_median, pop < 5, NA),
              loc_p90_b = replace(loc_p90, pop < 5, NA)) %>%
    ungroup()
  
  if (match.arg(average) == "locsperdraws_avgL_qInit") {
    
    ## Here 3 quantiles will be dealt with, as if they were locs
    
    L_loc_q <- draws_loc_q %>%
      tidybayes::gather_draws(L_loc[loc,tax]) %>% # state_init[loc,pop]
      group_by(tax, .draw, .iteration, .chain) %>%
      summarize(p10 = quantile(.value, prob = 0.1, type = 1, na.rm = T),
                median = quantile(.value, prob = 0.5, type = 1, na.rm = T),
                p90 = quantile(.value, prob = 0.9, type = 1, na.rm = T)) %>%
      pivot_longer(any_of(c("p10", "median", "p90")), names_to = "quantile", values_to = "L_loc") %>%
      mutate(loc = as.integer(factor(quantile, levels = c("p10", "median", "p90")))) %>%
      pivot_wider(id_cols = c(".draw", ".iteration", ".chain"), names_from = c("loc", "tax"), values_from = "L_loc", names_glue = "L_loc[{loc},{tax}]") %>% 
      as_draws_rvars()
    
    state_init_q <- Quantiles_init %>%
      pivot_longer(any_of(c("p10", "median", "p90")), names_to = "quantile", values_to = "state_init") %>%
      mutate(loc = as.integer(factor(quantile, levels = c("p10", "median", "p90")))) %>%
      # mutate(state_init = exp(state_init_log)) %>% ## !!!
      pivot_wider(id_cols = c(".draw", ".iteration", ".chain"), names_from = c("loc", "pop"), values_from = "state_init", names_glue = "state_init[{loc},{pop}]") %>% 
      as_draws_rvars()
    
    draws_loc_q <- list(state_init = state_init_q$state_init, L_loc = L_loc_q$L_loc)
    n_locs_q <- 3
  }
  
  if (match.arg(average) == "locsperdraws_all") {
    draws_loc <- draws_loc_avg
    n_locs <- n_locs_avg
  } else if (match.arg(average) == "locsperdraws_avgL") {
    L_loc_avg_rep <- draws_loc_avg$L_loc[rep(1, n_locs),]
    draws_loc$L_loc <- L_loc_avg_rep
  } else if (match.arg(average) == "locsperdraws_avgL_qInit") {
    L_loc_avg_rep <- draws_loc_avg$L_loc[rep(1, n_locs_q),]
    draws_loc <- draws_loc_q
    n_locs <- n_locs_q
  }
  
  
  ######-----------------------------------------------------
  ### Nested simulations: locs/draws
  ######-----------------------------------------------------
  
  ### Functions ---------------------------------------
  ## (ordered inner to outer)
  
  ##### iterateModel -------------------------
  #
  iterateModel <- function(initialstate,
                           pars,
                           times,
                           format = c("long", "matrix")) {
    
    b <- pars$b # Length n_species vector of basal area increment rates.
    c_a <- pars$c_a # Length n_species vector ## only used in model versions, where there are distinct parameters for competition on b an a
    c_b <- pars$c_b # Length n_species vector
    c_j <- pars$c_j # Length n_species vector
    g <- pars$g # Length n_species vector of transition rates.
    h <- pars$h # Length n_species vector of transition rates.
    
    l <- pars$l # Length n_species vector of input rates.
    # l <- pars$k # Length n_species vector of input rates.
    
    r <- pars$r # Length n_species vector of input rates
    s <- pars$s # Length n_species vector of shading rates
    
    ba_a_upper <- data_stan_priors$ba_a_upper
    ba_a_avg <- data_stan_priors$ba_a_avg
    
    ## Set the count variables
    n <- length(r) # no. of species
    
    times_intern <- 1:max(times)
    n_times <- length(times_intern)
    
    # Prepare a state matrix
    whichstate <- rep(1:3, each = n)
    State <- matrix(rep(initialstate, times = n_times), nrow = n_times, byrow = T)
    
    ## Here comes the model.
    for (t in 2:n_times) {
      
      ## States at t-1: *State*, and *State*
      J <- State[t-1,whichstate == 1]
      A <- State[t-1,whichstate == 2]
      B <- State[t-1,whichstate == 3]
      
      ## The total basal area of big trees
      BA <- A * ba_a_avg + B
      BA_sum <- sum(BA)
      
      J_trans <- r*BA + l + (J - g*J) # use rlnorm + u[ ,1]
      J_t <- J_trans * 1/(1 + c_j*sum(J) + s*BA_sum) # count of juveniles J
      # with all variables > 0, and 0 < g, m_j < 1
      
      A_trans <- g*J + (A - h*A) # + u[ ,2]
      A_t <-  A_trans * 1/(1 + c_a*BA_sum) # count of small adults
      # with all variables > 0, and 0 < h, m_a < 1
      
      A_ba <- A * h * ba_a_upper # Basal area of small adults A. Conversion by multiplication with basal area of State exit (based on upper dhh boundary of the class)
      B_trans <- A_ba + B # + u[ ,3]
      B_t <- (1+b)*B_trans * 1/(1 + c_b*BA_sum)  # basal area of big adults B
      ## b is the net basal area increment (including density-independent m) basically equivalent to a Ricker model, i.e. constant increment rate leading to exponential growth, negative density dependent limitation scaled with total BA_sum.
      
      State[t, ] <- c(J_t, A_t, B_t)
    }
    
    whichtimes <- which(times_intern %in% times)
    State <- State[whichtimes,]
    
    if(match.arg(format) == "long") State <- data.frame(abundance = c(State),
                                                        stage = factor(rep(c("J", "A", "B"), each = n * length(times)), levels = c("J", "A", "B")),
                                                        tax = rep(1:n, each = length(times)),
                                                        time = times)
    
    return(State)
  }
  
  
  ##### iterateModel_draws() ------------------
  iterateModel_draws <- function(locpars, pars, time, averageperlocs) {
    
    if (averageperlocs) {
      locpars <- lapply(locpars, mean, na.omit = T) ## will work with rvars, while retaining the data structure of a matrix, e.g.
      
    } else {
      locpars <- lapply(locpars, as_draws_matrix)
    }
    
    ## Exping in any case, but after averaging!
    name_locpars <- names(locpars)
    locpars <- lapply(name_locpars, function(n) { if(str_ends(n, "_log")) exp(locpars[[n]]) else locpars[[n]] })
    names(locpars) <- str_remove(name_locpars, "_log")
    
    ## Assign local parameters to draws, adopt manually if necessary
    pars <- lapply(1:length(pars), function(i) within(pars[[i]], { if("L_loc" %in% locparname) l <- c(locpars$L_loc[i,])
                                                                   if("B_log" %in% locparname) b <- c(locpars$B[i,])
                                                                   if("C_a_log" %in% locparname) c_a <- c(locpars$C_a[i,])
                                                                   if("C_b_log" %in% locparname) c_b <- c(locpars$C_b[i,])
                                                                   if("C_j_log" %in% locparname) c_j <- c(locpars$C_j[i,])
                                                                   if("G_log" %in% locparname) g <- c(locpars$G[i,])
                                                                   if("H_log" %in% locparname) h <- c(locpars$H[i,])
                                                                   if("R_log" %in% locparname) r <- c(locpars$R[i,])
                                                                   if("S_log" %in% locparname) s <- c(locpars$S[i,])
                                                                 }
                                                      ))
    
    return( lapply(1:length(pars), function(i) iterateModel(c(locpars$state_init[i,]), pars = pars[[i]], time)) )
    
  }
  
  ##### iterateLocs -------------------
  iterateLocs <- function(i, lp, p, t, avgperlocs) {
    iterateModel_draws(locpars = lapply(lp, function(x) x[i,]), pars = p, time = t, averageperlocs = avgperlocs)
  }
  
  
  ### Apply functions over par structure -------------------
  
  averageperlocs <- match.arg(average) == "drawsperlocs_all"
  
  sims <- future_sapply(1:n_locs, iterateLocs,
                        lp = draws_loc, p = pars, t = time, avg = averageperlocs, simplify = F) ## returns a nested list[locs, draws] of matrices[times]
  sims <- lapply(sims, bind_rows, .id = "draw")
  Sims <- bind_rows(sims, .id = "loc")
  
  Quantiles_init <- Quantiles_init %>% 
    mutate(pop = as.character(pop),
           stage = fct_recode(pop, "J" = "1", "J" = "2", "A" = "3", "A" = "4", "B" = "5", "B" = "6"),
           tax = as.integer(fct_recode(pop, "1" = "1", "2" = "2", "1" = "3", "2" = "4", "1" = "5", "2" = "6"))) %>%
    rename(draw = ".draw") %>%
    dplyr::select(-c(".chain", "pop", ".iteration"))
  
  Sims <- Sims %>%
    group_by(loc, draw) %>%
    mutate(diff = abundance - c(NA, abundance[1:(n()-1)]),
           absdiff = abs(diff),
           isconverged = replace_na(absdiff <= data_stan_priors$tolerance_fix, F),
           isflat = replace_na(absdiff <= data_stan_priors$tolerance_fix * 10, F),
           time_fix = first(time[isconverged]),
           time_flat = first(time[isflat]),
           state_fix = first(abundance[isconverged]),
           state_flat = first(abundance[isflat])) %>%
    ungroup() %>%
    mutate(draw = as.integer(draw)) %>%
    left_join(Quantiles_init, by = c("draw", "stage", "tax"))
  
  return(Sims)
}


## selectParnameEnvironmental --------------------------------
## Helper for excluding variables that are not in the fit.
selectParnameEnvironmental <- function(parname, Environmental) {
  
  parname_drawn <- unique(Environmental$.variable)
  notdrawn <- paste(setdiff(parname, parname_drawn), collapse = ", ")
  if (length(notdrawn) > 1) message("selectParnameEnvironmental(): The variables ", notdrawn, " are not in the posterior to be regressed against the environmental variables!")
  parname_environmental <- intersect(parname, parname_drawn)
  
  ## reorder and concatenate objects that are unknown to the order
  o <- c("L_loc", "L_log", "R_log", "C_j_log", "S_log",
         "G_log", "C_a_log",
         "H_log", "B_log", "C_b_log",
         "ba_init", "ba_fix", "ba_frac_init", "ba_frac_fix", "major_init", "major_fix")
  i <- match(o, parname)
  parname_o <- parname_environmental[i]
  parname_environmental <- union(parname_o, parname_environmental)
  names(parname_environmental) <- str_remove(parname_environmental, "_log")
  
  return(parname_environmental[!is.na(parname_environmental)])
}


## fitEnvironmental_gam --------------------------------
# Environmental <- tar_read("Environmental_env")
# Environmental <- tar_read("Diff_envko_env")
# parname <- tar_read("parname_environmental_ba_env")[4]
# envname <- tar_read("predictor_select")
# path  <- tar_read("dir_publish")
fitEnvironmental_gam <- function(Environmental, parname, envname = tar_read(predictor_select),
                                 taxon = c(1:2, 0), fam = c("gaussian", "binomial"), path = tar_read("dir_publish")) {
  
  parname <- selectParnameEnvironmental(parname, Environmental)

  taxon <- head(as.integer(taxon), 1)
  fam <- match.arg(fam)
  
  E <- Environmental %>% 
    filter(tax == taxon) %>%
    filter(.variable == parname) %>%
    rename(v = .value)
  
  splineformula <- paste0("v ~ ", "te(", paste(envname, collapse = ", "), ", k = c(5, 5))")
  
  ### mgcv
  fit <- gam(as.formula(splineformula), family = fam, data = E)
  
  if(!is.null(path)) {
    s <- summary(fit)
    textext <- itsadug::gamtabs(s, caption = "Summary of the thin plate spline fit for environmental ...", label = paste("tab:gam", parname, taxon, sep = "_"))
    dir.create(file.path(path, "fits_Environmental"), showWarnings = FALSE)
    cat(textext, file = file.path(path, "fits_Environmental", paste0("Environmental", "_", parname, "_", taxon, "_summary_gam.tex")), fill = T) %>% invisible()
  }
  
  attr(fit, "tax") <- taxon
  attr(fit, "par") <- parname
  
  return(fit)
}


## fitEnvironmental_glmnet --------------------------------
# Environmental <- tar_read("Environmental_env")
# parname <- tar_read("parname_environmental")[3]
# parname <- tar_read("parname_environmental_binomial")[1]
# envname <- tar_read("predictor_select")
# path  <- tar_read("dir_fit")
fitEnvironmental_glmnet <- function(Environmental, parname, envname = tar_read(predictor_select),
                                    taxon = c(1:2, 0), fam = c("gaussian", "binomial"), path = tar_read("dir_publish")) {
  
  parname <- selectParnameEnvironmental(parname, Environmental)
  
  taxon <- head(as.integer(taxon), 1)
  fam <- match.arg(fam)
  
  E <- Environmental %>% 
    filter(tax %in% taxon) %>%
    filter(.variable == parname) %>%
    rename(v = .value)
  
  if (fam == "binomial") {
    E$v <- as.factor(round(E$v))
  }
  
  formula <- paste0("~ 1 + ", paste0("poly(", envname, ", 2, raw = T)" , collapse = " + "))
  M <- model.matrix(as.formula(formula), data = E)
  
  ### glmnet
  fit <- cv.glmnet(x = M, y = E$v, family = fam, nlambda = 100, nfolds = 20, data = E,  parallel = TRUE) # plot(fit)
  
  
  if(!is.null(path)) {
    
    S <- coef(fit, s = fit$lambda.1se) %>% as.matrix()
    dir.create(file.path(path, "fits_Environmental"), showWarnings = FALSE)
    write.csv(S, file = file.path(path, "fits_Environmental", paste0("Environmental", "_", parname, "_", taxon, "_summary_glmnet.csv")))
    
  }
  
  attr(fit, "E") <- E
  # attr(fit, "M") <- M
  attr(fit, "tax") <- taxon
  attr(fit, "par") <- parname
  
  return(fit)
}


## fitEnvironmental_glm --------------------------------
# Environmental <- tar_read("Environmental_env")
# parname <- tar_read("parname_environmental")[3]
# parname <- tar_read("parname_environmental_binomial")[1]
# envname <- tar_read("predictor_select")
# path  <- tar_read("dir_fit")
fitEnvironmental_glm <- function(Environmental, parname, envname = tar_read(predictor_select),
                                 taxon = c(1:2, 0), fam = c("gaussian", "binomial"), path = tar_read("dir_publish")) {
  
  parname <- selectParnameEnvironmental(parname, Environmental)
  
  taxon <- head(as.integer(taxon), 1)
  fam <- match.arg(fam)
  
  E <- Environmental %>% 
    filter(tax %in% taxon) %>%
    filter(.variable == parname) %>%
    rename(v = .value)
  
  
  formula <- paste0("v ~ 1 + ", paste0("poly(", envname, ", 2, raw = T)" , collapse = " + "))
  
  ### glm
  ## glm, compared to unregularized glmnet, has the advantage to predict continuosly on [0, 1] when binomial
  fit <- glm(as.formula(formula), family = fam, data = E) # plot(fit)
  
  if(!is.null(path)) {
    dir.create(file.path(path, "fits_Environmental"), showWarnings = FALSE)
    stargazer(fit, out = file.path(path, "fits_Environmental", paste0("Environmental", "_", parname, "_", taxon, "_summary_glm.tex"))) %>% invisible()
  }
  
  attr(fit, "tax") <- taxon
  attr(fit, "par") <- parname
  
  return(fit)
}


## predictEnvironmental --------------------------------
# fit <- tar_read("fit_environmental_env")[[1]]
# envname <- tar_read("predictor_select")
# path  <- dir_fit <- tar_read("dir_publish")
# basename  <- basename_fit_env <-  tar_read("basename_fit_env")
# color  <- tar_read("twocolors")
# themefun  <- theme_fagus <-  tar_read("themefunction")
predictEnvironmental <- function(fit, envname,
                                 path = dir_fit, basename = basename_fit_env,
                                 color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  if(length(envname) != 2) stop("There are more or fewer than 2 environmental gradients.")
  
  taxon <- attr(fit, "tax")
  parname <- attr(fit, "par")
  
  res <- 400
  name_x <- envname[1]
  name_y <- envname[2]
  
  if (is(fit, "cv.glmnet")) {
    
    E <- attr(fit, "E")
    O <- cbind(x = E[[name_x]], y = E[[name_y]])
    range_x <- range(E[[name_x]], na.rm = T)
    range_y <- range(E[[name_y]], na.rm = T)
    
  } else if (is(fit, "gam")) { ## gams are also glms!
    
    O <- cbind(x = fit$model[[name_x]], y = fit$model[[name_y]])
    range_x <- range(fit$model[[name_x]], na.rm = T)
    range_y <- range(fit$model[[name_y]], na.rm = T)
    
  } else if (is(fit, "glm")) {
    
    O <- cbind(x = fit$data[[name_x]], y = fit$data[[name_y]])
    range_x <- range(fit$data[[name_x]], na.rm = T)
    range_y <- range(fit$data[[name_y]], na.rm = T)
    
  } else {
    
    warning("predictEnvironmental(): unknown class of fit for ", paste(parname, taxon), "!")
    return(NULL)
    
  }
  
  ch <- chull(O)
  Hullpoint <- O[c(ch, ch[1]),] ## close the polygon
  poly <- st_sfc(st_polygon(list(Hullpoint)))
  
  x <- seq(range_x[1], range_x[2], length.out = res)
  y <- seq(range_y[1], range_y[2], length.out = res)
  P <- expand.grid(x, y) %>% setNames(c(name_x, name_y))
  points <- st_sfc(lapply(1:nrow(P), function(i) st_point(c(as.matrix(P[i,])))))
  
  iscovered <- st_covered_by(points, poly, sparse = F)
  P_covered <- P[iscovered,]
  
  if (is(fit, "cv.glmnet")) {
    
    formula <- paste0("~ 1 + ", paste0("poly(", envname, ", 2, raw = T)" , collapse = " + "))
    X <- model.matrix(as.formula(formula), data = P_covered)
    p <- predict(fit, newx = X, type = "response", s = fit$lambda.1se) %>% c() # minimal lambda with the best prediction: lambda.min
    
  } else {
    p <- predict(fit, newdata = P_covered, type = "response") %>% c() # %>% matrix(nrow = res, ncol = res)
  }
  
  D <- cbind(P_covered, z = p)
  
  
  attr(D, "name_x") <- name_x
  attr(D, "name_y") <- name_y
  attr(D, "taxon") <- taxon
  attr(D, "parname") <- parname
  
  return(D)
}


# ————————————————————————————————————————————————————————————————————————————————— #
# Plot posterior         -----------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## plotStanfit --------------------------------
# stanfit  <- tar_read("stanfit")
# stanfit  <- tar_read("stanfit_test")
# stanfit  <- tar_read("stanfit_test_plotting")
# exclude <- tar_read("exclude")
# path  <- tar_read("dir_fit")
# basename  <- tar_read("basename_fit_test")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotStanfit <- function(stanfit, exclude, path, basename,
                        color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  usedmcmc <- "sample" == attr(stanfit, "stan_args")[[1]]$method
  
  # basename <- attr(stanfit, "model_name") %>%
  #   str_replace("-[1-9]-", "-x-")
  
  parname <- setdiff(stanfit@model_pars, exclude)
  parname_sansprior <- parname[!grepl("prior$", parname)]
  
  traceplot <- rstan::traceplot(stanfit, pars = parname_sansprior, include = T)
  # parallelplot_c <- bayesplot::mcmc_parcoord(stanfit, pars = vars(starts_with(c("c_", "s_"))))
  # parallelplot_others <- bayesplot::mcmc_parcoord(stanfit, pars = vars(!matches(c(exclude, "c_", "log_", "phi_", "lp_", "s_", "_prior"))))
  
  plots <- list(traceplot = traceplot) # parallelplot_c = parallelplot_c, parallelplot_others = parallelplot_others,
  
  mapply(function(p, n) ggsave(paste0(path, "/", basename, "_", n, ".png"), p, device = "png"), plots, names(plots))
  
  if(usedmcmc) {
    
    png(paste0(path, "/", basename, "_", "pairsplot", ".png"), width = 2600, height = 2600)
    pairs(stanfit, pars = c(parname_sansprior, "lp__"), include = T)
    dev.off()
    
  }
  
  return(NULL) # return(plots)
}


## plotParameters --------------------------------
# draws  <- tar_read("draws_env")
# exclude <- tar_read("exclude")
# parname <- tar_read("parname_plotorder")
# path  <- tar_read("dir_fit")
# basename  <- tar_read("basename_fit_env")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotParameters <- function(draws, parname, exclude = tar_read("exclude"), path, basename,
                           color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  extendedcolor <- c(color, "#555555") # add a third neutral colour for inspecific priors
  priorlinecolor <- c("#000000", "#000000")
  prioralpha <- c(1, 0.3)
  
  
  getRidgedata <- function(name, d = draws) {
    R <- bayesplot::mcmc_areas_data(d,
                                    pars = vars(any_of(paste0(name, c("[1]", "[2]", "_prior", "_prior[1]", "_prior[2]")))), # vars(starts_with(name, ignore.case = F)))
                                    point_est = c("median"),
                                    prob = 0.8)
    
    M <- bayesplot::mcmc_intervals_data(d,
                                        pars = vars(any_of(paste0(name, c("[1]", "[2]", "_prior", "_prior[1]", "_prior[2]"))))) # vars(starts_with(name, ignore.case = F)))
    
    R %<>%
      mutate(tax = str_extract(parameter, '\\[[12]\\]$')) %>%
      mutate(par = str_extract(parameter, '([a-z_])*')) %>%
      mutate(prior = str_ends(parameter, 'prior(\\[[12]\\])*')) %>%
      mutate(group = case_when(tax == '[1]' & prior ~ 'Fagus prior',
                               tax == '[2]' & prior ~ 'other prior',
                               tax == '[1]' & !prior ~ 'Fagus',
                               tax == '[2]' & !prior ~ 'other',
                               is.na(tax) & prior ~ 'prior')) %>%
      mutate(group = factor(group, levels = c('Fagus', 'other', 'Fagus prior', 'other prior', 'prior'), ordered = T)) %>%
      mutate(arrangement = sort(as.integer(group), decreasing = T)) %>%
      left_join(M[, c("parameter", "m")], by = "parameter")
    
    # group_by(parameter, tax) %>%
    # mutate(m_d = median(x, na.rm = T), max_d = max(scaled_density, na.rm = T)) %>%
    # ungroup()
    
    return(R)
  }
  
  
  plotRidges <- function(Data, plotlegend = FALSE) {
    
    Data$parameter <- fct_reorder(Data$parameter, as.integer(Data$group), .desc = T)
    
    ggplot(Data, aes(y = parameter, height = scaled_density, x = x, col = prior, fill = tax, alpha = prior)) +
      geom_density_ridges(stat = "identity", size = 0.8, rel_min_height = 0.01) +
      geom_segment(aes(x = m, xend = m, y = parameter, yend = as.integer(parameter) + scaled_density), color = "black", linetype = 3, size = 0.3) +
      scale_color_manual(values = priorlinecolor) +
      scale_fill_manual(values = extendedcolor) +
      scale_alpha_manual(values = prioralpha) +
      
      ggtitle(paste("log", str_remove(first(Data$par), "_log"))) +
      scale_y_discrete(labels = function(parameter) Data$group[match(parameter, Data$parameter)], expand = expansion(mult = c(0.1, 1))) +
      
      themefun() +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank()) +
      { if(!plotlegend) theme(legend.position="none") }
    
  }
  
  # plotRidges <- function(startswith, d = draws) {
  #   bayesplot::mcmc_areas_ridges(d, pars = vars(starts_with(startswith, ignore.case = F))) +
  #     themefun()
  # }
  
  parname_sansprior <- parname # parname[!grepl("prior$", parname)]
  ## parnamestart <- na.omit(unique(str_extract(parname_sansprior, "^[a-z]_[ljab]"))) ## Not used. Instead, full parname is used as matching character since it excludes vertex form parameters. Originally intended for later use in starts_with: Everything that starts with a small letter, and has the right index after that to be a meaningful parameter. (Small letter is important!)
  ## parname <- c(parname, paste0(parname, "_prior"))
  
  bayesplot::color_scheme_set("gray")
  ## areas for all parameters (the exclude might not be complete)
  # areasplot <- bayesplot::mcmc_areas(draws, area_method = "scaled height", pars = vars(!matches(c(exclude, "log_", "lp_", "prior")))) + themefun()
  
  ridgedata <- parallel::mclapply(parname_sansprior, getRidgedata, mc.cores = getOption("mc.cores", 18L))
  ridgeplots <- lapply(ridgedata, plotRidges)
  ridgeplotgrid <- cowplot::plot_grid(plotlist = ridgeplots,  align = "v")
  legendplot <- plotRidges(ridgedata[[1]], plotlegend = TRUE)
  
  plots <- list(ridgeplotgrid = ridgeplotgrid, ridge_legendplot = legendplot) # areasplot = areasplot
  
  mapply(function(p, n) ggsave(paste0(path, "/", basename, "_", n, ".pdf"), p, device = "pdf", height = 10, width = 12), plots, names(plots))
  message("Parameter plots complete.")
  
  return(NULL) # return(plots)
}


## plotPosterior --------------------------------
# varname  <- tar_read("parname_env_vertex_spread")
# cmdstanfit  <- tar_read("fit_env")
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotPosterior <- function(cmdstanfit, varname,
                          path = tar_read("dir_publish"), color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  varname_draws <- cmdstanfit$metadata()$stan_variables
  varname <- intersect(varname, varname_draws)
  
  if(length(varname) == 0) {
    warning("plotPosterior(): None of the provided names are variables in the fit.")
    return(NULL)
  }
  
  d <- cmdstanfit$draws(variables = varname)
  n_vars <- length(unique(attr(d, "dimnames")$variable))
  
  rhat <- rhat(cmdstanfit, pars = varname)
  # nutsparam <- nuts_params(d)
  
  # color_scheme_set("viridis")
  # plot <- mcmc_dens_chains(d, bw = "SJ") + theme_fagus()
  plot <- mcmc_areas(d, bw = "SJ", rhat = rhat, area_method = "equal height") + theme_fagus()
  
  ggsave(paste0(path, "/", basename_cmdstanfit, "_varname.", varname[1], ".pdf"), plot, device = "pdf", height = n_vars * 1.2, width = 6)
  
  return(plot)
}


## plotPredictions --------------------------------
# cmdstanfit  <- tar_read("priorsim_test")
# cmdstanfit  <- tar_read("fit_test")
# data_stan_priors <- tar_read("data_stan_priors")
# draws <- tar_read("draws_test") ## this is here as an option for plotting draw objects if the fit has NaNs in generated quantities
# path  <- tar_read("dir_fit")
plotPredictions <- function(cmdstanfit, data_stan_priors, draws = NULL, check = c("prior", "posterior"), path) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  data <- data_stan_priors$y
  Longdata <- attr(data_stan_priors, "Long")
  grp <- with(Longdata, interaction(as.integer(as.factor(obsid)), stage, substr(tax, 1, 1)))
  
  if(match.arg(check) == "prior") {
    
    if (is.null(draws)) {
      Sim <- cmdstanfit$draws(variables = "y_prior_sim", format = "draws_matrix")
    } else {
      Sim <- draws$y_prior_sim
    }
    
  } else if (match.arg(check) == "posterior") {
    
    if (is.null(draws)) {
      Sim <- cmdstanfit$draws(variables = "y_sim", format = "draws_matrix")
      Fix <- cmdstanfit$draws(variables = c("J_fix", "A_fix", "B_fix", "ba_fix")) %>%
        # subset_draws(draw = which(isconverged)) %>%
        as_draws_rvars() # ## For some reason, only extraction as array first and then as_draws_rvars() restores the desired data_structure!
      
    } else {
      Sim <- draws$y_hat_offset
      
      ## Fix <- 
      
    }
  }
  
  densplots <- list("predictions" = bayesplot::ppc_dens_overlay_grouped(log1p(data), log1p(Sim), group = grp))
  
  if (match.arg(check) == "posterior") {
    
    n_species <- 2
    convertVectorToDrawsList <- function(x) {
      lapply(1:length(x), function(i) as_draws_rvars(x[i]) )
    }
    
    fix_draws <- lapply(Fix, function(Rvar) lapply(1:n_species,
                                                   function(i) do.call(function(...) bind_draws(... , along = "draw"), convertVectorToDrawsList(Rvar[,i, drop = T]))
                                                   )
                        )
    fix_draws <- as_draws(lapply(fix_draws, function(f) do.call(cbind, lapply(f, function(l) l$x))))
    fix_draws <- thin_draws(fix_draws, thin = 10)
    
    M <- as_draws_matrix(fix_draws) ## enforce proper naming for plot methods
    fixplot <- bayesplot::mcmc_areas_ridges(log(M))
    
    densplots <- c(densplots, list("equilibria" = fixplot))
  }
  
  plotname <- paste(names(densplots), check, sep = "_")
  
  mapply(function(p, n) ggsave(paste0(path, "/", basename_cmdstanfit, "_", n, ".png"), p, device = "png", width = 15, height = 10),
         densplots, plotname)
  ## cowplot::plot_grid(densplot, fixdensplot, labels = c("States", "Equilibria"), ncol = 1) #  axis = "b", align = "h"
  
  return(densplots)
}


## plotSensitivity ------------------------------------------------------------
# cmdstanfit <- tar_read(fit_test)
# include <- tar_read()
# path  <- tar_read("dir_publish")
plotSensitivity <- function(cmdstanfit, include, measure = "cjs_dist", path) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  senssequence <- powerscale_sequence(cmdstanfit,
                                      variables = include,
                                      log_prior_fn = extract_log_prior, # require(priorsense)
                                      div_measure = measure)
  
  plot_powerscale <- powerscale_plot_dens(senssequence,
                                          variables = names(senssequence$base_draws)[1:(length(senssequence$base_draws)-3)]) ## These are the variable names in "include", but with indices.
  
  ggsave(paste0(path, "/", basename_cmdstanfit, "_", "sens_powerscale", ".pdf"), plot_powerscale, width = 42, height = 7)
  
  return(plot_powerscale)
}


## plotStates --------------------------------
# States <- tar_read("States")[[1]]
# allstatevars <- tar_read(basalareaname)
# path  <- tar_read("dir_publish")
# basename  <- tar_read("basename_fit")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotStates <- function(States,
                       allstatevars = basalareaname,
                       path, basename, color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  allstatevars <- intersect(as.character(unique(States$var)), allstatevars)
  
  statelabel <- c(ba_init = "Initial state",
                  ba_fix = "Equilibrium state",
                  ba_frac_init = "Initial fraction of the total basal area",
                  ba_frac_fix = "Equilibrium fraction of the total basal area",
                  
                  ba_fix_ko_b = "Equilibrium without b",
                  ba_fix_ko_s = "Equilibrium without s",
                  ba_fix_ko_2_b = "Equilibrium without b of others",
                  ba_fix_ko_2_s = "Equilibrium without s effect on others",
                  ba_fix_ko_b_l_r = "Equilibrium without the respective other",
                  ba_fix_ko_b_l_r_ko = "K.O. state",
                  
                  ba_fix_switch_b = "Equilibrium with switched b",
                  ba_fix_switch_c_a = "Equilibrium with switched c_A",
                  ba_fix_switch_c_b = "Equilibrium with switched c_B",
                  ba_fix_switch_c_j = "Equilibrium with switched c_J",
                  ba_fix_switch_b_c_b = "Equilibrium with switched b and c_B",
                  ba_fix_switch_b_c_a_c_b_h = "Equilibrium with switched overstory parameters",
                  ba_fix_switch_g = "Equilibrium with switched g",
                  ba_fix_switch_h = "Equilibrium with switched h",
                  ba_fix_switch_l = "Equilibrium with switched l",
                  ba_fix_switch_l_r = "Equilibrium with switched l and r",
                  ba_fix_switch_g_l_r_s = "Equilibrium with switched understory parameters",
                  ba_fix_switch_s = "Equilibrium with switched s")
  
  States <- States[!is.na(States$value),]
  
  T_major <- pivot_wider(States[1:6], names_from = "var", values_from = "value") %>%
    na.omit() %>% ## implicit NAs appear through completion by pivot_wider
    mutate(major_fix = as.logical(major_fix), major_init = as.logical(major_init))
  
  plot_major <- ggplot(T_major, aes(x = major_fix, y = ba_fix, col = tax, fill = tax)) +
    geom_violin(trim = T, col = "black", scale = "width") +
    
    ## scale_y_continuous(trans = "log10", n.breaks = 25) + # ggallin::pseudolog10_trans
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(base = 10, sides = "l", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) +
    themefun() +
    theme(axis.title.x = element_blank()) +
    theme(panel.grid.minor = element_blank()) + ## !!! remove the minor gridlines
    
    ggtitle("Equilibrium BA by taxon and whether Fagus ultimately has majority") +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color)
  # geom_jitter(position = position_jitter(0.2))
  
  #### When
  # whenvar <- c("ba_init", "ba_fix")
  # T_when <- filter(States, var %in% whenvar) %>% # filter(States, str_starts(var, "ba")) %>%
  #   rename(when = var) %>%
  #   mutate(when = factor(when, levels = whenvar)) %>%
  #   group_by(when, loc, draw) %>%
  #   mutate(diff_ba = value[tax == "Fagus"] - value[tax == "other"]) %>%
  #   ungroup()
  # 
  # plot_when <- ggplot(T_when, aes(x = tax, y = value, col = tax, fill = tax)) + # without facet wrap: ggplot(T_when, aes(x = when, y = value, col = tax, fill = tax))
  #   geom_violin(trim = T, col = "black", scale = "width") +
  #   
  #   # geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 4, fill = "transparent", scale = "width", data = T_when[T_when$is_loc_p10_draw,]) +
  #   geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 3, fill = "transparent", scale = "width", data = T_when[T_when$is_loc_median_draw,]) +
  #   # geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 2, fill = "transparent", scale = "width", data = T_when[T_when$is_loc_p90_draw,]) +
  #   
  #   ## scale_y_continuous(trans = "log10", n.breaks = 25) + # ggallin::pseudolog10_trans
  #   scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
  #                 labels = scales::trans_format("log10", scales::math_format(10^.x))
  #                 ) +
  #   annotation_logticks(base = 10, sides = "l", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) +
  #   themefun() +
  #   theme(axis.title.x = element_blank()) +
  #   theme(panel.grid.minor = element_blank()) + ## !!! remove the minor gridlines
  #   
  #   facet_grid(rows = . ~ when, labeller = labeller(when = statelabel)) +
  #   
  #   labs(y = "basal area [m^2 ha^-1]") +
  #   
  #   scale_color_manual(values = color) +
  #   scale_fill_manual(values = color)
  
  
  # Scatter_when <- States %>%
  #   filter(var %in% whenvar) %>% # filter(States, str_starts(var, "ba")) %>%
  #   filter(!is.na(value)) %>%
  #   rename(when = var) %>%
  #   mutate(when = factor(when, levels = whenvar)) %>%
  #   dplyr::select(tax, loc, value, when, draw) %>%
  #   pivot_wider(id_cols = c("draw", "when", "loc"), names_from = "tax", values_from = "value", names_prefix = "ba_")
  # 
  #   ## For adding density colours to points
  #   # group_by(when) %>%
  #   # mutate(denscol = densCols(x = log10(ba_other), y = log10(ba_Fagus),
  #   #                          nbin = 4000, colramp = colorRampPalette(c("#DEDEDE", "black"))))
  # 
  # plot_scatter_when <- ggplot(Scatter_when, aes(x = ba_other, y = ba_Fagus)) +
  #   
  #   geom_hex(bins = 100) +
  #   scale_fill_gradient(low = "#DDDDDD", high = "#000000", trans = "sqrt") +
  #   geom_abline(slope = 1, intercept = 0, linetype = 3) +
  #   # annotate(geom = 'text',  label = 'f(x) = x',
  #   #          x = diff(range(Scatter_when$ba_other)) * 0.001 + min(Scatter_when$ba_other), y = diff(range(Scatter_when$ba_Fagus)) * 0.001 + min(Scatter_when$ba_Fagus),
  #   #          size = 4, angle = 45) +
  #   
  #   ## For adding density colours directly to points, with col = denscol
  #   # geom_point(size = 0.1) + ## alpha = 0.1
  #   # scale_color_identity() +
  #   ## Other density viz:
  #   # geom_density_2d() + ## contour
  #   # geom_smooth(method='lm', formula = ba_Fagus ~ ba_other) +
  #   facet_grid(rows = . ~ when, labeller = labeller(when = statelabel)) +
  #   labs(y = "Fagus", x = "other", title = "Specific states basal area [m^2 ha^-1]") +
  #   
  #   ## scale_y_continuous(trans = "log10", n.breaks = 25) + # ggallin::pseudolog10_trans
  #   scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6),
  #                 labels = scales::trans_format("log10", scales::math_format(10^.x))
  #                 ) +
  #   scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6),
  #                 labels = scales::trans_format("log10", scales::math_format(10^.x))\
  #                 ) +
  #   annotation_logticks(base = 10, sides = "lb", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) +
  #   themefun() +
  #   theme(panel.grid.minor = element_blank()) + ## !!! remove the minor gridlines
  #   theme(legend.position = c(0.1, 0.65), legend.background = element_rect(fill = "transparent"))
  
  
  #### Main
  threevar <- c("ba_init", "ba_fix", "ba_fix_switch_s")
  T_main <- filter(States, var %in% threevar) %>%
    rename(when = var) %>%
    mutate(when = factor(when, levels = threevar)) %>%
    group_by(when, loc, draw) %>%
    mutate(diff_ba = value[tax == "Fagus"] - value[tax == "other"]) %>%
    ungroup()
  
  plot_main <- ggplot(T_main, aes(x = tax, y = value, col = tax, fill = tax)) +
    geom_violin(trim = T, col = "black", scale = "width") +
    
    # geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 4, fill = "transparent", scale = "width", data = T_main[T_main$is_loc_p10_draw,]) +
    geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 3, fill = "transparent", scale = "width", data = T_main[T_main$is_loc_median_draw,]) +
    # geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 2, fill = "transparent", scale = "width", data = T_main[T_main$is_loc_p90_draw,]) +
    
    ## scale_y_continuous(trans = "log10", n.breaks = 25) + # ggallin::pseudolog10_trans
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(base = 10, sides = "l", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) +
    themefun() +
    theme(axis.title.x = element_blank()) +
    theme(panel.grid.minor = element_blank()) + ## !!! remove the minor gridlines
    
    facet_grid(rows = . ~ when, labeller = labeller(when = statelabel)) +
    
    labs(y = "basal area [m^2 ha^-1]") +
    
    scale_color_manual(values = color) +
    scale_fill_manual(values = color)
  
  Scatter_main <- States %>%
    filter(var %in% threevar) %>% # filter(States, str_starts(var, "ba")) %>%
    filter(!is.na(value)) %>%
    rename(when = var) %>%
    mutate(when = factor(when, levels = threevar)) %>%
    dplyr::select(tax, loc, value, when, draw) %>%
    pivot_wider(id_cols = c("draw", "when", "loc"), names_from = "tax", values_from = "value", names_prefix = "ba_")
  
  plot_scatter_main <- ggplot(Scatter_main, aes(x = ba_other, y = ba_Fagus)) +
    geom_hex(bins = 100) +
    scale_fill_gradient(low = "#DDDDDD", high = "#000000", trans = "sqrt") +
    geom_abline(slope = 1, intercept = 0, linetype = 3) +
    facet_grid(rows = . ~ when, labeller = labeller(when = statelabel)) +
    labs(y = "Fagus", x = "other", title = "Specific states basal area [m^2 ha^-1]") +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(base = 10, sides = "lb", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) +
    themefun() +
    theme(panel.grid.minor = element_blank()) + ## !!! remove the minor gridlines
    theme(legend.position = c(0.1, 0.65), legend.background = element_rect(fill = "transparent"))
  
  
  #### Supplementary
  fourvar <- c("ba_fix", "ba_fix_switch_s", "ba_fix_c_j", "ba_fix_switch_l_r",
               "ba_fix_switch_g", "ba_fix_switch_b_c_b")
  T_supp <- filter(States, var %in% fourvar) %>%
    rename(when = var) %>%
    mutate(when = factor(when, levels = fourvar)) %>%
    group_by(when, loc, draw) %>%
    mutate(diff_ba = value[tax == "Fagus"] - value[tax == "other"]) %>%
    ungroup()
  
  plot_supp <- ggplot(T_supp, aes(x = tax, y = value, col = tax, fill = tax)) +
    geom_violin(trim = T, col = "black", scale = "width") +
    
    # geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 4, fill = "transparent", scale = "width", data = T_supp[T_supp$is_loc_p10_draw,]) +
    geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 3, fill = "transparent", scale = "width", data = T_supp[T_supp$is_loc_median_draw,]) +
    # geom_violin(aes(x = tax, y = value), trim = T, col = "black", linetype = 2, fill = "transparent", scale = "width", data = T_supp[T_supp$is_loc_p90_draw,]) +
    
    ## scale_y_continuous(trans = "log10", n.breaks = 25) + # ggallin::pseudolog10_trans
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(base = 10, sides = "l", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) +
    themefun() +
    theme(axis.title.x = element_blank()) +
    theme(panel.grid.minor = element_blank()) + ## !!! remove the minor gridlines
    
    facet_grid(rows = . ~ when, labeller = labeller(when = statelabel)) +
    
    labs(y = "basal area [m^2 ha^-1]") +
    
    scale_color_manual(values = color) +
    scale_fill_manual(values = color)
  
  Scatter_supp <- States %>%
    filter(var %in% fourvar) %>% # filter(States, str_starts(var, "ba")) %>%
    filter(!is.na(value)) %>%
    rename(when = var) %>%
    mutate(when = factor(when, levels = fourvar)) %>%
    dplyr::select(tax, loc, value, when, draw) %>%
    pivot_wider(id_cols = c("draw", "when", "loc"), names_from = "tax", values_from = "value", names_prefix = "ba_")
  
  plot_scatter_supp <- ggplot(Scatter_supp, aes(x = ba_other, y = ba_Fagus)) +
    geom_hex(bins = 100) +
    scale_fill_gradient(low = "#DDDDDD", high = "#000000", trans = "sqrt") +
    geom_abline(slope = 1, intercept = 0, linetype = 3) +
    facet_grid(rows = . ~ when, labeller = labeller(when = statelabel)) +
    labs(y = "Fagus", x = "other", title = "Specific states basal area [m^2 ha^-1]") +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(base = 10, sides = "lb", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) +
    themefun() +
    theme(panel.grid.minor = element_blank()) + ## !!! remove the minor gridlines
    theme(legend.position = c(0.1, 0.65), legend.background = element_rect(fill = "transparent"))
  
  
  #### All
  T_all <- filter(States, var %in% allstatevars) %>%
    rename(gq = var) %>%
    mutate(gq = factor(gq, levels = allstatevars))
  
  plot_all <- ggplot(T_all, aes(x = tax, y = value, col = tax, fill = tax)) +
    geom_violin(trim = T, col = "black", scale = "width") +
    facet_wrap(~ gq, labeller = labeller(gq = statelabel)) +
    ggtitle("BA") +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    
    ## scale_y_continuous(trans = "log10", n.breaks = 25) + # ggallin::pseudolog10_trans
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(base = 10, sides = "l", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) +
    themefun() +
    theme(panel.grid.minor = element_blank())## !!! remove the minor gridlines
  
  # plot_diff <- ggplot(T_when, aes(x = when, y = diff_ba)) +
  #   geom_violin(trim = FALSE, col = "black", scale = "width") +
  #   ggtitle("log_BA_Fagus - log_BA_other at equilibirum and at initial time") +
  #   scale_color_manual(values = color) +
  #   themefun()
  
  plots <- list(plot_states_major = plot_major,
                # plot_states_when = plot_when,
                # plot_states_scatter_when = plot_scatter_when,
                # plot_states_main = plot_main,
                # plot_states_scatter_main = plot_scatter_main,
                # plot_states_supp = plot_supp,
                # plot_states_scatter_supp = plot_scatter_supp,
                plot_states_all = plot_all) # plot_states_diff = plot_diff
  
  mapply(function(p, n) ggsave(paste0(path, "/", basename, "_", n, ".pdf"), p, device = "pdf", width = 10, height = 8),
         plots, names(plots))
  
  # stateplotgrid <- cowplot::plot_grid(plot_when + theme(legend.position = "none"), plot_scatter_when, labels = c("(A)", "(B)"),  align = "h", axis = "rl",  nrow = 2, rel_heights = c(1.4, 1))
  # ggsave(paste0(path, "/", basename, "_plot_states_combined", ".png"), stateplotgrid, device = "png", width = 8, height = 10)
  # ggsave(paste0(path, "/", basename, "_plot_states_combined", ".pdf"), stateplotgrid, device = "pdf", width = 8, height = 10)
  
  stateplotgrid_main <- cowplot::plot_grid(plot_main + theme(legend.position = "none"), plot_scatter_main, labels = c("(A)", "(B)"),  align = "h", axis = "rl",  nrow = 2, rel_heights = c(1.3 ,1))
  ggsave(paste0(path, "/", basename, "_plot_states_main", ".png"), stateplotgrid_main, device = "png", width = 10, height = 9)
  ggsave(paste0(path, "/", basename, "_plot_states_main", ".pdf"), stateplotgrid_main, device = "pdf", width = 10, height = 9)
  
  stateplotgrid_supp <- cowplot::plot_grid(plot_supp + theme(legend.position = "none"), plot_scatter_supp, labels = c("(A)", "(B)"),  align = "h", axis = "rl",  nrow = 2, rel_heights = c(1.3 ,1))
  ggsave(paste0(path, "/", basename, "_plot_states_supp", ".png"), stateplotgrid_supp, device = "png", width = 14, height = 8)
  ggsave(paste0(path, "/", basename, "_plot_states_supp", ".pdf"), stateplotgrid_supp, device = "pdf", width = 14, height = 8)
  
  return(NULL) # return(plots)
}



## plotConditional_resampling --------------------------------
## This resamples the posterior two times based on the frequencies of Fagus majority freq_major and 1-freq_major (referring to the frequencies of subpopulations).
## - The resampling was done with replacement proportional to the frequency of predominance.
## - both posteriors Fagus_major and others_major have the same number of samples, although the actual ratio might be different.
## - The sum of the two sample sets is not equal to the posterior.
##
## LaTeX caption: \caption{Correlations of two separate posterior parameter distributions, conditioned on whether Fagus (green) or others (yellow) more frequently are the major population (in basal area) at the equilibrium state. To get the two conditional posteriors, all parameter combinations (MCMC draws) were resampled with replacement, proportional to the frequency of the subpopulations where either Fagus or others were in the majority. E.g., one MCMC draw that would result in 70\% of the subpopulations having Fagus in the majority at equilibrium, would have a resampling probability of 0.7 for the conditional posterior for Fagus and vice versa. Note that due to separate resampling with replacements, the two posteriors have an equal number of samples, although the actual ratio of \emph{Fagus'} to \emph{others'} predominance is \textbf{X} (Table \textbf{X}).}
##
# parname  <- tar_read("parname_plotorder")
# cmdstanfit  <- tar_read("fit_test")
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotConditional_resampling <- function(cmdstanfit, parname, path,
                                       color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  isconverged <- cmdstanfit$draws(variables = "converged_fix", format = "draws_matrix") %>%
    apply(1, all)
  
  message("Dropped ", sum(!isconverged), " draw(s), because the simulations have not converged to fixpoint.")
  
  freq_major <- cmdstanfit$draws(variables = "major_fix", format = "draws_matrix") %>%
    subset_draws(draw = which(isconverged)) %>%
    rowMeans()
  
  if(all(freq_major == 1) | all(freq_major == 0)) {
    warning("Either one species has the majority in all draws in all clusters. Thus, resampling is not possible")
    return(NULL)
  }
  
  ## Compare random effects
  ## also: K_loc_log_raw * sigma_k_loc
  # draws_L_major <- cmdstanfit$draws(variables = "L_loc") %>%
  #   as_draws_rvars() %>%
  #   resample_draws(weights = freq_major)
  # draws_L_minor <- cmdstanfit$draws(variables = "L_loc") %>%
  #   as_draws_rvars() %>%
  #   resample_draws(weights = 1-freq_major)
  # 
  # diff <- draws_L_major[[1]] - draws_L_minor[[1]]
  
  renameAll <- function(draws, suffix = "major") {
    
    name <- posterior::variables(draws)
    name_new <- paste(name, suffix, sep = "_")
    posterior::variables(draws) <- name_new
    
    return(draws)
  }
  
  
  draws_par <- cmdstanfit$draws(variables = parname) %>%
    subset_draws(draw = which(isconverged)) %>%
    as_draws_rvars() ## For some reason, only extraction as array first and then as_draws_rvars() restores the desired data_structure!
  
  draws_par_weighted_major <- posterior::resample_draws(draws_par, weights = freq_major, method = "stratified") %>%
    renameAll(suffix = "major")
  
  draws_par_weighted_minor <- posterior::resample_draws(draws_par, weights = 1 - freq_major, method = "stratified") %>%
    renameAll(suffix = "minor")
  
  
  d <- posterior::bind_draws(draws_par_weighted_major, draws_par_weighted_minor)
  names_order <- names(d) %>% sort()
  d <- d[names_order]
  
  # if ("L_loc" %in% parname) {
  #   ## just a smaller subset for random effects
  #   d <- lapply(d, function(x) x[20:40,])
  # }
  
  ## Marginal plots
  d_1 <- lapply(d, function(i) i[1]) %>% as_draws_array()
  d_2 <- lapply(d, function(i) i[2]) %>% as_draws_array()
  
  plots_parameters_conditional <- list(
    Fagus.sylvatica = bayesplot::mcmc_areas_ridges(d_1) + themefun(),
    other = bayesplot::mcmc_areas_ridges(d_2) + themefun()
  )
  
  plotgrid <- cowplot::plot_grid(plotlist = plots_parameters_conditional, ncol = 2, labels = names(plots_parameters_conditional))
  ggsave(paste0(path, "/", basename_cmdstanfit, "_plot_conditional", ".pdf"), plotgrid, dev = "pdf", height = 20, width = 24)
  
  
  ## Pairs plot
  D <- d %>%
    gather_draws(`.*`[i], regex = T) %>%
    ungroup() %>%
    mutate(major = if_else(str_ends(.variable, "_major"), "Fagus_major", "other_major")) %>%
    mutate(tax = fct_recode(as.character(i), "Fagus" = "1", "other" = "2")) %>%
    mutate(parameter = str_extract(.variable, "([a-z]|c_.+)_log")) %>%
    mutate(parameter = factor(parameter, levels = parname)) %>%
    pivot_wider(id_cols = c(".draw", "major"), names_from = c("parameter", "tax"), values_from = ".value")
  
  ## Custom density for colorscale
  plotDensity <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_density(..., color = "black") +
      scale_fill_manual(values = color)
    }
  
  pairsplot <- ggpairs(D,
                       mapping = aes(col = major, fill = major),
                       columns = match(paste0(rep(parname, each = 2), c("_Fagus", "_other")), colnames(D)),
                       diag = list(continuous = plotDensity),
                       upper = list(continuous = wrap("cor", size = 3.3)),
                       lower = list(continuous = wrap("points", alpha = 0.1, size = 0.6))) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    themefun() +
    theme(panel.spacing = unit(0.1, "lines"))
  
  ggsave(paste0(path, "/", basename_cmdstanfit, "_pairs_conditional", ".png"), pairsplot, device = "png", height = 26, width = 26)
  
  return(c(plots_parameters_conditional, 'pairs' = pairsplot))
}


## plotConditional --------------------------------
## This plots a pairs correlation plot and does not resample as above,
## but simply colors the posterior samples by freq_major.
## - here, the posterior is the actual posterior.
# parname  <- tar_read("parname_plotorder")
# cmdstanfit  <- tar_read("fit")[[1]]
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotConditional <- function(cmdstanfit, parname, conditional = T, path, 
                            color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  isconverged <- cmdstanfit$draws(variables = "converged_fix", format = "draws_matrix") %>%
    apply(1, all)
  
  message("Dropped ", sum(!isconverged), " draw(s), because the simulations have not converged to fixpoint.")
  
  freq_major <- cmdstanfit$draws(variables = "major_fix", format = "draws_matrix") %>%
    subset_draws(draw = which(isconverged)) %>%
    rowMeans()
  
  ismajor <- freq_major > 0.5
  
  d <- cmdstanfit$draws(variables = parname) %>%
    subset_draws(draw = which(isconverged)) %>%
    as_draws_rvars() ## For some reason, only extraction as array first and then as_draws_rvars() restores the desired data_structure!
  
  
  ## Pairs plot
  D <- d %>%
    gather_draws(`.*`[i], regex = T) %>% suppressWarnings() %>% ## package tidyr warns about using deprecated gather_()
    ungroup() %>%
    mutate(tax = fct_recode(as.character(i), "Fagus" = "1", "other" = "2")) %>%
    mutate(parameter = str_extract(.variable, "([a-z]|c_.+)_log")) %>%
    mutate(parameter = factor(parameter, levels = parname)) %>%
    pivot_wider(id_cols = c(".draw"), names_from = c("parameter", "tax"), values_from = ".value") %>%
    bind_cols(freq_major = freq_major, major = if_else(ismajor, "Fagus", "other"))
  
  ## Custom density for colorscale
  plotDensity <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_density(..., color = "black") +
      scale_fill_manual(values = color)
  }
  
  plotPoints <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(mapping = aes(col = 1-freq_major), ..., alpha = 0.1, size = 0.6) + ## everything else but col gets inherited!
      scale_color_gradient(low = color[1], high = color[2], guide = "colourbar", aesthetics = c("colour"))
  }
  
  plotHex <- function(data, mapping, ...) {
    ggplot(data = data, mapping) +
      geom_hex(bins = 30) + ## everything else but col gets inherited!
      scale_fill_continuous(trans = "reverse")
  }
  
  plotStats <- function(data, mapping, ...) {
    
    printStats <- function(x, y) {
      m <- lm(y ~ x)
      slopetext <- paste(formatC(coef(m)[2], digits = 3, format = "f"), signif_stars(summary(m)$coefficients[2,4]))
      
      cor <- stats::cor.test(x, y)
      cortext <- formatC(as.numeric(cor$estimate), digits = 3, format = "f")
      cortext <- str_c(cortext, signif_stars(cor$p.value), "\n")
      
      return(c(cortext, slopetext))
    }
    
    gtext <- ggally_statistic(data, mapping, printStats, title = c("cor", "slope"))
    
    return(gtext)
  }
  
  
  pairsplot <- ggpairs(D,
                       mapping = aes(fill = !!if (conditional) D$major else NULL), ## discrete mapping for densities
                       columns = match(paste0(rep(parname, each = 2), c("_Fagus", "_other")), colnames(D)),
                       diag = list(continuous = plotDensity),
                       upper = list(continuous = plotStats), # wrap("cor", size = 3.3)
                       lower = list( continuous = if (conditional) plotPoints else plotHex ) ) + ## wrap("points", alpha = 0.1, size = 0.6)
    
    # scale_color_manual(values = color) +
    { if(conditional) scale_fill_manual(values = color)  else scale_fill_continuous(trans = "reverse") } +
    themefun() +
    theme(panel.spacing = unit(0.1, "lines"))
  
  
  ggsave(paste0(path, "/", basename_cmdstanfit, "_pairs_conditional", ".png"), pairsplot, device = "png", height = 26, width = 26)
  
  return(list('pairs' = pairsplot))
}


## plotPairs --------------------------------
# parname  <- tar_read("parname_plotorder")
# cmdstanfit  <- tar_read("fit_test")
# cmdstanfit  <- tar_read("fit")[[1]]
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotPairs <- function(cmdstanfit, parname,
                      path = tar_read("dir_publish"), color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  varname_draws <- cmdstanfit$metadata()$stan_variables
  parname <- intersect(parname, varname_draws)
  
  d <- cmdstanfit$draws(variables = parname)
  
  ## Extract info from variable names
  # D <- d %>%
  #   as_draws_rvars() %>%
  #   gather_draws(`.*`[i], regex = T) %>% suppressWarnings() %>% ## package tidyr warns about using deprecated gather_()
  #   ungroup() %>%
  #   mutate(tax = fct_recode(as.character(i), "Fagus" = "1", "other" = "2")) %>%
  #   mutate(parameter = str_extract(.variable, "([a-z]|c_.+)_log")) %>%
  #   mutate(parameter = factor(parameter, levels = parname)) %>%
  #   pivot_wider(id_cols = c(".draw"), names_from = c("parameter", "tax"), values_from = ".value")
  
  nutsparam <- nuts_params(cmdstanfit)
  # logposterior <- log_posterior(cmdstanfit)
  # mcmc_nuts_divergence(nutsparam, logposterior)
  
  color_scheme_set("viridis")
  pairsplot <- mcmc_pairs(d,
                          # pars = c("..."), transform = list(sigma = "log"),
                          diag_fun = "dens",
                          # diag_args = "",
                          off_diag_fun = "hex", # "scatter",
                          # off_diag_args = list(size = 0.3, alpha = 0.5), ## for "scatter"
                          # condition = pairs_condition(nuts = "lp__"),
                          # lp = logposterior,
                          np = nutsparam,
                          np_style = pairs_style_np(div_color = "red",
                                                    div_shape = 4,
                                                    div_size = 1,
                                                    div_alpha = 1,
                                                    td_color = "black",
                                                    td_shape = 3,
                                                    td_size = 1,
                                                    td_alpha = 1)) # + theme_fagus()
  
  ggsave(paste0(path, "/", basename_cmdstanfit, "_pairs_", parname[1], ".png"), pairsplot, device = "png", height = 28, width = 28)
  ggsave(paste0(path, "/", basename_cmdstanfit, "_pairs_", parname[1], ".pdf"), pairsplot, device = "pdf", height = 28, width = 28)
  
  return(list('pairs' = pairsplot))
}


## plotTrace --------------------------------
# parname  <- tar_read("parname_plotorder")
# cmdstanfit  <- tar_read("fit_test")
# cmdstanfit  <- tar_read("fit")[[1]]
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotTrace <- function(cmdstanfit, parname, path,
                      color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  d <- cmdstanfit$draws(variables = parname)
  
  color_scheme_set("viridis")
  traceplot <- mcmc_trace(d) + # # pars = c("..."), transformations = list(sigma = "log")
    theme_fagus()
  
  ggsave(paste0(path, "/", basename_cmdstanfit, "_trace", ".png"), traceplot, device = "png", height = 15, width = 20)
  
  return(list('trace' = traceplot))
}


## plotContributions --------------------------------
# parname  <- c(tar_read("parname_plotorder"), b_c_b = "b_c_b_log")
# contribname <- paste("sum_ko", rep(1:2, each = length(parname)), names(parname), "fix", sep = "_")
# cmdstanfit  <- tar_read("fit")[[1]]
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")

plotContributions <- function(cmdstanfit, parname, path, contribution = c("sum_ko", "sum_ko_avg", "sum_ko_prop", "sum_switch"), plotlog = FALSE,
                              color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  n_species <- 2 ## is also used in functions below
  parorder <- names(parname)
  parname <- str_remove(setdiff(parname, "k_log"), "_log")
  plotprop <- match.arg(contribution) == "sum_ko_prop"
  if (plotprop & plotlog) warning("Log-scale option was discarded because it does not make sense with proportions.")
  
  contribname <- if (match.arg(contribution) == "sum_ko_prop") { paste("sum_ko", rep(1:n_species, each = length(parname)), "prop", parname, "fix", sep = "_") }
  else if (match.arg(contribution) == "sum_ko") { paste("sum_ko", rep(1:n_species, each = length(parname)), parname, "fix", sep = "_") }
  else if (match.arg(contribution) == "sum_ko_avg") { paste("sum_ko", rep(1:n_species, each = length(parname)), parname, "fix_avg", sep = "_") }
  else if (match.arg(contribution) == "sum_switch") { paste("sum_switch", parname, "fix", sep = "_") }
  
  
  C <- cmdstanfit$draws(variables = contribname) %>%
    as_draws_rvars()
  
  convertVectorToDrawsList <- function(x) {
    lapply(1:length(x), function(i) as_draws_rvars(x[i]) )
  }
  
  fix_draws <- lapply(C, function(Rvar) lapply(1:n_species,
                                               function(i) do.call(function(...) bind_draws(... , along = "draw"), convertVectorToDrawsList(Rvar[,i, drop = T]))
                                               )
                      )
  fix_draws <- as_draws(lapply(fix_draws, function(f) do.call(cbind, lapply(f, function(l) l$x))))
  # fix_draws <- thin_draws(fix_draws, thin = 10)
  
  M <- as_draws_matrix(fix_draws) ## enforce proper naming for plot methods
  I <- bayesplot::mcmc_intervals_data(M, point_est = "median", prob = 0.5, prob_outer = 0.8) %>%
    mutate(p = parameter,
           parameter = str_extract(p, "(?<=_)(sKoEnvB|b_c_b|c_a|c_b|c_j|[bghlrs]{1})(?=_)"),
           kotax = suppressWarnings( fct_recode(str_extract(p, "(?<=_)(\\d)(?!=_)"), "Fagus" = "1", "other" = "2") ),
           tax = suppressWarnings( fct_recode(str_extract(p, "(\\d+)(?!.*\\d)"), "Fagus" = "1", "other" = "2") ), # the last number in the string
           reciprocal = as.character(kotax) != as.character(tax), # there might be different level sets
           stage = fct_collapse(parameter, "J" = c("c_j", "r", "l", "s", "sKoEnvB"), "A" = c("g", "c_a"), "B" = c("c_b", "b", "h", "b_c_b"))
           ) %>%
    mutate(stage = ordered(stage, c("J", "A", "B"))) %>%
    mutate(stagepos = as.numeric(as.character(fct_recode(stage, "1" = "J", "5.5" = "A", "7.5" = "B")))) %>%
    mutate(parameter = ordered(parameter, parorder)) %>%
    mutate(parameter = fct_reorder(parameter, as.numeric(stage))) %>%
    
    ## Letter Positions
    group_by(reciprocal) %>%
    mutate(xletterpos_h = max(hh) * 0.85, ## for !plotprop
           xletterpos_l = min(ll) * 1.01) %>% ## for plotprop
    
    arrange(stage, parameter)
  
  # plot_contributions <- bayesplot::mcmc_areas_ridges(M)
  
  pos <- position_nudge(y = (as.integer(I$tax) - 1.5) * 0.3)
  
  plot_contributions <- ggplot(I, aes(x = m, y = parameter, yend = parameter,
                                      color = tax, group = stage)) + 
    geom_linerange(aes(xmin = l, xmax = h), size = 2.6, position = pos) +
    # geom_segment(aes(x = l, xend = h), size = 3, lineend = "round", position = pos) + ## unfortunately the lineend
    geom_segment(aes(x = ll, xend = hh), size = 1.2, lineend = "round", position = pos) +
    geom_point(color = "black", position = pos, size = 1.7) +
    coord_flip() +
    geom_vline(xintercept = if (plotprop) 1 else 0, linetype = 3, size = 0.6, col = "#222222") +
    # geom_hline(yintercept = c(4.5, 6.5), linetype = 3, size = 0.6, col = "#222222") + ## lines before g/h
    geom_hline(yintercept = c(5, 7), linetype = 3, size = 0.6, col = "#222222") + ## lines through g/h
    
    { if (!plotprop) geom_text(aes(y = stagepos, x = xletterpos_h, label = stage), size = 9, col = "#222222") } +
    { if (plotprop) geom_text(aes(y = stagepos, x = xletterpos_l, label = stage), size = 9, col = "#222222") } +
    
    { if (match.arg(contribution) != "sum_switch") facet_wrap(~reciprocal, # ~kotax
                                                              scales = "free",
                                                              labeller = labeller(kotax = function(kotax) paste(kotax, "demographic rates"),
                                                                                  reciprocal = function(reciprocal) if_else(reciprocal == 'TRUE', "indirect", "direct"))) } +
    themefun() +
    scale_color_manual(values = color) +
    theme(axis.title.x = element_blank()) +
    
    ## Only for log-scale
    { if (!plotprop & plotlog) scale_x_continuous(trans = pseudo_log_trans(sigma = 0.1, base = 10)) } + ## breaks = c(-10^(1:3), 10^(1:3)) #ggallin::pseudolog10_trans, breaks = scales::trans_breaks("log10", function(x) 10^x, n = 8) , labels = scales::trans_format("log10", scales::math_format(10^.x)))# labels = scales::trans_format("log10", scales::math_format(10^.x))
    ## logticks do not seem work for the pseudo-log close to zero!
    # { if (!plotprop & plotlog) annotation_logticks(base = 10, sides = "l", scaled = T, short = unit(1, "mm"), mid = unit(2, "mm"), long = unit(2.5, "mm"), colour = "black", size = 0.25) } +
    # { if (!plotprop & plotlog) theme(panel.grid.minor = element_blank()) } + ## !!! remove the minor gridlines
    
    { if (plotprop) labs(x = "Average yearly increment in proportion to the total basal area increment [ ]", y = "Parameter", title = "Yearly propotional contributions to the basal") } +
    { if (!plotprop) labs(x = "Cumulated rate of basal area increment [m2 ha-1 yr-1]", y = "Parameter", title = "Contributions to the basal area") }
  
  ggsave(paste0(path, "/", basename_cmdstanfit, "_plot_contributions_", contribution, if(plotlog) "_log_" else "_", length(contribname), ".pdf"),
         plot_contributions, dev = "pdf", height = 8, width = 12)
  
  return(NULL) # return(plot_contributions)
}


## plotPredominant --------------------------------
# States  <- tar_read("States")[[1]]
# majorname  <- tar_read(majorname)
# path  <- tar_read("dir_publish")
# basename  <- tar_read("basename_fit")[[1]]
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotPredominant <- function(States, majorname,
                            path = tar_read("dir_publish"), basename, color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  States <- States[!is.na(States$value),]
  
  majorname <- intersect(unique(States$var), majorname)
  
  majorlabel <- c(major_init = "Initial state",
                  major_fix = "Equilibrium state",
                  
                  major_fix_switch_l = "switched l",
                  major_fix_switch_l_r = "switched l and r",
                  major_fix_switch_c_j = "switched c_J",
                  major_fix_switch_s = "switched s",
                  major_fix_switch_g = "switched g",
                  major_fix_switch_g_l_r_s = "switched understory parameters",
                  major_fix_switch_c_a = "switched c_A",
                  major_fix_switch_h = "switched h",
                  major_fix_switch_b = "switched b",
                  major_fix_switch_c_b = "switched c_B",
                  major_fix_switch_b_c_b = "switched b and c_B",
                  major_fix_switch_b_c_a_c_b_h = "switched overstory parameters")
  
  
  D <- States[1:6] %>%
    
    dplyr::filter(var %in% majorname) %>%
    
    mutate(value_ = 1 - value) %>%
    mutate(value = if_else(tax == "Fagus", value, value_)) %>% ## more efficient than doing the if_else inside groups
    mutate(stage = fct_collapse(var,
                                J = c("major_fix_switch_l_r", "major_fix_switch_c_j", "major_fix_switch_s", "major_fix_switch_g"),
                                A = c("major_fix_switch_c_a", "major_fix_switch_h"),
                                B = c("major_fix_switch_b_c_b"))
    ) %>%
    
    ## average over draws
    group_by(var, tax, stage) %>%
    summarize(value = mean(value), sd = sd(value), q10 = quantile(value, 0.1), q90 = quantile(value, 0.9)) %>%
    ungroup()
  
  plot_predominant <- ggplot(D, aes(x = factor(var, levels = rev(names(majorlabel))),
                                    y = value,
                                    group = stage,
                                    # ymin = if_else(tax == "Fagus", q10, as.numeric(NA)),
                                    # ymax = if_else(tax == "Fagus", q90, as.numeric(NA)),
                                    fill = tax)) + 
    geom_col(position = position_fill(reverse = TRUE), width = 0.4, color = "black", size = 0.5) +
    # geom_errorbar(width = .2, position = "identity") +
    geom_text(aes(label = if_else(value >= 0.07, paste0(sprintf("%.1f", value * 100), "% ", if_else(value >= 0.5, as.character(tax), "")),"")),
              position = position_fill(reverse = TRUE, vjust = 0.5), colour = "white", size = 4) +
    ylab("% predominant subpopulations") +
    scale_y_continuous(labels = scales::percent) +
    # scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    themefun() +
    coord_flip() +
    theme(axis.title.y = element_blank()) +
    theme(legend.position = c(-0.4, 0))
  
  
  ggsave(paste0(path, "/", basename, "_plot_predominance", ".pdf"),
         plot_predominant,
         dev = "pdf", height = 8, width = 5)
  
  return(plot_predominant)
}


## plotMarginal --------------------------------
# Marginal <- tar_read(Marginal_env)
# parname <- str_to_sentence(tar_read(parname_plotorder))
# basename  <- tar_read("basename_fit_env")
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotMarginal <- function(Marginal, parname,
                         path, basename, color = c("#208E50", "#FFC800"), themefun = theme_fagus) {

  M <- Marginal %>%
    filter(.variable %in% parname) %>%
    mutate(tax = fct_recode(as.character(tax), "Fagus" = "1", "other" = "2")) %>%
    group_by(.variable, tax) %>%
    mutate(medianofmedians = median(median))
  
  plot <- ggplot(M, aes(x = .value, group = interaction(.variable, tax), fill = tax)) + 
    geom_density(trim = T, alpha = 0.8) +
    scale_fill_manual(values = color) +
    geom_vline(aes(xintercept = medianofmedians, color = tax)) +
    scale_color_manual(values = color) +
    facet_grid(rows = vars(.variable), scales = "free") + # cols = vars(tax)
    themefun() +
    theme(axis.title.x = element_blank())

  ggsave(paste0(path, "/", basename, "_plot_marginal", ".pdf"),
         plot, dev = "pdf", height = 14, width = 6)
  
  return(NULL) # return(plot)
}


## plotEnvironmental --------------------------------
# surfaces <- tar_read(surface_environmental_env)
# basename <- tar_read(basename_fit_env)
# path <-  tar_read(dir_publish)
plotEnvironmental <- function(surfaces = surface_environmental_env, binaryname = "major_fix", commonscale = FALSE, removevar = NULL,
                              basename = tar_read("basename_fit_env"), path = tar_read("dir_publish"), color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  
  #### Handling of the binary border
  parname_surfaces <- sapply(surfaces, function(s) attr(s, "parname")) ## Will be used to name plots below
  i_binary <- which(binaryname == parname_surfaces)
  Binary <- if(isTRUE(i_binary >= 1)) surfaces[[i_binary]] else NULL
  Binary$z <- if(!is.null(Binary)) round(Binary$z) else NULL
  bnotnull <- !is.null(Binary)
  
  if(bnotnull) {
    
    Binary_1 <- dplyr::filter(Binary, z == 1) %>% ## is ordered, so that slicing should return some central location
      slice(round(nrow(.)*0.5)) %>%
      bind_cols(label = "Fagus predominant")
  
  } else{
    
    warning("plotEnvironmental(): surface with binaryname is not in provided surfaces.")
  }
  
  
  if (!isTRUE(commonscale)) {
    
    #### 1. CASE: !commonscale, i.e. differernt plots per variable on a plotgrid
    
    plotE <- function(D) {
      
      name_x <- attr(D, "name_x")
      name_y <- attr(D, "name_y")
      parname <- attr(D, "parname")
      taxon <- attr(D, "taxon")
      
      parstart_reverse <- c("S_", "C_j_", "C_a_", "C_b_",
                            "sum_ko_1_s_fix", "sum_ko_2_s_fix",
                            "sum_ko_1_c_j_fix", "sum_ko_2_c_j_fix",
                            "sum_ko_1_c_a_fix", "sum_ko_2_c_a_fix",
                            "sum_ko_1_c_b_fix", "sum_ko_2_c_b_fix")
      
      isanyreversepar <- str_starts(parname, fixed(parstart_reverse)) %>% any()
      
      direction <- if (isanyreversepar) 1 else -1
      
      plot <- ggplot(D, aes_string(x = name_x, y = name_y, z = "z")) +
        geom_raster(aes(fill = z)) +
        metR::geom_contour2(mapping = aes(z = z, label = round(..level.., 3)), col = "white") +
        # scale_color_manual(values = color) +
        scale_fill_viridis_c(direction = direction) +
        
        { if(bnotnull) geom_contour(mapping = aes_string(x = name_x, y = name_y, z = "z"),
                                    data = Binary, bins = 2, col = "black", linetype = 2, size = 0.8, inherit.aes = F) } + 
        { if(bnotnull) geom_text(data = Binary_1, mapping = aes_string(x = name_x, y = name_y, label = "label"), col = "black") } +
        
        themefun() +
        ggtitle(paste(parname, taxon)) +
        scale_y_reverse() ## invert water level scale, consistent with Ökogramm.
      
      return(plot)
    }
    
    plots <- lapply(surfaces, plotE)
    names(plots) <- parname_surfaces
    plotgrid <- cowplot::plot_grid(plotlist = plots, ncol = 2)
    
    ggsave(filename = paste0(path, "/", basename, "_plotgrid_environmental_", format(Sys.time(), "%H.%M.%S"), ".pdf"),
           plot = plotgrid, device = "pdf", width = 15, height = 4 * length(plots), limitsize = FALSE)
    
  } else {
    
    #### 2. CASE: commonscale, i.e. common plot with facets
    
    name_x <- attr(surfaces[[1]], "name_x")
    name_y <- attr(surfaces[[2]], "name_y")
    
    surfaces <- lapply(surfaces, function(S) bind_cols(S, tax = attr(S, "taxon"), variable = attr(S, "parname")))
    
    D <- bind_rows(surfaces, .id = "bindingid") %>%
      suppressMessages( mutate(tax = fct_recode(as.character(tax), "Fagus" = "1", "others" = "2", "both" = "0")) ) %>%
      dplyr::filter(!variable %in% removevar)
    
    scale_fill_div <- function(...) scale_fill_gradient2(low = color[1], mid = "white", high = color[2], midpoint = 0, ...)
    
    plots <- ggplot(D, aes_string(x = name_x, y = name_y, z = "z")) +
      geom_raster(aes(fill = z)) +
      metR::geom_contour2(mapping = aes(z = z, label = round(..level.., 3)), col = "white") +
      scale_fill_div() +
      
      { if(bnotnull) geom_contour(mapping = aes_string(x = name_x, y = name_y, z = "z"),
                                  data = Binary, bins = 2, col = "black", linetype = 2, size = 0.8, inherit.aes = F) } + 
      { if(bnotnull) geom_text(data = Binary_1, mapping = aes_string(x = name_x, y = name_y, label = "label"), col = "black") } +
      
      themefun() +
      scale_y_reverse() + ## invert water level scale, consistent with Ökogramm.
      { if( length(unique(D$tax)) > 1 ) facet_wrap(~ variable + tax, ncol = 2) else facet_wrap(~ variable, ncol = 2) }

    
    ggsave(filename = paste0(path, "/", basename, "_plot_environmental_commonscale_", format(Sys.time(), "%H.%M.%S"), ".pdf"),
           plot = plots, device = "pdf", width = 10, height = 2.2 * length(surfaces), limitsize = FALSE)
    
  }
  
  return(plots)
}


## plotPoly --------------------------------
# Surfaces <- tar_read(Surfaces_poly_env)
# Environmental <- tar_read(Environmental_env)
# basename <- tar_read(basename_fit_env)
# path <-  tar_read(dir_publish)
plotPoly <- function(Surfaces, Environmental = NULL,
                     basename = tar_read("basename_fit_env"), path = tar_read("dir_publish"),
                     color = c("#208E50", "#FFC800"), themefun = theme_fagus) {

  name_x <- attr(Surfaces, "name_x")
  name_y <- attr(Surfaces, "name_y")
  
  if (!is.null(Environmental)) {
    
    E <- Environmental %>%
      dplyr::filter(str_starts(.variable, "sum_ko_") & str_ends(.variable, "fix")) %>% 
      rename(parname = .variable) %>%
      group_by_at(c("tax", "parname", "loc", name_x, name_y)) %>%
      dplyr::summarize(value = mean(.value, na.rm = T), .groups = "drop") %>%
      mutate(isdirectcontribution = str_extract(parname, "[1-2]") == tax) %>%
      filter(isdirectcontribution) %>%
      mutate(taxon = fct_recode(as.character(tax), Fagus = "1", others = "2")) %>%
      mutate(parname = str_extract(parname, "(c_[jab]|[a-z])(?=_fix)")) %>%
      mutate(parname = paste0(parname, "_log")) %>%
      filter(parname %in% unique(Surfaces$parname))
    
    }
  
  scale_color_div <- function(...) scale_color_gradient2(low = color[1], mid = "white", high = color[2], midpoint = 0, ...)
  
  plot <- ggplot(Surfaces, aes_string(x = name_x, y = name_y)) +
    # geom_raster(aes(fill = z)) +
    
    { if (!is.null(Environmental)) geom_jitter(data = E, mapping = aes(color = value), width = 0.1, height = 0.2, alpha = 0.6) } + # width = 0.3, height = 0.3, size = 0.5
    { if (!is.null(Environmental)) scale_color_div() } +
    
    metR::geom_contour2(data = Surfaces, mapping = aes(z = z, label = round(..level.., 3)), # 
                        colour = "gray30", global.breaks = F, margin = unit(rep(4, 4), "pt"), label.placer = label_placer_flattest(), lineend = "round", skip = 1) +
    

    
    facet_grid(rows = vars(parname), cols = vars(taxon)) +
    themefun() +
    
    scale_y_reverse()

  ggsave(filename = paste0(path, "/", basename, "_plot_poly", ".pdf"),
         plot = plot, device = "pdf", width = 10, height = 40, limitsize = FALSE)
  
  return(NULL) # return(plot)
}


## plotBinary --------------------------------
# Environmental <- tar_read("Environmental_env")
# parname <- str_to_sentence(tar_read(parname_plotorder))
# parname <- tar_read(contribname_env)
# tar_load(fit_environmental_env_binomial)
# fit_bin <- fit_environmental_env_binomial[[sapply(fit_environmental_env_binomial, function(x) attr(x, "par") == "major_fix") %>% which()]]
# themefun  <- tar_read("themefunction")
# color  <- tar_read("twocolors")
# basename  <- tar_read("basename_fit_env")
# path  <- tar_read("dir_publish")
plotBinary <- function(Environmental, parname, fit_bin = NULL, binarythreshold = 0.5,
                       facetscale = c("free", "free_x", "free_y"), plotlog = F,
                       path = tar_read("dir_publish"), basename = tar_read("basename_fit_env"),
                       color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  
  binaryname <- if (is.null(fit_bin)) "major_fix" else attr(fit_bin, "par")
  facetscale <- match.arg(facetscale)
  
  B <- Environmental %>%
    ungroup() %>%
    dplyr::filter(.variable == binaryname) %>%
    rename(binary = .value)
    
  if (is.null(fit_bin)) {
    
    B %<>% dplyr::select(binary, loc, .draw)
  
  } else {
    
    ## replace the drawn variable with an environmental prediction, but make locs unique first
    B %<>% distinct(loc, .keep_all = T)
    
    if (is(fit_bin, "cv.glmnet")) {
     error("plotBinary(): fit of type cv.glmnet not implemented.") ## see above how to predict with cv.glmnet
    } else {
      p <- predict(fit_bin, newdata = B, type = "response") %>% c()
    }
    
    t <- as.integer(p > binarythreshold)
    
    B %<>%
      dplyr::select(loc) %>%
      bind_cols(binary = t, binary_pred = p)
  }

  Environmental %<>%
    ungroup() %>%
    dplyr::filter(.variable %in% parname) %>%
    
    dplyr::filter(.variable != binaryname) %>%
    
    { if (is.null(fit_bin)) left_join(., B, by = c("loc", ".draw")) else left_join(., B, by = "loc")  } %>%
    dplyr::filter(!is.na(binary))
    # rename(chain = .chain, iteration = .iteration, draw = .draw, variable = .variable, value = .value) %>%
  
  
  E <- Environmental %>%
    rename(value = .value) %>%
    mutate(binary = as.logical(binary)) %>%
    mutate(binary = if_else(tax == 1, binary, !binary)) %>%
    

    ## aggregate by draws by loc first, to reduce the variation to the environmental variation
    group_by(tax, .variable, binary, loc) %>%
    summarize(value = mean(value, na.rm = T)) %>%
    ungroup() %>%
    
    group_by(tax, .variable, binary) %>%
    summarize(ll = quantile(value, 0.1), l = quantile(value, 0.25),
              m = median(value),
              u = quantile(value, 0.75), uu = quantile(value, 0.9)) %>%
    ungroup() %>%
    
    mutate(stage = suppressWarnings( fct_collapse(.variable,
                                                  "J" = c("L_log", "R_log", "S_log", "C_j_log"),
                                                  "A" = c("G_log", "C_a_log"),
                                                  "B" = c("H_log", "B_log", "C_b_log"),
                                                  other_level = "JAB")
                                   )
           ) %>%
    
    mutate(tax = fct_recode(as.character(tax), Fagus = "1", others = "2")) %>% 
    rename(p = .variable) %>% 
    mutate(parameter = if_else(str_starts(p, "sum_"),
                               str_extract(p, "(?<=_)(sKoEnvB|b_c_b|c_a|c_b|c_j|[bghlrs]{1})(?=_)"),
                               p), ## distinguish whether contribution or parameter
           parameter = factor(parameter, levels = unique(c(parname, parameter, p))), # for the order
           kotax = suppressWarnings( fct_recode(str_extract(p, "(?<=_)(\\d)(?!=_)"), "Fagus" = "1", "others" = "2") ),
           # avg = str_detect(p, "_avg"),
           predominance = as.factor(if_else(binary, "wins", "looses")),
           reciprocal = as.character(kotax) != as.character(tax)) %>%
    mutate(stage = ordered(stage, c("J", "A", "B"))) %>%
    mutate(stagepos = as.numeric(as.character(fct_recode(stage, "1" = "J", "5.5" = "A", "7.5" = "B")))) %>%
    mutate(labelxpos = if_else(m > 0, ll - 0.2, uu + 0.2)) %>%
    filter(!Vectorize(isTRUE)(reciprocal)) ## removes indirect effect!
  
  
  pos <- position_dodge(width = 0.8)
  
  plot_binary <- ggplot(E, aes(x = m, y = tax, color = tax, group = predominance)) +  # y = parameter, 
    
    geom_linerange(aes(xmin = l, xmax = u), size = 2.6, position = pos) +
    geom_linerange(aes(xmin = ll, xmax = uu), size = 1.2, position = pos) +
    geom_point(color = "black", position = pos, size = 0.8) +
    geom_text(aes(x = labelxpos, label = predominance), position = pos) +
    coord_flip() +
    # geom_vline(xintercept = if (plotprop) 1 else 0, linetype = 3, size = 0.6, col = "#222222") +
    # geom_hline(yintercept = c(5, 7), linetype = 3, size = 0.6, col = "#222222") + ## lines through g/h
    
    facet_wrap(~parameter, scales = facetscale,
               labeller = labeller(tax = function(tax) paste(tax))) +
    
    { if (plotlog) scale_x_continuous(trans = pseudo_log_trans(sigma = 0.1, base = 10)) } +
    
    themefun() +
    scale_color_manual(values = color)
  
  ggsave(paste0(path, "/", basename, "_plot_binary_", parname[1], ".pdf"),
         plot_binary, dev = "pdf", height = 10, width = 10)
  
  return(plot_binary)
}


## plotTrajectories --------------------------------
# Trajectories <- tar_read("Trajectories_avg_test")
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotTrajectories <- function(Trajectories, thicker = FALSE, path, basename, plotpdf = FALSE,
                             color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  Trajectories %<>%
    group_by(loc, tax, stage, draw) %>%
    filter(any(isconverged)) %>%
    # filter(time >= 3) %>%
    filter(time < 3000) %>%
    # mutate(time_shifted = time - time_fix, time_log = log(time)) %>% ## shifts the time, so that trajectories are aligned by fixpoint
    mutate(grp = interaction(loc, tax, draw), tax = as.factor(tax)) %>%
    ungroup() %>%
    mutate(tax = fct_recode(as.character(tax), "Fagus" = "1", "other" = "2")) %>%
    
    ## Cut stage B above a certain value for pretty facet ylims
    filter(!(stage == "B" & abundance > 400))
  
  aes_lines <- aes(x = time, y = abundance, group = grp, col = tax)
  
  plot <- ggplot(Trajectories, aes_lines) +
    { if(thicker) geom_line(size = 0.5, alpha = 0.06) else geom_line(size = 0.2, alpha = 0.05) } +
    facet_wrap(~stage, scales = "free_y",
               labeller = labeller(stage = c(J = "J (count [ha^-1])",
                                             A = "A (count [ha^-1])",
                                             B = "B (basal area [m^2 ha^-1])"))) +
    coord_trans(y = "sqrt", x = "sqrt") + # coord_trans(y = "log2", x = "log2") + # 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 12)) +
    # scale_y_continuous(trans = "pseudo_log", n.breaks = 8) +
    scale_color_manual(values = color) +
    xlab("time (years)") +
    themefun()
  
  if(is.null(basename)) basename <- "Model"
  ggsave(paste0(path, "/", basename, "_trajectories", ".png"), plot, device = "png", height = 5.2, width = 12)
  if(plotpdf) ggsave(paste0(path, "/", basename, "_trajectories", ".pdf"), plot, device = "pdf", height = 5.2, width = 12)
  
  return(NULL) # return(plot)
}


## animateTrajectories --------------------------------
# plot_trajectories <- tar_read("plot_trajectories_avg_test")
# path  <- tar_read("dir_publish")
animateTrajectories <- function(plot_trajectories, path, basename) {
  
  animation <- plot_trajectories +
    geom_point(size = 0.01) + ## draws points at the lineends
    transition_reveal(time, range = c(0, 3000))
  
  if(is.null(basename)) basename <- "Model"
  gganimate::animate(animation, duration = 8, fps = 15, width = 1100, height = 500, renderer = ffmpeg_renderer(format = "mp4"))
  anim_save(paste0(path, "/", basename, "_trajectories_animation", ".mp4"))
  
  return(NULL) # return(animation)
}



# ————————————————————————————————————————————————————————————————————————————————— #
# Posterior tests         -----------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## testSensitivity ------------------------------------------------------------
# cmdstanfit <- tar_read(fit_test)
# include <- tar_read(parname)
# path  <- tar_read("dir_publish")
## For CJSdist, we consider an ad hoc threshold ≥ 0.05 to be indicative of sensitivity. For a normal distribution, this corresponds to the mean differing by approximately more than 0.3 standard deviations,
## or the standard deviation differing by a factor greater than approximately 0.3, when the power-scaling factor is changed by a factor of two.
testSensitivity <- function(cmdstanfit, include, measure = "cjs_dist", path) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  sensitivity <- powerscale_sensitivity(cmdstanfit,
                                        variables = include,
                                        log_prior_fn = extract_log_prior, # require(priorsense)
                                        div_measure = measure)
  
  write.csv(sensitivity[[1]], paste0(path, "/", basename_cmdstanfit, "_sensitivity.csv"))
  
  return(sensitivity)
}


