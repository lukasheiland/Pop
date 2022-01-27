# ————————————————————————————————————————————————————————————————————————————————— #
# Posterior verbs         -----------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## extract — read draws. Returns a draws, or a stanfit object.
## format — formats the extracted draws for later use.
## summarize — tabulating the posterior, Returns a data.framoid
## generate — generate quantities from the posterior.
## test — statistical tests on the posteriors.
## plot - plots. Side effect: file. Returns a plot object or a list of plot objects.



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
# stanfit  <- tar_read("stanfit")
# stanfit  <- tar_read("stanfit_test")
# helpers_exclude  <- tar_read("helpers_exclude")
extractDraws <- function(stanfit, exclude = helpers_exclude) {
  
  draws <- rstan::extract(stanfit, pars = exclude, include = F)
  
  return(draws)
}


## formatLoc --------------------------------
## helper for all variables that have structure array[N_locs] int/real
# cmdstanfit <- tar_read("fit_test")
formatLoc <- function(name, locmeans = FALSE, cmdstanfit_ = cmdstanfit, data_stan_priors_ = data_stan_priors) {
  n_locs <- data_stan_priors_$N_locs
  
  Draws <- cmdstanfit_$draws(variables = name, format = "draws_matrix")
  n_draws <- if (locmeans) 1 else nrow(Draws)
  n_i <- cmdstanfit_$metadata()$stan_variable_sizes[[name]][2]
  if(is.na(n_i)) n_i <- 1
  
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


## formatStates --------------------------------
# cmdstanfit <- tar_read("fit_test")
# data_stan_priors <- tar_read("data_stan_priors")
formatStates <- function(cmdstanfit, data_stan_priors) {
  
  statename <- c("major_init", "major_fix", "ba_init", "ba_fix", "ba_fix_ko_s")
  
  States <- lapply(statename, formatLoc, cmdstanfit_ = cmdstanfit, data_stan_priors_ = data_stan_priors)
  States <- lapply(States, function(S) if( length(unique(S$i)) == 1 ) bind_rows(S, within(S, {i <- 2})) else S )
  States %<>%
    bind_rows() %>%
    mutate(tax = factor(c("Fagus", "other")[i]))
  
  States$value[States$value == 9] <- NA
  
  return(States)
}



# ————————————————————————————————————————————————————————————————————————————————— #
# Summarize posterior         -----------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## summarizeFit --------------------------------
# cmdstanfit <- tar_read("fit_test")
# path <- tar_read("dir_publish")
summarizeFit <- function(cmdstanfit, exclude = NULL, path) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  allpar <- cmdstanfit$metadata()$stan_variables
  includepar <- setdiff(allpar, exclude)
  summary <- cmdstanfit$summary(includepar)
  
  write.csv(summary, paste0(path, "/", basename_cmdstanfit, "_summary.csv"))
  
  head(summary, 20) %>%
    as.data.frame() %>%
    print()
  
  drawconverged <- cmdstanfit$draws(variables = "converged_fix", format = "draws_matrix") %>%
    apply(1, all)
  
  if (all(drawconverged)) {
    message("All clusters in all draws have converged to the fix point, i.e. the population trajectories ran into an equilibrium.")
  } else {
    message("For ", sum(!drawconverged), " draws, not all of the clusters have converged to the fix point, i.e. not all population trajectories ran into an equilibrium.")
  }
  
  return(summary)
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
# cmdstanfit  <- tar_read("fit")
# data_stan_priors  <- tar_read("data_stan_priors")
# path  <- tar_read("dir_fit")
generateResiduals <- function(cmdstanfit, data_stan_priors, path) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")

  Sim <- cmdstanfit$draws(variables = "y_sim", format = "draws_matrix") %>% t()# matrix of observations simulated from the fitted model - row index for observations and colum index for simulations
  Sim[is.na(Sim)] <- 0
  y <- data_stan_priors$y
  y_hat <- cmdstanfit$draws(variables = "y_hat_rep", format = "draws_matrix") %>% apply(2, median, na.rm = T)
  y_hat[is.na(y_hat)] <- 0
  
  Longdata <- attr(data_stan_priors, "Long")
  grp <- with(Longdata, interaction(as.integer(as.factor(obsid)), stage, substr(tax, 1, 1)))
  
  residuals <- DHARMa::createDHARMa(simulatedResponse = Sim, observedResponse = y, fittedPredictedResponse = y_hat, integerResponse = T)

  png(paste0(path, "/", basename_cmdstanfit, "_", "DHARMa", ".png"), width = 1600, height = 1000)
  plot(residuals, quantreg = T, smoothScatter = F)
  dev.off()
  
  # residuals_grouped <- recalculateResiduals(residuals, group = grp)
  
  png(paste0(path, "/", basename_cmdstanfit, "_", "DHARMa_grouped", ".png"), width = 2200, height = 800)
  plot(residuals, form = grp, quantreg = T, smoothScatter = F)
  dev.off()
  
  return(residuals)
}


## generateTrajectories --------------------------------
# cmdstanfit <- tar_read("fit_test")
# parname <- tar_read("parname")
# data_stan_priors <- tar_read("data_stan_priors")
# locparname <- tar_read("parname_loc")

generateTrajectories <- function(cmdstanfit, data_stan_priors, parname, locparname = "state_init_log",
                                 time = c(seq(1, 491, by = 10), seq(500, 5000, by = 100)), thinstep = 1, usemean = FALSE) {
  
  parname <- setdiff(parname, c("phi_obs", "sigma_k_loc"))
  parname_sans_log <- gsub("_log$", "", parname)
  
  
  ## iterateModel --------------------------------
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
                                                        stage = rep(c("j", "a", "b"), each = n * length(times)),
                                                        tax = rep(1:n, each = length(times)),
                                                        time = times)
    
    return(State)
  }
  
  
  Draws <- cmdstanfit$draws(variables = c(parname, locparname), format = "draws_list") %>%
    thin_draws(thin = thinstep)
  
  # Makes pars to be a nested list/data.frame[[parameters]][[draws]]
  if(usemean) {
    
    M <- cmdstanfit$summary(variables = parname) %>%
      dplyr::select(variable, mean)
    parmeans <- sapply(parname, function(n) M$mean[str_starts(M$variable, n)], USE.NAMES = T, simplify = F)
    pars <- sapply(parname, function(n) { if(str_ends(n, "_log")) exp(parmeans[[n]]) else parmeans[[n]] }, USE.NAMES = T, simplify = F)
    names(pars) <- parname_sans_log
    pars <- list(pars)
    
  } else {
    
    D <- subset_draws(Draws, variable = parname) %>%
      as_draws_matrix() %>%
      as.data.frame() %>%
      dplyr::mutate(across(contains("_log["), exp))
    pars <- split(D, 1:nrow(D))
    pars <- lapply(pars, function(p) matrix(c(as.matrix(p)), byrow = F, nrow = data_stan_priors$N_species, dimnames = list(NULL, parname_sans_log)))
    pars <- lapply(pars, as.data.frame)
  }
  
  ## Nested simulations: locs/draws
  
  iterateModel_draws <- function(locpars, pars, time, usemean) {
    
    locpars <- lapply(locpars, function(x) if(usemean) mean(x) else as_draws_matrix(x))
    
    ## Assign local parameters to draws, adopt manually if necessary
    pars <- lapply(1:length(pars), function(i) within(pars[[i]], l <- c(locpars$L_loc[i,])))
    
    return( lapply(1:length(pars), function(i) iterateModel(c(locpars$state_init[i,]), pars = pars[[i]], time)) )
  }
  
  iterateLocs <- function(i, lp, p, t, um) {
    iterateModel_draws(locpars = lapply(lp, function(x) x[i,]), pars = p, time = t, usemean = um)
  }
  
  
  draws_loc <- subset_draws(Draws, variable = locparname) %>% # c("state_init_log", "L_loc")
    as_draws_rvars()
  draws_loc$state_init <- exp(draws_loc$state_init_log)
  n_locs <- data_stan_priors$N_locs
  
  sims <- future_sapply(1:n_locs, iterateLocs,
                        lp = draws_loc, p = pars, t = time, um = usemean, simplify = F) # a nested list[locs, draws] of matrices[times]
  sims <- lapply(sims, bind_rows, .id = "draw")
  sims <- bind_rows(sims, .id = "loc")
  
  sims %<>%
    group_by(loc, draw) %>%
    mutate(diff = abundance - c(NA, abundance[1:(n()-1)]),
           absdiff = abs(diff),
           isconverged = replace_na(absdiff <= data_stan_priors$tolerance_fix, F),
           time_fix = first(time[isconverged]),
           state_fix = first(abundance[isconverged])) %>%
    ungroup()
  
  return(sims)
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
  
  plotRidges <- function(startswith, fit = stanfit) {
    bayesplot::mcmc_areas_ridges(fit, pars = vars(starts_with(startswith, ignore.case = F))) +
      themefun()
  }
  
  usedmcmc <- "sample" == attr(stanfit, "stan_args")[[1]]$method
  
  # basename <- attr(stanfit, "model_name") %>%
  #   str_replace("-[1-9]-", "-x-")
  
  parname <- setdiff(stanfit@model_pars, exclude)
  parnamestart <- na.omit(unique(str_extract(parname, "^[a-z]_[ljab]"))) # Everything that starts with a small letter, and has the right index after that to be a meaningful parameter. (Small letter is important!)
  parname_sansprior <- parname[!grepl("prior$", parname)]
  
  
  traceplot <- rstan::traceplot(stanfit, pars = parname_sansprior, include = T)
  areasplot <- bayesplot::mcmc_areas(stanfit, area_method = "scaled height", pars = vars(!matches(c(exclude, "log_", "lp_", "prior")))) + themefun()
  ridgeplots <- parallel::mclapply(parnamestart, plotRidges, mc.cores = getOption("mc.cores", 7L))
  ridgeplotgrid <- cowplot::plot_grid(plotlist = ridgeplots)
  
  # parallelplot_c <- bayesplot::mcmc_parcoord(stanfit, pars = vars(starts_with(c("c_", "s_"))))
  # parallelplot_others <- bayesplot::mcmc_parcoord(stanfit, pars = vars(!matches(c(exclude, "c_", "log_", "phi_", "lp_", "s_", "_prior"))))
  
  plots <- list(traceplot = traceplot,
                ridgeplotgrid = ridgeplotgrid,
                areasplot = areasplot) # parallelplot_c = parallelplot_c, parallelplot_others = parallelplot_others,
  
  mapply(function(p, n) ggsave(paste0(path, "/", basename, "_", n, ".png"), p, device = "png"), plots, names(plots))
  
  if(usedmcmc) {
    
    png(paste0(path, "/", basename, "_", "pairsplot", ".png"), width = 2600, height = 2600)
    pairs(stanfit, pars = c(parname_sansprior, "lp__"), include = T)
    dev.off()
    
  }
  
  return(plots)
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
      Sim <- draws$y_hat_rep
      
      ## Fix <- 
      
    }
  }

  densplots <- list("predictions" = bayesplot::ppc_dens_overlay_grouped(log(data), log(Sim), group = grp))
  
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
# States <- tar_read("States_test")
# path  <- tar_read("dir_publish")
# basename  <- tar_read("basename_fit_test")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotStates <- function(States, path, basename, color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  T_major <- pivot_wider(States, names_from = "var") %>%
    mutate(major_fix = as.logical(major_fix), major_fix = as.logical(major_fix))

  plot_major <- ggplot(T_major, aes(x = tax, y = log(ba_fix), col = major_fix)) +
    geom_violin(trim = FALSE) +
    ggtitle("Equilibrium BA by taxon and whether Fagus ultimately has majority") +
    scale_color_manual(values = color) +
    themefun()
    # geom_jitter(position = position_jitter(0.2))
  
  
  T_when <- filter(States, str_starts(var, "ba")) %>%
    rename(when = var) %>%
    group_by(when, loc, draw) %>%
    mutate(diff_ba = value[tax == "Fagus"] - value[tax == "other"]) %>%
    ungroup()
  
  plot_when <- ggplot(T_when, aes(x = when, y = log(value), col = tax)) +
    geom_violin(trim = FALSE) +
    ggtitle("BA at equilibirum and at initial time") +
    scale_color_manual(values = color) +
    themefun()
  
  # plot_diff <- ggplot(T_when, aes(x = when, y = diff_ba)) +
  #   geom_violin(trim = FALSE) +
  #   ggtitle("log_BA_Fagus - log_BA_other at equilibirum and at initial time") +
  #   scale_color_manual(values = color) +
  #   themefun()
  
  plots <- list(plot_states_major = plot_major, plot_states_when = plot_when) # plot_states_diff = plot_diff
  
  mapply(function(p, n) ggsave(paste0(path, "/", basename, "_", n, ".png"), p, device = "png"), plots, names(plots))
  
  return(plots)
}



## plotConditional --------------------------------
# parname  <- tar_read("parname_sim")
# cmdstanfit  <- tar_read("fit_test")
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotConditional <- function(cmdstanfit, parname, path,
                            color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  isconverged <- cmdstanfit$draws(variables = "converged_fix", format = "draws_matrix") %>%
    apply(1, all)
  
  message("Dropped ", sum(!isconverged), " draw(s), because the simulations have not converged to fixpoint.")
  
  freq_major <- cmdstanfit$draws(variables = "major_fix", format = "draws_matrix") %>%
    subset_draws(draw = which(isconverged)) %>%
    rowMeans()
  
  
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
  
  draws_par_weighted_major <- posterior::resample_draws(draws_par, weights = freq_major, method = "simple") %>%
    renameAll(suffix = "major")
  
  draws_par_weighted_minor <- posterior::resample_draws(draws_par, weights = 1 - freq_major, method = "simple") %>%
    renameAll(suffix = "minor")
  
  
  d <- posterior::bind_draws(draws_par_weighted_major, draws_par_weighted_minor)
  names_order <- names(d) %>% sort()
  d <- d[names_order]
  
  if (parname %in% c("L_loc")) {
    ## just a smaller subset for random effects
    d <- lapply(d, function(x) x[20:40,])
  }
  
  d_1 <- lapply(d, function(i) i[1]) %>% as_draws_array()
  d_2 <- lapply(d, function(i) i[2]) %>% as_draws_array()
  
  plots_parameters_conditional <- list(
    Fagus.sylvatica = bayesplot::mcmc_areas_ridges(d_1) + themefun(),
    other = bayesplot::mcmc_areas_ridges(d_2) + themefun()
  )
  
  plotgrid <- cowplot::plot_grid(plotlist = plots_parameters_conditional, ncol = 2, labels = names(plots_parameters_conditional))
  ggsave(paste0(path, "/", basename_cmdstanfit, "_plot_conditional", ".png"), plotgrid, dev = "png", height = 20, width = 24)
  
  return(plots_parameters_conditional)
}


## plotContributions --------------------------------
# parname  <- tar_read("parname_sim")
# cmdstanfit  <- tar_read("fit_test")
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")

plotContributions <- function(cmdstanfit, parname, path, plotprop = FALSE,
                              color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  basename_cmdstanfit <- attr(cmdstanfit, "basename")
  
  parname <- str_remove(setdiff(parname, "k_log"), "_log")
  contribname <- if (plotprop) paste0("sum_ko_prop_", parname, "_fix") else paste0("sum_ko_", parname, "_fix")
  
  C <- cmdstanfit$draws(variables = contribname) %>%
    as_draws_rvars()
  
  n_species <- 2
  
  convertVectorToDrawsList <- function(x) {
    lapply(1:length(x), function(i) as_draws_rvars(x[i]) )
  }
  
  fix_draws <- lapply(C, function(Rvar) lapply(1:n_species,
                                               function(i) do.call(function(...) bind_draws(... , along = "draw"), convertVectorToDrawsList(Rvar[,i, drop = T]))
  )
  )
  fix_draws <- as_draws(lapply(fix_draws, function(f) do.call(cbind, lapply(f, function(l) l$x))))
  fix_draws <- thin_draws(fix_draws, thin = 10)
  
  M <- as_draws_matrix(fix_draws) ## enforce proper naming for plot methods
  I <- bayesplot::mcmc_intervals_data(M, point_est = "median", prob = 0.5, prob_outer = 0.9) %>%
    mutate(p = parameter,
           parameter = str_extract(p, "(?<=_)([bghlrs]{1}|c_a|c_b|c_j)(?=_)"),
           tax = fct_recode(str_extract(p, "(\\d+)(?!.*\\d)"), "Fagus sylvatica" = "1", "other" = "2"),
           stage = fct_collapse(parameter, "J" = c("c_j", "r", "l", "s"), "A" = c("g", "c_a"), "B" = c("c_b", "b", "h"))
    ) %>%
    mutate(stage = ordered(stage, c("J", "A", "B"))) %>%
    arrange(stage, parameter)
  
  
  # plot_contributions <- bayesplot::mcmc_areas_ridges(M)
  
  pos <- position_nudge(y = (as.integer(I$tax) - 1.5) * 0.2)
  
  plot_contributions <- ggplot(I, aes(x = m, y = fct_reorder(parameter, as.numeric(stage)), color = tax, group = stage)) + 
    geom_linerange(aes(xmin = l, xmax = h), size = 2, position = pos) +
    geom_linerange(aes(xmin = ll, xmax = hh), position = pos) +
    geom_point(color = "black", position = pos) +
    coord_flip() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = color) +
    themefun()
  
  ggsave(paste0(path, "/", basename_cmdstanfit, "_plot_contributions", if(plotprop) "_prop" else "", ".png"), plot_contributions, dev = "png", height = 30, width = 10)
  
  return(plot_contributions)
}


## plotTrajectories --------------------------------
# Trajectories <- tar_read("Trajectories_test")
# path  <- tar_read("dir_publish")
# color  <- tar_read("twocolors")
# themefun  <- tar_read("themefunction")
plotTrajectories <- function(Trajectories, thicker = FALSE, path, basename,
                             color = c("#208E50", "#FFC800"), themefun = theme_fagus) {
  
  Trajectories %<>%
    group_by(loc, tax, stage, draw) %>%
    filter(any(isconverged)) %>%
    mutate(time_shifted = time - time_fix) %>%
    mutate(grp = interaction(loc, tax, draw), tax = as.factor(tax)) %>%
    mutate(abundance_log = log(abundance))
  
  aes_lines <- aes(x = time_shifted, y = abundance_log, group = grp, col = tax)
  
  plot <- ggplot(Trajectories, aes_lines) +
    { if(thicker) geom_line(size = 0.5, alpha = 0.1) else geom_line(size = 0.1, alpha = 0.01) } +
    facet_wrap(~stage) +
    scale_color_manual(values = color) +
    themefun()
  
  if(is.null(basename)) basename <- "Model"
  ggsave(paste0(path, "/", basename, "_equilines", ".png"), plot, dev = "png", height = 14, width = 20)
  
  return(plot)
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


