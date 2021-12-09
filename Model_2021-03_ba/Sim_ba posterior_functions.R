# ————————————————————————————————————————————————————————————————————————————————— #
# Posterior simulations ------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## iterateModel --------------------------------
# cmdstanfit  <- tar_read("fit_test")
# data_stan_priors <- tar_read("data_stan_priors")
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
  r <- pars$r # Length n_species vector of input rates
  s <- pars$s # Length n_species vector of shading rates
  
  ba_a_upper <- pars$ba_a_upper
  ba_a_avg <- pars$ba_a_avg
  
  ## Set the count variables
  n <- length(r) # no. of species
  
  ## 
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


## simulateTrajectories --------------------------------
# cmdstanfit <- tar_read("fit_test_pq")
# parname <- tar_read("parname")
# data_stan_priors <- tar_read("data_stan_priors")

simulateTrajectories <- function(cmdstanfit, data_stan_priors, parname, time = seq(1, 501, by = 10), thinstep = 1, mean = FALSE) {
  
  parname <- setdiff(parname, c("phi_obs", "sigma_l"))
  parname_sans_log <- gsub("_log$", "", parname)
  
  Draws <- cmdstanfit$draws(variables = c(parname, "state_init_log"), format = "draws_list") %>%
    thin_draws(thin = thinstep)
  
  # Makes pars to be a nested list/data.frame[[parameters]][[draws]]
  if(mean) {

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
  
  State_init <- subset_draws(Draws, variable = "state_init_log") %>%
    as_draws_rvars() %$%
    state_init_log %>%
    exp()
  
  iterateModel_draws <- function(state_init, pars, time) {
    S <- as_draws_matrix(state_init)
    lapply(1:nrow(S), function(i) iterateModel(c(S[i,]), pars = pars[[i]], time))
  }
  
  # State_init <- State_init[10:11,]
  # time <- 1:300
  
  sims <- future.apply::future_apply(State_init, 1, iterateModel_draws, pars = pars, time = time, simplify = F) # a nested list[locs, draws] of matrices[times]
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


## plotTrajectories --------------------------------
# Trajectories <- tar_read("Trajectories")

plotTrajectories <- function(Trajectories) {
  
  Trajectories %<>%
    # filter(time > 1) %>%
    mutate(grp = interaction(loc, tax, draw), tax = as.factor(tax))
  
  aes_lines <- aes(x = time, y = abundance, group = grp, col = tax, alpha = 0.001)
  
  plot <- ggplot(Trajectories, aes_lines) +
    geom_line() +
    facet_wrap(~stage)
  
  return(plot)
}



# ————————————————————————————————————————————————————————————————————————————————— #
# Notes            ------------------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## One time series for playing ---------------------------------------------------------------
# S <- tar_read(summary_test)
# data <- tar_read(data_stan_priors)
# 
# pars <- list(
#   b = filter(S, str_starts(variable, "b_log\\["))$mean,
#   c_a = filter(S, str_starts(variable, "c_a_log\\["))$mean,
#   c_b = filter(S, str_starts(variable, "c_b_log\\["))$mean,
#   c_j = filter(S, str_starts(variable, "c_j_log\\["))$mean,
#   g = filter(S, str_starts(variable, "g_logit\\["))$mean,
#   h = filter(S, str_starts(variable, "h_logit\\["))$mean,
#   l = filter(S, str_starts(variable, "l_log\\["))$mean,
#   r = filter(S, str_starts(variable, "r_log\\["))$mean,
#   s = filter(S, str_starts(variable, "s_log\\["))$mean
#   )
# 
# pars_trans <- c(lapply(pars, exp), list(ba_a_upper = data$ba_a_upper, ba_a_avg = data$ba_a_avg, n_species = 2))
# times <- 2:200
# 
# Sim1 <- calcModel(rep(c(1000, 100, 100), each = 2), times = times, pars = pars_trans)
# 
# matplot(Sim1[,-1], type = "b", ylab="N",
#         pch = rep(c("J", "A", "B"), each = pars_trans$n_species),
#         col = 1:pars_trans$n_species, xlab = "time") # log='y'

