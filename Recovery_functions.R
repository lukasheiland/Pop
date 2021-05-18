# Library -----------------------------------------------------------------

library(GGally)


# Plot predicted (generated quantities) vs. true values -----------------------------------------
plotPredictedVsTrue <- function(draws, data, modelname = "ba-rect") {
  
  truepars <- attr(data, "pars")
  
  if(modelname == "ba-rect") {
    
    par(mfrow = c(2, 2))
        
    # predictions per life stages
    plot(c(draws$y_hat_rep[100,,3,,data$i_j]) ~ c(data$y[,3,,data$i_j]))
    abline(0, 1)
    plot(c(draws$y_hat_rep[100,,3,,data$i_a]) ~ c(data$y[,3,,data$i_a]))
    abline(0, 1)
    plot(c(draws$y_hat_rep[1,,3,,data$i_b]) ~ c(data$y[,3,,data$i_b]))
    abline(0, 1)
    
    # plot(c(draws$y_sim[1,,3,,data$i_j]) ~ c(data$y[,3,,data$i_j]))
    # abline(0, 1)
    # plot(c(draws$y_sim[1,,3,,data$i_a]) ~ c(data$y[,3,,data$i_a]))
    # abline(0, 1)
    # plot(c(draws$y_sim[1,,3,,data$i_b]) ~ c(data$y[,3,,data$i_b]))
    # abline(0, 1)
        
    par(mfrow = c(1, 1))
        
  }
  
  if(modelname %in% c("ba", "ba-rag", "ba-rag-ranef")) {
    R <- data.frame(y_hat_rep = draws$y_hat_rep[1,], y = data$y,
                    pops = data$pops[rep(data$rep_init2y0, each = data$n_reobs[1])])
    
    R %>% ggformula::gf_point(y_hat_rep ~ y) %>%
      gf_abline(gformula = NULL, slope = 1, intercept = 0)
    
    # plot(draws$state_init[1, data$rep_init2y0] ~ data$y0)
    
  }
  
}


#### Draws vs true ----------------------
plotEstimateVsTrue <- function(parname = "h",
                          simdata,
                          stanfit) {
  
  simpar <- attr(simdata, "pars")[[parname]]
  
  if(is.vector(simpar)) {
    
    Stanpar <- rstan::extract(stanfit, pars = parname)[[1]] %>%
      as.data.frame() %>%
      pivot_longer(cols = everything(), names_to = "species", values_to = "draw") %>%
      mutate(species = as.integer(as.factor(species))) %>%
      bind_cols(true = rep(simpar, length.out = nrow(.)))
    
    Stanpar %>%
      # sample_n(1000) %>%
      gf_point(draw ~ true, alpha = 0.3, size = 0.5, title = parname) %>%
      gf_abline(slope = 1, intercept = 0, gformula = NULL)
    
  } else if (str_starts(parname, "Beta")) {
    
    print("Not yet implemented.")
    
  } else {
    
    # Stanpar <- rstan::extract(rstandraws, pars = parname)[[1]] %>%
    #   as.data.frame() %>%
    #   pivot_longer(cols = everything(), names_to = "species", values_to = "draw") %>%
    #   mutate(species = as.integer(as.factor(species))) %>%
    #   bind_cols(true = rep(simpar, length.out = nrow(.)))
    # 
    # Stanpar %>%
    #   # sample_n(1000) %>%
    #   gf_point(draw ~ true | species, alpha = 0.3, size = 0.1) %>%
    #   gf_abline(slope = 1, intercept = 0, gformula = NULL)
    
  }
}



#### Parameter vs Env ----------------------
predictParameterInEnv <- function(B, x = seq(-1, 1, by = 0.01), n_beta = 2, species = 1) {
  E <- as.data.frame(replicate(n_beta, x))
  polyformula <- as.formula(paste("~", paste("poly(", colnames(E), ", 2)", collapse = "+")))
  X <- model.matrix(polyformula, data = E)
  X %*% B[,1]
}

# poly2 <- function(x, b) {b[1] + x*b[2] + x^2*b[3]}
# getVertex <- function(b) {-b[2]/(2*b[3])}

plotParameterInEnv <- function(betaname, x = Env[,1]) {
  truepars <- attr(data, "pars")
  Beta_true <- truepars[[betaname]]
  Beta_mean <- matrix(fit$summary(variables = betaname)$mean, data$N_beta, data$N_totalspecies)
  
  X <- cbind(true = predictParameterInEnv(Beta_true, x), meandraw = predictParameterInEnv(Beta_mean, x))
  matplot(x, X, col = c("blue", "black"), pch = c("T", "D"), type = "p")
}


# plotParameterInEnv("Beta_g")

#### Draws vs true ----------------------
plotStatePairs <- function(standata) {
  
  ba_a_avg <- pi * ((standata$dbh_lower_a + standata$dbh_lower_b)/2/2)^2 * 1e-6
  
  D <- standata$Long %>%
    mutate(id = interaction(time, plot, loc)) %>%
    select(stage, species, abundance, id, time) %>%
    pivot_wider(id_cols = c("id", "time"), names_from = c("stage", "species"), values_from = c("abundance"))
  
  plotlevel <- as.integer(as.factor(D$time))
  
  D %>% 
    mutate(ba_sum = (a_1*ba_a_avg + b_1) + (a_2*ba_a_avg + b_2)) %>%
    select(-c("id", "time")) %>%
    pairs(col = plotlevel, pch = plotlevel)
  
}




