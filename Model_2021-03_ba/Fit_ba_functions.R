

####  Do the fit ----------------------
drawSamples <- function(model, data, method = c("variational", "mcmc", "sim"), n_chains = 3, initfunc = 0) {
  
  if(match.arg(method) == "variational") {
    fit <- model$variational(data = data,
                             output_dir = "Fits.nosync",
                             init = initfunc,
                             eta = 0.001,
                             iter = 20**4) # convergence after 1500 iterations
    
  } else if (match.arg(method) == "mcmc") {
    
    fit <- model$sample(data = data,
                        output_dir = "Fits.nosync",
                        init = initfunc,
                        iter_warmup = 400, iter_sampling = 600,
                        adapt_delta = 0.99,
                        max_treedepth = 16,
                        chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))
    
  } else if (match.arg(method) == "sim") {
    
    fit <- model$sample(data = data,
                        output_dir = NULL,
                        init = initfunc,
                        iter_sampling = 500,
                        fixed_param = TRUE)
  }
  
  return(fit)
}
