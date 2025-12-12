library(loo)

run.all.models.1 <- function(model.list, jags.data, params.to.monitor, MCMC.params, Rhat.params){
  
  loo.out <- list()
  Rmax <- list() #vector(mode = "numeric", length = length(out.names))
  
  
  for (k in 1:length(model.list)){
    file.name.root <- unlist(str_split(unlist(str_split(model.list[[k]]$file.name, ".txt"))[1],
                                       pattern = "Model_JAGS_"))[2]
    
    MCMC.params$model.file = paste0("models/", model.list[[k]]$file.name)
    jags.data$X <- model.list[[k]]$Cov
    
    if (!file.exists(paste0("RData/", model.list[[k]]$out.file.name))){
      
      tic <- Sys.time()
      jm <- jags(data = jags.data,
                 #inits = inits,
                 parameters.to.save= params.to.monitor,
                 model.file = MCMC.params$model.file,
                 n.chains = MCMC.params$n.chains,
                 n.burnin = MCMC.params$n.burnin,
                 n.thin = MCMC.params$n.thin,
                 n.iter = MCMC.params$n.samples,
                 DIC = T, 
                 parallel=T)
      
      toc <- Sys.time()
      
      out.list <- list(jags.out = jm,
                       jags.data = jags.data,
                       Run.Date = tic,
                       Run.Time = toc - tic,
                       MCMC.params = MCMC.params)
      
      saveRDS(out.list, 
              file = paste0("RData/", model.list.FL[[k]]$out.file.name))
      
    } else {
      out.list <- readRDS(file = paste0("RData/", model.list.FL[[k]]$out.file.name))
    }
    
    Rmax[[k]] <- rank.normalized.R.hat(out.list$jags.out$samples, 
                                       params = Rhat.params, 
                                       MCMC.params = MCMC.params)
    
    loglik.mat <- out.list$jags.out$sims.list$loglik
    
    n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
    Reff <- relative_eff(exp(loglik.mat),
                         chain_id = rep(1:MCMC.params$n.chains,
                                        each = n.per.chain),
                         cores = MCMC.params$n.chains)
    
    loo.out <- rstanarm::loo(loglik.mat, 
                             r_eff = Reff, 
                             cores = MCMC.params$n.chains, 
                             k_threshold = 0.7)
    
    loo.out[[k]] <- list(Reff = Reff,
                         loo.out = loo.out)
    
  }
  
  return(list(Rmax = Rmax,
              loo.out = loo.out))
}

compute.LOOIC <- function(loglik.array, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  n.samples <- dim(loglik.array)[1]
  
  # Create an empty list to hold the log likelihood values
  #ll.list <- vector("list", length = n.samples)
  
  loglik.list <- lapply(1:n.samples, function(s) {
    # Extract the slice for the current sample 's'
    slice <- loglik.array[s, , , ]
    # Return the non-NA values from that slice
    slice[!is.na(slice)]
  })
  
  #loglik.list <- apply(loglik.array, 1, function(slice) slice[!is.na(slice)])
  
  # loop through each MCMC sample and collect log likelihood values
  # for (s in 1:n.samples){
  #   ll.cube <- loglik.array[s,,,]
  #   ll.list[[s]] <- ll.cube[!is.na(loglik.array[s,,,])]
  # }
  
  loglik.mat <- do.call(rbind, loglik.list)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = MCMC.params$n.chains)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = MCMC.params$n.chains, 
                           k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}


# Extracting posterior samples of deviance or any other variable from jags output:
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}


# Compute the "rank-normalized R-hat" by Vehtari et al. (2021) from jagsUI
# output.
# Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. Bayesian Analysis, 16(2), 667–718.
# https://doi.org/10.1214/20-BA1221
# 
# The first input is MCMC samples. If the jagsUI output is jm, this is jm$samples.
# The second input is a string of regular expression. This is a bit
# complicated. For example, to select "BF.Fixed" and all K parameters, which 
# are indexed, Use "^BF\\.Fixed|^K\\[" A '^' specifies that the following letter
# is the beginning of a string. '\\.' specifies a literal period, which needs
# to be "escaped" by two backslashes (\\). A square bracket needs to be escaped
# with two backslashes as well. The pipe (|) indicates 'or'.  
# 

rank.normalized.R.hat <- function(samples, params, MCMC.params){
  library(posterior)
  library(coda)
  
  col.names <- grep(params, varnames(samples), value = TRUE, perl = TRUE)
  subset.mcmc.samples <- samples[, col.names]
  
  subset.mcmc.array <- as_draws_array(subset.mcmc.samples, .nchains = MCMC.params$n.chains)
  
  rhat.values <- apply(subset.mcmc.array, MARGIN = 3, FUN = rhat)
  
  return(rhat.values)
}
