library(loo)


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
