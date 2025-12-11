# Analysis of leatherback turtle nesting data for Florida,
# Puerto Rico, and St. Croix. A collaboration with Kelly
# Stewart. Started several years ago but it went down in the
# priority list for both of us. It was picked up again at 
# the end of 2025. 
# 
# In this script, all models are run and their results are 
# compared using LOOIC. Convergence of MCMC is determined
# using the rank-normalized Rhat statistics. Results from 
# this script should be summarized elsewhere.
# 

rm(list=ls())
library(tidyverse)
library(ggplot2)
library(readr)
library(lubridate)
library(jagsUI)
library(bayesplot)
library(stringr)

source("Caribbean_Dc_fcns.R")

# MCMC setup
MCMC.params <- list(n.samples = 50000,
                    n.burnin = 30000,
                    n.thin = 5,
                    n.chains = 5)

n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin

# Get data 
col.def <- cols(ID = col_integer(),
                year = col_integer(),
                beach = col_character(),
                latitude = col_double(),
                distance = col_double(),
                days_week = col_integer(),
                days_year = col_integer(),
                nests = col_integer())

FL.nest.counts <- read_csv("data/FL_Sept2020.csv", 
                           col_types = col.def) %>% 
  mutate(beach_f = as.factor(toupper(beach)),
         dataset = "FL",
         ID2 = as.numeric(as.factor(ID)),
         lat_band = ifelse(latitude < 26, 25,
                           ifelse(latitude < 27, 26,
                                  ifelse(latitude < 28, 27,
                                         ifelse(latitude < 29, 28,
                                                ifelse(latitude < 30, 29,
                                                       ifelse(latitude < 31, 30)))))),
         beach_f = as.factor(toupper(beach)))


# parameters to monitor
# These are for the original model by Michelle Sims
# Each area was ran separately with area-specific covariates
parameters.to.monitor.FL <- c("a0", "a1", "beta", "B.hat", "a1pc",
                             "mu.a0", "mu.a1",
                             "rea0", "rea1",
                             "delta1", "delta2",
                             "sigma.e", "Sigma.B",
                             "sigma.a0", "sigma.a1",
                             "floridapc", "probP", "deviance",
                             "Devobs", "Devpred", "loglik")

# DAta for FL


# Define covariates for FL
Cov.FL <- c("distance", "days_week", "days_year")
X.FL <- list(cbind(log(FL.nest.counts$distance) -
                     median(log(FL.nest.counts$distance)),
                   FL.nest.counts$days_week -
                     median(FL.nest.counts$days_week),
                   FL.nest.counts$days_year -
                     median(FL.nest.counts$days_year)),
             cbind(log(FL.nest.counts$distance) -
                     median(log(FL.nest.counts$distance)),
                   FL.nest.counts$days_week -
                     median(FL.nest.counts$days_week)),
             cbind(log(FL.nest.counts$distance) -
                     median(log(FL.nest.counts$distance)),
                   FL.nest.counts$days_year -
                     median(FL.nest.counts$days_year)),
             cbind(FL.nest.counts$days_week -
                     median(FL.nest.counts$days_week),
                   FL.nest.counts$days_year -
                     median(FL.nest.counts$days_year)),
             log(FL.nest.counts$distance) -
               median(log(FL.nest.counts$distance)),
             FL.nest.counts$days_week - 
               median(FL.nest.counts$days_week),
             FL.nest.counts$days_year - 
               median(FL.nest.counts$days_year),
             0)

model.names.FL <- c("3Covs", rep("2Covs", times = 3), 
                 rep("1Cov", times = 3), "0Cov")

out.names.FL <- c("3Covs", "logD_DayWk", "logD_DayYr", "DayWk_DayYr",
               "logD", "DayWk", "DayYr", "0Cov") 

model.list.FL <- list()
c <- 1
for (k in 1:length(out.names.FL)){
  model.list.FL[[c]] <- c(ID = c, 
                         file.name = paste0("Model_JAGS_rSlope_rInt_",
                                            model.names.FL[k], ".txt"))
  c <- c + 1
  model.list.FL[[c]] <- c(ID = c, 
                         file.name = paste0("Model_JAGS_negbin_rSlope_rInt_",
                                            model.names.FL[k], ".txt"))
  c <- c + 1  
}

# Run all models for FL
loo.out.FL <- list()
Rmax.FL <- list() #vector(mode = "numeric", length = length(out.names))

for (k in 1:length(out.names.FL)){
  MCMC.params$model.file = paste0("models/Model_JAGS_rSlope_rInt_",
                                  model.names.FL[k], ".txt")
  FL.jags.data$X <- X.FL[[k]]
  
  if (!file.exists(paste0("RData/JAGS_out_rSlope_rInt_", 
                          out.names[k], "_FL.rds"))){
    
    tic <- Sys.time()
    jm <- jags(data = FL.jags.data,
               #inits = inits,
               parameters.to.save= parameters.to.monitor.FL,
               model.file = MCMC.params$model.file,
               n.chains = MCMC.params$n.chains,
               n.burnin = MCMC.params$n.burnin,
               n.thin = MCMC.params$n.thin,
               n.iter = MCMC.params$n.samples,
               DIC = T, 
               parallel=T)
    
    toc <- Sys.time()
    
    out.list <- list(jags.out = jm,
                     jags.data = FL.jags.data,
                     Run.Date = tic,
                     Run.Time = toc - tic,
                     MCMC.params = MCMC.params)
                     
    saveRDS(out.list, 
            file = paste0("RData/JAGS_out_rSlope_rInt_", 
                          out.names[k], "_FL.rds"))
    
  } else {
    out.list <- readRDS(file = paste0("RData/JAGS_out_rSlope_rInt_",
                                out.names[k], "_FL.rds"))
  }
  
  Rmax.FL[[k]] <- rank.normalized.R.hat(out.list$jags.out$samples, 
                                        params = params.to.monitor.1, 
                                        MCMC.params = MCMC.params)
  
  loo.out.FL[[k]] <- compute.LOOIC(loglik.array = jm$sims.list$loglik, 
                                   MCMC.params = MCMC.params)
  
}


# Covarates for STX
Cov.STX <- c("days_year")
X.STX <- list(STX.nest.counts$days_year -
                median(STX.nest.counts$days_year),
              0)

model.names.STX <- c("1Cov", "0Cov")
out.names.STX <- c("DayYr", "0Cov") 

model.list.STX <- list()
c <- 1
for (k in 1:length(out.names)){
  model.list.STX[[c]] <- c(ID = c, 
                          file.name = paste0("Model_JAGS_rSlope_rInt_",
                                             model.names.STX[k], ".txt"))
  c <- c + 1
  model.list.STX[[c]] <- c(ID = c, 
                          file.name = paste0("Model_JAGS_negbin_rSlope_rInt_",
                                             model.names.STX[k], ".txt"))
  c <- c + 1  
}


nest.counts <- read_csv("data/FL_PR_STX.csv", 
                        col_types = col.def) %>% 
  mutate(lat_band = ifelse(latitude < 26, 25,
                           ifelse(latitude < 27, 26,
                                  ifelse(latitude < 28, 27,
                                         ifelse(latitude < 29, 28,
                                                ifelse(latitude < 30, 29,
                                                       ifelse(latitude < 31, 30)))))),
         beach_f = as.factor(toupper(beach)))

# There are some skipped years, which needs to be filled in
nest.counts %>% 
  group_by(ID) %>%
  summarise(n.years = max(year) - min(year) + 1,
            year.1 = min(year),
            year.2 = max(year),
            lat = first(latitude)) -> summary.years

year.mat <- y <- matrix(nrow = nrow(summary.years), 
                        ncol = max(summary.years$n.years))

# find which data points are missing.
NA.idx <- vector(mode = "list", length = nrow(summary.years))

k <- 3
for (k in 1:nrow(summary.years)){
  nest.counts %>% 
    filter(ID == summary.years$ID[k]) %>%
    select(year, nests) %>% 
    mutate(seq.yr = year - min(year) + 1) -> tmp
  
  y[k,tmp$seq.yr] <- tmp$nests
  year.mat[k, tmp$seq.yr] <- tmp$year
  n.years.k <- summary.years$n.years[k]
  if (sum(is.na(year.mat[k, 1:n.years.k])) > 0)
    NA.idx[[k]] <- c(1:n.years.k)[is.na(year.mat[k, 1:n.years.k])]
}

# Create a vector of variable names for missing years
missing.y <- c()

for (k in 1:length(NA.idx)){
  if (length(NA.idx[[k]]) > 0){
    for (k1 in 1:length(NA.idx[[k]])){
      missing.y <- c(missing.y, paste0("y[", k, ",", NA.idx[[k]][k1], "]"))
    }
  }
}


# These are for my new models:
parameters.to.monitor.2 <- c(c("U", "mean.U1", "sigma.U1",
                               "mean.U2", "sigma.U2", "p",
                               "sigma.Q",  "loglik", "N"),
                             missing.y)

# Model files need to be in a list
models.list.2 <- list(ID = c(1:5),
                    file.names = c("Model_norm_norm_QsRs_2Us.txt",
                                   "Model_norm_Pois_Qs_2Us.txt",
                                   "Model_norm_Pois_Qs_2Us_lat.txt",
                                   "Model_norm_Pois_Qs_2Us_skip.txt",
                                   "Model_norm_Pois_2Us_skip_Nsum.txt"))





