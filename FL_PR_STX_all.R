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

# Define data columns 
col.def <- cols(ID = col_integer(),
                year = col_integer(),
                beach = col_character(),
                latitude = col_double(),
                distance = col_double(),
                days_week = col_integer(),
                days_year = col_integer(),
                nests = col_integer())

## Running models by Michelle Sims. Using each beach separately ##

# Parameters to monitor for convergence:
Rhat.params <- "^a0\\[|^a1\\[|^beta\\[|^mu.a0|^m.a1|^delta1|^delta2|^sigma."

################################ FL #######################################

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

FL.nest.counts %>% 
  group_by(ID2) %>%
  summarise(n = n()) -> FL.ns

FL.nest.counts %>% 
  select(ID2, ID, latitude, beach_f, year) %>% 
  group_by(ID) %>%
  summarise(ID2 = first(ID2),
            name = first(beach_f),
            latitude = first(latitude),
            median.yr = median(year),
            min.yr = min(year),
            n = n()) %>%
  mutate(latc = latitude - mean(latitude),
         latc2 = latc^2,
         beach = ID,
         ID2 = ID2,
         name = name) -> FL.lat.dat

median.yr.FL <- select(FL.lat.dat, c(ID, median.yr, min.yr)) %>%
  right_join(FL.nest.counts, by = "ID") %>% 
  select(ID, median.yr, min.yr, year)

# parameters to monitor
# These are for the original model by Michelle Sims
# Each area was ran separately with area-specific covariates
parameters.to.monitor.FL <- c("a0", "a1", "beta", "B.hat", "a1pc",
                             "mu.a0", "m.a1",
                             "rea0", "rea1",
                             "delta1", "delta2",
                             "sigma.e", 
                             "sigma.a0", "sigma.a1", "rho",
                             "loglik")

# Data for FL
jags.data.FL <- list(N = length(FL.nest.counts$nests),
                     nbeach = length(unique(FL.nest.counts$ID2)),
                     count = FL.nest.counts$nests,
                     beach = FL.nest.counts$ID2,
                     yearc = median.yr.FL$year - median.yr.FL$min.yr,
                     latc = FL.lat.dat$latc,
                     latc2 = FL.lat.dat$latc2)


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

model.names.FL <- c("3Covs", 
                    rep("2Covs", times = 3), 
                    rep("1Cov", times = 3), 
                    "0Cov")

out.names.FL <- c("3Covs", "logD_DayWk", "logD_DayYr", "DayWk_DayYr",
                  "logD", "DayWk", "DayYr", "0Cov") 

model.list.FL <- list()
c <- 1
for (k in 1:length(out.names.FL)){
  model.list.FL[[c]] <- list(ID = c, 
                             file.name = paste0("Model_JAGS_Pois_rSlope_rInt_",
                                                model.names.FL[k], ".txt"),
                             Cov = X.FL[[k]],
                             out.file.name = paste0("JAGS_out_Pois_rSlope_rInt_",
                                                    out.names.FL[k], "_FL.rds"))
  c <- c + 1
  # model.list.FL[[c]] <- list(ID = c, 
  #                            file.name = paste0("Model_JAGS_negbin_rSlope_rInt_",
  #                                               model.names.FL[k], ".txt"),
  #                            Cov = X.FL[[k]],
  #                            out.file.name = paste0("JAGS_out_negbin_rSlope_rInt_",
  #                                                   out.names.FL[k], "_FL.rds"))
  # c <- c + 1  
}

Rhat.params.FL <- "^a0\\[|^a1\\[|^beta\\[|^mu.a0|^m.a1|^delta1|^delta2|^sigma."

# Run all models for FL
out.FL <- run.all.models.1(model.list = model.list.FL,
                           jags.data = jags.data.FL,
                           params.to.monitor = parameters.to.monitor.FL,
                           MCMC.params = MCMC.params,
                           Rhat.params = Rhat.params.FL)

###################################### FL-2 #################################

# Also try running two-slope models:
model.list.FL <- list()
c <- 1
for (k in 1:length(out.names.FL)){
  model.list.FL[[c]] <- list(ID = c, 
                             file.name = paste0("Model_JAGS_Pois_r2Slopes_rInt_",
                                                model.names.FL[k], ".txt"),
                             Cov = X.FL[[k]],
                             out.file.name = paste0("JAGS_out_Pois_r2Slopes_rInt_",
                                                    out.names.FL[k], "_FL.rds"))
  c <- c + 1
  # model.list.FL[[c]] <- list(ID = c, 
  #                            file.name = paste0("Model_JAGS_negbin_rSlope_rInt_",
  #                                               model.names.FL[k], ".txt"),
  #                            Cov = X.FL[[k]],
  #                            out.file.name = paste0("JAGS_out_negbin_rSlope_rInt_",
  #                                                   out.names.FL[k], "_FL.rds"))
  # c <- c + 1  
}

parameters.to.monitor.FL <- c("a0.1", "a0.2",
                              "a1.1", "a1.2",
                              "beta", "B.hat", 
                              "mu.a0.1", "m.a1.1",
                              "mu.a0.2", "m.a1.2",
                              "rho.1", "rho.2",
                              "delta1", "delta2",
                              "sigma.e", 
                              "sigma.a0.1", "sigma.a1.1",
                              "sigma.a0.2", "sigma.a1.2",
                              "year.change",
                              "loglik")

Rhat.params.FL <- "^a0.1\\[|^a0.2\\[|^a1.1\\[^a1.2\\[|^beta\\[|^mu.a0|^m.a1|^delta1|^delta2|^sigma.|^rho.|year.change"

# Data for FL
jags.data.FL <- list(N = length(FL.nest.counts$nests),
                     nbeach = length(unique(FL.nest.counts$ID2)),
                     count = FL.nest.counts$nests,
                     beach = FL.nest.counts$ID2,
                     yearc = median.yr.FL$year - median.yr.FL$min.yr,
                     latc = FL.lat.dat$latc,
                     latc2 = FL.lat.dat$latc2,
                     minT = 23, maxT = 33)

# Run all models for FL
out.FL.2 <- run.all.models.1(model.list = model.list.FL,
                           jags.data = jags.data.FL,
                           params.to.monitor = parameters.to.monitor.FL,
                           MCMC.params = MCMC.params,
                           Rhat.params = Rhat.params.FL)

################################ STX #######################################
# Data - only one beach
STX.nest.counts <- read_csv("data/STX_Sept2020.csv", 
                            col_types = col.def) %>% 
  mutate(beach_f = as.factor(toupper(beach)),
         dataset = "STX",
         ID2 = as.numeric(as.factor(ID)))

STX.nest.counts %>% 
  group_by(ID2) %>%
  summarise(n = n()) -> STX.ns

STX.nest.counts %>% 
  select(ID2, ID, beach_f, year) %>% 
  group_by(ID) %>%
  summarise(ID2 = first(ID2),
            name = first(beach_f),
            median.yr = median(year),
            min.yr = min(year),
            n = n()) %>%
  mutate(beach = ID,
         ID2 = ID2,
         name = name) -> STX.dat

STX.dat %>%
  select(c(ID, median.yr, min.yr)) %>%
  right_join(STX.nest.counts, by = "ID") %>% 
  select(ID, median.yr, min.yr, year) -> median.yr.STX 

jags.data.STX <- list(N = length(STX.nest.counts$nests),
                      nbeach = length(unique(STX.nest.counts$ID2)),
                      count = STX.nest.counts$nests,
                      beach = STX.nest.counts$ID2,
                      yearc = median.yr.STX$year - median.yr.STX$min.yr)

# Covarates for STX
Cov.STX <- c("days_year")
X.STX <- list(STX.nest.counts$days_year -
                median(STX.nest.counts$days_year),
              0)

model.names.STX <- c("1Cov_NoLat", "0Cov_NoLat")
out.names.STX <- c("DayYr", "0Cov") 

model.list.STX <- list()
c <- 1
for (k in 1:length(out.names.STX)){
  model.list.STX[[c]] <- list(ID = c, 
                              file.name = paste0("Model_JAGS_Pois_",
                                                 model.names.STX[k], ".txt"),
                              Cov = X.STX[[k]],
                              out.file.name = paste0("JAGS_out_Pois_",
                                                     out.names.STX[k], "_STX.rds"))
  c <- c + 1
  # model.list.STX[[c]] <- list(ID = c, 
  #                             file.name = paste0("Model_JAGS_negbin_rSlope_rInt_",
  #                                            model.names.STX[k], ".txt"),
  #                             Cov = X.STX[[k]])
  # c <- c + 1  
}

parameters.to.monitor.STX <- c("a0", "a1", "a2", "beta", 
                               "sigma.e",
                               "epsilon",
                               "loglik")

Rhat.params.STX <- "^a0|^a1|^a2|^epsilon|^sigma."

# Run all models for STX
out.STX <- run.all.models.1(model.list = model.list.STX,
                            jags.data = jags.data.STX,
                            params.to.monitor = parameters.to.monitor.STX,
                            MCMC.params = MCMC.params,
                            Rhat.params = Rhat.params.STX)

###################################### STX-2 #################################

# Also try running two-slope models:
model.list.STX <- list()
c <- 1
for (k in 1:length(out.names.STX)){
  model.list.STX[[c]] <- list(ID = c, 
                             file.name = paste0("Model_JAGS_Pois_2Slopes_",
                                                model.names.STX[k], ".txt"),
                             Cov = X.STX[[k]],
                             out.file.name = paste0("JAGS_out_Pois_2Slopes_",
                                                    out.names.STX[k], "_STX.rds"))
  c <- c + 1
  # model.list.FL[[c]] <- list(ID = c, 
  #                            file.name = paste0("Model_JAGS_negbin_rSlope_rInt_",
  #                                               model.names.FL[k], ".txt"),
  #                            Cov = X.FL[[k]],
  #                            out.file.name = paste0("JAGS_out_negbin_rSlope_rInt_",
  #                                                   out.names.FL[k], "_FL.rds"))
  # c <- c + 1  
}

parameters.to.monitor.STX <- c("a0.1", "a0.2",
                               "a1.1", "a1.2",
                               "beta",  
                               "sigma.e", 
                               "year.change",
                               "loglik")

Rhat.params.STX <- "^a0.1|^a0.2|^a1.1^a1.2|^beta|^sigma.|year.change"

# Data for STX
jags.data.STX <- list(N = length(STX.nest.counts$nests),
                     nbeach = length(unique(STX.nest.counts$ID2)),
                     count = STX.nest.counts$nests,
                     beach = STX.nest.counts$ID2,
                     yearc = median.yr.STX$year - median.yr.STX$min.yr,
                     minT = 23, maxT = 33)

# Run all models for STX
out.STX.2 <- run.all.models.1(model.list = model.list.STX,
                             jags.data = jags.data.STX,
                             params.to.monitor = parameters.to.monitor.STX,
                             MCMC.params = MCMC.params,
                             Rhat.params = Rhat.params.STX)


################################ PR #######################################
# For PR, no covariates were used for now
PR.nest.counts <- read_csv("data/PR_Sept2020.csv", 
                            col_types = col.def) %>% 
  mutate(beach_f = as.factor(toupper(beach)),
         dataset = "PR",
         ID2 = as.numeric(as.factor(ID)))

PR.nest.counts %>% 
  group_by(ID2) %>%
  summarise(n = n()) -> PR.ns

PR.nest.counts %>% 
  select(ID2, ID, beach_f, year) %>% 
  group_by(ID) %>%
  summarise(ID2 = first(ID2),
            name = first(beach_f),
            median.yr = median(year),
            min.yr = min(year),
            n = n()) %>%
  mutate(beach = ID,
         ID2 = ID2,
         name = name) -> PR.lat.dat

PR.lat.dat %>%
  select(c(ID, median.yr, min.yr)) %>%
  right_join(PR.nest.counts, by = "ID") %>% 
  select(ID, median.yr, min.yr, year) -> median.yr.PR 

jags.data.PR <- list(N = length(PR.nest.counts$nests),
                     nbeach = length(unique(PR.nest.counts$ID2)),
                     count = PR.nest.counts$nests,
                     beach = PR.nest.counts$ID2,
                     yearc = median.yr.PR$year - median.yr.PR$min.yr)

Cov.PR <- c("days_year")
X.PR <- list(PR.nest.counts$days_year -
               median(PR.nest.counts$days_year),
             0)

model.names.PR <- c("1Cov_NoLat", "0Cov_NoLat")
out.names.PR <- c("DayYr", "0Cov") 

model.list.PR <- list()
c <- 1
for (k in 1:length(out.names.PR)){
  model.list.PR[[c]] <- list(ID = c, 
                             file.name = paste0("Model_JAGS_Pois_",
                                                 model.names.PR[k], ".txt"),
                              Cov = X.PR[[k]],
                              out.file.name = paste0("JAGS_out_Pois_",
                                                     out.names.PR[k], "_PR.rds"))
  c <- c + 1
}

parameters.to.monitor.PR <- c("a0", "a1", "a2", "beta", 
                              "sigma.e",
                              "epsilon",
                              "loglik")

Rhat.params.PR <- "^a0|^a1|^a2|^epsilon|^sigma."

# Run all models for PR
out.PR <- run.all.models.1(model.list = model.list.PR,
                            jags.data = jags.data.PR,
                            params.to.monitor = parameters.to.monitor.PR,
                            MCMC.params = MCMC.params,
                            Rhat.params = Rhat.params.PR)

###################################### PR-2 #################################

# Also try running two-slope models:
model.list.PR <- list()
c <- 1
for (k in 1:length(out.names.PR)){
  model.list.PR[[c]] <- list(ID = c, 
                              file.name = paste0("Model_JAGS_Pois_2Slopes_",
                                                 model.names.PR[k], ".txt"),
                              Cov = X.PR[[k]],
                              out.file.name = paste0("JAGS_out_Pois_2Slopes_",
                                                     out.names.PR[k], "_PR.rds"))
  c <- c + 1
}

parameters.to.monitor.PR <- c("a0.1", "a0.2",
                               "a1.1", "a1.2",
                               "beta",  
                               "sigma.e", 
                               "year.change",
                               "loglik")

Rhat.params.PR <- "^a0.1|^a0.2|^a1.1^a1.2|^beta|^sigma.|year.change"

# Data for PR
jags.data.PR <- list(N = length(PR.nest.counts$nests),
                      nbeach = length(unique(PR.nest.counts$ID2)),
                      count = PR.nest.counts$nests,
                      beach = PR.nest.counts$ID2,
                      yearc = median.yr.PR$year - median.yr.PR$min.yr,
                      minT = 23, maxT = 33)

# Run all models for PR
out.PR.2 <- run.all.models.1(model.list = model.list.PR,
                              jags.data = jags.data.PR,
                              params.to.monitor = parameters.to.monitor.PR,
                              MCMC.params = MCMC.params,
                              Rhat.params = Rhat.params.PR)




############## A new approach starts here ###################################

# Create a vector of variable names for missing years
missing.y <- c()

for (k in 1:length(NA.idx)){
  if (length(NA.idx[[k]]) > 0){
    for (k1 in 1:length(NA.idx[[k]])){
      missing.y <- c(missing.y, 
                     paste0("y[", k, ",", NA.idx[[k]][k1], "]"))
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





