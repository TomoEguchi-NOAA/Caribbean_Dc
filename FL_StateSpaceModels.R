# Analysis of leatherback turtle nesting data for Florida,
# using the new state-space models for debugging purposes.
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
                                                       ifelse(latitude < 31, 30)))))))

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

############## A new approach starts here ###################################

# for the Norm-Norm model, observed y+1 is in the log space.
# for the Norm-Pois models, observed y is a count
# "Model_norm_norm_QsRs_2Us.txt",

model.file.names <- c("Model_norm_Pois_Qs_2Us.txt",
                      "Model_norm_Pois_Qs_U_skip.txt",
                      "Model_norm_Pois_Qs_2Us_skip.txt",
                      "Model_norm_Pois_Qs_U.txt")
#                      "Model_norm_Pois_2Us_skip_Nsum.txt")

# Model files in a list for creating a summary list later
models.table.2 <- data.frame(ID = c(1:length(model.file.names)),
                             file.names = model.file.names)

######################## FL ######################################
# Data for FL0
FL.nest.counts %>% 
  group_by(ID) %>%
  summarise(n.years = max(year) - min(year) + 1,
            year.1 = min(year),
            year.2 = max(year),
            lat = first(latitude)) -> FL.summary.years

year.1 <- min(FL.nest.counts$year)

# Convert y, beach ID, and year into vectors
FL.year.vec <- as.vector(FL.nest.counts$year)
FL.y.vec <- as.vector(FL.nest.counts$nests)

FL.ID.vec <- as.vector(FL.nest.counts$ID) #$ rep(as.vector(FL.summary.years$ID),
                 #times = as.vector(FL.summary.years$n.years))

FL.lat.vec <- as.vector(FL.summary.years$lat)

# FL.year.idx.vec <- vector(mode = "numeric", length = length(FL.y.vec))
# FL.year.idx.vec[FL.year.vec < 2012] <- 1
# FL.year.idx.vec[FL.year.vec > 2011] <- 2

c <- 1
k <- 3
uniq.IDs <- unique(FL.summary.years$ID)
min.seq.yr <- high.low.1 <- vector(mode = "numeric", length = length(uniq.IDs))
FL.y.1 <- FL.ID.1 <- FL.years.1 <- vector(mode = "numeric", length = length(uniq.IDs))

FL.ID.vec.all <- FL.y.vec.all <- vector(mode = "numeric")
FL.years <- FL.year.seq.vec <- vector(mode = "numeric")
FL.year.1toT.vec <- vector(mode = "numeric")

k <- 3
for (k in 1:length(uniq.IDs)){
  FL.summary.years %>%
    filter(ID == uniq.IDs[k]) -> summary.1
  
  FL.nest.counts %>%
    filter(ID == uniq.IDs[k]) -> counts.1

  # Create a dataframe with sequential years, including years without data
  year.df <- data.frame(year = seq(min(counts.1$year), 
                                   max(counts.1$year)))
  
  # Create a dataframe with NAs in nest counts for unsampled years
  counts.1 %>%
    select(year, nests) %>% 
    right_join(y = year.df, by = "year") %>%
    arrange(by = year) %>%
    mutate(seq.yr = year - year.1 + 1,
           seq.1toT = seq(from = 1, 
                          to = (max(year) - min(year) + 1))) -> tmp

  FL.ID.vec.all[c:(c+summary.1$n.years-1)] <- uniq.IDs[k]
  FL.y.vec.all[c:(c+summary.1$n.years-1)] <- tmp$nests
  
  # convert sampling years to sequential years from the first year (year.1)  
  FL.year.seq.vec[c:(c+summary.1$n.years-1)] <- tmp$seq.yr
  FL.year.1toT.vec[c:(c+summary.1$n.years-1)] <- tmp$seq.1toT
  
  FL.years[c:(c+summary.1$n.years-1)] <- tmp$year
  
  min.seq.yr[k] <- min(tmp$seq.yr, na.rm = T)
  high.low.1[k] <- ifelse(min.seq.yr[k] %% 2 == 0, 1, 2)
  FL.y.1[k] <- counts.1$nests[1]
  FL.ID.1[k] <- uniq.IDs[k]
  
  FL.years.1[k] <- ifelse(tmp$year[1] < 2012, 1, 2)
  
  c <- c + nrow(tmp)
}

FL.year.1toT.1 <- FL.year.1toT.vec[FL.year.1toT.vec == 1]
FL.year.1toT.2plus <- FL.year.1toT.vec[FL.year.1toT.vec > 1]
FL.y.2plus <- FL.y.vec.all[FL.year.1toT.vec > 1]
FL.year.2plus <- FL.year.seq.vec[FL.year.1toT.vec > 1]
FL.ID.2plus <- FL.ID.vec.all[FL.year.1toT.vec > 1]
FL.years.2plus <- ifelse(FL.years[FL.year.1toT.vec > 1] < 2012, 1, 2)
high.low.2plus <- ifelse(FL.years[FL.year.1toT.vec > 1] %% 2 == 0, 1, 2)

# Need to find missing y index. 
FL.missing.y <- c()
NA.idx <- c(1:length(FL.y.2plus))[is.na(FL.y.2plus)]

for (k in 1:length(NA.idx)){
  FL.missing.y <- c(FL.missing.y, paste0("y.2[", NA.idx[k], "]"))
}

#
# These are for my new models:
parameters.to.monitor.2 <- c("U1", "U2",
                             "mean.U1", "sigma.U1",
                             "mean.U2", "sigma.U2", 
                             "p",
                             "sigma.Q",
                             "r",  "loglik", "N", "N0",
                             "b0.U1", "b1.U1",
                             "b0.U2", "b1.U2",
                             FL.missing.y)

Rhat.params <- "^U|^mean.|^sigma.|^b|^N"

out.list <- list()
k <- 1
for (k in 1:length(model.file.names)){
  MCMC.params$model.file <- paste0("models/", model.file.names[k]) 
  tmp.1 <- str_split(model.file.names[k], "Model_")[[1]][2]
  out.file.name <- paste0("RData/JAGS_out_", str_split(tmp.1, ".txt")[[1]][1], "_FL.rds") 
  if (str_detect(model.file.names[k], "norm_norm")){
    # For the normal likelihood, I have to add 1 to counts so that there is no
    # log(0) = -Inf, which causes problems (errors)
    jags.data.FL.2 <- list(ID.1 = FL.ID.1,
                           ID.2 = FL.ID.2plus,
                           n.ID = length(FL.y.1),
                           years.1 = FL.years.1,
                           years.2 = FL.years.2plus,
                           #year.1 = ,
                           year.2 = FL.year.1toT.2plus,
                           y.1 = log(FL.y.1 + 1),
                           y.2 = log(FL.y.2plus + 1),
                           n.2plus = length(FL.y.2plus))
      
      # n.beaches = nrow(FL.y),
      #                      n.years = FL.summary.years$n.years,
      #                      y = log(FL.y + 1),
      #                      years = FL.years)
      # 
  } else if (str_detect(model.file.names[k], "norm_Pois")){
    jags.data.FL.2 <- list(ID.1 = FL.ID.1,
                           ID.2 = FL.ID.2plus,
                           n.ID = length(FL.y.1),
                           years.1 = FL.years.1,
                           years.2 = FL.years.2plus,
                           #year.1 = min.seq.yr,
                           year.2 = FL.year.1toT.2plus,
                           y.1 = FL.y.1,
                           y.2 = FL.y.2plus,
                           n.2plus = length(FL.y.2plus),
                           lat = FL.lat.vec,
                           high.low.1 = high.low.1,
                           high.low.2plus = high.low.2plus)
    
  }
  
  if (!file.exists(out.file.name)){
    tic <- Sys.time()
    jm <- jags(data = jags.data.FL.2,
               #inits = inits,
               parameters.to.save= parameters.to.monitor.2,
               model.file = MCMC.params$model.file,
               n.chains = MCMC.params$n.chains,
               n.burnin = MCMC.params$n.burnin,
               n.thin = MCMC.params$n.thin,
               n.iter = MCMC.params$n.samples,
               DIC = T, 
               parallel=T)
    toc <- Sys.time() - tic
    
    Rmax <- rank.normalized.R.hat(jm$samples, 
                                  params = Rhat.params, 
                                  MCMC.params = MCMC.params)
    
    loglik.mat <- jm$sims.list$loglik
    
    n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
    Reff <- relative_eff(exp(loglik.mat),
                         chain_id = rep(1:MCMC.params$n.chains,
                                        each = n.per.chain),
                         cores = MCMC.params$n.chains)
    
    loo.out <- rstanarm::loo(loglik.mat, 
                             r_eff = Reff, 
                             cores = MCMC.params$n.chains, 
                             k_threshold = 0.7)
    
    out.list[[k]] <- list(jags.out = jm,
                          jags.data = jags.data.FL.2,
                          Run.Date = tic,
                          Run.Time = toc,
                          MCMC.params = MCMC.params,
                          Rmax = Rmax,
                          loo.out = loo.out,
                          parameters = parameters.to.monitor.2)
    
    saveRDS(out.list[[k]], 
            file = out.file.name)
  } else {
    
    out.list[[k]] <- readRDS(out.file.name)
  
    
  }
}

# check convergence and compare models
rank.Rhat <- lapply(out.list, 
                    FUN = function(x) x$Rmax)

big.rank.Rhat <- lapply(rank.Rhat,
                        FUN = function(x) x[x > 1.01])

n.big.rank.Rhat <- lapply(big.rank.Rhat,
                          FUN = function(x) length(x))

max.big.rank.Rhat <- lapply(rank.Rhat,
                            FUN = function(x) max(x))

# From rank-normalized Rhats, all models seem okay.

# Check goodness-of-fit
looic <- lapply(out.list,
                FUN = function(x) x$loo.out)

looic.table.list <- lapply(looic, FUN = function(x){
  out <- data.frame(looic = x$estimates["looic", "Estimate"],
                    SE = x$estimates["looic", "SE"])
  return(out)
})

looic.table <- do.call("rbind", looic.table.list) %>%
  rownames_to_column(var = "ID") %>%
  arrange(by = "looic") %>%
  mutate(dLOOIC = looic - min(looic))

# Models 1-3 are not that different
#best.model <- as.numeric(looic.table[1, "ID"])
best.model <- 4
out.best <- out.list[[best.model]]
rm(list = "out.list")

# Missing data posteriors
mcmc_hist(out.best$jags.out$samples, FL.missing.y[1:12], binwidth = 1)
mcmc_hist(out.best$jags.out$samples, FL.missing.y[13:24], binwidth = 1)
mcmc_hist(out.best$jags.out$samples, FL.missing.y[25:36], binwidth = 1)
mcmc_hist(out.best$jags.out$samples, FL.missing.y[37:48], binwidth = 1)
mcmc_hist(out.best$jags.out$samples, FL.missing.y[49:61], binwidth = 1)

if (str_detect(model.file.names[best.model], "2Us"))
  # Rate of change for before and after 2012
  mcmc_dens(out.best$jags.out$samples, c("mean.U1", "mean.U2"))

if (str_detect(model.file.names[best.model], "_U_"))
  # Rate of change for before and after 2012
  mcmc_dens(out.best$jags.out$samples, c("mean.U1"))


jags.data.df <- data.frame(year = c(FL.years[FL.year.1toT.vec == 1], 
                                    FL.years[FL.year.1toT.vec > 1]),
                           beach = c(FL.ID.1, FL.ID.2plus),
                           nest = c(FL.y.1, FL.y.2plus))

k <- 1
for (k in 1:nrow(FL.summary.years)){
  jags.data.df %>% 
    filter(beach == k) -> jags.data.1beach 
  
  p.1beach <- ggplot(jags.data.1beach) +
    geom_point(aes(x = year, y = nest))
  
}

# Look at the number per year as the sum of all beaches
jags.data.df %>%
  arrange(by = beach) %>%
  pivot_wider(id_cols = year, 
              id_expand = T,
              names_from = beach,
              values_from = nest) %>%
  column_to_rownames("year") -> nest.data.beach

nest.per.year <- data.frame(year = seq(from = min(FL.summary.years$year.1),
                                       to = max(FL.summary.years$year.2)),
                            n.nests = rowSums(nest.data.beach, na.rm = T),
                            n.beaches = rowSums(!is.na(nest.data.beach))) %>%
  mutate(nests.per.beach = n.nests/n.beaches,
         var = apply(nest.data.beach, FUN = var, MARGIN = 1, na.rm = T))

N.hats.sum <- data.frame(Year = as.numeric(rownames(nest.data.beach)),
                         N = colSums(out.best$jags.out$mean$N, na.rm = T),
                         low.N = colSums(out.best$jags.out$q2.5$N, na.rm = T),
                         high.N  = colSums(out.best$jags.out$q97.5$N, na.rm = T))

p.N.hats <- ggplot(N.hats.sum) +
  geom_ribbon(aes(x = Year, 
                  ymin = low.N, ymax = high.N),
              fill = "orange", alpha = 0.5) +
  geom_point(aes(x = Year, y = N)) +
  geom_path(aes(x = Year, y = N)) 
  #geom_errorbar(aes(x = Year, ymin = low.N, ymax = high.N))

