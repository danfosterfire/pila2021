library(here)
library(tidyverse)
library(posterior)
library(bayesplot)

#### observed plots ############################################################

# load mcmc results
posterior = readRDS(here::here('02-data',
                               '03-results',
                               'real_fits',
                               'posterior_draws.rds'))

# extract parameters
beta_s = 
  posterior %>%
  select(contains('beta_s')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

beta_g = 
  posterior %>%
  select(contains('beta_g')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

beta_f = 
  posterior %>%
  select(contains('beta_f')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

plotEffect_s = 
  posterior %>%
  select(contains('plotEffect_s')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

plotEffect_g = 
  posterior %>%
  select(contains('plotEffect_g')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

plotEffect_f = 
  posterior %>%
  select(contains('plotEffect_f')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

ecoEffect_s = 
  posterior %>%
  select(contains('ecoEffect_s')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

ecoEffect_g = 
  posterior %>%
  select(contains('ecoEffect_g')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

ecoEffect_f = 
  posterior %>%
  select(contains('ecoEffect_f')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

sigmaEpsilon_g = 
  posterior %>%
  summarise_all(median) %>%
  pull(sigmaEpsilon_g) %>%
  as.numeric()

size_metadata = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.rds')) %>%
  # convert to metric
  mutate(bin_midpoint = bin_midpoint * 0.0254,
         bin_lower = bin_lower * 0.0254,
         bin_upper = bin_upper * 0.0254,
         dbh_m.mean = dbh_in.mean * 0.0254) 

size_metadata$r = 
  c(readRDS(here::here('02-data',
                       '02-for_analysis',
                       'pila_training.rds'))$r,
    rep(0, times = 18))

plots.pila = 
  readRDS(here::here('02-data', '01-preprocessed', 'plot_data.rds'))%>%
  mutate(ba_scaled = as.numeric(scale(ba_ft2ac)),
         cwd_dep90_scaled = as.numeric(scale(cwd_departure90)),
         cwd_mean_scaled = as.numeric(scale(cwd_mean)),
         intercept = 1) %>%
  select(plot_id, lat, lon, ecosubcd, intercept, fire, wpbr, ba_scaled, 
         cwd_dep90_scaled,cwd_mean_scaled) %>%
  right_join(
    readRDS(here::here('02-data', '02-for_analysis', 'union_plots.rds'))
  ) %>%
  right_join(
    readRDS(here::here('02-data', '02-for_analysis', 'union_ecosubs.rds'))
  )


A_observed = 
    array(dim = list(nrow(size_metadata),
                     nrow(size_metadata),
                     nrow(plots.pila)),
          dimnames = list('class_to' = 1:nrow(size_metadata),
                          'class_from' = 1:nrow(size_metadata),
                          'plot' = 1:nrow(plots.pila)),
          data = 
            sapply(X = 1:nrow(plots.pila),
                   FUN = function(plot){
                     
                     # construct explanatory variable matrix for vital rate 
                     # functions for the current plot
                     X_sg = 
                       plots.pila %>%
                       slice(plot) %>%
                       expand(nesting(intercept, fire, wpbr, ba_scaled,
                                      cwd_dep90_scaled,cwd_mean_scaled),
                              dbh = size_metadata$dbh_m.mean) %>%
                       mutate(dbh_fire = dbh*fire,
                              dbh_wpbr = dbh*wpbr,
                              dbh_ba = dbh*ba_scaled,
                              dbh_cwd_dep90 = dbh*cwd_dep90_scaled,
                              dbh_cwd_mean = dbh*cwd_mean_scaled) %>%
                       select(intercept, dbh, fire, wpbr, ba_scaled,
                              cwd_dep90_scaled, cwd_mean_scaled, 
                              dbh_fire, dbh_wpbr, dbh_ba,
                              dbh_cwd_dep90, dbh_cwd_mean) %>%
                       as.matrix()
                     
                     X_f = 
                       plots.pila %>%
                       slice(plot) %>%
                       expand(nesting(intercept, fire, wpbr, ba_scaled,
                                      cwd_dep90_scaled,cwd_mean_scaled),
                              dbh = size_metadata$dbh_m.mean) %>%
                       mutate(dbh2 = dbh**2,
                              dbh_fire = dbh*fire,
                              dbh2_fire = dbh2*fire,
                              dbh_wpbr = dbh*wpbr,
                              dbh2_wpbr = dbh2*wpbr,
                              dbh_ba = dbh*ba_scaled,
                              dbh2_ba = dbh2*ba_scaled,
                              dbh_cwd_dep90 = dbh*cwd_dep90_scaled,
                              dbh2_cwd_dep90 = dbh2*cwd_dep90_scaled,
                              dbh_cwd_mean = dbh*cwd_mean_scaled,
                              dbh2_cwd_mean = dbh2*cwd_mean_scaled) %>%
                       select(intercept, dbh, dbh2, fire, wpbr, ba_scaled,
                              cwd_dep90_scaled, cwd_mean_scaled, 
                              dbh_fire,dbh2_fire, dbh_wpbr, dbh2_wpbr,dbh_ba,
                              dbh2_ba, dbh_cwd_dep90,dbh2_cwd_dep90, 
                              dbh_cwd_mean, dbh2_cwd_mean) %>%
                       as.matrix()
                     
                     
                     # calculate vector of survival probabilities for each 
                     # size class on this plot with this parameter draw
                     p = 
                       boot::inv.logit(as.numeric(X_sg %*% beta_s) +
                                         ecoEffect_s[plots.pila$ecosub.i[plot]]+
                                         plotEffect_s[plots.pila$plot_id.i][plot])
                     
                     # calculate vector of mean size at time 2 for each size 
                     # class on this plot with this parameter draw
                     mu = as.numeric(X_sg %*% beta_g)+
                       ecoEffect_g[plots.pila$ecosub.i[plot]]+
                       plotEffect_g[plots.pila$plot_id.i[plot]]
                     
                     # calculate vector of fecundity for each size class on this 
                     # plot with this parameter draw
                     f = 
                       exp(as.numeric(X_f %*% beta_f)+
                             ecoEffect_f[plots.pila$ecosub.i[plot]]+
                             plotEffect_f[plots.pila$plot_id.i[plot]])
                     
                     # loop over each "from" size class
                     sapply(X = 1:nrow(size_metadata),
                            FUN = function(class_from){
                              
                              # growth kernel from this size class into 
                              # each other size class, using the cumulative 
                              # density function as recommended by Doak et al. 2021
                              g = 
                                ((pnorm(size_metadata$bin_upper,
                                        mu[class_from],
                                        sigmaEpsilon_g) - 
                                    pnorm(size_metadata$bin_lower,
                                          mu[class_from],
                                          sigmaEpsilon_g))/
                                   (1-pnorm(0,
                                            mu[class_from],
                                            sigmaEpsilon_g)))
                              
                              # loop over every destination size class; 
                              # now vectorized
                              # calculate the transition kernel between the 
                              # current 'from' class and every 'to' class
                              transition_kern = 
                                # survival of each from class
                                (p[class_from] *
                                   # probability of growth from this class to 
                                   # every other class
                                   g)+(
                                     # total number of new recruits from this class 
                                     f[class_from]*
                                       
                                       # proportion of new recruits falling in the 
                                       # to class
                                       size_metadata$r)
                              
                              return(transition_kern)
                              
                            })
                     
                   }))

# save this so we don't have to rebuild it every time we knit
saveRDS(A_observed,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'A_observed.rds'))
