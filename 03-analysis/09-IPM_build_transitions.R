
library(here)
library(tidyverse)
library(posterior)
library(bayesplot)
library(foreach)
library(doParallel)

registerDoParallel()

# load mcmc results
posterior = readRDS(here::here('02-data',
                               '03-results',
                               'real_fits',
                               'posterior_draws.rds'))



size_metadata = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.rds')) %>%
  # convert to metric
  mutate(bin_midpoint = bin_midpoint * 0.0254,
         bin_lower = bin_lower * 0.0254,
         bin_upper = bin_upper * 0.0254,
         dbh_m.mean = dbh_in.mean * 0.0254) %>%
  
  # for testing, restrict to ~95th percentile observed size class
  filter(dbh_m.mean < 1.25) 
# I've played around with only using "smaller" 
# sizes in the IPM transition kernel, because I don't really believe some 
# of the transition rates for very large trees (>1.25m dbh), because I doubt 
# that the relationship between size and fecundity really keeps increasing 
# exponentially that far into the tail of DBH. The 
# parameter-medians-all-subplots results are not sensitive to the range of 
# sizes modeled in the transition matrix. The 
# posterior-params-hypothetical-subplots are slightly sensitive, in that the 
# posterior lambda distributions for the different scenarios do shift (in 
# particular, the posterior lambda distribution for burned subplots shifts 
# downward slightly), but the core findings (strong negative effect of fire 
# reducing lambda below 1, weaker negative effects of WPBR and high BA reducing 
# lambda to ~1, high positive effect of low BA, weak effects of drought and 
# dryness) remain unchanged. For simplicity and clarity when writing a paper, 
# i'd like the post-estimation IPM to have the same structure as the IPM used 
# to estimate fecundity, and the results are not sensitive to the decision 
# about what sizes to include in the IPM, so I'm leaving in the large size bins 
# for the final publication analysis.

size_metadata

size_metadata$r = 
  c(readRDS(here::here('02-data',
                       '02-for_analysis',
                       'pila_training.rds'))$r,
    rep(0, times = 8))

head(size_metadata)



#### observed fixed and random #################################################

subplots = 
  readRDS(here::here('02-data', '01-preprocessed', 'subplot_data.rds'))%>%
  mutate(ba_scaled = as.numeric(scale(ba_ft2ac)),
         cwd_dep90_scaled = as.numeric(scale(cwd_departure90)),
         cwd_mean_scaled = as.numeric(scale(cwd_mean)),
         intercept = 1) %>%
  select(plot_id, subp_id, lat, lon, ecosubcd, intercept, fire, wpbr, ba_scaled, 
         cwd_dep90_scaled,cwd_mean_scaled)



subplots.pila = 
  subplots %>%
  right_join(
    readRDS(here::here('02-data', '02-for_analysis', 'union_plots.rds'))
  ) %>%
  right_join(
    readRDS(here::here('02-data', '02-for_analysis', 'union_ecosubs.rds'))
  )


A_subplot = 
  foreach(draw = 1:nrow(posterior),
          .packages = c('tidyverse')) %dopar% {
            
            # get beta_s for the current draw
            beta_s = 
              posterior %>%
              select(contains('beta_s')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            beta_g = 
              posterior %>%
              select(contains('beta_g')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            beta_f = 
              posterior %>%
              select(contains('beta_f')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            plotEffect_s = 
              posterior %>%
              select(contains('plotEffect_s')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            plotEffect_g = 
              posterior %>%
              select(contains('plotEffect_g')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            plotEffect_f = 
              posterior %>%
              select(contains('plotEffect_f')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            ecoEffect_s = 
              posterior %>%
              select(contains('ecoEffect_s')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            ecoEffect_g = 
              posterior %>%
              select(contains('ecoEffect_g')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            ecoEffect_f = 
              posterior %>%
              select(contains('ecoEffect_f')) %>%
              slice(draw) %>%
              as.data.frame() %>%
              as.numeric()
            
            sigmaEpsilon_g = 
              posterior %>%
              slice(draw) %>%
              pull(sigmaEpsilon_g) %>%
              as.numeric()
            
            A = 
              array(dim = list(nrow(size_metadata),
                               nrow(size_metadata),
                               nrow(subplots.pila)),
                    dimnames = list('class_to' = 1:nrow(size_metadata),
                                    'class_from' = 1:nrow(size_metadata),
                                    'subplot' = 1:nrow(subplots.pila)),
                    data = 
                      sapply(X = 1:nrow(subplots.pila),
                             FUN = function(subplot){
                               
                               # construct explanatory variable matrix for survival 
                               # for teh current subplot
                               X = 
                                 subplots.pila %>%
                                 slice(subplot) %>%
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
                               
                               # calculate size_from length vector of survival 
                               # probabilities on this subplot with this parameter draw
                               p = 
                                 boot::inv.logit(as.numeric(X %*% beta_s) +
                                                   ecoEffect_s[subplots.pila$ecosub.i[subplot]]+
                                                   plotEffect_s[subplots.pila$plot_id.i][subplot])
                               
                               mu = as.numeric(X %*% beta_g)+
                                 ecoEffect_g[subplots.pila$ecosub.i[subplot]]+
                                 plotEffect_g[subplots.pila$plot_id.i[subplot]]
                               
                               f = 
                                 exp(as.numeric(X %*% beta_f)+
                                       ecoEffect_f[subplots.pila$ecosub.i[subplot]]+
                                       plotEffect_f[subplots.pila$plot_id.i[subplot]])
                               
                               sapply(X = 1:nrow(size_metadata),
                                      FUN = function(class_from){
                                        
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
                                        
                                        sapply(X = 1:nrow(size_metadata),
                                               FUN = function(class_to){
                                                 
                                                 # for testing
                                                 #paste0('to:',class_to,',from:',class_from,
                                                 #       ',subplot:',subplot)
                                                 transition_prob = 
                                                   # survival of each from class
                                                   (p[class_from] *
                                                      # prob of growth from to
                                                      g[class_to]) +
                                                   # number of new recruits
                                                   (f[class_from] *
                                                      size_metadata$r[class_to])
                                                 return(transition_prob)
                                                 
                                               })
                                        
                                      })
                               
                             }))
            
            return(A)
          }

A_subplot = 
  array(dim = c(nrow(size_metadata), # sizeclass to
                nrow(size_metadata), # sizeclass from
                nrow(subplots.pila), # subplots
                length(A_subplot)), # posterior draws
        dimnames = list('class_to' = 1:nrow(size_metadata),
                        'class_from' = 1:nrow(size_metadata),
                        'subplot' = 1:nrow(subplots.pila),
                        'draw' = 1:length(A_subplot)),
        data = 
          sapply(X = 1:length(A_subplot),
                 FUN = function(draw){
                   
                   return(A_subplot[[draw]])
                   
                 }))

saveRDS(A_subplot,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'subplot_As.rds'))

#### using hypothetical subplots (only fixed effects) ##########################

head(subplots)

hypothetical_subplots = 
  data.frame(intercept = rep(1, times = 9),
             fire = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
             wpbr = c(FALSE, FALSE, T, F, F, F, F, F, F),
             ba_scaled = c(0, 0, 0, -1, 1, 0, 0, 0, 0),
             cwd_dep90_scaled = c(0, 0, 0, 0, 0, -1, 1, 0, 0),
             cwd_mean_scaled = c(0, 0, 0, 0, 0, 0, 0, -1, 1),
             subp_id = 1:9,
             name = c('Undisturbed', 'Fire', 'WPBR', 'Low BA', 'High BA',
                      'Low Drought', 'High Drought', 'Wet Site', 'Dry Site'))

A_hypotheticals = 
  array(dim = c(nrow(size_metadata), # sizeclass to
                nrow(size_metadata), # sizeclass from
                nrow(hypothetical_subplots), # subplots
                nrow(posterior)), # posterior draws
        dimnames = list('class_to' = 1:nrow(size_metadata),
                        'class_from' = 1:nrow(size_metadata),
                        'subplot' = 1:nrow(hypothetical_subplots),
                        'draw' = 1:nrow(posterior)),
        data = 
          sapply(X = 1:nrow(posterior),
                 FUN = function(draw){
                   
                   # get beta_s for the current draw
                   beta_s = 
                     posterior %>%
                     select(contains('beta_s')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   beta_g = 
                     posterior %>%
                     select(contains('beta_g')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   beta_f = 
                     posterior %>%
                     select(contains('beta_f')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaEpsilon_g = 
                     posterior %>%
                     slice(draw) %>%
                     pull(sigmaEpsilon_g) %>%
                     as.numeric()
                   
                   sapply(X = 1:nrow(hypothetical_subplots),
                          FUN = function(subplot){
                            
                            # construct explanatory variable matrix for survival 
                            # for the current subplot
                            X = 
                              hypothetical_subplots %>%
                              slice(subplot) %>%
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
                            
                            # calculate size_from length vector of survival 
                            # probabilities on this subplot with this parameter draw
                            p = boot::inv.logit(as.numeric(X %*% beta_s))
                            
                            mu = as.numeric(X %*% beta_g)
                            
                            f = exp(as.numeric(X %*% beta_f))
                            
                            sapply(X = 1:nrow(size_metadata),
                                   FUN = function(class_from){
                                     
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
                                     
                                     sapply(X = 1:nrow(size_metadata),
                                            FUN = function(class_to){
                                              
                                              transition_prob = 
                                                # survival of each from class
                                                (p[class_from] *
                                                # prob of growth from to
                                                g[class_to]) +
                                                # number of new recruits
                                                (f[class_from] *
                                                 size_metadata$r[class_to])
                                              return(transition_prob)
                                              
                                              # for testing
                                              #paste0('d:',draw,'|s:',subplot,
                                              #  '|f:',class_from,'|t:',class_to)
                                            })
                                   })
                            })
                 }))


saveRDS(A_hypotheticals,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'hypothetical_As.rds'))
