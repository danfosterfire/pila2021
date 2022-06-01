
#### setup #####################################################################

library(here)
library(posterior)
library(tidyverse)

#### tph #######################################################################

tph_posterior = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'tph_fit.rds'))$draws() %>%
  as_draws_df()

tph_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'tph.rds'))

tph_post_predictive = 
  do.call('bind_rows',
          lapply(X = 1:nrow(tph_posterior),
                 FUN = function(d){
                   
                   beta = tph_posterior %>%
                     select(contains('beta')) %>%
                     slice(d) %>%
                     as.numeric()
                   
                   sigma_ecosub = 
                     tph_posterior$sigma_ecosub[d]
                   
                   sigma_epsilon = 
                     tph_posterior$sigma_epsilon[d]
                 
                   # random effect realizations
                   ecosub_effects = rnorm(n = tph_data$E,
                                          mean = 0, 
                                          sd = sigma_ecosub)
                   
                   # linear predictor
                   XB = tph_data$X %*% beta
                   mu = sapply(X = 1:tph_data$N,
                               FUN = function(i){
                                 XB[i] + ecosub_effects[tph_data$ecosub_id[i]]
                               })
                   
                   # responses
                   Y = rcauchy(n = tph_data$N,
                               location = mu, scale = sigma_epsilon)
                   
                   # results
                   results = tph_data$X %>%
                     as_tibble() %>%
                     mutate(delta_tph_real = tph_data$Y,
                            delta_tph_sim = Y,
                            draw = d) 
                   return(results)
                 }))

tph_post_predictive %>%
  group_by(intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled, cwd_mean_scaled,
           delta_tph_real) %>%
  summarise(delta_tph_sim.05 = quantile(delta_tph_sim, p = 0.05),
            delta_tph_sim.50 = median(delta_tph_sim),
            delta_tph_sim.95 = quantile(delta_tph_sim, p = 0.95)) %>%
  ungroup() %>%
  ggplot(aes(x = delta_tph_sim.50))+
  geom_point(aes(y = delta_tph_real))+
  geom_ribbon(aes(ymin = delta_tph_sim.05, ymax = delta_tph_sim.95),
              alpha = 0.25)+
  theme_minimal()+
  scale_y_continuous(limits = c(-100, 100))


tph_post_predictive %>%
  group_by(intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled, cwd_mean_scaled,
           delta_tph_real) %>%
  summarise(delta_tph_sim.05 = quantile(delta_tph_sim, p = 0.05),
            delta_tph_sim.50 = median(delta_tph_sim),
            delta_tph_sim.95 = quantile(delta_tph_sim, p = 0.95)) %>%
  ungroup() %>%
  mutate(inside_90 = delta_tph_real > delta_tph_sim.05 & 
           delta_tph_real < delta_tph_sim.95) %>%
  summary()


#### ba ########################################################################

ba_posterior = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'ba_fit.rds'))$draws() %>%
  as_draws_df()

ba_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'ba.rds'))

ba_post_predictive = 
  do.call('bind_rows',
          lapply(X = 1:nrow(ba_posterior),
                 FUN = function(d){
                   
                   beta = ba_posterior %>%
                     select(contains('beta')) %>%
                     slice(d) %>%
                     as.numeric()
                   
                   sigma_ecosub = 
                     ba_posterior$sigma_ecosub[d]
                   
                   sigma_epsilon = 
                     ba_posterior$sigma_epsilon[d]
                 
                   # random effect realizations
                   ecosub_effects = rnorm(n = ba_data$E,
                                          mean = 0, 
                                          sd = sigma_ecosub)
                   
                   # linear predictor
                   XB = ba_data$X %*% beta
                   mu = sapply(X = 1:ba_data$N,
                               FUN = function(i){
                                 XB[i] + ecosub_effects[ba_data$ecosub_id[i]]
                               })
                   
                   # responses
                   Y = rnorm(n = ba_data$N,
                               mean = mu, sd = sigma_epsilon)
                   
                   # results
                   results = ba_data$X %>%
                     as_tibble() %>%
                     mutate(delta_ba_real = ba_data$Y,
                            delta_ba_sim = Y,
                            draw = d) 
                   return(results)
                 }))

ba_post_predictive %>%
  group_by(intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled, cwd_mean_scaled,
           delta_ba_real) %>%
  summarise(delta_ba_sim.05 = quantile(delta_ba_sim, p = 0.05),
            delta_ba_sim.50 = median(delta_ba_sim),
            delta_ba_sim.95 = quantile(delta_ba_sim, p = 0.95)) %>%
  ungroup() %>%
  ggplot(aes(x = delta_ba_sim.50, y = delta_ba_real))+
  geom_point(size = 0)+
  geom_ribbon(aes(ymin = delta_ba_sim.05, ymax = delta_ba_sim.95),
              alpha = 0.25)+
  geom_smooth(method = 'lm')+
  geom_abline(intercept = 0, slope = 1, color = 'red')+
  theme_minimal()


ba_post_predictive %>%
  group_by(intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled, cwd_mean_scaled,
           delta_ba_real) %>%
  summarise(delta_ba_sim.05 = quantile(delta_ba_sim, p = 0.05),
            delta_ba_sim.50 = median(delta_ba_sim),
            delta_ba_sim.95 = quantile(delta_ba_sim, p = 0.95)) %>%
  ungroup() %>%
  mutate(inside_90 = delta_ba_real > delta_ba_sim.05 & 
           delta_ba_real < delta_ba_sim.95) %>%
  summary()
