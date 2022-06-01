
#### setup #####################################################################

library(here)
library(posterior)
library(tidyverse)

ba_posterior = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'ba_fit.rds'))$draws() %>%
  as_draws_df()



#### posterior v prior #########################################################

param_names = 
  data.frame(param = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]', 
                       'beta[5]', 'beta[6]'),
             name = c('Intercept', 'Fire', 'WPBR', 'BA', 'Drought', 'Site Dryness'),
             param_order = 1:6)
pretty_param_names = param_names$name
names(pretty_param_names) = param_names$param_order
beta_priorvpost = 
  ba_posterior %>%
  select('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]', 'beta[5]', 'beta[6]') %>%
  rowid_to_column('iter') %>%
  pivot_longer(cols = c(-iter),
               names_to = 'param',
               values_to = 'value') %>%
  left_join(param_names) %>%
  ggplot()+
  geom_density(aes(x = value))+
  geom_density(data = data.frame(prior = rnorm(n = 1000, mean = 0, sd = 10)),
               aes(x = prior), color = 'red')+
  facet_wrap(~param_order, 
             labeller = labeller(param_order = pretty_param_names),
             scales = 'free')

sigmas_priorvpost = 
  ba_posterior %>%
  select('sigma_epsilon', 'sigma_ecosub') %>%
  rowid_to_column('iter') %>%
  pivot_longer(cols = c(-iter),
               names_to = 'param',
               values_to = 'value') %>%
  ggplot()+
  geom_density(aes(x = value))+
  geom_density(data = data.frame(prior = truncnorm::rtruncnorm(n = 1000, mean = 0, sd = 10, a = 0)),
               aes(x = prior), color = 'red')+
  facet_wrap(~param)


#### scenarios #################################################################

hypothetical_plots = 
  data.frame(intercept = rep(1, times = 9),
             fire = c(0, 1, 0, 0, 0, 0, 0, 0, 0),
             wpbr = c(0, 0, 1, 0, 0, 0, 0, 0, 0),
             ba = c(0, 0, 0, -1, 1, 0, 0, 0, 0),
             drought = c(0, 0, 0, 0, 0, -1, 1, 0, 0),
             dryness = c(0, 0, 0, 0, 0, 0, 0, -1, 1),
             name = c('Unstressed', 'Fire', 'WPBR', 'Low BA', 'High BA',
                      'Low Drought', 'High Drought', 'Low Dryness', 'High Dryness'))

predict_mean_delta = 
  function(dataset){
    
    do.call('rbind',
            lapply(X = 1:4000,
                   FUN = function(d){
              
                     beta_d = ba_posterior %>%
                       select(contains('beta')) %>%
                       slice(d) %>%
                       as.numeric()
               
                     
                     mu = as.matrix(dataset[,1:6]) %*% beta_d 
                     
                     results = dataset
                     results$mu = mu
                     
                     return(results)
                     
            }))
    
    
  }

predicted_mean_deltas = 
  do.call('rbind', 
          lapply(X = 1:9,
                 FUN = function(i){
                   predict_mean_delta(hypothetical_plots[i,])
                 }))

predicted_mean_deltas %>%
  mutate(name = factor(name, levels = c('Unstressed', 
                                        'Fire', 'WPBR', 'Low BA', 'High BA',
                                        'Low Drought', 'High Drought',
                                        'Low Dryness', 'High Dryness'))) %>%
  ggplot(aes(x = name, y = mu))+
  geom_hline(yintercept = 0, col = 'black', lty = 2, lwd = 1)+
  geom_violin(fill = 'lightgrey', alpha = 0.75)+
  theme_minimal()+
  labs(y = 'Mean change in BA (m^2 / ha)', x = 'Scenario')


