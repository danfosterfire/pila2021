
library(here)
library(tidyverse)
library(posterior)
library(bayesplot)
library(foreach)
library(doParallel)


# load mcmc results
pila_fit = readRDS(here::here('02-data', '03-results', 'real_fits', 'usgs', 'pila.rds'))

pila_data = readRDS(here::here('02-data', '02-for_analysis', 'usgs', 'pila_data.rds'))


size_metadata = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'usgs',
                     'size_metadata.rds'))  %>%
  # convert to meters
  mutate(lower = lower * 0.01,
         upper = upper * 0.01,
         dbh_med = dbh_med * 0.01)


head(size_metadata)

posterior = as_draws_df(pila_fit$draws())


names(posterior)

#### all parameter estimates, all plots, fixed and random ######################


A_array = 
  array(dim = c(nrow(size_metadata), # sizeclass to
                nrow(size_metadata), # sizeclass from
                pila_data$P, # plots
                4000), # posterior draws
        dimnames = list('class_to' = 1:nrow(size_metadata),
                        'class_from' = 1:nrow(size_metadata),
                        'plot' = 1:pila_data$P,
                        'draw' = 1:4000),
        data = 
          sapply(X = 1:4000,
                 FUN = function(draw){
                   
                   # get beta_s for the current draw
                   beta_s = 
                     posterior %>%
                     select(contains('beta_s')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   plotEffects_s = 
                     posterior %>%
                     select(contains('plotEffect_s')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   beta_g = 
                     posterior %>%
                     select(contains('beta_g')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   
                   plotEffects_g = 
                     posterior %>%
                     select(contains('plotEffect_g')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   beta_f = 
                     posterior %>%
                     select(contains('beta_f')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   
                   plotEffects_f = 
                     posterior %>%
                     select(contains('plotEffect_f')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaEpsilon_g = 
                     posterior %>%
                     slice(draw) %>%
                     pull(sigmaEpsilon_g) %>%
                     as.numeric()
                   
                   mu_r = posterior$mu_r[draw]
                   
                   sigmaEpsilon_r = posterior$sigmaEpsilon_r[draw]
                   
                   # recruit size kernel
                    r = ((pnorm(log(size_metadata$upper*100),
                                mu_r,
                                sigmaEpsilon_r)-
                           pnorm(log(size_metadata$lower*100),
                                 mu_r,
                                 sigmaEpsilon_r)))
                   
                   sapply(X = 1:pila_data$P,
                          FUN = function(plot){
                            
                            # construct explanatory variable matrix for survival 
                            # for the current subplot
                            X = 
                              matrix(nrow = pila_data$M_f, ncol = pila_data$K,
                                     byrow = FALSE,
                                     data = c(rep(1, times = pila_data$M_f),
                                              size_metadata$dbh_med))
                            
                            # calculate size_from length vector of survival 
                            # probabilities on this subplot with this parameter draw
                            p = 
                              boot::inv.logit(as.numeric(X %*% beta_s))+
                              plotEffects_s[plot]
                            
                            mu_size = 
                              as.numeric(X %*% beta_g)+
                              plotEffects_g[plot]
                            
                            f = 
                              exp(as.numeric(X %*% beta_f)+
                                    plotEffects_f[plot])
                            
                            sapply(X = 1:nrow(size_metadata),
                                   FUN = function(class_from){
                                     
                                     g = 
                                       ((pnorm(size_metadata$upper,
                                               mu_size[class_from],
                                               sigmaEpsilon_g) - 
                                           pnorm(size_metadata$lower,
                                                 mu_size[class_from],
                                                 sigmaEpsilon_g))/
                                          (1-pnorm(0,
                                                   mu_size[class_from],
                                                   sigmaEpsilon_g)))
                                     
                                     sapply(X = 1:nrow(size_metadata),
                                            FUN = function(class_to){
                                              
                                              transition_prob = 
                                                # survival of each from class
                                                (p[class_from] *
                                                   # prob of growth from to
                                                   g[class_to]) +
                                                # number of new recruits
                                                (f[class_from] *r[class_to])
                                              return(transition_prob)
                                              
                                              # for testing
                                              #paste0('d:',draw,'|s:',subplot,
                                              #  '|f:',class_from,'|t:',class_to)
                                            })
                                   })
                          })
                 }))

ipm_results = 
  readRDS(here::here('02-data', 
                     '02-for_analysis',
                     'usgs',
                     'plots_spatial.rds')) %>%
  filter(!is.na(plot_id.i)) %>%
  expand(nesting(plot, areaha, lat, lon, plot_id.i),
         draw = 1:4000)

ipm_results$lambda =
  sapply(X = 1:4000,
         FUN = function(draw){
           sapply(X = 1:pila_data$P,
                  FUN = function(plot){
                    A = A_array[,,plot,draw]
                    lambda = max(as.numeric(eigen(A)$values))
                    return(lambda)
                  })
         }) %>%
  as.numeric()



ipm_density = 
  ggplot(data = ipm_results,
       aes(x = lambda, color = plot))+
  geom_density()+
  theme_minimal()+
  geom_vline(xintercept = 1, lty = 2, lwd = 1)+
  labs(color = 'Plot', y = 'Posterior probability density', x = 'Lambda')

ggsave(ipm_density,
       filename = here::here('04-communication',
                             'figures',
                             'report',
                             'ipm_density.png'),
       height = 4, width = 6.5, units = 'in')

head(ipm_results)

ipm_ci = 
  ipm_results %>%
  group_by(plot) %>%
  summarise(lambda.05 = quantile(lambda, probs = c(0.05)),
            lambda.25 = quantile(lambda, probs = c(0.25)),
            lambda.50 = median(lambda),
            lambda.75 = quantile(lambda, probs = 0.75),
            lambda.95 = quantile(lambda, probs = 0.95)) %>%
  ungroup() %>%
  ggplot()+
  geom_segment(aes(y = plot, yend = plot, x = lambda.05, xend = lambda.95),
               lwd = 1, color = 'slategray3')+
  geom_segment(aes(y = plot, yend = plot, x = lambda.25, xend = lambda.75),
               lwd = 2, color = 'slategray4')+
  geom_point(aes(x = lambda.50, y = plot), size = 5)+
  theme_minimal()+
  labs(x = 'Lambda', y = 'Plot')

ggsave(ipm_ci,
       filename = here::here('04-communication',
                             'figures',
                             'report',
                             'ipm_ci.png'),
       height = 4, width = 6.5, units = 'in')

saveRDS(ipm_results,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'usgs',
                   'ipm_results.rds'))
