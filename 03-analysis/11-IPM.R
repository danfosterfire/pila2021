
library(here)
library(tidyverse)
library(posterior)
library(bayesplot)
library(foreach)
library(doParallel)


# load mcmc results
pila_fit = readRDS(here::here('02-data', '03-results', 'real_fits', 'pila.rds'))

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

head(size_metadata)

subplots = 
  readRDS(here::here('02-data', '01-preprocessed', 'subplot_data.rds'))%>%
  mutate(ba_scaled = as.numeric(scale(ba_ft2ac)),
         cwd_dep90_scaled = as.numeric(scale(cwd_departure90)),
         cwd_mean_scaled = as.numeric(scale(cwd_mean)),
         intercept = 1) %>%
  select(plot_id, subp_id, lat, lon, ecosubcd, intercept, fire, wpbr, ba_scaled, 
         cwd_dep90_scaled,cwd_mean_scaled)

posterior = as_draws_df(pila_fit$draws())

posterior %>%
  select(contains('beta_s')) %>%
  slice(1) %>%
  as.data.frame() %>%
  as.numeric()

explan_data =  
  subplots %>%
  expand(nesting(subp_id, lat, lon, ecosubcd, intercept, fire, wpbr,
                 ba_scaled, cwd_dep90_scaled, cwd_mean_scaled),
         nesting(size_metadata %>%
                   select(dbh_m.mean))) 



#### only pila subplots, median parameter estimates, fixed and random ##########

names(subplots)

subplots.pila = 
  subplots %>%
  right_join(
    readRDS(here::here('02-data', '02-for_analysis', 'union_plots.rds'))
  ) %>%
  right_join(
    readRDS(here::here('02-data', '02-for_analysis', 'union_ecosubs.rds'))
  )


# get beta_s for the current draw
beta_s.med = 
  posterior %>%
  select(contains('beta_s')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

beta_g.med = 
  posterior %>%
  select(contains('beta_g')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

beta_f.med = 
  posterior %>%
  select(contains('beta_f')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

plotEffect_s.med = 
  posterior %>%
  select(contains('plotEffect_s')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

plotEffect_g.med = 
  posterior %>%
  select(contains('plotEffect_g')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

plotEffect_f.med = 
  posterior %>%
  select(contains('plotEffect_f')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

ecoEffect_s.med = 
  posterior %>%
  select(contains('ecoEffect_s')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

ecoEffect_g.med = 
  posterior %>%
  select(contains('ecoEffect_g')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

ecoEffect_f.med = 
  posterior %>%
  select(contains('ecoEffect_f')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()


sigmaEpsilon_g.med = 
  posterior %>%
  summarise_all(median) %>%
  pull(sigmaEpsilon_g) %>%
  as.numeric()
            
subplot_lambdas.med = 
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
           boot::inv.logit(as.numeric(X %*% beta_s.med) +
                             ecoEffect_s.med[subplots.pila$ecosub.i[subplot]]+
                             plotEffect_s.med[subplots.pila$plot_id.i][subplot])
         
         mu = as.numeric(X %*% beta_g.med)+
           ecoEffect_g.med[subplots.pila$ecosub.i[subplot]]+
           plotEffect_g.med[subplots.pila$plot_id.i[subplot]]
         
         f = 
           exp(as.numeric(X %*% beta_f.med)+
                 ecoEffect_f.med[subplots.pila$ecosub.i[subplot]]+
                 plotEffect_f.med[subplots.pila$plot_id.i[subplot]])
         
         A.subplot = 
           sapply(X = 1:nrow(size_metadata),
                FUN = function(class_from){
                  
                  g = 
                    ((pnorm(size_metadata$bin_upper,
                            mu[class_from],
                            sigmaEpsilon_g.med) - 
                        pnorm(size_metadata$bin_lower,
                              mu[class_from],
                              sigmaEpsilon_g.med))/
                       (1-pnorm(0,
                                mu[class_from],
                                sigmaEpsilon_g.med)))
                  
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
         lambda.subplot = max(as.numeric(eigen(A.subplot)$values))
         return(lambda.subplot)
       })
 

subplots.pila$lambda_postmed = 
  subplot_lambdas.med

postmed_lambda_distribution = 
  ggplot(data = subplots.pila,
       aes(x = lambda_postmed))+
  geom_density()+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 2.5))+
  geom_vline(xintercept = 1, color = 'grey', lty = 2, lwd = 1)+
  labs(y = 'Density', x = 'Lambda')

postmed_lambda_distribution_ppt = 
  ggplot(data = subplots.pila,
         aes(x = lambda_postmed))+
  geom_histogram()+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 2.5))+
  geom_vline(xintercept = 1, color = 'grey', lty = 2, lwd = 1)+
  labs(y = 'N subplots', x = 'Lambda')+
  theme(text = element_text(size = 18))

ggsave(postmed_lambda_distribution_ppt,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'subplot_lambdas.png'),
       height = 4, width = 6.5, units = 'in')

summary(subplots.pila$lambda_postmed)

# what proportion of plots is the pop predicted to decline on
length(subplots.pila$lambda_postmed[subplots.pila$lambda_postmed<1])/
  length(subplots.pila$lambda_postmed)

ggsave(postmed_lambda_distribution,
       filename = here::here('04-communication',
                  'figures',
                  'manuscript',
                  'postmed_lambda_distribution.png'),
       height = 4, width = 6.5, units = 'in')

postmed_lambda_distribution

library(sf)
library(scales)
subplots.pila.sf = 
  subplots.pila %>%
  st_as_sf(coords = c('lon', 'lat'),
           crs = 4269)

ggplot(data = 
         subplots.pila.sf,
       aes(color = lambda_postmed))+
  geom_sf(size = 1)+
  scale_color_viridis_c(limits = c(0.9,1.5), oob = scales::squish)+
  theme_minimal()




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
                4000), # posterior draws
        dimnames = list('class_to' = 1:nrow(size_metadata),
                        'class_from' = 1:nrow(size_metadata),
                        'subplot' = 1:nrow(hypothetical_subplots),
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



hypothetical_lambdas = 
  hypothetical_subplots %>%
  expand(nesting(subp_id, name, intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled, 
                 cwd_mean_scaled),
         data.frame(draw = 1:4000))


hypothetical_lambdas$lambda = 
  sapply(X = 1:nrow(hypothetical_subplots),
         FUN = function(subplot){
           sapply(X = 1:4000,
                  FUN = function(draw){
                    # paste0('s:',subplot,'d:',draw) for testing
                    max(as.numeric(eigen(A_hypotheticals[,,subplot,draw])$values))
                  })
         }) %>%
  as.numeric()


saveRDS(hypothetical_lambdas,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'hypothetical_lambdas.rds'))

hypothetical_lambdas = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'hypothetical_lambdas.rds'))

pretty_names = 
  hypothetical_subplots$name
names(pretty_names) = hypothetical_subplots$subp_id

head(hypothetical_lambdas)

hypothetical_lambda_summary = 
  hypothetical_lambdas %>%
  group_by(subp_id, name) %>%
  summarise(lambda.med = median(lambda),
            lambda.05 = quantile(lambda, probs = 0.05),
            lambda.95 = quantile(lambda, probs = 0.95))

hypothetical_lambda_summary

write.csv(hypothetical_lambda_summary,
          here::here('04-communication', 
                     'tables',
                     'hypothetical_lambdas_summary.csv'),
          row.names = FALSE)

hypothetical_lambdas_plot = 
  ggplot(data = 
         hypothetical_lambdas,
       aes(x = lambda))+
  geom_density(lwd = 1)+
  geom_vline(xintercept = 1, color = 'grey', lty = 2, lwd = 1)+
  theme_minimal()+
  facet_grid(subp_id~., scales = 'free_y',
             labeller = labeller(subp_id = pretty_names))+
  scale_x_continuous(limits = c(0.5,3))+
  coord_cartesian(xlim = c(0.9, 1.5))+
  theme(axis.text.y = element_blank())+
  labs(x = 'Lambda', y = 'Posterior Density')

hypothetical_lambdas_plot

ggsave(hypothetical_lambdas_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'hypotheticals_lambda_post.png'),
       height = 7.5, width = 4, units = 'in')

hypothetical_lambdas_ppt = 
  ggplot(data = 
         hypothetical_lambdas,
       aes(x = lambda))+
  geom_density(lwd = 1)+
  geom_vline(xintercept = 1, color = 'grey', lty = 2, lwd = 1)+
  theme_minimal()+
  facet_grid(.~subp_id, scales = 'free_x',
             labeller = labeller(subp_id = pretty_names))+
  scale_x_continuous(limits = c(0.5,3))+
  theme(axis.text.x = element_blank(),
        text = element_text(size = 16))+
  labs(x = 'Lambda', y = 'Posterior Density')+
  coord_flip(xlim = c(0.9, 1.4))

hypothetical_lambdas_ppt

ggsave(hypothetical_lambdas_ppt,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'hypotheticals_lambda_post.png'),
       height = 4, width = 11, units = 'in')



