library(here)
library(tidyverse)
library(posterior)
library(bayesplot)
library(data.table)

#### build transitions #########################################################

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
  readRDS(here::here('02-data',
                       '02-for_analysis',
                       'pila_training.rds'))$r

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
                     
                     X = 
                       plots.pila %>%
                       slice(plot) %>%
                       tidyr::expand(nesting(intercept, fire, wpbr, ba_scaled,
                                      cwd_dep90_scaled,cwd_mean_scaled),
                              dbh = size_metadata$bin_midpoint) %>%
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
                       boot::inv.logit(as.numeric(X %*% beta_s) +
                                         ecoEffect_s[plots.pila$ecosub.i[plot]]+
                                         plotEffect_s[plots.pila$plot_id.i][plot])
                     
                     # calculate vector of mean size at time 2 for each size 
                     # class on this plot with this parameter draw
                     mu = as.numeric(X %*% beta_g)+
                       ecoEffect_g[plots.pila$ecosub.i[plot]]+
                       plotEffect_g[plots.pila$plot_id.i[plot]]
                     
                     # calculate vector of fecundity for each size class on this 
                     # plot with this parameter draw
                     f = 
                       exp(as.numeric(X %*% beta_f)+
                             ecoEffect_f[plots.pila$ecosub.i[plot]])
                     
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

# save this so we don't have to rebuild it every time
saveRDS(A_observed,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'A_observed.rds'))

A_observed = readRDS(here::here('02-data',
                                '03-results',
                                'real_fits',
                                'A_observed.rds'))


transitions.df = 
  as.data.table(A_observed) %>% as_tibble()

head(transitions.df)

transitions.df %>%
  group_by(class_to, class_from) %>%
  summarise(value = median(value)) %>%
  ungroup() %>%
  left_join(size_metadata %>%
              select(class_to = bin_id, 
                     class_to_midpoint = bin_midpoint) %>%
              mutate(class_to = as.character(class_to))) %>%
  left_join(size_metadata %>%
              select(class_from = bin_id, 
                     class_from_midpoint = bin_midpoint) %>%
              mutate(class_from = as.character(class_from))) %>%
  ggplot(aes(x = class_from_midpoint, y = class_to_midpoint, fill = value))+
  geom_tile()+
  scale_fill_viridis_c()

#### lambda ####################################################################


lambda_observed = 
  sapply(X = 1:nrow(plots.pila),
         FUN = function(plot){
           A_plot = A_observed[,,plot]
           lambda_plot = max(as.numeric(Re(eigen(A_plot)$values)))
           return(lambda_plot)
         })

summary(lambda_observed)

lambda_observed_plot = 
  ggplot(data.frame(lambda_observed),
       aes(x = lambda_observed))+
  geom_histogram(fill = 'lightgrey')+
  geom_vline(xintercept = 1, color = 'red', lty = 2, lwd = 1)+
  theme_minimal()+
  labs(y = 'Number of plots', x = 'Asymptotic Population Growth Rate')

lambda_observed_plot

ggsave(lambda_observed_plot,
       filename = here::here('04-communication',
                  'figures',
                  'manuscript',
                  'lambda_observed.png'),
       height = 4.5, width = 5.5, units = 'in')

# proportion of plots where lambda < 1
length(lambda_observed[lambda_observed<1])/length(lambda_observed)

#### ssd #######################################################################

# stable size distribution
ssd_observed = 
  matrix(nrow = nrow(size_metadata),
         ncol = nrow(plots.pila),
         byrow = FALSE,
         data = 
           sapply(X = 1:nrow(plots.pila),
                  FUN = function(plot){
                    A.plot = A_observed[,,plot]
                    # from supplamentory materials for merow et al 2014 
                    # "On using integral projection models..."
                    w.eigen = Re(eigen(A.plot)$vectors[,1])
                    ssd = w.eigen / sum(w.eigen)
                    return(ssd)
                  }))

ssd_observed.df = 
  expand.grid('plot' = 1:nrow(plots.pila),
              'sizeclass' = 1:nrow(size_metadata)) %>%
  left_join(size_metadata %>%
              mutate(bin_midpoint_cm = bin_midpoint*100) %>%
              select(sizeclass = bin_id, bin_midpoint_cm))

ssd_observed.df$class_proportion = 
  
  as.numeric(
    sapply(X = 1:nrow(size_metadata),
           FUN = function(sizeclass){
             
             return(ssd_observed[sizeclass,])
             
           })
  )

# stable size distribution is inverse J not surprising
ssd_observed.df %>%
  group_by(bin_midpoint_cm, sizeclass) %>%
  summarise(prop.med = median(class_proportion),
            prop.05 = quantile(class_proportion, 0.05),
            prop.95 = quantile(class_proportion, 0.95)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint_cm))+
  geom_point(aes(y = prop.med), size = 3)+
  geom_errorbar(aes(ymin = prop.05, ymax = prop.95),
                width = 1)+
  theme_minimal()

# stable size distribution is inverse J not surprising
ssd_observed.df %>%
  group_by(bin_midpoint_cm, sizeclass) %>%
  summarise(prop.med = median(class_proportion),
            prop.05 = quantile(class_proportion, 0.05),
            prop.95 = quantile(class_proportion, 0.95)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint_cm))+
  geom_point(aes(y = prop.med), size = 3)+
  geom_errorbar(aes(ymin = prop.05, ymax = prop.95),
                width = 1)+
  theme_minimal()+
  scale_y_log10()

summary(ssd_observed.df)

#### reproductive value ########################################################

repro_observed = 
  matrix(nrow = nrow(size_metadata),
         ncol = nrow(plots.pila),
         byrow = FALSE,
         data = 
           sapply(X = 1:nrow(plots.pila),
                  FUN = function(plot){
                    A.plot = A_observed[,,plot]
                    # from supplementary materials for merow et al 2014 
                    # "On using integral projection models..."
                    v.eigen = Re(eigen(t(A.plot))$vectors[,1])
                    rv = v.eigen / v.eigen[1]
                    return(rv)
                  }),
         dimnames = list(size_metadata$bin_midpoint,
                         plots.pila$plot_id))

repro_observed.df = 
  expand.grid('plot' = 1:nrow(plots.pila),
              'sizeclass' = 1:nrow(size_metadata)) %>%
  left_join(size_metadata %>%
              select(sizeclass = bin_id, bin_midpoint))

repro_observed.df$reproductive_value = 
  
  as.numeric(
    sapply(X = 1:nrow(size_metadata),
           FUN = function(sizeclass){
             
             return(repro_observed[sizeclass,])
             
           })
  )

repro_observed.df %>% summary()

repro_observed.df %>%
  
  group_by(bin_midpoint, sizeclass) %>%
  
  # there's a couple of NA plots for reproductive value, looks like cases 
  # where numerical errors are resulting in a divide by zero when going from 
  # v.eigen to reproductive value? 
  summarise(repr.med = median(reproductive_value, na.rm = TRUE),
            
            # there's a couple of NAs in h
            repr.05 = quantile(reproductive_value, 0.25, na.rm = TRUE),
            repr.95 = quantile(reproductive_value, 0.75, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint))+
  geom_point(aes(y = repr.med), size = 3)+
  geom_errorbar(aes(ymin = repr.05, ymax = repr.95),
                width = 0.01)+
  theme_minimal()

length(repro_observed.df$reproductive_value[repro_observed.df$reproductive_value<0])

repro_observed.df %>% filter(reproductive_value<0) %>% pull(plot) %>% unique()

plots.pila %>% 
  filter(is.element(plot_id.i, c(257, 546, 1005, 971, 634))) %>%
  print(width = Inf)

# looking at the bad plots
A.plot = A_observed[,,634]
# they all have mixed signs and small magnitudes in the left eigenvector
v.eigen = Re(eigen(t(A.plot))$vectors[,1])
rv = v.eigen / v.eigen[1]

#### sensitivity and elasticity ################################################
v.dot.w_observed = 
  
  sapply(X = 1:nrow(plots.pila),
         FUN = function(plot){
           sum(ssd_observed[,plot] * repro_observed[,plot])*0.127
         })

sens_observed = 
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(plots.pila)),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'plot' = 1:nrow(plots.pila)),
        data = 
          sapply(X = 1:nrow(plots.pila),
                 FUN = function(plot){
                   outer(repro_observed[,plot], ssd_observed[,plot])/
                     v.dot.w_observed[plot]
                 }))


elas_observed = 
  
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(plots.pila)),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'plot' = 1:nrow(plots.pila)),
        data = 
          sapply(X = 1:nrow(plots.pila),
                 FUN = function(plot){
                   matrix(as.vector(sens_observed[,,plot])*
                            as.vector(A_observed[,,plot])/
                            lambda_observed[plot])
                 }))


sens_elas_observed.df = 
  expand.grid(size_to = size_metadata$bin_id,
              size_from = size_metadata$bin_id) %>%
  left_join(size_metadata %>%
              select(size_to = bin_id, 
                     bin_midpoint_to_m = bin_midpoint)) %>%
  left_join(size_metadata %>%
              select(size_from = bin_id,
                     bin_midpoint_from_m = bin_midpoint)) %>%
  expand(nesting(size_to, size_from, bin_midpoint_to_m, bin_midpoint_from_m),
         plot = 1:nrow(plots.pila)) %>%
  
  arrange(plot, size_to, size_from)

sens_elas_observed.df$sensitivity = 
  
  sapply(X = 1:nrow(plots.pila),
         FUN = function(plot){
           
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    
                    sapply(X = 1:nrow(size_metadata),
                           FUN = function(size_from){
                             
                             sens_observed[size_to, size_from, plot]
                             
                           })
                    
                  })
           
         }) %>%
  as.vector()

sens_elas_observed.df$elasticity = 
  
  sapply(X = 1:nrow(plots.pila),
         FUN = function(plot){
           
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    
                    sapply(X = 1:nrow(size_metadata),
                           FUN = function(size_from){
                             
                             elas_observed[size_to, size_from, plot]
                             
                           })
                    
                  })
           
         }) %>%
  as.vector()

sens_elas_observed.df %>%
  group_by(size_to, size_from, bin_midpoint_to_m, bin_midpoint_from_m) %>%
  summarise(sensitivity = median(sensitivity, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint_from_m,
           y = bin_midpoint_to_m,
           fill = sensitivity))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c()

sens_elas_observed.df %>%
  group_by(size_to, size_from, bin_midpoint_to_m, bin_midpoint_from_m) %>%
  summarise(elasticity = median(elasticity, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint_from_m,
           y = bin_midpoint_to_m,
           fill = elasticity))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c()


