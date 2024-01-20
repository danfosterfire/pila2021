library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

set.seed(110819)

#### survival ##################################################################

fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'surv_fit.rds'))

samples = as_draws_df(fitted_model$draws())

pila_training = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_mort_training.rds'))

# simulating one data set per draw

samples.beta_s = samples %>% select(contains('beta'))
samples.ecoEffect_s = samples %>% select(contains('effect_ecosub')) %>% as.data.frame()
samples.plotEffect_s = samples %>% select(contains('effect_plot')) %>% as.data.frame()


samples.plotEffect_s %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, group = parameter))+
  geom_line(stat = 'density', alpha = 0.4)+
  theme_minimal()

fitted_model$summary(variables = c('beta', 'sigma_plot', 'sigma_ecosub'))

beta_s.med = 
  samples.beta_s %>%
  summarise_all(median) %>%
  as.numeric()

eco_effects = 
  samples.ecoEffect_s %>%
  summarise_all(median) %>%
  as.numeric()

plot_effects = 
  samples.plotEffect_s %>%
  summarise_all(median) %>%
  as.numeric()

logitp.med = 
  as.numeric(pila_training$X %*% beta_s.med)+
  as.numeric(eco_effects)[pila_training$ecosub_id]+
  as.numeric(plot_effects)[pila_training$plot_id]

median_results = 
  pila_training$X %>%
  as_tibble(x = .)

median_results$logitp = logitp.med
median_results$p = boot::inv.logit(logitp.med)
median_results$surv_true = pila_training$surv

median_results = 
  median_results %>%
  mutate(r = dense_rank(p),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20)))

median_results %>%
  ggplot()+
  geom_point(aes(x = r, y = p))+
  geom_jitter(aes(x = r, y = surv_true),
              color = 'red',
              size = 0,
              height = 0.1, width = 0)+
  geom_point(data = 
               median_results %>%
               group_by(r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             color = 'blue',
             aes(x = r, y = surv_true))+
  theme_minimal()

intermediate = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   
                   eco_effects = 
                     rnorm(n = pila_training$E,
                           mean = 0,
                           sd = samples$sigma_ecosub[i])
                   
                   plot_effects = 
                     rnorm(n = pila_training$P,
                           mean = 0,
                           sd = samples$sigma_plot)
                   
                   XB = as.numeric(pila_training$X %*% beta_s)
                   realized_eco = as.numeric(samples.ecoEffect_s[i,])[pila_training$ecosub_id]
                   realized_plot = as.numeric(samples.plotEffect_s[i,])[pila_training$plot_id]
                   
                   logitp = XB  + realized_eco + realized_plot
                     as.numeric(pila_training$X %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_training$ecosub_id]+
                     as.numeric(samples.plotEffect_s[i,])[pila_training$plot_id]
                      eco_effects[pila_training$ecosub_id]+
                     plot_effects[pila_training$plot_id]
                     
                     
                   p = boot::inv.logit(logitp)
                   surv_sim = rbinom(n = pila_training$N, 
                                 size = 1, 
                                 prob = p)
                   
                   result = 
                     pila_training$X %>% 
                     as_tibble(x = .) %>%
                     rownames_to_column('tree_id')
                   result$XB = XB
                   result$realized_eco = realized_eco
                   result$realized_plot = realized_plot
                   result$logitp = logitp
                   result$p = p
                   result$surv_sim = surv_sim
                   result$surv_true = pila_training$surv
                   result$iter = i
                   return(result)
                 })) 

intermediate %>%
  group_by(tree_id, surv_true) %>%
  summarise(p.mean = mean(p)) %>%
  ungroup() %>%
  mutate(cutoff = ifelse(p.mean>=0.85, 'above',
                         ifelse(p.mean>=0.75 & p.mean < 0.85,
                                'around', 
                                'below'))) %>%
  left_join(intermediate) %>%
  ggplot(aes(x = realized_plot, group = tree_id))+
  geom_line(stat = 'density', alpha = 0.2)+
  facet_grid(cutoff ~.)


intermediate %>%
  group_by(tree_id, surv_true) %>%
  summarise(XB.mean = mean(XB),
            eco.mean = mean(realized_eco),
            plot.mean = mean(realized_plot),
            p.mean = mean(p)) %>%
  ungroup() %>%
  ggplot(aes(x = plot.mean, y = p.mean))+
  geom_point(size = 0)+
  theme_minimal()

intermediate %>%
  group_by(tree_id, surv_true) %>%
  summarise(XB.mean = mean(XB),
            p.mean = mean(p)) %>%
  ungroup() %>%
  mutate(r = dense_rank(XB.mean),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20))) %>%
  ggplot()+
  geom_point(aes(x = r, y = boot::inv.logit(XB.mean)))+
  geom_point(data = 
               intermediate %>%
               group_by(tree_id, surv_true) %>%
               summarise(XB.mean = mean(XB), p.mean = mean(p)) %>%
               ungroup() %>%
               mutate(r = dense_rank(XB.mean),
                      r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20))) %>%
               group_by(., r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             aes(x = r, y = surv_true),
             color = 'blue')

mort_retrodictions = 
  intermediate %>%
  group_by(tree_id, surv_true) %>%
  summarise(p.50 = quantile(p, 0.5),
            p.mean = mean(p),
            p.025 = quantile(p, 0.025),
            p.975 = quantile(p, 0.975)) %>%
  ungroup() %>%
  mutate(r = dense_rank(p.mean),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20)))




intermediate %>%
  group_by(tree_id, dbh_m.init, surv_true) %>%
  summarise(p.50 = quantile(p, 0.5),
            p.mean = mean(p),
            p.025 = quantile(p, 0.025),
            p.975 = quantile(p, 0.975)) %>%
  ungroup() %>%
  mutate(dbh_class = cut(dbh_m.init, breaks = seq(from = 0, to = 2.54, by = 0.0254),
                         labels = FALSE)) %>%
  group_by(dbh_class) %>%
  summarise(p.real = mean(surv_true),
            p.50 = mean(p.50),
            p.mean = mean(p.mean),
            p.025 = mean(p.025),
            p.975 = mean(p.975)) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(dbh_class)))+
  geom_point(aes(y = p.real), color = 'red')+
  geom_point(aes(y = p.mean), color = 'black')+
  geom_errorbar(aes(ymin = p.025, ymax = p.975), width = 0.2)
  
head(mort_retrodictions)
surv_retrodictions_plot = 
  ggplot(
  data = mort_retrodictions,
  aes(x = r, y = surv_true))+
  geom_jitter(height = 0.1, width = 0, size = 1, color = 'red')+
  theme_minimal()+
  geom_point(data = mort_retrodictions,
             aes(x = r, y = p.mean))+
  geom_ribbon(data = mort_retrodictions,
              aes(x = r, ymin = p.025, ymax = p.975, y = p.50),
              alpha = 0.2)+
  geom_point(data = 
               mort_retrodictions %>% 
               group_by(r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             aes(x = r, y = surv_true),
             color = 'blue', pch = 4)+
  labs(x = 'Rank (of predicted survival probability)',
       y = 'Survival')

surv_retrodictions_plot

head(intermediate)

head(mort_retrodictions)

mort_retrodictions %>%
  filter(r_bin == '(538,605]') %>%
  left_join(intermediate %>%
              select(tree_id,
                     dbh_m.init, dbh_m2.init, fire, wpbr, ba_scaled, cwd_dep90_scaled, cwd_mean_scaled)) %>%
  group_by(tree_id, dbh_m.init, surv_true) %>%
  summarise() %>%
  ungroup() %>%
  ggplot(aes(x = dbh_m.init, y = surv_true))+
  geom_jitter(width = 0, height = 0.1)+
  geom_smooth(method = 'glm',
              method.arg = list(family = 'binomial'),
              formula = y~x+I(x**2))

mort_retrodictions %>%
  filter(r_bin == '(538,605]') %>%
  left_join(intermediate %>%
              select(tree_id,
                     dbh_m.init, dbh_m2.init, fire, wpbr, ba_scaled, cwd_dep90_scaled, cwd_mean_scaled))
  


# looks pretty good
ggsave(surv_retrodictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'retrodictions_s.png'),
       height = 4.5, width = 6, units = 'in')

#### growth ####################################################################


fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'growth_fit.rds'))

samples = as_draws_df(fitted_model$draws())

pila_training = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_growth_training.rds'))



samples.beta_g = samples %>% select(contains('beta'))
samples.ecoEffect_g = samples %>% select(contains('effect_ecosub')) %>% as.data.frame()
samples.plotEffect_g = samples %>% select(contains('effect_plot')) %>% as.data.frame()
growth_retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_g = as.numeric(samples.beta_g[i,])
                   mu = 
                     as.numeric(pila_training$X %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_training$ecosub_id]+
                     as.numeric(samples.plotEffect_g[i,])[pila_training$plot_id]
                   size1_sim = truncnorm::rtruncnorm(n = pila_training$N,
                                                     a = 0.0254,
                                                     mean = mu,
                                                     sd = samples$sigma_epsilon[i])
                   
                   result = 
                     pila_training$X %>% 
                     as_tibble(x = .) %>%
                     rownames_to_column('tree_id')
                   result$mu = mu
                   result$size1_sim = size1_sim
                   result$size1_true = pila_training$size1
                   result$iter = i
                   return(result)
                 }))

growth_retrodictions_plot = 
  ggplot(data = 
         growth_retrodictions %>%
         group_by(tree_id,size1_true) %>%
         summarise(size1_sim.50 = quantile(size1_sim, 0.5),
                   size1_sim.025 = quantile(size1_sim, 0.025),
                   size1_sim.975 = quantile(size1_sim, 0.975)) %>%
         ungroup() %>%
         mutate(r = dense_rank(size1_sim.50)))+
  geom_point(aes(x = r, y = size1_true),size = 0)+
  theme_minimal()+
  geom_line(aes(x = r, y = size1_sim.50), color = 'blue')+
  geom_ribbon(aes(x = r, ymin = size1_sim.025, ymax = size1_sim.975, y = size1_sim.50),
              alpha = 0.2)+
  labs(x = 'Rank (simulated size)', y = 'Size at remeasure')

growth_retrodictions_plot

# model is slightly overpredicting size of the smallest trees, missing some variation
# in size

growth_retrodictions_plot2 = 
  ggplot(data = 
         growth_retrodictions %>%
         group_by(tree_id,size1_true) %>%
         summarise(size1_sim.50 = quantile(size1_sim, 0.5),
                   size1_sim.025 = quantile(size1_sim, 0.025),
                   size1_sim.975 = quantile(size1_sim, 0.975)) %>%
         ungroup() %>%
         mutate(r = dense_rank(size1_sim.50)))+
  geom_point(aes(x = size1_sim.50, y = size1_true),size = 0)+
  theme_minimal()+
  geom_ribbon(aes(x = size1_sim.50, ymin = size1_sim.025, ymax = size1_sim.975, y = size1_sim.50),
              alpha = 0.2)+
  geom_abline(intercept = 0, slope = 1, color = 'blue')+
  labs(x = 'Median posterior predicted DBH (m)',
       y = 'True DBH (m)')

growth_retrodictions_plot2

growth_retrodictions_plot2+
  coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.25))+
  scale_y_continuous(breaks = seq(from = 0, to = 0.254, by = 0.0254))+
  scale_x_continuous(breaks = seq(from = 0, to = 2.54, by = 0.0254))

ggsave(growth_retrodictions_plot2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'retrodictions_g.png'),
       height = 4.5, width = 6, units = 'in')

#### recruitment ###############################################################

fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'fecd_fit.rds'))

samples = as_draws_df(fitted_model$draws())

pila_training = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_fecd_training.rds'))

samples.beta_f = samples %>% select(contains('beta'))

recruitment.retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_f = as.numeric(samples.beta_f[i,])
                   
                   # to get mean, .05, and .95th fecundities, run this line by line
                   # with beta_f = mean(samples.beta_f), beta_f = quantile(samples.beta_f, 0.05), etc.
                   
                   logf = 
                     as.numeric(pila_training$X %*% beta_f)
                   
                   f = 
                     matrix(nrow = pila_training$M, ncol = pila_training$P, 
                            byrow = FALSE,
                            data = 
                              sapply(X = 1:pila_training$P,
                                     FUN = function(plot){
                                       sapply(X = 1:pila_training$M,
                                              FUN = function(sizeclass_from){
                                                exp(logf[sizeclass_from+(pila_training$M*(plot-1))])
                                              })
                                     }))
                   nprime = 
                     sapply(X = 1:pila_training$P,
                            FUN = function(plot){
                              sum(f[,plot])
                            })
                   
                   cprime_pred = 
                     rnbinom(n = length(nprime),
                             mu = nprime*pila_training$a,
                             size = samples$kappa[i])
                   
                   result = 
                     data.frame(plot = 1:pila_training$P) %>%
                     as_tibble(x = .)
                   
                   result$density_pred = 
                     nprime %>%
                     as.numeric()
                   result$count_sim = 
                     cprime_pred
                   result$count_true = 
                     pila_training$cprime
                   result$iter = i
                   return(result)
                 }))

head(recruitment.retrodictions)

recr_retrodictions_plot = 
  recruitment.retrodictions %>%
  pivot_longer(cols = c('count_sim', 'count_true'),
               names_to = 'source', values_to = 'count') %>%
  ggplot(aes(x = count, fill = source))+
  geom_bar(position = position_dodge())+
  scale_x_continuous(limits = c(-1, 10))+
  theme_minimal()+
  labs(x = 'Count of new recreuits', y = 'Count of simulated plots')

recr_retrodictions_plot
# looks good

recr_retrodictions_plot2 = 
  recruitment.retrodictions %>%
  group_by(plot, count_true) %>%
  summarise(density_pred.50 = quantile(density_pred,probs = 0.5),
            count_sim.975 = quantile(count_sim, probs = 0.975),
            count_sim.025 = quantile(count_sim, probs = 0.025),
            count_sim.50 = quantile(count_sim, probs = 0.5)) %>%
  ungroup() %>%
  ggplot(aes(x = density_pred.50, y = count_sim.50))+
  geom_abline(intercept = 0, slope = 1, color = 'blue')+
  geom_errorbar(aes(ymin = count_sim.025, ymax = count_sim.975),
              alpha = 0.2)+
  geom_point(aes(y = count_true))+
  theme_minimal()


recr_retrodictions_plot2

# looks good
ggsave(recr_retrodictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'retrodictions_r.png'))

