library(here)
library(tidyverse)

fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'real_fits',
                                  'pila.rds'))

pila_training = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_training.rds'))

#### survival ##################################################################

samples = as_draws_df(fitted_model$draws())

# simulating one data set per draw

samples.beta_s = samples %>% select(contains('beta_s'))
samples.ecoEffect_s = samples %>% select(contains('ecoEffect_s')) %>% as.data.frame()
samples.plotEffect_s = samples %>% select(contains('plotEffect_s')) %>% as.data.frame()

mort_retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   logitp = 
                     as.numeric(pila_training$X_s %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_training$ecosub_s]+
                     as.numeric(samples.plotEffect_s[i,])[pila_training$plotid_s]
                   p = boot::inv.logit(logitp)
                   surv_sim = rbinom(n = pila_training$N_s, 
                                 size = 1, 
                                 prob = p)
                   
                   result = 
                     pila_training$X_s %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$logitp = logitp
                   result$p = p
                   result$surv_sim = surv_sim
                   result$surv_true = pila_training$surv
                   result$iter = i
                   return(result)
                 })) %>%
  group_by(tree_id, surv_true) %>%
  summarise(p.50 = quantile(p, 0.5),
            surv_sim = mean(surv_sim),
            p.025 = quantile(p, 0.025),
            p.975 = quantile(p, 0.975)) %>%
  ungroup() %>%
  mutate(r = dense_rank(p.50),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20)))

head(mort_retrodictions)
ggplot(
  data = mort_retrodictions,
  aes(x = r, y = surv_true))+
  geom_jitter(height = 0.1, width = 0, size = 0, color = 'red')+
  theme_minimal()+
  geom_point(data = mort_retrodictions %>%
               group_by(r_bin) %>%
               summarise(r = mean(r),
                         p.50 = mean(p.50)) %>%
               ungroup(),
             aes(x = r, y = p.50))+
  geom_ribbon(data = mort_retrodictions %>%
                group_by(r_bin) %>%
                summarise(r = mean(r),
                          p.50 = mean(p.50),
                          p.025 = mean(p.025),
                          p.975 = mean(p.975)) %>%
                ungroup(),
              aes(x = r, ymin = p.025, ymax = p.975, y = p.50),
              alpha = 0.2)+
  geom_point(data = 
               mort_retrodictions %>% 
               group_by(r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             aes(x = r, y = surv_true),
             color = 'blue', pch = 4)


# looks pretty good


#### growth ####################################################################

samples.beta_g = samples %>% select(contains('beta_g'))
samples.ecoEffect_g = samples %>% select(contains('ecoEffect_g')) %>% as.data.frame()
samples.plotEffect_g = samples %>% select(contains('plotEffect_g')) %>% as.data.frame()
growth_retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_g = as.numeric(samples.beta_g[i,])
                   mu = 
                     as.numeric(pila_training$X_g %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_training$ecosub_g]+
                     as.numeric(samples.plotEffect_g[i,])[pila_training$plotid_g]
                   size1_sim = truncnorm::rtruncnorm(n = pila_training$N_g,
                                                     a = 0,
                                                     mean = mu,
                                                     sd = samples$sigmaEpsilon_g[i])
                   
                   result = 
                     pila_training$X_g %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$mu = mu
                   result$size1_sim = size1_sim
                   result$size1_true = pila_training$size1_g
                   result$iter = i
                   return(result)
                 }))

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
  labs(x = 'rank (simulated size)', y = 'Size at remeasure')

# model is slightly overpredicting size of the smallest trees, missing some variation
# in size

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
  geom_abline(intercept = 0, slope = 1, color = 'blue')

