library(here)
library(tidyverse)

fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'real_fits',
                                  'pila.rds'))

pila_validation = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_validation.rds'))

#### survival ##################################################################

samples = as_draws_df(fitted_model$draws())

# simulating one data set per draw

samples.beta_s = samples %>% select(contains('beta_s'))
samples.ecoEffect_s = samples %>% select(contains('ecoEffect_s')) %>% as.data.frame()
mort_predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   logitp = 
                     as.numeric(pila_validation$X_s %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_validation$ecosub_s]
                   p = boot::inv.logit(logitp)
                   surv_sim = rbinom(n = pila_validation$N_s, 
                                 size = 1, 
                                 prob = p)
                   
                   result = 
                     pila_validation$X_s %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$logitp = logitp
                   result$p = p
                   result$surv_sim = surv_sim
                   result$surv_true = pila_validation$surv
                   result$iter = i
                   return(result)
                 }))  %>%
  group_by(tree_id, surv_true) %>%
  summarise(p.mid = mean(p),
            surv_sim = mean(surv_sim)) %>%
  ungroup() %>%
  mutate(r = dense_rank(p.mid),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 10)))

ggplot(
  data = mort_predictions,
  aes(x = r, y = surv_true))+
  geom_jitter(height = 0.1, width = 0, size = 0, color = 'red')+
  theme_minimal()+
  geom_point(data = mort_predictions %>%
               group_by(r_bin) %>%
               summarise(r = mean(r),
                         surv_sim = sum(surv_sim)/n()) %>%
               ungroup(),
             aes(x = r, y = surv_sim))+
  geom_ribbon(data = mort_predictions %>%
                group_by(r_bin) %>%
                summarise(r = mean(r),
                          q50 = qbinom(p = 0.5,
                                       size = n(),
                                       prob = mean(p.mid))/n(),
                          q05 = qbinom(p = 0.025,
                                       size = n(),
                                       prob = mean(p.mid))/n(),
                          q95 = qbinom(p = 0.95,
                                       size = n(),
                                       prob = mean(p.mid))/n()) %>%
                ungroup(),
              aes(x = r, ymin = q05, ymax = q95, y = q50),
              alpha = 0.2)+
  geom_point(data = 
               mort_predictions %>% 
               group_by(r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             aes(x = r, y = surv_true),
             color = 'blue', pch = 4)

# underestimating the variability in survival

#### growth ####################################################################

samples.beta_g = samples %>% select(contains('beta_g'))
samples.ecoEffect_g = samples %>% select(contains('ecoEffect_g')) %>% as.data.frame()
growth_predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_g = as.numeric(samples.beta_g[i,])
                   mu = 
                     as.numeric(pila_validation$X_g %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_validation$ecosub_g]
                   size1_sim = truncnorm::rtruncnorm(n = pila_validation$N_g,
                                                     a = 0,
                                                     mean = mu,
                                                     sd = samples$sigmaEpsilon_g[i])
                   
                   result = 
                     pila_validation$X_g %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$mu = mu
                   result$size1_sim = size1_sim
                   result$size1_true = pila_validation$size1_g
                   result$iter = i
                   return(result)
                 }))

ggplot(data = 
         growth_predictions %>%
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
         growth_predictions %>%
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

# looks great