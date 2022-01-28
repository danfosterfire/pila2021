library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'real_fits',
                                  'usgs',
                                  'pila.rds'))

pila_data = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'usgs',
                                   'pila_data.rds'))

#### survival ##################################################################

samples = as_draws_df(fitted_model$draws())

# simulating one data set per draw

samples.beta_s = samples %>% select(contains('beta_s'))
samples.plotEffect_s = samples %>% select(contains('plotEffect_s')) %>% as.data.frame()

mort_retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   logitp = 
                     as.numeric(pila_data$X_s %*% beta_s)+
                     as.numeric(samples.plotEffect_s[i,])[pila_data$plotid_s]
                   p = boot::inv.logit(logitp)
                   surv_sim = rbinom(n = pila_data$N_s, 
                                 size = 1, 
                                 prob = p)
                   
                   result = 
                     pila_data$X_s %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$logitp = logitp
                   result$p = p
                   result$surv_sim = surv_sim
                   result$surv_true = pila_data$surv_s
                   result$iter = i
                   return(result)
                 })) %>%
  group_by(tree_id, surv_true) %>%
  summarise(p.50 = quantile(p, 0.5),
            p.mean = mean(p),
            p.025 = quantile(p, 0.025),
            p.975 = quantile(p, 0.975)) %>%
  ungroup() %>%
  mutate(r = dense_rank(p.mean),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20)))

head(mort_retrodictions)
surv_retrodictions_plot = 
  ggplot(
  data = mort_retrodictions,
  aes(x = r, y = surv_true))+
  geom_jitter(height = 0.1, width = 0, size = 0, color = 'red')+
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
             color = 'blue', pch = 4)

surv_retrodictions_plot
# this one looks amazing... why does the FIA one look so bad by comparison?
ggsave(surv_retrodictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'report',
                             'usgs_retrodictions_s.png'))

#### growth ####################################################################

samples.beta_g = samples %>% select(contains('beta_g'))
samples.plotEffect_g = samples %>% select(contains('plotEffect_g')) %>% as.data.frame()
growth_retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_g = as.numeric(samples.beta_g[i,])
                   mu = 
                     as.numeric(pila_data$X_g %*% beta_g)+
                     as.numeric(samples.plotEffect_g[i,])[pila_data$plotid_g]
                   size1_sim = truncnorm::rtruncnorm(n = pila_data$N_g,
                                                     a = 0,
                                                     mean = mu,
                                                     sd = samples$sigmaEpsilon_g[i])
                   
                   result = 
                     pila_data$X_g %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$mu = mu
                   result$size1_sim = size1_sim
                   result$size1_true = pila_data$dbhT_g
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
  labs(x = 'rank (simulated size)', y = 'Size at remeasure')

growth_retrodictions_plot

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
  geom_abline(intercept = 0, slope = 1, color = 'blue')

growth_retrodictions_plot2 # looks good

ggsave(growth_retrodictions_plot2,
       filename = here::here('04-communication',
                             'figures',
                             'report',
                             'ugss_retrodictions_g.png'))

#### fecundity #################################################################


samples.beta_f = samples %>% select(contains('beta_f'))
samples.plotEffect_f = samples %>% select(contains('plotEffect_f'))

recruitment.retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_f = as.numeric(samples.beta_f[i,])
                   
                   logf = 
                     as.numeric(pila_data$X_f %*% beta_f)+
                     as.numeric(samples.plotEffect_f[i,])[pila_data$plotid_f]
                   
                   A = 
                     matrix(nrow = pila_data$C_f,
                            ncol = pila_data$M_f,
                            byrow = TRUE,
                            data = sapply(X = 1:pila_data$C_f,
                                          FUN = function(plot){
                                            sapply(X = 1:pila_data$M_f,
                                                   FUN = function(size_from){
                                                     exp(logf[size_from+(pila_data$M_f*(plot-1))])
                                                   })
                                          }))
                   
                   nprime = 
                     sapply(X = 1:pila_data$C_f,
                            FUN = function(plot){
                              A[plot,] %*% pila_data$n[,plot]
                            })
                   
                   
                   density_pred = 
                     nprime * pila_data$a
                     
                   count_sim = rnbinom(n = length(density_pred),
                                         mu = density_pred,
                                         size = samples$kappa_f[i]) %>%
                     as.integer()
                   
                   result = 
                     data.frame(census_id = 1:pila_data$C_f,
                                density_pred = density_pred,
                                count_sim = count_sim,
                                count_true = as.integer(pila_data$cprime_f)) %>%
                     as_tibble()
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
  theme_minimal()

recr_retrodictions_plot

# looks good
recr_retrodictions_plot2 = 
  recruitment.retrodictions %>%
  group_by(census_id, count_true) %>%
  summarise(density_pred.50 = quantile(density_pred,probs = 0.5),
            count_sim.975 = quantile(count_sim, probs = 0.975),
            count_sim.025 = quantile(count_sim, probs = 0.025),
            count_sim.50 = quantile(count_sim, probs = 0.5)) %>%
  ungroup() %>%
  ggplot(aes(x = density_pred.50, y = count_sim.50))+
  geom_abline(intercept = 0, slope = 1, color = 'blue')+
  geom_ribbon(aes(ymin = count_sim.025, ymax = count_sim.975),
              alpha = 0.2)+
  geom_point(aes(y = count_true))+
  theme_minimal()

recr_retrodictions_plot2

# looks good
ggsave(recr_retrodictions_plot2,
       filename = here::here('04-communication',
                             'figures',
                             'report',
                             'usgs_retrodictions_r.png'))


#### recruit size ##############################################################


recruit_size.retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(i){
                   
                   mu = samples$mu_r[i]
                   sigma = samples$sigmaEpsilon_r[i]
                   
                   sim_logsize = rnorm(n = length(pila_data$logdbh_r),
                                       mean = mu, sd = sigma)
                   
                   results = 
                     data.frame(sim_logsize = sim_logsize,
                                true_logsize = pila_data$logdbh_r,
                                iter = i)
                   
                   return(results)
                 }))

ggplot(data = 
         recruit_size.retrodictions %>%
         pivot_longer(cols = c('sim_logsize', 'true_logsize'),
                      names_to = 'source', values_to = 'value'),
       aes(x = value, fill = source))+
  geom_histogram(position = position_dodge())

recrsize_retrodictions_plot = 
  ggplot(data = 
         recruit_size.retrodictions %>%
         pivot_longer(cols = c('sim_logsize', 'true_logsize'),
                      names_to = 'source', values_to = 'value'),
       aes(x = exp(value), fill = source))+
  geom_histogram(position = position_dodge())+
  theme_minimal()

ggsave(recrsize_retrodictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'report',
                             'usgs_retrodictions_rsize.png'))
