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

pila_validation = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_mort_validation.rds'))

# simulating one data set per draw

samples.beta_s = samples %>% select(contains('beta'))
samples.ecoEffect_s = samples %>% select(contains('effect_ecosub')) %>% as.data.frame()
samples.plotEffect_s = samples %>% select(contains('effect_plot')) %>% as.data.frame()


mort_predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   logitp = 
                     as.numeric(pila_validation$X %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_validation$ecosub_id]+
                     as.numeric(samples.plotEffect_s[i,])[pila_validation$plot_id]
                   p = boot::inv.logit(logitp)
                   surv_sim = rbinom(n = pila_validation$N, 
                                     size = 1, 
                                     prob = p)
                   
                   result = 
                     pila_validation$X %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$logitp = logitp
                   result$p = p
                   result$surv_sim = surv_sim
                   result$surv_true = pila_validation$surv
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

do.call('bind_rows',
        lapply(X = 1:nrow(samples),
               FUN = function(i){
                 beta_s = as.numeric(samples.beta_s[i,])
                 logitp = 
                   as.numeric(pila_validation$X %*% beta_s)+
                   as.numeric(samples.ecoEffect_s[i,])[pila_validation$ecosub_id]+
                   as.numeric(samples.plotEffect_s[i,])[pila_validation$plot_id]
                 p = boot::inv.logit(logitp)
                 surv_sim = rbinom(n = pila_validation$N, 
                                   size = 1, 
                                   prob = p)
                 
                 result = 
                   pila_validation$X %>% 
                   as_tibble() %>%
                   rownames_to_column('tree_id')
                 result$logitp = logitp
                 result$p = p
                 result$surv_sim = surv_sim
                 result$surv_true = pila_validation$surv
                 result$iter = i
                 return(result)
               })) %>%
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

head(mort_predictions)
surv_predictions_plot = 
  ggplot(
    data = mort_predictions,
    aes(x = r, y = surv_true))+
  geom_jitter(height = 0.1, width = 0, size = 1, color = 'red')+
  theme_minimal()+
  geom_point(data = mort_predictions,
             aes(x = r, y = p.mean))+
  geom_ribbon(data = mort_predictions,
              aes(x = r, ymin = p.025, ymax = p.975, y = p.50),
              alpha = 0.2)+
  geom_point(data = 
               mort_predictions %>% 
               group_by(r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             aes(x = r, y = surv_true),
             color = 'blue', pch = 4)


surv_predictions_plot
# looks pretty good
ggsave(surv_predictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'predictions_s.png'),
       height = 4.5, width = 6, units = 'in')

#### growth ####################################################################


fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'growth_fit.rds'))

samples = as_draws_df(fitted_model$draws())

pila_validation = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_growth_validation.rds'))



samples.beta_g = samples %>% select(contains('beta'))
samples.ecoEffect_g = samples %>% select(contains('effect_ecosub')) %>% as.data.frame()
samples.plotEffect_g = samples %>% select(contains('effect_plot')) %>% as.data.frame()
growth_predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_g = as.numeric(samples.beta_g[i,])
                   mu = 
                     as.numeric(pila_validation$X %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_validation$ecosub_id]+
                     as.numeric(samples.plotEffect_g[i,])[pila_validation$plot_id]
                   size1_sim = truncnorm::rtruncnorm(n = pila_validation$N,
                                                     a = 0.0254,
                                                     mean = mu,
                                                     sd = samples$sigma_epsilon[i])
                   
                   result = 
                     pila_validation$X %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$mu = mu
                   result$size1_sim = size1_sim
                   result$size1_true = pila_validation$size1
                   result$iter = i
                   return(result)
                 }))

growth_predictions_plot = 
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

growth_predictions_plot

# model is slightly overpredicting size of the smallest trees, missing some variation
# in size

growth_predictions_plot2 = 
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

growth_predictions_plot2

growth_predictions_plot2+
  coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.25))+
  scale_y_continuous(breaks = seq(from = 0, to = 0.254, by = 0.0254))+
  scale_x_continuous(breaks = seq(from = 0, to = 2.54, by = 0.0254))

ggsave(growth_predictions_plot2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'predictions_g.png'),
       height = 4.5, width = 6, units = 'in')

#### recruitment ###############################################################

fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'fecd_fit.rds'))

samples = as_draws_df(fitted_model$draws())

pila_validation = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_fecd_validation.rds'))

samples.beta_f = samples %>% select(contains('beta'))

recruitment.predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_f = as.numeric(samples.beta_f[i,])
                   
                   logf = 
                     as.numeric(pila_validation$X %*% beta_f)
                   
                   f = 
                     matrix(nrow = pila_validation$M, ncol = pila_validation$P, 
                            byrow = FALSE,
                            data = 
                              sapply(X = 1:pila_validation$P,
                                     FUN = function(plot){
                                       sapply(X = 1:pila_validation$M,
                                              FUN = function(sizeclass_from){
                                                exp(logf[sizeclass_from+(pila_validation$M*(plot-1))])
                                              })
                                     }))
                   nprime = 
                     sapply(X = 1:pila_validation$P,
                            FUN = function(plot){
                              sum(f[,plot])
                            })
                   
                   cprime_pred = 
                     rnbinom(n = length(nprime),
                             mu = nprime*pila_validation$a,
                             size = samples$kappa[i])
                   
                   result = 
                     data.frame(plot = 1:pila_validation$P) %>%
                     as_tibble()
                   
                   result$density_pred = 
                     nprime %>%
                     as.numeric()
                   result$count_sim = 
                     cprime_pred
                   result$count_true = 
                     pila_validation$cprime
                   result$iter = i
                   return(result)
                 }))

head(recruitment.predictions)

recr_predictions_plot = 
  recruitment.predictions %>%
  pivot_longer(cols = c('count_sim', 'count_true'),
               names_to = 'source', values_to = 'count') %>%
  ggplot(aes(x = count, fill = source))+
  geom_bar(position = position_dodge())+
  #scale_x_continuous(limits = c(-1, 10))+
  theme_minimal()

recr_predictions_plot
# looks good

recr_predictions_plot2 = 
  recruitment.predictions %>%
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


recr_predictions_plot2

# looks good
ggsave(recr_predictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'predictions_r.png'))

