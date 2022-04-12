library(here)
library(tidyverse)
library(posterior)
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
samples.plotEffect_s = samples %>% select(contains('plotEffect_s')) %>% as.data.frame()
mort_predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   logitp = 
                     as.numeric(pila_validation$X_s %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_validation$ecosub_s]+
                     as.numeric(samples.plotEffect_s[i,])[pila_validation$plotid_s]
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
  summarise(p.50 = quantile(p, 0.5),
            surv_sim = mean(surv_sim),
            p.mean = mean(p),
            p.025 = quantile(p, 0.025),
            p.975 = quantile(p, 0.975)) %>%
  ungroup() %>%
  mutate(r = dense_rank(p.mean),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 10)))

head(mort_predictions)
surv_predictons_plot = 
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

surv_predictons_plot

ggsave(surv_predictons_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'predictions_s.png'))

#### growth ####################################################################

samples.beta_g = samples %>% select(contains('beta_g'))
samples.ecoEffect_g = samples %>% select(contains('ecoEffect_g')) %>% as.data.frame()
samples.plotEffect_g = samples %>% select(contains('plotEffect_g')) %>% as.data.frame()
growth_predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_g = as.numeric(samples.beta_g[i,])
                   mu = 
                     as.numeric(pila_validation$X_g %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_validation$ecosub_g]+
                     as.numeric(samples.plotEffect_g[i,])[pila_validation$plotid_g]
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
  geom_point(aes(x = size1_sim.50, y = size1_true),size = 1)+
  theme_minimal()+
  geom_ribbon(aes(x = size1_sim.50, ymin = size1_sim.025, ymax = size1_sim.975, y = size1_sim.50),
              alpha = 0.2)+
  geom_abline(intercept = 0, slope = 1, color = 'blue')

growth_predictions_plot2

# looks great
ggsave(growth_predictions_plot2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'predictions_g.png'))
#### recruitment ###############################################################


samples.beta_f = samples %>% select(contains('beta_f'))
samples.ecoEffect_f = samples %>% select(contains('ecoEFfect_f'))
samples.plotEffect_f = samples %>% select(contains('plotEffect_f'))

recruitment.predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_f = as.numeric(samples.beta_f[i,])
                   beta_g = as.numeric(samples.beta_g[i,])
                   beta_s = as.numeric(samples.beta_s[i,])
                   
                   logf = 
                     as.numeric(pila_validation$X_r %*% beta_f)+
                     as.numeric(samples.ecoEffect_f[i,])[pila_validation$ecosub_r]+
                     as.numeric(samples.plotEffect_f[i,])[pila_validation$plotid_r]
                   
                   mu = 
                     as.numeric(pila_validation$X_r %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_validation$ecosub_r]+
                     as.numeric(samples.plotEffect_g[i,])[pila_validation$plotid_r]
                   
                   logitp = 
                     as.numeric(pila_validation$X_r %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_validation$ecosub_r]+
                     as.numeric(samples.plotEffect_s[i,])[pila_validation$plotid_r]
                   
                   s = 
                     matrix(nrow = 1, ncol = pila_validation$S_r, byrow = TRUE,
                            data = 
                              sapply(X = 1:pila_validation$S_r,
                                     FUN = function(subplot){
                                       boot::inv.logit(logitp[1+(pila_validation$M_r*(subplot-1))])
                                     }))
                   
                   g = 
                     matrix(nrow = 2, ncol = pila_validation$S_r, byrow = FALSE,
                            data = 
                              sapply(X = 1:pila_validation$S_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:2,
                                              FUN = function(sizeclass_to){
                                                (pnorm(pila_validation$u_bounds[sizeclass_to],
                                                      mean = mu[1+(pila_validation$M_r*(subplot-1))],
                                                      sd = samples$sigmaEpsilon_g[i])-
                                                   pnorm(pila_validation$l_bounds[sizeclass_to],
                                                         mean = mu[1+(pila_validation$M_r*(subplot-1))],
                                                         sd = samples$sigmaEpsilon_g[i])) / 
                                                  (1-pnorm(0,
                                                           mean = mu[1+(pila_validation$M_r*(subplot-1))],
                                                           sd = samples$sigmaEpsilon_g[i]))
                                              })
                                     }))
                   
                   # growth kernel from sizeclass 1
                   growKern = 
                     matrix(nrow = pila_validation$S_r, ncol = 2, byrow = TRUE,
                            data = 
                              sapply(X = 1:pila_validation$S_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:2,
                                              FUN = function(sizeclass_to){
                                                s[1,subplot]*g[sizeclass_to,subplot]
                                              })
                                     }))
                   
                   f = 
                     matrix(nrow = pila_validation$M_r, ncol = pila_validation$S_r, 
                            byrow = FALSE,
                            data = 
                              sapply(X = 1:pila_validation$S_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:pila_validation$M_r,
                                              FUN = function(sizeclass_from){
                                                exp(logf[sizeclass_from+(pila_validation$M_r*(subplot-1))])
                                              })
                                     }))
                   
                   recKern = 
                     array(dim = 
                             c(pila_validation$S_r,
                               2,
                               pila_validation$M_r),
                           dimnames = 
                             list(subplot = 1:pila_validation$S_r,
                                  sizeclass_to = 1:2,
                                  sizeclass_from = 1:pila_validation$M_r),
                           data = 
                             sapply(X = 1:pila_validation$M_r,
                                    FUN = function(sizeclass_from){
                                      sapply(X = 1:2,
                                             FUN = function(sizeclass_to){
                                               sapply(X = 1:pila_validation$S_r,
                                                      FUN = function(subplot){
                                                        pila_validation$r[sizeclass_to]*
                                                          f[sizeclass_from,subplot]
                                                      })
                                             })
                                    }))
                   
                   A = 
                     array(dim = 
                             c(pila_validation$S_r,
                               2,
                               pila_validation$M_r),
                           dimnames = 
                             list(subplot = 1:pila_validation$S_r,
                                  sizeclass_to = 1:2,
                                  sizeclass_from = 1:pila_validation$M_r),
                           data = 
                             sapply(X = 1:pila_validation$M_r,
                                    FUN = function(sizeclass_from){
                                      sapply(X = 1:2,
                                             FUN = function(sizeclass_to){
                                               sapply(X = 1:pila_validation$S_r,
                                                      FUN = function(subplot){
                                                        if(sizeclass_from == 1){
                                                          growKern[subplot, sizeclass_to]+
                                                            recKern[subplot, sizeclass_to, sizeclass_from]
                                                        } else {
                                                          recKern[subplot,sizeclass_to,sizeclass_from]
                                                        }
                                                      })
                                             })
                                    }))
                   
                   nprime = 
                     matrix(nrow = pila_validation$S_r, ncol = 2, byrow = TRUE,
                            data = 
                              sapply(X = 1:pila_validation$S_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:2,
                                              FUN = function(sizeclass){
                                                as.numeric(A[subplot,sizeclass,] %*% 
                                                             pila_validation$n[subplot,])
                                              })
                                     }))
                   
                   cprime_pred = 
                     matrix(nrow = 2, ncol = pila_validation$S_r, byrow = FALSE,
                            data = 
                              sapply(X = 1:pila_validation$S_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:2,
                                              FUN = function(sizeclass){
                                                rnbinom(n = 1,
                                                        mu = 
                                                          nprime[subplot,sizeclass]*
                                                          pila_validation$a[sizeclass],
                                                        size = samples$kappa_r[i])
                                              })
                                     }))
                   
                   result = 
                     expand.grid(subplot = 1:pila_validation$S_r,
                                 sizeclass = 1:2) %>%
                     as_tibble() %>%
                     arrange(subplot,sizeclass)
                   
                   result$density_pred = 
                     sapply(X = 1:pila_validation$S_r,
                            FUN = function(subplot){
                              sapply(X = 1:2,
                                     FUN = function(sizeclass){
                                       nprime[subplot,sizeclass]*pila_validation$a[sizeclass]
                                     })
                            }) %>%
                     as.numeric()
                   result$count_sim = 
                     sapply(X = 1:pila_validation$S_r,
                            FUN = function(subplot){
                              sapply(X = 1:2,
                                     FUN = function(sizeclass){
                                       cprime_pred[sizeclass,subplot]
                                     })
                            }) %>%
                     as.integer()
                   result$count_true = 
                     sapply(X = 1:pila_validation$S_r,
                            FUN = function(subplot){
                              sapply(X = 1:2,
                                     FUN = function(sizeclass){
                                       pila_validation$cprime[sizeclass,subplot]
                                     })
                            }) %>%
                     as.integer()
                   result$iter = i
                   return(result)
                 }))

head(recruitment.predictions)

recruitment.predictions %>%
  pivot_longer(cols = c('count_sim', 'count_true'),
               names_to = 'source', values_to = 'count') %>%
  ggplot(aes(x = count, fill = source))+
  geom_bar(position = position_dodge())+
  scale_x_continuous(limits = c(-1, 10))+
  theme_minimal()

# looks good

recruitment_predictions_plot = 
  recruitment.predictions %>%
  group_by(subplot,sizeclass, count_true) %>%
  summarise(density_pred.50 = quantile(density_pred,probs = 0.5),
            count_sim.975 = quantile(count_sim, probs = 0.975),
            count_sim.025 = quantile(count_sim, probs = 0.025),
            count_sim.50 = quantile(count_sim, probs = 0.5)) %>%
  ungroup() %>%
  ggplot(aes(x = density_pred.50, y = count_true))+
  geom_abline(intercept = 0, slope = 1, color = 'blue')+
  geom_ribbon(aes(ymin = count_sim.025, ymax = count_sim.975),
              alpha = 0.2)+
  geom_point(aes(y = count_true))+
  theme_minimal()

recruitment_predictions_plot
# looks good
ggsave(recruitment_predictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'predictions_r.png'))

