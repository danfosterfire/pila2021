library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

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


mort_retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   logitp = 
                     as.numeric(pila_training$X %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_training$ecosub_id]+
                     as.numeric(samples.plotEffect_s[i,])[pila_training$plot_id]
                   p = boot::inv.logit(logitp)
                   surv_sim = rbinom(n = pila_training$N, 
                                 size = 1, 
                                 prob = p)
                   
                   result = 
                     pila_training$X %>% 
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
                   as.numeric(pila_training$X %*% beta_s)+
                   as.numeric(samples.ecoEffect_s[i,])[pila_training$ecosub_id]+
                   as.numeric(samples.plotEffect_s[i,])[pila_training$plot_id]
                 p = boot::inv.logit(logitp)
                 surv_sim = rbinom(n = pila_training$N, 
                                   size = 1, 
                                   prob = p)
                 
                 result = 
                   pila_training$X %>% 
                   as_tibble() %>%
                   rownames_to_column('tree_id')
                 result$logitp = logitp
                 result$p = p
                 result$surv_sim = surv_sim
                 result$surv_true = pila_training$surv
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
# looks pretty good
ggsave(surv_retrodictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'retrodictions_s.png'))

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
                                                     a = 0,
                                                     mean = mu,
                                                     sd = samples$sigma_epsilon[i])
                   
                   result = 
                     pila_training$X %>% 
                     as_tibble() %>%
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
  labs(x = 'rank (simulated size)', y = 'Size at remeasure')

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
  geom_abline(intercept = 0, slope = 1, color = 'blue')

growth_retrodictions_plot2

growth_retrodictions_plot2+
  coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.25))+
  scale_y_continuous(breaks = seq(from = 0, to = 0.254, by = 0.0254))+
  scale_x_continuous(breaks = seq(from = 0, to = 2.54, by = 0.0254))

ggsave(growth_retrodictions_plot2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'retrodictions_g.png'))

#### recruitment ###############################################################


samples.beta_f = samples %>% select(contains('beta_f'))
samples.ecoEffect_f = samples %>% select(contains('ecoEFfect_f'))

recruitment.retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_f = as.numeric(samples.beta_f[i,])
                   beta_g = as.numeric(samples.beta_g[i,])
                   beta_s = as.numeric(samples.beta_s[i,])
                   
                   logf = 
                     as.numeric(pila_training$X_r %*% beta_f)+
                     as.numeric(samples.ecoEffect_f[i,])[pila_training$ecosub_r]
                   
                   mu = 
                     as.numeric(pila_training$X_rg %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_training$ecosub_r]+
                     as.numeric(samples.plotEffect_g[i,])[pila_training$plotid_r]
                   
                   logitp = 
                     as.numeric(pila_training$X_r %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_training$ecosub_r]+
                     as.numeric(samples.plotEffect_s[i,])[pila_training$plotid_r]
                   
                   s = 
                     matrix(nrow = 1, ncol = pila_training$P_r, byrow = TRUE,
                            data = 
                              sapply(X = 1:pila_training$P_r,
                                     FUN = function(subplot){
                                       boot::inv.logit(logitp[1+(pila_training$M_r*(subplot-1))])
                                     }))
                   
                   g = 
                     matrix(nrow = pila_training$max_recr_class, ncol = pila_training$P_r, byrow = FALSE,
                            data = 
                              sapply(X = 1:pila_training$P_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:pila_training$max_recr_class,
                                              FUN = function(sizeclass_to){
                                                (pnorm(pila_training$u_bounds[sizeclass_to],
                                                      mean = mu[1+(pila_training$M_r*(subplot-1))],
                                                      sd = samples$sigmaEpsilon_g[i])-
                                                   pnorm(pila_training$l_bounds[sizeclass_to],
                                                         mean = mu[1+(pila_training$M_r*(subplot-1))],
                                                         sd = samples$sigmaEpsilon_g[i])) / 
                                                  (1-pnorm(0,
                                                           mean = mu[1+(pila_training$M_r*(subplot-1))],
                                                           sd = samples$sigmaEpsilon_g[i]))
                                              })
                                     }))
                   
                   # growth kernel from sizeclass 1
                   growKern = 
                     matrix(nrow = pila_training$P_r, ncol = 2, byrow = TRUE,
                            data = 
                              sapply(X = 1:pila_training$P_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:pila_training$max_recr_class,
                                              FUN = function(sizeclass_to){
                                                s[1,subplot]*g[sizeclass_to,subplot]
                                              })
                                     }))
                   
                   f = 
                     matrix(nrow = pila_training$M_r, ncol = pila_training$P_r, 
                            byrow = FALSE,
                            data = 
                              sapply(X = 1:pila_training$P_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:pila_training$M_r,
                                              FUN = function(sizeclass_from){
                                                exp(logf[sizeclass_from+(pila_training$M_r*(subplot-1))])
                                              })
                                     }))
                   
                   recKern = 
                     array(dim = 
                             c(pila_training$P_r,
                               pila_training$max_recr_class,
                               pila_training$M_r),
                           dimnames = 
                             list(subplot = 1:pila_training$P_r,
                                  sizeclass_to = 1:pila_training$max_recr_class,
                                  sizeclass_from = 1:pila_training$M_r),
                           data = 
                             sapply(X = 1:pila_training$M_r,
                                    FUN = function(sizeclass_from){
                                      sapply(X = 1:pila_training$max_recr_class,
                                             FUN = function(sizeclass_to){
                                               sapply(X = 1:pila_training$P_r,
                                                      FUN = function(subplot){
                                                        pila_training$r[sizeclass_to]*
                                                          f[sizeclass_from,subplot]
                                                      })
                                             })
                                    }))
                   
                   A = 
                     array(dim = 
                             c(pila_training$P_r,
                               pila_training$max_recr_class,
                               pila_training$M_r),
                           dimnames = 
                             list(subplot = 1:pila_training$P_r,
                                  sizeclass_to = 1:pila_training$max_recr_class,
                                  sizeclass_from = 1:pila_training$M_r),
                           data = 
                             sapply(X = 1:pila_training$M_r,
                                    FUN = function(sizeclass_from){
                                      sapply(X = 1:pila_training$max_recr_class,
                                             FUN = function(sizeclass_to){
                                               sapply(X = 1:pila_training$P_r,
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
                     matrix(nrow = pila_training$P_r, ncol = pila_training$max_recr_class, byrow = TRUE,
                            data = 
                              sapply(X = 1:pila_training$P_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:pila_training$max_recr_class,
                                              FUN = function(sizeclass){
                                                as.numeric(A[subplot,sizeclass,] %*% 
                                                             pila_training$n[subplot,])
                                              })
                                     }))
                   
                   cprime_pred = 
                     matrix(nrow = pila_training$max_recr_class, ncol = pila_training$P_r, byrow = FALSE,
                            data = 
                              sapply(X = 1:pila_training$P_r,
                                     FUN = function(subplot){
                                       sapply(X = 1:pila_training$max_recr_class,
                                              FUN = function(sizeclass){
                                                rnbinom(n = 1,
                                                        mu = 
                                                          nprime[subplot,sizeclass]*
                                                          pila_training$a[sizeclass],
                                                        size = samples$kappa_r[i])
                                              })
                                     }))
                   
                   result = 
                     expand.grid(subplot = 1:pila_training$P_r,
                                 sizeclass = 1:pila_training$max_recr_class) %>%
                     as_tibble() %>%
                     arrange(subplot,sizeclass)
                   
                   result$density_pred = 
                     sapply(X = 1:pila_training$P_r,
                            FUN = function(subplot){
                              sapply(X = 1:pila_training$max_recr_class,
                                     FUN = function(sizeclass){
                                       nprime[subplot,sizeclass]*pila_training$a[sizeclass]
                                     })
                            }) %>%
                     as.numeric()
                   result$count_sim = 
                     sapply(X = 1:pila_training$P_r,
                            FUN = function(subplot){
                              sapply(X = 1:pila_training$max_recr_class,
                                     FUN = function(sizeclass){
                                       cprime_pred[sizeclass,subplot]
                                     })
                            }) %>%
                     as.integer()
                   result$count_true = 
                     sapply(X = 1:pila_training$P_r,
                            FUN = function(subplot){
                              sapply(X = 1:pila_training$max_recr_class,
                                     FUN = function(sizeclass){
                                       pila_training$cprime[sizeclass,subplot]
                                     })
                            }) %>%
                     as.integer()
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
  group_by(subplot,sizeclass, count_true) %>%
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
                             'manuscript',
                             'retrodictions_r.png'))

