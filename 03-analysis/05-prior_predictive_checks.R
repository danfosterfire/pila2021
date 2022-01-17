
# i'm doing this because my first crack at the model, with very weakly informative 
# priors matching those used in the shriver paper, did not perform well in 
# terms of parameter recovery. In the shriver paper, they only check a single 
# set of parameters. You are really supposed to check multiple "true" parameter 
# sets, and according to the lit I'm seeing on SBC these true values should be 
# drawn from the prior distribution(s) in order for the posterior credible 
# intervals to be interpretable as something like a CI. 
# With some sets of parameters, the model recovers them 
# well, but for many parameter sets at least some true values fall outside the 
# posterior credible intervals, and in some cases the sampler runs into 
# divergence trouble or is very slow to converge. I think the issue is that 
# with too-flexible priors we can get "True" parameter values that are 
# wildly unrealistic, and they are conflicting with the greater-than-0 constraint 
# in the truncated-normal growth model. In reality I don't think the 
# weakly informative prior distributions are plausible - in particular I expect 
# the true value of betaSize_g (the coefficient for initial measurement size 
# determining followup measurement size) to be near 1, and negative values of 
# XB_size (the linear predictor for size.remeasure from fixed effects) are 
# not plausible either. With weird paramter values the truncated normal, where 
# all the responses are positive, becomes difficult to fit. The reading I've done on this 
# (the stan user's guide and 
# https://betanalpha.github.io/assets/case_studies/principled_bayesian_workflow.html)
# suggests doing prior predictive checks: basically drawing true parameter values 
# from the prior, simulating responses using those parameter values, and plotting 
# the resulting distributions of responses. The priors should be calibrated to 
# give plausible responses (where plausible is defined using prior knowledge 
# external to the data I'll use to fit the posterior). So, eg., size at remeasurement 
# shouldn't include values greater than a few meters, and the relationship between 
# size at initial measurement and size at remeasurement should be plausible. 
# We want to avoid the bulk of the survival model priors falling in regions where 
# observed survival is 0 or 100%. That sort of thing. 

library(here)
library(tidyverse)

pila_data = readRDS(here::here('02-data', '02-for_analysis', 'pila_training.rds'))


#### function defs #############################################################

#### priors for growth model ###################################################

# pull parameters from prior
beta_g = matrix(nrow = pila_data$K, ncol = 1000, byrow = FALSE,
                data = sapply(X = 1:1000,
                              FUN = function(i){
                                rnorm(n = pila_data$K, mean = 0, sd = 0.1)
                              }))
# special priors for intercept and initial size
beta_g[1,] = rnorm(n = 1000, mean = 0.15, sd = 0.1)
beta_g[2,] = rnorm(n = 1000, mean = 1, sd = 0.25)

sigmaPlot_g = truncnorm::rtruncnorm(n = 1000, mean = 0, sd = 0.25, a = 0)
sigmaEco_g = truncnorm::rtruncnorm(n = 1000, mean = 0, sd = 0.25, a = 0)
sigmaEpsilon_g = truncnorm::rtruncnorm(n = 1000, mean = 0, sd = 0.25, a = 0)

# simulate responses
growth_sims = 
  do.call('bind_rows',
    lapply(X = 1:1000,
         FUN = function(i){
           XB =
             as.numeric(as.matrix(pila_data$X_g) %*% beta_g[,i])
           ecoEffects = rnorm(n = pila_data$E, mean = 0, sd = sigmaEco_g[i])
           plotEffects = rnorm(n = pila_data$P, mean = 0, sd = sigmaPlot_g[i])
           
           y = truncnorm::rtruncnorm(n = nrow(pila_data$X_g),
                                     mean = 
                                       XB+
                                       ecoEffects[pila_data$ecosub_g]+
                                       plotEffects[pila_data$plotid_g],
                                     sd = sigmaEpsilon_g[i],
                                     a = 0)
           results = 
             pila_data$X_g %>%
             as_tibble()
           results$beta1_g = beta_g[1,i]
           results$beta2_g = beta_g[2,i]
           results$beta3_g = beta_g[3,i]
           results$beta4_g = beta_g[4,i]
           results$beta5_g = beta_g[5,i]
           results$beta6_g = beta_g[6,i]
           results$beta7_g = beta_g[7,i]
           results$beta8_g = beta_g[8,i]
           results$beta9_g = beta_g[9,i]
           results$beta10_g = beta_g[10,i]
           results$beta11_g = beta_g[11,i]
           results$beta12_g = beta_g[12,i]
           
           results$sigmaEco_g = sigmaEco_g[i]
           results$sigmaPlot_g = sigmaPlot_g[i]
           results$sigmaEpsilon_g = sigmaEpsilon_g[i]
           results$size1_g = y
           results$sim = i
           return(results)
         }))

# plot results

## y distributions
ggplot(data = growth_sims,
       aes(x = size1_g))+
  geom_histogram()

ggplot(data = growth_sims,
       aes(x = size1_g))+
  geom_histogram()+
  geom_vline(xintercept = quantile(growth_sims$size1_g, probs = 0.5), color = 'blue')+
  geom_vline(xintercept = quantile(growth_sims$size1_g, probs = 0.99), color = 'red')

## y distributions vs parameter values

lapply(X = 1:pila_data$K,
       FUN = function(k){
         d = as.data.frame(growth_sims)
         d[,paste0('beta_bin')] = 
           cut(d[,paste0('beta',k,'_g')],
               breaks = seq(from = min(as.numeric(d[,paste0('beta',k,'_g')]))-0.01,
                            to = max(as.numeric(d[,paste0('beta', k, '_g')]))+0.01,
                            length.out = 10))
         ggplot(data = d,
                aes(y = size1_g, x = beta_bin))+
           geom_boxplot()+
           labs(x = paste0('beta', k))
       })



growth_sims %>%
  mutate(sigmaEpsilon_bin = 
           cut(sigmaEpsilon_g, 
               breaks = 
                 seq(from = min(sigmaEpsilon_g)-0.1, to = max(sigmaEpsilon_g)+0.1, length.out = 10))) %>%
  ggplot(aes(x = sigmaEpsilon_bin, y = size1_g))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = -90))

growth_sims %>%
  mutate(size0_bin = 
           cut(dbh_m.init,
               breaks = seq(from = min(dbh_m.init)-0.01, to = max(dbh_m.init)+0.01, length.out =10))) %>%
  ggplot(aes(x = size0_bin, y = size1_g))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = -90))

#### priors for survival model #################################################

# pull parameters from prior
# pull parameters from prior
beta_s = matrix(nrow = pila_data$K, ncol = 1000, byrow = FALSE,
                data = sapply(X = 1:1000,
                              FUN = function(i){
                                rnorm(n = pila_data$K, mean = 0, sd = 1)
                              }))

sigmaEco_s = truncnorm::rtruncnorm(n = 1000, mean = 0, sd = 1, a = 0)
sigmaPlot_s = truncnorm::rtruncnorm(n = 1000, mean = 0, sd = 1, a = 0)

# simulate responses
growth_sims = 
  do.call('bind_rows',
    lapply(X = 1:1000,
         FUN = function(i){
           XB = as.numeric(as.matrix(pila_data$X_s) %*% beta_s[,i])
           ecoEffects = rnorm(n = pila_data$E, mean = 0, sd = sigmaEco_s[i])
           plotEffects = rnorm(n = pila_data$P, mean = 0, sd = sigmaPlot_s[i])
           y = rbinom(n = nrow(pila_data$X_s),
                      size = 1,
                      prob = 
                        boot::inv.logit(XB+
                                          ecoEffects[pila_data$ecosub_s]+
                                          plotEffects(pila_data$plotid_s)))
           results = 
             pila_data$X_s %>%
             as_tibble()
           results$beta1_s = beta_s[1,i]
           results$beta2_s = beta_s[2,i]
           results$beta3_s = beta_s[3,i]
           results$beta4_s = beta_s[4,i]
           results$beta5_s = beta_s[5,i]
           results$beta6_s = beta_s[6,i]
           results$beta7_s = beta_s[7,i]
           results$beta8_s = beta_s[8,i]
           results$beta9_s = beta_s[9,i]
           results$beta10_s = beta_s[10,i]
           results$beta11_s = beta_s[11,i]
           results$beta12_s = beta_s[12,i]
           results$sigmaEco_s = sigmaEco_s[i]
           results$sigmaPlot_s = sigmaPlot_s[i]
           results$surv = y
           results$sim = i
           return(results)
         }))

# plot results
growth_sims %>%
  ggplot(aes(x = as.factor(surv)))+
  geom_bar()


lapply(X = 1:pila_data$K,
       FUN = function(k){
         d = as.data.frame(growth_sims)
         d[,paste0('beta_bin')] = 
           cut(d[,paste0('beta',k,'_s')],
               breaks = seq(from = min(as.numeric(d[,paste0('beta',k,'_s')]))-0.01,
                            to = max(as.numeric(d[,paste0('beta', k, '_s')]))+0.01,
                            length.out = 10))
         ggplot(data = d,
                aes(fill = as.factor(surv), x = beta_bin))+
           geom_bar(position = position_fill())+
           labs(x = paste0('beta', k))
       })


# only works for two parameters
growth_sims %>%
  mutate(dbh_bin = 
           cut(dbh_m.init,
               breaks = seq(from = min(dbh_m.init)-0.01, to = max(dbh_m.init)+0.01, length.out =10)),
         beta2_bin = 
           cut(beta2_s, 
               breaks = seq(from = min(beta2_s)-0.01, 
                            to = max(beta2_s)+0.01,
                            length.out = 10)),
         beta1_bin = 
           cut(beta1_s, 
               breaks = seq(from = min(beta1_s)-0.01, 
                            to = max(beta1_s)+0.01,
                            length.out = 10))) %>%
  ggplot(aes(x = dbh_bin, fill = as.factor(surv)))+
  geom_bar(position = position_fill())+
  theme(axis.text.x = element_text(angle = -90))+
  facet_grid(beta1_bin~beta2_bin)

#### priors for fecundity model ################################################

# pull parameters from prior

