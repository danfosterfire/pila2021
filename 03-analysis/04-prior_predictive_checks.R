
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

pila_data = readRDS(here::here('02-data', '02-for_analysis', 'pila_data.rds'))


#### function defs #############################################################

#### priors for growth model ###################################################

# pull parameters from prior

beta0_g = rnorm(n = 1000, mean = 0, sd = 0.25)
betaSize_g = rnorm(n = 1000, mean = 1, sd = 0.25)

sigmaEpsilon_g = truncnorm::rtruncnorm(n = 1000, mean = 0, sd = 0.1, a = 0)

# simulate responses
growth_sims = 
  do.call('bind_rows',
    lapply(X = 1:1000,
         FUN = function(i){
           XB = 
             pila_data$X_g[,'intercept']*beta0_g[i] + 
             pila_data$X_g[,'dbh_m.init']*betaSize_g[i]
           y = truncnorm::rtruncnorm(n = nrow(pila_data$X_g),
                                     mean = XB,
                                     sd = sigmaEpsilon_g[i],
                                     a = 0)
           results = 
             pila_data$X_g %>%
             as_tibble()
           results$beta0_g = beta0_g[i]
           results$betaSize_g = betaSize_g[i]
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
ggplot(data = 
         growth_sims %>%
         mutate(beta0_bin = 
                  cut(beta0_g, breaks = quantile(beta0_g, probs = seq(0, 1, 0.1), include.lowest = TRUE)),
                betaSize_bin = 
                  cut(betaSize_g, breaks = quantile(betaSize_g, probs = seq(0, 1, 0.1), include.lowest = TRUE))) %>%
         group_by(beta0_bin, betaSize_bin) %>%
         summarise(size1_g = quantile(size1_g, prob = 0.99)),
       aes(fill = size1_g,
           x = beta0_bin, y = betaSize_bin))+
  geom_tile()+
  scale_fill_viridis_c()

growth_sims %>%
  mutate(beta0_bin = 
           cut(beta0_g, 
               breaks = 
                 seq(from = min(beta0_g)-0.1, to = max(beta0_g)+0.1, length.out = 10))) %>%
  ggplot(aes(x = beta0_bin, y = size1_g))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = -90))

growth_sims %>%
  mutate(betaSize_bin = 
           cut(betaSize_g, 
               breaks = 
                 seq(from = min(betaSize_g)-0.1, to = max(betaSize_g)+0.1, length.out = 10))) %>%
  ggplot(aes(x = betaSize_bin, y = size1_g))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = -90))

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
beta0_s = rnorm(n = 1000, mean = 0, sd = 1)
betaSize_s = rnorm(n = 1000, mean = 0, sd = 1)

# simulate responses
growth_sims = 
  do.call('bind_rows',
    lapply(X = 1:1000,
         FUN = function(i){
           XB = 
             pila_data$X_s[,'intercept']*beta0_s[i] + 
             pila_data$X_s[,'dbh_m.init']*betaSize_s[i]
           y = rbinom(n = nrow(pila_data$X_s),
                      size = 1,
                      prob = boot::inv.logit(XB))
           results = 
             pila_data$X_s %>%
             as_tibble()
           results$beta0_s = beta0_s[i]
           results$betaSize_s = betaSize_s[i]
           results$surv = y
           results$sim = i
           return(results)
         }))

# plot results
growth_sims %>%
  ggplot(aes(x = as.factor(surv)))+
  geom_bar()

growth_sims %>%
  mutate(beta0_bin = 
           cut(beta0_s, 
               breaks = seq(from = min(beta0_s)-0.01, 
                            to = max(beta0_s)+0.01,
                            length.out = 10))) %>%
  ggplot(aes(x = beta0_bin, fill = as.factor(surv)))+
  geom_bar(position = position_fill())


growth_sims %>%
  mutate(betaSize_bin = 
           cut(betaSize_s, 
               breaks = seq(from = min(betaSize_s)-0.01, 
                            to = max(betaSize_s)+0.01,
                            length.out = 10))) %>%
  ggplot(aes(x = betaSize_bin, fill = as.factor(surv)))+
  geom_bar(position = position_fill())

growth_sims %>%
  mutate(size0_bin = 
           cut(dbh_m.init,
               breaks = seq(from = min(dbh_m.init)-0.01, to = max(dbh_m.init)+0.01, length.out =10)),
         betaSize_bin = 
           cut(betaSize_s, 
               breaks = seq(from = min(betaSize_s)-0.01, 
                            to = max(betaSize_s)+0.01,
                            length.out = 10)),
         beta0_bin = 
           cut(beta0_s, 
               breaks = seq(from = min(beta0_s)-0.01, 
                            to = max(beta0_s)+0.01,
                            length.out = 10))) %>%
  ggplot(aes(x = size0_bin, fill = as.factor(surv)))+
  geom_bar(position = position_fill())+
  theme(axis.text.x = element_text(angle = -90))+
  facet_grid(beta0_bin~betaSize_bin)

#### priors for fecundity model ################################################
