
#### setup #####################################################################

#### setup #####################################################################

# load packages
library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)


# load model data
pila_data = 
  readRDS(here::here('02-data',
                   '02-for_analysis',
                   'pila_data.rds'))


#### build model and run sampler ###############################################

pila_fit.model = 
  cmdstan_model(here::here('03-analysis', 'model.stan'))



pila_fit.samples = 
  pila_fit.model$sample(data = pila_data,
                        init = 
                          
                          # the default initialization scheme (random draws in the 
                          # range -2:2) is not working well with the truncated 
                          # normal distribution for the growth model; samplers 
                          # are often rejecting more 
                          # than 100 sets of initial parameter values and then 
                          # giving up. Setting initial value for the fixed 
                          # effect of size0 to 1 helps by nudging the initial 
                          # state towards plausibly-positive means for the 
                          # size1 size. The values of 0 for other fixed effects 
                          # and 1 for the variances are fairly arbitrary defaults.
                          # the model takes a while to get in gear, so sampling is 
                        # slow at the start as most parameter proposals get 
                        # rejected for predicting a negative mean size1, but it 
                        # eventually wanders into a good region of parameter 
                        # space and samples well.
                          list(list(beta_s = 
                                      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    sigmaPlot_s = 1,
                                    beta_g = 
                                      c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    sigmaPlot_g = 1,
                                    sigmaEpsilon_g = 1),
                               list(beta_s = 
                                      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    sigmaPlot_s = 1,
                                    beta_g = 
                                      c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    sigmaPlot_g = 1,
                                    sigmaEpsilon_g = 1),
                               list(beta_s = 
                                      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    sigmaPlot_s = 1,
                                    beta_g = 
                                      c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    sigmaPlot_g = 1,
                                    sigmaEpsilon_g = 1),
                               list(beta_s = 
                                      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    sigmaPlot_s = 1,
                                    beta_g = 
                                      c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    sigmaPlot_g = 1,
                                    sigmaEpsilon_g = 1)),
                        parallel_chains = 4)


#### model diagnostics #########################################################

pila_fit.samples$summary(c('beta_s', 'sigmaPlot_s', 'beta_g', 'sigmaPlot_g', 'sigmaEpsilon_g')) %>%
  print(n = Inf)

mcmc_trace(pila_fit.samples$draws(variables = 
                                    c('beta_s', 'sigmaPlot_s', 
                                      'beta_g', 'sigmaPlot_g', 'sigmaEpsilon_g')))
mcmc_dens_overlay(pila_fit.samples$draws(variables = 
                                    c('beta_s', 'sigmaPlot_s', 
                                      'beta_g', 'sigmaPlot_g', 'sigmaEpsilon_g')))

pila_fit.samples$cmdstan_diagnose()

#### save model results ########################################################

pila_fit.samples$save_object(here::here('02-data', '03-results', 'pila_fit_mcmc.rds'))

