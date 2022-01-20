#### setup #####################################################################

# load packages
library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)


# load model data
pila_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'pila_training.rds'))


#### build model and run sampler ###############################################

stan_model = 
  cmdstan_model(here::here('03-analysis', 'model.stan'))


fitted_model = 
  stan_model$sample(
    data = pila_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results', 'real_fits'),
    output_basename = 'pila',
    seed = 110819)



#### model diagnostics #########################################################


fitted_model$summary(c('beta_s',  'sigmaEco_s', 'sigmaPlot_s',
                       'beta_g', 'sigmaEco_g', 'sigmaPlot_g', 'sigmaEpsilon_g',
                       'beta_f', 'sigmaEco_f', 'sigmaPlot_f', 'kappa_r')) %>% 
  print(n = Inf)

# I think it's ok for there to be some correlation in the parameter esitmates 
# between the intercept and the size effect? the parameters are only weakly 
# identified but I don't have a good a priori reason to want to set beta_size 
# to 1
mcmc_pairs(fitted_model$draws(),
           pars = c('beta_g[1]', 'beta_g[2]','sigmaEpsilon_g', 'sigmaEco_g', 'sigmaPlot_g'))

mcmc_pairs(fitted_model$draws(),
           pars = c('beta_s[1]', 'beta_s[2]', 'sigmaEco_s', 'sigmaPlot_s'))

mcmc_pairs(fitted_model$draws(),
           pars = c('beta_f[1]', 'beta_f[2]', 'sigmaEco_f', 'sigmaPlot_f', 'kappa_r'))

mcmc_pairs(fitted_model$draws(),
           pars = c('beta_s[1]', 'beta_g[1]', 'beta_f[1]'))

mcmc_dens_overlay(fitted_model$draws(variables = 
                                       c('beta_s', 'sigmaEco_s','sigmaPlot_s')))

mcmc_dens_overlay(fitted_model$draws(variables = 
                                       c('beta_g', 'sigmaEco_g','sigmaPlot_s', 'sigmaEpsilon_g')))

mcmc_dens_overlay(fitted_model$draws(variables = 
                                       c('beta_f', 'sigmaEco_f','sigmaPlot_f', 'kappa_r')))

fitted_model$cmdstan_diagnose()

fitted_model$save_object(here::here('02-data', '03-results', 'real_fits', 'pila.rds'))
