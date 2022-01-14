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
  cmdstan_model(here::here('03-analysis', 'model_simple.stan'))


fitted_model = 
  stan_model$sample(
    data = pila_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results', 'real_fits'),
    output_basename = 'pila')



#### model diagnostics #########################################################


pila_fit.samples_noRanEff$summary(c('beta_s', 'beta_g', 'sigmaEpsilon_g',
                                    'beta_f', 'kappa_r')) %>%
  print(n = Inf)

lapply(X = c('beta_s', 'beta_g', 'sigmaEpsilon_g',
             'beta_f', 'kappa_r'),
       FUN = function(v){
         mcmc_trace(pila_fit.samples_noRanEff$draws(variables = v))})

lapply(X = c('beta_s', 'beta_g', 'sigmaEpsilon_g',
             'beta_f', 'kappa_r'),
       FUN = function(v){
         mcmc_dens_overlay(pila_fit.samples_noRanEff$draws(variables = v))})


mcmc_pairs(pila_fit.samples_noRanEff$draws(),
           pars = c('beta_f[1]', 'kappa_r', 'sigmaEpsilon_g'),
           off_diag_fun = 'hex', 
           np = nuts_params(pila_fit.samples_noRanEff),
           condition = pairs_condition(nuts = 'accept_stat__'))