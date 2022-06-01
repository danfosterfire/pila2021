#### setup #####################################################################

# load packages
library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)


# load model data
tph_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'tph.rds'))

ba_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'ba.rds'))


#### tph #######################################################################

# build model
tph_model = cmdstan_model(here::here('03-analysis', 'tph_model.stan'))

# estimate parameters
tph_fit = 
  tph_model$sample(
    data = tph_data,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results', 'real_fits'),
    output_basename = 'tph',
    seed = 110819
      )

# model diagnostics
tph_fit$summary(c('beta', 'sigma_ecosub', 'sigma_epsilon'))

mcmc_pairs(tph_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]',
                    'beta[5]', 'beta[6]', 'sigma_ecosub', 'sigma_epsilon'))


mcmc_dens_overlay(tph_fit$draws(variables = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]',
                    'beta[5]', 'beta[6]', 'sigma_ecosub', 'sigma_epsilon')))

tph_fit$cmdstan_diagnose()

tph_fit$save_object(here::here('02-data',
                               '03-results',
                               'real_fits',
                               'tph_fit.rds'))

#### BA ########################################################################
# build model
ba_model = cmdstan_model(here::here('03-analysis', 'ba_model.stan'))

# estimate parameters
ba_fit = 
  ba_model$sample(
    data = ba_data,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results', 'real_fits'),
    output_basename = 'ba',
    seed = 110819
      )

# model diagnostics
ba_fit$summary(c('beta', 'sigma_ecosub', 'sigma_epsilon'))

mcmc_pairs(ba_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]',
                    'beta[5]', 'beta[6]', 'sigma_ecosub', 'sigma_epsilon'))


mcmc_dens_overlay(ba_fit$draws(variables = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]',
                    'beta[5]', 'beta[6]', 'sigma_ecosub', 'sigma_epsilon')))

ba_fit$cmdstan_diagnose()

ba_fit$save_object(here::here('02-data',
                               '03-results',
                               'real_fits',
                               'ba_fit.rds'))
