#### setup #####################################################################

# load packages
library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)



surv_data = readRDS(here::here('02-data',
                               '02-for_analysis',
                               'surv_data_pila_usgs.rds'))

growth_data = readRDS(here::here('02-data',
                                 '02-for_analysis',
                                 'growth_data_pila_usgs.rds'))

fecd_data = readRDS(here::here('02-data',
                               '02-for_analysis',
                               'fecd_data_pila_usgs.rds'))


#### survival ##################################################################

# build model
surv_model = cmdstan_model(here::here('03-analysis', 'surv_model_usgs.stan'))

# run sampler
surv_fit = 
  surv_model$sample(
    data = surv_data,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results', 'real_fits', 'usgs'),
    output_basename = 'pila_surv',
    seed = 110819
  )

# model diagnostics
surv_fit$cmdstan_diagnose()

surv_fit$summary()

mcmc_pairs(surv_fit$draws(variables = c('beta', 'sigma_plot')))

mcmc_dens_overlay(surv_fit$draws(variables = c('beta', 'sigma_plot')))

surv_fit$save_object(here::here('02-data', 
                                '03-results',
                                'real_fits',
                                'usgs',
                                'surv_fit_usgs.rds'))

surv_posterior = as_draws_df(surv_fit$draws())

saveRDS(surv_posterior,
        here::here('02-data', 
                   '03-results',
                   'real_fits',
                   'usgs',
                   'surv_posterior_usgs.rds'))


#### growth ####################################################################

# build model
growth_model = cmdstan_model(here::here('03-analysis', 'growth_model_usgs.stan'))

# run sampler
growth_fit = 
  growth_model$sample(
    data = growth_data,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results', 'real_fits', 'usgs'),
    output_basename = 'pila_growth',
    seed = 110819
  )

# model diagnostics
growth_fit$cmdstan_diagnose()

growth_fit$summary()

mcmc_pairs(growth_fit$draws(variables = c('beta', 'sigma_plot')))

mcmc_dens_overlay(growth_fit$draws(variables = c('beta', 'sigma_plot')))

growth_fit$save_object(here::here('02-data', 
                                '03-results',
                                'real_fits',
                                'usgs',
                                'growth_fit_usgs.rds'))

growth_posterior = as_draws_df(growth_fit$draws())

saveRDS(growth_posterior,
        here::here('02-data', 
                   '03-results',
                   'real_fits',
                   'usgs',
                   'growth_posterior_usgs.rds'))


#### survival ##################################################################

# build model
fecd_model = cmdstan_model(here::here('03-analysis', 'fecd_model_usgs.stan'))

# run sampler
fecd_fit = 
  fecd_model$sample(
    data = fecd_data,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results', 'real_fits', 'usgs'),
    output_basename = 'pila_fecd',
    seed = 110819,
    adapt_delta = 0.95
  )

# model diagnostics
fecd_fit$cmdstan_diagnose()

fecd_fit$summary()

mcmc_pairs(fecd_fit$draws(variables = c('beta', 'sigma_plot')))

mcmc_dens_overlay(fecd_fit$draws(variables = c('beta', 'sigma_plot')))

fecd_fit$save_object(here::here('02-data', 
                                '03-results',
                                'real_fits',
                                'usgs',
                                'fecd_fit_usgs.rds'))

fecd_posterior = as_draws_df(fecd_fit$draws())

saveRDS(fecd_posterior,
        here::here('02-data', 
                   '03-results',
                   'real_fits',
                   'usgs',
                   'fecd_posterior_usgs.rds'))



