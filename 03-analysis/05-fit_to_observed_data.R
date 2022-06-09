#### setup #####################################################################

# load packages
library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)


# load model data
pila_surv_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'pila_mort_training.rds'))

pila_growth_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'pila_growth_training.rds'))

pila_fecd_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'fecd_data_pila_usgs.rds'))

#### survival ##################################################################

# build model
surv_model = cmdstan_model(here::here('03-analysis', 'surv_model.stan'))

# run sampler
surv_fit = 
  surv_model$sample(
    data = pila_surv_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results'),
    output_basename = 'pila_surv',
    seed = 110819,
    adapt_delta = 0.8
  )

# check summary
surv_fit$summary()

# check diagnostics
surv_fit$cmdstan_diagnose()

# parameter pair plots
mcmc_pairs(surv_fit$draws(variables = c('beta[1]', 'beta[2]', 'beta[3]', 'sigma_plot', 'sigma_ecosub')))

# posterior density plots
mcmc_dens_overlay(surv_fit$draws(variables = c('beta', 'sigma_plot', 'sigma_ecosub')))

# save fitted model
surv_fit$save_object(here::here('02-data', '03-results', 'surv_fit.rds'))

# save posterior df
surv_posterior = as_draws_df(surv_fit$draws())

saveRDS(surv_posterior, here::here('02-data', '03-results', 'surv_post.rds'))

#### growth ####################################################################

# build model
growth_model = cmdstan_model(here::here('03-analysis', 'growth_model.stan'))

# run sampler
growth_fit = 
  growth_model$sample(
    data = pila_growth_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results'),
    output_basename = 'pila_growth',
    seed = 112188,
    adapt_delta = 0.8,
  )

# check summary
growth_fit$summary()

# check diagnostics
growth_fit$cmdstan_diagnose()

# parameter pair plots
mcmc_pairs(growth_fit$draws(variables = c('beta[1]', 'beta[2]', 'sigma_plot', 
                                          'sigma_ecosub', 'sigma_epsilon')))

# posterior density plots
mcmc_dens_overlay(growth_fit$draws(variables = c('beta', 'sigma_plot', 'sigma_ecosub')))

# save fitted model
growth_fit$save_object(here::here('02-data', '03-results', 'growth_fit.rds'))

# save posterior df
growth_posterior = as_draws_df(growth_fit$draws())

saveRDS(growth_posterior, here::here('02-data', '03-results', 'growth_post.rds'))

#### fecundity ###############################################################

# build model
fecd_model = cmdstan_model(here::here('03-analysis', 'fecd_model.stan'))

# run sampler
fecd_fit = 
  fecd_model$sample(
    data = pila_fecd_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results'),
    output_basename = 'pila_fecd',
    seed = 110819,
    adapt_delta = 0.8,
  )

# check summary
fecd_fit$summary()

# check diagnostics
fecd_fit$cmdstan_diagnose()

# parameter pair plots
mcmc_pairs(fecd_fit$draws(variables = c('beta[1]', 'beta[2]', 'beta[3]', 'sigma_ecosub')))

# posterior density plots
mcmc_dens_overlay(fecd_fit$draws(variables = c('beta', 'sigma_ecosub')))

# save fitted model
fecd_fit$save_object(here::here('02-data', '03-results', 'fecd_fit.rds'))

# save posterior df
fecd_posterior = as_draws_df(fecd_fit$draws())

saveRDS(fecd_posterior, here::here('02-data', '03-results', 'fecd_post.rds'))



#### build model and run sampler ###############################################



stan_model = 
  cmdstan_model(here::here('03-analysis', 'model.stan'))


fitted_model = 
  stan_model$sample(
    data = pila_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results', 'real_fits'),
    output_basename = 'pila_quadgrow',
    seed = 020188,
    adapt_delta = 0.8)



#### model diagnostics #########################################################


fitted_model$summary(c('beta_s',  'sigmaEco_s', 'sigmaPlot_s',
                       'beta_g', 'sigmaEco_g', 'sigmaPlot_g', 'sigmaEpsilon_g',
                       'beta_f', 'sigmaEco_f', 'kappa_r')) %>% 
  print(n = Inf)

# I think it's ok for there to be some correlation in the parameter esitmates 
# between the intercept and the size effect? the parameters are only weakly 
# identified but I don't have a good a priori reason to want to set beta_size 
# to 1
mcmc_pairs(fitted_model$draws(),
           pars = c('beta_g[1]', 'beta_g[2]','sigmaEpsilon_g', 'sigmaEco_g', 'sigmaPlot_g'))


mcmc_pairs(fitted_model$draws(),
           pars = c('beta_s[1]', 'beta_s[2]', 'beta_s[3]', 'sigmaEco_s', 'sigmaPlot_s'))

mcmc_pairs(fitted_model$draws(),
           pars = c('beta_f[1]', 'beta_f[2]', 'sigmaEco_f',  'kappa_r'))

mcmc_pairs(fitted_model$draws(),
           pars = c('beta_s[1]', 'beta_g[1]', 'beta_f[1]'))

mcmc_dens_overlay(fitted_model$draws(variables = 
                                       c('beta_s', 'sigmaEco_s','sigmaPlot_s')))

mcmc_dens_overlay(fitted_model$draws(variables = 
                                       c('beta_g', 'sigmaEco_g','sigmaPlot_s', 'sigmaEpsilon_g')))

mcmc_dens_overlay(fitted_model$draws(variables = 
                                       c('beta_f', 'sigmaEco_f', 'kappa_r')))

fitted_model$cmdstan_diagnose()

fitted_model$save_object(here::here('02-data', '03-results', 'real_fits', 'pila.rds'))

# save the posterior dataframe directly for portability
posterior = as_draws_df(fitted_model$draws())
saveRDS(posterior,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'posterior_draws.rds'))
