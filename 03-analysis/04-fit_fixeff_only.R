
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

pila_fit.model_noRanEff = 
  cmdstan_model(here::here('03-analysis', 'model_noRanEff.stan'))




# HAVING TROUBLE WITH THE RECRUITMENT MODEL; THE VARIABLE "r" KEEPS GETTING 
# SET TO A VECTOR OF NANS INSTEAD OF A VECTOR OF REALS. LOOKING AT THE 
# PARAMETER VALUES THIS IS HAPPENING WHEN NU IS A LARGE NEGATIVE NUMBER AND UPSILON IS 
# VERY SMALL; IT LOOKS LIKE THE SAMPLER CATCHES ON THAT THOSE VALUES SUCK AND STOPS 
# PROPOSING THEM
pila_fit.samples_noRanEff = 
  pila_fit.model_noRanEff$sample(data = pila_data,
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
                                  beta_g = 
                                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaEpsilon_g = 1,
                                  beta_f = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
                             list(beta_s = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  beta_g = 
                                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaEpsilon_g = 1,
                                  beta_f = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
                             list(beta_s = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  beta_g = 
                                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaEpsilon_g = 1,
                                  beta_f = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
                             list(beta_s = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  beta_g = 
                                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaEpsilon_g = 1,
                                  beta_f = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))),
                        parallel_chains = 4,
                        output_dir = here::here('02-data',
                                                '03-results'),
                        output_basename = 'noRanEff')


#### model diagnostics #########################################################


pila_fit.samples$summary(c('beta_s', 'sigmaPlot_s', 'sigmaEco_s', 'beta_g', 
                           'sigmaPlot_g', 'sigmaEco_g', 'sigmaEpsilon_g',
                           'beta_f', 'sigmaPlot_f', 'sigmaEco_f',
                           'upsilon', 'nu', 'kappa_r')) %>%
  print(n = Inf)

mcmc_trace(pila_fit.samples$draws(variables = 
                                    c('beta_s', 'sigmaPlot_s', 'sigmaEco_s', 'beta_g', 
                           'sigmaPlot_g', 'sigmaEco_g', 'sigmaEpsilon_g',
                           'beta_f', 'sigmaPlot_f', 'sigmaEco_f',
                           'upsilon', 'nu', 'kappa_r')))
mcmc_dens_overlay(pila_fit.samples$draws(variables = 
                                    c('beta_s', 'sigmaPlot_s', 'sigmaEco_s', 'beta_g', 
                           'sigmaPlot_g', 'sigmaEco_g', 'sigmaEpsilon_g',
                           'beta_f', 'sigmaPlot_f', 'sigmaEco_f',
                           'upsilon', 'nu', 'kappa_r')))

pila_fit.samples$cmdstan_diagnose()

pila_fit.samples_noRanEff = 
  as_cmdstan_fit(c(here::here('02-data', '03-results', 'noRanEff-1.csv'),
                 here::here('02-data', '03-results', 'noRanEff-2.csv'),
                 here::here('02-data', '03-results', 'noRanEff-3.csv'),
                 here::here('02-data', '03-results', 'noRanEff-4.csv')))

mcmc_pairs(pila_fit.samples_noRanEff$draws(),
           pars = c('beta_f[1]', 'upsilon', 'nu', 'kappa_r', 'sigmaEpsilon_g'),
           off_diag_fun = 'hex', 
           np = nuts_params(pila_fit.samples_noRanEff),
           condition = pairs_condition(nuts = 'accept_stat__'))


#### check for patterns in the residuals #######################################

# which random effects need to be included?

#### save model results ########################################################

pila_fit.samples$save_object(here::here('02-data', '03-results', 'pila_fit_mcmc_fullran.rds'))

