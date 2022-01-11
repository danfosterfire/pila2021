
# following shriver, doing three types of model validation:
# 1) simulate data using known parameters, make sure these are recovered when 
# model fit to simulated data
# 2) within-sample posterior predictive p-values
# 3) out-of-sample predictive checks


#### setup #####################################################################

library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)

# load observed data
pila_data = readRDS(here::here('02-data', '02-for_analysis', 'pila_data.rds'))

# compile the SBC code
SBC_model = cmdstan_model(here::here('03-analysis', 'SBC_model.stan'))


#### simulate, estimate, rank ##################################################

sim_est_rank = 
  function(id){
    
    # load the observed data and the stan model
    sim_data = readRDS(here::here('02-data', '02-for_analysis', 'pila_data.rds'))
    stan_model = cmdstan_model(here::here('03-analysis','model_noRanEff.stan'))
    
    # pull some true parameter values (these are normally pulled from the 
    # prior distribution, but the prior for beta_g[2] (starting size) is 
    # unrealistically wide and causes problems with the truncated normal 
    # response size1, so pull it from a fairly narrow prior around 1. I don't
    # include this narrow prior in the parameter estimation code because it 
    # would be a hassle to split out the fixeff for size from the rest of the 
    # fixed effects)
    true_beta_s = rnorm(n = 14, mean = 0, sd = 5)
    true_beta_g = rnorm(n = 14, mean = 0, sd = 5)
    true_beta_g[2] = rnorm(n = 1, mean = 1, sd = 0.1)
    true_sigmaEpsilon_g = rcauchy(n = 1, location = 0, scale = 5)
    
    # simulate some fake random effects
    
    # simulate some fake responses
    sim_surv = 
      rbinom(n = sim_data$N_s,
             prob = boot::inv.logit(as.numeric(sim_data$X_s %*% true_beta_s)),
             size = 1)
    sim_size1 = 
      truncnorm::rtruncnorm(n = sim_data$N_g,
                            a = 0,
                            mean = as.numeric(sim_data$X_g %*% true_beta_g),
                            sd = true_sigmaEpsilon_g)
    
    sim_data$surv = sim_surv
    sim_data$size1_g = sim_size1
  }

#### check uniformity of each parameter ########################################