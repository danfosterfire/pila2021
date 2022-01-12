
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


#### simulate, estimate, rank ##################################################

sim_est_rank = 
  function(id){
    
    library(here)
    library(tidyverse)
    library(cmdstanr)
    library(bayesplot)
    library(posterior)
    
    # load the observed data and the stan model
    sim_data = readRDS(here::here('02-data', '02-for_analysis', 'pila_data.rds'))
    stan_model = cmdstanr::cmdstan_model(here::here('03-analysis','model_simple.stan'))
    
    # pull some true parameter values (these are normally pulled from the 
    # prior distribution, but the prior for beta_g[2] (starting size) is 
    # unrealistically wide and causes problems with the truncated normal 
    # response size1, so pull it from a fairly narrow prior around 1. I don't
    # include this narrow prior in the parameter estimation code because it 
    # would be a hassle to split out the fixeff for size from the rest of the 
    # fixed effects)
    true_beta_s = rnorm(n = 14, mean = 0, sd = 5)
    true_beta_g = rnorm(n = 14, mean = 0, sd = 5)
    true_sigmaEpsilon_g = 0
    while (true_sigmaEpsilon_g <= 0){
      true_sigmaEpsilon_g = rcauchy(n = 1, location = 0, scale = 5)
    }
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
    
    # get posterior parameter distributions
    sim_fit = 
      stan_model$sample(
        data = sim_data,
        init = 
          list(
            list(beta_s = rep(0, times = 14),
                 beta_g = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 sigmaEpsilon_g = 1),
            list(beta_s = rep(0, times = 14),
                 beta_g = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 sigmaEpsilon_g = 1),
            list(beta_s = rep(0, times = 14),
                 beta_g = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 sigmaEpsilon_g = 1),
            list(beta_s = rep(0, times = 14),
                 beta_g = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 sigmaEpsilon_g = 1)),
        parallel_chains = 1,
        output_dir = here::here('02-data', '03-results', 'sim_fits'),
        output_basename = paste0('simfit_', id),
        thin = 4,
        iter_warmup = 1000,
        iter_sampling = 1000)
    
    sim_fit$save_object(here::here('02-data', '03-results', 'sim_fits',
                                   paste0('simfit_',id,'.rds')))
    
    # compare posterior samples against true parameter values
    parameter_ranks = 
      
      # start with the (thinned) posterior samples
      as_draws_df(sim_fit$draws()) %>%
      
      # create indicator columns comparing posterior samples to true values
      mutate(
        Rbeta01_s = as.integer(`beta_s[1]` < true_beta_s[1]),
        Rbeta02_s = as.integer(`beta_s[2]` < true_beta_s[2]),
        Rbeta03_s = as.integer(`beta_s[3]` < true_beta_s[3]),
        Rbeta04_s = as.integer(`beta_s[4]` < true_beta_s[4]),
        Rbeta05_s = as.integer(`beta_s[5]` < true_beta_s[5]),
        Rbeta06_s = as.integer(`beta_s[6]` < true_beta_s[6]),
        Rbeta07_s = as.integer(`beta_s[7]` < true_beta_s[7]),
        Rbeta08_s = as.integer(`beta_s[8]` < true_beta_s[8]),
        Rbeta09_s = as.integer(`beta_s[9]` < true_beta_s[9]),
        Rbeta10_s = as.integer(`beta_s[10]` < true_beta_s[10]),
        Rbeta11_s = as.integer(`beta_s[11]` < true_beta_s[11]),
        Rbeta12_s = as.integer(`beta_s[12]` < true_beta_s[12]),
        Rbeta13_s = as.integer(`beta_s[13]` < true_beta_s[13]),
        Rbeta14_s = as.integer(`beta_s[14]` < true_beta_s[14]),
        
        Rbeta01_g = as.integer(`beta_g[1]` < true_beta_g[1]),
        Rbeta02_g = as.integer(`beta_g[2]` < true_beta_g[2]),
        Rbeta03_g = as.integer(`beta_g[3]` < true_beta_g[3]),
        Rbeta04_g = as.integer(`beta_g[4]` < true_beta_g[4]),
        Rbeta05_g = as.integer(`beta_g[5]` < true_beta_g[5]),
        Rbeta06_g = as.integer(`beta_g[6]` < true_beta_g[6]),
        Rbeta07_g = as.integer(`beta_g[7]` < true_beta_g[7]),
        Rbeta08_g = as.integer(`beta_g[8]` < true_beta_g[8]),
        Rbeta09_g = as.integer(`beta_g[9]` < true_beta_g[9]),
        Rbeta10_g = as.integer(`beta_g[10]` < true_beta_g[10]),
        Rbeta11_g = as.integer(`beta_g[11]` < true_beta_g[11]),
        Rbeta12_g = as.integer(`beta_g[12]` < true_beta_g[12]),
        Rbeta13_g = as.integer(`beta_g[13]` < true_beta_g[13]),
        Rbeta14_g = as.integer(`beta_g[14]` < true_beta_g[14]),
        
        RsigmaEpsilon_g = as.integer(sigmaEpsilon_g < true_sigmaEpsilon_g)
      ) %>%
      
      select(contains(c('Rbeta','RsigmaEpsilon'))) %>%
      
      # sum the indicators to get the rank for each variable
      summarise_all(sum)
    
    return(parameter_ranks)
    
    
  }

library(foreach)
library(doParallel)

registerDoParallel(12)

batch_param_ranks = 
  foreach(i = 1:100) %dopar% {
    sim_est_rank(i)
  }

saveRDS(batch_param_ranks,
        here::here('02-data','03-results', 'sim_fits', 'batch_param_ranks.rds'))

stopImplicitCluster()

#### check uniformity of each parameter ########################################

batch_param_ranks = 
  do.call('bind_rows',
          batch_param_ranks)

batch_param_ranks %>%
  rowid_to_column('run') %>%
  pivot_longer(cols = c(-run), names_to = 'param', values_to = 'rank') %>%
  ggplot(aes(x = rank))+
  geom_histogram(scale = 'free')+
  facet_wrap(~param)+
  theme_minimal()
