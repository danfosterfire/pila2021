
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
library(foreach)
library(doParallel)

#### function defs #############################################################


simulate_data = 
  function(id){
    
    # load in the observed X data
    sim_data = readRDS(here::here('02-data', '02-for_analysis', 'pila_data.rds'))
    
    # pull some true parameter values (these are normally pulled from the 
    # prior distribution, but the prior for beta_g[2] (starting size) is 
    # unrealistically wide and causes problems with the truncated normal 
    # response size1, so pull it from a fairly narrow prior around 1.
    true_beta_s = rnorm(n = sim_data$K, mean = 0, sd = 1)
    true_beta_g = rnorm(n = sim_data$K, mean = 0, sd = 0.25)
    true_beta_g[1] = rnorm(n = 1, mean = 0.15, sd = 0.1)
    true_beta_g[2] = rnorm(n = 1, mean = 1, sd = 0.25)
    true_sigmaEpsilon_g = truncnorm::rtruncnorm(n = 1, mean = 0, sd = 0.25, a = 0)
    
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
    
    results = list('id' = id,
                   'true_beta_s' = true_beta_s,
                   'true_beta_g' = true_beta_g,
                   'true_sigmaEpsilon_g' = true_sigmaEpsilon_g,
                   'sim_data' = sim_data)
    
    return(results)
  }



estimate_parameters = 
  function(id){
    library(here)
    library(cmdstanr)
    
    # load the data and the stan model
    sim_data = simulated_data[[id]]$sim_data
    stan_model = cmdstanr::cmdstan_model(here::here('03-analysis', 'model_simple.stan'))
    
    # fit the model
    sim_fit = 
      stan_model$sample(
        data = sim_data,
        parallel_chains = 1,
        output_dir = here::here('02-data', '03-results', 'sim_fits'),
        output_basename = paste0('simfit_', id),
        thin = 4,
        iter_warmup = 1000,
        iter_sampling = 1000)
    
    sim_fit$save_object(here::here('02-data', '03-results', 'sim_fits',
                                   paste0('simfit_',id,'.rds')))
    
    print(paste0('successful run ', id))
    return(paste0('Successful run ', id))
  }


rank_parameters = 
  function(id){
    sim_fit = readRDS(here::here('02-data',
                                 '03-results',
                                 'sim_fits',
                                 paste0('simfit_',id,'.rds')))
    
    # compare posterior samples against true parameter values
    parameter_ranks = 
      
      # start with the (thinned) posterior samples
      as_draws_df(sim_fit$draws()) %>%
      
      # create indicator columns comparing posterior samples to true values
      mutate(
        Rbeta01_s = as.integer(`beta_s[1]` < simulated_data[[id]]$true_beta_s[1]),
        Rbeta02_s = as.integer(`beta_s[2]` < simulated_data[[id]]$true_beta_s[2]),
        
        Rbeta01_g = as.integer(`beta_g[1]` < simulated_data[[id]]$true_beta_g[1]),
        Rbeta02_g = as.integer(`beta_g[2]` < simulated_data[[id]]$true_beta_g[2]),
        RsigmaEpsilon_g = as.integer(sigmaEpsilon_g < simulated_data[[id]]$true_sigmaEpsilon_g)
      ) %>%
      
      select(contains(c('Rbeta','RsigmaEpsilon'))) %>%
      
      # sum the indicators to get the rank for each variable
      summarise_all(sum)
    
    return(parameter_ranks)
    
  }



#### simulate, estimate, rank ##################################################

set.seed(110819)
simulated_data = 
  lapply(X = 1:100,
         FUN = function(i){
           simulate_data(i)
         })

registerDoParallel(cores = 12)

foreach(i = 1:100) %dopar% {
    estimate_parameters(i)
  }

stopImplicitCluster()


batch_param_ranks = 
  do.call('bind_rows',
          lapply(X = c(1:100),
       FUN = function(i){
         rank_parameters(i)
       }))



#### check uniformity of each parameter ########################################
batch_param_ranks %>%
  rowid_to_column('run') %>%
  pivot_longer(cols = c(-run), names_to = 'param', values_to = 'rank') %>%
  mutate(rank_sim = sample(1:1000, size = nrow(.), replace = TRUE)) %>%
  mutate(rank_bin = cut(rank, breaks = seq(from = 0, to = 1001, by = 100), include.lowest = TRUE)) %>%
  ggplot(aes(x = rank_bin))+
  geom_bar()+
  facet_wrap(~param, scales = 'free')+
  theme_minimal()

#### check model diagnostics ###################################################

model_diagnostics = 
  sapply(X = c(1:100),
       FUN = function(i){
         sim_fit = readRDS(here::here('02-data', '03-results', 'sim_fits', 
                                      paste0('simfit_', i, '.rds')))
         
         sim_fit$cmdstan_diagnose()$stdout
       })

model_diagnostics[1]
c(1:100)[!str_detect(model_diagnostics, pattern = 'Processing complete, no problems detected.')]


#### scratch ###################################################################

# fit the real data and compare the posteriors against the priors to see how 
# informative these priors are
pila_data = readRDS(here::here('02-data', '02-for_analysis', 'pila_data.rds'))

stan_model = cmdstan_model(here::here('03-analysis', 'model_simple.stan'))

pila_fit = stan_model$sample(data = pila_data,
                             parallel_chains = 4)

pila_fit$cmdstan_diagnose()

samples = as_draws_df(pila_fit$draws())

names(samples)

ggplot()+
  geom_density(data = samples,
               aes(x = `beta_s[1]`),
               color = 'blue',
               lwd = 1)+
  geom_density(data = data.frame(prior = rnorm(n = 4000, mean = 0, sd = 1)),
               aes(x = prior),
               color = 'red', lty = 2, lwd = 1)+
  theme_minimal()

ggplot()+
  geom_density(data = samples,
               aes(x = `beta_s[2]`),
               color = 'blue',
               lwd = 1)+
  geom_density(data = data.frame(prior = rnorm(n = 4000, mean = 0, sd = 1)),
               aes(x = prior),
               color = 'red', lty = 2, lwd = 1)+
  theme_minimal()

ggplot()+
  geom_density(data = samples,
               aes(x = `beta_g[1]`),
               color = 'blue',
               lwd = 1)+
  geom_density(data = data.frame(prior = rnorm(n = 4000, mean = 0.15, sd = 0.1)),
               aes(x = prior),
               color = 'red', lty = 2, lwd = 1)+
  theme_minimal()

ggplot()+
  geom_density(data = samples,
               aes(x = `beta_g[2]`),
               color = 'blue',
               lwd = 1)+
  geom_density(data = data.frame(prior = rnorm(n = 4000, mean = 1, sd = 0.25)),
               aes(x = prior),
               color = 'red', lty = 2, lwd = 1)+
  theme_minimal()+
  scale_x_continuous(limits = c(0.95, 1.05))

ggplot()+
  geom_density(data = samples,
               aes(x = sigmaEpsilon_g),
               color = 'blue',
               lwd = 1)+
  geom_density(data = data.frame(prior = truncnorm::rtruncnorm(n = 4000, mean = 0, sd = 0.25, a = 0)),
               aes(x = prior),
               color = 'red', lty = 2, lwd = 1)+
  theme_minimal()


