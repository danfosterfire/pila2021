
# following shriver, doing three types of model validation:
# 1) simulate data using known parameters, make sure these are recovered when 
# model fit to simulated data -> the way shriver et al did this is actually not 
# valid; you're not supposed to simulate a single data set and fit it, you're 
# supposed to simulate a whole ensemble of datasets and fit them. You're supposed 
# to select the priors for this simulation from the prior distribution. I can 
# see how this is a good way to check the model, but I'm starting to wonder if 
# i'm actually doing anything useful here. The truncated normal distribution for
# growth is finnickey enough that it's easy to generate pathological datasets
# (I've already carefully calibrated the priors for intercept and size0), and 
# i'm not seeing an obvious way to calibrate the priors for additional 
# parameters. I don't feel like i'm actually improving the model, just 
# carefully tailoring the priors to make sure that in a 100 datasets I dont
# have any pathological ones, which is getting difficult to do with more 
# parameters. I've also bumped up against the prior restrictions I really 
# think I can impose apriori just based on domain knowledge. I suspect that 
# restrictiong the SD of the additional beta parameters would make the 
# simulated datasets more consistently nice, but just based off of 
# domain knowledge I don't feel like I have a good reason to restrict them 
# further. An SD for effect of fire/disease/drought/BA of ~1cm/year seems 
# very reasonable to me. 
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
    true_beta_g = rnorm(n = sim_data$K, mean = 0, sd = 0.1)
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
        Rbeta03_s = as.integer(`beta_s[3]` < simulated_data[[id]]$true_beta_s[3]),
        Rbeta04_s = as.integer(`beta_s[4]` < simulated_data[[id]]$true_beta_s[4]),
        Rbeta05_s = as.integer(`beta_s[5]` < simulated_data[[id]]$true_beta_s[5]),
        Rbeta06_s = as.integer(`beta_s[6]` < simulated_data[[id]]$true_beta_s[6]),
        Rbeta07_s = as.integer(`beta_s[7]` < simulated_data[[id]]$true_beta_s[7]),
        
        Rbeta01_g = as.integer(`beta_g[1]` < simulated_data[[id]]$true_beta_g[1]),
        Rbeta02_g = as.integer(`beta_g[2]` < simulated_data[[id]]$true_beta_g[2]),
        Rbeta03_g = as.integer(`beta_g[3]` < simulated_data[[id]]$true_beta_g[3]),
        Rbeta04_g = as.integer(`beta_g[4]` < simulated_data[[id]]$true_beta_g[4]),
        Rbeta05_g = as.integer(`beta_g[5]` < simulated_data[[id]]$true_beta_g[5]),
        Rbeta06_g = as.integer(`beta_g[6]` < simulated_data[[id]]$true_beta_g[6]),
        Rbeta07_g = as.integer(`beta_g[7]` < simulated_data[[id]]$true_beta_g[7]),
        
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
  geom_hline(yintercept = qbinom(0.005, size = 100, prob = 1/10))+
  geom_hline(yintercept = qbinom(0.995, size = 100, prob = 1/10))+
  facet_wrap(~param, scales = 'free')+
  theme_minimal()

#### check model diagnostics ###################################################

model_diagnostics = 
  sapply(X = c(1:100),
       FUN = function(i){
         sim_fit = readRDS(here::here('02-data', '03-results', 'sim_fits', 
                                      paste0('simfit_', i, '.rds')))
         
         return(sim_fit$cmdstan_diagnose()$stdout)
       })

# check these to see which sets of parameters were causing the trouble
bad_batches = 
  c(1:100)[!str_detect(model_diagnostics, pattern = 'Processing complete, no problems detected.')]

matrix(ncol = simulated_data[[1]]$sim_data$K,
       nrow = 100,
       byrow = TRUE,
       data = sapply(X = 1:100,
                     FUN = function(i){
                       simulated_data[[i]]$true_beta_s
                     }),
       dimnames = list(c(1:100), 
                       c('beta_1','beta_2','beta_3', 'beta_4', 'beta_5', 
                         'beta_6', 'beta_7'))) %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),
               names_to = 'param', values_to = 'value') %>%
  ggplot(data = .,
         aes(x = value))+
  geom_histogram()+
  theme_minimal()+
  facet_wrap(~param)+
  geom_vline(
    data =matrix(ncol = simulated_data[[1]]$sim_data$K,
                 nrow = length(bad_batches),
                 byrow = TRUE,
                 data = sapply(X = bad_batches,
                               FUN = function(i){
                                 simulated_data[[i]]$true_beta_s
                               }),
                 dimnames = list(bad_batches, 
                                 c('beta_1','beta_2','beta_3', 'beta_4', 'beta_5', 
                                   'beta_6', 'beta_7'))) %>%
      as_tibble() %>%
      pivot_longer(cols = everything(),
                   names_to = 'param', values_to = 'value'),
    aes(xintercept = value), color = 'red')


matrix(ncol = simulated_data[[1]]$sim_data$K,
       nrow = 100,
       byrow = TRUE,
       data = sapply(X = 1:100,
                     FUN = function(i){
                       simulated_data[[i]]$true_beta_g
                     }),
       dimnames = list(c(1:100), 
                       c('beta_1','beta_2','beta_3', 'beta_4', 'beta_5', 
                         'beta_6', 'beta_7'))) %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),
               names_to = 'param', values_to = 'value') %>%
  ggplot(data = .,
         aes(x = value))+
  geom_histogram()+
  theme_minimal()+
  facet_wrap(~param, scales = 'free')+
  geom_vline(
    data =matrix(ncol = simulated_data[[1]]$sim_data$K,
                 nrow =length(bad_batches),
                 byrow = TRUE,
                 data = sapply(X = bad_batches,
                               FUN = function(i){
                                 simulated_data[[i]]$true_beta_g
                               }),
                 dimnames = list(bad_batches, 
                                 c('beta_1','beta_2','beta_3', 'beta_4', 'beta_5', 
                                   'beta_6', 'beta_7'))) %>%
      as_tibble() %>%
      pivot_longer(cols = everything(),
                   names_to = 'param', values_to = 'value'),
    aes(xintercept = value), color = 'red')

matrix(ncol = simulated_data[[1]]$sim_data$K,
                 nrow =length(bad_batches),
                 byrow = TRUE,
                 data = sapply(X = bad_batches,
                               FUN = function(i){
                                 simulated_data[[i]]$true_beta_g
                               }),
                 dimnames = list(bad_batches, 
                                 c('beta_1','beta_2','beta_3', 'beta_4', 'beta_5', 
                                   'beta_6', 'beta_7')))

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


