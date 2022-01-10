
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

# load mcmc results
pila_fit = readRDS(here::here('02-data', '03-results', 'pila_fit_mcmc_fullran.rds'))

samples.real = as_draws_df(pila_fit$draws())

samples.real %>% summarise_all(median)
#### simulate data  ############################################################

# DEFINE PARAMETERS
true_params = 
  list(
    beta_s = samples.real %>% select(contains('beta_s')) %>% summarise_all(median) %>% as.numeric(),
    sigmaPlot_s = samples.real %>% select(contains('sigmaPlot_s')) %>% summarise_all(median) %>% as.numeric(),
    sigmaEco_s = samples.real %>% select(contains('sigmaEco_s')) %>% summarise_all(median) %>% as.numeric(),
    
    beta_g = samples.real %>% select(contains('beta_g')) %>% summarise_all(median) %>% as.numeric(),
    sigmaPlot_g = samples.real %>% select(contains('sigmaPlot_g')) %>% summarise_all(median) %>% as.numeric(),
    sigmaEco_g = samples.real %>% select(contains('sigmaEco_g')) %>% summarise_all(median) %>% as.numeric(),
    sigmaEpsilon_g = samples.real %>% select(contains('sigmaEpsilon_g')) %>% summarise_all(median) %>% as.numeric(),
    
    beta_f = samples.real %>% select(contains('beta_f')) %>% summarise_all(median) %>% as.numeric(),
    sigmaPlot_f = samples.real %>% select(contains('sigmaPlot_f')) %>% summarise_all(median) %>% as.numeric(),
    sigmaEco_f = samples.real %>% select(contains('sigmaEco_f')) %>% summarise_all(median) %>% as.numeric(),
    kappa_r = samples.real %>% select(contains('kappa_r')) %>% summarise_all(median) %>% as.numeric()
    #kappa_r = 1
  )


# BUILD EXPLANATORY VARIABLES
names(pila_data)

sim_explanatory = 
  list(N_s = pila_data$N_s, # number of individuals observed for survival (live at t = 0)
       P_s = pila_data$P_s, # number of plots in survival data
       E_s = pila_data$E_s, # number of ecoregions in survival data
       plotid_s = pila_data$plotid_s, # plotids vector
       ecosub_s = pila_data$ecosub_s, # ecosub id vector
       X_s = pila_data$X_s, # fixed effect covariates; using the actual values here 
                            # to get the most accurate distributions
       
       N_g = pila_data$N_g, # number of individuals observed for growth (live at t=0 and t=1)
       P_g = pila_data$P_g, # number of plots in growth data
       E_g = pila_data$E_g, # number of ecosubs in growth data
       plotid_g = pila_data$plotid_g, # plot ids
       ecosub_g = pila_data$ecosub_g, # ecosub ids
       X_g = pila_data$X_g, # fixeff covariates
       
       N_r = pila_data$N_r, # number of rows in recruitment data; 1 row per 
                            # sizeclass:subplot
       S_r = pila_data$S_r, # number of subplots in recruitment data
       P_f = pila_data$P_f, # number of plots in the fecundity model (plots 
                            # which are in both the recr and surv models)
       E_f = pila_data$E_f, # see above
       X_r = pila_data$X_r, # fixeff covariates for recruitment data, used for 
                            # the fecundity model as well as predicting growth 
                            # and survival for each DBH class; importantly 
                            # this is ordered increasing by sizeclass within each subplot
       plotid_sr = pila_data$plotid_sr, # ids linking rows in recruitment data to survival plot raneffects
       ecosub_sr = pila_data$ecosub_sr, # " linking to survival ecoregion raneffs
       plotid_gr = pila_data$plotid_gr, # " linking to growth plot raneffs
       ecosub_gr = pila_data$ecosub_gr, # " linking to growth ecosub raneffs
       plotid_fr = pila_data$plotid_fr, # " linking to fecundity plot raneffs
       ecosub_fr = pila_data$ecosub_fr, # " linking to fecudnity ecosub raneffs
       M_r = pila_data$M_r, # number of distinct size classes
       u_bounds = pila_data$u_bounds, # size class upper bounds
       l_bounds = pila_data$l_bounds, # size class lower bounds
       midpoints = pila_data$midpoints, # size class midpoints
       a = pila_data$a, # plot areas for each size class
       n = pila_data$n, # number of individuals in each size class at time t = 0,
       r = pila_data$r # new recruit size kernel
       )

# SIMULATE RANDOM EFFECTS AND RESPONSES

plotEffects_s = rnorm(n = sim_explanatory$P_s, mean = 0, sd = true_params$sigmaPlot_s)
ecoEffects_s = rnorm(n = sim_explanatory$E_s, mean = 0, sd = true_params$sigmaEco_s)
plotEffects_g = rnorm(n = sim_explanatory$P_g, mean = 0, sd = true_params$sigmaPlot_g)
ecoEffects_g = rnorm(n = sim_explanatory$E_g, mean = 0, sd = true_params$sigmaEco_g)
plotEffects_f = rnorm(n = sim_explanatory$P_f, mean = 0, sd = true_params$sigmaPlot_f)
ecoEffects_f = rnorm(n = sim_explanatory$E_f, mean = 0, sd = true_params$sigmaEco_f)

logitp_s =  
  as.numeric(as.matrix(sim_explanatory$X_s) %*% true_params$beta_s) +
  plotEffects_s[sim_explanatory$plotid_s] +
  ecoEffects_s[sim_explanatory$ecosub_s]
p_s = boot::inv.logit(logitp_s)
surv = rbinom(n = sim_explanatory$N_s, size = 1, prob = p_s)

mu_g = 
  as.numeric(as.matrix(sim_explanatory$X_g) %*% true_params$beta_g) +
  plotEffects_g[sim_explanatory$plotid_g] +
  ecoEffects_g[sim_explanatory$ecosub_g]
size1_g = truncnorm::rtruncnorm(n = sim_explanatory$N_g, a = 0, mean = mu_g, sd = true_params$sigmaEpsilon_g)

logf = 
  as.numeric(as.matrix(sim_explanatory$X_r) %*% true_params$beta_f) +
  plotEffects_f[sim_explanatory$plotid_fr] +
  ecoEffects_f[sim_explanatory$ecosub_fr]

f = matrix(nrow = sim_explanatory$M_r, ncol = sim_explanatory$S_r,
           data = exp(logf), byrow = FALSE)

logits = 
  as.numeric(as.matrix(sim_explanatory$X_r) %*% true_params$beta_s) +
  plotEffects_s[sim_explanatory$plotid_sr] + 
  ecoEffects_s[sim_explanatory$ecosub_sr] 
s = boot::inv.logit(logits)

mu_gr = 
  as.numeric(as.matrix(sim_explanatory$X_r) %*% true_params$beta_g)+
  plotEffects_g[sim_explanatory$plotid_gr] +
  ecoEffects_g[sim_explanatory$ecosub_gr]

g = array(dim = c(2,
                  1,
                  sim_explanatory$S_r), 
          dimnames = list(toclass = sim_explanatory$midpoints[1:2], 
                          fromclass = 1,
                          subplot = 1:sim_explanatory$S_r),
          data = 
            sapply(X = 1:sim_explanatory$S_r,
                   FUN = function(d){
                     sapply(X = 1,
                            FUN = function(j){
                              sapply(X = 1:sim_explanatory$M_r,
                                     FUN = function(h){
                                        #paste0('to',h,';from',j,';on',d) #testing
                                       (pnorm(q = sim_explanatory$u_bounds[h],
                                             mean = mu_gr[j+(sim_explanatory$M_r * (d-1))],
                                             sd = true_params$sigmaEpsilon_g) - 
                                         pnorm(q = sim_explanatory$l_bounds[h],
                                               mean = mu_gr[j+sim_explanatory$M_r*(d-1)],
                                               sd = true_params$sigmaEpsilon_g))/
                                         (1-pnorm(q = 0, 
                                                  mean = mu_gr[j+(sim_explanatory$M_r*(d-1))],
                                                  sd = true_params$sigmaEpsilon_g))
                                     })
                            })
                   }))

A = array(dim = c(2,
                   sim_explanatory$M_r,
                   sim_explanatory$S_r),
          dimnames = list(toclass = sim_explanatory$midpoints[1:2],
                          fromclass = sim_explanatory$midpoints,
                          subplot = 1:sim_explanatory$S_r),
          data = sapply(X = 1:sim_explanatory$S_r,
                        FUN = function(d){
                          sapply(X = 1:sim_explanatory$M_r,
                                 FUN = function(j){
                                   sapply(X = 1:2,
                                          FUN = function(h){
                                            
                                            g[h,1,d]*s[1+(sim_explanatory$M_r*(d-1))] +
                                              sim_explanatory$r[h]*f[j,d]
                                          })
                                 })
                        }))

nprime = matrix(nrow = 2, 
                ncol = sim_explanatory$S_r,
                byrow = FALSE,
                data = sapply(X = 1:sim_explanatory$S_r,
                              FUN = function(d){
                                as.numeric(A[,,d] %*% sim_explanatory$n[d,])
                              }))

cprime = 
  matrix(nrow = 2, 
         ncol = sim_explanatory$S_r,
         byrow = FALSE,
         data = sapply(X = 1:sim_explanatory$S_r,
                       FUN = function(d){
                         rnbinom(mu = nprime[,d]*sim_explanatory$a,
                                 size = true_params$kappa_r,
                                 n = 2)
                       }))

#### try to recover parameters #################################################

sim_data = 
  sim_explanatory

sim_data$surv = surv
sim_data$size1_g = size1_g
sim_data$cprime = cprime

saveRDS(sim_data,
        here::here('02-data',
                   '02-for_analysis',
                   'sim_data.rds'))

sim_fit.model = cmdstan_model(here::here('03-analysis', 'model.stan'))

sim_fit.samples = 
  sim_fit.model$sample(data = sim_data,
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
                                  sigmaPlot_s = 1,
                                  sigmaEco_s = 1,
                                  beta_g = 
                                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_g = 1,
                                  sigmaEco_g = 1,
                                  sigmaEpsilon_g = 1,
                                  beta_f = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_f = 1,
                                  sigmaEco_f = 1),
                             list(beta_s = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_s = 1,
                                  sigmaEco_s = 1,
                                  beta_g = 
                                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_g = 1,
                                  sigmaEco_g = 1,
                                  sigmaEpsilon_g = 1,
                                  beta_f = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_f = 1,
                                  sigmaEco_f = 1),
                             list(beta_s = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_s = 1,
                                  sigmaEco_s = 1,
                                  beta_g = 
                                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_g = 1,
                                  sigmaEco_g = 1,
                                  sigmaEpsilon_g = 1,
                                  beta_f = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_f = 1,
                                  sigmaEco_f = 1),
                             list(beta_s = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_s = 1,
                                  sigmaEco_s = 1,
                                  beta_g = 
                                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_g = 1,
                                  sigmaEco_g = 1,
                                  sigmaEpsilon_g = 1,
                                  beta_f = 
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                  sigmaPlot_f = 1,
                                  sigmaEco_f = 1)),
                        parallel_chains = 4,
                        output_dir = here::here('02-data',
                                                '03-results'),
                        output_basename = 'sim_model')


sim_fit.samples$save_object(here::here('02-data', '03-results', 'sim_fit.rds'))



#### compare results to true values ############################################

sim_fit.samples$summary(c('beta_s', 'sigmaPlot_s', 'sigmaEco_s', 'beta_g', 
                           'sigmaPlot_g', 'sigmaEco_g', 'sigmaEpsilon_g',
                           'beta_f', 'sigmaPlot_f', 'sigmaEco_f','kappa_r')) %>%
  print(n = Inf)

mcmc_trace(sim_fit.samples$draws(variables = 
                                    c('beta_s', 'sigmaPlot_s', 'sigmaEco_s', 'beta_g', 
                                      'sigmaPlot_g', 'sigmaEco_g', 'sigmaEpsilon_g',
                                      'beta_f', 'sigmaPlot_f', 'sigmaEco_f','kappa_r')))


mcmc_dens_overlay(sim_fit.samples$draws(variables = 
                                           c('beta_s', 'sigmaPlot_s', 'sigmaEco_s', 'beta_g', 
                                             'sigmaPlot_g', 'sigmaEco_g', 'sigmaEpsilon_g',
                                             'beta_f', 'sigmaPlot_f', 'sigmaEco_f','kappa_r')))

sim_fit.samples$cmdstan_diagnose()

samples = as_draws_df(sim_fit.samples$draws(variables = 
                                               c('beta_s', 'sigmaPlot_s', 'sigmaEco_s', 'beta_g', 
                                             'sigmaPlot_g', 'sigmaEco_g', 'sigmaEpsilon_g',
                                             'beta_f', 'sigmaPlot_f', 'sigmaEco_f', 'kappa_r')))

samples_all = as_draws_df(sim_fit.samples$draws())

samples_all %>%
  select(contains('ecoEffect_f')) %>%
  summarise_all(median) %>%
  pivot_longer(cols = everything()) %>%
  pull(value)

ecoEffects_f

samples_all %>%
  select(contains('ecoEffect_g')) %>%
  summarise_all(median) %>%
  pivot_longer(cols = everything()) %>%
  pull(value)

ecoEffects_g

head(samples)

# beta_s
lapply(X = 1:14,
       FUN = function(i){
         
         param_name = paste0('beta_s[',i,']')
         
         bounds95 = quantile(pull(samples[,param_name]), probs = c(0.025, 0.975))
         
         ggplot(data = samples,
                aes(x = .data[[param_name]]))+
           geom_histogram()+
           geom_vline(xintercept = true_params$beta_s[i], color = 'red')+
           geom_vline(xintercept = bounds95[1], color = 'blue')+
           geom_vline(xintercept = bounds95[2], color = 'blue')+
           theme_minimal()
       })

# beta_g
lapply(X = 1:14,
       FUN = function(i){
         
         param_name = paste0('beta_g[',i,']')
         
         bounds95 = quantile(pull(samples[,param_name]), probs = c(0.025, 0.975))
         
         ggplot(data = samples,
                aes(x = .data[[param_name]]))+
           geom_histogram()+
           geom_vline(xintercept = true_params$beta_g[i], color = 'red')+
           geom_vline(xintercept = bounds95[1], color = 'blue')+
           geom_vline(xintercept = bounds95[2], color = 'blue')+
           theme_minimal()
       })

# beta_f
lapply(X = 1:14,
       FUN = function(i){
         
         param_name = paste0('beta_f[',i,']')
         
         bounds95 = quantile(pull(samples[,param_name]), probs = c(0.025, 0.975))
         
         ggplot(data = samples,
                aes(x = .data[[param_name]]))+
           geom_histogram()+
           geom_vline(xintercept = true_params$beta_f[i], color = 'red')+
           geom_vline(xintercept = bounds95[1], color = 'blue')+
           geom_vline(xintercept = bounds95[2], color = 'blue')+
           theme_minimal()
       })

lapply(X = c('sigmaPlot_s', 'sigmaPlot_g', 'sigmaPlot_f',
             'sigmaEco_s', 'sigmaEco_g', 'sigmaEco_f',
             'sigmaEpsilon_g', 'kappa_r'),
       FUN = function(p){
         bounds95 = quantile(pull(samples[,p]), probs = c(0.025, 0.975))
         ggplot(data = samples,
                aes(x = .data[[p]]))+
           geom_histogram()+
           geom_vline(xintercept = true_params[[p]], color = 'red')+
           geom_vline(xintercept = bounds95[1], color = 'blue')+
           geom_vline(xintercept = bounds95[2], color = 'blue')+
           theme_minimal()
       })




# beta estimates wonky because the raneff estimates are bad?
