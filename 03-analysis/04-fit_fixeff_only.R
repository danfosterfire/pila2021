
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

#pila_fit.samples_noRanEff = 
#  as_cmdstan_fit(c(here::here('02-data', '03-results', 'noRanEff-1.csv'),
#                   here::here('02-data', '03-results', 'noRanEff-2.csv'),
#                   here::here('02-data', '03-results', 'noRanEff-3.csv'),
#                   here::here('02-data', '03-results', 'noRanEff-4.csv')))

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


#### check for patterns in the residuals #######################################

# which random effects need to be included?
samples = as_draws_df(pila_fit.samples_noRanEff$draws())

# load subplot data to get coordinates for plotting
mort_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'mort_data.rds')) %>%
  filter(species=='PILA')%>%
  left_join(readRDS(here::here('02-data', '01-preprocessed', 'subplot_data.rds')) %>%
              select(plot_id, subp_id, lat, lon, ecosubcd)) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)),
         ecosub.i = as.integer(factor(ecosubcd)))

# again, loading to get coordinates for plotting
growth_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'growth_data.rds')) %>%
  filter(species=='PILA')%>%
  left_join(readRDS(here::here('02-data', '01-preprocessed', 'subplot_data.rds')) %>%
              select(plot_id, subp_id, lat, lon, ecosubcd)) %>%
  # construct data for size:stressor interactions
  mutate(plot_id.i = as.integer(factor(plot_id)),
         ecosub.i = as.integer(factor(ecosubcd)))


# extract fixed effect coefficients for survival model
betas_surv = select(samples, contains('beta_s'))

# get the expected probability and residuals for each observed datapoint
beta_s = 
  betas_surv %>%
  summarise_all(median) %>%
  as.numeric()
survival_predictions = pila_data$X_s
survival_predictions$plot_id = pila_data$plotid_s
survival_predictions$logitp = as.numeric(as.matrix(pila_data$X_s) %*% beta_s)
survival_predictions$p = boot::inv.logit(survival_predictions$logitp)
survival_predictions$observed = pila_data$surv
survival_predictions$pearson_resids = 
  (survival_predictions$observed - survival_predictions$p) / 
  sqrt(survival_predictions$p * (1-survival_predictions$p))

ggplot(data = survival_predictions,
       aes(y = p, x = as.factor(observed)))+
  geom_boxplot()+
  coord_flip()+
  theme_minimal()

ggplot(data = survival_predictions,
       aes(x = p, y = observed))+
  geom_jitter(size = 0, height = 0.1, width = 0, alpha = 0.1)+
  theme_minimal()+
  geom_smooth(method = 'glm', formula = y~x, method.args = list(family = 'binomial'))

ggplot(data = survival_predictions,
       aes(x = p, y = as.factor(observed)))+
  geom_violin()+
  theme_minimal() 



# also want to plot the residuals, mostly across space to see if we need to 
# add a spatial random effect like the shriver paper did
survival_predictions %>%
  left_join(mort_data.pila %>%
              select(plot_id.i, lat, lon) %>%
              group_by(plot_id.i, lat, lon) %>%
              summarise() %>%
              ungroup(),
            by = c('plot_id' =  'plot_id.i')) %>%
  group_by(plot_id, lat, lon) %>%
  summarise(pearson_resids = mean(pearson_resids)) %>%
  ungroup() %>%
  ggplot(data = .,
         aes(x = lon, y = lat, color = pearson_resids))+
  geom_point(size = 1)+
  theme_minimal()+
  coord_sf()+
  scale_color_viridis_c()
# this looks good

ggplot(data = 
           survival_predictions %>%
           filter(plot_id>800&plot_id<900),
         aes(x = as.factor(plot_id), y = pearson_resids))+
  geom_point(alpha = 0.5)
# not sure what to make of this tbh, probably not the most helpful diagnostic

survival_predictions %>%
  left_join(mort_data.pila %>%
              select(plot_id.i, lat, lon) %>%
              group_by(plot_id.i, lat, lon) %>%
              summarise() %>%
              ungroup(),
            by = c('plot_id' =  'plot_id.i')) %>%
  group_by(plot_id, lat, lon) %>%
  summarise(p = mean(p)) %>%
  ungroup() %>%
  ggplot(data = .,
         aes(x = lon, y = lat, color = p))+
  geom_point(size = 1)+
  theme_minimal()+
  coord_sf()+
  scale_color_viridis_c()

# extract fixed effect coefficients for survival model
betas_grow = select(samples, contains('beta_g'))

# get the expected probability and residuals for each observed datapoint
beta_g = 
  betas_grow %>%
  summarise_all(median) %>%
  as.numeric()
growth_predictions = pila_data$X_g
growth_predictions$plot_id = pila_data$plotid_g
growth_predictions$mu = as.numeric(as.matrix(pila_data$X_g) %*% beta_g)
growth_predictions$observed = pila_data$size1_g
growth_predictions$resids = 
  (growth_predictions$observed - growth_predictions$mu)

ggplot(data = growth_predictions,
       aes(x = mu, y = observed))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, color = 'red')+
  geom_smooth(method = 'lm')+
  theme_minimal()



# also want to plot the residuals, mostly across space to see if we need to 
# add a spatial random effect like the shriver paper did
growth_predictions %>%
  left_join(growth_data.pila %>%
              select(plot_id.i, lat, lon) %>%
              group_by(plot_id.i, lat, lon) %>%
              summarise() %>%
              ungroup(),
            by = c('plot_id' =  'plot_id.i')) %>%
  group_by(plot_id, lat, lon) %>%
  summarise(resids = mean(resids)) %>%
  ungroup() %>%
  filter(resids<1 & resids>-1) %>%
  ggplot(data = .,
         aes(x = lon, y = lat, color = resids))+
  geom_point(size = 1)+
  theme_minimal()+
  coord_sf()+
  scale_color_viridis_c()
# needs an ecoregion effect

ggplot(data = 
           growth_predictions %>%
           filter(plot_id>800&plot_id<900),
         aes(x = as.factor(plot_id), y = resids))+
  geom_point(alpha = 0.5)
# also wants a plot random effect



