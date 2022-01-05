
library(here)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(tidyverse)

# load observed data
pila_data = readRDS(here::here('02-data', '02-for_analysis', 'pila_data.rds'))

# load mcmc results
pila_fit = readRDS(here::here('02-data', '03-results', 'pila_fit_mcmc_fullran.rds'))

samples = as_draws_df(pila_fit$draws())

# load subplot data to get coordinates for plotting
mort_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'mort_data.rds')) %>%
  filter(species=='PILA')%>%
  left_join(readRDS(here::here('02-data', '01-preprocessed', 'subplot_data.rds')) %>%
              select(plot_id, subp_id, lat, lon, ecosubcd)) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)))

# again, loading to get coordinates for plotting
growth_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'growth_data.rds')) %>%
  filter(species=='PILA')%>%
  left_join(readRDS(here::here('02-data', '01-preprocessed', 'subplot_data.rds')) %>%
              select(plot_id, subp_id, lat, lon, ecosubcd)) %>%
  # construct data for size:stressor interactions
  mutate(plot_id.i = as.integer(factor(plot_id)))

#### survival predictions from observed data ###################################

# extract fixed effect coefficients for survival model
betas_surv = select(samples, contains('beta_s'))
# extract plot level random effects for survival model
plotEffects_surv = select(samples, contains('plotEffect_s'))
# extract ecoregion level random effects for survival model
ecoEffects_surv = select(samples, contains('ecoEffect_s'))

# get the expected probability and residuals for each observed datapoint
survival_predictions = 
  do.call('rbind',
          lapply(X = sample(x = 1:4000, size = 100, replace = FALSE),
                 FUN = function(i){
                   
                   beta_i = as.numeric(betas_surv[i,])
                   plotEffects_i = as.numeric(plotEffects_surv[i,])
                   ecoEffects_i = as.numeric(ecoEffects_surv[i,])
                   
                   X = pila_data$X_s
                   plot_ids = pila_data$plotid_s
                   ecosub_ids = pila_data$ecosub_s
                   
                   predictions = X
                   predictions$plot_id = plot_ids
                   predictions$ecosub_id = ecosub_ids
                   predictions$logitp = 
                     as.numeric(as.matrix(X) %*% beta_i) + 
                     plotEffects_i[plot_ids] +
                     ecoEffects_i[ecosub_ids]
                   predictions$p = boot::inv.logit(predictions$logitp)
                   predictions$observed = pila_data$surv
                   predictions$pearson_resids = 
                     (predictions$observed-predictions$p) / 
                     sqrt(predictions$p * (1-predictions$p))
                   return(predictions)
                   
                 }))

head(survival_predictions)

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

# survival model is very good at predicting survival (actual survivors tend to 
# have very high values for p) but not as good at predicting death (trees that 
# died have p values all over the place); plotting the basic glm looks pretty 
# good; looks like the model slightly underpredicts survival and mortality at 
# intermediate levels of p; it's a little under-confident, which is ok.


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

#### plotting survival fixed effects ###########################################

predict_survival_fixeff = 
  function(dataset){
    do.call('rbind',
            lapply(X = 1:4000,
                   FUN = function(i){
                     
                     beta_i = as.numeric(betas_surv[i,])
                     
                     X = dataset[,c('intercept', 'dbh_in.init', 'fire', 
                                    'insects', 'disease', 'ba_scaled', 
                                    'cwd_dep90_scaled', 'cwd_mean_scaled', 
                                    'dbh_fire', 'dbh_insects', 'dbh_disease',
                                    'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')]
                     
                     predictions = X
                     predictions$logitp = 
                       as.numeric(as.matrix(X) %*% beta_i)
                     predictions$p = boot::inv.logit(predictions$logitp)
                     return(predictions)
                     
                   }))
    
  }

# first, "main effects" (other than the focal variable, disturbances held at false 
# and continuous variables held at 0 (other than size, they are scaled) 

## size main effect (no disturbances, continuous variables held at 0, which is 
## the subplot-level mean (so like the average subplot))
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_survival_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(p.median = median(p),
            p.025 = quantile(p, probs = 0.025),
            p.25 = quantile(p, probs = 0.25),
            p.75 = quantile(p, probs = 0.75),
            p.975 = quantile(p, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init))+
  geom_ribbon(aes(ymin = p.025, ymax = p.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = p.25, ymax = p.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = p.median), lwd = 1.5)+
  scale_y_continuous(limits = c(0, 1))+
  theme_minimal()

# size X fire
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = c(0,1),
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_survival_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(p.median = median(p),
            p.025 = quantile(p, probs = 0.025),
            p.25 = quantile(p, probs = 0.25),
            p.75 = quantile(p, probs = 0.75),
            p.975 = quantile(p, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(fire)))+
  geom_ribbon(aes(ymin = p.025, ymax = p.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = p.25, ymax = p.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = p.median), lwd = 1.5)+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  scale_y_continuous(limits = c(0, 1))

# size X insects
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = c(0,1),
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_survival_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(p.median = median(p),
            p.025 = quantile(p, probs = 0.025),
            p.25 = quantile(p, probs = 0.25),
            p.75 = quantile(p, probs = 0.75),
            p.975 = quantile(p, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(insects)))+
  geom_ribbon(aes(ymin = p.025, ymax = p.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = p.25, ymax = p.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = p.median), lwd = 1.5)+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  scale_y_continuous(limits = c(0, 1))

# size X disease
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = c(0,1),
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_survival_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(p.median = median(p),
            p.025 = quantile(p, probs = 0.025),
            p.25 = quantile(p, probs = 0.25),
            p.75 = quantile(p, probs = 0.75),
            p.975 = quantile(p, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(disease)))+
  geom_ribbon(aes(ymin = p.025, ymax = p.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = p.25, ymax = p.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = p.median), lwd = 1.5)+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  scale_y_continuous(limits = c(0, 1))

# size X BA
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = c(-1, 0, 1),
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_survival_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(p.median = median(p),
            p.025 = quantile(p, probs = 0.025),
            p.25 = quantile(p, probs = 0.25),
            p.75 = quantile(p, probs = 0.75),
            p.975 = quantile(p, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(ba_scaled)))+
  geom_ribbon(aes(ymin = p.025, ymax = p.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = p.25, ymax = p.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = p.median), lwd = 1.5)+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  scale_y_continuous(limits = c(0, 1))

# size X drought
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = c(-1, 0, 1),
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_survival_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(p.median = median(p),
            p.025 = quantile(p, probs = 0.025),
            p.25 = quantile(p, probs = 0.25),
            p.75 = quantile(p, probs = 0.75),
            p.975 = quantile(p, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(cwd_dep90_scaled)))+
  geom_ribbon(aes(ymin = p.025, ymax = p.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = p.25, ymax = p.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = p.median), lwd = 1.5)+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  scale_y_continuous(limits = c(0, 1))

# size X dryness
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = c(-1, 0, 1)) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_survival_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(p.median = median(p),
            p.025 = quantile(p, probs = 0.025),
            p.25 = quantile(p, probs = 0.25),
            p.75 = quantile(p, probs = 0.75),
            p.975 = quantile(p, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(cwd_mean_scaled)))+
  geom_ribbon(aes(ymin = p.025, ymax = p.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = p.25, ymax = p.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = p.median), lwd = 1.5)+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  scale_y_continuous(limits = c(0, 1))


#### growth predictions from observed data #####################################


# extract fixed effect coefficients for survival model
betas_grow = select(samples, contains('beta_g'))
# extract plot level random effects for survival model
plotEffects_grow = select(samples, contains('plotEffect_g'))
ecoEffects_grow = select(samples, contains('ecoEffect_g'))

# get the expected probability and residuals for each observed datapoint
growth_predictions = 
  do.call('rbind',
          lapply(X = sample(x = 1:4000, size = 100, replace = FALSE),
                 FUN = function(i){
                   
                   beta_i = as.numeric(betas_grow[i,])
                   plotEffects_i = as.numeric(plotEffects_grow[i,])
                   ecoEffects_i = as.numeric(ecoEffects_grow[i,])
                   
                   X = pila_data$X_g
                   plot_ids = pila_data$plotid_g
                   ecosub_ids = pila_data$ecosub_g
                   
                   predictions = X
                   predictions$plot_id = plot_ids
                   predictions$ecosub_id = pila_data$ecosub_g
                   predictions$mu = 
                     as.numeric(as.matrix(X) %*% beta_i) + 
                     plotEffects_i[plot_ids] + 
                     ecoEffects_i[ecosub_ids]
                   predictions$observed = pila_data$size1_g
                   predictions$resids = 
                     (predictions$observed-predictions$mu)
                   return(predictions)
                 }))

head(growth_predictions)

ggplot(data = growth_predictions,
       aes(y = observed, x = mu))+
  geom_point(size = 0, alpha = 0.1)+
  theme_minimal()+
  geom_smooth(method = 'lm')+
  geom_abline(intercept = 0, slope = 1, color = 'red')

# growth model looks good

# also want to plot the residuals, mostly across space to see if we need to 
# add a spatial random effect like the shriver paper did
growth_predictions %>%
  left_join(growth_data.pila %>%
              select(plot_id.i, lat, lon, ecosubcd) %>%
              group_by(plot_id.i, lat, lon, ecosubcd) %>%
              summarise() %>%
              ungroup(),
            by = c('plot_id' =  'plot_id.i')) %>%
  group_by(plot_id, lat, lon, ecosubcd) %>%
  summarise(resids = mean(resids)) %>%
  ungroup() %>%
  filter(resids < 1 & resids > -1) %>%
  ggplot(data = .,
         aes(x = lon, y = lat, color = resids))+
  geom_point(alpha = 0.7)+
  theme_minimal()+
  coord_sf()+
  scale_color_viridis_c()
# this looks good; needed the ecoregion random effect

growth_predictions %>%
  left_join(growth_data.pila %>%
              select(plot_id.i, lat, lon) %>%
              group_by(plot_id.i, lat, lon) %>%
              summarise() %>%
              ungroup(),
            by = c('plot_id' =  'plot_id.i')) %>%
  group_by(plot_id, lat, lon) %>%
  summarise(mu = mean(mu)) %>%
  ungroup() %>%
  ggplot(data = .,
         aes(x = lon, y = lat, color = mu))+
  geom_point(size = 1)+
  theme_minimal()+
  coord_sf()+
  scale_color_viridis_c()

#### fixed effects for growth model ############################################


predict_size2_fixeff = 
  function(dataset){
    do.call('rbind',
            lapply(X = 1:4000,
                   FUN = function(i){
                     
                     beta_i = as.numeric(betas_grow[i,])
                     
                     X = dataset[,c('intercept', 'dbh_in.init', 'fire', 
                                    'insects', 'disease', 'ba_scaled', 
                                    'cwd_dep90_scaled', 'cwd_mean_scaled', 
                                    'dbh_fire', 'dbh_insects', 'dbh_disease',
                                    'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')]
                     
                     predictions = X
                     predictions$XB = 
                       as.numeric(as.matrix(X) %*% beta_i)
                     predictions$growth = 
                       predictions$XB - predictions$dbh_in.init
                     return(predictions)
                     
                   }))
    
  }

# first, "main effects" (other than the focal variable, disturbances held at false 
# and continuous variables held at 0 (other than size, they are scaled) 

## size main effect (no disturbances, continuous variables held at 0, which is 
## the subplot-level mean (so like the average subplot))
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_size2_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(growth.median = median(growth),
            growth.025 = quantile(growth, probs = 0.025),
            growth.25 = quantile(growth, probs = 0.25),
            growth.75 = quantile(growth, probs = 0.75),
            growth.975 = quantile(growth, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init))+
  geom_ribbon(aes(ymin = growth.025, ymax = growth.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = growth.25, ymax = growth.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = growth.median), lwd = 1.5)+
  theme_minimal()

# size X fire
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = c(0,1),
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_size2_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(growth.median = median(growth),
            growth.025 = quantile(growth, probs = 0.025),
            growth.25 = quantile(growth, probs = 0.25),
            growth.75 = quantile(growth, probs = 0.75),
            growth.975 = quantile(growth, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(fire)))+
  geom_ribbon(aes(ymin = growth.025, ymax = growth.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = growth.25, ymax = growth.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = growth.median), lwd = 1.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  theme_minimal()

expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = c(0,1),
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_size2_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(growth.median = median(growth),
            growth.025 = quantile(growth, probs = 0.025),
            growth.25 = quantile(growth, probs = 0.25),
            growth.75 = quantile(growth, probs = 0.75),
            growth.975 = quantile(growth, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(insects)))+
  geom_ribbon(aes(ymin = growth.025, ymax = growth.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = growth.25, ymax = growth.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = growth.median), lwd = 1.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  theme_minimal()

# size X disease
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = c(0,1),
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_size2_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(growth.median = median(growth),
            growth.025 = quantile(growth, probs = 0.025),
            growth.25 = quantile(growth, probs = 0.25),
            growth.75 = quantile(growth, probs = 0.75),
            growth.975 = quantile(growth, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(disease)))+
  geom_ribbon(aes(ymin = growth.025, ymax = growth.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = growth.25, ymax = growth.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = growth.median), lwd = 1.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  theme_minimal()

# size X basal area
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = c(-1, 0, 1),
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_size2_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(growth.median = median(growth),
            growth.025 = quantile(growth, probs = 0.025),
            growth.25 = quantile(growth, probs = 0.25),
            growth.75 = quantile(growth, probs = 0.75),
            growth.975 = quantile(growth, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(ba_scaled)))+
  geom_ribbon(aes(ymin = growth.025, ymax = growth.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = growth.25, ymax = growth.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = growth.median), lwd = 1.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  theme_minimal()


expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = c(-1, 0, 1),
            'cwd_mean_scaled' = 0) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_size2_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(growth.median = median(growth),
            growth.025 = quantile(growth, probs = 0.025),
            growth.25 = quantile(growth, probs = 0.25),
            growth.75 = quantile(growth, probs = 0.75),
            growth.975 = quantile(growth, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(cwd_dep90_scaled)))+
  geom_ribbon(aes(ymin = growth.025, ymax = growth.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = growth.25, ymax = growth.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = growth.median), lwd = 1.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  theme_minimal()

# size X dryness
expand.grid('intercept' = 1,
            'dbh_in.init' = seq(from = 0.5, to = 49.5, by = 1),
            'fire' = 0,
            'insects' = 0,
            'disease' = 0,
            'ba_scaled' = 0,
            'cwd_dep90_scaled' = 0,
            'cwd_mean_scaled' = c(-1, 0, 1)) %>%
  mutate(dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled) %>%
  predict_size2_fixeff() %>%
  group_by(intercept, dbh_in.init, fire, insects, disease, ba_scaled, 
           cwd_dep90_scaled, cwd_mean_scaled, dbh_fire, dbh_insects, dbh_disease,
           dbh_ba, dbh_cwd90, dbh_cwdmean) %>%
  summarise(growth.median = median(growth),
            growth.025 = quantile(growth, probs = 0.025),
            growth.25 = quantile(growth, probs = 0.25),
            growth.75 = quantile(growth, probs = 0.75),
            growth.975 = quantile(growth, probs = 0.975)) %>%
  ggplot(data = .,
         aes(x = dbh_in.init, fill = as.factor(cwd_mean_scaled)))+
  geom_ribbon(aes(ymin = growth.025, ymax = growth.975),
              color = NA, alpha = 0.1)+
  geom_ribbon(aes(ymin = growth.25, ymax = growth.75),
              color = NA, alpha = 0.25)+
  geom_line(aes(y = growth.median), lwd = 1.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85)+
  theme_minimal()

