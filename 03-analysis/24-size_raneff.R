
# load packages
library(here)
library(tidyverse)

set.seed(110819)


# a row for every subplot, columns for the subplot-level covariates like 
# disturbance data, basal area, and drought; includes all plots where 
# pila was present at initial or followup
plot_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'plot_data.rds')) %>%
  mutate(ba_scaled = as.numeric(scale(ba_ft2ac)),
         cwd_dep90_scaled = as.numeric(scale(cwd_departure90)),
         cwd_mean_scaled = as.numeric(scale(cwd_mean)),
         intercept = 1)




# a row for every individual tagged tree which was alive at initial and 
# remeasurement, and columns for the size of the tree at each time; 
# joined to the subplot data to get covariates for each individual tree 
# and interactions of those with size
growth_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'growth_data.rds')) %>%
  filter(species=='PILA')%>%
  left_join(plot_data %>%
              select(plot_id, plot_id, ecosubcd, 
                     intercept,
                     fire, insects, disease, wpbr,
                     ba_scaled, cwd_dep90_scaled, cwd_mean_scaled)) %>%
  
  # get DBH in meters, which will be on a nicer (close to 0-1, similar to other 
  # covariates) scale
  mutate(dbh_m.init = dbh_in.init*0.0254,
         dbh_m.re = dbh_in.re * 0.0254,
         dbh_m2.init = dbh_m.init**2,
         dbh_m3.init = dbh_m.init**3,
         dbh_fire = dbh_m.init*fire,
         dbh2_fire = dbh_m2.init*fire,
         dbh3_fire = dbh_m3.init*fire,
         dbh_insects = dbh_m.init*insects,
         dbh_disease = dbh_m.init*disease,
         dbh2_disease = dbh_m2.init*disease,
         dbh3_disease = dbh_m3.init*disease,
         dbh_ba = dbh_m.init*ba_scaled,
         dbh2_ba = dbh_m2.init*ba_scaled,
         dbh3_ba = dbh_m3.init*ba_scaled,
         dbh_cwd90 = dbh_m.init*cwd_dep90_scaled,
         dbh2_cwd90 = dbh_m2.init*cwd_dep90_scaled,
         dbh3_cwd90 = dbh_m3.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_m.init*cwd_mean_scaled,
         dbh2_cwdmean = dbh_m2.init*cwd_mean_scaled,
         dbh3_cwdmean = dbh_m3.init*cwd_mean_scaled,
         dbh_wpbr = dbh_m.init*wpbr,
         dbh2_wpbr = dbh_m2.init*wpbr,
         dbh3_wpbr = dbh_m3.init*wpbr)


# a row for every individual tagged tree which was alive at the initial 
# measurement, and columns indicating its survival status at remeasurement, 
# plus covariate columns for plot-level data and interactions with size
mort_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'mort_data.rds')) %>%
  filter(species=='PILA')%>%
  left_join(plot_data %>%
              select(plot_id, plot_id, ecosubcd,
                     intercept,
                     fire, insects, disease, wpbr,
                     ba_scaled, cwd_dep90_scaled, cwd_mean_scaled)) %>%
  
  # convert to meters scale
  mutate(
    dbh_m.init = dbh_in.init*0.0254,
    dbh_m2.init = dbh_m.init**2,
    dbh_m3.init = dbh_m.init**3,
    dbh_fire = dbh_m.init*fire,
    dbh2_fire = dbh_m2.init*fire,
    dbh3_fire = dbh_m3.init*fire,
    dbh_insects = dbh_m.init*insects,
    dbh_disease = dbh_m.init*disease,
    dbh2_disease = dbh_m2.init*disease,
    dbh3_disease = dbh_m3.init*disease,
    dbh_ba = dbh_m.init*ba_scaled,
    dbh2_ba = dbh_m2.init*ba_scaled,
    dbh3_ba = dbh_m3.init*ba_scaled,
    dbh_cwd90 = dbh_m.init*cwd_dep90_scaled,
    dbh2_cwd90 = dbh_m2.init*cwd_dep90_scaled,
    dbh3_cwd90 = dbh_m3.init*cwd_dep90_scaled,
    dbh_cwdmean = dbh_m.init*cwd_mean_scaled,
    dbh2_cwdmean = dbh_m2.init*cwd_mean_scaled,
    dbh3_cwdmean = dbh_m3.init*cwd_mean_scaled,
    dbh_wpbr = dbh_m.init*wpbr,
    dbh2_wpbr = dbh_m2.init*wpbr,
    dbh3_wpbr = dbh_m3.init*wpbr)

# a row for each unique combination of plot:species:size class, for 
# 1 inch size bins from 1.5-99.5"; "tpa_unadj.init" and "tpa_unadj.re" give 
# the area-adjusted density of stems in each species:size 
# bin for each plot at the initial and remeasurement. dbh_class gives the 
# ID of each size bin (integer equal to the open lower bound of the bin, so 
# stems 1"<=dbh<2" go in class 1, stems with dbh = 2.5" go in bin 2, etc.)
sizedist_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'sizedist_data.rds')) %>%
  filter(species=='PILA')

# a row for each size class, for 1 inch size bins from 1.5-99.5"; gives 
# metadata about each class including the midpoints, upper and lower bounds, 
# and the sampling area in the FIA design.
size_metadata = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.rds'))

mort_data.pila %>% filter(plot_id =='41-1-39-79537') %>% print(width = Inf)
growth_data.pila %>% filter(plot_id == '41-1-39-79537') %>% print(width = Inf)
sizedist_data.pila %>% filter(plot_id == '41-1-39-79537')

# a row for each plot, with only plots where PILA was present at 
# initial measurement

# a row for each size class and plot, with only plots where PILA was present 
# at initial measurement
recr_data.pila = 
  plot_data %>%
  filter(is.element(plot_id, mort_data.pila$plot_id)) %>%
  
  # filter to only subplots where both initial and remeasurement had manual 
  # greater than or equal to 2
  tidyr::expand(nesting(plt_cn, prev_plt_cn, plot_id, elev_ft, lat, lon,
                        ecosubcd, invdate.re, inv_manual.re, macro_break, fire,
                        insects, disease, cutting, invdate.init, inv_manual.init, 
                        ba_ft2ac, cwd_departure90, cwd_mean, ba_scaled, cwd_dep90_scaled,
                        cwd_mean_scaled, intercept, wpbr),
                dbh_in.init = size_metadata$bin_midpoint) %>%
  mutate(dbh_class = cut(dbh_in.init,
                         breaks = seq(from = 1, to= 100, by = 1),
                         labels = FALSE,
                         right = FALSE),
         dbh_m.init = dbh_in.init*0.0254,
         dbh_m2.init = dbh_m.init**2,
         dbh_m3.init = dbh_m.init**3,
         dbh_fire = dbh_m.init*fire,
         dbh2_fire = dbh_m2.init*fire,
         dbh3_fire = dbh_m3.init*fire,
         dbh_insects = dbh_m.init*insects,
         dbh_disease = dbh_m.init*disease,
         dbh2_disease = dbh_m2.init*disease,
         dbh3_disease = dbh_m3.init*disease,
         dbh_ba = dbh_m.init*ba_scaled,
         dbh2_ba = dbh_m2.init*ba_scaled,
         dbh3_ba = dbh_m3.init*ba_scaled,
         dbh_cwd90 = dbh_m.init*cwd_dep90_scaled,
         dbh2_cwd90 = dbh_m2.init*cwd_dep90_scaled,
         dbh3_cwd90 = dbh_m3.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_m.init*cwd_mean_scaled,
         dbh2_cwdmean = dbh_m2.init*cwd_mean_scaled,
         dbh3_cwdmean = dbh_m3.init*cwd_mean_scaled,
         dbh_wpbr = dbh_m.init*wpbr,
         dbh2_wpbr = dbh_m2.init*wpbr,
         dbh3_wpbr = dbh_m3.init*wpbr) %>%
  # add in the observed counts 
  left_join(readRDS(here::here('02-data',
                               '01-preprocessed',
                               'recruits_data.rds')) %>%
              filter(species=='PILA')) %>%
  
  left_join(sizedist_data.pila %>%
              select(plot_id, species, dbh_class, tpa_unadj.init)) %>%
  
  # there are 26 plots which, despite being included in the mort data (ie, 
  # have PILA which are alive at initial measure), 
  # have 0 TPA of pila in any size class at initial measurement. Upon investigating,
  # I think these are all cases where the species ID of the tree was changed 
  # from something else to PILA, so the tree number shows alive at initial 
  # and remasurement and the species at remeasurement is PILA (ensuring the 
  # subplot is included in the growth and mortality data frames) but 0 trees 
  # in PILA size bins at initial measurement. solution is to just drop the 
  # offending subplots from the recruitment dataset, where the all 0s 
  # were screwing up the IPM model estimates for seedling density.
filter(!is.element(plot_id,
                   group_by(., plot_id) %>%
                     summarise(tpa_unadj.init = sum(tpa_unadj.init)) %>%
                     ungroup() %>%
                     filter(tpa_unadj.init==0) %>%
                     pull(plot_id))) %>%
  
  # order by size and subplot
  arrange(plot_id, dbh_in.init) %>%
  rename(recruit_count = count) 

recr_data.pila %>%
  pull(plot_id) %>%
  unique() %>%
  length()



#### asssign plot and ecoregion indices ########################################

# going to use a single unified index across all the datasets, which cleans 
# up the code substantially and *i think* should process ok?
union_plots = 
  growth_data.pila %>%
  select(plot_id) %>%
  bind_rows(mort_data.pila %>%
              select(plot_id)) %>%
  bind_rows(recr_data.pila %>%
              select(plot_id)) %>%
  group_by(plot_id) %>%
  summarise() %>%
  ungroup() %>%
  arrange(plot_id) %>%
  mutate(plot_id.i = as.integer(factor(plot_id)))

union_ecosubs = 
  growth_data.pila %>%
  select(ecosubcd) %>%
  bind_rows(mort_data.pila %>%
              select(ecosubcd)) %>%
  bind_rows(recr_data.pila %>%
              select(ecosubcd))%>%
  group_by(ecosubcd) %>%
  summarise() %>%
  ungroup() %>%
  arrange(ecosubcd)%>%
  mutate(ecosub.i = as.integer(factor(ecosubcd)))

growth_data.pila = 
  growth_data.pila %>%
  left_join(union_plots) %>%
  left_join(union_ecosubs)

mort_data.pila = 
  mort_data.pila %>%
  left_join(union_plots) %>%
  left_join(union_ecosubs) 

recr_data.pila = 
  recr_data.pila %>%
  left_join(union_plots) %>%
  left_join(union_ecosubs)


#### get distribution of inventory intervals ###################################

plot_data %>%
  filter(is.element(plot_id, union_plots$plot_id)) %>%
  mutate(interval_days = invdate.re - invdate.init,
         interval_years = as.numeric(interval_days) / 365) %>%
  pull(interval_years) %>%
  summary()


plot_data %>%
  filter(is.element(plot_id, union_plots$plot_id)) %>%
  mutate(interval_days = invdate.re - invdate.init,
         interval_years = as.numeric(interval_days) / 365) %>%
  pull(interval_years) %>%
  quantile(., probs = c(0.1, 0.9))


#### split into training and validation data ###################################

validation_plots = 
  recr_data.pila %>%
  group_by(plot_id) %>%
  summarise() %>%
  ungroup() %>%
  sample_frac(size = 0.1, replace = FALSE) %>%
  pull(plot_id)


growth_data.pila_training = 
  growth_data.pila %>%
  filter(!is.element(plot_id, validation_plots))

growth_data.pila_validation = 
  growth_data.pila %>%
  filter(is.element(plot_id, validation_plots))

mort_data.pila_training = 
  mort_data.pila %>%
  filter(!is.element(plot_id, validation_plots))

mort_data.pila_validation = 
  mort_data.pila %>%
  filter(is.element(plot_id, validation_plots))

recr_data.pila_training = 
  recr_data.pila %>%
  filter(!is.element(plot_id, validation_plots))

recr_data.pila_validation = 
  recr_data.pila %>%
  filter(is.element(plot_id, validation_plots))

# a row for each unique combination of plot:species:size class, for 1 
# inch size bins from 0.5-9.5"; "count" gives the number of untagged 
# trees on the plot in the species:size bin at the remeasurement (ie 
# ingrowth)

untagged_data.pila_training = 
  recr_data.pila_training %>%
  arrange(plot_id, dbh_in.init)

untagged_data.pila_validation = 
  recr_data.pila_validation %>%
  arrange(plot_id, dbh_in.init)

#### prepare mortality training and validation data ############################

pila_mort_training = 
  list(K = 3,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(mort_data.pila_training),
       surv = as.integer(mort_data.pila_training$survived),
       X = 
         as.matrix(mort_data.pila_training[,c('intercept', 'dbh_m.init', 
                                              'dbh_m2.init')]),
       plot_id = mort_data.pila_training$plot_id.i,
       ecosub_id = mort_data.pila_training$ecosub.i)


pila_mort_validation = 
  list(K = 3,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(mort_data.pila_validation),
       surv = as.integer(mort_data.pila_validation$survived),
       X = 
         as.matrix(mort_data.pila_validation[,c('intercept', 'dbh_m.init', 
                                                'dbh_m2.init')]),
       plot_id = mort_data.pila_validation$plot_id.i,
       ecosub_id = mort_data.pila_validation$ecosub.i)



#### prepare growth training and validation ####################################


pila_growth_training = 
  list(K = 3,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(growth_data.pila_training),
       size1 = growth_data.pila_training$dbh_m.re,
       X = 
         as.matrix(growth_data.pila_training[,c('intercept', 'dbh_m.init', 
                                                'dbh_m2.init')]),
       plot_id = growth_data.pila_training$plot_id.i,
       ecosub_id = growth_data.pila_training$ecosub.i)


pila_growth_validation = 
  list(K = 3,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(growth_data.pila_validation),
       size1 = growth_data.pila_validation$dbh_m.re,
       X = 
         as.matrix(growth_data.pila_validation[,c('intercept', 'dbh_m.init', 
                                                  'dbh_m2.init')]),
       plot_id = growth_data.pila_validation$plot_id.i,
       ecosub_id = growth_data.pila_validation$ecosub.i)


#### prepare fecudnity training and validation #################################

pila_fecd_training = 
  list(K = 1,
       N = nrow(recr_data.pila_training),
       E = nrow(union_ecosubs),
       P = length(unique(recr_data.pila_training$plot_id.i)),
       M = max(size_metadata$bin_id),
       X = 
         as.matrix(recr_data.pila_training[,c('intercept')]),
       plot_id = recr_data.pila_training$plot_id.i,
       ecosub_id = recr_data.pila_training$ecosub.i,
       cprime = 
         recr_data.pila_training %>% 
         group_by(plot_id, recruit_count) %>% 
         summarise() %>%
         ungroup() %>%
         arrange(plot_id) %>%
         pull(recruit_count),
       a = size_metadata$plot_area_ac[1],
       n = matrix(ncol = length(unique(recr_data.pila_training$plot_id.i)),
                  nrow = nrow(size_metadata),
                  data = recr_data.pila_training$tpa_unadj.init,
                  byrow = FALSE))


pila_fecd_validation = 
  list(K = 3,
       N = nrow(recr_data.pila_validation),
       E = nrow(union_ecosubs),
       P = length(unique(recr_data.pila_validation$plot_id.i)),
       M = max(size_metadata$bin_id),
       X = 
         as.matrix(recr_data.pila_validation[,c('intercept')]),
       plot_id = recr_data.pila_validation$plot_id.i,
       ecosub_id = recr_data.pila_validation$ecosub.i,
       cprime = 
         recr_data.pila_validation %>% 
         group_by(plot_id, recruit_count) %>% 
         summarise() %>%
         ungroup() %>%
         arrange(plot_id) %>%
         pull(recruit_count),
       a = size_metadata$plot_area_ac[1],
       n = matrix(ncol = length(unique(recr_data.pila_validation$plot_id.i)),
                  nrow = nrow(size_metadata),
                  data = recr_data.pila_validation$tpa_unadj.init,
                  byrow = FALSE))



#### write results #############################################################

saveRDS(pila_mort_training, 
        here::here('02-data', '02-for_analysis', 'pila_mort_training_size.rds'))

saveRDS(pila_mort_validation,
        here::here('02-data', '02-for_analysis', 'pila_mort_validation_size.rds'))


saveRDS(pila_growth_training, 
        here::here('02-data', '02-for_analysis', 'pila_growth_training_size.rds'))

saveRDS(pila_growth_validation,
        here::here('02-data', '02-for_analysis', 'pila_growth_validation_size.rds'))


saveRDS(pila_fecd_training, 
        here::here('02-data', '02-for_analysis', 'pila_fecd_training_size.rds'))

saveRDS(pila_fecd_validation,
        here::here('02-data', '02-for_analysis', 'pila_fecd_validation_size.rds'))


# these are useful for linking back up with external data later:
saveRDS(union_plots,
        here::here('02-data', '02-for_analysis', 'union_plots_size.rds'))
saveRDS(union_ecosubs,
        here::here('02-data', '02-for_analysis', 'union_ecosubs_size.rds'))

#### setup #####################################################################

# load packages
library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)

# for reproducibility, use the version of stan we used
set_cmdstan_path("C:/Users/dfoster/Documents/.cmdstanr/cmdstan-2.28.2")

# load model data
pila_surv_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'pila_mort_training_size.rds'))

pila_growth_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'pila_growth_training_size.rds'))

pila_fecd_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'pila_fecd_training_size.rds'))

#### survival ##################################################################

# build model
surv_model = cmdstan_model(here::here('03-analysis', 'surv_model.stan'))

# run sampler
surv_fit = 
  surv_model$sample(
    data = pila_surv_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results'),
    output_basename = 'pila_surv_size',
    seed = 112188,
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
surv_fit$save_object(here::here('02-data', '03-results', 'surv_fit_size.rds'))

# save posterior df
surv_posterior = as_draws_df(surv_fit$draws())

saveRDS(surv_posterior, here::here('02-data', '03-results', 'surv_post_size.rds'))

#### growth ####################################################################

library(posterior)
library(bayesplot)

# build model
growth_model = cmdstan_model(here::here('03-analysis', 'growth_model.stan'))

# run sampler
growth_fit = 
  growth_model$sample(
    data = pila_growth_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results'),
    output_basename = 'pila_growth_size',
    seed = 112188,
    adapt_delta = 0.9,
    init = 
      list(
        list('beta' = c(0, 1, 0)),
        list('beta' = c(0, 1, 0)),
        list('beta' = c(0, 1, 0)),
        list('beta' = c(0, 1, 0))
      )
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
growth_fit$save_object(here::here('02-data', '03-results', 'growth_fit_size.rds'))

# save posterior df
growth_posterior = as_draws_df(growth_fit$draws())

saveRDS(growth_posterior, here::here('02-data', '03-results', 'growth_post_size.rds'))

#### fecundity ###############################################################

# build model
fecd_model = cmdstan_model(here::here('03-analysis', 'fecd_model.stan'))

# run sampler
fecd_fit = 
  fecd_model$sample(
    data = pila_fecd_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results'),
    output_basename = 'pila_fecd_size',
    seed = 110819,
    adapt_delta = 0.8,
  )

# check summary
fecd_fit$summary()

# check diagnostics
fecd_fit$cmdstan_diagnose()


# posterior density plots
mcmc_dens_overlay(fecd_fit$draws(variables = c('beta')))

# save fitted model
fecd_fit$save_object(here::here('02-data', '03-results', 'fecd_fit_size.rds'))

# save posterior df
fecd_posterior = as_draws_df(fecd_fit$draws())

saveRDS(fecd_posterior, here::here('02-data', '03-results', 'fecd_post_size.rds'))

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

set.seed(110819)

#### survival ##################################################################

fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'surv_fit_size.rds'))

samples = as_draws_df(fitted_model$draws())

pila_training = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_mort_training_size.rds'))

# simulating one data set per draw

samples.beta_s = samples %>% select(contains('beta'))
samples.ecoEffect_s = samples %>% select(contains('effect_ecosub')) %>% as.data.frame()
samples.plotEffect_s = samples %>% select(contains('effect_plot')) %>% as.data.frame()


samples.plotEffect_s %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, group = parameter))+
  geom_line(stat = 'density', alpha = 0.4)+
  theme_minimal()

fitted_model$summary(variables = c('beta', 'sigma_plot', 'sigma_ecosub'))

beta_s.med = 
  samples.beta_s %>%
  summarise_all(median) %>%
  as.numeric()

eco_effects = 
  samples.ecoEffect_s %>%
  summarise_all(median) %>%
  as.numeric()

plot_effects = 
  samples.plotEffect_s %>%
  summarise_all(median) %>%
  as.numeric()

logitp.med = 
  as.numeric(pila_training$X %*% beta_s.med)+
  as.numeric(eco_effects)[pila_training$ecosub_id]+
  as.numeric(plot_effects)[pila_training$plot_id]

median_results = 
  pila_training$X %>%
  as_tibble(x = .)

median_results$logitp = logitp.med
median_results$p = boot::inv.logit(logitp.med)
median_results$surv_true = pila_training$surv

median_results = 
  median_results %>%
  mutate(r = dense_rank(p),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20)))

median_results %>%
  ggplot()+
  geom_point(aes(x = r, y = p))+
  geom_jitter(aes(x = r, y = surv_true),
              color = 'red',
              size = 0,
              height = 0.1, width = 0)+
  geom_point(data = 
               median_results %>%
               group_by(r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             color = 'blue',
             aes(x = r, y = surv_true))+
  theme_minimal()

intermediate = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   
                   eco_effects = 
                     rnorm(n = pila_training$E,
                           mean = 0,
                           sd = samples$sigma_ecosub[i])
                   
                   plot_effects = 
                     rnorm(n = pila_training$P,
                           mean = 0,
                           sd = samples$sigma_plot)
                   
                   XB = as.numeric(pila_training$X %*% beta_s)
                   realized_eco = as.numeric(samples.ecoEffect_s[i,])[pila_training$ecosub_id]
                   realized_plot = as.numeric(samples.plotEffect_s[i,])[pila_training$plot_id]
                   
                   logitp = XB  + realized_eco + realized_plot
                   as.numeric(pila_training$X %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_training$ecosub_id]+
                     as.numeric(samples.plotEffect_s[i,])[pila_training$plot_id]
                   eco_effects[pila_training$ecosub_id]+
                     plot_effects[pila_training$plot_id]
                   
                   
                   p = boot::inv.logit(logitp)
                   surv_sim = rbinom(n = pila_training$N, 
                                     size = 1, 
                                     prob = p)
                   
                   result = 
                     pila_training$X %>% 
                     as_tibble(x = .) %>%
                     rownames_to_column('tree_id')
                   result$XB = XB
                   result$realized_eco = realized_eco
                   result$realized_plot = realized_plot
                   result$logitp = logitp
                   result$p = p
                   result$surv_sim = surv_sim
                   result$surv_true = pila_training$surv
                   result$iter = i
                   return(result)
                 })) 


mort_retrodictions = 
  intermediate %>%
  group_by(tree_id, surv_true) %>%
  summarise(p.50 = quantile(p, 0.5),
            p.mean = mean(p),
            p.025 = quantile(p, 0.025),
            p.975 = quantile(p, 0.975)) %>%
  ungroup() %>%
  mutate(r = dense_rank(p.mean),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20)))

head(mort_retrodictions)
surv_retrodictions_plot = 
  ggplot(
    data = mort_retrodictions,
    aes(x = r, y = surv_true))+
  geom_jitter(height = 0.1, width = 0, size = 1, color = 'red')+
  theme_minimal()+
  geom_point(data = mort_retrodictions,
             aes(x = r, y = p.mean))+
  geom_ribbon(data = mort_retrodictions,
              aes(x = r, ymin = p.025, ymax = p.975, y = p.50),
              alpha = 0.2)+
  geom_point(data = 
               mort_retrodictions %>% 
               group_by(r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             aes(x = r, y = surv_true),
             color = 'blue', pch = 4)

surv_retrodictions_plot

# looks pretty good
ggsave(surv_retrodictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'retrodictions_s_size.png'),
       height = 4.5, width = 6, units = 'in')

#### growth ####################################################################


fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'growth_fit_size.rds'))

samples = as_draws_df(fitted_model$draws())

pila_training = readRDS(here::here('02-data',
                                   '02-for_analysis',
                                   'pila_growth_training_size.rds'))



samples.beta_g = samples %>% select(contains('beta'))
samples.ecoEffect_g = samples %>% select(contains('effect_ecosub')) %>% as.data.frame()
samples.plotEffect_g = samples %>% select(contains('effect_plot')) %>% as.data.frame()
growth_retrodictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_g = as.numeric(samples.beta_g[i,])
                   mu = 
                     as.numeric(pila_training$X %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_training$ecosub_id]+
                     as.numeric(samples.plotEffect_g[i,])[pila_training$plot_id]
                   size1_sim = truncnorm::rtruncnorm(n = pila_training$N,
                                                     a = 0.0254,
                                                     mean = mu,
                                                     sd = samples$sigma_epsilon[i])
                   
                   result = 
                     pila_training$X %>% 
                     as_tibble(x = .) %>%
                     rownames_to_column('tree_id')
                   result$mu = mu
                   result$size1_sim = size1_sim
                   result$size1_true = pila_training$size1
                   result$iter = i
                   return(result)
                 }))

growth_retrodictions_plot = 
  ggplot(data = 
           growth_retrodictions %>%
           group_by(tree_id,size1_true) %>%
           summarise(size1_sim.50 = quantile(size1_sim, 0.5),
                     size1_sim.025 = quantile(size1_sim, 0.025),
                     size1_sim.975 = quantile(size1_sim, 0.975)) %>%
           ungroup() %>%
           mutate(r = dense_rank(size1_sim.50)))+
  geom_point(aes(x = r, y = size1_true),size = 0)+
  theme_minimal()+
  geom_line(aes(x = r, y = size1_sim.50), color = 'blue')+
  geom_ribbon(aes(x = r, ymin = size1_sim.025, ymax = size1_sim.975, y = size1_sim.50),
              alpha = 0.2)+
  labs(x = 'rank (simulated size)', y = 'Size at remeasure')

growth_retrodictions_plot

# model is slightly overpredicting size of the smallest trees, missing some variation
# in size

growth_retrodictions_plot2 = 
  ggplot(data = 
           growth_retrodictions %>%
           group_by(tree_id,size1_true) %>%
           summarise(size1_sim.50 = quantile(size1_sim, 0.5),
                     size1_sim.025 = quantile(size1_sim, 0.025),
                     size1_sim.975 = quantile(size1_sim, 0.975)) %>%
           ungroup() %>%
           mutate(r = dense_rank(size1_sim.50)))+
  geom_point(aes(x = size1_sim.50, y = size1_true),size = 0)+
  theme_minimal()+
  geom_ribbon(aes(x = size1_sim.50, ymin = size1_sim.025, ymax = size1_sim.975, y = size1_sim.50),
              alpha = 0.2)+
  geom_abline(intercept = 0, slope = 1, color = 'blue')

growth_retrodictions_plot2

growth_retrodictions_plot2+
  coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.25))+
  scale_y_continuous(breaks = seq(from = 0, to = 0.254, by = 0.0254))+
  scale_x_continuous(breaks = seq(from = 0, to = 2.54, by = 0.0254))

ggsave(growth_retrodictions_plot2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'retrodictions_g_size.png'),
       height = 4.5, width = 6, units = 'in')

#### recruitment ###############################################################
#skipped, this is slow and exactly the same model/data

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

set.seed(110819)

#### survival ##################################################################

fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'surv_fit_size.rds'))

samples = as_draws_df(fitted_model$draws())

pila_validation = readRDS(here::here('02-data',
                                     '02-for_analysis',
                                     'pila_mort_validation_size.rds'))

# simulating one data set per draw

samples.beta_s = samples %>% select(contains('beta'))
samples.ecoEffect_s = samples %>% select(contains('effect_ecosub')) %>% as.data.frame()
samples.plotEffect_s = samples %>% select(contains('effect_plot')) %>% as.data.frame()


mort_predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_s = as.numeric(samples.beta_s[i,])
                   logitp = 
                     as.numeric(pila_validation$X %*% beta_s)+
                     as.numeric(samples.ecoEffect_s[i,])[pila_validation$ecosub_id]+
                     as.numeric(samples.plotEffect_s[i,])[pila_validation$plot_id]
                   p = boot::inv.logit(logitp)
                   surv_sim = rbinom(n = pila_validation$N, 
                                     size = 1, 
                                     prob = p)
                   
                   result = 
                     pila_validation$X %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$logitp = logitp
                   result$p = p
                   result$surv_sim = surv_sim
                   result$surv_true = pila_validation$surv
                   result$iter = i
                   return(result)
                 })) %>%
  group_by(tree_id, surv_true) %>%
  summarise(p.50 = quantile(p, 0.5),
            p.mean = mean(p),
            p.025 = quantile(p, 0.025),
            p.975 = quantile(p, 0.975)) %>%
  ungroup() %>%
  mutate(r = dense_rank(p.mean),
         r_bin = cut(r, breaks = seq(0, nrow(.)+1, length.out = 20)))

do.call('bind_rows',
        lapply(X = 1:nrow(samples),
               FUN = function(i){
                 beta_s = as.numeric(samples.beta_s[i,])
                 logitp = 
                   as.numeric(pila_validation$X %*% beta_s)+
                   as.numeric(samples.ecoEffect_s[i,])[pila_validation$ecosub_id]+
                   as.numeric(samples.plotEffect_s[i,])[pila_validation$plot_id]
                 p = boot::inv.logit(logitp)
                 surv_sim = rbinom(n = pila_validation$N, 
                                   size = 1, 
                                   prob = p)
                 
                 result = 
                   pila_validation$X %>% 
                   as_tibble() %>%
                   rownames_to_column('tree_id')
                 result$logitp = logitp
                 result$p = p
                 result$surv_sim = surv_sim
                 result$surv_true = pila_validation$surv
                 result$iter = i
                 return(result)
               })) %>%
  group_by(tree_id, dbh_m.init, surv_true) %>%
  summarise(p.50 = quantile(p, 0.5),
            p.mean = mean(p),
            p.025 = quantile(p, 0.025),
            p.975 = quantile(p, 0.975)) %>%
  ungroup() %>%
  mutate(dbh_class = cut(dbh_m.init, breaks = seq(from = 0, to = 2.54, by = 0.0254),
                         labels = FALSE)) %>%
  group_by(dbh_class) %>%
  summarise(p.real = mean(surv_true),
            p.50 = mean(p.50),
            p.mean = mean(p.mean),
            p.025 = mean(p.025),
            p.975 = mean(p.975)) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(dbh_class)))+
  geom_point(aes(y = p.real), color = 'red')+
  geom_point(aes(y = p.mean), color = 'black')+
  geom_errorbar(aes(ymin = p.025, ymax = p.975), width = 0.2)

head(mort_predictions)
surv_predictions_plot = 
  ggplot(
    data = mort_predictions,
    aes(x = r, y = surv_true))+
  geom_jitter(height = 0.1, width = 0, size = 1, color = 'red')+
  theme_minimal()+
  geom_point(data = mort_predictions,
             aes(x = r, y = p.mean))+
  geom_ribbon(data = mort_predictions,
              aes(x = r, ymin = p.025, ymax = p.975, y = p.50),
              alpha = 0.2)+
  geom_point(data = 
               mort_predictions %>% 
               group_by(r_bin) %>%
               summarise(surv_true = sum(surv_true)/n(),
                         r = mean(r)) %>%
               ungroup(),
             aes(x = r, y = surv_true),
             color = 'blue', pch = 4)


surv_predictions_plot
# looks pretty good
ggsave(surv_predictions_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'predictions_s_size.png'),
       height = 4.5, width = 6, units = 'in')

#### growth ####################################################################


fitted_model = readRDS(here::here('02-data',
                                  '03-results',
                                  'growth_fit_size.rds'))

samples = as_draws_df(fitted_model$draws())

pila_validation = readRDS(here::here('02-data',
                                     '02-for_analysis',
                                     'pila_growth_validation_size.rds'))



samples.beta_g = samples %>% select(contains('beta'))
samples.ecoEffect_g = samples %>% select(contains('effect_ecosub')) %>% as.data.frame()
samples.plotEffect_g = samples %>% select(contains('effect_plot')) %>% as.data.frame()
growth_predictions = 
  do.call('bind_rows',
          lapply(X = 1:nrow(samples),
                 FUN = function(i){
                   beta_g = as.numeric(samples.beta_g[i,])
                   mu = 
                     as.numeric(pila_validation$X %*% beta_g)+
                     as.numeric(samples.ecoEffect_g[i,])[pila_validation$ecosub_id]+
                     as.numeric(samples.plotEffect_g[i,])[pila_validation$plot_id]
                   size1_sim = truncnorm::rtruncnorm(n = pila_validation$N,
                                                     a = 0.0254,
                                                     mean = mu,
                                                     sd = samples$sigma_epsilon[i])
                   
                   result = 
                     pila_validation$X %>% 
                     as_tibble() %>%
                     rownames_to_column('tree_id')
                   result$mu = mu
                   result$size1_sim = size1_sim
                   result$size1_true = pila_validation$size1
                   result$iter = i
                   return(result)
                 }))

growth_predictions_plot = 
  ggplot(data = 
           growth_predictions %>%
           group_by(tree_id,size1_true) %>%
           summarise(size1_sim.50 = quantile(size1_sim, 0.5),
                     size1_sim.025 = quantile(size1_sim, 0.025),
                     size1_sim.975 = quantile(size1_sim, 0.975)) %>%
           ungroup() %>%
           mutate(r = dense_rank(size1_sim.50)))+
  geom_point(aes(x = r, y = size1_true),size = 0)+
  theme_minimal()+
  geom_line(aes(x = r, y = size1_sim.50), color = 'blue')+
  geom_ribbon(aes(x = r, ymin = size1_sim.025, ymax = size1_sim.975, y = size1_sim.50),
              alpha = 0.2)+
  labs(x = 'rank (simulated size)', y = 'Size at remeasure')

growth_predictions_plot

# model is slightly overpredicting size of the smallest trees, missing some variation
# in size

growth_predictions_plot2 = 
  ggplot(data = 
           growth_predictions %>%
           group_by(tree_id,size1_true) %>%
           summarise(size1_sim.50 = quantile(size1_sim, 0.5),
                     size1_sim.025 = quantile(size1_sim, 0.025),
                     size1_sim.975 = quantile(size1_sim, 0.975)) %>%
           ungroup() %>%
           mutate(r = dense_rank(size1_sim.50)))+
  geom_point(aes(x = size1_sim.50, y = size1_true),size = 0)+
  theme_minimal()+
  geom_ribbon(aes(x = size1_sim.50, ymin = size1_sim.025, ymax = size1_sim.975, y = size1_sim.50),
              alpha = 0.2)+
  geom_abline(intercept = 0, slope = 1, color = 'blue')

growth_predictions_plot2

growth_predictions_plot2+
  coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.25))+
  scale_y_continuous(breaks = seq(from = 0, to = 0.254, by = 0.0254))+
  scale_x_continuous(breaks = seq(from = 0, to = 2.54, by = 0.0254))

ggsave(growth_predictions_plot2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'predictions_g_size.png'),
       height = 4.5, width = 6, units = 'in')

#### IPM #######################################################################

library(here)
library(tidyverse)
library(posterior)
library(bayesplot)
library(data.table)

#### build transitions #########################################################

# load mcmc results
surv_post = readRDS(here::here('02-data', 
                               '03-results',
                               'surv_post_size.rds'))

growth_post = readRDS(here::here('02-data',
                                 '03-results',
                                 'growth_post_size.rds'))

fecd_post = readRDS(here::here('02-data',
                               '03-results',
                               'fecd_post_size.rds'))

# extract parameters
beta_s = 
  surv_post %>%
  select(contains('beta')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

beta_g = 
  growth_post %>%
  select(contains('beta')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

beta_f = 
  fecd_post %>%
  select(contains('beta')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

#plotEffect_s = 
#  surv_post %>%
#  select(contains('effect_plot')) %>%
#  summarise_all(median) %>%
#  as.data.frame() %>%
#  as.numeric()

#plotEffect_g = 
#  growth_post %>%
#  select(contains('effect_plot')) %>%
#  summarise_all(median) %>%
#  as.data.frame() %>%
#  as.numeric()


ecoEffect_s = 
  surv_post %>%
  select(contains('effect_ecosub')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()

ecoEffect_g = 
  growth_post %>%
  select(contains('effect_ecosub')) %>%
  summarise_all(median) %>%
  as.data.frame() %>%
  as.numeric()


sigmaEpsilon_g = 
  growth_post %>%
  summarise_all(median) %>%
  pull(sigma_epsilon) %>%
  as.numeric()

size_metadata = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.rds')) %>%
  # convert to metric
  mutate(bin_midpoint = bin_midpoint * 0.0254,
         bin_lower = bin_lower * 0.0254,
         bin_upper = bin_upper * 0.0254) 

head(size_metadata)

ecoregions = 
  readRDS(here::here('02-data', '02-for_analysis', 'union_ecosubs.rds'))


A_observed = 
    array(dim = list(nrow(size_metadata),
                     nrow(size_metadata),
                     nrow(ecoregions)),
          dimnames = list('class_to' = 1:nrow(size_metadata),
                          'class_from' = 1:nrow(size_metadata),
                          'ecoregion' = 1:nrow(ecoregions)),
          data = 
            sapply(X = 1:nrow(ecoregions),
                   FUN = function(eco){
                     
                     print(paste0('Working on ecoregion ', eco))
                     
                     # construct explanatory variable matrix for vital rate 
                     # functions for the current plot
                     X_g = 
                       ecoregions %>%
                       slice(eco) %>%
                       expand(intercept = 1,
                              dbh = size_metadata$bin_midpoint) %>%
                       mutate(dbh2 = dbh**2) %>%
                       select(intercept, dbh, dbh2) %>%
                       as.matrix()
                     
                     X_s = 
                       ecoregions %>%
                       slice(eco) %>%
                       expand(intercept = 1,
                              dbh = size_metadata$bin_midpoint) %>%
                       mutate(dbh2 = dbh**2) %>%
                       select(intercept, dbh, dbh2) %>%
                       as.matrix()
                     
                     X_f = 
                       ecoregions %>%
                       slice(eco) %>%
                       expand(intercept = 1,
                              dbh = size_metadata$bin_midpoint) %>%
                       select(intercept) %>%
                       as.matrix()
                     
                     # calculate vector of survival probabilities for each 
                     # size class on this plot with this parameter draw
                     p = 
                       boot::inv.logit(as.numeric(X_s %*% beta_s) +
                                         ecoEffect_s[ecoregions$ecosub.i[eco]])
                     
                     # calculate vector of mean size at time 2 for each size 
                     # class on this plot with this parameter draw
                     mu = as.numeric(X_g %*% beta_g)+
                       ecoEffect_g[ecoregions$ecosub.i[eco]]
                     
                     # calculate vector of fecundity for each size class on this 
                     # plot with this parameter draw
                     f = 
                       exp(as.numeric(X_f %*% beta_f))
                     
                     # loop over each "from" size class
                     sapply(X = 1:nrow(size_metadata),
                            FUN = function(class_from){
                              
                              # growth kernel from this size class into 
                              # each other size class, using the cumulative 
                              # density function as recommended by Doak et al. 2021
                              g = 
                                ((pnorm(size_metadata$bin_upper,
                                        mu[class_from],
                                        sigmaEpsilon_g) - 
                                    pnorm(size_metadata$bin_lower,
                                          mu[class_from],
                                          sigmaEpsilon_g))/
                                   (1-pnorm(0.0254,
                                            mu[class_from],
                                            sigmaEpsilon_g)))
                              
                              # loop over every destination size class; 
                              # now vectorized
                              # calculate the transition kernel between the 
                              # current 'from' class and every 'to' class
                              transition_kern = 
                                # survival of each from class
                                (p[class_from] *
                                   # probability of growth from this class to 
                                   # every other class
                                   g)+(
                                     # total number of new recruits from this class 
                                     f[class_from]*
                                       
                                       # proportion of new recruits falling in the 
                                       # to class
                                       size_metadata$r)
                              
                              return(transition_kern)
                              
                            })
                     
                   }))

# save this so we don't have to rebuild it every time
saveRDS(A_observed,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'A_observed_size_eco.rds'))

A_observed = readRDS(here::here('02-data',
                                '03-results',
                                'real_fits',
                                'A_observed_size_eco.rds'))


transitions.df = 
  as.data.table(A_observed) %>% as_tibble()

head(transitions.df)

transitions.df %>%
  group_by(class_to, class_from) %>%
  summarise(value = median(value)) %>%
  ungroup() %>%
  left_join(size_metadata %>%
              select(class_to = bin_id, 
                     class_to_midpoint = bin_midpoint) %>%
              mutate(class_to = as.character(class_to))) %>%
  left_join(size_metadata %>%
              select(class_from = bin_id, 
                     class_from_midpoint = bin_midpoint) %>%
              mutate(class_from = as.character(class_from))) %>%
  ggplot(aes(x = class_from_midpoint, y = class_to_midpoint, fill = value))+
  geom_tile()+
  scale_fill_viridis_c()

#### lambda ####################################################################


lambda_observed = 
  sapply(X = 1:nrow(ecoregions),
         FUN = function(eco){
           A_eco = A_observed[,,eco]
           lambda_eco = max(as.numeric(Re(eigen(A_eco)$values)))
           return(lambda_eco)
         })

summary(lambda_observed)

results = 
  ecoregions %>%
  bind_cols(tibble(lambda = lambda_observed))

library(sf)

ecosub.sf = 
  st_read(here::here('02-data',
                     '00-source',
                     'usfs',
                     'S_USA.EcomapSubsections',
                     'S_USA.EcomapSubsections.shp')) %>%
  st_make_valid() %>%
  select(MAP_UNIT_S, MAP_UNIT_N) %>%
  left_join(
    results,
    by = c('MAP_UNIT_S' = 'ecosubcd')
  )

head(ecosub.sf)

library(ggspatial)

palette_breaks = 
  c(-Inf, 0.9, 0.95, 0.99, 1.01, 1.05, 1.1, Inf)

palette_labels = 
  c('Less than 0.90', '0.90 to 0.95', '0.95 to 0.99', 
    '0.99 to 1.01','1.01 to 1.05', '1.05 to 1.10', 'Greater than 1.10')

palette_colors = 
  RColorBrewer::brewer.pal(n = length(palette_breaks)-1,
                           name = 'RdBu')
  
# spectral
#palette_colors = "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#FFFFBF" "#E6F598" "#ABDDA4" "#66C2A5" "#3288BD"

# RdBu
palette_colors = 
  c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#67A9CF", "#2166AC")

ecosub.sf$lambda_binned = 
  cut(ecosub.sf$lambda,
      breaks = palette_breaks,
      labels = palette_labels)

ecosub.sf$lambda_binned = 
  factor(ecosub.sf$lambda_binned,
         levels = rev(levels(ecosub.sf$lambda_binned)))

ecosub_lambda_map = 
  ggplot(data = ecosub.sf %>% filter(!is.na(lambda_binned)),
       aes(fill = lambda_binned))+
  geom_sf(data = 
            USAboundaries::us_states(),
          fill = NA)+
  geom_sf(color = 'black', linewidth = 0.1)+
  theme_minimal()+
  scale_fill_manual(values = rev(palette_colors), na.value = 'transparent', drop = FALSE)+
  labs(fill = 'Ecoregion lambda')+
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = 'white', color = 'black'))+
  coord_sf(xlim = c(390000, 1100000), ylim = c(3740000, 5020000),
           crs = "EPSG:26910")+
  annotation_scale()

ggsave(ecosub_lambda_map,
       filename = here::here('04-communication', 'figures', 'manuscript', 'ecosub_lambda_map.png'),
       height = 7, width = 4.5)


