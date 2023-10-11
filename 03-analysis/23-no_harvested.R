
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


#### remove harvested plots ####################################################

harvested_plots = 
  mort_data.pila %>%
  filter(tree_status.re=='harvested') %>%
  pull(plot_id)

growth_data.pila = 
  growth_data.pila %>%
  filter(!is.element(plot_id, harvested_plots))


mort_data.pila = 
  mort_data.pila %>%
  filter(!is.element(plot_id, harvested_plots))

plot_data = 
  plot_data%>%
  filter(!is.element(plot_id, harvested_plots))

recr_data.pila = 
  recr_data.pila %>%
  filter(!is.element(plot_id, harvested_plots))

sizedist_data.pila = 
  sizedist_data.pila %>%
  filter(!is.element(plot_id, harvested_plots))
  


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
  list(K = 18,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(mort_data.pila_training),
       surv = as.integer(mort_data.pila_training$survived),
       X = 
         as.matrix(mort_data.pila_training[,c('intercept', 'dbh_m.init', 
                                              'dbh_m2.init', 'fire', 
                                              'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                              'cwd_mean_scaled', 'dbh_fire', 'dbh2_fire',
                                              'dbh_wpbr','dbh2_wpbr',
                                              'dbh_ba', 'dbh2_ba', 'dbh_cwd90','dbh2_cwd90',
                                              'dbh_cwdmean', 'dbh2_cwdmean')]),
       plot_id = mort_data.pila_training$plot_id.i,
       ecosub_id = mort_data.pila_training$ecosub.i)


pila_mort_validation = 
  list(K = 18,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(mort_data.pila_validation),
       surv = as.integer(mort_data.pila_validation$survived),
       X = 
         as.matrix(mort_data.pila_validation[,c('intercept', 'dbh_m.init', 
                                                'dbh_m2.init', 'fire', 
                                                'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                                'cwd_mean_scaled', 'dbh_fire', 'dbh2_fire',
                                                'dbh_wpbr','dbh2_wpbr',
                                                'dbh_ba', 'dbh2_ba', 'dbh_cwd90','dbh2_cwd90',
                                                'dbh_cwdmean', 'dbh2_cwdmean')]),
       plot_id = mort_data.pila_validation$plot_id.i,
       ecosub_id = mort_data.pila_validation$ecosub.i)



#### prepare growth training and validation ####################################


pila_growth_training = 
  list(K = 18,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(growth_data.pila_training),
       size1 = growth_data.pila_training$dbh_m.re,
       X = 
         as.matrix(growth_data.pila_training[,c('intercept', 'dbh_m.init', 
                                                'dbh_m2.init', 'fire', 
                                                'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                                'cwd_mean_scaled', 'dbh_fire', 'dbh2_fire',
                                                'dbh_wpbr','dbh2_wpbr',
                                                'dbh_ba', 'dbh2_ba', 'dbh_cwd90','dbh2_cwd90',
                                                'dbh_cwdmean', 'dbh2_cwdmean')]),
       plot_id = growth_data.pila_training$plot_id.i,
       ecosub_id = growth_data.pila_training$ecosub.i)


pila_growth_validation = 
  list(K = 18,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(growth_data.pila_validation),
       size1 = growth_data.pila_validation$dbh_m.re,
       X = 
         as.matrix(growth_data.pila_validation[,c('intercept', 'dbh_m.init', 
                                                  'dbh_m2.init', 'fire', 
                                                  'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                                  'cwd_mean_scaled', 'dbh_fire', 'dbh2_fire',
                                                  'dbh_wpbr','dbh2_wpbr',
                                                  'dbh_ba', 'dbh2_ba', 'dbh_cwd90','dbh2_cwd90',
                                                  'dbh_cwdmean', 'dbh2_cwdmean')]),
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
        here::here('02-data', '02-for_analysis', 'pila_mort_training_nh.rds'))

saveRDS(pila_mort_validation,
        here::here('02-data', '02-for_analysis', 'pila_mort_validation_nh.rds'))


saveRDS(pila_growth_training, 
        here::here('02-data', '02-for_analysis', 'pila_growth_training_nh.rds'))

saveRDS(pila_growth_validation,
        here::here('02-data', '02-for_analysis', 'pila_growth_validation_nh.rds'))


saveRDS(pila_fecd_training, 
        here::here('02-data', '02-for_analysis', 'pila_fecd_training_nh.rds'))

saveRDS(pila_fecd_validation,
        here::here('02-data', '02-for_analysis', 'pila_fecd_validation_nh.rds'))


# these are useful for linking back up with external data later:
saveRDS(union_plots,
        here::here('02-data', '02-for_analysis', 'union_plots_nh.rds'))
saveRDS(union_ecosubs,
        here::here('02-data', '02-for_analysis', 'union_ecosubs_nh.rds'))

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
                     'pila_mort_training_nh.rds'))

pila_growth_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'pila_growth_training_nh.rds'))

pila_fecd_training = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'pila_fecd_training_nh.rds'))

#### survival ##################################################################

# build model
surv_model = cmdstan_model(here::here('03-analysis', 'surv_model.stan'))

# run sampler
surv_fit = 
  surv_model$sample(
    data = pila_surv_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results'),
    output_basename = 'pila_surv_nh',
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
surv_fit$save_object(here::here('02-data', '03-results', 'surv_fit_nh.rds'))

# save posterior df
surv_posterior = as_draws_df(surv_fit$draws())

saveRDS(surv_posterior, here::here('02-data', '03-results', 'surv_post_nh.rds'))

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
    output_basename = 'pila_growth_nh',
    seed = 112188,
    adapt_delta = 0.9,
    init = 
      list(
        list('beta' = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
        list('beta' = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
        list('beta' = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
        list('beta' = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
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
growth_fit$save_object(here::here('02-data', '03-results', 'growth_fit_nh.rds'))

# save posterior df
growth_posterior = as_draws_df(growth_fit$draws())

saveRDS(growth_posterior, here::here('02-data', '03-results', 'growth_post_nh.rds'))

#### fecundity ###############################################################

# build model
fecd_model = cmdstan_model(here::here('03-analysis', 'fecd_model.stan'))

# run sampler
fecd_fit = 
  fecd_model$sample(
    data = pila_fecd_training,
    parallel_chains = 4,
    output_dir = here::here('02-data', '03-results'),
    output_basename = 'pila_fecd_nh',
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
fecd_fit$save_object(here::here('02-data', '03-results', 'fecd_fit_nh.rds'))

# save posterior df
fecd_posterior = as_draws_df(fecd_fit$draws())

saveRDS(fecd_posterior, here::here('02-data', '03-results', 'fecd_post_nh.rds'))

library(here)
library(tidyverse)
library(posterior)
library(bayesplot)


#### setup #####################################################################


hypothetical_plots = 
  data.frame(intercept = rep(1, times = 9),
             fire = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
             wpbr = c(FALSE, FALSE, T, F, F, F, F, F, F),
             ba_scaled = c(0, 0, 0, -1, 1, 0, 0, 0, 0),
             cwd_dep90_scaled = c(0, 0, 0, 0, 0, -1, 1, 0, 0),
             cwd_mean_scaled = c(0, 0, 0, 0, 0, 0, 0, -1, 1),
             subp_id = 1:9,
             name = c('Undisturbed', 'Fire', 'WPBR', 'Low BA', 'High BA',
                      'Low Drought', 'High Drought', 'Wet Site', 'Dry Site'))

# load mcmc results
posterior_s = readRDS(here::here('02-data',
                                 '03-results',
                                 'surv_post_nh.rds'))

posterior_g = readRDS(here::here('02-data',
                                 '03-results',
                                 'growth_post_nh.rds'))

posterior_f = readRDS(here::here('02-data',
                                 '03-results',
                                 'fecd_post_nh.rds'))




size_metadata = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.rds')) %>%
  # convert to metric
  mutate(bin_midpoint = bin_midpoint * 0.0254,
         bin_lower = bin_lower * 0.0254,
         bin_upper = bin_upper * 0.0254) 


#### transitions ###############################################################


A_hypotheticals = 
  array(dim = c(nrow(size_metadata), # sizeclass to
                nrow(size_metadata), # sizeclass from
                nrow(hypothetical_plots), # plots
                nrow(posterior_s)), # posterior draws
        dimnames = list('class_to' = 1:nrow(size_metadata),
                        'class_from' = 1:nrow(size_metadata),
                        'plot' = 1:nrow(hypothetical_plots),
                        'draw' = 1:nrow(posterior_s)),
        data = 
          sapply(X = 1:nrow(posterior_s),
                 FUN = function(draw){
                   
                   print(paste0('Working on draw ', draw))
                   
                   # get beta_s for the current draw
                   beta_s = 
                     posterior_s %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   beta_g = 
                     posterior_g %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   beta_f = 
                     posterior_f %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaEpsilon_g = 
                     posterior_g %>%
                     slice(draw) %>%
                     pull(sigma_epsilon) %>%
                     as.numeric()
                   
                   sapply(X = 1:nrow(hypothetical_plots),
                          FUN = function(plot){
                            
                            # construct explanatory variable matrix for survival 
                            # for the current plot
                            X_g = 
                              hypothetical_plots %>%
                              slice(plot) %>%
                              expand(nesting(intercept, fire, wpbr, ba_scaled,
                                             cwd_dep90_scaled,cwd_mean_scaled),
                                     dbh = size_metadata$bin_midpoint) %>%
                              mutate(dbh2 = dbh**2,
                                     dbh_fire = dbh*fire,
                                     dbh2_fire = dbh2*fire,
                                     dbh_wpbr = dbh*wpbr,
                                     dbh2_wpbr = dbh2*wpbr,
                                     dbh_ba = dbh*ba_scaled,
                                     dbh2_ba = dbh2*ba_scaled,
                                     dbh_cwd_dep90 = dbh*cwd_dep90_scaled,
                                     dbh2_cwd_dep90 = dbh2*cwd_dep90_scaled,
                                     dbh_cwd_mean = dbh*cwd_mean_scaled,
                                     dbh2_cwd_mean = dbh2*cwd_mean_scaled) %>%
                              select(intercept, dbh, dbh2, fire, wpbr, ba_scaled,
                                     cwd_dep90_scaled, cwd_mean_scaled, 
                                     dbh_fire, dbh2_fire, dbh_wpbr, dbh2_wpbr, 
                                     dbh_ba, dbh2_ba,
                                     dbh_cwd_dep90, dbh2_cwd_dep90,
                                     dbh_cwd_mean, dbh2_cwd_mean) %>%
                              as.matrix()
                            
                            
                            X_s = 
                              hypothetical_plots %>%
                              slice(plot) %>%
                              expand(nesting(intercept, fire, wpbr, ba_scaled,
                                             cwd_dep90_scaled,cwd_mean_scaled),
                                     dbh = size_metadata$bin_midpoint) %>%
                              mutate(dbh2 = dbh**2,
                                     dbh_fire = dbh*fire,
                                     dbh2_fire = dbh2*fire,
                                     dbh_wpbr = dbh*wpbr,
                                     dbh2_wpbr = dbh2*wpbr,
                                     dbh_ba = dbh*ba_scaled,
                                     dbh2_ba = dbh2*ba_scaled,
                                     dbh_cwd_dep90 = dbh*cwd_dep90_scaled,
                                     dbh2_cwd_dep90 = dbh2*cwd_dep90_scaled,
                                     dbh_cwd_mean = dbh*cwd_mean_scaled,
                                     dbh2_cwd_mean = dbh2*cwd_mean_scaled) %>%
                              select(intercept, dbh, dbh2, fire, wpbr, ba_scaled,
                                     cwd_dep90_scaled, cwd_mean_scaled, 
                                     dbh_fire, dbh2_fire, dbh_wpbr, dbh2_wpbr, 
                                     dbh_ba, dbh2_ba,
                                     dbh_cwd_dep90, dbh2_cwd_dep90,
                                     dbh_cwd_mean, dbh2_cwd_mean) %>%
                              as.matrix()
                            
                            X_f = 
                              hypothetical_plots %>%
                              slice(plot) %>%
                              expand(nesting(intercept, fire, wpbr, ba_scaled,
                                             cwd_dep90_scaled,cwd_mean_scaled),
                                     dbh = size_metadata$bin_midpoint) %>%
                              mutate(dbh2 = dbh**2,
                                     dbh_fire = dbh*fire,
                                     dbh2_fire = dbh2*fire,
                                     dbh_wpbr = dbh*wpbr,
                                     dbh2_wpbr = dbh2*wpbr,
                                     dbh_ba = dbh*ba_scaled,
                                     dbh2_ba = dbh2*ba_scaled,
                                     dbh_cwd_dep90 = dbh*cwd_dep90_scaled,
                                     dbh2_cwd_dep90 = dbh2*cwd_dep90_scaled,
                                     dbh_cwd_mean = dbh*cwd_mean_scaled,
                                     dbh2_cwd_mean = dbh2*cwd_mean_scaled) %>%
                              select(intercept) %>%
                              as.matrix()
                            
                            
                            # calculate size_from length vector of survival 
                            # probabilities on this plot with this parameter draw
                            p = boot::inv.logit(as.numeric(X_s %*% beta_s))
                            
                            mu = as.numeric(X_g %*% beta_g)
                            
                            f = exp(as.numeric(X_f %*% beta_f))
                            
                            sapply(X = 1:nrow(size_metadata),
                                   FUN = function(class_from){
                                     
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
                                     
                                     transition_prob = 
                                       # survival of each from class
                                       (p[class_from] *
                                          # prob of growth from to
                                          g) +
                                       # number of new recruits
                                       (f[class_from] *
                                          size_metadata$r)
                                     
                                     return(transition_prob)
                                   })
                          })
                 }))


saveRDS(A_hypotheticals,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'hypothetical_As_nh.rds'))

#### lambda ####################################################################

A_hypotheticals = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'hypothetical_As_nh.rds'))

hypothetical_lambdas_full = 
  hypothetical_plots %>%
  expand(nesting(subp_id, name, intercept, fire, wpbr, ba_scaled, 
                 cwd_dep90_scaled, cwd_mean_scaled),
         data.frame(draw = 1:4000))


hypothetical_lambdas_full$lambda = 
  sapply(X = 1:nrow(hypothetical_plots),
         FUN = function(plot){
           sapply(X = 1:4000,
                  FUN = function(draw){
                    # paste0('s:',plot,'d:',draw) for testing
                    max(as.numeric(Re(eigen(A_hypotheticals[,,plot,draw])$values)))
                  })
         }) %>%
  as.numeric()


pretty_names = 
  hypothetical_plots$name
names(pretty_names) = hypothetical_plots$subp_id

hypothetical_lambdas_plot = 
  ggplot(data = 
           hypothetical_lambdas_full,
         aes(x = lambda))+
  geom_density(lwd = 1, fill = 'lightgrey')+
  geom_vline(xintercept = 1, color = 'red', lty = 2, lwd = 1)+
  theme_minimal()+
  facet_grid(subp_id~., scales = 'free_y',
             labeller = labeller(subp_id = pretty_names))+
  theme(axis.text.y = element_blank())+
  labs(x = 'Asymptotic Population Growth Rate', y = 'Posterior Density')

hypothetical_lambdas_plot

ggsave(hypothetical_lambdas_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'hypotheticals_lambda_post_nh.png'),
       height = 7.5, width = 4, units = 'in')


hypothetical_lambda_summary = 
  hypothetical_lambdas_full %>%
  group_by(subp_id, name) %>%
  summarise(lambda.med = median(lambda),
            lambda.05 = quantile(lambda, probs = 0.05),
            lambda.95 = quantile(lambda, probs = 0.95))

hypothetical_lambda_summary

write.csv(hypothetical_lambda_summary,
          here::here('04-communication',
                     'tables',
                     'hypothetical_lambdas_summary_nh.csv'))


