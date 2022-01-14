
# load packages
library(here)
library(tidyverse)

set.seed(110819)

# a row for every subplot, columns for the subplot-level covariates like 
# disturbance data, basal area, and drought
subplot_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'subplot_data.rds')) %>%
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
  left_join(subplot_data %>%
              select(plot_id, subp_id, ecosubcd, 
                     intercept,
                     fire, insects, disease, 
                     ba_scaled, cwd_dep90_scaled, cwd_mean_scaled)) %>%
  
  # get DBH in meters, which will be on a nicer (close to 0-1, similar to other 
  # covariates) scale
  mutate(dbh_m.init = dbh_in.init*0.0254,
         dbh_m.re = dbh_in.re * 0.0254) %>%
  # construct data for size:stressor interactions
  mutate(dbh_fire = dbh_m.init*fire,
         dbh_insects = dbh_m.init*insects,
         dbh_disease = dbh_m.init*disease,
         dbh_ba = dbh_m.init*ba_scaled,
         dbh_cwd90 = dbh_m.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_m.init*cwd_mean_scaled)


# a row for every individual tagged tree which was alive at the initial 
# measurement, and columns indicating its survival status at remeasurement, 
# plus covariate columns for subplot-level data and interactions with size
mort_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'mort_data.rds')) %>%
  filter(species=='PILA')%>%
  left_join(subplot_data %>%
              select(plot_id, subp_id, ecosubcd,
                     intercept,
                     fire, insects, disease, 
                     ba_scaled, cwd_dep90_scaled, cwd_mean_scaled)) %>%
  
  # convert to meters scale
  mutate(dbh_m.init = dbh_in.init * 0.0254) %>%
  # construct data for size:stressor interactions
  mutate(dbh_fire = dbh_m.init*fire,
         dbh_insects = dbh_m.init*insects,
         dbh_disease = dbh_m.init*disease,
         dbh_ba = dbh_m.init*ba_scaled,
         dbh_cwd90 = dbh_m.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_m.init*cwd_mean_scaled)

# a row for each unique combination of subplot:species:size class, for 
# 1 inch size bins from 0.5-99.5"; "tpa_unadj.init" and "tpa_unadj.re" give 
# the area-adjusted density of stems in each species:size 
# bin for each subplot at the initial and remeasurement. dbh_class gives the 
# ID of each size bin (integer equal to the closed upper bound of the bin, so 
# stems 0"<=dbh<1" go in class 1, stems with dbh = 1.5" go in bin 2, etc.)
sizedist_data.pila = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'sizedist_data.rds')) %>%
  filter(species=='PILA')

# a row for each size class, for 1 inch size bins from 0.5-99.5"; gives 
# metadata about each class including the midpoints, upper and lower bounds, 
# and the sampling area in the FIA design.
size_metadata = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.rds'))

mort_data.pila %>% filter(subp_id =='41-1-39-79537-2') %>% print(width = Inf)
growth_data.pila %>% filter(subp_id == '41-1-39-79537-2') %>% print(width = Inf)
sizedist_data.pila %>% filter(subp_id == '41-1-39-79537-2')
# a row for each size class and subplot, with only subplots which are in 
# BOTH the growth and mortality datasets
recr_data.pila = 
  subplot_data %>%
  filter(is.element(subp_id, mort_data.pila$subp_id)&
           is.element(subp_id, growth_data.pila$subp_id)) %>%  
  # filter to only subplots where both initial and remeasurement had manual 
  # greater than or equal to 2
  filter(inv_manual.init >= 2.0 & inv_manual.re >= 2.0) %>%
  expand(nesting(plt_cn, prev_plt_cn, plot_id, subp_id, elev_ft, lat, lon,
                 ecosubcd, invdate.re, inv_manual.re, macro_break, fire,
                 insects, disease, cutting, invdate.init, inv_manual.init, 
                 ba_ft2ac, cwd_departure90, cwd_mean, ba_scaled, cwd_dep90_scaled,
                 cwd_mean_scaled, intercept),
         dbh_in.init = size_metadata$dbh_in.mean) %>%
  mutate(dbh_class = cut(dbh_in.init,
                         breaks = seq(from = 0, to= 100, by = 5),
                         labels = FALSE,
                         right = FALSE),
         dbh_m.init = dbh_in.init*0.0254,
         dbh_fire = dbh_m.init*fire,
         dbh_insects = dbh_m.init*insects,
         dbh_disease = dbh_m.init*disease,
         dbh_ba = dbh_m.init*ba_scaled,
         dbh_cwd90 = dbh_m.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_m.init*cwd_mean_scaled) %>%
  # add in the observed counts 
  left_join(readRDS(here::here('02-data',
                               '01-preprocessed',
                               'untagged_data.rds')) %>%
              filter(species=='PILA')) %>%
  
  left_join(sizedist_data.pila %>%
              select(subp_id, species, dbh_class, tpa_unadj.init)) %>%
  
  # there are 22 subplots which, despite being included in the mortality and 
  # growth data (ie, have PILA which are alive at initial and remeasurement), 
  # have 0 TPA of pila in any size class at initial measurement. Upon investigating,
  # I think these are all cases where the species ID of the tree was changed 
  # from something else to PILA, so the tree number shows alive at initial 
  # and remasurement and the species at remeasurement is PILA (ensuring the 
  # subplot is included in the growth and mortality data frames) but 0 trees 
  # in PILA size bins at initial measurement. solution is to just drop the 
  # offending subplots from the recruitment dataset, where the all 0s 
  # were screwing up the IPM model estimates for seedling density.
  filter(!is.element(subp_id,
                       group_by(., subp_id) %>%
                       summarise(tpa_unadj.init = sum(tpa_unadj.init)) %>%
                       ungroup() %>%
                       filter(tpa_unadj.init==0) %>%
                       pull(subp_id))) %>%
  
  # order by size and subplot
  arrange(subp_id, dbh_in.init) %>%
  rename(untagged_count = count)







# estimating the mean and variance parameters for the recruitment size kernel
# is not working well; with the normalization to sum to 1 the mean and variance 
# become unidentifiable, which I think is also causing the divergent transitions
# in the sampler. Because these aren't really the parameters of interest anyways, 
# I'm going to avoid trying to estimate them in the model and just supply the 
# recruitment size kernel as data. 
r = 
  # start with all the untagged data
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'untagged_data.rds')) %>%
  filter(species=='PILA'&is.element(subp_id, recr_data.pila$subp_id)) %>%
  
  # keep only the first two size classes (assume everything >10" dbh isn't 
  # really a new recruit, just a missed tree)
  filter(dbh_class <= 2) %>%
  
  # right join to the total number of new recruits on each subplot, keeping 
  # only subplots with at least 1 new recruit
  right_join(
      # get the total number of untagged trees in the first two size classes (new recxruits) 
    # on each subplot
    readRDS(here::here('02-data',
                       '01-preprocessed',
                       'untagged_data.rds')) %>%
      filter(species=='PILA'&is.element(subp_id, recr_data.pila$subp_id))  %>%
      filter(dbh_class <= 2) %>%
      group_by(subp_id) %>% 
      summarise(total_count = sum(count)) %>%
      ungroup() %>%
      # keep only subplots with at least 1 new recruit
      filter(total_count > 0)
  ) %>%
  
  # get the proportion of total new recruits in each of the first two size classes
  mutate(p_total = count / total_count) %>%
  
  # get the average distribution across all subplots
  group_by(dbh_class) %>%
  summarise(p_total = mean(p_total)) %>%
  ungroup() %>%
  
  pull(p_total)

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

# a row for each unique combination of subplot:species:size class, for 1 
# inch size bins from 0.5-9.5"; "count" gives the number of untagged 
# trees on the subplot in the species:size bin at the remeasurement (ie 
# ingrowth)

untagged_data.pila_training = 
  recr_data.pila_training %>%
  filter(dbh_class <= 2) %>%
  arrange(subp_id, dbh_in.init)

untagged_data.pila_validation = 
  recr_data.pila_validation %>%
  filter(dbh_class <= 2) %>%
  arrange(subp_id, dbh_in.init)

#### prepare training and validation data ######################################



pila_training = 
  list(
    # number of fixef parameters
    K = 12,
    # number of plots and ecoregions
    P = nrow(union_plots),
    E = nrow(union_ecosubs),
    
    # survival data
    N_s = nrow(mort_data.pila_training),
       surv = as.integer(mort_data.pila_training$survived),
       plotid_s = mort_data.pila_training$plot_id.i,
       ecosub_s = mort_data.pila_training$ecosub.i,
       X_s = 
         as.matrix(mort_data.pila_training[,c('intercept', 'dbh_m.init', 'fire', 
                                              'disease','ba_scaled', 'cwd_dep90_scaled', 
                                              'cwd_mean_scaled', 'dbh_fire', 'dbh_disease',
                                              'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')]),
    
    # growth data 
       N_g = nrow(growth_data.pila_training),
       size1_g = growth_data.pila_training$dbh_m.re,
       plotid_g = growth_data.pila_training$plot_id.i,
       ecosub_g = growth_data.pila_training$ecosub.i,
       X_g = 
         as.matrix(growth_data.pila_training[,c('intercept', 'dbh_m.init', 'fire', 
                                                'disease','ba_scaled', 'cwd_dep90_scaled', 
                                                'cwd_mean_scaled', 'dbh_fire', 'dbh_disease',
                                                'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')]),
    
    # recruitment data
    N_r = nrow(recr_data.pila_training),
    S_r = length(unique(recr_data.pila_training$subp_id)),
    X_r = 
      as.matrix(recr_data.pila_training[,c('intercept', 'dbh_m.init', 'fire', 
                                     'disease','ba_scaled', 'cwd_dep90_scaled', 
                                     'cwd_mean_scaled', 'dbh_fire', 
                                     'dbh_disease', 'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')]),
    plotid_r = recr_data.pila_training$plot_id.i,
    ecosub_r = recr_data.pila_training$ecosub.i,
    M_r = nrow(size_metadata),
    u_bounds = size_metadata$bin_upper*0.0254,
    l_bounds = size_metadata$bin_lower*0.0254,
    a = size_metadata$plot_area_ac[1:2],
    cprime = matrix(ncol = length(unique(recr_data.pila_training$subp_id)),
                    nrow = 2,
                    data = untagged_data.pila_training$untagged_count,
                    byrow = FALSE),
    n = matrix(nrow = length(unique(recr_data.pila_training$subp_id)),
               ncol = nrow(size_metadata),
               data = recr_data.pila_training$tpa_unadj.init,
               byrow = TRUE),
    r = r)

pila_validation = 
  list(
    # number of fixef parameters
    K = 12,
    # number of plots and ecoregions
    P = nrow(union_plots),
    E = nrow(union_ecosubs),
    
    # survival data
    N_s = nrow(mort_data.pila_validation),
       surv = as.integer(mort_data.pila_validation$survived),
       plotid_s = mort_data.pila_validation$plot_id.i,
       ecosub_s = mort_data.pila_validation$ecosub.i,
       X_s = 
         as.matrix(mort_data.pila_validation[,c('intercept', 'dbh_m.init', 'fire', 
                                     'disease','ba_scaled', 'cwd_dep90_scaled', 
                                     'cwd_mean_scaled', 'dbh_fire', 'dbh_disease',
                                     'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')]),
    
    # growth data 
       N_g = nrow(growth_data.pila_validation),
       size1_g = growth_data.pila_validation$dbh_m.re,
       plotid_g = growth_data.pila_validation$plot_id.i,
       ecosub_g = growth_data.pila_validation$ecosub.i,
       X_g = 
         as.matrix(growth_data.pila_validation[,c('intercept', 'dbh_m.init', 'fire', 
                                                  'disease','ba_scaled', 'cwd_dep90_scaled', 
                                                  'cwd_mean_scaled', 'dbh_fire', 'dbh_disease',
                                                  'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')]),
    
    # recruitment data
    N_r = nrow(recr_data.pila_validation),
    S_r = length(unique(recr_data.pila_validation$subp_id)),
    X_r = 
      as.matrix(recr_data.pila_validation[,c('intercept', 'dbh_m.init', 'fire', 
                                     'disease','ba_scaled', 'cwd_dep90_scaled', 
                                     'cwd_mean_scaled', 'dbh_fire', 
                                     'dbh_disease', 'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')]),
    plotid_r = recr_data.pila_validation$plot_id.i,
    ecosub_r = recr_data.pila_validation$ecosub.i,
    M_r = nrow(size_metadata),
    u_bounds = size_metadata$bin_upper*0.0254,
    l_bounds = size_metadata$bin_lower*0.0254,
    a = size_metadata$plot_area_ac[1:2],
    cprime = matrix(ncol = length(unique(recr_data.pila_validation$subp_id)),
                    nrow = 2,
                    data = untagged_data.pila_validation$untagged_count,
                    byrow = FALSE),
    n = matrix(nrow = length(unique(recr_data.pila_validation$subp_id)),
               ncol = nrow(size_metadata),
               data = recr_data.pila_validation$tpa_unadj.init,
               byrow = TRUE),
    r = r)


#### write results #############################################################

saveRDS(pila_training, 
        here::here('02-data', '02-for_analysis', 'pila_training.rds'))

saveRDS(pila_validation,
        here::here('02-data', '02-for_analysis', 'pila_validation.rds'))
