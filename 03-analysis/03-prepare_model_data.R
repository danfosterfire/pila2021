
# load packages
library(here)
library(tidyverse)


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
  # construct data for size:stressor interactions
  mutate(plot_id.i = as.integer(factor(plot_id)),
         ecosub.i = as.integer(factor(ecosubcd)),
         dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled)


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
  # construct data for size:stressor interactions
  mutate(plot_id.i = as.integer(factor(plot_id)),
         ecosub.i = as.integer(factor(ecosubcd)),
         dbh_fire = dbh_in.init*fire,
         dbh_insects = dbh_in.init*insects,
         dbh_disease = dbh_in.init*disease,
         dbh_ba = dbh_in.init*ba_scaled,
         dbh_cwd90 = dbh_in.init*cwd_dep90_scaled,
         dbh_cwdmean = dbh_in.init*cwd_mean_scaled)

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

# a row for each size class, for 1 inch size bins from 0.5-9.5"; gives 
# metadata about each class including the midpoints, upper and lower bounds, 
# and the sampling area in the FIA design. Only includes the smallest 10 
# size bins, which are modeled as the response in the recruitment submodel.
size_metadata = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.rds'))

# a row for each unique combination of subplot:species:size class, for 1 
# inch size bins from 0.5-9.5"; "count" gives the number of untagged 
# trees on the subplot in the species:size bin at the remeasurement (ie 
# ingrowth)
untagged_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'untagged_data.rds'))


#### prepare model data ########################################################
pila_data = 
  list(N_s = nrow(mort_data.pila),
       P_s = length(unique(mort_data.pila$plot_id)),
       E_s = length(unique(mort_data.pila$ecosubcd)),
       surv = as.integer(mort_data.pila$survived),
       plotid_s = mort_data.pila$plot_id.i,
       ecosub_s = mort_data.pila$ecosub.i,
       X_s = 
         mort_data.pila[,c('intercept', 'dbh_in.init', 'fire', 'insects', 
                           'disease', 'ba_scaled', 'cwd_dep90_scaled', 
                           'cwd_mean_scaled', 'dbh_fire', 'dbh_insects',
                           'dbh_disease', 'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')],
       N_g = nrow(growth_data.pila),
       P_g = length(unique(growth_data.pila$plot_id)),
       E_g = length(unique(growth_data.pila$ecosubcd)),
       size1_g = growth_data.pila$dbh_in.re,
       plotid_g = growth_data.pila$plot_id.i,
       ecosub_g = growth_data.pila$ecosub.i,
       X_g = 
         growth_data.pila[,c('intercept', 'dbh_in.init', 'fire', 'insects', 
                             'disease', 'ba_scaled', 'cwd_dep90_scaled', 
                             'cwd_mean_scaled', 'dbh_fire', 'dbh_insects',
                             'dbh_disease', 'dbh_ba', 'dbh_cwd90', 'dbh_cwdmean')])


#### write results #############################################################

saveRDS(pila_data, 
        here::here('02-data', '02-for_analysis', 'pila_data.rds'))

