
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
<<<<<<< HEAD
  filter(inv_manual.init >= 2.0 & inv_manual.re >= 2.0) %>%
    tidyr::expand(nesting(plt_cn, prev_plt_cn, plot_id, elev_ft, lat, lon,
=======
  tidyr::expand(nesting(plt_cn, prev_plt_cn, plot_id, elev_ft, lat, lon,
>>>>>>> aa9ed2400f79a4e56bf84bf92919f936f2c934e6
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


<<<<<<< HEAD
summary(predict(tph_lme4_fit))

=======
pila_growth_training = 
  list(K = 12,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(growth_data.pila_training),
       size1 = growth_data.pila_training$dbh_m.re,
       X = 
         as.matrix(growth_data.pila_training[,c('intercept', 'dbh_m.init', 
                                           'fire', 
                                     'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                     'cwd_mean_scaled', 'dbh_fire', 
                                     'dbh_wpbr',
                                     'dbh_ba', 'dbh_cwd90',
                                     'dbh_cwdmean')]),
       plot_id = growth_data.pila_training$plot_id.i,
       ecosub_id = growth_data.pila_training$ecosub.i)


pila_growth_validation = 
  list(K = 12,
       P = nrow(union_plots),
       E = nrow(union_ecosubs),
       N = nrow(growth_data.pila_validation),
       size1 = growth_data.pila_validation$dbh_m.re,
       X = 
         as.matrix(growth_data.pila_validation[,c('intercept', 'dbh_m.init', 
                                                  'fire', 
                                                  'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                                  'cwd_mean_scaled', 'dbh_fire', 
                                                  'dbh_wpbr',
                                                  'dbh_ba', 'dbh_cwd90',
                                                  'dbh_cwdmean')]),
       plot_id = growth_data.pila_validation$plot_id.i,
       ecosub_id = growth_data.pila_validation$ecosub.i)


#### prepare fecudnity training and validation #################################

pila_fecd_training = 
  list(K = 3,
       N = nrow(recr_data.pila_training),
       E = nrow(union_ecosubs),
       P = length(unique(recr_data.pila_training$plot_id.i)),
       M = max(size_metadata$bin_id),
       X = 
          as.matrix(recr_data.pila_training[,c('intercept', 'dbh_m.init', 
                                               'dbh_m2.init')]),
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
          as.matrix(recr_data.pila_validation[,c('intercept', 'dbh_m.init', 
                                               'dbh_m2.init')]),
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


#### tph and BA ################################################################
>>>>>>> aa9ed2400f79a4e56bf84bf92919f936f2c934e6

plot_data %>%
  select(plot_id,
         `2001-2009` = pila_ba_m2ha.init,
         `2010-2019` = pila_ba_m2ha.re) %>%
  pivot_longer(cols = c(`2001-2009`, `2010-2019`),
               names_to = 'timestep',
               values_to = 'pila_ba_m2ha') %>%
  group_by(timestep) %>%
  summarise(mean_ba_m2ha = mean(pila_ba_m2ha),
            se_ba_m2ha = sd(pila_ba_m2ha)/sqrt(n()),
            count = n()) %>%
  ungroup()

basal_area_plot = 
  plot_data %>%
  select(plot_id,
         `2001-2009` = pila_ba_m2ha.init,
         `2010-2019` = pila_ba_m2ha.re) %>%
  pivot_longer(cols = c(`2001-2009`, `2010-2019`),
               names_to = 'timestep',
               values_to = 'pila_ba_m2ha') %>%
  group_by(timestep) %>%
  summarise(mean_ba_m2ha = mean(pila_ba_m2ha),
            se_ba_m2ha = sd(pila_ba_m2ha)/sqrt(n())) %>%
  ungroup() %>%
  ggplot(aes(x = timestep, y = mean_ba_m2ha))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = mean_ba_m2ha-se_ba_m2ha,
                    ymax = mean_ba_m2ha+se_ba_m2ha),
                width = 0.25)+
  theme_minimal()+
  labs(x = 'Period', y = 'Basal Area (m^2 / ha)')

basal_area_plot



sizedist_data.pila %>%
  group_by(plot_id) %>%
  summarise(tph.init = sum(tpa_unadj.init/0.404686),
            tph.re = sum(tpa_unadj.re/0.404686)) %>%
  ungroup() %>%
  summarise(mean_tph.init = mean(tph.init),
            se_tph.init = sd(tph.init)/sqrt(n()),
            mean_tph.re = mean(tph.re),
            se_tph.re = sd(tph.re)/sqrt(n()),
            count = n())

tph_plot = 
  sizedist_data.pila %>%
  group_by(plot_id) %>%
  summarise(tph.init = sum(tpa_unadj.init/0.404686),
            tph.re = sum(tpa_unadj.re/0.404686)) %>%
  ungroup() %>%
  pivot_longer(cols = c(tph.init, tph.re),
               names_to = 'Period',
               values_to = 'tph') %>%
  mutate(Period = ifelse(Period=='tph.init',
                         '2001-2009',
                         '2010-2019')) %>%
  group_by(Period) %>%
  summarise(mean_tph = mean(tph),
            se_tph = sd(tph)/sqrt(n())) %>%
  ungroup() %>%
  ggplot(aes(x = Period, y = mean_tph))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = mean_tph-(se_tph),
                    ymax = mean_tph+(se_tph)),
                width = 0.25)+
  theme_minimal()+
  labs(y = 'Stem Density (trees / ha)',
       x = 'Period')
  
tph_plot



abundance_plot = 
  cowplot::plot_grid(basal_area_plot, tph_plot, ncol = 2)

abundance_plot

ggsave(abundance_plot,
       filename = 
         here::here('04-communication',
                  'figures',
                  'manuscript',
                  'abundance.png'),
       height = 4, width = 6.5, units = 'in')



# size distribution of PILA
sizedist_data.pila %>%
  group_by(dbh_class) %>%
  summarise(tpa_unadj.init = mean(tpa_unadj.init),
            tpa_unadj.re = mean(tpa_unadj.re)) %>%
  ungroup() %>%
  mutate(total_tpa_unadj.init = 
           sizedist_data.pila %>%
              group_by(plt_cn) %>%
              summarise(tpa_unadj.init = sum(tpa_unadj.init),
                        tpa_unadj.re = sum(tpa_unadj.re)) %>%
              ungroup() %>%
              summarise(total_tpa_unadj.init = mean(tpa_unadj.init),
                        total_tpa_unadj.re = mean(tpa_unadj.re)) %>%
           pull(total_tpa_unadj.init),
         total_tpa_unadj.re = 
           sizedist_data.pila %>%
              group_by(plt_cn) %>%
              summarise(tpa_unadj.init = sum(tpa_unadj.init),
                        tpa_unadj.re = sum(tpa_unadj.re)) %>%
              ungroup() %>%
              summarise(total_tpa_unadj.init = mean(tpa_unadj.init),
                        total_tpa_unadj.re = mean(tpa_unadj.re)) %>%
           pull(total_tpa_unadj.re)) %>%
  mutate(proportion.init = tpa_unadj.init / total_tpa_unadj.init,
         proportion.re = tpa_unadj.re / total_tpa_unadj.re) %>%
  ggplot(aes(x = dbh_class, y = proportion.init))+
  geom_point()+
  scale_y_log10()


#### write results #############################################################

saveRDS(pila_mort_training, 
        here::here('02-data', '02-for_analysis', 'pila_mort_training.rds'))

saveRDS(pila_mort_validation,
        here::here('02-data', '02-for_analysis', 'pila_mort_validation.rds'))


saveRDS(pila_growth_training, 
        here::here('02-data', '02-for_analysis', 'pila_growth_training.rds'))

saveRDS(pila_growth_validation,
        here::here('02-data', '02-for_analysis', 'pila_growth_validation.rds'))


saveRDS(pila_fecd_training, 
        here::here('02-data', '02-for_analysis', 'pila_fecd_training.rds'))

saveRDS(pila_fecd_validation,
        here::here('02-data', '02-for_analysis', 'pila_fecd_validation.rds'))


# these are useful for linking back up with external data later:
saveRDS(union_plots,
        here::here('02-data', '02-for_analysis', 'union_plots.rds'))
saveRDS(union_ecosubs,
        here::here('02-data', '02-for_analysis', 'union_ecosubs.rds'))

