
# load packages
library(here)
library(tidyverse)

set.seed(110819)

# use this to control how many bins new recruits can fall into
max_recr_class = 1

# a row for every subplot, columns for the subplot-level covariates like 
# disturbance data, basal area, and drought
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
  
  #  scale height to ~ the 0-1 range
  mutate(
    height_scaled.init = height_ft.init / 270,
    height_scaled.re = height_ft.re / 270,
    height2 = height_scaled.init**2,
    height_scaled.re = height_ft.re / 270,
    height_fire = height_scaled.init*fire,
    height2_fire = height2 * fire,
    height_wpbr = height_scaled.init*wpbr,
    height2_wpbr = height2*wpbr,
    height_ba = height_scaled.init*ba_scaled,
    height2_ba = height2*ba_scaled,
    height_cwd90 = height_scaled.init*cwd_dep90_scaled,
    height2_cwd90 = height2*cwd_dep90_scaled,
    height_cwdmean = height_scaled.init*cwd_mean_scaled,
    height2_cwdmean = height2*cwd_mean_scaled)


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
  
  #  scale height to ~ the 0-1 range
  mutate(
    height_scaled.init = height_ft.init / 270,
    height_scaled.re = height_ft.re / 270,
    height2 = height_scaled.init**2,
    height_scaled.re = height_ft.re / 270,
    height_fire = height_scaled.init*fire,
    height2_fire = height2 * fire,
    height_wpbr = height_scaled.init*wpbr,
    height2_wpbr = height2*wpbr,
    height_ba = height_scaled.init*ba_scaled,
    height2_ba = height2*ba_scaled,
    height_cwd90 = height_scaled.init*cwd_dep90_scaled,
    height2_cwd90 = height2*cwd_dep90_scaled,
    height_cwdmean = height_scaled.init*cwd_mean_scaled,
    height2_cwdmean = height2*cwd_mean_scaled)

# a row for each unique combination of plot:species:size class, for 
# 1 inch size bins from 0.5-99.5"; "tpa_unadj.init" and "tpa_unadj.re" give 
# the area-adjusted density of stems in each species:size 
# bin for each plot at the initial and remeasurement. dbh_class gives the 
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

head(size_metadata)


mort_data.pila %>% filter(plot_id =='41-1-39-79537') %>% print(width = Inf)
growth_data.pila %>% filter(plot_id == '41-1-39-79537') %>% print(width = Inf)
sizedist_data.pila %>% filter(plot_id == '41-1-39-79537')
# a row for each size class and plot, with only plots which are in 
# BOTH the growth and mortality datasets
recr_data.pila  = 
  plot_data %>%
  filter(is.element(plot_id, mort_data.pila$plot_id)&
           is.element(plot_id, growth_data.pila$plot_id)) %>%
  # filter to only subplots where both initial and remeasurement had manual 
  # greater than or equal to 2
  filter(inv_manual.init >= 2.0 & inv_manual.re >= 2.0) %>%
    tidyr::expand(nesting(plt_cn, prev_plt_cn, plot_id, elev_ft, lat, lon,
                 ecosubcd, invdate.re, inv_manual.re, macro_break, fire,
                 insects, disease, cutting, invdate.init, inv_manual.init, 
                 ba_ft2ac, cwd_departure90, cwd_mean, ba_scaled, cwd_dep90_scaled,
                 cwd_mean_scaled, intercept, wpbr),
         height_ft.init = size_metadata$bin_midpoint) %>%
  mutate(height_class = cut(height_ft.init,
                         breaks = seq(from = 0, to= 270, by = 5),
                         labels = FALSE,
                         right = FALSE),
  #  scale height to ~ the 0-1 range
    height_scaled.init = height_ft.init / 270,
    height2 = height_scaled.init**2,
    height_fire = height_scaled.init*fire,
    height2_fire = height2 * fire,
    height_wpbr = height_scaled.init*wpbr,
    height2_wpbr = height2*wpbr,
    height_ba = height_scaled.init*ba_scaled,
    height2_ba = height2*ba_scaled,
    height_cwd90 = height_scaled.init*cwd_dep90_scaled,
    height2_cwd90 = height2*cwd_dep90_scaled,
    height_cwdmean = height_scaled.init*cwd_mean_scaled,
    height2_cwdmean = height2*cwd_mean_scaled) %>%
  # add in the observed counts 
  left_join(readRDS(here::here('02-data',
                               '01-preprocessed',
                               'untagged_data.rds')) %>%
              filter(species=='PILA')) %>%
  
  left_join(sizedist_data.pila %>%
              select(plot_id, species, height_class, tpa_unadj.init)) %>%
  
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
  filter(!is.element(plot_id,
                       group_by(., plot_id) %>%
                       summarise(tpa_unadj.init = sum(tpa_unadj.init)) %>%
                       ungroup() %>%
                       filter(tpa_unadj.init==0) %>%
                       pull(plot_id))) %>%
  
  # order by size and subplot
  arrange(plot_id, height_ft.init) %>%
  rename(untagged_count = count) 
  
recr_data.pila %>%
  pull(plot_id) %>%
  unique() %>%
  length()







# estimating the mean and variance parameters for the recruitment size kernel
# is not working well; with the normalization to sum to 1 the mean and variance 
# become unidentifiable, which I think is also causing the divergent transitions
# in the sampler. Because these aren't really the parameters of interest anyways, 
# I'm going to avoid trying to estimate them in the model and just supply the 
# recruitment size kernel as data. 
# Based on York et al. 2004, 90% of sugar pine seedlings fall in the 
# height range ~0.55-1.4m after 5 years of growth. New recruits have 
# on average 5 years to grow, and at most 10 years. Assume that  100% fall in the 
# DBH range 0-5" after 10 years of growth. 
r = c(1, rep(0, times = 53))
r
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
  filter(height_class <= max_recr_class) %>%
  arrange(plot_id, height_ft.init)

untagged_data.pila_validation = 
  recr_data.pila_validation %>%
  filter(height_class <= max_recr_class) %>%
  arrange(plot_id, height_ft.init)

#### prepare training and validation data ######################################



pila_training = 
  list(
    # number of fixef parameters
    K = 18,
    K_g = 18,
    # number of plots and ecoregions
    P = nrow(union_plots),
    E = nrow(union_ecosubs),
    
    # survival data
    N_s = nrow(mort_data.pila_training),
       surv = as.integer(mort_data.pila_training$survived),
       plotid_s = mort_data.pila_training$plot_id.i,
       ecosub_s = mort_data.pila_training$ecosub.i,
       X_s = 
         as.matrix(mort_data.pila_training[,c('intercept', 'height_scaled.init', 
                                           'height2', 'fire', 
                                     'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                     'cwd_mean_scaled', 'height_fire', 'height2_fire',
                                     'height_wpbr','height2_wpbr',
                                     'height_ba', 'height2_ba', 'height_cwd90','height2_cwd90',
                                     'height_cwdmean', 'height2_cwdmean')]),
    
    # growth data 
       N_g = nrow(growth_data.pila_training),
       size1_g = growth_data.pila_training$height_scaled.re,
       plotid_g = growth_data.pila_training$plot_id.i,
       ecosub_g = growth_data.pila_training$ecosub.i,
       X_g = 
         as.matrix(growth_data.pila_training[,c('intercept', 'height_scaled.init', 
                                                 'height2', 'fire', 
                                                 'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                                 'cwd_mean_scaled', 'height_fire', 'height2_fire',
                                                 'height_wpbr','height2_wpbr',
                                                 'height_ba', 'height2_ba', 'height_cwd90','height2_cwd90',
                                                 'height_cwdmean', 'height2_cwdmean')]),
    
    # recruitment data
    max_recr_class = max_recr_class,
    N_r = nrow(recr_data.pila_training),
    P_r = length(unique(recr_data.pila_training$plot_id)),
    X_r = 
      as.matrix(recr_data.pila_training[,c('intercept', 'height_scaled.init', 
                                           'height2', 'fire', 
                                           'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                           'cwd_mean_scaled', 'height_fire', 'height2_fire',
                                           'height_wpbr','height2_wpbr',
                                           'height_ba', 'height2_ba', 'height_cwd90','height2_cwd90',
                                           'height_cwdmean', 'height2_cwdmean')]),
    X_rg = 
      as.matrix(recr_data.pila_training[,c('intercept', 'height_scaled.init', 
                                           'height2', 'fire', 
                                           'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                           'cwd_mean_scaled', 'height_fire', 'height2_fire',
                                           'height_wpbr','height2_wpbr',
                                           'height_ba', 'height2_ba', 'height_cwd90','height2_cwd90',
                                           'height_cwdmean', 'height2_cwdmean')]),
    plotid_r = recr_data.pila_training$plot_id.i,
    ecosub_r = recr_data.pila_training$ecosub.i,
    M_r = nrow(size_metadata),
    u_bounds = size_metadata$bin_upper/270,
    l_bounds = size_metadata$bin_lower/270,
    a = size_metadata$plot_area_ac,
    cprime = matrix(ncol = length(unique(recr_data.pila_training$plot_id)),
                    nrow = max_recr_class,
                    data = untagged_data.pila_training$untagged_count,
                    byrow = FALSE),
    n = matrix(nrow = length(unique(recr_data.pila_training$plot_id)),
               ncol = nrow(size_metadata),
               data = recr_data.pila_training$tpa_unadj.init,
               byrow = TRUE),
    r = r)

pila_validation = 
  list(
    # number of fixef parameters
    K = 18,
    K_g = 18,
    # number of plots and ecoregions
    P = nrow(union_plots),
    E = nrow(union_ecosubs),
    
    # survival data
    N_s = nrow(mort_data.pila_validation),
    surv = as.integer(mort_data.pila_validation$survived),
    plotid_s = mort_data.pila_validation$plot_id.i,
    ecosub_s = mort_data.pila_validation$ecosub.i,
    X_s = 
      as.matrix(mort_data.pila_validation[,c('intercept', 'height_scaled.init', 
                                           'height2', 'fire', 
                                           'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                           'cwd_mean_scaled', 'height_fire', 'height2_fire',
                                           'height_wpbr','height2_wpbr',
                                           'height_ba', 'height2_ba', 'height_cwd90','height2_cwd90',
                                           'height_cwdmean', 'height2_cwdmean')]),
    
    # growth data 
    N_g = nrow(growth_data.pila_validation),
    size1_g = growth_data.pila_validation$height_scaled.re,
    plotid_g = growth_data.pila_validation$plot_id.i,
    ecosub_g = growth_data.pila_validation$ecosub.i,
    X_g = 
      as.matrix(growth_data.pila_validation[,c('intercept', 'height_scaled.init', 
                                             'height2', 'fire', 
                                             'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                             'cwd_mean_scaled', 'height_fire', 'height2_fire',
                                             'height_wpbr','height2_wpbr',
                                             'height_ba', 'height2_ba', 'height_cwd90','height2_cwd90',
                                             'height_cwdmean', 'height2_cwdmean')]),
    
    # recruitment data
    max_recr_class = max_recr_class,
    N_r = nrow(recr_data.pila_validation),
    P_r = length(unique(recr_data.pila_validation$plot_id)),
    X_r = 
      as.matrix(recr_data.pila_validation[,c('intercept', 'height_scaled.init', 
                                           'height2', 'fire', 
                                           'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                           'cwd_mean_scaled', 'height_fire', 'height2_fire',
                                           'height_wpbr','height2_wpbr',
                                           'height_ba', 'height2_ba', 'height_cwd90','height2_cwd90',
                                           'height_cwdmean', 'height2_cwdmean')]),
    X_rg = 
      as.matrix(recr_data.pila_validation[,c('intercept', 'height_scaled.init', 
                                           'height2', 'fire', 
                                           'wpbr','ba_scaled', 'cwd_dep90_scaled', 
                                           'cwd_mean_scaled', 'height_fire', 'height2_fire',
                                           'height_wpbr','height2_wpbr',
                                           'height_ba', 'height2_ba', 'height_cwd90','height2_cwd90',
                                           'height_cwdmean', 'height2_cwdmean')]),
    plotid_r = recr_data.pila_validation$plot_id.i,
    ecosub_r = recr_data.pila_validation$ecosub.i,
    M_r = nrow(size_metadata),
    u_bounds = size_metadata$bin_upper/270,
    l_bounds = size_metadata$bin_lower/270,
    a = size_metadata$plot_area_ac,
    cprime = matrix(ncol = length(unique(recr_data.pila_validation$plot_id)),
                    nrow = max_recr_class,
                    data = untagged_data.pila_validation$untagged_count,
                    byrow = FALSE),
    n = matrix(nrow = length(unique(recr_data.pila_validation$plot_id)),
               ncol = nrow(size_metadata),
               data = recr_data.pila_validation$tpa_unadj.init,
               byrow = TRUE),
    r = r)


#### tph and BA ################################################################

tph = 
  list(N = nrow(plot_data),
       E = length(unique(plot_data$ecosubcd)),
       Y = plot_data$pila_tph.re-plot_data$pila_tph.init,
       X = plot_data[,c('intercept', 'fire', 'wpbr', 'ba_scaled',
                        'cwd_dep90_scaled', 'cwd_mean_scaled')] %>%
         as.matrix(),
       ecosub_id = as.integer(factor(plot_data$ecosubcd)))

ba = 
  list(N = nrow(plot_data),
       E = length(unique(plot_data$ecosubcd)),
       Y = plot_data$pila_ba_m2ha.re-plot_data$pila_ba_m2ha.init,
       X = plot_data[,c('intercept', 'fire', 'wpbr', 'ba_scaled',
                        'cwd_dep90_scaled', 'cwd_mean_scaled')] %>%
         as.matrix(),
       ecosub_id = as.integer(factor(plot_data$ecosubcd)))

tph.df = 
  tph$X %>%
  as.data.frame() %>%
  mutate(y = tph$Y,
         ecosub = as.factor(tph$ecosub_id))

library(lme4)

tph_lme4_fit = lme4::lmer(data = tph.df,
                      formula = y ~ fire + wpbr + ba_scaled + cwd_dep90_scaled + 
                        cwd_mean_scaled + (1|ecosub))

summary(tph_lme4_fit)

summary(predict(tph_lme4_fit))


plot_data %>% summary()
plot_data %>%
  #filter(inv_manual.init >= 2.0 & inv_manual.re >= 2.0) %>%
  #filter(fire == FALSE) %>%
  select(plot_id,
         `2004-2009` = pila_ba_m2ha.init,
         `2014-2019` = pila_ba_m2ha.re) %>%
  pivot_longer(cols = c(`2004-2009`, `2014-2019`),
               names_to = 'timestep',
               values_to = 'pila_ba_m2ha') %>%
  group_by(timestep) %>%
  summarise(mean_ba_m2ha = mean(pila_ba_m2ha),
            se_ba_m2ha = sd(pila_ba_m2ha)/sqrt(n())) %>%
  ungroup()

basal_area_plot = 
  plot_data %>%
  filter(inv_manual.init >= 2.0 & inv_manual.re >= 2.0) %>%
  select(plot_id,
         `2004-2009` = pila_ba_m2ha.init,
         `2014-2019` = pila_ba_m2ha.re) %>%
  pivot_longer(cols = c(`2004-2009`, `2014-2019`),
               names_to = 'timestep',
               values_to = 'pila_ba_m2ha') %>%
  group_by(timestep) %>%
  summarise(mean_ba_m2ha = mean(pila_ba_m2ha),
            se_ba_m2ha = sd(pila_ba_m2ha)/sqrt(n())) %>%
  ungroup() %>%
  ggplot(aes(x = timestep, y = mean_ba_m2ha))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = mean_ba_m2ha-2*se_ba_m2ha,
                    ymax = mean_ba_m2ha+2*se_ba_m2ha),
                width = 0.25)+
  theme_minimal()+
  labs(x = 'Period', y = 'Basal Area (m^2/ha)')

basal_area_plot


tph_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'tph_data.rds'))

sizedist_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'sizedist_data.rds'))

sizedist_data %>%
  filter(species=='PILA') %>%
  group_by(plot_id) %>%
  summarise(tph.init = sum(tpa_unadj.init/0.404686),
            tph.re = sum(tpa_unadj.re/0.404686)) %>%
  ungroup() %>%
  summarise(mean_tph.init = mean(tph.init),
            se_tph.init = sd(tph.init)/sqrt(n()),
            mean_tph.re = mean(tph.re),
            se_tph.re = sd(tph.re)/sqrt(n()))

tph_plot = 
  sizedist_data %>%
  filter(species=='PILA') %>%
  group_by(plot_id) %>%
  summarise(tph.init = sum(tpa_unadj.init/0.404686),
            tph.re = sum(tpa_unadj.re/0.404686)) %>%
  ungroup() %>%
  pivot_longer(cols = c(tph.init, tph.re),
               names_to = 'Period',
               values_to = 'tph') %>%
  mutate(Period = ifelse(Period=='tph.init',
                         '2004-2009',
                         '2014-2019')) %>%
  group_by(Period) %>%
  summarise(mean_tph = mean(tph),
            se_tph = sd(tph)/sqrt(n())) %>%
  ungroup() %>%
  ggplot(aes(x = Period, y = mean_tph))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = mean_tph-(2*se_tph),
                    ymax = mean_tph+(2*se_tph)),
                width = 0.25)+
  theme_minimal()+
  labs(y = 'Stem Density (trees / ha)',
       x = 'Period')
  
  sizedist_data %>%
  filter(species=='PILA') %>%
  group_by(plot_id) %>%
  summarise(tph.init = sum(tpa_unadj.init/0.404686),
            tph.re = sum(tpa_unadj.re/0.404686)) %>%
  ungroup() %>%
  pivot_longer(cols = c(tph.init, tph.re),
               names_to = 'Period',
               values_to = 'tph') %>%
  mutate(Period = ifelse(Period=='tph.init',
                         '2004-2009',
                         '2014-2019')) %>%
  group_by(Period) %>%
  summarise(mean_tph = mean(tph),
            se_tph = sd(tph)/sqrt(n())) %>%
  ungroup()

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

# proportional and absolute changes in TPA
readRDS(here::here('02-data',
                   '01-preprocessed',
                   'tph_data.rds')) %>%
  filter(species=='PILA' & is.element(plot_id, plot_data$plot_id)) %>%
  group_by(plot_id) %>%
  summarise(pila_tph_unadj.init = sum(tpa_unadj.init)/0.404686,
            pila_tph_unadj.re = sum(tpa_unadj.re)/0.404686) %>%
  ungroup() %>%
  left_join(
    readRDS(here::here('02-data',
                     '01-preprocessed',
                     'tph_data.rds')) %>%
    group_by(plot_id) %>%
    summarise(total_tph_unadj.init = sum(tpa_unadj.init)/0.404686,
              total_tph_unadj.re = sum(tpa_unadj.re)/0.404686) %>%
    ungroup()
  ) %>%
  mutate(pila_prop.init = pila_tph_unadj.init / total_tph_unadj.init,
         pila_prop.re = pila_tph_unadj.re / total_tph_unadj.re) %>%
  summarise(mean_pila_tph_re = mean(pila_tph_unadj.re),
            mean_pila_tph_init = mean(pila_tph_unadj.init),
            mean_pila_prop_re = mean(pila_prop.re,na.rm = TRUE),
            mean_pila_prop_init = mean(pila_prop.init, na.rm = TRUE),
            se_pila_tph_re = sd(pila_tph_unadj.re)/sqrt(n()),
            se_pila_tph_init = sd(pila_tph_unadj.init)/sqrt(n()),
            se_pila_prop_re = sd(pila_prop.re,na.rm = TRUE)/sqrt(n()),
            se_pila_prop_init = sd(pila_prop.init, na.rm = TRUE)/sqrt(n()),
            )



# proportional and absolute changes in BA


# size distribution of PILA
sizedist_data.pila %>%
  group_by(height_class) %>%
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
  ggplot(aes(x = height_class, y = proportion.init))+
  geom_point()+
  scale_y_log10()

# proportion of plot total BA and TPH

#### write results #############################################################

saveRDS(pila_training, 
        here::here('02-data', '02-for_analysis', 'pila_training.rds'))

saveRDS(pila_validation,
        here::here('02-data', '02-for_analysis', 'pila_validation.rds'))

# these are useful for linking back up with external data later:
saveRDS(union_plots,
        here::here('02-data', '02-for_analysis', 'union_plots.rds'))
saveRDS(union_ecosubs,
        here::here('02-data', '02-for_analysis', 'union_ecosubs.rds'))

saveRDS(tph,
        here::here('02-data', 
                   '02-for_analysis',
                   'tph.rds'))

saveRDS(ba,
        here::here('02-data',
                   '02-for_analysis',
                   'ba.rds'))


#### scratch ###################################################################

# plots in growth data but not mortality data? why?
unique(growth_data.pila$plot_id)[!is.element(unique(growth_data.pila$plot_id),
                                             mort_data.pila$plot_id)] 
