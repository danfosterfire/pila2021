
#### setup #####################################################################

library(here)
library(sf)
library(tidyverse)

treelist = readRDS(here::here('02-data',
                              '01-preprocessed',
                              'fia_treelist.rds'))

context = readRDS(here::here('02-data',
                             '01-preprocessed',
                             'fia_context.rds'))

seedlings = readRDS(here::here('02-data',
                               '01-preprocessed',
                               'fia_seedlings.rds'))

pila_range.sf = 
  st_read(here::here('02-data',
                     '01-preprocessed',
                     'pila_range_map.shp'))

#### quality control ###########################################################

# trees got duplicated at some point, check to see that tree CNs are unique 
# here, and that there are at most 2 instances of a tree_id
nrow(treelist) == length(unique(treelist$tre_cn))
treelist %>%
  group_by(tree_id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n>2)
# there are 205 instances where a tree appears more than twice; probably not a
# problem here because I still am including plots which were inventoried more than 
# twice


#### aggregate context data to subplots ########################################

# when a tree is present, I can tie a specific condition to a specific tree, 
# but this isn't possible when the tree is missing at either the initial or 
# followup. However, because subplots are constant over time (I thinK????)
# I can use a tree's subplot when it is present to ID it when it's not.

names(context)

unique(c(context$disturb1, context$disturb2, context$disturb3))
unique(c(context$treat1, context$treat2, context$treat3))
unique(c(context$treatpnw1, context$treatpnw2, context$treatpnw3))
subplots = 
  context %>%
  
  # treat 9999 disturb years as NA
  mutate(
    disturb1_yr = ifelse(disturb1_yr==9999,
                         NA,
                         disturb1_yr),
    disturb2_yr = ifelse(disturb2_yr==9999,
                         NA,
                         disturb2_yr),
    disturb3_yr = ifelse(disturb3_yr==9999,
                         NA,
                         disturb3_yr)
  ) %>%
  
  # treat NA disturbance and treatment codes as 'none'
  mutate(
    disturb1 = ifelse(is.na(disturb1), 'none', disturb1),
    disturb2 = ifelse(is.na(disturb2), 'none', disturb2),
    disturb3 = ifelse(is.na(disturb3), 'none', disturb3),
    treat1 = ifelse(is.na(treat1), 'none', treat1),
    treat2 = ifelse(is.na(treat2), 'none', treat2),
    treat3 = ifelse(is.na(treat3), 'none', treat3),
    treatpnw1 = ifelse(is.na(treatpnw1), 'none', treatpnw1),
    treatpnw2 = ifelse(is.na(treatpnw2), 'none', treatpnw2),
    treatpnw3 = ifelse(is.na(treatpnw3), 'none', treatpnw3),
  ) %>%
  
  # scan the disturbance columns and turn them into relevant logicals
  mutate(
    
    # first, scan disturb1, 2, and 3 separately for fire, and fill in fire year
    # disturb_yr is well covered for fire, not so much for others
    fire_1 = is.element(disturb1, c('fire', 'fire (crown)', 'fire (ground)')),
    crownfire_1 = disturb1=='fire (crown)',
    fireyr_1 = ifelse(fire_1, disturb1_yr, NA),
    fire_2 = is.element(disturb2, c('fire', 'fire (crown)', 'fire (ground)')),
    crownfire_2 = disturb2=='fire (crown)',
    fireyr_2 = ifelse(fire_2, disturb2_yr, NA),
    fire_3 = is.element(disturb3, c('fire', 'fire (crown)', 'fire (ground)')),
    crownfire_3 = disturb3=='fire (crown)',
    fireyr_3 = ifelse(fire_3, disturb3_yr, NA),
    
    fire = fire_1|fire_2|fire_3,
    fire_yr = ifelse(fire_1,
                     fireyr_1,
                     ifelse(fire_2,
                            fireyr_2,
                            ifelse(fire_3,
                                   fireyr_3,
                                   NA))),
    
    insects = 
      is.element(disturb1, c('insects (trees)', 'insects')) | 
      is.element(disturb2, c('insects (trees)', 'insects')) | 
      is.element(disturb3, c('insects (trees)', 'insects')),
    
    disease = 
      is.element(disturb1, c('disease (trees)', 'disease')) | 
      is.element(disturb2, c('disease (trees)', 'disease')) | 
      is.element(disturb3, c('disease (trees)', 'disease')),
    
    drought = 
      disturb1=='drought'|disturb2=='drought'|disturb3=='drought',
    
    suppression = 
      disturb1=='vegetation'|disturb2=='vegetation'|disturb3=='vegetation',
    
    other = 
      (is.element(disturb1, 
                  c('weather (other)', 'human (other)', 'unknown / other',
                    'animal', 'geologic', 'insects (understory)', 
                    'disease (understory)')) | treat1=='other'|treatpnw2=='other')| 
      (is.element(disturb2, 
                  c('weather (other)', 'human (other)', 'unknown / other',
                    'animal', 'geologic', 'insects (understory)', 
                    'disease (understory)')) | treat2=='other' | treatpnw2=='other')|
      (is.element(disturb3, 
                  c('weather (other)', 'human (other)', 'unknown / other',
                    'animal', 'geologic', 'insects (understory)', 
                    'disease (understory)')) | treat3=='other'|treatpnw3=='other'),
    
    cutting = 
      treat1=='cutting'|treat2=='cutting'|treat3=='cutting'|
      treatpnw1=='cutting'|treatpnw2=='cutting'|treatpnw3=='cutting',
    
    siteprep = 
      treat1=='siteprep'|treat2=='siteprep'|treat3=='siteprep',
    
    artificialregen = 
      treat1=='regen (artificial)'|treat2=='regen (artificial)'|treat3=='regen (artificial)',
    
    naturalregen = 
      treat1=='regen (natural)'|treat2=='regen (natural)'|treat3=='regen (natural)'
    ) %>%
  
  # aggregate the conds together to subplots
  group_by(plt_cn, state_id, plot_id, subp_id, inv_year_nominal,
           inv_date, inv_kind, inv_design, inv_manual, inv_macrobrk, inv_plotstatus,
           inv_nonsamp, inv_subpstatus, prev_plt_cn, 
           elev_ft, lat, lon, ecosub_cd, topographic) %>%
  summarise(
    fire = any(fire),
    fire_yr = max(fire_yr), 
    # if there were multiple fires, flag the oldest, 
    # because we'll set trees after that as burned
    insects = any(insects),
    disease = any(disease),
    drought = any(drought),
    suppression = any(suppression),
    cutting = any(cutting),
    siteprep = any(siteprep),
    artificialregen = any(artificialregen),
    naturalregen = any(naturalregen)
    
  ) %>%
  ungroup()




#### filtering subplots ########################################################


# aspatial filtering first
subplots = 
  subplots %>%
  
  # within CA or OR
  filter(is.element(state_id, c(6, 41))) %>%
  
  # invyr in correct range; drops the phase 3 "off subpanel"
  filter(inv_year_nominal >= 2000 & inv_year_nominal <= 2021)
  
  
# lat,lon within PILA range polygon
plt_cn_within_range = 
  
  # create spatial points
  st_as_sf(subplots,
           coords = c('lon', 'lat'),
           crs = 4269) %>%
  
  # intersect them with the range polygon
  st_intersection(st_transform(pila_range.sf, st_crs(.))) %>%
  
  # extract the plt_cns 
  pull(plt_cn)

subplots = 
  subplots %>%
  filter(is.element(plt_cn, plt_cn_within_range))


#### data exploration ##########################################################

names(subplots)

lapply(names(subplots)[c(c(7, 10), c(18:27))],
       FUN = function(response){
         ggplot(data = subplots,
                aes(x = as.character(inv_year_nominal),
                    fill = .data[[response]]))+
           geom_bar()+
           theme_minimal()
       })

ggplot(data = subplots,
       aes(x = as.character(inv_year_nominal),
           fill = is.na(prev_plt_cn)))+
  geom_bar()+
  theme_minimal()

#### joining other spatial data ################################################


#### pivot across time and join trees with subplot data ########################

# need to combine the trees data with the subplots data, and also pivot the 
# initial and remeasurement data wide; not sure whether to pivot the trees and 
# subplots data separately, and then join? or join them, then pivot? should be 
# easier to fill in the subplot data for missing trees, now that I've merged 
# it down to subplots instead of conds. Can compare matching by plt_cn:prev_plt_cn
# and tre_cn:prev_tre_cn against matching explititly on tree_id and subplot_id.
# hopefully this will resolve the missing trees enough to help untangle that 
# problem.

head(subplots)
head(treelist)

# pivot first, then join:
subplots_wide = 
  
  # start with the subplots table
  subplots %>%
  
  # select only the observations that were followup measurements made on or after 2011
  filter(inv_kind=='national_remeasure') %>%
  
  # left join with the matching observations
  left_join(subplots %>%
              filter(inv_kind != 'national_remeasure'),
            by = c('prev_plt_cn' = 'plt_cn',
                   'state_id' = 'state_id',
                   'plot_id' = 'plot_id', 
                   'subp_id' = 'subp_id'),
            suffix = c('.re', '.init'))

# these data points should not have changed
subplots_wide %>%
  filter((elev_ft.re!=elev_ft.init)|
           (lat.re!=lat.init)|
           (lon.re!=lon.init)|
           (ecosub_cd.re!=ecosub_cd.init)|
           (topographic.re!=topographic.init)) %>%
  dplyr::select(lat.re, lat.init, lon.re, lon.init,
         ecosub_cd.re, ecosub_cd.init, topographic.re, topographic.init) %>%
  print(n = Inf)

subplots_wide %>%
  filter(inv_design.re != inv_design.init)

# topographic position descriptions are not entirely consistent; the others 
# we can merge
# pivot first, then join:
subplots_wide = 
  
  # start with the subplots table
  subplots %>%
  
  # select only the observations that were followup measurements made on or after 2011
  filter(inv_kind=='national_remeasure') %>%
  
  # left join with the matching observations
  left_join(subplots %>%
              filter(inv_kind != 'national_remeasure'),
            by = c('prev_plt_cn' = 'plt_cn',
                   'state_id' = 'state_id',
                   'plot_id' = 'plot_id', 
                   'subp_id' = 'subp_id',
                   'elev_ft' = 'elev_ft',
                   'lat' = 'lat',
                   'lon' = 'lon',
                   'ecosub_cd' = 'ecosub_cd',
                   'inv_design' = 'inv_design'),
            suffix = c('.re', '.init'))

names(subplots_wide)

ggplot(data = subplots_wide,
       aes(x = as.character(inv_year_nominal.re),
           fill = as.character(insects.re)))+
  geom_bar()

ggplot(data = subplots_wide,
       aes(fill = is.na(prev_plt_cn.init),
           x = as.character(inv_year_nominal.init)))+
  geom_bar()


names(treelist)

treelist_wide = 
  treelist %>%
  
  # filter to only trees which are in a remeasurement subplots
  filter(is.element(plt_cn, subplots_wide$plt_cn)) %>%
  
  # full join with trees which are in an initial subplots; note that this 
  # includes trees which are ont 'initial' subplots which have a prev_plt_cn
  full_join(treelist %>%
              filter(is.element(plt_cn, subplots_wide$prev_plt_cn)),
            by = c('prev_tre_cn' = 'tre_cn'),
            suffix = c('.re', '.init'))

names(treelist_wide)


# these columns should match up too, or one of them should be NA
treelist_wide %>%
  filter(!(is.na(state_id.re)|is.na(state_id.init)|state_id.re==state_id.init) |
           !(is.na(plot_id.re)|is.na(plot_id.init)|plot_id.re==plot_id.init)|
           !(is.na(subp_id.re)|is.na(subp_id.init)|subp_id.re==subp_id.init)|
           !(is.na(tree_id.re)|is.na(tree_id.init)|tree_id.re==tree_id.init)) %>%
  dplyr::select(state_id.re, state_id.init, plot_id.re, plot_id.init,
         subp_id.re, subp_id.init, tree_id.re, tree_id.init,
         spp.re, spp.init, dbh_in.re, dbh_in.init) %>%
  print(n = Inf, width = Inf)

names(treelist_wide)

ggplot(data = treelist,
       aes(x = as.character(inv_year_nominal),
           fill = as.character(prev_status)))+
  geom_bar()
# these are data columns, not ID colums, which should match in theory but probably 
# wont always in reality, and i should just drop them and trust my tree_cn merging instead...
treelist_wide %>%
  mutate(prev_status.re = 
           ifelse(prev_status.re==1, 'live',
                  ifelse(prev_status.re==2, 'dead',
                         NA))) %>%
  filter(!(is.na(tree_status.init)|is.na(prev_status.re)|tree_status.init==prev_status.re)|
           !(is.na(prev_cond.re)|is.na(cond_id.init)|prev_cond.re==cond_id.init)) %>%
  dplyr::select(tree_status.init, tree_status.re, prev_status.init, prev_status.re,
         cond_id.init, cond_id.re, prev_cond.init, prev_cond.re) %>%
  print(width = Inf)

# huh they do always match (at least for trees with two observations)

# looks like there are a few exceptions (5 out of 140k); these are rare and 
# the trees themselves seem to match, so I'm going to treat these as typos 
# and trust the tre_cn=prev_tre_cn matching 
names(treelist)
combined_wide.present_at_remeasure = 
  
  # first handle the trees which were present at remeasurement 
  # (which may or may not have been present initially)
  subplots_wide %>%
  
  # join in the remeasurement trees, creating a table with 1 
  # row per tree; we are dropping rows for any subplots
  # which didn't have any remeasurement trees
  inner_join(treelist %>%
              dplyr::select(-state_id, -plot_id, -subp_id),
            by = c('plt_cn' = 'plt_cn',
                   # i checked and this always matches too, not losing any here
                   'inv_year_nominal.re' = 'inv_year_nominal'),
            suffix = c('.subp', '.tree')) %>%
    
    # now, join the initial trees which DO have a followup observation (ie, !is.na(prev_tre_cn)
    # for the followup observation currently in the table); this doesn't add 
    # any rows
    left_join(treelist %>%
                dplyr::select(-state_id, -plot_id, -subp_id),
              by = c('prev_plt_cn' = 'plt_cn',
                     'prev_tre_cn' = 'tre_cn',
                     'inv_year_nominal.init' = 'inv_year_nominal'),
              suffix = c('.re', '.init'))
  
names(combined_wide.present_at_remeasure)
names(subplots_wide)
combined_wide.not_present_at_remeasure = 

  # now, handle the trees which were present initially but did not have 
  # any followup data
  # again start with the full subplots table
  subplots_wide %>%
    
    # inner join it to the subset of the treelist which hasn't appeared already 
    # in the table; this keeps only subplots where some of the trees 
  # were not present during the remeasurement, and keeps only trees which were not present in the 
  # remeasurement, but are included in the set of valid subplot plt_cn s 
    inner_join(treelist %>%
                 dplyr::select(-state_id, -plot_id, -subp_id) %>%
                 filter(!is.element(tre_cn, 
                                    combined_wide.present_at_remeasure$tre_cn)&
                          !is.element(tre_cn, 
                                      combined_wide.present_at_remeasure$prev_tre_cn)),
               by = c('prev_plt_cn' = 'plt_cn',
                      'inv_year_nominal.init' = 'inv_year_nominal'),
               suffix = c('.subp', '.tree')) %>%
  # need to rename all the tree data columns as being '.init'
  dplyr::select(plt_cn, state_id, plot_id, subp_id,
         inv_year_nominal.re, inv_date.re, inv_kind.re, inv_design, 
         inv_macrobrk.re, inv_plotstatus.re, inv_nonsamp.re, inv_subpstatus.re,
         prev_plt_cn, elev_ft, lat, lon, ecosub_cd, topographic.re, fire.re, fire_yr.re,
         insects.re, disease.re, drought.re, suppression.re, cutting.re, 
         siteprep.re, artificialregen.re, naturalregen.re,
         inv_year_nominal.init, inv_date.init, inv_kind.init, 
         inv_macrobrk.init, inv_plotstatus.init, inv_nonsamp.init, inv_subpstatus.init,
         prev_plt_cn.init, topographic.init, fire.init, fire_yr.init, 
         insects.init, disease.init, drought.init, suppression.init,
         cutting.init, siteprep.init, artificialregen.init, naturalregen.init,
         prev_tre_cn = tre_cn, cond_id.init = cond_id, tree_id.init = tree_id,
         spp.init = spp, tree_status.init = tree_status, dbh_in.init = dbh_in, 
         height_ft.init = height_ft, cclass.init = cclass, crownrat.init = crownrat, 
         damage1.init = damage1, damage2.init = damage2, damage3.init = damage3,
         damage1_pnw.init = damage1_pnw, damage2_pnw.init = damage2_pnw, 
         damage3_pnw.init = damage3_pnw, death_cause.init = death_cause, 
         death_year.init = death_year, tpa_unadj.init = tpa_unadj, 
         prev_tre_cn.init = prev_tre_cn, prev_status.init = prev_status,
         prev_cond.init = prev_cond, prev_subc.init = prev_subc, 
         reconcile.init = reconcile, snag_missing.init = snag_missing)

names(combined_wide.present_at_remeasure)
names(combined_wide.not_present_at_remeasure)

names(combined_wide.not_present_at_remeasure)[!is.element(names(combined_wide.not_present_at_remeasure),
                                                          names(combined_wide.present_at_remeasure))]


combined_wide = 
  # bind rows correctly makes all the .re tree data in "not present at remeasure"  
  # NA; a subplot with no trees at the remeasure will show up in 
  # the "not_present_at_remeasure" table, and a subplot with no trees 
  # at the initial observation will show up in the "present_at_remeasure"
  # table; a subplot with trees at both times will likewise show up in the 
  # present at remeasure table, and a subplot with no trees at either time
  # will not show up
  bind_rows(combined_wide.present_at_remeasure,
            combined_wide.not_present_at_remeasure)

ggplot(data = 
         combined_wide,
       aes(x = inv_year_nominal.re,
           fill = tree_status.re))+
  geom_bar()

ggplot(data = 
         combined_wide %>% filter(is.na(tree_status.re)),
       aes(x = as.character(inv_year_nominal.init),
           fill = inv_plotstatus.re))+
  geom_bar()

ggplot(data = 
         combined_wide,
       aes(x = as.character(inv_year_nominal.init),
           fill = inv_plotstatus.re))+
  geom_bar()+
  facet_wrap(~inv_plotstatus.init, scales = 'free_y')


ggplot(data = 
         combined_wide,
       aes(x = as.character(inv_year_nominal.init),
           fill = inv_subpstatus.re))+
  geom_bar()+
  facet_wrap(~inv_subpstatus.init, scales = 'free_y')


ggplot(data = 
         combined_wide %>% filter(inv_subpstatus.re=='nonsampled'),
       aes(x = as.character(inv_year_nominal.init),
           fill = inv_nonsamp.re))+
  geom_bar()+
  facet_wrap(~inv_subpstatus.init, scales = 'free_y')

ggplot(data = 
         combined_wide,
       aes(x = as.character(inv_year_nominal.re),
           fill = inv_subpstatus.init))+
  geom_bar()+
  facet_wrap(~inv_subpstatus.re)

ggplot(data = 
         combined_wide %>% filter(is.na(tree_status.re)),
       aes(x = as.character(inv_year_nominal.init),
           fill = inv_nonsamp.re))+
  geom_bar()

# it is suspicious that the number of disappearing trees varies so much 
# from year to year
ggplot(data = 
         combined_wide %>% filter(is.na(tree_status.re)),
       aes(x = as.character(inv_year_nominal.init),
           fill = tree_status.init))+
  geom_bar()
# but they are almost all trees which were dead, harvested, or out of sample 
# in the initial measurement, or from 2001 plots which were skipped and 
# never got a followup

# of the trees which were alive at the initial measurement but went missing,
# the vast majority were on plots which were sampled in 2001 (or 2002) and then skpped
# during the remeasurement (~1700 trees); those initially sampled in
# 2004-2009 were all on plots which were 
# resampled but had no forest condition present at either initial or remeasure 
# (<50 trees per year)
ggplot(data = 
         combined_wide %>% 
         filter(is.na(tree_status.re)&
                  tree_status.init=='live'),
       aes(x = as.character(inv_year_nominal.init),
           fill = inv_nonsamp.re))+
  geom_bar()
ggplot(data = 
         combined_wide %>%
         filter(is.na(tree_status.re)&
                  tree_status.init=='live'),
       aes(x = as.character(inv_year_nominal.init),
           fill = inv_plotstatus.re))+
  geom_bar()
ggplot(data = 
         combined_wide,
       aes(x = as.character(inv_year_nominal.init),
           fill = tree_status.init))+
  geom_bar()
#### final filtering ###########################################################

combined_wide.filtered = 
  combined_wide %>%
  # trees which were dead, cut, or out of sample at the initial measurement don't provide any useful 
  # info about recruitment, mortality, or growth rates, so ditch them
  filter(!is.element(tree_status.init,
                     c('dead', 'harvested', 'outofsample'))) %>%
  
  # we're also not interested in any trees on plots (or subplots) which were nonsampled at 
  # either the initial or remeasurement;
  # this might introduce small bias into the dataset (this is what the FIA programs 
  # complicated stratification scheme attempts to adjust for), because plots 
  # may be nonsampled due to hazard or denied access, both of which would plausibly
  # correlate with forest processes and structure (e.g. steep topography, 
  #private landownership). However, I don't have a good way to correct for the bias 
  # (the FIA stratification system is super opaque and I'm not sure how it 
  # interacts with the ways i've filtered and categorized the data). In any case, 
  # I expect that any bias would be very minor, for two reasons. First, 
  # trees which were on subplots or plots which were nonsampled at either 
  # initial or followup make up less than 4% of the dataset pre-filtering. 
  # Second, nonsampled status (e.g. due to steep terrain or private ownership) 
  # would have have a strong **interaction** with forest health variables 
  # and disturbances for dropping nonsampled plots to have a major effect on 
  # parameter estimates for the effect of forest health / disturbances on 
  # demographic rates. Such an interaction is plausible, but IMO would be 
  # weaker than the correlation between nonsampled status and the 
  # forest health / disturbance variables themselves, and probalby weak relative 
  # to the main effects of the forest health / disturbance variables. If nonsampled 
  # status has a strong correlation with the forest health / disturbance 
  # variables themselves (or other important variables I'm not capturing), it 
  # could conceiably bias our esimates of the overall demographic rates, which 
  # is a more realistic problem. However, with only 4% of the sample falling 
  # into the 'nonsampled' bin, the effect on demographic rates would have to 
  # be very strong to cause a major bias in "population-wide" rates estimated 
  # from only the sampled plots/subplots. 
  filter(inv_plotstatus.init!='nonsampled'&
           inv_plotstatus.re!='nonsampled'&
           inv_subpstatus.init!='nonsampled'&
           inv_subpstatus.re!='nonsampled') %>%
  
  # we want to exclude any trees which were out of sample at the followup 
  # measurement
  filter(tree_status.re != 'outofsample' | is.na(tree_status.re)) %>%
  
  # and we want to exclude any trees which should have been in the sample at 
  # the initial measurement (reconcile.re=='missed dead', 'missed_live', )
  filter(!is.element(reconcile.re, c('missed dead', 'missed live'))) %>%
  
  # and i'm not sure what else procedural change means when its coded 
  # as 'live' and not NA at the followup, but get rid of those
  filter(!is.element(reconcile.re, c('missing (moved)', 
                                     'previous error',
                                     'procedural change',
                                     'shrank',
                                     'area nonsampled')))


#### seedlings #################################################################

# when plot manual < 2.0, they didn't actually count more than 6 seedlings, 
# instead recording NA for treecount and 6 for treecount_calc. Rather than 
# try to incorporate the seedlings counts from 2001-2003, i'm just 
# going to ditch them
head(seedlings)
head(subplots)
head(treelist)
seedlings_combined = 
  seedlings %>%
  as_tibble() %>%
  
  # keep only observations on included subplots
  right_join(select(subplots, plt_cn, state_id, plot_id, subp_id)) %>%
  
  select(plt_cn, state_id, plot_id, subp_id, spp, treecount_calc) %>%
  
  # get complete combinations of plt_cn and spp, fill in the missings with 0
  complete(nesting(plt_cn, state_id, plot_id, subp_id), spp) %>%
  mutate(treecount_calc = ifelse(is.na(treecount_calc), 0, treecount_calc)) %>%
  
  right_join(subplots) %>%
  
  # ditch observations with manual < 2.0
  filter(inv_manual >= 2.0)


#### quality control ###########################################################


ggplot(data = combined_wide.filtered %>% filter(is.na(tree_status.init)),
       aes(x = as.character(inv_year_nominal.re),
           fill = reconcile.re))+
  geom_bar()+
  facet_grid(tree_status.re~., scales = 'free_y')

ggplot(data = combined_wide.filtered,
       aes(x = as.character(inv_year_nominal.re),
           fill = tree_status.re))+
  geom_bar()+
  facet_grid(reconcile.re~., scales = 'free_y')


names(combined_wide.filtered)
unique(combined_wide.filtered$tree_status.re)
ggplot(data = combined_wide.filtered,
       aes(x = as.character(inv_year_nominal.re),
           fill = tree_status.re))+
  geom_bar()

ggplot(data = combined_wide.filtered,
       aes(x = as.character(inv_year_nominal.re),
           fill = reconcile.re))+
  geom_bar()+
  facet_grid(tree_status.re~., scales = 'free_y')
ggplot(data = combined_wide.filtered,
       aes(x = as.character(inv_year_nominal.re),
           fill = reconcile.re))+
  geom_bar()+
  facet_grid(tree_status.init~., scales = 'free_y')
ggplot(data = combined_wide.filtered %>% filter(tree_status.re=='outofsample'),
       aes(x = as.character(inv_year_nominal.re),
           fill = tree_status.init))+
  geom_bar()

ggplot(data = combined_wide.filtered %>% filter(tree_status.re=='dead'),
       aes(x = as.character(inv_year_nominal.re),
           fill = snag_missing.re))


# does the breakpoint diameter always match up with the tpa right? yes
combined_wide.filtered %>%
  mutate(diam_class.init = 
           ifelse(dbh_in.init<24,
                  '<24',
                  ifelse(dbh_in.init<30,
                         '24-30',
                         '>30'))) %>%
  group_by(inv_macrobrk.init, diam_class.init) %>%
  summarise(unique(round(tpa_unadj.init, 2))) %>%
  ungroup()
combined_wide.filtered %>%
  mutate(diam_class.re = 
           ifelse(dbh_in.re<24,
                  '<24',
                  ifelse(dbh_in.re<30,
                         '24-30',
                         '>30'))) %>%
  group_by(inv_macrobrk.re, diam_class.re) %>%
  summarise(unique(round(tpa_unadj.re, 2))) %>%
  ungroup()


ggplot(data = combined_wide.filtered,
       aes(x = as.character(inv_year_nominal.init),
           fill = damage1.init))+
  geom_bar()

unique(combined_wide.filtered$damage1.init)

#### finesse post-join columns  ################################################

combined_wide.filtered = 
  combined_wide.filtered %>%
  dplyr::select(state_id, plot_id, subp_id, tree_id.init, tree_id.re,
         plt_cn.init = prev_plt_cn, plt_cn.re = plt_cn, 
         tre_cn.init = prev_tre_cn, tre_cn.re = tre_cn,
         inv_year_nominal.init, inv_year_nominal.re,
         inv_date.init, inv_date.re, 
         elev_ft, lat, lon, ecosub_cd,
         subp_fire.init = fire.init, subp_fire.re = fire.re, 
         subp_fireyr.init = fire_yr.init, subp_fireyr.re = fire_yr.re,
         subp_insects.init = insects.init, subp_insects.re = insects.re,
         subp_disease.init = disease.init, subp_disease.re = disease.re,
         subp_drought.init = drought.init, subp_drought.re = drought.re,
         subp_supp.init = suppression.init, supb_supp.re = suppression.re,
         subp_cutting.init = cutting.init, subp_cutting.re = cutting.re,
         subp_siteprep.init = siteprep.init, subp_siteprep.re = siteprep.re,
         subp_artificialregen.init = artificialregen.init, 
         subp_artificialregen.re = artificialregen.re,
         subp_naturalregen.init = naturalregen.init,
         subp_naturalregen.re = naturalregen.re,
         spp.init, spp.re, tree_status.init, tree_status.re,
         dbh_in.init, dbh_in.re, height_ft.init, height_ft.re,
         cclass.init, cclass.re, crownrat.init, crownrat.re,
         damage1.init, damage2.init, damage3.init, 
         damage1_pnw.init, damage2_pnw.init, damage3_pnw.init,
         damage1.re, damage2.re, damage3.re, 
         damage1_pnw.re, damage2_pnw.re, damage3_pnw.re,
         death_cause.re, death_year.re, tpa_unadj.init, tpa_unadj.re) %>%
  
  # assume NA means 'none'; this isn't really true because the damageX and 
  # damageX_pnw columns have different temporal ranges of coverage (2013-on and 
  # 2000-2013, respectively, with some temporal structure in whether crews tended to 
  # explicitly record 'none', or just leave it blank); but it makes the code 
  # below easier and doesn't affect our interpretation
  mutate(
    damage1.init = ifelse(is.na(damage1.init), 'none', damage1.init),
    damage2.init = ifelse(is.na(damage2.init), 'none', damage2.init),
    damage3.init = ifelse(is.na(damage3.init), 'none', damage3.init),
    damage1_pnw.init = ifelse(is.na(damage1_pnw.init), 'none', damage1_pnw.init),
    damage2_pnw.init = ifelse(is.na(damage2_pnw.init), 'none', damage2_pnw.init),
    damage3_pnw.init = ifelse(is.na(damage3_pnw.init), 'none', damage3_pnw.init),
    damage1.re = ifelse(is.na(damage1.re), 'none', damage1.re),
    damage2.re = ifelse(is.na(damage2.re), 'none', damage2.re),
    damage3.re = ifelse(is.na(damage3.re), 'none', damage3.re),
    damage1_pnw.re = ifelse(is.na(damage1_pnw.re), 'none', damage1_pnw.re),
    damage2_pnw.re = ifelse(is.na(damage2_pnw.re), 'none', damage2_pnw.re),
    damage3_pnw.re = ifelse(is.na(damage3_pnw.re), 'none', damage3_pnw.re)) %>%
  
  # turn these damage columns into binary flags
  mutate(tree_fire.init = 
           damage1.init=='fire'|damage2.init=='fire'|damage3.init=='fire'|
           damage1_pnw.init=='fire'|damage2_pnw.init=='fire'|damage3_pnw.init=='fire',
         tree_disease.init = 
           damage1.init=='disease'|damage2.init=='disease'|damage3.init=='disease'|
           damage1_pnw.init=='disease'|damage2_pnw.init=='disease'|damage3_pnw.init=='disease',
         tree_insects.init = 
           damage1.init=='insects'|damage2.init=='insects'|damage3.init=='insects'|
           damage1_pnw.init=='insects'|damage2_pnw.init=='insects'|damage3_pnw.init=='insects',
         tree_drought.init = 
           damage1.init=='drought'|damage2.init=='drought'|damage3.init=='drought'|
           damage1_pnw.init=='drought'|damage2_pnw.init=='drought'|damage3_pnw.init=='drought',
         tree_suppression.init = 
           damage1.init=='suppression'|damage2.init=='suppression'|damage3.init=='suppression'|
           damage1_pnw.init=='suppression'|damage2_pnw.init=='suppression'|damage3_pnw.init=='suppression',
         tree_cutting.init = 
           damage1.init=='harvest'|damage2.init=='harvest'|damage3.init=='harvest'|
           damage1_pnw.init=='harvest'|damage2_pnw.init=='harvest'|damage3_pnw.init=='harvest',
         tree_other.init = 
           is.element(damage1.init, c('abiotic', 'animal', 'other')) |
           is.element(damage2.init, c('abiotic', 'animal', 'other')) |
           is.element(damage3.init, c('abiotic', 'animal', 'other')) |
           is.element(damage1_pnw.init, c('abiotic', 'animal', 'other'))|
           is.element(damage2_pnw.init, c('abiotic', 'animal', 'other'))|
           is.element(damage3_pnw.init, c('abiotic', 'animal', 'other')),
         
         tree_fire.re = 
           damage1.re=='fire'|damage2.re=='fire'|damage3.re=='fire'|
           damage1_pnw.re=='fire'|damage2_pnw.re=='fire'|damage3_pnw.re=='fire',
         tree_disease.re = 
           damage1.re=='disease'|damage2.re=='disease'|damage3.re=='disease'|
           damage1_pnw.re=='disease'|damage2_pnw.re=='disease'|damage3_pnw.re=='disease',
         tree_insects.re = 
           damage1.re=='insects'|damage2.re=='insects'|damage3.re=='insects'|
           damage1_pnw.re=='insects'|damage2_pnw.re=='insects'|damage3_pnw.re=='insects',
         tree_drought.re = 
           damage1.re=='drought'|damage2.re=='drought'|damage3.re=='drought'|
           damage1_pnw.re=='drought'|damage2_pnw.re=='drought'|damage3_pnw.re=='drought',
         tree_suppression.re = 
           damage1.re=='suppression'|damage2.re=='suppression'|damage3.re=='suppression'|
           damage1_pnw.re=='suppression'|damage2_pnw.re=='suppression'|damage3_pnw.re=='suppression',
         tree_cutting.re = 
           damage1.re=='harvest'|damage2.re=='harvest'|damage3.re=='harvest'|
           damage1_pnw.re=='harvest'|damage2_pnw.re=='harvest'|damage3_pnw.re=='harvest',
         tree_other.re = 
           is.element(damage1.re, c('abiotic', 'animal', 'other')) |
           is.element(damage2.re, c('abiotic', 'animal', 'other')) |
           is.element(damage3.re, c('abiotic', 'animal', 'other')) |
           is.element(damage1_pnw.re, c('abiotic', 'animal', 'other'))|
           is.element(damage2_pnw.re, c('abiotic', 'animal', 'other'))|
           is.element(damage3_pnw.re, c('abiotic', 'animal', 'other'))) %>%
    
    # keep the binary flags, get rid of the damage columns
    dplyr::select(-damage1.init, -damage2.init, -damage3.init, 
           -damage1_pnw.init, -damage2_pnw.init, -damage3_pnw.init,
           -damage1.re, -damage2.re, -damage3.re, 
           -damage1_pnw.re, -damage2_pnw.re, -damage3_pnw.re) %>%
  
  # treat state code as character
  mutate(state_id = as.character(state_id))

summary(combined_wide.filtered)

#### derived columns ###########################################################

# can we just use crown class as a proxy for competition? looks like yes
ggplot(data = combined_wide.filtered,
       aes(x = inv_year_nominal.init, fill = cclass.init))+
  geom_bar()+
  facet_wrap(~tree_status.init)

ggplot(data = combined_wide.filtered,
       aes(x = inv_year_nominal.re, fill = cclass.re))+
  geom_bar()+
  facet_wrap(~tree_status.re)


#### quality control ###########################################################

# checking for duplicates
nrow(combined_wide.filtered)
length(unique(combined_wide.filtered$tre_cn.init))
length(unique(combined_wide.filtered$tre_cn.re))

length(unique(treelist$tre_cn))
head(treelist_wide)
length(unique(treelist_wide$tre_cn))
nrow(treelist_wide %>% filter(!is.na(tre_cn)))
#### write results #############################################################
write.csv(combined_wide.filtered,
          here::here('02-data',
                     '02-for_analysis',
                     'FIA_censuses.csv'),
          row.names = FALSE)
saveRDS(combined_wide.filtered,
        here::here('02-data',
                   '02-for_analysis',
                   'FIA_censuses.rds'))

write.csv(seedlings_combined,
          here::here('02-data',
                     '02-for_analysis',
                     'FIA_seedlings.csv'),
          row.names = FALSE)
saveRDS(seedlings_combined,
        here::here('02-data',
                   '02-for_analysis',
                   'FIA_seedlings.rds'))
