library(here)
library(sf)
library(tidyverse)
library(RODBC)
library(terra)

db = 
  odbcConnectAccess2007(here::here('02-data', 
                               '00-source',
                               'usgs',
                               'DemogDATABASE2021_10_08.accdb'))



usgs = 
  sqlFetch(db, 'allplotarchive') %>%
  # CIRQUE is too much of a headache to be worth a plot with no PILA, PINELOCO 
  # and PINEBUVI too new to have a growth model, EMRIDGE has a bunch of non-ingrowth
  # new trees added in
  filter(PLOT != 'CIRQUE' & PLOT != 'PINELOCO' & PLOT != 'EMRIDGE' & PLOT != 'PINEBUVI')
  

#### individual data ###########################################################

# tree_id | spp | plot | year_0 | dbh_0 | status_0 | year_1 | dbh_1 | status_1 | 
#  fire | wpbr | cwd | BA
individual_data = 
  usgs %>%
  mutate(tree_id = paste(PLOT, SUBPLOT, TAGNUMBER, sep = '-')) %>%
  
  # there are a number of instances where a census was split across multiple years
  # in every case, the majority of the inventory took place in the earliest 
  # year, which is also the year that is on 5-year census cycle. Want to set the min 
  # year as the canonical year of each census. Fudging the years in this way 
  # allows us to assume an even 5 year census interval across the growth, 
  # mortality, and recruitment submodels, which is important for being able to 
  # integrate them into an IPM later. The way we're modelling growth (really, 
  # "size at second census" with a term for initial size) doesn't allow for 
  # flexibly including a term for teh census duration, as far as I can see.
  # specifically: 
    ## GIBBS 1996/1997, LOLOG 1984/1985, SUPILA 1982/1983, SURIP 1982/1983
    ## YOHOPIPO 1996/1997, BBBPIPO 2013 and 2014 -> 2012, CCRPIPO 2016 -> 2015
    ## EMSLOPE 2009->2008, EMSLOPE 2014->2013, FFS2BURN 2012->2011
    ## FFS6BURN 2013->2012, FFS7CONTROL 2012->2011, GIBBS 2014->2013 
    ## HOMEPICO 2011->2010, LMCC 2014-> 2013, LOGPIJE 2016->2015, LOLOG 2014->2012
    ## LOTHOR 2014->2013, SURIP 2015->2014, WTABMA 2014->2013
  # first remove the offending year columns
  select(-YEAR1, -YEAR2, -YEAR3, -YEAR4, -YEAR5, -YEAR6, -YEAR7, -YEAR8) %>%
  # join in a fixed version of the years
  left_join(
    usgs %>%
      # fix janky years in FFS6BURN
      mutate(YEAR5 = ifelse(PLOT=='FFS6BURN'&YEAR5==2013,2021,YEAR5)) %>%
      group_by(PLOT) %>%
      summarise(YEAR1 = min(YEAR1, na.rm = TRUE),
                YEAR2 = min(YEAR2, na.rm = TRUE),
                YEAR3 = min(YEAR3, na.rm = TRUE),
                YEAR4 = min(YEAR4, na.rm = TRUE),
                YEAR5 = min(YEAR5, na.rm = TRUE),
                YEAR6 = min(YEAR6, na.rm = TRUE),
                YEAR7 = min(YEAR7, na.rm = TRUE),
                YEAR8 = min(YEAR8, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(
        YEAR1 = ifelse(is.infinite(YEAR1), NA, YEAR1),
        YEAR2 = ifelse(is.infinite(YEAR2), NA, YEAR2),
        YEAR3 = ifelse(is.infinite(YEAR3), NA, YEAR3),
        YEAR4 = ifelse(is.infinite(YEAR4), NA, YEAR4),
        YEAR5 = ifelse(is.infinite(YEAR5), NA, YEAR5),
        YEAR6 = ifelse(is.infinite(YEAR6), NA, YEAR6),
        YEAR7 = ifelse(is.infinite(YEAR7), NA, YEAR7),
        YEAR8 = ifelse(is.infinite(YEAR8), NA, YEAR8)
      )
  ) %>%
  
  select(plot = PLOT, tree_id, spp = SppCode, yearin = IngrowthYear, 
         yearfirst = YearFirstRecorded, 
         yeardead = MortalityYear, dbh_recr = DBH0, dbh_mort = DBHMort,
         dbh_1 = DBH1, dbh_2 = DBH2, dbh_3 = DBH3, dbh_4 = DBH4, dbh_5 = DBH5,
         dbh_6 = DBH6, dbh_7 = DBH7, dbh_8 = DBH8,
         year_1 = YEAR1, year_2 = YEAR2, year_3 = YEAR3, year_4 = YEAR4, 
         year_5 = YEAR5, year_6 = YEAR6, year_7 = YEAR7, year_8 = YEAR8) %>%
  
  # pivot longer, then wider, to get a row per observation, 
  # with DBH and year at time of observation
  pivot_longer(cols = c(dbh_1, dbh_2, dbh_3, dbh_4, dbh_5, dbh_6, dbh_7, dbh_8,
                        year_1, year_2, year_3, year_4, year_5, year_6, year_7, 
                        year_8),
               names_to = c('variable', 'obs'), names_sep = '_') %>%
  # there's some dubplicate rows, so aggregate these 
  group_by(plot, tree_id, spp, yearin, yearfirst, yeardead, dbh_recr, dbh_mort, 
           variable, obs) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup()  %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  
  # ditch rows with no DBH or year
  filter(!(is.nan(dbh)&is.nan(year))) %>%
  
  # replace -99 with NA
  mutate(yearin = ifelse(yearin<1980, NA, yearin),
         yearfirst = ifelse(yearfirst < 1980, NA, yearfirst),
         yeardead = ifelse(yeardead<1980, NA, yeardead),
         dbh_recr = ifelse(dbh_recr<0, NA, dbh_recr),
         dbh_mort = ifelse(dbh_mort<0, NA, dbh_mort),
         dbh = ifelse(dbh<0, NA, dbh),
         year = ifelse(year<0, NA, year),
         obs = as.integer(obs)) %>%
  
  # ditch rows with unrealistic (typo) DBHs
  filter(is.na(dbh) | (dbh > 0 & !(dbh>300 & spp != 'SEGI'))) %>%
  
  # add a status flag
  mutate(status = ifelse(is.na(yeardead), 'live',
                         ifelse(year<yeardead, 
                                'live',
                                'dead'))) %>%
  
  # ditch live trees with missing DBHs
  filter(!(is.na(dbh)&status=='live')) %>%
  
  # get the dbh and status at the end of each census interval
  select(plot, tree_id, spp, obs_0 = obs, year, status, dbh) %>%
  mutate(obs_1 = obs_0 + 1) %>% 
  left_join(select(., plot, tree_id, spp, obs_0, year, status, dbh),
            by = c('plot' = 'plot',
                   'tree_id' = 'tree_id',
                   'spp' = 'spp', 
                   'obs_1' = 'obs_0'),
            suffix = c('_0', '_T')) %>%
  
  # get rid of trees which are dead at the start of a census, or which haven't 
  # had a followup census
  filter(status_0=='live' & !is.na(year_T)) 

#### size metadata #############################################################

size_metadata = 
  individual_data %>%
  mutate(dbh_class = cut(dbh_0,
                         breaks = c(seq(from = 0, to = 200, by = 10),1000),
                         labels = FALSE)) %>%
  group_by(dbh_class) %>%
  summarise(dbh_med = median(dbh_0)) %>%
  ungroup() %>%
  mutate(lower = seq(from = 0, to = 200, by = 10),
         upper = c(seq(from = 10, to = 200, by = 10), 1000))

#### complete crossing of plot and year ########################################


census_years = 
  usgs %>%
  
  # fix janky years in FFS6BURN
  mutate(YEAR5 = ifelse(PLOT=='FFS6BURN'&YEAR5==2013,2021,YEAR5)) %>%

  group_by(PLOT) %>%
  summarise(YEAR1 = min(YEAR1, na.rm = TRUE),
            YEAR2 = min(YEAR2, na.rm = TRUE),
            YEAR3 = min(YEAR3, na.rm = TRUE),
            YEAR4 = min(YEAR4, na.rm = TRUE),
            YEAR5 = min(YEAR5, na.rm = TRUE),
            YEAR6 = min(YEAR6, na.rm = TRUE),
            YEAR7 = min(YEAR7, na.rm = TRUE),
            YEAR8 = min(YEAR8, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    YEAR1 = ifelse(is.infinite(YEAR1), NA, YEAR1),
    YEAR2 = ifelse(is.infinite(YEAR2), NA, YEAR2),
    YEAR3 = ifelse(is.infinite(YEAR3), NA, YEAR3),
    YEAR4 = ifelse(is.infinite(YEAR4), NA, YEAR4),
    YEAR5 = ifelse(is.infinite(YEAR5), NA, YEAR5),
    YEAR6 = ifelse(is.infinite(YEAR6), NA, YEAR6),
    YEAR7 = ifelse(is.infinite(YEAR7), NA, YEAR7),
    YEAR8 = ifelse(is.infinite(YEAR8), NA, YEAR8)
  ) %>%
  pivot_longer(cols = c('YEAR1', 'YEAR2', 'YEAR3', 'YEAR4', 'YEAR5', 'YEAR6',
                        'YEAR7', 'YEAR8'),
               names_to = 'census', values_to = 'year') %>%
  mutate(census = as.integer(gsub(census, pattern = 'YEAR', replacement = '')),
         census_1 = census+1) %>%
  select(plot = PLOT, census, census_1, year) %>%
  left_join(select(., plot, census, year),
            by = c('plot' = 'plot',
                   'census_1' = 'census'),
            suffix = c('_0', '_T')) %>%
  filter(!is.na(year_0)&!is.na(year_T)) %>%
  select(-census_1)


#### adding plot level spatial data ############################################

plots_spatial = 
  st_read(here::here('02-data',
                     '00-source',
                     'usgs',
                     'Boundary Polygons.gdb'),
          layer = 'allplotscombined') %>%
  st_transform(crs = 'EPSG:4326') %>%
  st_centroid() %>%
  select(plot = PLOT, Shape, areaha) %>%
  mutate(lat = st_coordinates(.)[,'Y'],
         lon = st_coordinates(.)[,'X']) %>%
  as.data.frame() %>%
  select(-Shape) %>%
  as_tibble()

plots_spatial = 
  st_read(here::here('02-data',
                     '00-source',
                     'usgs',
                     'Boundary Polygons.gdb'),
          layer = 'allplotscombined') %>%
  st_transform(crs = 'EPSG:4326') %>%
  st_centroid() %>%
  select(plot = PLOT, Shape, areaha) %>%
  mutate(lat = st_coordinates(.)[,'Y'],
         lon = st_coordinates(.)[,'X']) %>%
  as.data.frame() %>%
  select(-Shape) %>%
  as_tibble()





#### plot:census-level explanatory variables ###################################

head(census_years)


census_data = 
  census_years %>%
  
  # join in the spatial information
  left_join(plots_spatial) %>%
  
  # fire: did a fire occur on the plot during the census interval?
  left_join(# per adrian not flagging LMCC as burned because only a small portion of the 
            # plot was affected
            data.frame(plot = c('CRCRPIPO', 'FFS2BURN', 'FFS5BURN', 'FFS6BURN',
                                'LOTHAR', 'UPTHAR', 'YOHOPIPO'),
                       fire_year1 = c(2009, 2001, 2001, 2001, 1990, 1990, 2007),
                       fire_year2 = c(NA, NA, NA, NA, 2004, 2004, 2011)),
            by = c('plot' = 'plot')) %>%
  mutate(burned_0T = 
           ifelse((!is.na(fire_year1)&
                     fire_year1>=year_0 &
                     fire_year1<year_T)|
                    (!is.na(fire_year2)&
                       fire_year2>=year_0 &
                       fire_year2<year_T),
                  TRUE,
                  FALSE)) %>%
  select(-fire_year1, -fire_year2) %>%
  
  # wpbr: has a tree ever died of WPBR on the plot?
  mutate(wpbr = 
           ifelse(is.element(plot,
                             usgs %>%
                                group_by(PLOT) %>%
                                summarise(wpbr = 
                                            any(Proximate==61, na.rm = TRUE)|
                                            any(Primary==61, na.rm = TRUE)) %>%
                                ungroup() %>%
                                filter(wpbr) %>%
                                pull(PLOT)),
                  TRUE,
                  FALSE)) %>%
  # basal area: plot-level live BA at the beginning of the census
  left_join(individual_data %>%
              filter(status_0 == 'live') %>%
              left_join(plots_spatial %>%
                          select(plot, areaha),
                        by = c('plot' = 'plot')) %>%
              group_by(plot, year_0) %>%
              summarise(ba_m2ha = 
                          sum(pi*((dbh_0/200)**2) / areaha)) %>%
              ungroup(),
            by = c('plot' = 'plot',
                   'year_0' = 'year_0'))
  
# drought: 90th percentile (scaled) CWD departure over the census interval
plots_bbox = 
  list('x_min' = min(plots_spatial$lon)-1,
       'x_max' = max(plots_spatial$lon)+1,
       'y_min' = min(plots_spatial$lat)-1,
       'y_max' = max(plots_spatial$lat)+1)

cwd_growseason_means = 
  lapply(X = 1982:2020,
         FUN = function(y){
           
           cwd_year = 
             rast(here::here('02-data',
                             '00-source',
                             'terraclimate',
                             paste0('TerraClimate_def_',y,'.nc')))
           
           cwd_year = 
             crop(cwd_year,
                  c(plots_bbox$x_min, plots_bbox$x_max,
                    plots_bbox$y_min, plots_bbox$y_max))
           
           cwd_year = 
             mean(cwd_year[[5:10]])
           
           return(cwd_year)
            
         }) %>%
  rast()

names(cwd_growseason_means) = 
  paste0('growseasonmean_',as.character(1982:2020))

plot(cwd_growseason_means)

# get departure from "normal" (1982-2020 mean) CWD for each 
# year on each location
cwd_departure = 
  cwd_growseason_means - mean(cwd_growseason_means)

names(cwd_departure) = paste0('cwddeparture_', as.character(1982:2020))

plot(cwd_departure)

# extract the departure for each year at each plot location
cwd_departures = 
  extract(cwd_departure, census_data[,c('lon', 'lat')]) %>%
  bind_cols('plot' = census_data$plot,
            'year_0' = census_data$year_0,
            'year_T' = census_data$year_T)

cwd_departure_span = 
  sapply(X = 1:nrow(cwd_departures),
         FUN = function(i){
           year_range = seq(from = cwd_departures$year_0[i],
                            to = cwd_departures$year_T[i]-1,
                            by = 1)
           
           columns_to_select = 
             paste0('cwddeparture_', year_range)
           
           values_in_span = 
             as.numeric(cwd_departures[i,columns_to_select])
           
           # return the 90th percentile CWD departure
           return(as.numeric(quantile(values_in_span, probs = 0.9)))
           
         })


census_data$cwd_departure90 = cwd_departure_span


head(census_data)

census_data = 
  census_data %>%
  mutate(ba_scaled = as.numeric(scale(ba_m2ha)),
         cwd_scaled = as.numeric(scale(cwd_departure90)))

#### complete crossing of plot, year, species and size for fecundity ###########

# want one row per species:plot:year:dbh_class, with columns for explanatory variables and 
# the response, a count of the number of new recruits in that dbh class
# note that we only have annual recruitment inventories (which are recorded 
# as ingrowthyear in the trees table) starting in 1999
plot_year_spp_sizeclass = 
  
  # get all the combinations of plot, year, species, and size class
  census_data %>%
  expand(nesting(plot, census, year_0, year_T, areaha, lat, lon, 
                 burned_0T, wpbr, ba_m2ha, cwd_departure90),
         spp = unique(individual_data$spp)) %>%
  expand(nesting(plot, census, year_0, year_T, areaha, lat, lon, burned_0T,
                 wpbr, ba_m2ha, cwd_departure90, spp),
         dbh_class = 1:21) %>%
  
  # add in the median dbh for each size class
  left_join(size_metadata %>%
              select(dbh_class, dbh_med)) %>%
  
  arrange(plot, year_0, spp, dbh_class)

head(plot_year_spp_sizeclass)

#### add size class tallies to fecundity data ##################################
fecundity_data = 
  plot_year_spp_sizeclass %>%
  left_join(individual_data %>%
              mutate(dbh_class_0 = 
                       cut(dbh_0, 
                           breaks = c(seq(from = 0, to = 200, by = 20), 1000),
                           labels = FALSE)) %>%
              left_join(plots_spatial %>% 
                          select(plot, areaha)) %>%
              mutate(tph_0 = 1 / areaha) %>%
              group_by(plot, spp, year_0, dbh_class_0) %>%
              summarise(tph_0 = sum(tph_0)) %>%
              ungroup(),
            by = c('plot' = 'plot',
                   'spp' = 'spp',
                   'year_0' = 'year_0',
                   'dbh_class' = 'dbh_class_0')) %>%
  mutate(tph_0 = ifelse(is.na(tph_0), 0, tph_0)) %>%
  arrange(plot, year_0, spp, dbh_class)

#### new recruit sizes #########################################################

head(individual_data)

individual_recruits = 
  
  # get the first occurence of each tree; I'm ignoring the "ingrowth year" 
  # information (the annual recruitment surveys) because they only started in 
  # 1999, and I'm basing everything off the 5 year censuses anyways. This does 
  # include a bunch of cases where non-ingrowth new trees were added 
  # (specifically:   EMRIDGE 1999-2004 (many trees), EMSLOPE 1998 (24 trees), 
  # GIBBS 2008 subplot 3, and 27 (28, and 17 trees), Gibbs 2009 33 trees in 
  # subplot 23);  given the small number of trees involved 
  # for plots other than EMRIDGE, I'm going to filter these by size to get only 
  # the plausible (initial DBH < 10cm) new recruits
  individual_data %>%
  
  # first, ditch trees which were present at the initial census of each plot
  filter(!is.element(tree_id,
                     individual_data %>%
                     left_join(census_years %>%
                                 group_by(plot) %>%
                                 summarise(initial_year = min(year_0)) %>%
                                 ungroup()) %>%
                       filter(year_0 == initial_year) %>%
                       pull(tree_id))) %>%
  
  group_by(plot, tree_id, spp) %>%
  summarise(year_recr = min(year_0),
            dbh_recr = min(dbh_0)) %>%
  ungroup() %>%
  
  # exclude "recruits" with > 10 cm DBH at time of "recruitment"
  filter(dbh_recr < 10)


#### count of new recruits per census ##########################################

recruit_counts = 
  census_years %>%
  expand(nesting(plot, census, year_0, year_T),
         spp = unique(individual_data$spp)) %>%
  left_join(
    individual_recruits %>%
      group_by(plot, spp, year_recr) %>%
      summarise(count = n()) %>%
      ungroup(),
    by = c('plot' = 'plot',
           'year_T' = 'year_recr',
           'spp' = 'spp')
  ) %>%
  mutate(count = ifelse(is.na(count),0,count)) %>%
  arrange(plot, year_0, spp)
  

head(recruit_counts)


#### filter to sugar pine ######################################################

individual_data.pila = 
  individual_data %>%
  filter(spp == 'PILA')

recruit_counts.pila = 
  recruit_counts %>%
  mutate(pc = paste0(plot,'-',census)) %>%
  filter(spp == 'PILA' &
           is.element(pc,
                      individual_data.pila %>%
                        mutate(pc = paste0(plot, '-', obs_0)) %>%
                        pull(pc)))

individual_recruits.pila = 
  individual_recruits %>%
  filter(spp == 'PILA')

fecundity_data.pila = 
  fecundity_data %>%
  mutate(pc = paste0(plot,'-',census)) %>%
  filter(spp == 'PILA' &
           is.element(pc,
                      individual_data.pila %>%
                        mutate(pc = paste0(plot, '-', obs_0)) %>%
                        pull(pc)))


#### add indices ###############################################################

unique_plots.pila = 
  individual_data.pila %>%
  group_by(plot) %>%
  summarise() %>%
  mutate(plot_id.i = as.integer(factor(plot)))

unique_trees.pila = 
  individual_data.pila %>%
  group_by(tree_id) %>%
  summarise() %>%
  mutate(tree_id.i = as.integer(factor(tree_id)))

plots_spatial = 
  plots_spatial %>%
  left_join(unique_plots.pila)

individual_data.pila = 
  individual_data.pila %>%
  left_join(unique_plots.pila) %>%
  left_join(unique_trees.pila)

recruit_counts.pila = 
  recruit_counts.pila %>%
  left_join(unique_plots.pila)%>% 
  arrange(plot, census)

individual_recruits.pila = 
  individual_recruits.pila %>%
  left_join(unique_plots.pila)

fecundity_data.pila = 
  fecundity_data.pila %>%
  left_join(unique_plots.pila)  %>%
  arrange(plot, census, dbh_class)

#### survival data #############################################################

survival_data.pila = 
  individual_data.pila %>%
  left_join(census_data,
            by = c('plot' = 'plot',
                   'year_0' = 'year_0',
                   'year_T' = 'year_T',
                   'obs_0' = 'census')) %>%
  filter(burned_0T == FALSE) %>%
  mutate(intercept = 1,
         year_scaled = year_0 - 1982,
         dbh_0.m = dbh_0 * 0.01,
         dbh_fire = dbh_0.m * burned_0T,
         dbh_wpbr = dbh_0.m * wpbr,
         dbh_ba = dbh_0.m * ba_scaled,
         dbh_cwd = dbh_0.m * cwd_scaled,
         dbh_year = dbh_0.m * year_scaled) %>%
  select(plot, plot_id.i, tree_id, spp, obs_0, year_0, status_0, 
         year_T, status_T, intercept, dbh_0.m, year_scaled, fire = burned_0T, wpbr,
         ba_scaled, cwd_scaled, dbh_fire, dbh_wpbr, dbh_ba, dbh_cwd, dbh_year)



#### growth data ###############################################################

growth_data.pila = 
  individual_data.pila %>%
  filter(status_0 == 'live' & status_T == 'live') %>%
  left_join(census_data,
            by = c('plot' = 'plot',
                   'year_0' = 'year_0',
                   'year_T' = 'year_T',
                   'obs_0' = 'census')) %>%
  filter(burned_0T == FALSE) %>%
  mutate(intercept = 1,
         year_scaled = year_0-1982,
         dbh_0.m = dbh_0 * 0.01,
         dbh_T.m = dbh_T * 0.01,
         dbh_fire = dbh_0.m * burned_0T,
         dbh_wpbr = dbh_0.m * wpbr,
         dbh_ba = dbh_0.m * ba_scaled,
         dbh_cwd = dbh_0.m * cwd_scaled,
         dbh_year = dbh_0.m * year_scaled) %>%
  select(plot, plot_id.i, tree_id, spp, obs_0, year_0, status_0, 
         year_T, status_T, dbh_T.m, intercept, dbh_0.m, year_scaled, fire = burned_0T, wpbr,
         ba_scaled, cwd_scaled, dbh_fire, dbh_wpbr, dbh_ba, dbh_cwd, dbh_year)

head(growth_data.pila)

ggplot(data = survival_data.pila,
       aes(x = fire))+
  geom_bar()

#### fecundity data ############################################################

head(fecundity_data.pila)
fecundity_data.pila = 
  fecundity_data.pila %>%
  filter(burned_0T == FALSE) %>%
  mutate(intercept = 1,
         dbh_0.m = dbh_med*0.01,
         year_scaled = year_0 - 1982,
         dbh_year = dbh_0.m * year_scaled)


head(recruit_counts)

recruit_counts.pila = 
  recruit_counts.pila %>%
  left_join(census_data,
            by = c('plot' = 'plot',
                   'year_0' = 'year_0',
                   'year_T' = 'year_T',
                   'census' = 'census')) %>%
  filter(burned_0T == FALSE)

head(individual_recruits.pila)


individual_recruits.pila = 
  individual_recruits.pila %>%
  left_join(census_data,
            by = c('plot' = 'plot',
                   'year_recr' = 'year_0')) %>%
  filter(burned_0T == FALSE)

#### build model data ##########################################################

pila_data = 
  list(
    K = 2,
    P = nrow(unique_plots.pila),
    
    # survival data
    N_s = nrow(survival_data.pila),
    surv_s = as.integer(survival_data.pila$status_T=='live'),
    X_s = survival_data.pila[,c('intercept', 'dbh_0.m')] %>%
      as.matrix(),
    plotid_s = survival_data.pila$plot_id.i,
    
    # growth data
    N_g = nrow(growth_data.pila),
    dbhT_g = growth_data.pila$dbh_T.m,
    X_g = 
      growth_data.pila[,c('intercept', 'dbh_0.m')] %>% 
      as.matrix(),
    plotid_g = growth_data.pila$plot_id.i,
    
    # fecundity data
    N_f = nrow(fecundity_data.pila),
    C_f = nrow(fecundity_data.pila %>%
                 group_by(plot, census) %>%
                 summarise() %>%
                 ungroup()),
    M_f = nrow(size_metadata),
    cprime_f = recruit_counts.pila$count,
    n = matrix(nrow = nrow(size_metadata),
               ncol = nrow(fecundity_data.pila %>%
                             group_by(plot,census) %>%
                             summarise() %>%
                             ungroup()),
               data = fecundity_data.pila$tph_0,
               byrow = FALSE),
    X_f = 
      fecundity_data.pila[,c('intercept', 'dbh_0.m')] %>% 
      as.matrix(),
    plotid_f = fecundity_data.pila$plot_id.i,
    a = 
      fecundity_data.pila %>%
      group_by(plot, census, areaha) %>%
      summarise() %>%
      ungroup() %>%
      pull(areaha),
    
    # recruit size data
    N_r = nrow(individual_recruits.pila),
    logdbh_r = log(individual_recruits.pila$dbh_recr))


saveRDS(pila_data,
        here::here('02-data', '02-for_analysis', 'usgs', 'pila_data.rds'))

# this is useful for plotting later
saveRDS(plots_spatial,
        here::here('02-data', '02-for_analysis', 'usgs', 'plots_spatial.rds'))
