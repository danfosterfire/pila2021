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
  


#### core census data ##########################################################

# tree_id | spp | plot | year_0 | dbh_0 | status_0 | year_1 | dbh_1 | status_1 | 
#  fire | wpbr | cwd | BA
census_data = 
  usgs %>%
  mutate(tree_id = paste(PLOT, SUBPLOT, TAGNUMBER, sep = '-')) %>%
  select(plot = PLOT, tree_id, spp = SppCode, yearin = IngrowthYear, 
         yearfirst = YearFirstRecorded, 
         yeardead = MortalityYear, dbh_recr = DBH0, dbh_mort = DBHMort,
         dbh_1 = DBH1, dbh_2 = DBH2, dbh_3 = DBH3, dbh_4 = DBH4, dbh_5 = DBH5,
         dbh_6 = DBH6, dbh_7 = DBH7, dbh_8 = DBH8,
         year_1 = YEAR1, year_2 = YEAR2, year_3 = YEAR3, year_4 = YEAR4, 
         year_5 = YEAR5, year_6 = YEAR6, year_7 = YEAR7, year_8 = YEAR8) %>%
  
  # pivot longer, then wider, to get a row per observation, with DBH and year at time of observation
  pivot_longer(cols = c(dbh_1, dbh_2, dbh_3, dbh_4, dbh_5, dbh_6, dbh_7, dbh_8,
                        year_1, year_2, year_3, year_4, year_5, year_6, year_7, year_8),
               names_to = c('variable', 'obs'), names_sep = '_') %>%
  # there's some dubplicate rows, so aggregate these 
  group_by(plot, tree_id, spp, yearin, yearfirst, yeardead, dbh_recr, dbh_mort, variable, obs) %>%
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
  
  # add in rows for the recruitment observations, which occur in inyear, 
  # have DBH dbh_recr, and occur at observation min(obs)-1
  bind_rows(filter(., !is.na(dbh) & !is.na(year)) %>%
              group_by(plot, tree_id, spp, yearin, yearfirst, yeardead, dbh_recr, dbh_mort) %>%
              summarise(obs_min = min(obs), year_min = min(year)) %>%
              ungroup() %>%
              filter(!is.na(dbh_recr) & year_min !=yearin) %>%
              mutate(obs = obs_min-1,
                     dbh = dbh_recr,
                     year = yearin) %>%
              select(plot, tree_id, spp, yearin, yearfirst, yeardead, dbh_recr, dbh_mort,
                     obs, dbh, year)) %>%
  
  # add in rows for the final observation,l which occur in yeardead, have 
  # DBH dbh_mort, and occur at observation (max(obs)+1)
  bind_rows(filter(., !is.na(dbh) & !is.na(year)) %>%
              group_by(plot, tree_id, spp, yearin, yearfirst, yeardead, dbh_recr, dbh_mort) %>%
              summarise(obs_max = max(obs), year_max = max(year)) %>%
              ungroup() %>%
              filter(!is.na(dbh_mort) & year_max != yeardead) %>%
              mutate(obs = obs_max+1,
                     dbh = dbh_mort,
                     year = yeardead) %>%
              select(plot, tree_id, spp, yearin, yearfirst, yeardead, dbh_recr, dbh_mort,
                     obs, dbh, year)) %>%
  
  # ditch rows with unrealistic (typo) DBHs and missing DBHs
  filter(dbh > 0 & !(dbh>300 & spp != 'SEGI') & !is.na(dbh)) %>%
  mutate(status = ifelse(is.na(yeardead), 'live',
                         ifelse(year<yeardead, 
                                'live',
                                'dead'))) %>%
  
  # get the dbh and status at the end of each census interval
  select(plot, tree_id, spp, obs_0 = obs, year, status, dbh) %>%
  mutate(obs_1 = obs_0 + 1) %>% 
  left_join(select(., plot, tree_id, spp, obs_0, year, status, dbh),
            by = c('plot' = 'plot',
                   'tree_id' = 'tree_id',
                   'spp' = 'spp', 
                   'obs_1' = 'obs_0'),
            suffix = c('_0', '_1')) %>%
  
  # get rid of trees which are dead at the start of a census, or which haven't 
  # had a followup census
  filter(status_0=='live' & !is.na(year_1)) %>%
  
  # theres a handful of trees where it looks like the years were recorded in the 
  # wrong columns; assume everything else is OK and just swap the years when 
  # year_1 < year_0
  mutate(year_swap = year_0,
         swapped_years = ifelse(year_0 > year_1, TRUE, FALSE),
         year_0 = ifelse(swapped_years, year_1, year_0),
         year_1 = ifelse(swapped_years, year_swap, year_1)) %>%
  select(-year_swap, -swapped_years)


# have followup census for snags that died after the last true census, want to 
# get rid of them; warnings are for cases where there are no non-NA values
census_years = 
  usgs %>%
  group_by(PLOT) %>%
  summarise(YEAR1 = min(YEAR1, na.rm = TRUE),
            YEAR2 = max(YEAR2, na.rm = TRUE),
            YEAR3 = max(YEAR3, na.rm = TRUE),
            YEAR4 = max(YEAR4, na.rm = TRUE),
            YEAR5 = max(YEAR5, na.rm = TRUE),
            YEAR6 = max(YEAR6, na.rm = TRUE),
            YEAR7 = max(YEAR7, na.rm = TRUE),
            YEAR8 = max(YEAR8, na.rm = TRUE)) %>%
  ungroup()

census_years$latest_census = 
  sapply(X = 1:nrow(census_years),
         FUN = function(i){
           max(as.numeric(census_years[i,c('YEAR1', 'YEAR2', 'YEAR3', 'YEAR4',
                                           'YEAR5', 'YEAR6', 'YEAR7', 'YEAR8')]),
               na.rm = TRUE)
         })


census_data = 
  census_data %>%
  left_join(census_years %>%
              select(plot = PLOT, latest_census)) %>%
  filter(year_1 <= latest_census) %>%
  select(-latest_census)



#### fire ######################################################################



fire_years = 
  # per adrian not flagging LMCC as burned because only a small portion of the 
  # plot was affected
  data.frame(plot = c('CRCRPIPO', 'FFS2BURN', 'FFS5BURN', 'FFS6BURN',
                      'LOTHAR', 'UPTHAR', 'YOHOPIPO'),
             fire_year1 = c(2009, 2001, 2001, 2001, 1990, 1990, 2007),
             fire_year2 = c(NA, NA, NA, NA, 2004, 2004, 2011))

fire_years

head(census_data)

census_data = 
  census_data %>%
  left_join(fire_years) %>%
  mutate(burned_01 = 
           ifelse((!is.na(fire_year1)&
                     fire_year1>=year_0 &
                     fire_year1<year_1)|
                    (!is.na(fire_year2)&
                       fire_year2>=year_0 &
                       fire_year2<year_1),
                  TRUE,
                  FALSE)) %>%
  select(-fire_year1, -fire_year2)


#### plot:year level basal area ################################################

plot_metadata = 
  census_data %>%
  group_by(plot) %>%
  summarise(year_min = min(year_0), year_max = max(year_0)) %>%
  ungroup() %>%
  left_join(sf::st_read(here::here('02-data',
                            '00-source',
                            'usgs',
                            'Boundary Polygons.gdb'),
              layer = 'allplotscombined') %>%
              as.data.frame() %>%
              select(plot = PLOT, area_ha = areaha))

plot_ba = 
  do.call('bind_rows',
          lapply(X = 1:nrow(plot_metadata),
                 FUN = function(p){
                   year_min = plot_metadata$year_min[p]
                   year_max = plot_metadata$year_max[p]
                   area_ha = plot_metadata$area_ha[p]
                   
                   current_plot_ba = 
                     lapply(X = year_min:year_max,
                          FUN = function(y){
                            
                            current_year_ba = 
                              census_data %>%
                              filter(status_0=='live' & status_1=='live' & 
                                       year_0 <= y & year_1 >= y & 
                                       plot == plot_metadata$plot[p]) %>%
                              group_by(plot) %>%
                              summarise(ba_m2ha = 
                                          sum(pi*((dbh_0/200)**2) / plot_metadata$area_ha[p])) %>%
                              ungroup() %>%
                              mutate(year = y)
                            
                            
                            return(current_year_ba)
                          })
                   
                   return(current_plot_ba)
                 }))

head(census_data)

head(plot_ba)

census_data = 
  census_data %>%
  left_join(plot_ba %>%
            select(plot, bam2ha_0 = ba_m2ha, year_0 = year))

#### wpbr ######################################################################

head(usgs)
census_data = 
  census_data %>%
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
                  FALSE))




#### cwd #######################################################################


censuses_spatial = 
  st_read(here::here('02-data',
                     '00-source',
                     'usgs',
                     'Boundary Polygons.gdb'),
          layer = 'allplotscombined') %>%
  st_transform(crs = 'EPSG:4326') %>%
  st_centroid() %>%
  select(plot = PLOT, Shape) %>%
  right_join(census_data %>%
               group_by(plot, year_0, year_1) %>%
               summarise() %>%
               ungroup()) %>%
  mutate(lat = st_coordinates(.)[,'Y'],
         lon = st_coordinates(.)[,'X']) %>%
  as.data.frame() %>%
  select(-Shape) %>%
  as_tibble()

head(censuses_spatial)

plots_bbox = 
  list('x_min' = min(censuses_spatial$lon)-1,
       'x_max' = max(censuses_spatial$lon)+1,
       'y_min' = min(censuses_spatial$lat)-1,
       'y_max' = max(censuses_spatial$lat)+1)


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


cwd_departures = 
  extract(cwd_departure, censuses_spatial[,c('lon', 'lat')]) %>%
  bind_cols('plot' = censuses_spatial$plot,
            'year_0' = censuses_spatial$year_0,
            'year_1' = censuses_spatial$year_1)

cwd_departure_span = 
  sapply(X = 1:nrow(cwd_departures),
         FUN = function(i){
           
           year_range = seq(from = cwd_departures$year_0[i],
                            to = cwd_departures$year_1[i]-1,
                            by = 1)
           
           columns_to_select = 
             paste0('cwddeparture_', year_range)
           
           values_in_span = 
             as.numeric(cwd_departures[i,columns_to_select])
           
           # return the 90th percentile CWD departure
           return(as.numeric(quantile(values_in_span, probs = 0.9)))
           
         })

head(cwd_departure_span)

censuses_spatial$cwd_departure90 = cwd_departure_span

census_data = 
  census_data %>%
  left_join(censuses_spatial %>%
              select(plot, year_0, year_1, cwd_departure90))

head(census_data)

#### filter to just PILA #######################################################

census_data.pila = 
  census_data %>%
  filter(spp == 'PILA')

#### add indexing columns ######################################################

census_data.pila = 
  census_data.pila %>%
  mutate(plot.i = as.integer(factor(plot)),
         tree_id.i = as.integer(factor(tree_id)))


#### assemble survival data ####################################################

survival_data.pila = 
  census_data.pila %>%
  filter(status_0 == 'live') %>%
  mutate(intercept = 1)

#### assemble growth data ######################################################

growth_data.pila = 
  census_data.pila %>%
  filter(status_0 == 'live' & status_1 == 'live') %>%
  mutate(intercept = 1)


#### assemble fecundity data ###################################################


# want one row per species:plot:year:dbh_class, with columns for explanatory variables and 
# the response, a count of the number of new recruits in that dbh class
# note that we only have annual recruitment inventories (which are recorded 
# as ingrowthyear in the trees table) starting in 1999
plot_year_spp_sizeclass = 
  do.call('bind_rows',
          lapply(X = 1:nrow(plot_metadata),
                 FUN = function(p){
                   year_min = 
                     max(plot_metadata$year_min[p], 1999)
                   year_max = plot_metadata$year_max[p]
                   
                   lapply(X = year_min:year_max,
                          FUN = function(y){
                            
                            expand.grid(plot = plot_metadata$plot[p],
                                       year = y, 
                                       spp = unique(usgs$SppCode)) %>%
                              as_tibble() %>%
                              mutate(plot = as.character(plot),
                                     spp = as.character(spp))
                          })
                 }))  %>%
  expand(nesting(plot, year, spp),
         dbh_class = 1:21)

new_recruits = 
  usgs %>%
  mutate(tree_id = paste(PLOT, SUBPLOT, TAGNUMBER, sep = '-'))  %>%
  filter(!is.na(DBH0) & DBH0>=0 & !is.na(IngrowthYear) & IngrowthYear>=1982) %>%
  filter(DBH0 <= 7.5) %>% # 7.4 is the 99.9th percentile of DBH0; assume values above this are jank 
  select(plot = PLOT, tree_id, spp = SppCode, yearin = IngrowthYear,
         dbh_recr = DBH0) %>%
  as_tibble() %>%
  mutate(dbh_class = 
           cut(dbh_recr, breaks = c(seq(0, 200, by = 10), 1000), labels = FALSE))

recr_data = 
  plot_year_spp_sizeclass %>%
  left_join(
    new_recruits %>%
      select(plot, spp, year = yearin, dbh_class) %>%
      group_by(plot, spp, year, dbh_class) %>%
      summarise(n_recruits = n()) %>%
      ungroup()
  ) %>%
  mutate(n_recruits = ifelse(is.na(n_recruits), 0, n_recruits))


head(recr_data)


recruit_data = 
  plot_year_spp %>%
  left_join(new_recruits %>%
              group_by(plot, spp, yearin) %>%
              summarise(n_recruits = n()) %>%
              ungroup(),
            by = c('plot' = 'plot',
                   'spp' = 'spp',
                   'year' = 'yearin')) %>%
  mutate(n_recruits = ifelse(is.na(n_recruits), 0, n_recruits)) %>%
  left_join(plot_metadata %>%
              select(plot, area_ha)) %>%
  
  # join in BA data
  left_join(plot_ba) %>%
  
  # add WPBR flag
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
  
  # add BURNED5yr
  left_join(fire_years) %>%
  mutate(burned5 = ifelse((!is.na(fire_year1)&year==fire_year1+5)|
                             (!is.na(fire_year2)&year==fire_year2+5),
                           TRUE,
                           FALSE)) %>%
  select(-fire_year1, -fire_year2) %>%
  
  # add cwd
  left_join(cwd_departures %>%
              select(-year_0, -year_1) %>%
              pivot_longer(cols = starts_with('cwddeparture'), 
                           names_to = c('variable', 'year'),
                           values_to = 'cwd_departure', 
                           names_sep = '_') %>%
              mutate(year = as.integer(year)) %>%
              select(plot, year, cwd_departure),
            by = c('plot' = 'plot',
                   'year' = 'year'))
         
head(recruit_data)     



recruit_data.pila = 
  recruits_tally %>%
  filter(spp == 'PILA')

head(recruits_tally)


#### assemble size distribution data ###########################################
# size distribution data: need a row for each plot:year:spp:sizeclass, and 
# column giving the area-standardized rates of occurence (TPH) of trees
tallies_by_size_class = 
  do.call('bind_rows',
          lapply(X = 1:nrow(plot_metadata),
                 FUN = function(p){
                   year_min = 
                     max(plot_metadata$year_min[p], 1999)
                   year_max = plot_metadata$year_max[p]
                   area_ha = plot_metadata$area_ha[p]
                   
                   current_plot_sizedists = 
                     lapply(X = year_min:year_max,
                            FUN = function(y){
                              
                              current_year_sizedist = 
                                census_data %>%
                                filter(status_0=='live' & status_1=='live' & 
                                         year_0 <= y & year_1 >= y & 
                                         plot == plot_metadata$plot[p]) %>%
                                mutate(dbh_class = 
                                         cut(dbh_0,
                                             breaks = c(seq(from = 0, to = 200, by = 10),1000),
                                             labels = FALSE)) %>%
                                group_by(plot, spp, dbh_class) %>%
                                summarise(tph = n() / plot_metadata$area_ha[p]) %>%
                                ungroup() %>%
                                mutate(year = y)
                              
                              
                              return(current_year_sizedist)
                            })
                   
                   return(current_plot_sizedists)
                 }))

head(size_density_data)

head(tallies_by_size_class)

size_density_data = 
  size_density_data %>%
  left_join(tallies_by_size_class) %>%
  mutate(tph = ifelse(is.na(tph), 0, tph))


size_density_data.pila = 
  size_density_data %>%
  filter(spp=='PILA')

#### assemble recruit size data ################################################

head(new_recruits)

new_recruits.pila = 
  new_recruits %>%
  filter(spp == 'PILA')

hist(log(new_recruits.pila$dbh_recr))

#### assemble size metadata ####################################################


size_metadata = 
  census_data %>%
  mutate(dbh_class = cut(dbh_0,
                         breaks = c(seq(from = 0, to = 200, by = 10),1000),
                         labels = FALSE)) %>%
  group_by(dbh_class) %>%
  summarise(dbh_med = median(dbh_0)) %>%
  ungroup() %>%
  mutate(lower = seq(from = 0, to = 200, by = 10),
         upper = c(seq(from = 10, to = 200, by = 10), 1000))


#### assemble model data #######################################################

usgs_model_data = 
  list(P = length(unique(census_data.pila$plot)),
       
       # survival submodel
       N_s = nrow(survival_data.pila),
       surv = survival_data.pila$status_1=='live',
       X_s = 
         survival_data.pila[,c('intercept', 'burned_01', 'wpbr', 'bam2ha_0', 
                               'cwd_departure90')] %>%
         as.matrix(),
       plotid_s = survival_data.pila$plot.i,
       
       # growth model 
       N_g = nrow(growth_data.pila),
       size1_g = growth_data.pila$dbh_1,
       X_g = growth_data.pila[,c('intercept', 'burned_01', 'wpbr', 'bam2ha_0',
                                 'cwd_departure90')] %>%
         as.matrix(),
       plotid_g = growth_data.pila$plot.i
       
       # recruitment model
       ) 


#### write results #############################################################



#### scratch ###################################################################

censuses = 
  usgs %>%
  group_by(PLOT) %>%
  # there are a number of instances where a census was split across multiple years
  # in every case, the majority of the inventory took place in the earliest 
  # year, which is also the year that is on 5-year census cycle, set the min 
  # year as the canonical year of each census
  # specifically: 
    ## GIBBS 1996/1997, 
    ## LOLOG 1984/1985
    ## SUPILA 1982/1983
    ## SURIP 1982/1983
    ## YOHOPIPO 1996/1997
    ## BBBPIPO 2013 and 2014 -> 2012
    ## CCRPIPO 2016 -> 2015
    ## EMSLOPE 2009->2008
    ## EMSLOPE 2014->2013
    ## FFS2BURN 2012->2011
    ## FFS6BURN 2013->2011
    ## FFS7CONTROL 2012->2011
    ## GIBBS 2014->2013 
    ## HOMEPICO 2011->2010
    ## LMCC 2014-> 2013
    ## LOGPIJE 2016->2015
    ## LOLOG 2014->2012
    ## LOTHOR 2014->2013
    ## SURIP 2015->2014
    ## WTABMA 2014->2013
  summarise(YEAR1 = min(YEAR1, na.rm = TRUE),
            YEAR2 = min(YEAR2, na.rm = TRUE),
            YEAR3 = min(YEAR3, na.rm = TRUE),
            YEAR4 = min(YEAR4, na.rm = TRUE),
            YEAR5 = min(YEAR5, na.rm = TRUE),
            YEAR6 = min(YEAR6, na.rm = TRUE),
            YEAR7 = min(YEAR7, na.rm = TRUE),
            YEAR8 = min(YEAR8, na.rm = TRUE)) %>%
  ungroup()

