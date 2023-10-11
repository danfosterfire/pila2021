
#### setup #####################################################################


library(here)
library(tidyverse)
library(sf)
library(DBI)
library(RSQLite)
library(terra)

#### get data used for models ##################################################

model_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'plot_data.rds')) %>%
  mutate(ba_scaled = as.numeric(scale(ba_ft2ac)),
         cwd_dep90_scaled = as.numeric(scale(cwd_departure90)),
         cwd_mean_scaled = as.numeric(scale(cwd_mean)),
         intercept = 1)



#### get extra fiadb info ######################################################

# sqlite dbs downloaded from 
# https://apps.fs.usda.gov/fia/datamart/datamart_sqlite.html
# on 12/14/2021

fiadb.ca = 
  dbConnect(RSQLite::SQLite(),
            here::here('02-data',
                       '00-source',
                       'fia',
                       'FIADB_CA.db'))

fiadb.nv = 
  dbConnect(RSQLite::SQLite(),
            here::here('02-data',
                       '00-source',
                       'fia',
                       'FIADB_NV.db'))
  
fiadb.or = 
  dbConnect(RSQLite::SQLite(),
            here::here('02-data',
                       '00-source',
                       'fia',
                       'FIADB_OR.db'))

fiadb.wa = 
  dbConnect(RSQLite::SQLite(),
            here::here('02-data',
                       '00-source',
                       'fia',
                       'FIADB_WA.db'))

# they're not all the same which is super annoying
dbListTables(fiadb.ca) == dbListTables(fiadb.nv) &
  dbListTables(fiadb.ca) == dbListTables(fiadb.or) &
  dbListTables(fiadb.ca) == dbListTables(fiadb.wa)

dbListTables(fiadb.ca)

dbReadTable(fiadb.ca, 'REF_FIADB_VERSION')

fia = 
  lapply(X =  c('COND', 'PLOT', 'REF_FOREST_TYPE',
                'REF_FOREST_TYPE_GROUP', 'REF_SPECIES', 'SEEDLING', 'SUBPLOT', 
                'SUBP_COND', 'TREE'),
         FUN = function(tname){
           dbReadTable(fiadb.ca, tname) %>%
             bind_rows(dbReadTable(fiadb.nv, tname)) %>%
             bind_rows(dbReadTable(fiadb.or, tname)) %>%
             bind_rows(dbReadTable(fiadb.wa, tname)) %>%
             as_tibble()
         })

names(fia) = c('COND', 'PLOT', 'REF_FOREST_TYPE',
               'REF_FOREST_TYPE_GROUP', 'REF_SPECIES', 'SEEDLING', 'SUBPLOT', 
               'SUBP_COND', 'TREE')
         
dbDisconnect(fiadb.ca)
dbDisconnect(fiadb.or)
dbDisconnect(fiadb.nv)
dbDisconnect(fiadb.wa)

# one row per observation of a subplot, with columns indicating the disturbances 
# present at that observation, the spatial location, etc.
# the goal is to have a table with one row per observation of a condition (stand 
# on a subplot at a specific inventory date)
# context
subplots = 
  
  # start with the plot table, with one row per plot observation
  fia$PLOT %>%
  as_tibble() %>%
  rename(PLT_CN = CN) %>%
  select(PLT_CN, PREV_PLT_CN, INVYR, STATECD, UNITCD, COUNTYCD, PLOT, 
         PLOT_STATUS_CD, PLOT_NONSAMPLE_REASN_CD, MEASYEAR, 
         MEASMON, MEASDAY, KINDCD, DESIGNCD, MANUAL, MACRO_BREAKPOINT_DIA,
         LAT, LON, ELEV, ECOSUBCD) %>%
  
  # join in the subplot info, with one row per subplot observation
  left_join(fia$SUBPLOT %>%
              rename(SUBP_CN = CN) %>%
              select(SUBP_CN, PLT_CN, 
                     STATECD, UNITCD, COUNTYCD, PLOT, SUBP),
            by = c('PLT_CN' = 'PLT_CN',
                   'STATECD' = 'STATECD',
                   'UNITCD' = 'UNITCD',
                   'COUNTYCD' = 'COUNTYCD',
                   'PLOT' = 'PLOT')) %>%
  
  # join the subplot conditions; this creates one row per condition:subplot:plot:time
  left_join(fia$SUBP_COND %>%
              rename(SUBPCOND_CN = CN) %>%
              select(PLT_CN, SUBPCOND_CN, STATECD, UNITCD,
                     COUNTYCD, PLOT, SUBP, CONDID),
            by = c('PLT_CN' = 'PLT_CN', 
                   'STATECD' = 'STATECD',
                   'UNITCD' = 'UNITCD',
                   'COUNTYCD' = 'COUNTYCD',
                   'PLOT' = 'PLOT',
                   'SUBP' = 'SUBP')) %>%
  
  # join the condition details
  left_join(fia$COND %>%
              rename(CND_CN = CN) %>%
              select(PLT_CN, CND_CN, STATECD, UNITCD, 
                     COUNTYCD, PLOT, CONDID, OWNGRPCD, 
                     FORTYPCD, CONDPROP_UNADJ,
                     DSTRBCD1, DSTRBCD2, DSTRBCD3,
                     TRTCD1, TRTCD2, TRTCD3),
            by = c('PLT_CN' = 'PLT_CN',
                   'STATECD' = 'STATECD',
                   'UNITCD' = 'UNITCD',
                   'COUNTYCD' = 'COUNTYCD',
                   'PLOT' = 'PLOT',
                   'CONDID' = 'CONDID')) %>%
  
  filter(PLT_CN %in% model_data$plt_cn)


# how many fortypcalc cds per plot? how many owngrpcds per plot?
swapping_info = 
  subplots %>%
  group_by(PLT_CN, STATECD, COUNTYCD) %>%
  summarise(PRIVATE = any(OWNGRPCD==40),
            FORTYPCD = last(FORTYPCD, order_by = CONDPROP_UNADJ, na_rm = TRUE)) %>%
  ungroup()

swapping_info


#### true data #################################################################

true_data = 
  model_data %>%
  select(plt_cn, plot_id, lat_true = lat, lon_true = lon, 
         cwd_departure_90_true = cwd_departure90, 
         cwd_mean_true = cwd_mean,
         invdate.init, invdate.re) %>%
  left_join(swapping_info,
            by = c('plt_cn' = 'PLT_CN'))


head(true_data)

summary(true_data)

280/941


nrow(swapping_info %>% filter(PRIVATE))

nrow(swapping_info)

#### get 'observed' cwd ########################################################

# fuzz / swap one dataset
fuzz_and_swap = function(dataset){
  
  # first, attempt to swap all of the private plots
  private = 
    dataset %>% 
    # just the private plots
    filter(PRIVATE) %>%
    
    # group by state, county, forest type
    group_by(STATECD, COUNTYCD, FORTYPCD) %>%
    # assign an id
    mutate(swap_id = row_number()) %>%
    ungroup() %>%
    # get the partner id
    mutate(partner_id = 
             case_when(swap_id%%2==0~swap_id-1,
                       swap_id%%2==1~swap_id+1),
           pair_id = paste0(STATECD,'-',COUNTYCD,'-',FORTYPCD,'-',floor((swap_id+1)/2))) %>%
    
    # left join the partner's coordinates
    left_join(y = select(., lat_swapped = lat_true, lon_swapped = lon_true, 
                         STATECD, COUNTYCD, FORTYPCD, pair_id, partner_id),
              by = 
                c('STATECD' = 'STATECD',
                  'COUNTYCD' = 'COUNTYCD',
                  'FORTYPCD' = 'FORTYPCD',
                  'pair_id' = 'pair_id',
                  'swap_id' = 'partner_id'))

  n_private = nrow(private)
  
  pairs_to_swap = 
    private %>%
    filter(!is.na(lat_swapped)) %>%
    distinct(pair_id) %>%
    slice_sample(n = floor(0.1*n_private)) %>%
    pull(pair_id)
  
  # then, select only the ones where swapping was actually successful (ie, it 
  # had a partner), and keep
  # a number of them equal to 20% of n_private
  private_swapped = private %>% filter(pair_id %in% pairs_to_swap)
  
  # dothe swapping step, for plots that got swapped, or leave them in place
  # for the plots that didnt
  swapped = 
    dataset %>%
    left_join(private_swapped %>%
                select(plt_cn, lat_swapped, lon_swapped),
              by = c('plt_cn' = 'plt_cn')) %>%
    mutate(lat_afterwap = ifelse(is.na(lat_swapped), lat_true, lat_swapped),
           lon_afterswap = ifelse(is.na(lon_swapped), lon_true, lon_swapped))
  
  # for each point, get a random bearing and normally distributed distance to 
  # fuzz
  fuzzed = 
    swapped %>%
    # perturbing by a half normal with sd = 804, so 84 percent of plots perturbed 
    # "within 0.5 mile for most plots" and 98% of plots perterbued "up to 1.0 mile 
    # on a small subset of them"
    mutate(distance_m = abs(rnorm(n = nrow(.), sd = 804, mean = 0)),
           distance_m = ifelse(distance_m>1608, 1608, distance_m),
           bearing_deg = runif(n = nrow(.), min = 0, max = 360)) %>%
    
    rowwise() %>%
    mutate(
           lat_afterfuzz = 
             geosphere::destPoint(p = c(lon_afterswap, lat_afterwap),
                                  b = bearing_deg,
                                  d = distance_m)[,'lat'],
           
           lon_afterfuzz = 
             geosphere::destPoint(p = c(lon_afterswap, lat_afterwap),
                                  b = bearing_deg,
                                  d = distance_m)[,'lon'])
  
  result = 
    fuzzed %>%
    select(plt_cn, lat_fs= lat_afterfuzz, lon_fs = lon_afterfuzz)
  
  return(result)
  
  }


extract_cwd_metrics = 
  function(plot_data){
    
    plot_data.sf = 
      plot_data %>%
      sf::st_as_sf(coords = c('lon_fs', 'lat_fs'),
                   crs = 'EPSG:4269')
    
    plots_bbox = 
      list('lat_min' = min(plot_data$lat_fs)-1,
           'lat_max' = max(plot_data$lat_fs)+1,
           'lon_min' = min(plot_data$lon_fs)-1,
           'lon_max' = max(plot_data$lon_fs)+1)
    
    cwd_growseason_means = 
      lapply(X = 2000:2020,
             FUN = function(y){
               
               cwd_year = 
                 rast(here::here('02-data',
                                 '00-source',
                                 'terraclimate',
                                 paste0('TerraClimate_def_',y,'.nc')))
               
               cwd_year = 
                 crop(cwd_year,
                      c(plots_bbox$lon_min, plots_bbox$lon_max,
                        plots_bbox$lat_min, plots_bbox$lat_max))
               
               cwd_year = 
                 mean(cwd_year[[5:10]])
               
               return(cwd_year)
               
             }) %>%
      rast()
    
    names(cwd_growseason_means) = 
      paste0('growseasonmean_',as.character(2000:2020))
    
    # get departure from "normal" (20 year mean) CWD for each 
    # year on each location
    cwd_departure = 
      cwd_growseason_means - mean(cwd_growseason_means)
    
    names(cwd_departure) = paste0('cwddeparture_', as.character(2000:2020))
    
    cwd_departures = 
      extract(cwd_departure, plot_data[,c('lon_fs', 'lat_fs')]) %>%
      bind_cols('plot_id' = plot_data$plot_id,
                'year_begin' = lubridate::year(plot_data$invdate.init),
                'year_end' = lubridate::year(plot_data$invdate.re))
    
    cwd_departure_span = 
      sapply(X = 1:nrow(cwd_departures),
             FUN = function(i){
               
               year_range = seq(from = cwd_departures$year_begin[i],
                                to = cwd_departures$year_end[i],
                                by = 1)
               
               columns_to_select = 
                 paste0('cwddeparture_', year_range)
               
               values_in_span = 
                 as.numeric(cwd_departures[i,columns_to_select])
               
               # return the 90th percentile CWD departure
               return(as.numeric(quantile(values_in_span, probs = 0.9)))
               
             })
    
    plot_data$cwd_departure90_fs = cwd_departure_span
    
    plot_data$cwd_mean_fs = 
      extract(mean(cwd_growseason_means),
              plot_data[,c('lon_fs', 'lat_fs')])$mean
  
    return(plot_data)
      
  }





#### compare true with fuzzed/swapped CWD ######################################

set.seed(110819)

# fuzz/swap one dataset
single_fs = 
  true_data %>%
  left_join(
    fuzz_and_swap(dataset = true_data),
    by = c('plt_cn' = 'plt_cn')) %>%
  extract_cwd_metrics() %>%
  
  mutate(cwd_mean_fs_scaled = as.numeric(scale(cwd_mean_fs)),
         cwd_mean_true_scaled = as.numeric(scale(cwd_mean_true)),
         cwd_departure90_fs_scaled = as.numeric(scale(cwd_departure90_fs)),
         cwd_departure90_true_scaled = as.numeric(scale(cwd_departure_90_true)))


ggplot(single_fs,
       aes(x = cwd_mean_fs_scaled, y = cwd_mean_true_scaled))+
  geom_point(alpha = 0.5)+
  theme_minimal()

ggplot(single_fs,
       aes(x = cwd_departure90_fs_scaled, y = cwd_departure90_true_scaled))+
  geom_point(alpha = 0.5)+
  theme_minimal()

ggplot(single_fs,
       aes(x= cwd_mean_fs_scaled-cwd_mean_true_scaled))+
  geom_histogram()+
  theme_bw()

ggplot(single_fs,
       aes(x= cwd_departure90_fs_scaled-cwd_departure90_true_scaled))+
  geom_histogram()+
  theme_bw()


#### do this 1000 times ########################################################

fuzz_swap_one = 
  function(id){
    
    print(paste0('Working on ID ', id))
    
    single_fs = 
      true_data %>%
      left_join(
        fuzz_and_swap(dataset = true_data),
        by = c('plt_cn' = 'plt_cn')) %>%
      extract_cwd_metrics() %>%
      
      mutate(cwd_mean_fs_scaled = as.numeric(scale(cwd_mean_fs)),
             cwd_mean_true_scaled = as.numeric(scale(cwd_mean_true)),
             cwd_departure90_fs_scaled = as.numeric(scale(cwd_departure90_fs)),
             cwd_departure90_true_scaled = as.numeric(scale(cwd_departure_90_true)),
             simid = id)
    
    return(single_fs)
  }

fuzzed_swapped_1k = 
  do.call('bind_rows',
          lapply(X = 1:1000,
                 function(x){fuzz_swap_one(x)}))

results = 
  fuzzed_swapped_1k %>%
  mutate(error_cwd_mean = cwd_mean_fs_scaled-cwd_mean_true_scaled,
         error_cwd_departure90 = cwd_departure90_fs_scaled-cwd_departure90_true_scaled)

summary(results %>% select(error_cwd_mean, error_cwd_departure90))


results %>%
  select(error_cwd_mean, error_cwd_departure90) %>%
  summarise_all(function(x){mean(abs(x))})

results %>%
  select(error_cwd_mean, error_cwd_departure90) %>%
  summarise_all(function(x){quantile(abs(x), c(0.95))})

cwd_mean_plot = 
  ggplot(results,
       aes(x = error_cwd_mean))+
  geom_histogram()+
  theme_bw()+
  labs(x = 'Error introduced by fuzzing/swapping\n(scaled units of mean CWD)',
       y = 'Count of simulated plots')

ggsave(plot = cwd_mean_plot,
       filename = here::here('04-communication', 'figures', 'manuscript', 'fuzzswap_cwd_mean.png'),
       height = 5, width = 5)

cwd_dep90_plot = 
  ggplot(results,
         aes(x = error_cwd_departure90))+
  geom_histogram()+
  theme_bw()+
  labs(x = 'Error introduced by fuzzing/swapping\n(scaled units of 90th percentile CWD departure)',
       y = 'Count of simulated plots')

ggsave(plot = cwd_dep90_plot,
       filename = here::here('04-communication', 'figures', 'manuscript', 'fuzzswap_cwd_dep90.png'),
       height = 5, width = 5)


error_plot = 
  results %>%
  pivot_longer(cols = c('error_cwd_mean', 'error_cwd_departure90'),
               names_to = 'variable',
               values_to = 'value') %>%
  mutate(explanatory_variable = case_when(variable=='error_cwd_mean'~'Site Dryness\n(mean CWD)',
                                          variable=='error_cwd_departure90'~'Drought\n(90th percentile CWD departure)')) %>%
  ggplot(data = .,
         aes(x=value))+
  geom_histogram()+
  theme_bw()+
  facet_grid(.~explanatory_variable)+
  labs(x = 'Error introduced by fuzzing/swapping (standard deviations)',
       y = 'Count of simulated plots')

error_plot
ggsave(plot = error_plot,
       filename = here::here('04-communication', 'figures', 'manuscript', 'fuzzswap_both.png'),
       height = 5, width = 7)
