#### setup #####################################################################

library(here)
library(raster)
library(sf)
library(tidyverse)
library(spData)
library(units)
library(USAboundaries)
library(ggspatial)
library(cowplot)
library(DBI)
library(RSQLite)


# CA, OR borders
aoi = 
  spData::us_states %>%
  filter(is.element(NAME, c('California', 'Oregon', 'Nevada', 'Washington'))) %>%
  summarise() %>%
  st_buffer(dist = units::as_units(100, 'meters'))

# pila range from Live tree species basal area of the contiguous United STates (2000-2009)
# see https://www.fs.usda.gov/rds/archive/catalog/RDS-2013-0013
# Wilson, Barry Tyler; Lister, Andrew J.; Riemann, Rachel I.; Griffith, 
# Douglas M. 2013. Live tree species basal area of the contiguous United 
# States (2000-2009). Newtown Square, PA: USDA Forest Service, 
# Rocky Mountain Research Station. https://doi.org/10.2737/RDS-2013-0013
pila_range = 
  raster::raster(here::here('02-data',
                            '00-source',
                            'wilson2013',
                            's117.img')) %>%
  # clip, aggregate to bigger cells, and then reproject
  raster::crop(x = ., y = st_transform(aoi, crs(.))) %>%
  
  # 1 degree of latitude is ~ 69.2 miles, so to get pixels ~2 miles across (twice 
  # the distance of the FIA "fuzzing") we want pixels to be (2/69.2) =~ 0.029 decimal degrees, 
  # so a factor of 14
  raster::aggregate(fact = 14) %>%
  raster::projectRaster(from = ., crs = st_crs(aoi)$input)

pila_range

plot(pila_range)

#### optimize BA threshold #####################################################

# FIA locations in CA/OR/WA/NV, with a code indicating whether PILA is 
# present or absent

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

# START HERE; GET A TABLE WITH ONE ROW PER PLOT AND A COLUMN INDICATING WHETEHTER 
# THERE WAS A LIVE PILA AT ANY SAMPLING OF THE PLOT
pila_presabs = 
  fia$PLOT %>%
  left_join(fia$TREE %>%
              select(PLT_CN, STATUSCD, SPCD), 
            by = c('CN' = 'PLT_CN')) %>%
  mutate(live_pila = ifelse(STATUSCD==1&SPCD==117, TRUE, FALSE),
         plot_id = paste(STATECD, UNITCD, COUNTYCD, PLOT, 
                         sep = '-')) %>%
  select(plot_id, lat = LAT, lon = LON, live_pila) %>%
  group_by(plot_id, lat, lon) %>%
  summarise(live_pila = any(live_pila, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(lat))


pila_presabs$predicted_ba = 
  raster::extract(pila_range, pila_presabs[,c('lon', 'lat')])


# get the accuracy implied by a specific threshold
test_threshold = 
  function(params){
    
    threshold = exp(params[1])
    
    # use the threshold to predict presence
    pila_presabs = 
      pila_presabs %>%
      
      mutate(
        predicted_presence = ifelse(predicted_ba>=threshold, TRUE, FALSE),
        true_positive = 
          ifelse(predicted_presence==TRUE & live_pila==TRUE, TRUE, FALSE),
        true_negative = 
          ifelse(predicted_presence==FALSE & live_pila == FALSE, TRUE, FALSE),
        false_positive = 
          ifelse(predicted_presence==TRUE & live_pila == FALSE, TRUE, FALSE),
        false_negative = 
          ifelse(predicted_presence == FALSE & live_pila==TRUE, TRUE, FALSE)
      )
    
    summary_stats = 
      pila_presabs %>%
      summarise(
        accuracy = sum(true_positive+true_negative)/n(),
        misclass = sum(false_negative+false_positive)/n(),
        tpr = sum(true_positive)/sum(true_positive+false_negative),
        precision = sum(true_positive)/sum(true_positive+false_positive),
        tnr = sum(true_negative)/sum(true_negative+false_positive),
        fpr = sum(false_positive)/sum(true_negative+false_positive)
      )
    
    return(summary_stats$accuracy[1])
  }

# optimize the threshold value
optim_results = 
  optim(par = c(1), 
      fn = test_threshold, 
      control = list('fnscale' = -1),
      method = 'Brent',
      lower = -10, 
      upper = 10)

optim_results

test_raster = 
  pila_range>=exp(optim_results$par)

plot(test_raster)


# filter to only cells where the predicted BA/ac is greater than or equal to 2
# this gives a level of inclusivity which is about on par with the little polygons
# filtering to to include anything > 0 gives too inclusive an area, including many 
# plots where PILA presence isn't really plausible up in WA
pila_range[pila_range<2] = NA 

# turn the raster into a multipart polygon
pila_range.sf = 
  raster::rasterToPolygons(x = pila_range>=2) %>%
  sf::st_as_sf() %>%
  mutate(area_m2 = st_area(.)) %>%
  summarise(area_m2 = sum(area_m2)) %>%
  mutate(area_km2 = set_units(area_m2, 'km^2')) %>%
  st_intersection(aoi) %>%
  dplyr::select(area_km2)

plot(pila_range.sf)

test = spData::world



overview_map = 
  ggplot(data = 
         spData::world %>%
         filter(continent=='North America'))+
  geom_sf(fill = NA, lwd = 1)+
  geom_rect(xmin = 390000, xmax = 1100000, ymin = 3740000, ymax = 5020000,
            fill = NA, color = 'red', lwd = 2)+
  coord_sf(crs = "EPSG:26910",
           xlim = c(200000, 5000000), ylim = c(1000000, 8000000))+
  theme_minimal()+
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

overview_map

range_map = 
  ggplot(data = pila_range.sf)+
  geom_sf(data = 
            us_states(),
          fill = NA, lwd = 1)+
  geom_sf(color = 'black', fill = '#20A387FF')+
  theme_minimal()+
  coord_sf(xlim = c(390000, 1100000), ylim = c(3740000, 5020000),
           crs = "EPSG:26910")+
  annotation_scale()

range_map


fire_perims = 
  st_read(here::here('02-data',
                     '00-source',
                     'usfs',
                     'S_USA.FinalFirePerimeter',
                     'S_USA.FinalFirePerimeter.shp')) %>%
  filter(FIREYEAR>=2011)

head(fire_perims)

range_map_fire = 
  ggplot(data = pila_range.sf)+
  geom_sf(data = 
            us_states(),
          fill = NA, lwd = 1)+
  geom_sf(color = 'black', fill = '#20A387FF')+
  geom_sf(data = fire_perims, color = NA, fill = 'red')+
  theme_minimal()+
  coord_sf(xlim = c(390000, 1100000), ylim = c(3740000, 5020000),
           crs = "EPSG:26910")+
  annotation_scale()

range_map_fire

pila_rangemap = 
  ggdraw()+
  draw_plot(range_map)+
  draw_plot(overview_map,
            x = 0.48, y = 0.6, width = 0.3, height = 0.3)

ggsave(plot = pila_rangemap,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'pila_range_map.png'))


st_write(pila_range.sf,
         here::here('02-data',
                    '01-preprocessed',
                    'pila_range_map.shp'),
         delete_dsn=TRUE)

