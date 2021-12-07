#### setup #####################################################################

library(here)
library(raster)
library(sf)
library(tidyverse)
library(spData)
library(units)


# CA, OR borders
aoi = 
  spData::us_states %>%
  filter(is.element(NAME, c('California', 'Oregon'))) %>%
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
  raster::aggregate(fact = 12) %>%
  raster::projectRaster(from = ., crs = st_crs(aoi)$input)

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

st_write(pila_range.sf,
         here::here('02-data',
                    '01-preprocessed',
                    'pila_range_map.shp'),
         delete_dsn=TRUE)

