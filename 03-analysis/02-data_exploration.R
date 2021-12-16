
library(here)
library(tidyverse)
library(ggplot2)


growth_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'growth_data.rds'))

mort_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'mort_data.rds'))

sizedist_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'sizedist_data.rds'))

subplot_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'subplot_data.rds'))

#### growth data ###############################################################

growth_data %>%
  print(width = Inf)

unique(growth_data$tree_status.init)
unique(growth_data$tree_status.re)

nrow(growth_data) == length(unique(growth_data$tre_cn))

ggplot(data = growth_data,
       aes(x = dbh_in.init))+
  geom_histogram()

ggplot(data = growth_data,
       aes(x = dbh_in.re))+
  geom_histogram()

growth_data %>% filter(is.na(dbh_in.re)) %>% print(width = Inf)

summary(growth_data$dbh_in.re)

ggplot(data = growth_data,
       aes(x = dbh_in.re-dbh_in.init))+
  geom_histogram()

ggplot(data = growth_data,
       aes(x = species))+
  geom_bar()

growth_data %>% filter(dbh_in.re-dbh_in.init < -50) %>% print(width = Inf)

ggplot(data = growth_data,
       aes(x = dbh_in.init, color = species))+
  geom_density()+
  theme_minimal()

ggplot(data = growth_data,
       aes(x = dbh_in.re, color = species))+
  geom_density()+
  theme_minimal()

ggplot(data = growth_data,
       aes(x = dbh_in.re-dbh_in.init, color = species))+
  geom_density()+
  theme_minimal()+
  scale_x_continuous(limits = c(-5, 10))



#### mortality data ############################################################

unique(mort_data$tree_status.init)
unique(mort_data$tree_status.re)

nrow(mort_data) == length(unique(mort_data$tre_cn))

ggplot(data = mort_data,
       aes(x = survived))+
  geom_bar()

ggplot(data = mort_data,
       aes(x = tree_status.re))+
  geom_bar()

ggplot(data = mort_data,
       aes(x = species))+
  geom_bar()

ggplot(data = mort_data,
       aes(x = species, fill = survived))+
  geom_bar()

ggplot(data = mort_data,
       aes(x = species, fill = survived))+
  geom_bar(position = position_fill())

#### sizedist data #############################################################

head(sizedist_data)

summary(sizedist_data)

ggplot(data = 
         sizedist_data %>%
         group_by(subp_id, dbh_class) %>%
         summarise(tpa_unadj.init = sum(tpa_unadj.init)) %>%
         ungroup(),
       aes(x = factor(dbh_class), y = tpa_unadj.init))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 100))

# this LOOKS like a lot of zeros, but the bins are so small that most of them 
# are empty the vast majority of the time. spot-checking individual plots confirms 
# that trees are indeed showing up in the size distribution bins correctly


#### subplot data ##############################################################

head(subplot_data)

lapply(X = 
         c('elev_ft', 'lat', 'lon', 'ba_ft2ac', 'max_cwd_departure'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = .data[[v]]))+
           geom_histogram()
       })


lapply(X = 
         c('fire', 'insects', 'disease', 'cutting'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = .data[[v]]))+
           geom_bar()
       })

# plot variables across space
lapply(X = 
         c('fire', 'insects', 'disease', 'ba_ft2ac', 'max_cwd_departure'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = lon, y = lat, color = .data[[v]]))+
           geom_point(alpha = 0.5, size = 0)+
           coord_sf()+
           theme_minimal()
       })

# plot variables across time
lapply(X = 
         c('fire', 'insects', 'disease'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = factor(lubridate::year(invdate.re)),
                    fill = .data[[v]]))+
           geom_bar()
       })

lapply(X = 
         c('ba_ft2ac', 'max_cwd_departure'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = factor(lubridate::year(invdate.re)),
                    y = .data[[v]]))+
           geom_boxplot()
       }) # a little confusing; keep in mind that cwd_departure is over the previous 10 years, not 
# the year shown

# plot them across eachother
ggplot(data = subplot_data,
       aes(x = fire, y = ba_ft2ac))+
  geom_boxplot() # huh interesting looks like higher BA plots were more likely to burn

ggplot(data = subplot_data,
       aes(x = fire, y = max_cwd_departure))+
  geom_boxplot() # not much relationship between drought and fire

ggplot(data = subplot_data,
       aes(x = insects, y = ba_ft2ac))+
  geom_boxplot() # higher ba more likely to have insects?

ggplot(data = subplot_data,
       aes(x = insects, y = max_cwd_departure))+
  geom_boxplot() # higher drought more likely to have insects

ggplot(data = subplot_data,
       aes(x = disease, y = ba_ft2ac))+
  geom_boxplot() # higher BA more likley to have disease

ggplot(data = subplot_data,
       aes(x = disease, y = max_cwd_departure))+
  geom_boxplot() # drought decreases probability of disease

ggplot(data = subplot_data,
       aes(x = ba_ft2ac, y = max_cwd_departure))+
  geom_point()+
  geom_smooth(method = 'lm') # very weak positive relationship between BA and CWD





#### notes #####################################################################

# issues to go clean up in the data import script
#  - 17 rows in growth data with NA remeasure DBH; probably just drop them
#  - a few crazy growth deltas (ignore? thats what perry would want probably)
# - think about what to do with the "cutting" subplot flag, which is all FALSE
#   - consider a better measure of CWD than max; maybe 90th percentile? or max of a moving window?