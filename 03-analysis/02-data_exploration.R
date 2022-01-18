
library(here)
library(tidyverse)
library(ggplot2)


growth_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'growth_data.rds'))

mort_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'mort_data.rds'))

sizedist_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'sizedist_data.rds'))

subplot_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
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
         c('elev_ft', 'lat', 'lon', 'ba_ft2ac', 'cwd_departure90', 'cwd_mean'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = .data[[v]]))+
           geom_histogram()
       })


lapply(X = 
         c('fire', 'insects', 'disease', 'cutting', 'wpbr'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = .data[[v]]))+
           geom_bar()
       })

# plot variables across space
lapply(X = 
         c('fire', 'insects', 'disease','wpbr', 'ba_ft2ac', 'cwd_departure90', 'cwd_mean'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = lon, y = lat, color = .data[[v]]))+
           geom_point(alpha = 0.25)+
           coord_sf()+
           theme_minimal()
       })

# plot variables across time
lapply(X = 
         c('fire', 'insects', 'disease', 'wpbr'),
       FUN = function(v){
         ggplot(data = subplot_data,
                aes(x = factor(lubridate::year(invdate.re)),
                    fill = .data[[v]]))+
           geom_bar()
       })

lapply(X = 
         c('ba_ft2ac', 'cwd_departure90'),
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
       aes(x = fire, y = cwd_departure90))+
  geom_boxplot() # fire more likely with drought

ggplot(data = subplot_data,
       aes(x = insects, y = ba_ft2ac))+
  geom_boxplot() # higher ba more likely to have insects?

ggplot(data = subplot_data,
       aes(x = insects, y = cwd_departure90))+
  geom_boxplot() # not much effect of drought on insects 
# (there was with max departure instead of 90th)

ggplot(data = subplot_data,
       aes(x = disease, y = ba_ft2ac))+
  geom_boxplot() # higher BA more likley to have disease

ggplot(data = subplot_data,
       aes(x = disease, y = cwd_departure90))+
  geom_boxplot() # drought decreases probability of disease 
# (detection? or diseases actually spread more under damp conditions?)

ggplot(data = subplot_data,
       aes(x = ba_ft2ac, y = cwd_departure90))+
  geom_point()+
  geom_smooth(method = 'lm') # none

ggplot(data = subplot_data,
       aes(y = ba_ft2ac, x = cwd_mean))+
  geom_point()+
  geom_smooth(method = 'lm') # ba is lower at drier sites

ggplot(data = subplot_data,
       aes(x = fire, y = cwd_mean))+
  geom_boxplot() # not much relationship between dryness and fire

ggplot(data = subplot_data,
       aes(x = insects, y = cwd_mean))+
  geom_boxplot() # more insects at drier sites

ggplot(data = subplot_data,
       aes(x = disease, y= cwd_mean))+
  geom_boxplot() # more disease at drier sites

ggplot(data = subplot_data,
       aes(x = fire, fill = insects))+
  geom_bar(position = position_fill()) # fire doesnt affect the p(insects)

ggplot(data = subplot_data,
       aes(x = insects, fill = fire))+
  geom_bar(position = position_fill()) # insects dont affect p(fire)

ggplot(data = subplot_data,
       aes(x = fire, fill = disease))+
  geom_bar(position = position_fill()) # disease more likely when fire = false

ggplot(data = subplot_data,
       aes(x = disease, fill = fire))+
  geom_bar(position = position_fill()) # fire less likely when diease = TRUE

ggplot(data = subplot_data,
       aes(x = insects, fill = disease))+
  geom_bar(position = position_fill()) # disease way more likely when insects true




#### notes #####################################################################

# issues to go clean up in the data import script
#  - 17 rows in growth data with NA remeasure DBH; probably just drop them (done)
#  - a few crazy growth deltas (leave them in? thats what perry would want probably; done)
# - think about what to do with the "cutting" subplot flag, which is all FALSE 
#    (leave it in, ditch from model)
#   - consider a better measure of CWD than max; maybe 90th percentile? or max 
#     of a moving window? (90th percentile departure + mean cwd as separate variables)