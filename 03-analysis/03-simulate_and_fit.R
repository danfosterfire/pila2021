
library(here)
library(tidyverse)
library(truncnorm)
library(nimble)

fia_censuses = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'fia_censuses.rds'))

fia.pila = 
  fia_censuses %>%
  filter(spp.re=='PILA')


hist(fia.pila$dbh_in.init)

unique(fia_censuses$tpa_unadj.re)

ggplot(data = fia_censuses,
       aes(x = dbh_in.init))+
  geom_histogram()+
  facet_wrap(~tpa_unadj.init,
             scales = 'free')



#### simulate data #############################################################

N = 1000

subsample = fia_censuses %>% sample_n(size = N)

# vector of sizes at initial census
z_it = subsample$dbh_in.init

sigma = 12 # standard deviation of remeasurement dbh



# coefficient for initial size affecting followup size
alpha_z = mean(fia_censuses$dbh_in.re / fia_censuses$db, na.rm = TRUE)
# shouldnt this have a quadratic effect?

# coefficients for spatial covariates
b_z 



# i length vector of mean remeasurement dbhs
#u_it1 = 


