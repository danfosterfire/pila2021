
library(here)
library(tidyverse)
library(readxl)

seedlings = readxl::read_excel(here::here('02-data','00-source',
                                          'york','2009 compiled.xls'))

head(seedlings)

seedlings %>%
  filter(STATUS=='L' & SPECIES=='SP') %>%
  mutate(dbh_class = cut(`DBH(cm)`, breaks = seq(0, 50.8, 2.54), 
                         labels = FALSE,
                         include.lowest = TRUE)) %>%
  group_by(dbh_class) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(total = seedlings %>%
           filter(STATUS=='L'&SPECIES=='SP') %>%
           nrow()) %>%
  mutate(prop = count / total) %>%
  print()

seedlings %>%
  filter(is.na(`DBH(cm)`) & STATUS=='L')

# so 1/10 of the seedlings get 10 years to grow, 1/10 get 9 years, ... and 1/10 get 1 year

