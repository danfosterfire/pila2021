# an earlier version of this code included insects and their interaction with DBH
# as potential explanatory variables. These were dropped because this script 
# revealed too few observations of "insects" to reliably include them in the 
# model. Conceptually, the paper will treat insects and drought as 
# confounded, which they are to a large extent in reality.

#### setup #####################################################################

library(here)
library(tidyverse)
library(ggplot2)


# load data
mort_list = readRDS(here::here('02-data',
                               '02-for_analysis',
                               'pila_mort_training.rds'))

growth_list = readRDS(here::here('02-data',
                                 '02-for_analysis',
                                 'pila_growth_training.rds'))

recr_list = readRDS(here::here('02-data',
                               '02-for_analysis',
                               'pila_fecd_training.rds'))


#### growth data ###############################################################

growth_data = 
  growth_list$X %>%
  as_tibble()

growth_data$size1 = growth_list$size1


names(growth_data)
# plot X distributions
lapply(X = c('dbh_m.init','ba_scaled',
             'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = .data[[v]]))+
           geom_histogram()
       })
lapply(X = c('fire', 'wpbr'),
       FUN = function(v){
  ggplot(data = growth_data,
         aes(x = .data[[v]]))+
    geom_bar()
})


# plot Y distribution
response_growth = 
  ggplot(data = growth_data,
       aes(x = size1))+
  geom_histogram()+
  theme_minimal()+
  theme(text = element_text(size = 14))+
  labs(y = 'N trees', x = 'DBH at remeasurement (cm)')

response_growth

ggsave(response_growth,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'response_growth.png'),
       height = 2.5, width = 5.5, units = 'in')

# plot XX distributions
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v_i){
         lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
                FUN = function(v_j){
                  ggplot(data = growth_data,
                         aes(x = .data[[v_i]], y = .data[[v_j]]))+
                    geom_point()
                })
       })
# looks fine
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v_i){
         lapply(X = c('fire', 'wpbr'),
                FUN = function(v_j){
                  ggplot(data = growth_data,
                         aes(x = as.factor(.data[[v_j]]), y = .data[[v_i]]))+
                    geom_violin()+
                    geom_jitter(width = 0.1, size = 0, height = 0)+
                    geom_boxplot(width = 0.1, position= position_nudge(x = -0.25))
                })
       })
ggplot(data = growth_data,
       aes(x = as.factor(fire), fill = as.factor(wpbr)))+
  geom_bar(position = position_fill())
# fire is less common on wpbr plots; this should be fine I have a lot of plots with either 
# wpbr or fire, and a few with both


# plot XY distributions
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = .data[[v]], y = size1))+
           geom_point()
       })
lapply(X = c('fire', 'wpbr'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = as.factor(.data[[v]]), y = size1))+
                    geom_violin()+
                    geom_jitter(width = 0.1, size = 0, height = 0)+
                    geom_boxplot(width = 0.1, position= position_nudge(x = -0.25))
       })

# plot XXY distributions for size with other vars
lapply(X = c('ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = dbh_m.init, 
                    y = size1, 
                    color = cut(.data[[v]],
                               breaks = quantile(.data[[v]], probs = c(0, 0.33, 0.66, 1)),
                               include_lowest = TRUE)))+
           geom_point()+
           geom_smooth(method = 'lm')+
           theme(legend.position = 'bottom')
       })
lapply(X = c('fire', 'wpbr'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = dbh_m.init,
                    y = size1,
                    color = as.factor(.data[[v]])))+
           geom_point()+
           geom_smooth(method = 'lm')+
           theme(legend.position = 'bottom')
       })

#### this one ##################################################################
ggplot(data = growth_data,
       aes(x = dbh_m.init, y = size1-dbh_m.init))+
  geom_point(size = 1)+
  geom_smooth(method = 'lm',
              formula = y~x+I(x**2))+
  theme_minimal()

#### survival data #############################################################

surv_data = 
  mort_list$X %>%
  as_tibble()

surv_data$surv = mort_list$surv

names(surv_data)
# plot X distributions
lapply(X = c('dbh_m.init','ba_scaled',
             'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(x = .data[[v]]))+
           geom_histogram()
       })

explanatory_size = 
  ggplot(data = surv_data,
         aes(x = dbh_m.init))+
  geom_histogram()+
  theme_minimal()+
  theme(text = element_text(size = 14))+
  labs(x = 'DBH (m)', y = 'N trees')

explanatory_size

ggsave(explanatory_size,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'explanatory_size.png'),
       height = 4.5, width = 5.5, units = 'in')

explanatory_ba = 
  ggplot(data = surv_data,
         aes(x = ba_scaled))+
  geom_histogram()+
  theme_minimal()+
  theme(text = element_text(size = 14))+
  labs(x = 'Basal area (scaled)', y = 'N trees')

explanatory_ba

ggsave(explanatory_ba,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'explanatory_ba.png'),
       height = 4.5, width = 5.5, units = 'in')

explanatory_drought = 
  ggplot(data = surv_data,
         aes(x = cwd_dep90_scaled))+
  geom_histogram()+
  theme_minimal()+
  theme(text = element_text(size = 14))+
  labs(x = 'Drought (scaled CWD departure)', y = 'N trees')

ggsave(explanatory_drought,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'explanatory_drought.png'),
       height = 4.5, width = 5.5, units = 'in')

explanatory_dryness = 
  ggplot(data = surv_data,
         aes(x = cwd_mean_scaled))+
  geom_histogram()+
  theme_minimal()+
  theme(text = element_text(size = 14))+
  labs(x = 'Site dryness (scaled mean CWD)', y = 'N trees')

ggsave(explanatory_dryness,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'explanatory_dryness.png'),
       height = 4.5, width = 5.5, units = 'in')




lapply(X = c('fire', 'wpbr'),
       FUN = function(v){
  ggplot(data = surv_data,
         aes(x = .data[[v]]))+
    geom_bar()
})

explanatory_fire = 
  ggplot(data = surv_data,
         aes(x = as.factor(fire)))+
  geom_bar()+
  theme_minimal()+
  labs(y = 'N trees', x = 'Fire')+
  theme(text = element_text(size = 14))

ggsave(explanatory_fire,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'explanatory_fire.png'),
       height = 4.5, width = 5.5, units = 'in')

explanatory_wpbr = 
  ggplot(data = surv_data,
         aes(x = as.factor(wpbr)))+
  geom_bar()+
  theme_minimal()+
  labs(y = 'N trees', x = 'WPBR')+
  theme(text = element_text(size = 14))

ggsave(explanatory_wpbr,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'explanatory_wpbr.png'),
       height = 4.5, width = 5.5, units = 'in')


# too few instances of insects = 1

# plot Y distribution
response_surv = 
  ggplot(data = surv_data,
       aes(x = as.factor(surv)))+
  geom_bar()+
  theme_minimal()+
  theme(text = element_text(size = 14))+
  labs(y = 'N trees', x = 'Survived')

response_surv

ggsave(response_surv,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'response_surv.png'),
       height = 2.5, width = 5.5, units = 'in')

# plot XX distributions
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v_i){
         lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
                FUN = function(v_j){
                  ggplot(data = surv_data,
                         aes(x = .data[[v_i]], y = .data[[v_j]]))+
                    geom_point()
                })
       })
# looks fine
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v_i){
         lapply(X = c('fire', 'wpbr'),
                FUN = function(v_j){
                  ggplot(data = surv_data,
                         aes(x = as.factor(.data[[v_j]]), y = .data[[v_i]]))+
                    geom_violin()+
                    geom_jitter(width = 0.1, size = 0, height = 0)+
                    geom_boxplot(width = 0.1, position= position_nudge(x = -0.25))
                })
       })
ggplot(data = surv_data,
       aes(x = as.factor(fire), fill = as.factor(wpbr)))+
  geom_bar(position = position_fill())
# fire is less common on wpbr plots; this should be fine I have a lot of plots with either 
# wpbr or fire, and a few with both


# plot XY distributions
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(y = .data[[v]], x = as.factor(surv)))+
           geom_violin()+
           geom_jitter(width = 0.1, size = 0, height = 0)+
           geom_boxplot(width = 0.1, position = position_nudge(x = -0.25))
       })

#### this one ###################################################
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(x = .data[[v]], y = surv))+
           geom_jitter(width = 0, size = 0, height = 0.1)+
           geom_smooth(method = 'glm', formula = y~x+I(x**2),
                       method.args = list(family = 'binomial'))
       })

surv_data$ecosub.i = mort_list$ecosub_s

surv_data %>%
  filter(dbh_m.init>2) %>%
  print(width = Inf)

lapply(X = c('fire', 'wpbr'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(x = as.factor(.data[[v]]), fill = as.factor(surv)))+
                    geom_bar(position = position_fill())
       })

# plot XXY distributions for size with other vars
lapply(X = c('ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(x = dbh_m.init, 
                    y = surv, 
                    color = cut(.data[[v]],
                               breaks = c(-10, -1, 1, 10),
                               include_lowest = TRUE)))+
           geom_jitter(width = 0, height = 0.1, size = 0)+
           geom_smooth(method = 'glm', method.args = list(family = 'binomial'))+
           theme(legend.position = 'bottom')
       })
lapply(X = c('fire', 'wpbr'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(x = dbh_m.init,
                    y = surv,
                    color = as.factor(.data[[v]])))+
           geom_jitter(width = 0, height = 0.1, size = 0)+
           geom_smooth(method = 'glm', method.args = list(family = 'binomial'))+
           theme(legend.position = 'bottom')
       })

#### recr data ################################################################

# plot y distribution
hist(recr_list$cprime) # oof i am pretty skeptical this is gonna work



# check total TPA on each subplot looks OK
recr_tpa = 
  sapply(X = 1:ncol(recr_list$n),
         FUN = function(plot){
           subplot_tpa = sum(recr_list$n[,plot])
           return(subplot_tpa)
         })

ncol(recr_list$n)

hist(recr_tpa)
summary(recr_tpa/0.404686)


# start with the size distribution data for all 1221 plots where pila was 
# present at initial or followup
comparison_tph = 
  readRDS(here::here('02-data',
           '01-preprocessed',
           'sizedist_data.rds')) %>%
  # keep only the 967 plots which are in the recruitment training set
  filter(
    is.element(plot_id,
                 readRDS(here::here('02-data',
                   '02-for_analysis',
                   'union_plots.rds')) %>%
                 filter(is.element(plot_id.i,
                                   recr_list$plot_id)) %>%
                 pull(plot_id))
    ) %>%
  filter(species=='PILA') %>%
  group_by(plot_id) %>%
  summarise(tph.init = sum(tpa_unadj.init/0.404686)) %>%
  ungroup()  %>%
  pull(tph.init)

round(recr_tpa/0.404686,2)==round(comparison_tph,2)



recr_list$a


lapply(X = colnames(recr_list$X),
       FUN = function(n){
         recr_list$X %>%
           as.data.frame() %>%
           ggplot(aes(x = .data[[n]]))+
           geom_histogram()
       })
