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
pila_data = readRDS(here::here('02-data',
                               '02-for_analysis',
                               'pila_training.rds'))

#### growth data ###############################################################

growth_data = 
  pila_data$X_g %>%
  as_tibble()

growth_data$size1_g = pila_data$size1_g

growth_data %>%
ggplot(aes(x = dbh_m.init, y = size1_g-dbh_m.init))+
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x + I(x**2))

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

# too few instances of insects = 1

# plot Y distribution
response_growth = 
  ggplot(data = growth_data,
       aes(x = size1_g))+
  geom_histogram()+
  theme_minimal()+
  theme(text = element_text(size = 14))+
  labs(y = 'N trees', x = 'DBH at remeasurement (m)')

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
                aes(x = .data[[v]], y = size1_g))+
           geom_point()
       })
lapply(X = c('fire', 'wpbr'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = as.factor(.data[[v]]), y = size1_g))+
                    geom_violin()+
                    geom_jitter(width = 0.1, size = 0, height = 0)+
                    geom_boxplot(width = 0.1, position= position_nudge(x = -0.25))
       })

# plot XXY distributions for size with other vars
lapply(X = c('ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = dbh_m.init, 
                    y = size1_g, 
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
                    y = size1_g,
                    color = as.factor(.data[[v]])))+
           geom_point()+
           geom_smooth(method = 'lm')+
           theme(legend.position = 'bottom')
       })


#### survival data #############################################################

surv_data = 
  pila_data$X_s %>%
  as_tibble()

surv_data$surv = pila_data$surv

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

lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(x = .data[[v]], y = surv))+
           geom_jitter(width = 0, size = 0, height = 0.1)+
           geom_smooth(method = 'lm', formula = y~x+I(x**2))
       })

surv_data$ecosub.i = pila_data$ecosub_s

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
pila_data$cprime

hist(pila_data$cprime)


# check total TPA on each subplot looks OK
recr_tpa = 
  sapply(X = 1:nrow(pila_data$n),
         FUN = function(subplot){
           subplot_tpa = sum(pila_data$n[subplot,])
           return(subplot_tpa)
         })

hist(recr_tpa)
summary(recr_tpa/0.404686)

readRDS(here::here('02-data',
           '01-preprocessed',
           'sizedist_data.rds')) %>%
  filter(
    is.element(plot_id,
                 readRDS(here::here('02-data',
                   '02-for_analysis',
                   'union_plots.rds')) %>%
                 filter(is.element(plot_id.i,
                                   pila_data$plotid_r)) %>%
                 pull(plot_id))
    ) %>%
  filter(species=='PILA') %>%
  group_by(plot_id) %>%
  summarise(tph.init = sum(tpa_unadj.init/0.404686)) %>%
  ungroup() %>%
  summary()





# check TPA of bigger trees on each subplots
sapply(X = 1:nrow(pila_data$n),
       FUN = function(subplot){
         sum(pila_data$n[subplot,2:100])
       }) %>%
  summary()

pila_data$r

pila_data$a

pila_data$M_r

lapply(X = colnames(pila_data$X_r),
       FUN = function(n){
         pila_data$X_r %>%
           as.data.frame() %>%
           ggplot(aes(x = .data[[n]]))+
           geom_histogram()
       })

#### tph and BA ################################################################

tph_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'tph.rds'))

ba_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'ba.rds'))

ba_tph_df = 
  tph_data$X %>%
  as_tibble() %>%
  mutate(delta_tph = tph_data$Y,
         delta_ba = ba_data$Y,
         ecosub_id = tph_data$ecosub_id)

# y distribution
ggplot(data = ba_tph_df,
       aes(x = delta_tph))+
  geom_histogram()

ggplot(data = ba_tph_df,
       aes(x = delta_tph))+
  geom_histogram()+
  scale_x_continuous(limits = c(-250, 250))

ggplot(data = ba_tph_df,
       aes(x = delta_ba))+
  geom_histogram()

# x distributions
ggplot(data = ba_tph_df,
       aes(x = fire))+
  geom_bar()

ggplot(data = ba_tph_df,
       aes(x = wpbr))+
  geom_bar()

ggplot(data = ba_tph_df,
       aes(x = ba_scaled))+
  geom_histogram()

ggplot(data = ba_tph_df,
       aes(x = cwd_dep90_scaled))+
  geom_histogram()

ggplot(data = ba_tph_df,
       aes(x = cwd_mean_scaled))+
  geom_histogram()

# xx relationships
ggplot(data = ba_tph_df,
       aes(x = ba_scaled, color = as.factor(fire)))+
  geom_density()

ggplot(data = ba_tph_df,
       aes(x = cwd_dep90_scaled, color = as.factor(fire)))+
  geom_density()
ggplot(data = ba_tph_df,
       aes(x = cwd_mean_scaled, color = as.factor(fire)))+
  geom_density()

ggplot(data = ba_tph_df,
       aes(x = ba_scaled, color = as.factor(wpbr)))+
  geom_density()

ggplot(data = ba_tph_df,
       aes(x = cwd_dep90_scaled, color = as.factor(wpbr)))+
  geom_density()
ggplot(data = ba_tph_df,
       aes(x = cwd_mean_scaled, color = as.factor(wpbr)))+
  geom_density()

ggplot(data = ba_tph_df,
       aes(x = ba_scaled, y = cwd_dep90_scaled))+
  geom_point()

ggplot(data = ba_tph_df,
       aes(x = ba_scaled, y = cwd_mean_scaled))+
  geom_point()

ggplot(data = ba_tph_df,
       aes(x = cwd_dep90_scaled, y = cwd_mean_scaled))+
  geom_point()

# xy relationships
ggplot(data = ba_tph_df,
       aes(x = as.factor(fire), y = delta_tph))+
  geom_boxplot()+
  scale_y_continuous(limits = c(-250, 250))

ggplot(data = ba_tph_df,
       aes(x = as.factor(wpbr), y = delta_tph))+
  geom_boxplot()+
  scale_y_continuous(limits = c(-250, 250))

ggplot(data = ba_tph_df,
       aes(x = ba_scaled, y = delta_tph))+
  geom_point()+
  geom_smooth(method = 'lm')

ggplot(data = ba_tph_df,
       aes(x = cwd_dep90_scaled, y = delta_tph))+
  geom_point()+
  geom_smooth(method = 'lm')

ggplot(data = ba_tph_df,
       aes(x = cwd_mean_scaled, y = delta_tph))+
  geom_point()+
  geom_smooth(method = 'lm')

ggplot(data = ba_tph_df,
       aes(x = as.factor(fire), y = delta_ba))+
  geom_boxplot()

ggplot(data = ba_tph_df,
       aes(x = as.factor(wpbr), y = delta_ba))+
  geom_boxplot()

ggplot(data = ba_tph_df,
       aes(x = ba_scaled, y = delta_ba))+
  geom_point()+
  geom_smooth(method = 'lm')

ggplot(data = ba_tph_df,
       aes(x = cwd_dep90_scaled, y = delta_ba))+
  geom_point()+
  geom_smooth(method = 'lm')

ggplot(data = ba_tph_df,
       aes(x = cwd_mean_scaled, y = delta_ba))+
  geom_point()+
  geom_smooth(method = 'lm')



