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

names(growth_data)
# plot X distributions
lapply(X = c('dbh_m.init','ba_scaled',
             'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = .data[[v]]))+
           geom_histogram()
       })
lapply(X = c('fire', 'insects', 'disease'),
       FUN = function(v){
  ggplot(data = growth_data,
         aes(x = .data[[v]]))+
    geom_bar()
})

# too few instances of insects = 1

# plot Y distribution
ggplot(data = growth_data,
       aes(x = size1_g))+
  geom_histogram()

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
         lapply(X = c('fire', 'disease'),
                FUN = function(v_j){
                  ggplot(data = growth_data,
                         aes(x = as.factor(.data[[v_j]]), y = .data[[v_i]]))+
                    geom_violin()+
                    geom_jitter(width = 0.1, size = 0, height = 0)+
                    geom_boxplot(width = 0.1, position= position_nudge(x = -0.25))
                })
       })
ggplot(data = growth_data,
       aes(x = as.factor(fire), fill = as.factor(disease)))+
  geom_bar(position = position_fill())
# fire is less common on disease plots; this should be fine I have a lot of plots with either 
# disease or fire, and a few with both


# plot XY distributions
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = growth_data,
                aes(x = .data[[v]], y = size1_g))+
           geom_point()
       })
lapply(X = c('fire', 'disease'),
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
lapply(X = c('fire', 'disease'),
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
lapply(X = c('fire', 'insects', 'disease'),
       FUN = function(v){
  ggplot(data = surv_data,
         aes(x = .data[[v]]))+
    geom_bar()
})

# too few instances of insects = 1

# plot Y distribution
ggplot(data = surv_data,
       aes(x = as.factor(surv)))+
  geom_bar()

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
         lapply(X = c('fire', 'disease'),
                FUN = function(v_j){
                  ggplot(data = surv_data,
                         aes(x = as.factor(.data[[v_j]]), y = .data[[v_i]]))+
                    geom_violin()+
                    geom_jitter(width = 0.1, size = 0, height = 0)+
                    geom_boxplot(width = 0.1, position= position_nudge(x = -0.25))
                })
       })
ggplot(data = surv_data,
       aes(x = as.factor(fire), fill = as.factor(disease)))+
  geom_bar(position = position_fill())
# fire is less common on disease plots; this should be fine I have a lot of plots with either 
# disease or fire, and a few with both


# plot XY distributions
lapply(X = c('dbh_m.init', 'ba_scaled', 'cwd_dep90_scaled', 'cwd_mean_scaled'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(y = .data[[v]], x = as.factor(surv)))+
           geom_violin()+
           geom_jitter(width = 0.1, size = 0, height = 0)+
           geom_boxplot(width = 0.1, position = position_nudge(x = -0.25))
       })
lapply(X = c('fire', 'disease'),
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
                               breaks = quantile(.data[[v]], probs = c(0, 0.33, 0.66, 1)),
                               include_lowest = TRUE)))+
           geom_jitter(width = 0, height = 0.1, size = 0)+
           geom_smooth(method = 'glm', method.args = list(family = 'binomial'))+
           theme(legend.position = 'bottom')
       })
lapply(X = c('fire', 'disease'),
       FUN = function(v){
         ggplot(data = surv_data,
                aes(x = dbh_m.init,
                    y = surv,
                    color = as.factor(.data[[v]])))+
           geom_jitter(width = 0, height = 0.1, size = 0)+
           geom_smooth(method = 'glm', method.args = list(family = 'binomial'))+
           theme(legend.position = 'bottom')
       })

