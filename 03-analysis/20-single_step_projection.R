
#### setup #####################################################################

library(here)
library(tidyverse)


#### get initial size distribution #############################################

tph_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'tph_data.rds'))


sizedist_data = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'sizedist_data.rds'))

tph_vector = 

  sizedist_data %>%
      filter(species=='PILA') %>%
      group_by(dbh_class) %>%
      summarise(mean_tpa.init = mean(tpa_unadj.init),
                se_tpa.init = sd(tpa_unadj.init)/sqrt(n())) %>%
      ungroup() %>%
      mutate(source = 'subset') %>%
  arrange(dbh_class) %>%
  mutate(tph = mean_tpa.init / 0.404686) %>%
  pull(tph)


#### load transition matrices ##################################################

A_hypotheticals = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'hypothetical_As.rds'))

A_observed = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'A_observed.rds'))

#### multiply transition matrices by initial size distribution #################

nprime_hypotheticals = 
  array(dim = list(dim(A_hypotheticals)[1],
                   dim(A_hypotheticals)[3],
                   dim(A_hypotheticals)[4]),
        dimnames = list(dbh_class = 1:99,
                        plot = 1:9,
                        draw = 1:4000),
        data = 
          sapply(X = 1:4000,
                 FUN = function(draw){
                   
                   sapply(X = 1:9,
                          FUN = function(p){
                            
                            post_distribution = 
                              A_hypotheticals[,,p,draw] %*%
                              matrix(nrow = 99, ncol = 1, data = tph_vector)
                            
                            return(post_distribution)
                            
                          })
                 }))



#### plot hypotheticals post-size distribution #################################

hypothetical_results_df = 
  expand.grid(plot = 1:9,
              draw = 1:4000,
              dbh_class = 1:99) %>%
  left_join(
    data.frame(dbh_class = 1:99,
               tph_initial = tph_vector)
  ) %>%
  arrange(plot, draw, dbh_class)

hypothetical_results_df$tph_post = 
  sapply(X = 1:9,
         FUN = function(p){
           sapply(X = 1:4000,
                  FUN = function(draw){
                    nprime_hypotheticals[,p,draw]
                  })
         }) %>%
  as.vector()


tph_figure = 
  hypothetical_results_df %>%
  group_by(plot, draw) %>%
  summarise(tph_initial = sum(tph_initial),
            tph_post = sum(tph_post),
            tph_delta = sum(tph_post-tph_initial)) %>%
  ungroup() %>%
  left_join(
    data.frame(plot = 1:9,
             name = c('Undisturbed', 'Fire', 'WPBR', 'Low BA', 'High BA',
                      'Low Drought', 'High Drought', 'Wet Site', 'Dry Site'))
  ) %>%  
  mutate(name = factor(name, levels = 
                         c('Undisturbed', 'Fire', 'WPBR', 'Low BA', 'High BA',
                           'Low Drought', 'High Drought', 'Wet Site', 'Dry Site'))) %>%
  ggplot(aes(x = tph_post))+
  geom_density(fill = 'lightgrey')+
  facet_grid(name~.)+
  geom_vline(aes(xintercept = tph_initial),
             lty = 2, lwd = 1, col = 'red')+
  theme_minimal()+
  theme(axis.text.y = element_blank())+
  labs(y = 'Posterior Density',
       x = 'Stem density (trees / ha)')+
  scale_x_continuous(limits = c(0, 75))

tph_figure

ggsave(tph_figure,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'projected_tph.png'),
       height = 7.5, width = 4, units = 'in')

ba_figure = 
  hypothetical_results_df %>%
  left_join(
    data.frame(dbh_class = 1:100,
               dbh_m = seq(from = 0.0127, to = 2.5273, by = 0.0254)) %>%
      mutate(ba_m2 = pi*((dbh_m/2)**2))
  ) %>%
  mutate(ba_m2ha.initial = ba_m2*tph_initial,
         ba_m2ha.post = ba_m2*tph_post) %>%
  group_by(plot, draw) %>%
  summarise(tph_initial = sum(tph_initial),
            tph_post = sum(tph_post),
            ba_initial = sum(ba_m2ha.initial),
            ba_post = sum(ba_m2ha.post)) %>%
  ungroup() %>%
  left_join(
    data.frame(plot = 1:9,
             name = c('Undisturbed', 'Fire', 'WPBR', 'Low BA', 'High BA',
                      'Low Drought', 'High Drought', 'Wet Site', 'Dry Site'))
  ) %>%  
  mutate(name = factor(name, levels = 
                         c('Undisturbed', 'Fire', 'WPBR', 'Low BA', 'High BA',
                           'Low Drought', 'High Drought', 'Wet Site', 'Dry Site'))) %>%
  
  ggplot(aes(x = ba_post))+
  geom_density(fill = 'lightgrey')+
  facet_grid(name~.)+
  geom_vline(aes(xintercept = ba_initial),
             lty = 2, lwd = 1, col = 'red')+
  theme_minimal()+
  theme(axis.text.y = element_blank())+
  labs(y = 'Posterior Density',
       x = 'Basal Area (m^2 / ha)')+
  scale_x_continuous(limits = c(0, 7))

ba_figure

ggsave(ba_figure,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'projected_ba.png'),
       height = 7.5, width = 4, units = 'in')

#### checking random effects distribution ######################################

posterior = readRDS(here::here('02-data', 
                               '03-results',
                               'real_fits',
                               'posterior_draws.rds'))

posterior %>%
  select(contains('ecoEffect_f')) %>%
  summarise_all(median) %>%
  as.numeric() %>%
  hist()


#### multiply transition matrix by initial size distribution: observed plots ###

# start with the size distribution data
sizedist_to_project = 
  sizedist_data %>%
  
  # keep only the PILA
  filter(species=='PILA') %>%
  
  # inner join the plot indices
  inner_join(
    readRDS(here::here('02-data', 
                       '02-for_analysis',
                       'union_plots.rds'))
  ) %>%
  arrange(plot_id.i, dbh_class)
  

sizedist_to_project$nprime = 
  sapply(X = unique(sizedist_to_project$plot_id.i),
         FUN = function(p){
           post_distribution = 
             A_observed[,,p] %*%
             matrix(nrow = 99, ncol = 1,
                    data = sizedist_to_project[sizedist_to_project$plot_id.i==p,]$tpa_unadj.init)
           
           return(post_distribution)
         }) %>%
  as.vector()





#### plot observed plots post distribution #####################################


sizedist_to_project %>%
  group_by(plot_id) %>%
  summarise(tph_init = sum(tpa_unadj.init / 0.404686),
            tph_re = sum(tpa_unadj.re / 0.404686),
            tph_predicted = sum(nprime / 0.404686)) %>%
  ungroup() %>%
  summarise(tph_init.mean = mean(tph_init),
            tph_init.se = sd(tph_init)/sqrt(n()),
            tph_re.mean = mean(tph_re),
            tph_re.se = sd(tph_re)/sqrt(n()),
            tph_predicted.mean = mean(tph_predicted),
            tph_predicted.se = sd(tph_predicted)/sqrt(n()))

sizedist_to_project %>%
  group_by(plot_id) %>%
  summarise(tph_init = sum(tpa_unadj.init / 0.404686),
            tph_re = sum(tpa_unadj.re / 0.404686),
            tph_predicted = sum(nprime / 0.404686)) %>%
  ungroup() %>%
  ggplot(aes(x = tph_predicted, y = tph_re))+
  geom_point(size = 0, alpha = 0.5)+
  geom_abline(intercept = 0, slope = 1, color = 'red')+
  geom_smooth(method = 'lm')+
  theme_minimal()
  #coord_cartesian(xlim = c(0, 2500), ylim = c(0, 2500))

sizedist_to_project = 
  sizedist_to_project %>%
  left_join(
    data.frame(dbh_class = 1:100,
               dbh_m = seq(from = 0.0127, to = 2.5273, by = 0.0254)) %>%
      mutate(ba_m2 = pi*((dbh_m/2)**2))
  ) %>%
  mutate(tph_init = tpa_unadj.init/0.404686,
         tph_re = tpa_unadj.re / 0.404686,
         tph_pred = nprime / 0.404686) %>%
  mutate(ba_m2ha_init = ba_m2*tph_init,
         ba_m2ha_re = ba_m2*tph_re,
         ba_m2ha_pred = ba_m2*tph_pred)

sizedist_to_project %>%
  group_by(plot_id) %>%
  summarise(ba_init = sum(ba_m2ha_init),
            ba_re = sum(ba_m2ha_re),
            ba_pred = sum(ba_m2ha_pred)) %>%
  ungroup() %>%
  summarise(ba_init.mean = mean(ba_init),
            ba_init.se = sd(ba_init)/sqrt(n()),
            ba_re.mean = mean(ba_re),
            ba_re.se = sd(ba_re)/sqrt(n()),
            ba_pred.mean = mean(ba_pred),
            ba_pred.se = sd(ba_pred)/sqrt(n()))

head(sizedist_to_project)

sizedist_to_project %>%
  group_by(dbh_class) %>%
  summarise(tph_init = mean(tph_init),
            tph_re = mean(tph_re),
            tph_pred = mean(tph_pred),
            ba_init = mean(ba_m2ha_init),
            ba_re = mean(ba_m2ha_re),
            ba_pred = mean(ba_m2ha_pred)) %>%
  ungroup() %>%
  pivot_longer(cols = c(-dbh_class),
               names_to = c('stat', 'timestep'),
               values_to = 'value',
               names_sep = '_') %>%
  filter(stat == 'tph') %>%
  ggplot(aes(x = dbh_class, y = value, color = timestep))+
  #geom_point()+
  theme_minimal()+
  geom_line(lwd = 1)


#### scratch ###################################################################

observed_results_df = 
  expand.grid(plot = unique(sizedist_to_project$plot_id.i),
              dbh_class = 1:100) %>%
  left_join(
    sizedist_to_project 
  ) %>%
  arrange(plot, dbh_class)

sizedist_to_project$nprime = 
  sapply(X = unique(sizedist_to_project$plot_id.i))


observed_results_df$tph_post = 
  sapply(X = 1:1102,
         FUN = function(p){
           nprime_observed[,p]
         }) %>%
  as.vector()


head(observed_results_df)

observed_results_df %>%
  group_by(plot) %>%
  summarise(tph_initial = sum(tph_initial),
            tph_post = sum(tph_post)) %>%
  ungroup() %>%
  ggplot(aes(x = tph_post))+
  geom_density(fill = 'lightgrey')+
  geom_vline(aes(xintercept = tph_initial),
             lty = 2, lwd = 1, col = 'red')+
  theme_minimal()+
  theme(axis.text.y = element_blank())+
  labs(y = 'Posterior Density',
       x = 'Stem density (trees / ha)')

observed_results_df %>%
  left_join(
    data.frame(dbh_class = 1:100,
               dbh_m = seq(from = 0.0127, to = 2.5273, by = 0.0254)) %>%
      mutate(ba_m2 = pi*((dbh_m/2)**2))
  ) %>%
  mutate(ba_m2ha.initial = ba_m2*tph_initial,
         ba_m2ha.post = ba_m2*tph_post) %>%
  group_by(plot) %>%
  summarise(tph_initial = sum(tph_initial),
            tph_post = sum(tph_post),
            ba_initial = sum(ba_m2ha.initial),
            ba_post = sum(ba_m2ha.post)) %>%
  ungroup() %>%
  pull(ba_post) %>%
  mean()


  
ggplot(aes(x = ba_post))+
  geom_density(fill = 'lightgrey')+
  geom_vline(aes(xintercept = ba_initial),
             lty = 2, lwd = 1, col = 'red')+
  theme_minimal()+
  theme(axis.text.y = element_blank())+
  labs(y = 'Posterior Density',
       x = 'Basal Area (m^2 / ha)')

ba_figure

ggsave(ba_figure,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'projected_ba.png'),
       height = 7.5, width = 4, units = 'in')

