
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


#### multiply transition matrices by initial size distribution #################

nprime_hypotheticals = 
  array(dim = list(dim(A_hypotheticals)[1],
                   dim(A_hypotheticals)[3],
                   dim(A_hypotheticals)[4]),
        dimnames = list(dbh_class = 1:100,
                        plot = 1:9,
                        draw = 1:4000),
        data = 
          sapply(X = 1:4000,
                 FUN = function(draw){
                   
                   sapply(X = 1:9,
                          FUN = function(p){
                            
                            post_distribution = 
                              A_hypotheticals[,,p,draw] %*%
                              matrix(nrow = 100, ncol = 1, data = tph_vector)
                            
                            return(post_distribution)
                            
                          })
                 }))



#### plot post-size distribution ###############################################

hypothetical_results_df = 
  expand.grid(plot = 1:9,
              draw = 1:4000,
              dbh_class = 1:100) %>%
  left_join(
    data.frame(dbh_class = 1:100,
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
       x = 'Stem density (trees / ha)')

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
       x = 'Basal Area (m^2 / ha)')

ba_figure

ggsave(ba_figure,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'projected_ba.png'),
       height = 7.5, width = 4, units = 'in')

