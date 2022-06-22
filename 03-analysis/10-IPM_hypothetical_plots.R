library(here)
library(tidyverse)
library(posterior)
library(bayesplot)


#### setup #####################################################################

hypothetical_plots = 
  data.frame(intercept = rep(1, times = 9),
             fire = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
             wpbr = c(FALSE, FALSE, T, F, F, F, F, F, F),
             ba_scaled = c(0, 0, 0, -1, 1, 0, 0, 0, 0),
             cwd_dep90_scaled = c(0, 0, 0, 0, 0, -1, 1, 0, 0),
             cwd_mean_scaled = c(0, 0, 0, 0, 0, 0, 0, -1, 1),
             subp_id = 1:9,
             name = c('Undisturbed', 'Fire', 'WPBR', 'Low BA', 'High BA',
                      'Low Drought', 'High Drought', 'Wet Site', 'Dry Site'))

# load mcmc results
posterior_s = readRDS(here::here('02-data',
                               '03-results',
                               'surv_post.rds'))

posterior_g = readRDS(here::here('02-data',
                               '03-results',
                               'growth_post.rds'))

posterior_f = readRDS(here::here('02-data',
                               '03-results',
                               'fecd_post.rds'))




size_metadata = 
  readRDS(here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.rds')) %>%
  # convert to metric
  mutate(bin_midpoint = bin_midpoint * 0.0254,
         bin_lower = bin_lower * 0.0254,
         bin_upper = bin_upper * 0.0254) 


#### transitions ###############################################################


A_hypotheticals = 
  array(dim = c(nrow(size_metadata), # sizeclass to
                nrow(size_metadata), # sizeclass from
                nrow(hypothetical_plots), # plots
                nrow(posterior_s)), # posterior draws
        dimnames = list('class_to' = 1:nrow(size_metadata),
                        'class_from' = 1:nrow(size_metadata),
                        'plot' = 1:nrow(hypothetical_plots),
                        'draw' = 1:nrow(posterior_s)),
        data = 
          sapply(X = 1:nrow(posterior_s),
                 FUN = function(draw){
                   
                   # get beta_s for the current draw
                   beta_s = 
                     posterior_s %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   beta_g = 
                     posterior_g %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   beta_f = 
                     posterior_f %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   sigmaEpsilon_g = 
                     posterior_g %>%
                     slice(draw) %>%
                     pull(sigma_epsilon) %>%
                     as.numeric()
                   
                   sapply(X = 1:nrow(hypothetical_plots),
                          FUN = function(plot){
                            
                            # construct explanatory variable matrix for survival 
                            # for the current plot
                            X_g = 
                              hypothetical_plots %>%
                              slice(plot) %>%
                              expand(nesting(intercept, fire, wpbr, ba_scaled,
                                             cwd_dep90_scaled,cwd_mean_scaled),
                                     dbh = size_metadata$bin_midpoint) %>%
                              mutate(dbh2 = dbh**2,
                                     dbh_fire = dbh*fire,
                                     dbh2_fire = dbh2*fire,
                                     dbh_wpbr = dbh*wpbr,
                                     dbh2_wpbr = dbh2*wpbr,
                                     dbh_ba = dbh*ba_scaled,
                                     dbh2_ba = dbh2*ba_scaled,
                                     dbh_cwd_dep90 = dbh*cwd_dep90_scaled,
                                     dbh2_cwd_dep90 = dbh2*cwd_dep90_scaled,
                                     dbh_cwd_mean = dbh*cwd_mean_scaled,
                                     dbh2_cwd_mean = dbh2*cwd_mean_scaled) %>%
                              select(intercept, dbh, dbh2, fire, wpbr, ba_scaled,
                                     cwd_dep90_scaled, cwd_mean_scaled, 
                                     dbh_fire, dbh2_fire, dbh_wpbr, dbh2_wpbr, 
                                     dbh_ba, dbh2_ba,
                                     dbh_cwd_dep90, dbh2_cwd_dep90,
                                     dbh_cwd_mean, dbh2_cwd_mean) %>%
                              as.matrix()
                            
                            
                            X_s = 
                              hypothetical_plots %>%
                              slice(plot) %>%
                              expand(nesting(intercept, fire, wpbr, ba_scaled,
                                             cwd_dep90_scaled,cwd_mean_scaled),
                                     dbh = size_metadata$bin_midpoint) %>%
                              mutate(dbh2 = dbh**2,
                                     dbh_fire = dbh*fire,
                                     dbh2_fire = dbh2*fire,
                                     dbh_wpbr = dbh*wpbr,
                                     dbh2_wpbr = dbh2*wpbr,
                                     dbh_ba = dbh*ba_scaled,
                                     dbh2_ba = dbh2*ba_scaled,
                                     dbh_cwd_dep90 = dbh*cwd_dep90_scaled,
                                     dbh2_cwd_dep90 = dbh2*cwd_dep90_scaled,
                                     dbh_cwd_mean = dbh*cwd_mean_scaled,
                                     dbh2_cwd_mean = dbh2*cwd_mean_scaled) %>%
                              select(intercept, dbh, dbh2, fire, wpbr, ba_scaled,
                                     cwd_dep90_scaled, cwd_mean_scaled, 
                                     dbh_fire, dbh2_fire, dbh_wpbr, dbh2_wpbr, 
                                     dbh_ba, dbh2_ba,
                                     dbh_cwd_dep90, dbh2_cwd_dep90,
                                     dbh_cwd_mean, dbh2_cwd_mean) %>%
                              as.matrix()
                            
                            X_f = 
                              hypothetical_plots %>%
                              slice(plot) %>%
                              expand(nesting(intercept, fire, wpbr, ba_scaled,
                                             cwd_dep90_scaled,cwd_mean_scaled),
                                     dbh = size_metadata$bin_midpoint) %>%
                              mutate(dbh2 = dbh**2,
                                     dbh_fire = dbh*fire,
                                     dbh2_fire = dbh2*fire,
                                     dbh_wpbr = dbh*wpbr,
                                     dbh2_wpbr = dbh2*wpbr,
                                     dbh_ba = dbh*ba_scaled,
                                     dbh2_ba = dbh2*ba_scaled,
                                     dbh_cwd_dep90 = dbh*cwd_dep90_scaled,
                                     dbh2_cwd_dep90 = dbh2*cwd_dep90_scaled,
                                     dbh_cwd_mean = dbh*cwd_mean_scaled,
                                     dbh2_cwd_mean = dbh2*cwd_mean_scaled) %>%
                              select(intercept) %>%
                              as.matrix()
                            
                            
                            # calculate size_from length vector of survival 
                            # probabilities on this plot with this parameter draw
                            p = boot::inv.logit(as.numeric(X_s %*% beta_s))
                            
                            mu = as.numeric(X_g %*% beta_g)
                            
                            f = exp(as.numeric(X_f %*% beta_f))
                            
                            sapply(X = 1:nrow(size_metadata),
                                   FUN = function(class_from){
                                     
                                     g = 
                                       ((pnorm(size_metadata$bin_upper,
                                               mu[class_from],
                                               sigmaEpsilon_g) - 
                                           pnorm(size_metadata$bin_lower,
                                                 mu[class_from],
                                                 sigmaEpsilon_g))/
                                          (1-pnorm(0.0254,
                                                mu[class_from],
                                                sigmaEpsilon_g)))
                                     
                                     transition_prob = 
                                       # survival of each from class
                                       (p[class_from] *
                                          # prob of growth from to
                                          g) +
                                       # number of new recruits
                                       (f[class_from] *
                                          size_metadata$r)
                                     
                                     return(transition_prob)
                                   })
                            })
                 }))


saveRDS(A_hypotheticals,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'hypothetical_As.rds'))

#### lambda ####################################################################

A_hypotheticals = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'hypothetical_As.rds'))

hypothetical_lambdas_full = 
  hypothetical_plots %>%
  expand(nesting(subp_id, name, intercept, fire, wpbr, ba_scaled, 
                 cwd_dep90_scaled, cwd_mean_scaled),
         data.frame(draw = 1:4000))


hypothetical_lambdas_full$lambda = 
  sapply(X = 1:nrow(hypothetical_plots),
         FUN = function(plot){
           sapply(X = 1:4000,
                  FUN = function(draw){
                    # paste0('s:',plot,'d:',draw) for testing
                    max(as.numeric(Re(eigen(A_hypotheticals[,,plot,draw])$values)))
                  })
         }) %>%
  as.numeric()


pretty_names = 
  hypothetical_plots$name
names(pretty_names) = hypothetical_plots$subp_id

hypothetical_lambdas_plot = 
  ggplot(data = 
           hypothetical_lambdas_full,
         aes(x = lambda))+
  geom_density(lwd = 1, fill = 'lightgrey')+
  geom_vline(xintercept = 1, color = 'red', lty = 2, lwd = 1)+
  theme_minimal()+
  facet_grid(subp_id~., scales = 'free_y',
             labeller = labeller(subp_id = pretty_names))+
  theme(axis.text.y = element_blank())+
  labs(x = 'Asymptotic Population Growth Rate', y = 'Posterior Density')

hypothetical_lambdas_plot

ggsave(hypothetical_lambdas_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'hypotheticals_lambda_post.png'),
       height = 7.5, width = 4, units = 'in')

hypothetical_lambdas_plot.ppt = 
  ggplot(data = 
           hypothetical_lambdas_full,
         aes(x = lambda))+
  geom_density( fill = 'lightgrey')+
  geom_vline(xintercept = 1, color = 'red', lty = 2, lwd = 1)+
  theme_minimal()+
  facet_grid(.~subp_id, scales = 'free_x',
             labeller = labeller(subp_id = pretty_names))+
  theme(axis.text.x = element_blank(),
        text= element_text(size = 16))+
  labs(x = 'Asymptotic Population Growth Rate', y = 'Posterior Density')+
  coord_flip()

hypothetical_lambdas_plot.ppt


ggsave(hypothetical_lambdas_plot.ppt,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'hypotheticals_lambda_post.png'),
       height = 4, width = 11, units = 'in')

hypothetical_lambda_summary = 
  hypothetical_lambdas_full %>%
  group_by(subp_id, name) %>%
  summarise(lambda.med = median(lambda),
            lambda.05 = quantile(lambda, probs = 0.05),
            lambda.95 = quantile(lambda, probs = 0.95))

hypothetical_lambda_summary

write.csv(hypothetical_lambda_summary,
          here::here('04-communication',
                     'tables',
                     'hypothetical_lambdas_summary.csv'))

#### ssd #######################################################################

hypothetical_ssd = 
  array(dim = c(nrow(size_metadata),
                nrow(hypothetical_plots),
                4000),
        dimnames = 
          list('sizeclass' = 1:nrow(size_metadata),
               'plot' = 1:nrow(hypothetical_plots),
               'draw' = 1:4000),
        data = 
          sapply(X = 1:4000,
                 FUN = function(draw){
                   sapply(X = 1:nrow(hypothetical_plots),
                          FUN = function(plot){
                            A = A_hypotheticals[,,plot,draw]
                            w.eigen = Re(eigen(A)$vectors[,1])
                            ssd = w.eigen / sum(w.eigen)
                            return(ssd)
                          })
                 }))


hypothetical_ssd.df = 
  hypothetical_plots %>%
  expand(nesting(intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled,
                 cwd_mean_scaled, subp_id, name),
         sizeclass = 1:nrow(size_metadata),
         draw = 1:4000) %>%
  left_join(size_metadata %>%
              select(sizeclass = bin_id, bin_midpoint)) %>%
  arrange(subp_id, sizeclass, draw)

hypothetical_ssd.df$class_proportion = 
  
  as.numeric(
    sapply(X = 1:nrow(hypothetical_plots),
           FUN = function(plot){
             
             as.numeric(
               sapply(X = 1:nrow(size_metadata),
                      FUN = function(sizeclass){
                        
                        return(hypothetical_ssd[sizeclass,plot,])
                      })
             )
             
           }))

# Fire makes the SSD really for the very rare bigger size classes, fewer 
# mid-to-large trees and more superlarge 
hypothetical_ssd.df %>%
  group_by(bin_midpoint, sizeclass, name, subp_id) %>%
  summarise(prop.med = median(class_proportion),
            prop.05 = quantile(class_proportion, 0.05),
            prop.95 = quantile(class_proportion, 0.95)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint, color = name, fill = name))+
  #geom_point(aes(y = prop.med), size = 1)+
  geom_line(aes(y = prop.med))+
  geom_ribbon(aes(ymin = prop.05, ymax = prop.95), alpha = 0.25)+
  theme_minimal()+
  scale_y_log10()+
  facet_wrap(~name)

summary(hypothetical_ssd.df)

#### reproductive value ########################################################


hypothetical_repro = 
  array(dim = c(nrow(size_metadata),
                nrow(hypothetical_plots),
                4000),
        dimnames = 
          list('sizeclass' = 1:nrow(size_metadata),
               'plot' = 1:nrow(hypothetical_plots),
               'draw' = 1:4000),
        data = 
          sapply(X = 1:4000,
                 FUN = function(draw){
                   sapply(X = 1:nrow(hypothetical_plots),
                          FUN = function(plot){
                            A = A_hypotheticals[,,plot,draw]
                            v.eigen = Re(eigen(t(A))$vectors[,1])
                            rv = v.eigen / v.eigen[1]
                            return(rv)
                          })
                 }))


hypothetical_repro.df = 
  hypothetical_plots %>%
  expand(nesting(intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled,
                 cwd_mean_scaled, subp_id, name),
         sizeclass = 1:nrow(size_metadata),
         draw = 1:4000) %>%
  left_join(size_metadata %>%
              select(sizeclass = bin_id, bin_midpoint)) %>%
  arrange(subp_id, sizeclass, draw)

hypothetical_repro.df$repro = 
  
  as.numeric(
    sapply(X = 1:nrow(hypothetical_plots),
           FUN = function(plot){
             
             as.numeric(
               sapply(X = 1:nrow(size_metadata),
                      FUN = function(sizeclass){
                        
                        return(hypothetical_repro[sizeclass,plot,])
                      })
             )
             
           }))

# Again these reproductive values are wonky, esp for burned plots. Something 
# about the transition matrix for burned plots is weird.
hypothetical_repro.df %>%
  group_by(bin_midpoint, sizeclass, name, subp_id) %>%
  #filter(subp_id != 2) %>%
  summarise(repro.med = median(repro),
            repro.05 = quantile(repro, 0.05, na.rm = TRUE),
            repro.95 = quantile(repro, 0.95, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint, color = name, fill = name))+
  #geom_point(aes(y = prop.med), size = 1)+
  geom_line(aes(y = repro.med))+
  geom_ribbon(aes(ymin = repro.05, ymax = repro.95), alpha = 0.25)+
  theme_minimal()+
  facet_wrap(~name, scales = 'free_y')


summary(hypothetical_repro.df)

#### senstivity and elasticity #################################################

hypothetical_lambdas.matrix = 
  matrix(nrow = nrow(hypothetical_plots),
         ncol = 4000,
         data = 
           sapply(X = 1:nrow(hypothetical_plots),
                  FUN = function(plot){
                    sapply(X = 1:4000,
                           FUN = function(draw){
                             max(Re(eigen(A_hypotheticals[,,plot,draw])$values))
                           })
                  }))

# sensitivity and elasticity
hypothetical_vdotw = 
  matrix(nrow = nrow(hypothetical_plots),
         ncol = 4000,
         data = 
           sapply(X = 1:nrow(hypothetical_plots),
                  FUN = function(plot){
                    sapply(X = 1:4000,
                           FUN = function(draw){
                             sum(hypothetical_ssd[,plot,draw] *
                                   hypothetical_repro[,plot,draw])*0.127
                           })
                  }))


hypothetical_sens = 
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(hypothetical_plots),
                   1000),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'plot' = 1:nrow(hypothetical_plots),
               'draw' = 1:1000),
        data = 
          
          sapply(X = 1:1000,
                 FUN = function(draw){
                   sapply(X = 1:nrow(hypothetical_plots),
                          FUN = function(plot){
                            outer(hypothetical_repro[,plot,draw], 
                                  hypothetical_ssd[,plot,draw])/
                              hypothetical_vdotw[plot,draw]
                          })
                 })
  )

hypothetical_elas = 
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(hypothetical_plots),
                   1000),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'plot' = 1:nrow(hypothetical_plots),
               'draw' = 1:1000),
        data = 
          
          sapply(X = 1:1000,
                 FUN = function(draw){
                   sapply(X = 1:nrow(hypothetical_plots),
                          FUN = function(plot){
                            
                            matrix(as.vector(hypothetical_sens[,,plot,draw])*
                                     as.vector(A_hypotheticals[,,plot,draw])/
                                     hypothetical_lambdas.matrix[plot,draw])
                            
                            
                          })
                 })
  )



hypothetical_sens_elas.df  = 
  expand.grid(size_to = size_metadata$bin_id,
              size_from = size_metadata$bin_id) %>%
  left_join(size_metadata %>%
              select(size_to = bin_id, 
                     bin_midpoint_to_m = bin_midpoint)) %>%
  left_join(size_metadata %>%
              select(size_from = bin_id,
                     bin_midpoint_from_m = bin_midpoint)) %>%
  expand(nesting(size_to, size_from, bin_midpoint_to_m, bin_midpoint_from_m),
         subp_id = 1:9,
         draw = 1:1000) %>%
  arrange(subp_id, size_to, size_from, draw)


hypothetical_sens_elas.df$sensitivity = 
  sapply(X = 1:nrow(hypothetical_plots),
         FUN = function(plot){
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    
                    sapply(X = 1:nrow(size_metadata),
                           FUN = function(size_from){
                             hypothetical_sens[size_to, size_from, plot,]
                           })
                  })
         }) %>%
  as.vector()


hypothetical_sens_elas.df$elasticity = 
  sapply(X = 1:nrow(hypothetical_plots),
         FUN = function(plot){
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    
                    sapply(X = 1:nrow(size_metadata),
                           FUN = function(size_from){
                             hypothetical_elas[size_to, size_from, plot,]
                           })
                  })
         }) %>%
  as.vector()


sensitivity_plots = 
  lapply(X = 1:9,
         FUN = function(s){
                       
              hypothetical_sens_elas.df %>%
              group_by(size_to, size_from, bin_midpoint_to_m, bin_midpoint_from_m,
                       subp_id) %>%
              summarise(sensitivity = median(sensitivity, na.rm = TRUE),
                        elasticity = median(elasticity, na.rm = TRUE)) %>%
              ungroup() %>%
              left_join(hypothetical_plots %>%
                          select(subp_id, name)) %>%
              filter(subp_id==s) %>%
              ggplot(aes(x = bin_midpoint_from_m,
                         y = bin_midpoint_to_m,
                         fill = sensitivity))+
              geom_tile()+
              coord_fixed()+
              theme_minimal()+
              scale_fill_viridis_c()+
                labs(title = hypothetical_plots$name[s])
         })
  
sensitivity_plots  

elasticity_plots = 
  lapply(X = 1:9,
         FUN = function(s){
                       
              hypothetical_sens_elas.df %>%
              group_by(size_to, size_from, bin_midpoint_to_m, bin_midpoint_from_m,
                       subp_id) %>%
              summarise(sensitivity = median(sensitivity, na.rm = TRUE),
                        elasticity = median(elasticity, na.rm = TRUE)) %>%
              ungroup() %>%
              left_join(hypothetical_plots %>%
                          select(subp_id, name)) %>%
              filter(subp_id==s) %>%
              ggplot(aes(x = bin_midpoint_from_m,
                         y = bin_midpoint_to_m,
                         fill = elasticity))+
              geom_tile()+
              coord_fixed()+
              theme_minimal()+
              scale_fill_viridis_c()+
                labs(title = hypothetical_plots$name[s])
         })
  
elasticity_plots  


