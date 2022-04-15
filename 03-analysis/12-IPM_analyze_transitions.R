



#### subplot tranasition matrices ############################################## 



subplot_transitions.med = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'subplot_As.rds'))

#### subplot lambda ############################################################
subplot_lambdas.med = 
  sapply(X = 1:nrow(subplots.pila),
         FUN = function(subplot){
           
           A.subplot = subplot_transitions.med[,,subplot]
           
           lambda.subplot = max(as.numeric(Re(eigen(A.subplot)$values)))
           return(lambda.subplot)
         })


# lambda: asymptotic population growth rate
subplots.pila$lambda_postmed = 
  subplot_lambdas.med

postmed_lambda_distribution = 
  ggplot(data = subplots.pila,
         aes(x = lambda_postmed))+
  geom_density()+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 2.5))+
  geom_vline(xintercept = 1, color = 'grey', lty = 2, lwd = 1)+
  labs(y = 'Density', x = 'Lambda')

postmed_lambda_distribution

postmed_lambda_distribution_ppt = 
  ggplot(data = subplots.pila,
         aes(x = lambda_postmed))+
  geom_histogram()+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 2.5))+
  geom_vline(xintercept = 1, color = 'grey', lty = 2, lwd = 1)+
  labs(y = 'N subplots', x = 'Lambda')+
  theme(text = element_text(size = 18))

ggsave(postmed_lambda_distribution_ppt,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'subplot_lambdas.png'),
       height = 4, width = 6.5, units = 'in')

summary(subplots.pila$lambda_postmed)


# what proportion of plots is the pop predicted to decline on
length(subplots.pila$lambda_postmed[subplots.pila$lambda_postmed<1])/
  length(subplots.pila$lambda_postmed)

ggsave(postmed_lambda_distribution,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'postmed_lambda_distribution.png'),
       height = 4, width = 6.5, units = 'in')

postmed_lambda_distribution

#### subplot stable size distribution ##########################################

# stable size distribution
subplot_ssd.med = 
  matrix(nrow = nrow(size_metadata),
         ncol = nrow(subplots.pila),
         byrow = FALSE,
         data = 
           sapply(X = 1:nrow(subplots.pila),
                  FUN = function(subplot){
                    A.subplot = subplot_transitions.med[,,subplot]
                    # from supplamentory materials for merow et al 2014 
                    # "On using integral projection models..."
                    w.eigen = Re(eigen(A.subplot)$vectors[,1])
                    ssd = w.eigen / sum(w.eigen)
                    return(ssd)
                  }))

subplot_ssd.df = 
  expand.grid('subplot' = 1:nrow(subplots.pila),
              'sizeclass' = 1:nrow(size_metadata)) %>%
  left_join(size_metadata %>%
              mutate(bin_midpoint_cm = bin_midpoint*100) %>%
              select(sizeclass = bin_id, bin_midpoint_cm))

subplot_ssd.df$class_proportion = 
  
  as.numeric(
    sapply(X = 1:nrow(size_metadata),
           FUN = function(sizeclass){
             
             return(subplot_ssd.med[sizeclass,])
             
           })
  )

# stable size distribution is inverse J not surprising
subplot_ssd.df %>%
  group_by(bin_midpoint_cm, sizeclass) %>%
  summarise(prop.med = median(class_proportion),
            prop.05 = quantile(class_proportion, 0.05),
            prop.95 = quantile(class_proportion, 0.95)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint_cm))+
  geom_point(aes(y = prop.med), size = 3)+
  geom_errorbar(aes(ymin = prop.05, ymax = prop.95),
                width = 1)+
  theme_minimal()

#### subplot reproductive value ################################################
# reproductive value
subplot_repro.med = 
  matrix(nrow = nrow(size_metadata),
         ncol = nrow(subplots.pila),
         byrow = FALSE,
         data = 
           sapply(X = 1:nrow(subplots.pila),
                  FUN = function(subplot){
                    A.subplot = subplot_transitions.med[,,subplot]
                    # from supplementary materials for merow et al 2014 
                    # "On using integral projection models..."
                    v.eigen = Re(eigen(t(A.subplot))$vectors[,1])
                    
                    rv = v.eigen / v.eigen[1]
                    return(rv)
                  }))


subplot_repro.df = 
  expand.grid('subplot' = 1:nrow(subplots.pila),
              'sizeclass' = 1:nrow(size_metadata)) %>%
  left_join(size_metadata %>%
              mutate(bin_midpoint_cm = bin_midpoint*100) %>%
              select(sizeclass = bin_id, bin_midpoint_cm))

subplot_repro.df$reproductive_value = 
  
  as.numeric(
    sapply(X = 1:nrow(size_metadata),
           FUN = function(sizeclass){
             
             return(subplot_repro.med[sizeclass,])
             
           })
  )


subplot_repro.df %>%
  filter(!is.element(subplot, c(1150, 3000))) %>%
  summary()

# stable size distribution is inverse J not surprising
subplot_repro.df %>%
  
  group_by(bin_midpoint_cm, sizeclass) %>%
  
  # there's a couple of NA subplots for reproductive value, looks like cases 
  # where numerical errors are resulting in a divide by zero when going from 
  # v.eigen to reproductive value? 
  summarise(repr.med = median(reproductive_value, na.rm = TRUE),
            
            # there's a couple of NAs in h
            repr.05 = quantile(reproductive_value, 0.25, na.rm = TRUE),
            repr.95 = quantile(reproductive_value, 0.75, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint_cm))+
  geom_point(aes(y = repr.med), size = 3)+
  #geom_errorbar(aes(ymin = repr.05, ymax = repr.95),
  #              width = 1)+
  theme_minimal()


#### subplot sensitivity and elasticity ########################################
# sensitivity and elasticity
v.dot.w = 
  
  sapply(X = 1:nrow(subplots.pila),
         FUN = function(subplot){
           sum(subplot_ssd.med[,subplot] * subplot_repro.med[,subplot])*0.127
         })


test = outer(subplot_repro.med[,1], subplot_ssd.med[,1])/v.dot.w[1]


sens = 
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(subplots.pila)),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'subplot' = 1:nrow(subplots.pila)),
        data = 
          sapply(X = 1:nrow(subplots.pila),
                 FUN = function(subplot){
                   outer(subplot_repro.med[,subplot], subplot_ssd.med[,subplot])/
                     v.dot.w[subplot]
                 }))

sens.agg = 
  matrix(nrow = nrow(size_metadata),
         ncol = nrow(size_metadata),
         byrow = FALSE,
         data = 
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_from){
                    sapply(X = 1:nrow(size_metadata),
                           FUN = function(size_to){
                             median(as.vector(sens[size_to,size_from,]),
                                    na.rm= TRUE)
                           })
                  }))

sens.df = 
  expand.grid(size_to = size_metadata$bin_id,
              size_from = size_metadata$bin_id) %>%
  left_join(size_metadata %>%
              select(size_to = bin_id, 
                     bin_midpoint_to_m = bin_midpoint)) %>%
  left_join(size_metadata %>%
              select(size_from = bin_id,
                     bin_midpoint_from_m = bin_midpoint)) 

# looks like its transitions from the biggest classes into the smallest classes
# that matter most; edit: well now I'm confused, was I swapping the rows and 
# columns before or am I swapping them now? this doesn't look biologically 
# reasonable, but it does match the demo sensitivity in the merow paper
# where its from the small class into the big class that has the highest 
# sensitivity, now I'm wondering if the merow code has a bug?
fields::image.plot(size_metadata$bin_midpoint,
                   size_metadata$bin_midpoint,
                   t(sens[,,4]),
                   xlab = 'Size (t)', ylab = 'Size (t+1)')

sens.df$sensitivity = 
  sapply(X = 1:nrow(size_metadata),
         FUN = function(size_from){
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    sens.agg[size_to, size_from]
                  })
         }) %>%
  as.vector()

ggplot(sens.df,
       aes(x = bin_midpoint_from_m,
           y = bin_midpoint_to_m,
           fill = sensitivity))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c()


sens = 
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(subplots.pila)),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'subplot' = 1:nrow(subplots.pila)),
        data = 
          sapply(X = 1:nrow(subplots.pila),
                 FUN = function(subplot){
                   outer(subplot_repro.med[,subplot], subplot_ssd.med[,subplot])/
                     v.dot.w[subplot]
                 }))


elas = 
  
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(subplots.pila)),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'subplot' = 1:nrow(subplots.pila)),
        data = 
          sapply(X = 1:nrow(subplots.pila),
                 FUN = function(subplot){
                   matrix(as.vector(sens[,,subplot])*
                            as.vector(subplot_transitions.med[,,subplot])/
                            subplot_lambdas.med[subplot])
                 }))

elas.agg = 
  matrix(nrow = nrow(size_metadata),
         ncol = nrow(size_metadata),
         byrow = FALSE,
         data = 
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    sapply(X = 1:nrow(size_metadata),
                           FUN = function(size_from){
                             median(as.vector(elas[size_to,size_from,]),
                                    na.rm= TRUE)
                           })
                  }))

elas.df = 
  expand.grid(size_to = size_metadata$bin_id,
              size_from = size_metadata$bin_id) %>%
  left_join(size_metadata %>%
              select(size_to = bin_id, 
                     bin_midpoint_to_m = bin_midpoint)) %>%
  left_join(size_metadata %>%
              select(size_from = bin_id,
                     bin_midpoint_from_m = bin_midpoint)) 

# looks like its transitions from the biggest classes into the smallest classes
# that matter most
fields::image.plot(size_metadata$bin_midpoint,
                   size_metadata$bin_midpoint,
                   t(elas.agg),
                   xlab = 'Size (t)', ylab = 'Size (t+1)')

elas.df$elasticity = 
  sapply(X = 1:nrow(size_metadata),
         FUN = function(size_from){
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    elas.agg[size_to, size_from]
                  })
         }) %>%
  as.vector()

ggplot(elas.df,
       aes(x = bin_midpoint_from_m,
           y = bin_midpoint_to_m,
           fill = elasticity))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c()





v = Re(eigen(t(subplot_transitions.med[,,49]))$vectors[,1])



test = subplots.pila

bad_subplots = 
  subplot_repro.df %>%
  filter(reproductive_value<0) %>%
  group_by(subplot) %>%
  summarise() %>%
  ungroup() %>%
  pull(subplot)
test$subp_number = 1:nrow(subplots.pila)
head(test)
test$bad_subplot = 
  sapply(X = test$subp_number,
         FUN = function(s){
           is.element(s, bad_subplots)
         })

test %>%
  filter(bad_subplot) %>%
  print(width = Inf, n = Inf)

test %>% 
  print(width = Inf, n = Inf)

test %>%
  ggplot(aes(x = bad_subplot, y = lambda_postmed))+
  geom_boxplot()+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 2))

#### mapping lambda ############################################################

library(sf)
library(scales)
subplots.pila.sf = 
  subplots.pila %>%
  st_as_sf(coords = c('lon', 'lat'),
           crs = 4269)

ggplot(data = 
         subplots.pila.sf,
       aes(color = lambda_postmed))+
  geom_sf(size = 1)+
  scale_color_viridis_c(limits = c(0.9,1.5), oob = scales::squish)+
  theme_minimal()

ggplot(data = 
         subplots.pila.sf %>%
         mutate(lambda_binary = lambda_postmed>=1),
       aes(color = lambda_binary))+
  geom_sf(size = 1, alpha = 0.5)+
  #scale_color_viridis_d(begin = 0.25, end = 0.85, option = 'C')+
  theme_minimal()




A_hypotheticals = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'hypothetical_As.rds'))

#### hypothetical lambda #######################################################

hypothetical_lambdas.matrix = 
  matrix(nrow = nrow(hypothetical_subplots),
         ncol = 4000,
         data = 
           sapply(X = 1:nrow(hypothetical_subplots),
                  FUN = function(subplot){
                    sapply(X = 1:4000,
                           FUN = function(draw){
                             max(Re(eigen(A_hypotheticals[,,subplot,draw])$values))
                           })
                  }))

hypothetical_lambdas = 
  hypothetical_subplots %>%
  expand(nesting(subp_id, name, intercept, fire, wpbr, ba_scaled, 
                 cwd_dep90_scaled, cwd_mean_scaled),
         data.frame(draw = 1:4000))


hypothetical_lambdas$lambda = 
  sapply(X = 1:nrow(hypothetical_subplots),
         FUN = function(subplot){
           sapply(X = 1:4000,
                  FUN = function(draw){
                    # paste0('s:',subplot,'d:',draw) for testing
                    max(as.numeric(eigen(A_hypotheticals[,,subplot,draw])$values))
                  })
         }) %>%
  as.numeric()

saveRDS(hypothetical_lambdas,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'hypothetical_lambdas.rds'))

hypothetical_lambdas = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'hypothetical_lambdas.rds'))



pretty_names = 
  hypothetical_subplots$name
names(pretty_names) = hypothetical_subplots$subp_id

head(hypothetical_lambdas)

hypothetical_lambda_summary = 
  hypothetical_lambdas %>%
  group_by(subp_id, name) %>%
  summarise(lambda.med = median(lambda),
            lambda.05 = quantile(lambda, probs = 0.05),
            lambda.95 = quantile(lambda, probs = 0.95))

hypothetical_lambda_summary

write.csv(hypothetical_lambda_summary,
          here::here('04-communication', 
                     'tables',
                     'hypothetical_lambdas_summary.csv'),
          row.names = FALSE)

hypothetical_lambdas_plot = 
  ggplot(data = 
           hypothetical_lambdas,
         aes(x = lambda))+
  geom_density(lwd = 1)+
  geom_vline(xintercept = 1, color = 'grey', lty = 2, lwd = 1)+
  theme_minimal()+
  facet_grid(subp_id~., scales = 'free_y',
             labeller = labeller(subp_id = pretty_names))+
  scale_x_continuous(limits = c(0.5,3))+
  coord_cartesian(xlim = c(0.9, 1.5))+
  theme(axis.text.y = element_blank())+
  labs(x = 'Lambda', y = 'Posterior Density')

hypothetical_lambdas_plot

ggsave(hypothetical_lambdas_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'hypotheticals_lambda_post.png'),
       height = 7.5, width = 4, units = 'in')

hypothetical_lambdas_ppt = 
  ggplot(data = 
           hypothetical_lambdas,
         aes(x = lambda))+
  geom_density(lwd = 1)+
  geom_vline(xintercept = 1, color = 'grey', lty = 2, lwd = 1)+
  theme_minimal()+
  facet_grid(.~subp_id, scales = 'free_x',
             labeller = labeller(subp_id = pretty_names))+
  scale_x_continuous(limits = c(0.5,3))+
  theme(axis.text.x = element_blank(),
        text = element_text(size = 16))+
  labs(x = 'Lambda', y = 'Posterior Density')+
  coord_flip(xlim = c(0.9, 1.4))

hypothetical_lambdas_ppt

ggsave(hypothetical_lambdas_ppt,
       filename = here::here('04-communication',
                             'figures',
                             'powerpoint',
                             'hypotheticals_lambda_post.png'),
       height = 4, width = 11, units = 'in')


#### hypothetical stable size distribution #####################################
dim(A_hypotheticals)
hypothetical_ssd = 
  array(dim = c(nrow(size_metadata),
                nrow(hypothetical_subplots),
                4000),
        dimnames = 
          list('sizeclass' = 1:nrow(size_metadata),
               'subplot' = 1:nrow(hypothetical_subplots),
               'draw' = 1:4000),
        data = 
          sapply(X = 1:4000,
                 FUN = function(draw){
                   sapply(X = 1:nrow(hypothetical_subplots),
                          FUN = function(subplot){
                            A = A_hypotheticals[,,subplot,draw]
                            w.eigen = Re(eigen(A)$vectors[,1])
                            ssd = w.eigen / sum(w.eigen)
                            return(ssd)
                          })
                 }))


hypothetical_ssd.df = 
  hypothetical_subplots %>%
  expand(nesting(intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled,
                 cwd_mean_scaled, subp_id, name),
         sizeclass = 1:nrow(size_metadata),
         draw = 1:4000) %>%
  left_join(size_metadata %>%
              mutate(bin_midpoint_cm = bin_midpoint*100) %>%
              select(sizeclass = bin_id, bin_midpoint_cm)) %>%
  arrange(subp_id, sizeclass, draw)

head(hypothetical_ssd.df)

hypothetical_ssd.df$class_proportion = 
  
  as.numeric(
    sapply(X = 1:nrow(hypothetical_subplots),
           FUN = function(subplot){
             
             as.numeric(
               sapply(X = 1:nrow(size_metadata),
                      FUN = function(sizeclass){
                        
                        return(hypothetical_ssd[sizeclass,subplot,])
                      })
             )
             
           }))

# Fire makes the SSD really for the very rare bigger size classes, fewer 
# mid-to-large trees and more superlarge 
hypothetical_ssd.df %>%
  group_by(bin_midpoint_cm, sizeclass, name, subp_id) %>%
  summarise(prop.med = median(class_proportion),
            prop.05 = quantile(class_proportion, 0.05),
            prop.95 = quantile(class_proportion, 0.95)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint_cm, color = name, fill = name))+
  #geom_point(aes(y = prop.med), size = 1)+
  geom_line(aes(y = prop.med))+
  geom_ribbon(aes(ymin = prop.05, ymax = prop.95), alpha = 0.25)+
  theme_minimal()+
  scale_y_log10()+
  facet_wrap(~name)

#### hypothetical reproductive value ###########################################

hypothetical_repro = 
  array(dim = c(nrow(size_metadata),
                nrow(hypothetical_subplots),
                4000),
        dimnames = 
          list('sizeclass' = 1:nrow(size_metadata),
               'subplot' = 1:nrow(hypothetical_subplots),
               'draw' = 1:4000),
        data = 
          sapply(X = 1:4000,
                 FUN = function(draw){
                   sapply(X = 1:nrow(hypothetical_subplots),
                          FUN = function(subplot){
                            A = A_hypotheticals[,,subplot,draw]
                            v.eigen = Re(eigen(t(A))$vectors[,1])
                            rv = v.eigen / v.eigen[1]
                            return(rv)
                          })
                 }))


hypothetical_repro.df = 
  hypothetical_subplots %>%
  expand(nesting(intercept, fire, wpbr, ba_scaled, cwd_dep90_scaled,
                 cwd_mean_scaled, subp_id, name),
         sizeclass = 1:nrow(size_metadata),
         draw = 1:4000) %>%
  left_join(size_metadata %>%
              mutate(bin_midpoint_cm = bin_midpoint*100) %>%
              select(sizeclass = bin_id, bin_midpoint_cm)) %>%
  arrange(subp_id, sizeclass, draw)

hypothetical_repro.df$repro = 
  
  as.numeric(
    sapply(X = 1:nrow(hypothetical_subplots),
           FUN = function(subplot){
             
             as.numeric(
               sapply(X = 1:nrow(size_metadata),
                      FUN = function(sizeclass){
                        
                        return(hypothetical_repro[sizeclass,subplot,])
                      })
             )
             
           }))

# Again these reproductive values are wonky, esp for burned subplots. Something 
# about the transition matrix for burned subplots is weird.
hypothetical_repro.df %>%
  group_by(bin_midpoint_cm, sizeclass, name, subp_id) %>%
  #filter(subp_id != 2) %>%
  summarise(repro.med = median(repro),
            repro.05 = quantile(repro, 0.05, na.rm = TRUE),
            repro.95 = quantile(repro, 0.95, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = bin_midpoint_cm, color = name, fill = name))+
  #geom_point(aes(y = prop.med), size = 1)+
  geom_line(aes(y = repro.med))+
  #geom_ribbon(aes(ymin = repro.05, ymax = repro.95), alpha = 0.25)+
  theme_minimal()+
  facet_wrap(~name, scales = 'free_y')

# ok this kind of makes sense: under disturbances that really reduce the 
# survival of small stems (WPBR and esp fire) the reproductive value of 
# big stems relative to little ones is magnified, because so few little ones 
# survive to become real contributors to reproduction; still think theres
# some numerical instability or smth causing the credible interval bounds for 
# reproductive value to go wonky for fire, all the others look fine


#### hypothetical sensitivity and elasticity ###################################

# sensitivity and elasticity
hypothetical_vdotw = 
  matrix(nrow = nrow(hypothetical_subplots),
         ncol = 4000,
         data = 
           sapply(X = 1:nrow(hypothetical_subplots),
                  FUN = function(subplot){
                    sapply(X = 1:4000,
                           FUN = function(draw){
                             sum(hypothetical_ssd[,subplot,draw] *
                                   hypothetical_repro[,subplot,draw])*0.127
                           })
                  }))


hypothetical_sens = 
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(hypothetical_subplots),
                   4000),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'subplot' = 1:nrow(hypothetical_subplots),
               'draw' = 1:4000),
        data = 
          
          sapply(X = 1:4000,
                 FUN = function(draw){
                   sapply(X = 1:nrow(hypothetical_subplots),
                          FUN = function(subplot){
                            outer(hypothetical_repro[,subplot,draw], 
                                  hypothetical_ssd[,subplot,draw])/
                              hypothetical_vdotw[subplot,draw]
                          })
                 })
  )

hypothetical_elas = 
  array(dim = list(nrow(size_metadata),
                   nrow(size_metadata),
                   nrow(hypothetical_subplots),
                   4000),
        dimnames = 
          list('size_to' = 1:nrow(size_metadata),
               'size_from' = 1:nrow(size_metadata),
               'subplot' = 1:nrow(hypothetical_subplots),
               'draw' = 1:4000),
        data = 
          
          sapply(X = 1:4000,
                 FUN = function(draw){
                   sapply(X = 1:nrow(hypothetical_subplots),
                          FUN = function(subplot){
                            
                            matrix(as.vector(hypothetical_sens[,,subplot,draw])*
                                     as.vector(A_hypotheticals[,,subplot,draw])/
                                     hypothetical_lambdas.matrix[subplot,draw])
                            
                            
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
         draw = 1:4000) %>%
  arrange(subp_id, size_to, size_from, draw)


hypothetical_sens_elas.df$sensitivity = 
  sapply(X = 1:nrow(hypothetical_subplots),
         FUN = function(subplot){
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    
                    sapply(X = 1:nrow(size_metadata),
                           FUN = function(size_from){
                             hypothetical_sens[size_to, size_from, subplot,]
                           })
                  })
         }) %>%
  as.vector()


hypothetical_sens_elas.df$elasticity = 
  sapply(X = 1:nrow(hypothetical_subplots),
         FUN = function(subplot){
           sapply(X = 1:nrow(size_metadata),
                  FUN = function(size_to){
                    
                    sapply(X = 1:nrow(size_metadata),
                           FUN = function(size_from){
                             hypothetical_elas[size_to, size_from, subplot,]
                           })
                  })
         }) %>%
  as.vector()

hypothetical_sens_elas.df %>%
  group_by(size_to, size_from, bin_midpoint_to_m, bin_midpoint_from_m,
           subp_id) %>%
  summarise(sensitivity = median(sensitivity, na.rm = TRUE),
            elasticity = median(elasticity, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(hypothetical_subplots %>%
              select(subp_id, name)) %>%
  filter(subp_id==2) %>%
  ggplot(aes(x = bin_midpoint_from_m,
             y = bin_midpoint_to_m,
             fill = sensitivity))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c()


hypothetical_sens_elas.df %>%
  group_by(size_to, size_from, bin_midpoint_to_m, bin_midpoint_from_m,
           subp_id) %>%
  summarise(sensitivity = median(sensitivity, na.rm = TRUE),
            elasticity = median(elasticity, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(hypothetical_subplots %>%
              select(subp_id, name)) %>%
  filter(subp_id==3) %>%
  ggplot(aes(x = bin_midpoint_from_m,
             y = bin_midpoint_to_m,
             fill = elasticity))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()+
  scale_fill_viridis_c(limits = c(0,1))

fields::image.plot(size_metadata$bin_midpoint,
                   size_metadata$bin_midpoint,
                   t(A_hypotheticals[,,2,1000]),
                   xlab = 'Size (t)', ylab = 'Size (t+1)')


