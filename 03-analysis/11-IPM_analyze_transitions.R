library(here)
library(tidyverse)



#### subplot tranasition matrices ############################################## 



subplot_As = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'subplot_As.rds'))

#### subplot lambda ############################################################


subplot_df = 
  expand.grid(draw = 1:100,
         subplot = 1:dim(subplot_As)[3]) %>%
  as_tibble() %>%
  arrange(draw, subplot)

subplot_df$lambda = 
  sapply(X = 1:100,
         FUN = function(draw){
           sapply(X = 1:dim(subplot_As)[3],
                  FUN = function(subplot){
                    A.subplot = subplot_As[,,subplot,draw]
                    lambda.subplot = max(as.numeric(Re(eigen(A.subplot)$values)))
                    return(lambda.subplot)
                  })
         }) %>%
  as.vector()

postmed_lambda_distribution = 
  ggplot(data = subplot_df,
         aes(x = log(lambda), group = draw))+
  geom_line(stat = 'density', alpha = 0.1)+
  theme_minimal()+
  #scale_x_continuous(limits = c(-1, 1))+
  geom_vline(xintercept = 0, color = 'grey', lty = 2, lwd = 1)+
  labs(y = 'Density', x = 'Lambda')

postmed_lambda_distribution

subplot_df %>%
  group_by(draw) %>%
  summarise(lambda.med = median(lambda, na.rm = TRUE),
            lambda.05 = quantile(lambda, probs = 0.05),
            lambda.95 = quantile(lambda, probs = 0.95)) %>%
  ungroup() %>%
  ggplot(aes(x = lambda.med))+
  geom_density()


#### parameter sensitivities ###################################################

head(subplot_df)

posterior = readRDS(here::here('02-data',
                               '03-results',
                               'real_fits',
                               'posterior_draws.rds'))

posterior %>%
  select(contains(match = c('beta', 'sigma')))

ipm_results = 
  subplot_df %>%
  left_join(posterior %>%
              select(contains(match = c('beta', 'sigma'))) %>%
              mutate(draw = posterior$.draw),
            by = c('draw' = 'draw'))

saveRDS(ipm_results,
        here::here('02-data',
                   '03-results',
                   'real_fits',
                   'real_subplot_lambda_all_draws.rds'))

ipm_results = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'real_subplot_lambda_all_draws.rds'))

head(ipm_results)


beta_s_sensitivities = 
  lapply(X = 1:12,
         FUN = function(p){
           
           param_name = paste0('beta_s[',p,']')
           
           ipm_results %>%
              mutate(log_lambda = log(lambda)) %>%
              group_by(draw) %>%
              summarise(llambda05 = quantile(log_lambda, 0.05),
                        llambda25 = quantile(log_lambda, 0.25),
                        llambdamed = median(log_lambda),
                        llambda75 = quantile(log_lambda, 0.75),
                        llambda95 = quantile(log_lambda, 0.95)) %>%
              ungroup() %>%
              left_join(posterior %>%
                          select(contains(match = c('beta', 'sigma'))) %>%
                          mutate(draw = posterior$.draw),
                        by = c('draw' = 'draw')) %>%
              ggplot(aes(x = .data[[param_name]]))+
              geom_ribbon(aes(ymin = llambda05, ymax = llambda95), alpha = 0.25)+
              geom_ribbon(aes(ymin = llambda25, ymax = llambda75), alpha = 0.5)+
              geom_line(aes(y = llambdamed))+
              geom_smooth(data = 
                            ipm_results %>%
                            mutate(log_lambda = log(lambda)),
                          aes(x = .data[[param_name]], y = log_lambda),
                          method = 'lm')+
              theme_minimal()
         })

beta_s_sensitivities

beta_g_sensitivities = 
  lapply(X = 1:12,
         FUN = function(p){
           
           param_name = paste0('beta_g[',p,']')
           
           ipm_results %>%
              mutate(log_lambda = log(lambda)) %>%
              group_by(draw) %>%
              summarise(llambda05 = quantile(log_lambda, 0.05),
                        llambda25 = quantile(log_lambda, 0.25),
                        llambdamed = median(log_lambda),
                        llambda75 = quantile(log_lambda, 0.75),
                        llambda95 = quantile(log_lambda, 0.95)) %>%
              ungroup() %>%
              left_join(posterior %>%
                          select(contains(match = c('beta', 'sigma'))) %>%
                          mutate(draw = posterior$.draw),
                        by = c('draw' = 'draw')) %>%
              ggplot(aes(x = .data[[param_name]]))+
              geom_ribbon(aes(ymin = llambda05, ymax = llambda95), alpha = 0.25)+
              geom_ribbon(aes(ymin = llambda25, ymax = llambda75), alpha = 0.5)+
              geom_line(aes(y = llambdamed))+
              geom_smooth(data = 
                            ipm_results %>%
                            mutate(log_lambda = log(lambda)),
                          aes(x = .data[[param_name]], y = log_lambda),
                          method = 'lm')+
              theme_minimal()
         })

beta_g_sensitivities

beta_f_sensitivities = 
  lapply(X = 1:12,
         FUN = function(p){
           
           param_name = paste0('beta_f[',p,']')
           
           ipm_results %>%
              mutate(log_lambda = log(lambda)) %>%
              group_by(draw) %>%
              summarise(llambda05 = quantile(log_lambda, 0.05),
                        llambda25 = quantile(log_lambda, 0.25),
                        llambdamed = median(log_lambda),
                        llambda75 = quantile(log_lambda, 0.75),
                        llambda95 = quantile(log_lambda, 0.95)) %>%
              ungroup() %>%
              left_join(posterior %>%
                          select(contains(match = c('beta', 'sigma'))) %>%
                          mutate(draw = posterior$.draw),
                        by = c('draw' = 'draw')) %>%
              ggplot(aes(x = .data[[param_name]]))+
              geom_ribbon(aes(ymin = llambda05, ymax = llambda95), alpha = 0.25)+
              geom_ribbon(aes(ymin = llambda25, ymax = llambda75), alpha = 0.5)+
              geom_line(aes(y = llambdamed))+
              geom_smooth(data = 
                            ipm_results %>%
                            mutate(log_lambda = log(lambda)),
                          aes(x = .data[[param_name]], y = log_lambda),
                          method = 'lm')+
              theme_minimal()
         })

beta_f_sensitivities

head(ipm_results)

ipm_results_long = 
  ipm_results %>%
  mutate(log_lambda = log(lambda)) %>%
  pivot_longer(cols = c('beta_s[1]', 'beta_s[2]', 'beta_s[3]', 'beta_s[4]',
                        'beta_s[5]', 'beta_s[6]', 'beta_s[7]', 'beta_s[8]',
                        'beta_s[9]', 'beta_s[10]', 'beta_s[11]', 'beta_s[12]',
                        'beta_g[1]', 'beta_g[2]', 'beta_g[3]', 'beta_g[4]',
                        'beta_g[5]', 'beta_g[6]', 'beta_g[7]', 'beta_g[8]',
                        'beta_g[9]', 'beta_g[10]', 'beta_g[11]', 'beta_g[12]',
                        'beta_f[1]', 'beta_f[2]', 'beta_f[3]', 'beta_f[4]',
                        'beta_f[5]', 'beta_f[6]', 'beta_f[7]', 'beta_f[8]',
                        'beta_f[9]', 'beta_f[10]', 'beta_f[11]', 'beta_f[12]'),
               names_to = 'parameter',
               values_to = 'parameter_value') %>%
  select(draw, subplot, log_lambda, parameter, parameter_value)

head(ipm_results_long)

fixeff_sensitivities = 
  ipm_results_long %>%
  group_by(parameter) %>%
  summarise() %>%
  ungroup()

library(lme4)
fixeff_sensitivity_fits = 
  
  lapply(X = fixeff_sensitivities$parameter,
         FUN = function(p){
           lm_data = 
             ipm_results_long %>%
             filter(parameter==p) %>%
             mutate(parameter_value_s = scale(parameter_value))
           
           lm_fit = lme4::lmer(data = lm_data,
                               log_lambda ~ parameter_value_s + (1|subplot))
           
           return(lm_fit)
         })

fixeff_sensitivities$beta0 = 
  sapply(X = fixeff_sensitivity_fits,
         FUN = function(fit){
           fit@beta[1]
         })

fixeff_sensitivities$beta1 = 
  sapply(X = fixeff_sensitivity_fits,
         FUN = function(fit){
           fit@beta[2]
         })

fixeff_sensitivities$beta1_se = 
  sapply(X = fixeff_sensitivity_fits,
         FUN = function(fit){
           summary(fit)$coefficients[2,2]
         })


beta_names =   
  data.frame(
      param = c('beta_s[1]', 'beta_s[2]', 'beta_s[3]', 'beta_s[4]', 'beta_s[5]',
                'beta_s[6]', 'beta_s[7]', 'beta_s[8]', 'beta_s[9]', 'beta_s[10]',
                'beta_s[11]', 'beta_s[12]',
                'beta_g[1]', 'beta_g[2]', 'beta_g[3]', 'beta_g[4]', 'beta_g[5]',
                'beta_g[6]', 'beta_g[7]', 'beta_g[8]', 'beta_g[9]', 'beta_g[10]',
                'beta_g[11]', 'beta_g[12]',
                'beta_f[1]', 'beta_f[2]', 'beta_f[3]', 'beta_f[4]', 'beta_f[5]',
                'beta_f[6]', 'beta_f[7]', 'beta_f[8]', 'beta_f[9]', 'beta_f[10]',
                'beta_f[11]', 'beta_f[12]'),
      pretty_name = 
          c('Surv-Intercept','Surv-DBH (m)','Surv-Fire','Surv-WPBR','Surv-Basal Area',
            'Surv-Drought','Surv-Site Dryness','Surv-DBH x Fire','Surv-DBH x WPBR',
            'Surv-DBH x BA','Surv-DBH x Drought','Surv-DBH x Dryness',
            'Grow-Intercept','Grow-DBH (m)','Grow-Fire','Grow-WPBR','Grow-Basal Area',
            'Grow-Drought','Grow-Site Dryness','Grow-DBH x Fire','Grow-DBH x WPBR',
            'Grow-DBH x BA','Grow-DBH x Drought','Grow-DBH x Dryness',
            'Fec-Intercept','Fec-DBH (m)','Fec-Fire','Fec-WPBR','Fec-Basal Area',
            'Fec-Drought','Fec-Site Dryness','Fec-DBH x Fire','Fec-DBH x WPBR',
            'Fec-DBH x BA','Fec-DBH x Drought','Fec-DBH x Dryness'))


fixeff_sensitivities %>%
  left_join(beta_names, by = c('parameter' = 'param')) %>%
  ggplot(aes(y = reorder(pretty_name, abs(beta1))))+
  geom_point(aes(x = beta1))+
  geom_errorbarh(aes(xmin = beta1-(2*beta1_se),
                    xmax = beta1+(2*beta1_se)))


# this is interesting; 
# x axis shows the change in lambda associated with increasing the value of 
# some parameter by 1 SD of the posterior distribution
# so parameters which are highly influential *and/or highly uncertain* have 
# large values
# largest values are for fecundity intercept and fecundity size effect,
# then interaction of size and dryness on fecundity, effect of dryness on 
# growth, interaction of size and BA on survival



#### subplot stable size distribution ##########################################

# stable size distribution
subplot_ssd.med = 
  matrix(nrow = nrow(size_metadata),
         ncol = nrow(subplots.pila),
         byrow = FALSE,
         data = 
           sapply(X = 1:nrow(subplots.pila),
                  FUN = function(subplot){
                    A.subplot = subplot_As[,,subplot]
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
                    A.subplot = subplot_As[,,subplot]
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
                            as.vector(subplot_As[,,subplot])/
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





v = Re(eigen(t(subplot_As[,,49]))$vectors[,1])



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
# Fecival of small stems (WPBR and esp fire) the reproductive value of 
# big stems relative to little ones is magnified, because so few little ones 
# Fecive to become real contributors to reproduction; still think theres
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


