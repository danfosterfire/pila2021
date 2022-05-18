data {
  int<lower=1> K; // number of fixed effect covariates, incl intercept and interactions
  int<lower=1> P; // number of unique plots
  int<lower=1> E; // number of unique ecoregion subsections
  
  // surivval data
  int<lower=0> N_s; // number of individual trees tracked for survival
  int surv[N_s]; // 1 if individual i survived from t = 0 to t = 1, 0 otherwise
  int plotid_s[N_s]; // plot ID indices
  int ecosub_s[N_s]; // ecoregion subsection indices
  matrix[N_s,K] X_s; // design matrix for fixed effects
  
  // growth data
  int<lower=0> N_g; // number of individual trees tracked for growth
  real size1_g[N_g]; // sizes at time t+1
  int plotid_g[N_g]; // plot id indices
  int ecosub_g[N_g]; // ecoregion subsection indices
  matrix[N_g,K] X_g; // design matrix for fixed effects, including intercept
  
  // recruitment data
  int<lower=0> N_r; // M_r (20) times the number of unique plots for the 
  // recruitment submodel
  int<lower=0> P_r; // number of unique plots for the recruitment submodel
  matrix[N_r, K] X_r; // fixeff explanatory variables for each sizeclass:plot 
  // combination, used to predict growth and survival for the recruitment 
  // submodel sizeclasses; instead of observed sizes here we have the 
  // class_mean_dbhs for each size class; MUST BE ORDERED plot (slow) sizeclass (fast)
  int plotid_r[N_r]; // plot id indices for survival random effects in IPM
  int ecosub_r[N_r]; // ecosub indices for survival random effects in IPM
  int<lower=0> M_r; // number of modeled size classes for the recruitment submodel
  vector[M_r] u_bounds; // upper bounds of size bins
  vector[M_r] l_bounds; //lower bounds of size bins
  vector[2] a; // plot area for each of the size classes; here only for the 
  // smallest 2 size classes (only need these, and the macro breakpoint diameter is not consistent 
  // across all the plots)
  int cprime[2,P_r]; // counts of untagged trees in each size class on each plot
  vector<lower=0>[M_r] n[P_r]; // vector of area-standardaized rates of occurence for 
  // each of size class on each of P_r plots at time t = 0
  vector<lower=0>[2] r; // recruitment size kernel
  
}


parameters {
  // survival
  vector[K] beta_s; // fixeff coefficients, including intercept, for surv model
  real<lower=0> sigmaPlot_s; // standard deviation of plot random effect
  real<lower=0> sigmaEco_s; // sddev of ecoregion subsection effect
  vector[P] zPlot_s; // standard normal deviates for plot random effect
  vector[E] zEco_s; // sdnormal deviates for ecoregion random effect
  
  // growth 
  vector[K] beta_g; // fixeff coeffs, including intercept, for growth model
  real<lower=0> sigmaPlot_g; // std deviation of plot random effect
  vector[P] zPlot_g; // std normal deviates for plot random effect
  real<lower=0> sigmaEco_g; // std dev of ecoregion random effect
  vector[E] zEco_g; // stdnormal deviates
  real<lower=0> sigmaEpsilon_g; // std deviation of residuals for growth model
  
  // recruitment
  real<lower=0> kappa_r; // neg binomial dispersion parameter for recruitment
  vector[K] beta_f; // fixeff coefficients including intercept for fecundity model
  real<lower=0> sigmaPlot_f;
  vector[P] zPlot_f;
  real<lower=0> sigmaEco_f;
  vector[E] zEco_f;
}

transformed parameters {
  vector[P] plotEffect_s;
  vector[P] plotEffect_g;
  vector[P] plotEffect_f;
  vector[E] ecoEffect_s;
  vector[E] ecoEffect_g;
  vector[E] ecoEffect_f;
  
  plotEffect_s = sigmaPlot_s * zPlot_s;
  plotEffect_g = sigmaPlot_g * zPlot_g;
  plotEffect_f = sigmaPlot_f * zPlot_f;
  ecoEffect_s = sigmaEco_s * zEco_s;
  ecoEffect_g = sigmaEco_g * zEco_g;
  ecoEffect_f = sigmaEco_f * zEco_f;
}

model {
  
  // variable declarations
  // // survival
  vector[N_s] logitp_s;
  vector[N_s] XB_s;
  
  // // growth 
  vector[N_g] mu_g;
  vector[N_g] XB_g;
  
  // // recruitment
  vector[N_r] XBP_r; // linear predictor for sizeclass class_mean_dbhs on each plot
  vector[N_r] XBg_r; // linear predictor for sizeclass class_mean_dbhs on each plot
  vector[N_r] XBf_r;
  vector[N_r] mu_gr;
  vector[N_r] logitp_sr;
  vector[2] nprime[P_r]; // area-standardized occurence rates at time t+1 on 
  // each plot in each size class (only the smallest 2 classes)
  matrix[2,M_r] A[P_r]; // final IPM transition kernel describing how column 
  // size classes at time 0 lead to row size classes at time 1 (only smallest 2 classes)
  matrix[2,M_r] recKern[P_r]; // recruitment kernel (into the smallest 2 classes)
  matrix[2,M_r] growKern[P_r]; // growth*survival transitions (into the smallest 2 classes)
  matrix[M_r,1] g[P_r]; // growth kernel of transitions from size class to size class on each plot
                        // (only from the smallest class into the others)
  matrix[M_r,P_r] s; // survival rates on each plot for each sizeclass 
  matrix[M_r,P_r] f; // fecundity rates by size class and plots
  vector[N_r] logf; 
  
  // fixed effects
  XB_s = X_s * beta_s;
  XB_g = X_g * beta_g;
  XBP_r = X_r * beta_s;
  XBg_r = X_r * beta_g;
  XBf_r = X_r * beta_f;
  
  
  // linear predictor for survival model
  for (i in 1:N_s){
    logitp_s[i] = XB_s[i] + plotEffect_s[plotid_s[i]] + ecoEffect_s[ecosub_s[i]];
  }
  
  // linear predictor for growth (size2) model
  for (i in 1:N_g){
    mu_g[i] = XB_g[i] + plotEffect_g[plotid_g[i]] + ecoEffect_g[ecosub_g[i]];
  }
  
  // linear predictors for gtrowth and survival on the recruitment IPM
  for (i in 1:N_r){
    logitp_sr[i] = XBP_r[i] + plotEffect_s[plotid_r[i]]+ ecoEffect_s[ecosub_r[i]];
    mu_gr[i] = XBg_r[i] + plotEffect_g[plotid_r[i]] + ecoEffect_g[ecosub_r[i]];
    logf[i] = XBf_r[i] + plotEffect_f[plotid_r[i]] + ecoEffect_f[ecosub_r[i]];
  }
  

  //print("r(2): ", r);
  // IPM model for recruitment; loop over all the plots
  for (plot in 1:P_r){
    
          
    // expected survival in the smallest size class
    s[1,plot] = inv_logit(logitp_sr[1+(M_r*(plot-1))]);
    
    // in the shriver code, this loops over all combinations of size classes,
    // which isn't necessary because we're only using growth from the smallest 
    // size class into the smallest two size classes;  set "sizeclass_from" to 
    // be 1 everywhere and sizeclass_to to be 1:2
    for (sizeclass_to in 1:2){
        // equation 12
        g[plot, sizeclass_to, 1] = 
         (normal_cdf(u_bounds[sizeclass_to]| mu_gr[1+(M_r*(plot-1))], sigmaEpsilon_g) - 
          normal_cdf(l_bounds[sizeclass_to]| mu_gr[1+(M_r*(plot-1))], sigmaEpsilon_g)) / 
          (1-normal_cdf(0| mu_gr[1+(M_r*(plot-1))], sigmaEpsilon_g));

      // growth kernel is the product of growth into each size class from 
      // the smallest size class, and 
      // survival in the smallest size class // note: shriver's code has 
      // s[sizeclass_to, plot], which is also unnesessary, because only 
      // growskern[, , sizeclass_from=1] actually gets used
      growKern[plot,sizeclass_to,1] = g[plot, sizeclass_to,1]*s[1,plot];
    }
    
    // loop over all the size classes
    for (sizeclass in 1:M_r){
    
      // expected fecundity in each size class
      f[sizeclass,plot] = exp(logf[sizeclass+(M_r*(plot-1))]);
    
      // recruitment kernel (new recruits per existing adult divvied into the 
      // recruitment size kernel)
      recKern[plot,1:2,sizeclass] = r * f[sizeclass,plot];
    }
    
    // transition kernel is the sum of growth of existing small individuals 
    // plus new recruits from the smallest size class (<1" dbh) (eq 11)
    A[plot,,1] = recKern[plot,1:2,1]+growKern[plot,1:2,1];
    A[plot,,2:20] = recKern[plot,1:2,2:M_r]; // or just recruitment from bigger classes
    
    nprime[plot,] = A[plot,,]*n[plot,]; // eqation 10 in shriver et al; density at t=1 is 
    // the transition kernel matrix multiplied by the density at time t = 0
  }
  
  
  // priors
  beta_s ~ normal(0, 5);
  beta_g ~ normal(0, 5);
  beta_f ~ normal(0, 5);
  sigmaPlot_s ~ normal(0, 5);
  sigmaPlot_g ~ normal(0, 5);
  sigmaPlot_f ~ normal(0, 5);
  sigmaEco_s ~ normal(0, 5);
  sigmaEco_g ~ normal(0, 5);
  sigmaEco_f ~ normal(0, 5);
  sigmaEpsilon_g ~ normal(0, 5);
  kappa_r ~ cauchy(0,5);
  
  // random effect realizations
  zPlot_s ~ std_normal();
  zPlot_g ~ std_normal();
  zPlot_f ~ std_normal();
  zEco_s ~ std_normal();
  zEco_g ~ std_normal();
  zEco_f ~ std_normal();
  
 
  // likelihoods
  surv ~ bernoulli_logit(logitp_s); // survival
  for (i in 1:N_g){
    size1_g[i] ~ normal(mu_g[i], sigmaEpsilon_g) T[0,]; // size at time t+1 (growth)
    }
  for (plot in 1:P_r){
    for (sizeclass in 1:2){
        cprime[sizeclass,plot] ~ 
          neg_binomial_2(nprime[plot,sizeclass]*a[sizeclass], kappa_r);
    }
  }
   
}

