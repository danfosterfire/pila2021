
data {
  // surivval data
  int<lower=0> N_s; // number of individual trees tracked for survival
  //int<lower=0> P_s; // number of unique plots
  //int<lower=0> E_s; // number of unique ecoregion subsections
  int surv[N_s]; // 1 if individual i survived from t = 0 to t = 1, 0 otherwise
  //int plotid_s[N_s]; // plot ID indices
  //int ecosub_s[N_s]; // ecoregion subsection indices
  matrix[N_s,14] X_s; // design matrix for fixed effects, not including size
  
  // growth data
  int<lower=0> N_g; // number of individual trees tracked for growth
  //int<lower=0> P_g;
  //int<lower=0> E_g;
  real size1_g[N_g]; // sizes at time t+1
  //int plotid_g[N_g]; // plot id indices
  //int ecosub_g[N_g]; // ecoregion subsection indices
  matrix[N_g,14] X_g; // design matrix for fixed effects, including intercept
  
  // recruitment data
  int<lower=0> N_r; // M_r (20) times the number of unique plots for the 
  // recruitment submodel
  //int<lower=0> P_f; // the number of unique plots for recruitment submodel (plots which 
  // occur in both the growth and survival datasets)
  int<lower=0> S_r; // number of unique subplots for the recruitment submodel
  //int<lower=0> E_f; // number of unique ecoregions for the recruitment submodel 
  // (ecoregions which occur in both the growth and survival datasets)
  matrix[N_r, 14] X_r; // fixeff explanatory variables for each sizeclass:subplot 
  // combination, used to predict growth and survival for the recruitment 
  // submodel sizeclasses; instead of observed sizes here we have the 
  // midpoints for each size class; MUST BE ORDERED subplot (slow) sizeclass (fast)
  //int plotid_sr[N_r]; // plot id indices for survival random effects in IPM
  //int ecosub_sr[N_r]; // ecosub indices for survival random effects in IPM
  //int plotid_gr[N_r]; // plot id indices for growth random effects in IPM
  //int ecosub_gr[N_r]; // ecosub indices for growth random effects in IPM
  //int plotid_fr[N_r]; // plot id indices for fecundity random effects in IPM
  //int ecosub_fr[N_r]; // ecosub indices for fecundity random effects in IPM
  int<lower=0> M_r; // number of modeled size classes for the recruitment submodel
  vector[M_r] u_bounds; // upper bounds of size bins
  vector[M_r] l_bounds; //lower bounds of size bins
  vector[M_r] midpoints; // midpoints of size bins
  vector[2] a; // plot area for each of the size classes; here only for the 
  // smallest 2 size classes (only need these, and the macro breakpoint diameter is not consistent 
  // across all the subplots)
  int cprime[M_r,S_r]; // counts of untagged trees in each size class on each subplot
  vector<lower=0>[M_r] n[S_r]; // vector of area-standardaized rates of occurence for 
  // each of size class on each of S_r subplots at time t = 0
  vector<lower=0>[2] r; 
}


parameters {
  // survival
  vector[14] beta_s; // fixeff coefficients, including intercept, for surv model
  
  // growth 
  vector[14] beta_g; // fixeff coeffs, including intercept, for growth model
  real<lower=0> sigmaEpsilon_g; // std deviation of residuals for growth model
  
  // recruitment
  real<lower=0> nu; // mean of recruitment size kernel
  real<lower=0> upsilon; // standard deviation of recruitment size kernel
  real<lower=0> kappa_r; // neg binomial dispersion parameter for recruitment
  vector[14] beta_f; // fixeff coefficients including intercept for fecundity model
}


model {
  
  // local variables
  vector[N_g] mu_g;  // linear predictor for mean size at time 1 (growth data)
  vector[N_s] logitp_s; // linear predictor for survival prob. (surv data)
  vector[N_r] mu_gr; // linear predictor for mean size at time 1 (recr data)
  vector[N_r] logitp_sr; // linear predictor for survival prob (recr data)
  
  vector[2] nprime[S_r]; // area-standardized occurence rates at time t+1 on 
  // each subplot in each size class
  matrix[2,M_r] A[S_r]; // final IPM transition kernel describing how column 
  // size classes at time 0 lead to row size classes at time 1
  matrix[2,M_r] recKern[S_r]; // recruitment kernel 
  matrix[2,M_r] growKern[S_r]; // growth*survival transitions 
  matrix[M_r,M_r] g[S_r]; // growth kernel of transitions from size class to size class on each plot
  matrix[M_r,S_r] s; // survival rates on each subplot for each of the 2 smallest 
  // size classes
  matrix[M_r,S_r] f; // fecundity rates by size class and subplots
  vector[N_r] logf; // linear predictor of fecundity rates
  
  // fixed effects
  logitp_s = X_s * beta_s;
  mu_g = X_g * beta_g;
  logitp_sr = X_r * beta_s;
  mu_gr = X_r * beta_g;
  logf = X_r * beta_f;
  
  
  //print("r(2): ", r);
  // IPM model for recruitment; loop over all the subplots
  for (subplot in 1:S_r){
    
    // loop over all combinations of size classes
    for (sizeclass_1 in 1:M_r){
      for (sizeclass_2 in 1:M_r){
        // equation 12
        g[subplot, sizeclass_2, sizeclass_1] = 
         (normal_cdf(u_bounds[sizeclass_2]| mu_gr[sizeclass_1+(M_r*(subplot-1))], sigmaEpsilon_g) - 
          normal_cdf(l_bounds[sizeclass_2]| mu_gr[sizeclass_1+(M_r*(subplot-1))], sigmaEpsilon_g)) / 
          (1-normal_cdf(0| mu_gr[sizeclass_1+(M_r*(subplot-1))], sigmaEpsilon_g));
      }
    }
    
    // loop over all the size classes
    for (sizeclass in 1:M_r){
    
      // expected fecundity in each size class
      f[sizeclass,subplot] = exp(logf[sizeclass+(M_r*(subplot-1))]);
    
      // expected survival in each size class
      s[sizeclass,subplot] = 
        inv_logit(logitp_sr[sizeclass+(M_r*(subplot-1))]); 
      
      // growth kernel is the product of growth into each size class and 
      // survival in that size class (eq 11.1)
      growKern[subplot,1:2,sizeclass] = g[subplot,1:2,sizeclass] * s[sizeclass,subplot];
      
      // recruitment kernel 
      recKern[subplot,1:2,sizeclass] = r * f[sizeclass,subplot];
    }
    
    // transition kernel is the sum of growth of existing small individuals 
    // plus new recruits for the smallest size class (<1" dbh) (eq 11)
    A[subplot,,1] = recKern[subplot,1:2,1]+growKern[subplot,1:2,1];
    A[subplot,,2:M_r] = recKern[subplot,1:2,2:M_r]; // or just recruitment for bigger classes
    
    nprime[subplot,] = A[subplot,,]*n[subplot,]; // eqation 10 in shriver et al; density at t=1 is 
    // the transition kernel multiplied by the density at time t = 0
  }
  
  
  // priors
  beta_s ~ normal(0, 5);
  beta_g ~ normal(0, 5);
  sigmaEpsilon_g ~ normal(0, 5);
  kappa_r ~ cauchy(0,5);
  nu ~ normal(0,5);
  upsilon~cauchy(0,1);
  beta_f ~ normal(0,5);
  
 
  // likelihoods
  surv ~ bernoulli_logit(logitp_s); // survival
  for (i in 1:N_g){
    size1_g[i] ~ normal(mu_g[i], sigmaEpsilon_g) T[0,]; // size at time t+1 (growth)
    }
  for (subplot in 1:S_r){
    for (sizeclass in 1:2){
        cprime[sizeclass,subplot] ~ 
          neg_binomial_2(nprime[subplot,sizeclass]*a[sizeclass], kappa_r);
    }
  }
   
}

