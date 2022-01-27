data {
  int<lower=1> K; // number of fixed effect parameters
  int<lower=1> P; // number of unique plots
  
  // survival data
  int<lower=0> N_s; // number of individual observations (tree:census) in the 
                    // survival model
  int<lower=0, upper=1> surv_s[N_s]; // observed status at time T for each observation,
                                // 1 if alive, 0 if dead
  matrix[N_s, K] X_s; // fixed effects covariates matrix for survival model
  int<lower=1, upper = P> plotid_s[N_s]; // plot indices for survival observations
  
  // growth data
  int<lower=0> N_g; // number of individual observations (tree:census) in the 
                    // growth model
  real<lower=0> dbhT_g[N_g]; // observed size at time T for each observation in meters
  matrix[N_g, K] X_g; // fixed effects covariates matrix for growth model
  int<lower=1, upper=P> plotid_g[N_g]; // plot indices for growth observations
  
  // fecundity data
  int<lower=0> N_f; // number of plot:year:sizeclass combinations
  int<lower=0> C_f; // number of unique censuses (plot:year)
  int<lower=0> M_f; // number of unique size classes
  int<lower=0> cprime_f[C_f]; // observed tally of new recruits
  matrix<lower=0>[M_f,C_f] n; // observed tph of existing individuals 
  matrix[N_f,K] X_f; // fixed effects covariates for fecundity
  int<lower=1,upper=P> plotid_f[N_f]; // plot indices
  vector[C_f] a; // plot areas for each census
  
  // recruit size data
  int<lower=0> N_r; // number of individual new recruits
  real logdbh_r[N_r]; // log(dbh) of new recruits
}


parameters {
  // survival
  vector[K] beta_s; // fixed effects coefficients 
  real<lower=0> sigmaPlot_s; // random effect SD
  vector[P] zPlot_s; // standard normal deviates for random effectrepeat
  
  // growth
  vector[K] beta_g; // fixed effects coeffs
  real<lower=0> sigmaPlot_g; // random effect SD
  vector[P] zPlot_g; // std normal deviates for random effect
  real<lower=0> sigmaEpsilon_g; // residual sd
  
  // fecundity
  vector[K] beta_f; // fixed effects coeffs
  real<lower=0> sigmaPlot_f; // random effect sd
  vector[P] zPlot_f; // std normal deviates for raneff
  real<lower=0> kappa_f; // neg binomial dispersion
  
  // recruitment
  real mu_r; // mean of log(dbh) of new recruits
  real<lower=0> sigmaEpsilon_r; // sd of log(dbh) of new recruits
}

transformed parameters {
  vector[P] plotEffect_s;
  vector[P] plotEffect_g;
  vector[P] plotEffect_f;
  
  plotEffect_s = sigmaPlot_s * zPlot_s;
  plotEffect_g = sigmaPlot_g * zPlot_g;
  plotEffect_f = sigmaPlot_f * zPlot_f;
}

model {
  // variable declarations
  // // survival
  vector[N_s] XB_s;
  vector[N_s] logitp_s;
  // // growth
  vector[N_g] XB_g;
  vector[N_g] mu_g;
  // // fecundity
  vector[N_f] XB_f;
  vector[N_f] logf;
  vector[C_f] nprime; // area standardized occurence rates of new recruits
  matrix[C_f,M_f] A; // transition "matrices" (1xM) for fecundity to recruits, one per census 
  
  // linear predictors
  // // survival
  XB_s = X_s * beta_s;
  for (i in 1:N_s){
    logitp_s[i] = XB_s[i] + plotEffect_s[plotid_s[i]];
  }
  // // growth
  XB_g = X_g * beta_g;
  for (i in 1:N_g){
    mu_g[i] = XB_g[i]+plotEffect_g[plotid_g[i]];
  }
  // // fecundity
  XB_f = X_f * beta_f;
  for (i in 1:N_f){
    logf[i] = XB_f[i]+plotEffect_f[plotid_f[i]];
  }
  for (census in 1:C_f){
    for (sizeclass_from in 1:M_f){
      A[census,sizeclass_from] = exp(logf[sizeclass_from+(M_f*(census-1))]);
    }
    
    nprime[census] = A[census,]*n[,census];
  }
  
  
  // priors
  beta_s ~ normal(0, 5);
  beta_g ~ normal(0,5);
  beta_f ~ normal(0,5);
  sigmaPlot_s ~ normal(0, 5);
  sigmaPlot_g ~ normal(0,5);
  sigmaPlot_f ~ normal(0,5);
  sigmaEpsilon_g ~ normal(0,5);
  kappa_f ~ cauchy(0,5);
  
  // random effect realizations
  zPlot_s ~ std_normal();
  zPlot_g ~ std_normal();
  zPlot_f ~ std_normal();
  
  // likelihoods
  // // survival
  surv_s ~ bernoulli_logit(logitp_s);
  // // growth
  for (i in 1:N_g){
    dbhT_g[i] ~ normal(mu_g[i], sigmaEpsilon_g) T[0,];
  }
  // // fecundity
  for (i in 1:C_f){
    cprime_f[i] ~ neg_binomial_2(nprime[i]*a[i], kappa_f);
  }
  // // recruitment
  logdbh_r ~ normal(mu_r, sigmaEpsilon_r);
}