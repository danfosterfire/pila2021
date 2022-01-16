
data {
  int<lower=1> K; // number of fixeff covariates, including intercept
  int<lower=1> E; // number of ecoregions
  int<lower=1> P; // number of plots
  
  // survival data
  int<lower=0> N_s; // number of responses for survival data
  matrix[N_s,K] X_s; // fixeff covariates for survival data
  int<lower=0, upper=1> surv[N_s]; // observed survival
  int<lower=1, upper = E> ecosub_s[N_s]; // ecoregion indices
  int<lower=1, upper = P> plotid_s[N_s]; // plot indices
  
  // growth data
  int<lower=0> N_g; // number of responses for growth data
  matrix[N_g,K] X_g; // fixeff covariates for growth data
  real<lower=0> size1_g[N_g]; // size
  int<lower=1, upper = E> ecosub_g[N_g];
  int<lower=1, upper = P> plotid_g[N_g]; // plot indices
}


parameters {
  // survival
  vector[K] beta_s;
  real<lower=0> sigmaEco_s;
  vector[E] zEco_s;
  real<lower=0> sigmaPlot_s;
  vector[P] zPlot_s;
  
  // growth 
  vector[K] beta_g;
  real<lower=0> sigmaEpsilon_g;
  real<lower=0> sigmaEco_g;
  vector[E] zEco_g;
  real<lower=0> sigmaPlot_g;
  vector[P] zPlot_g;
}

transformed parameters {
  vector[E] ecoEffect_s;
  vector[E] ecoEffect_g;
  vector[P] plotEffect_s;
  vector[P] plotEffect_g;
  
  ecoEffect_s = sigmaEco_s * zEco_s;
  ecoEffect_g = sigmaEco_g * zEco_g;
  plotEffect_s = sigmaPlot_s*zPlot_s;
  plotEffect_g = sigmaPlot_g * zPlot_g;
}

model {
  // variable declarations
  // // survival
  vector[N_s] XB_s;
  // // growth
  vector[N_g] XB_g;
  
  // priors
  // // survival
  beta_s ~ normal(0,1);
  sigmaEco_s ~ normal(0,1);
  sigmaPlot_s ~ normal(0,1);
  // // growth
  beta_g[1] ~ normal(0.15, 0.1); // special prior for the intercept, which should be mostly positive vals near 0
  beta_g[2] ~ normal(1, 0.25);// special prior for the effect of initial size, which should be near 1
  for (i in 3:K){
    beta_g[i] ~ normal(0, 0.1);
  }
  sigmaEpsilon_g ~ normal(0,0.25);
  sigmaEco_g ~ normal(0, 0.25); 
  sigmaPlot_g ~ normal(0, 0.25);
  // random effect realizations
  // // survival
  zEco_s ~ std_normal();
  zPlot_s ~ std_normal();
  // // growth
  zEco_g ~ std_normal();
  zPlot_g ~ std_normal();
  
  // linear predictors
  // // survival 
  XB_s = X_s * beta_s;
  // // growth
  XB_g = X_g * beta_g;
  
  // likelihoods
  // // survival 
  for (i in 1:N_s){
    surv[i] ~ bernoulli_logit(XB_s[i]+ecoEffect_s[ecosub_s[i]]+plotEffect_s[plotid_s[i]]);
    //surv[i] ~ bernoulli_logit(XB_s[i]+ecoEffect_s[ecosub_s[i]]);
  }
  // // growth
  for (i in 1:N_g){
    size1_g[i] ~ normal(XB_g[i]+ecoEffect_g[ecosub_g[i]]+plotEffect_g[plotid_g[i]], 
                        sigmaEpsilon_g)T[0,];
  }
}
