
data {
  int<lower=1> K; // number of fixeff covariates, including intercept
  
  // survival data
  int<lower=0> N_s; // number of responses for survival data
  matrix[N_s,K] X_s; // fixeff covariates for survival data
  int<lower=0, upper=1> surv[N_s]; // observed survival
  
  // growth data
  int<lower=0> N_g; // number of responses for growth data
  matrix[N_g,K] X_g; // fixeff covariates for growth data
  real<lower=0> size1_g[N_g]; // size
  
}


parameters {
  vector[K] beta_s;
  vector[K] beta_g;
  real<lower=0> sigmaEpsilon_g;
}

transformed parameters {
  
}

model {
  // variable declarations
  vector[N_s] XB_s;
  vector[N_g] XB_g;
  
  // priors
  beta_s ~ normal(0,1);
  beta_g[1] ~ normal(0.15, 0.1); // special prior for the intercept, which should be mostly positive vals near 0
  beta_g[2] ~ normal(1, 0.25);// special prior for the effect of initial size, which should be near 1
  for (i in 3:K){
    beta_g[i] ~ normal(0, 0.1);
  }
  
  sigmaEpsilon_g ~ normal(0,0.25);
  
  // random effect realizations
  
  // linear predictors
  XB_s = X_s * beta_s;
  XB_g = X_g * beta_g;
  
  // likelihoods
  surv ~ bernoulli_logit(XB_s);
  for (i in 1:N_g){
    size1_g[i] ~ normal(XB_g[i], sigmaEpsilon_g)T[0,];
  }
}
