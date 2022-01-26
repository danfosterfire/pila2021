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
}

parameters {
  // survival
  vector[K] beta_s; // fixed effects coefficients 
  real<lower=0> sigmaPlot_s; // random effect SD
  vector[P] zPlot_s; // standard normal deviates for random effect
}

transformed parameters {
  vector[P] plotEffect_s;
  
  plotEffect_s = sigmaPlot_s * zPlot_s;
  
}

model {
  // variable declarations
  // // survival
  vector[N_s] XB_s;
  vector[N_s] logitp_s;
  
  // linear predictors
  // // survival
  XB_s = X_s * beta_s;
  for (i in 1:N_s){
    logitp_s[i] = XB_s[i] + plotEffect_s[plotid_s[i]];
  }
  
  // priors
  beta_s ~ normal(0, 5);
  sigmaPlot_s ~ normal(0, 5);
  
  // random effect realizations
  zPlot_s ~ std_normal();
  
  // likelihoods
  // // survival
  surv_s ~ bernoulli_logit(logitp_s);
  
}