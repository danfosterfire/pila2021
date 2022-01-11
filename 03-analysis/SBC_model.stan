functions {
  // from the stan user's guide
  real normal_lb_rng(real mu, real sigma, real lb) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = (p_lb < 1) ? uniform_rng(p_lb, 1) : 1 ;
    real y = mu + sigma * inv_Phi(u);
    return y;
    }
}

data {
  int<lower=1> K; // number of fixeff covariates, including intercept
  
  // survival data
  int<lower=0> N_s; // number of responses for survival data
  matrix[N_s,K] X_s; // fixeff covariates for survival data
  
  // growth data
  int<lower=0> N_g; // number of responses for growth data
  matrix[N_g,K] X_g; // fixeff covariates for growth data
  
}

transformed data {
  
  // declare simulation parameters and responses
  vector[K] true_beta_s;
  vector[K] true_beta_g;
  real<lower=0> true_sigmaEpsilon_g = normal_lb_rng(0, 5, 0);
  
  vector[N_s] true_XB_s;
  vector[N_g] true_XB_g;
  
  int<lower=0, upper=1> sim_surv[N_s];
  real<lower=0> sim_size1_g[N_g];
  
  // draw simulation parameter values from the prior distributions
  for (i in 1:K){
    true_beta_s[i] = normal_rng(0,5);
    true_beta_g[i] = normal_rng(0,5);
  }
  
  // assemble linear predictor terms
  true_XB_s = X_s * true_beta_s;
  true_XB_g = X_g * true_beta_g;
  
  // simulate responses
  for (i in 1:N_s){
    sim_surv[i] = bernoulli_logit_rng(true_XB_s[i]);
  }
  for (i in 1:N_g){
    sim_size1_g[i] = normal_lb_rng(true_XB_g[i], true_sigmaEpsilon_g, 0);
  }

}

parameters {
  vector[K] hat_beta_s;
  vector[K] hat_beta_g;
  real<lower=0> hat_sigmaEpsilon_g;
}

transformed parameters {
  
}

model {
  // variable declarations
  vector[N_s] hat_XB_s;
  vector[N_g] hat_XB_g;
  
  // priors
  hat_beta_s ~ normal(0,5);
  hat_beta_g ~ normal(0,5);
  hat_sigmaEpsilon_g ~ normal(0,5);
  
  // random effect realizations
  
  // linear predictors
  hat_XB_s = X_s * hat_beta_s;
  hat_XB_g = X_g * hat_beta_g;
  
  // likelihoods
  sim_surv ~ bernoulli_logit(hat_XB_s);
  for (i in 1:N_g){
    sim_size1_g[i] ~ normal(hat_XB_g[i], hat_sigmaEpsilon_g);
  }
}

generated quantities {
  int<lower=0, upper=1> IltTrue_beta_s[K];
  int<lower=0, upper=1> IltTrue_beta_g[K];
  int<lower=0, upper=1> IltTrue_sigmaEpsilon_g;
    
  for (i in 1:K){
    IltTrue_beta_s[i] = hat_beta_s[i] < true_beta_s[i];
    IltTrue_beta_g[i] = hat_beta_g[i] < true_beta_g[i];
  }
  
  IltTrue_sigmaEpsilon_g = hat_sigmaEpsilon_g < true_sigmaEpsilon_g;
    
}