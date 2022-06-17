
data {
  int<lower=1> K; 
  int<lower=1> P;
  int<lower=1> N;
  
  array[N] int<lower=0,upper=1> surv;
  matrix[N,K] X;
  array[N] int<lower=1,upper=P> plot_id;
}

parameters {
  vector[K] beta;
  real<lower=0> sigma_plot;
  vector[P] z_plot;
}

transformed parameters {
  vector[P] effect_plot;
  
  effect_plot = sigma_plot * z_plot;
}

model {
  
  // variable declarations
  vector[N] XB;
  vector[N] logit_p;
  
  // priors
  beta ~ normal(0, 5);
  sigma_plot ~ normal(0, 5);
  
  // random effect realizations
  z_plot ~ std_normal();
  
  // linear predictor
  XB = X * beta;
  for (i in 1:N){
    logit_p[i] = XB[i]  + effect_plot[plot_id[i]];
  }
  
  // likelihood
  surv ~ bernoulli_logit(logit_p);
}