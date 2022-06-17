
data {
  int<lower=1> K; 
  int<lower=1> P;
  int<lower=1> N;
  
  array[N] real<lower=0> size1;
  matrix[N,K] X;
  array[N] int<lower=1,upper=P> plot_id;
}

parameters {
  vector[K] beta;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_epsilon;
  vector[P] z_plot;
}

transformed parameters {
  vector[P] effect_plot;
  
  effect_plot = sigma_plot * z_plot;
}

model {
  
  // variable declarations
  vector[N] XB;
  vector[N] mu;
  
  // priors
  beta ~ normal(0, 5);
  sigma_plot ~ normal(0, 5);
  sigma_epsilon ~ normal(0, 5);
  
  // random effect realizations
  z_plot ~ std_normal();
  
  // linear predictor
  XB = X * beta;
  for (i in 1:N){
    mu[i] = XB[i]  + effect_plot[plot_id[i]];
  }
  
  // likelihood
  for (i in 1:N){
    size1[i] ~ normal(mu[i], sigma_epsilon) T[0,];
  }
}