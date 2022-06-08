
data {
  int<lower=1> K; // number of fixeff params
  //int<lower=1> P; // number of plots
  int<lower=1> N; // number of plot:sizeclass:census combinations
  int<lower=1> M; // number of unique size classes
  int<lower=1> C; // number of unique censuses
  
  array[C] int<lower=0> cprime; // response; observed tally of new recrs
  matrix[N,K] X; // fixeff covariates for plot:sizeclass combos
  //array[N] int<lower=1,upper=P> plot_id; // plot indices
  matrix<lower=0>[M,C] n; // observed TPA of existing individuals for each size,census
  vector[C] a; // plot area searched for new recruits in hectares
  
  
}

parameters {
  vector[K] beta;
  //real<lower=0> sigma_plot;
  real<lower=0> kappa;
  //vector[P] z_plot;
}

transformed parameters {
  //vector[P] effect_plot;
  
  //effect_plot = sigma_plot * z_plot;
}

model {
  
  // variable declarations
  vector[N] XB; // fixed effects times coeffs
  vector[N] logf; // linear predictor for fecundity
  vector[C] nprime; // predicted density of new recruits
  
  // priors
  beta ~ normal(0, 5);
  //sigma_plot ~ normal(0, 5);
  kappa ~ cauchy(0, 5);
  
  // random effect realizations
  //z_plot ~ std_normal();
  
  // linear predictor
  XB = X * beta;
  for (i in 1:N){
    //logf[i] = XB[i] + effect_plot[plot_id[i]];
    logf[i] = XB[i];
  }
  
  for (i in 1:C){
    nprime[i] = sum(exp(logf[ (1+(M*(i-1))):(M+(M*(i-1))) ] ));
  }
  
  // likelihood
  for (i in 1:C){
    cprime[i] ~ neg_binomial_2(nprime[i]*a[i], kappa);
  }
}