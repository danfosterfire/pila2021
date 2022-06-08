
data {
  int<lower=1> K; // number of fixeff params
  int<lower=1> E; // number of ecoregions
  int<lower=1> P; // number of plots
  int<lower=1> N; // number of plot:sizeclass combinations
  int<lower=1> M; // number of unique size classes
  
  array[P] int<lower=0> cprime; // response; observed tally of new recrs
  matrix[N,K] X; // fixeff covariates for plot:sizeclass combos
  array[N] int<lower=1,upper=E> ecosub_id; // ecoregion indices
  matrix<lower=0>[M,P] n; // observed TPA of existing individuals
  real<lower=0> a; // plot area searched for new recruits in acres
  
  
}

parameters {
  vector[K] beta;
  real<lower=0> sigma_ecosub;
  real<lower=0> kappa;
  vector[E] z_ecosub;
}

transformed parameters {
  vector[E] effect_ecosub;
  
  effect_ecosub = sigma_ecosub * z_ecosub;
}

model {
  
  // variable declarations
  vector[N] XB; // fixed effects times coeffs
  vector[N] logf; // linear predictor for fecundity
  vector[P] nprime; // predicted density of new recruits
  
  // priors
  beta ~ normal(0, 5);
  sigma_ecosub ~ normal(0, 5);
  kappa ~ cauchy(0, 5);
  
  // random effect realizations
  z_ecosub ~ std_normal();
  
  // linear predictor
  XB = X * beta;
  for (i in 1:N){
    logf[i] = XB[i] + effect_ecosub[ecosub_id[i]];
  }
  
  for (i in 1:P){
    nprime[i] = sum(exp(logf[ (1+(M*(i-1))):(M+(M*(i-1))) ] ));
  }
  
  // likelihood
  cprime ~ neg_binomial_2(nprime*a, kappa);
}