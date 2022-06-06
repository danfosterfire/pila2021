// model for delta TPH, using a cauchy response
data {
  int<lower=0> N;
  int<lower=0> E;
  vector[N] Y;
  matrix[N,6] X;
  array[N] int ecosub_id;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[6] beta;
  real<lower=0> sigma_epsilon;
  real<lower = 0> sigma_ecosub;
  vector[E] z_ecosub;
}

transformed parameters {
  vector[E] ecosub_effects;
  
  ecosub_effects = sigma_ecosub * z_ecosub;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // variable declarations
  vector[N] XB;
  
  // priors
  beta ~ normal(0, 10);
  sigma_ecosub ~ normal(0, 10);
  sigma_epsilon ~ normal(0, 10);
  
  // random effect realizations
  z_ecosub ~ std_normal();
  
  // linear predictor
  XB = X * beta;
  
  for (i in 1:N){
    Y[i] ~ normal(XB[i]+ecosub_effects[ecosub_id[i]], sigma_epsilon);
  }
}

