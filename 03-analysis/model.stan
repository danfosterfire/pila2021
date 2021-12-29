data {
  int N_s; // number of individual trees tracked for survival
  int P_s; // number of unique plots
  int E_s; // number of unique ecoregion subsections
  int surv[N_s]; // 1 if individual i survived from t = 0 to t = 1, 0 otherwise
  int plotid_s[N_s]; // plot ID indices
  int ecosub_s[N_s]; // ecoregion subsection indices
  matrix[N_s,14] X_s; // design matrix for fixed effects, not including size
  
  
  int N_g; // number of individual trees tracked for growth
  int P_g; // number of unique plots in growth data
  int E_g; // number of unique ecoregion subsections
  real size1_g[N_g]; // sizes at time t+1
  
  int plotid_g[N_g]; // plot id indices
  int ecosub_g[N_g]; // ecoregion subsection indices
  matrix[N_g,14] X_g; // design matrix for fixed effects, including intercept
}


parameters {
  vector[14] beta_s; // fixeff coefficients, including intercept, for surv model
  real<lower=0> sigmaPlot_s; // standard deviation of plot random effect
  real<lower=0> sigmaEco_s; // sddev of ecoregion subsection effect
  vector[P_s] zPlot_s; // standard normal deviates for plot random effect
  vector[E_s] zEco_s; // sdnormal deviates for ecoregion random effect
  
  vector[14] beta_g; // fixeff coeffs, including intercept, for growth model
  real<lower=0> sigmaPlot_g; // std deviation of plot random effect
  vector[P_g] zPlot_g; // std normal deviates for plot random effect
  real<lower=0> sigmaEco_g; // std dev of ecoregion random effect
  vector[E_g] zEco_g; // stdnormal deviates
  real<lower=0> sigmaEpsilon_g; // std deviation of residuals for growth model
}

transformed parameters {
  vector[P_s] plotEffect_s;
  vector[P_g] plotEffect_g;
  vector[E_s] ecoEffect_s;
  vector[E_g] ecoEffect_g;
  
  plotEffect_s = sigmaPlot_s * zPlot_s;
  plotEffect_g = sigmaPlot_g * zPlot_g;
  ecoEffect_s = sigmaEco_s * zEco_s;
  ecoEffect_g = sigmaEco_g * zEco_g;
}

model {
  
  // linear predictors
  vector[N_g] mu_g;
  vector[N_s] logitp_s;
  vector[N_s] XB_s;
  vector[N_g] XB_g;
  
  XB_s = X_s * beta_s;
  XB_g = X_g * beta_g;
  
  for (i in 1:N_s){
    logitp_s[i] = XB_s[i] + plotEffect_s[plotid_s[i]] + ecoEffect_s[ecosub_s[i]];
  }
  
  for (i in 1:N_g){
    mu_g[i] = XB_g[i] + plotEffect_g[plotid_g[i]] + ecoEffect_g[ecosub_g[i]];
  }
  // priors
  beta_s ~ normal(0, 5);
  beta_g ~ normal(0, 5);
  sigmaPlot_s ~ normal(0, 5);
  sigmaPlot_g ~ normal(0, 5);
  sigmaEco_s ~ normal(0, 5);
  sigmaEco_g ~ normal(0, 5);
  sigmaEpsilon_g ~ normal(0, 5);
  
  // random effect realizations
  zPlot_s ~ std_normal();
  zPlot_g ~ std_normal();
  zEco_s ~ std_normal();
  zEco_g ~ std_normal();
  
  // likelihood
  surv ~ bernoulli_logit(logitp_s);
  for (i in 1:N_g){
    //print("size: ", size1_g[i]);
    //print("mu_g: ", mu_g[i]);
    //print("lpdf pre", target());
    size1_g[i] ~ normal(mu_g[i], sigmaEpsilon_g) T[0,];
    //print("target post", target());
    };

}

