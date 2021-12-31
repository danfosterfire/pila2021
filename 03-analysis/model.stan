
functions {

// Define probability density function; (just copying shriver's code here, not 
// sure why we need to define a custom pdf function instead of just using 
// the out of the box one?)
vector myNormal_pdf(vector y, real mu, real sigma) {
return (1 / sqrt(2*pi()*sigma^2 )) *  exp(-square(y- mu) / (2 * sigma^2));
}
}

data {
  // surivval data
  int<lower=0> N_s; // number of individual trees tracked for survival
  int<lower=0> P_s; // number of unique plots
  int<lower=0> E_s; // number of unique ecoregion subsections
  int surv[N_s]; // 1 if individual i survived from t = 0 to t = 1, 0 otherwise
  int plotid_s[N_s]; // plot ID indices
  int ecosub_s[N_s]; // ecoregion subsection indices
  matrix[N_s,14] X_s; // design matrix for fixed effects, not including size
  
  // growth data
  int<lower=0> N_g; // number of individual trees tracked for growth
  int<lower=0> P_g;
  int<lower=0> E_g;
  real size1_g[N_g]; // sizes at time t+1
  int plotid_g[N_g]; // plot id indices
  int ecosub_g[N_g]; // ecoregion subsection indices
  matrix[N_g,14] X_g; // design matrix for fixed effects, including intercept
  
  // recruitment data
  int<lower=0> N_r; // M_r (100) times the number of unique plots for the 
  // recruitment submodel
  int<lower=0> P_f; // the number of unique plots for recruitment submodel (plots which 
  // occur in both the growth and survival datasets)
  int<lower=0> S_r; // number of unique subplots for the recruitment submodel
  int<lower=0> E_f; // number of unique ecoregions for the recruitment submodel 
  // (ecoregions which occur in both the growth and survival datasets)
  matrix[N_r, 14] X_r; // fixeff explanatory variables for each sizeclass:subplot 
  // combination, used to predict growth and survival for the recruitment 
  // submodel sizeclasses; instead of observed sizes here we have the 
  // midpoints for each size class; MUST BE ORDERED subplot (slow) sizeclass (fast)
  int plotid_sr[N_r]; // plot id indices for survival random effects in IPM
  int ecosub_sr[N_r]; // ecosub indices for survival random effects in IPM
  int plotid_gr[N_r]; // plot id indices for growth random effects in IPM
  int ecosub_gr[N_r]; // ecosub indices for growth random effects in IPM
  int plotid_fr[N_r]; // plot id indices for fecundity random effects in IPM
  int ecosub_fr[N_r]; // ecosub indices for fecundity random effects in IPM
  int<lower=0> M_r; // number of modeled size classes for the recruitment submodel
  vector[M_r] u_bounds; // upper bounds of size bins
  vector[M_r] l_bounds; //lower bounds of size bins
  vector[M_r] midpoints; // midpoints of size bins
  vector[10] a; // plot area for each of the size classes; here only for the 
  // smallest 10 size classes (only need these, and the macro breakpoint diameter is not consistent 
  // across all the subplots)
  int cprime[M_r,S_r]; // counts of untagged trees in each size class on each subplot
  vector<lower=0>[M_r] n[S_r]; // vector of area-standardaized rates of occurence for 
  // each of size class on each of S_r subplots at time t = 0
}


parameters {
  // survival
  vector[14] beta_s; // fixeff coefficients, including intercept, for surv model
  real<lower=0> sigmaPlot_s; // standard deviation of plot random effect
  real<lower=0> sigmaEco_s; // sddev of ecoregion subsection effect
  vector[P_s] zPlot_s; // standard normal deviates for plot random effect
  vector[E_s] zEco_s; // sdnormal deviates for ecoregion random effect
  
  // growth 
  vector[14] beta_g; // fixeff coeffs, including intercept, for growth model
  real<lower=0> sigmaPlot_g; // std deviation of plot random effect
  vector[P_g] zPlot_g; // std normal deviates for plot random effect
  real<lower=0> sigmaEco_g; // std dev of ecoregion random effect
  vector[E_g] zEco_g; // stdnormal deviates
  real<lower=0> sigmaEpsilon_g; // std deviation of residuals for growth model
  
  // recruitment
  real nu; // mean of recruitment size kernel
  real<lower=0> upsilon; // standard deviation of recruitment size kernel
  real<lower=0> kappa_r; // neg binomial dispersion parameter for recruitment
  vector[14] beta_f; // fixeff coefficients including intercept for fecundity model
  real<lower=0> sigmaPlot_f;
  vector[P_f] zPlot_f;
  real<lower=0> sigmaEco_f;
  vector[E_f] zEco_f;
}

transformed parameters {
  vector[P_s] plotEffect_s;
  vector[P_g] plotEffect_g;
  vector[E_s] ecoEffect_s;
  vector[E_g] ecoEffect_g;
  vector[P_f] plotEffect_f;
  vector[E_f] ecoEffect_f;
  
  plotEffect_s = sigmaPlot_s * zPlot_s;
  plotEffect_g = sigmaPlot_g * zPlot_g;
  ecoEffect_s = sigmaEco_s * zEco_s;
  ecoEffect_g = sigmaEco_g * zEco_g;
  plotEffect_f = sigmaPlot_f * zPlot_f;
  ecoEffect_f = sigmaEco_f * zPlot_f;
}

model {
  
  // local variables
  vector[N_g] mu_g; 
  vector[N_s] logitp_s;
  vector[N_s] XB_s;
  vector[N_g] XB_g;
  vector[N_r] XBs_r; // linear predictor for smallest sizeclass midpoints on each subplot
  vector[N_r] XBg_r; // linear predictor for smallest sizeclass midpoints on each subplot
  vector[N_r] XBf_r;
  vector[N_r] mu_gr;
  vector[N_r] logitp_sr;
  
  vector[10] nprime[S_r]; // area-standardized occurence rates at time t+1 on 
  // each subplot in each size class
  matrix[10,M_r] A[S_r]; // final IPM transition kernel describing how column 
  // size classes at time 0 lead to row size classes at time 1
  matrix[10,M_r] recKern[S_r]; // recruitment kernel 
  matrix[10,M_r] growKern[S_r]; // growth*survival transitions 
  matrix[M_r,M_r] g[S_r]; // growth kernel of transitions from size class to size class on each plot
  matrix[M_r,S_r] s; // survival rates on each subplot for each of the 10 smallest 
  // size classes
  matrix[M_r,S_r] f; // fecundity rates by size class and subplots
  vector[N_r] logf;
  vector[10] r; // vector of recruitment rates into size classes
  // SETTING NU AND UPSILON
  //real nu = ;
  //real upsilon = ;
  
  // fixed effects
  XB_s = X_s * beta_s;
  XB_g = X_g * beta_g;
  XBs_r = X_r * beta_s;
  XBg_r = X_r * beta_g;
  XBf_r = X_r * beta_f;
  
  
  // linear predictor for survival model
  for (i in 1:N_s){
    logitp_s[i] = XB_s[i] + plotEffect_s[plotid_s[i]] + ecoEffect_s[ecosub_s[i]];
  }
  
  // linear predictor for growth (size2) model
  for (i in 1:N_g){
    mu_g[i] = XB_g[i] + plotEffect_g[plotid_g[i]] + ecoEffect_g[ecosub_g[i]];
  }
  
  // linear predictors for gtrowth and survival on the recruitment IPM
  for (i in 1:N_r){
    logitp_sr[i] = XBs_r[i] + plotEffect_s[plotid_sr[i]]+ ecoEffect_s[ecosub_sr[i]];
    mu_gr[i] = XBg_r[i] + plotEffect_g[plotid_gr[i]] + ecoEffect_g[ecosub_gr[i]];
    logf[i] = XBf_r[i] + plotEffect_f[plotid_fr[i]] + ecoEffect_f[ecosub_fr[i]];
  }
  
  // probability of recruits growing into size class h (equation 15)
  r[:] = 
    myNormal_pdf(midpoints[1:10], nu, upsilon) / sum(myNormal_pdf(midpoints[1:10], nu, upsilon));
  //r = exp(normal_lpdf(midpoints[1:10]|nu, upsilon)) / sum(exp(normal_lpdf(midpoints[1:10]|nu, upsilon)));

  // IPM model for recruitment; loop over all the subplots
  for (subplot in 1:S_r){
    
    // loop over all combinations of size classes
    for (sizeclass_1 in 1:M_r){
      for (sizeclass_2 in 1:M_r){
        // equation 12
        g[subplot, sizeclass_2, sizeclass_1] = 
         (normal_cdf(u_bounds[sizeclass_2]| mu_gr[sizeclass_1+(10*(subplot-1))], sigmaEpsilon_g) - 
          normal_cdf(l_bounds[sizeclass_2]| mu_gr[sizeclass_1+(10*(subplot-1))], sigmaEpsilon_g)) / 
          (1-normal_cdf(0| mu_gr[sizeclass_1+(10*(subplot-1))], sigmaEpsilon_g));
      }
    }
    
    // loop over all the size classes
    for (sizeclass in 1:M_r){
    
      // expected fecundity in each size class
      f[sizeclass,subplot] = exp(logf[sizeclass+(M_r*(subplot-1))]);
    
      // expected survival in each size class
      s[sizeclass,subplot] = 
        inv_logit(logitp_sr[sizeclass+(M_r*(subplot-1))]); 
      
      // growth kernel is the product of growth into each size class and 
      // survival in that size class (eq 11.1)
      growKern[subplot,1:10,sizeclass] = g[subplot,1:10,sizeclass] * s[sizeclass,subplot];
      
      // recruitment kernel 
      recKern[subplot,1:10,sizeclass] = r * f[sizeclass,subplot];
    }
    
    // transition kernel is the sum of growth of existing small individuals 
    // plus new recruits for the smallest size class (<1" dbh) (eq 11)
    A[subplot,,1] = recKern[subplot,1:10,1]+growKern[subplot,1:10,1];
    A[subplot,,2:100] = recKern[subplot,1:10,2:M_r]; // or just recruitment for bigger classes
    
    nprime[subplot,] = A[subplot,,]*n[subplot,]; // eqation 10 in shriver et al; density at t=1 is 
    // the transition kernel multiplied by the density at time t = 0
  }
  
  
  // priors
  beta_s ~ normal(0, 5);
  beta_g ~ normal(0, 5);
  sigmaPlot_s ~ normal(0, 5);
  sigmaPlot_g ~ normal(0, 5);
  sigmaEco_s ~ normal(0, 5);
  sigmaEco_g ~ normal(0, 5);
  sigmaEpsilon_g ~ normal(0, 5);
  kappa_r ~ cauchy(0,5);
  nu ~ normal(1,1);
  upsilon~cauchy(0,1);
  beta_f ~ normal(0,5);
  sigmaPlot_f ~ normal(0,5);
  sigmaEco_f ~ normal(0,5);
  
  // random effect realizations
  zPlot_s ~ std_normal();
  zPlot_g ~ std_normal();
  zEco_s ~ std_normal();
  zEco_g ~ std_normal();
  zPlot_f ~ std_normal();
  zEco_f ~ std_normal();
  
  // likelihoods
  surv ~ bernoulli_logit(logitp_s); // survival
  for (i in 1:N_g){
    size1_g[i] ~ normal(mu_g[i], sigmaEpsilon_g) T[0,]; // size at time t+1 (growth)
    }
  for (subplot in 1:S_r){
    //print("************Starting subplot: ", subplot);
    for (sizeclass in 1:10){
      //print("nprime: ", nprime[subplot,sizeclass]);
      //print("kappa_r: ", kappa_r);
      //if (is_nan(nprime[subplot,sizeclass]))
        //print("Found nan nprime; A: ", A[subplot,,sizeclass], "********************** n: ", n[subplot,sizeclass]);
      //  {print("Found nan nprime; sizeclass: ", sizeclass, 
      //        "   subplot: ", subplot, 
      //        "    recKern: ", recKern[subplot,1:10,sizeclass], 
      //        "     growKern: ", growKern[subplot,1:10, sizeclass]);
      //  print("found nan nprime; r: ", r, "      f[sizeclass,subplot]: ", f[sizeclass,subplot]);
      //  print("nu: ", nu, "   upsilon: ", upsilon, "   midpoints: ",midpoints[1:10], "   r: ", r);} 
        //print("    recKern: ", recKern[subplot,1:10,sizeclass]);
        //{print(nprime[subplot,sizeclass]);
        //print("A: ", A);
        //print("n: ", n);}
      //else if (nprime[subplot,sizeclass]==0)
      //  {print(nprime[subplot,sizeclass]);
      //  print("A: ", A);
      //  print("n: ", n);}
      //else
        cprime[sizeclass,subplot] ~ 
          neg_binomial_2(nprime[subplot,sizeclass]*a[sizeclass], kappa_r);
    }
  }
   
}

