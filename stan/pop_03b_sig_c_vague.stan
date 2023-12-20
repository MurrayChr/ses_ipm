// Population model with multivariate-normal prior on logit-scale vital rates
// and Poisson-lognormal 'Random Effects' prior on immigration.
// First year 1986 (not 1983 as in the capture-recapture data and hmm models )

// Here the count error variance sigma_c is estimated with a vague 
// normal(0,20) prior (cf. fixed to ten in pop_01_logit_mvn.stan)

functions {
  // moment-match lognormal to binomial
  row_vector get_lnorm_pars_binom(real log_N, real p){
    real mu;
    real sig2;
    row_vector[2] pars;
    if ((p>1)||(p<0)){
      reject("p must be in [0,1], but p = ",p);
    }
    else {
      mu = 1.5*(log_N + log(p)) - 0.5*log_sum_exp(log_N + log(p), log(1-p));
      sig2 = log_sum_exp(log_N + log(p), log(1-p)) - log_N - log(p);
    }
    pars = [mu, sqrt(sig2)];
    return pars;
  }
  // moment-match lognormal to poisson
  row_vector get_lnorm_pars_pois(real lambda){
    real mu;
    real sig2;
    row_vector[2] pars;
    if (0 >= lambda) {
      reject("lambda must be positive, but lambda = ",lambda);
    }
    else {
      mu = 1.5*log(lambda) - 0.5*log(lambda + 1);
      sig2 = log(lambda + 1) - log(lambda);
    }
    pars = [mu, sqrt(sig2)];
    return pars;
  }
}

data {
  int<lower=2> T;                         // number of years
  vector<lower=0>[T] y_count;             // count data
  
  // parameters for vital rate priors
  int<lower=2> n_pars;                       // number of vital rate parameters
  vector[n_pars] mean_logit_vr;            // means of logit-scale vital rate parameters
  matrix[n_pars, n_pars] cov_logit_vr;       // covariance matrix of logit-scale vital rates
  // indices
  array[T-1] int inds_s0;
  array[T-1] int inds_sN;
  array[T-1] int inds_sB;
  array[T-1] int inds_f3;
  array[T-1] int inds_f4;
  array[T-1] int inds_bb;
  array[T-1] int inds_nb;
  
  // mean and sd parameters for initial stage size priors
  array[2] real<lower=0> Pb1_1_pars;
  array[2] real<lower=0> Pb2_1_pars;
  array[2] real<lower=0> Pb3_1_pars;
  array[2] real<lower=0> Pb4_1_pars;
  array[2] real<lower=0> Br_1_pars;
  array[2] real<lower=0> Nb_1_pars;
}

transformed data {
  matrix[n_pars, n_pars] L = cholesky_decompose(cov_logit_vr);
}

parameters {
  // logit-scale vital rates
  vector[7*(T-1)] logit_vr;

  // survival to weaning probability
  vector<lower=0, upper=1>[T-1] sW;

  // initial log pop stage sizes
  real log_Pb1_1;
  real log_Pb2_1;
  real log_Pb3_1;
  real log_Pb4_1;
  real log_Br_1;
  real log_Nb_1;
  
  // raw log pop stage sizes for Pb1, ..., Pb4, In
  // for non-centered parametrisation
  vector[T-1] z_log_Pb1;
  vector[T-1] z_log_Pb2;       
  vector[T-1] z_log_Pb3;      
  vector[T-1] z_log_Pb4;
  vector[T] z_log_In;
  
  // raw log summands for Br and Nb stage sizes
  // for non-centered parmetrisation
  vector[T-1] z_log_Pb2Br;    
  vector[T-1] z_log_Pb3Br;
  vector[T-1] z_log_Pb4Br;
  vector[T-1] z_log_BrBr;
  vector[T-1] z_log_InBr;
  vector[T-1] z_log_NbBr;
  vector[T-1] z_log_BrNb;
  vector[T-1] z_log_InNb;
  vector[T-1] z_log_NbNb;
  
  // 'poisson-lognormal' prior parameters on immigration
  vector[T] z_log_lambda;      // raw "poisson parameter" lambda on log scale
  real mean_log_lambda;          // mean for hyperprior on poisson parameter 
  real<lower=0> sd_log_lambda;   // std.dev for hyperprior on poisson parameter
  
  // count error variance
  real<lower=0> sig_c;
}

transformed parameters {
  // survival probabilities
  vector<lower=0,upper=1>[T-1] s0 = inv_logit(logit_vr[inds_s0]);
  vector<lower=0,upper=1>[T-1] sN = inv_logit(logit_vr[inds_sN]);
  vector<lower=0,upper=1>[T-1] sB = inv_logit(logit_vr[inds_sB]);
  // conditional state transition probabilities
  vector<lower=0,upper=1>[T-1] f3 = inv_logit(logit_vr[inds_f3]);    
  vector<lower=0,upper=1>[T-1] f4 = inv_logit(logit_vr[inds_f4]);   
  vector<lower=0,upper=1>[T-1] nb = inv_logit(logit_vr[inds_nb]);     
  vector<lower=0,upper=1>[T-1] bb = inv_logit(logit_vr[inds_bb]);     
  
  //// declare variables for the population model
  
  // expected number of immigrants in each year
  vector<lower=0>[T] lambda;
  
  // log pop stage sizes 
  vector[T] log_Pb1;         // one-year-old prebreeders
  vector[T] log_Pb2;         // two-year-old prebreeders
  vector[T] log_Pb3;         // three-year-old prebreeders   
  vector[T] log_Pb4;         // four-year-old prebreeders
  vector[T] log_Br;          // all breeders (excluding new immigrants)
  vector[T] log_In;          // newly arrived breeding immigrants
  vector[T] log_Nb;          // all nonbreeders 
  
  // log summands for Br and Nb stage sizes
  vector[T-1] log_Pb2Br;    // three-year-old first-time breeders
  vector[T-1] log_Pb3Br;    // four-year-old first-time breeders
  vector[T-1] log_Pb4Br;    // five-year-old first-time breeders
  vector[T-1] log_BrBr;     // surviving breeders transition to breeders
  vector[T-1] log_InBr;     // surviving immigrants transition to breeders
  vector[T-1] log_NbBr;     // surviving nonbreeders transition to breeders
  vector[T-1] log_BrNb;     // surviving breeders transition to nonbreeders
  vector[T-1] log_InNb;     // surviving immigrants transition to nonbreeders
  vector[T-1] log_NbNb;     // surviving nonbreeders transition to nonbreeders
    
  // lognormal mean and sd parameters for Pb1, ..., Pb4, In
  array[T-1] row_vector[2] pars_log_Pb1;
  array[T-1] row_vector[2] pars_log_Pb2;
  array[T-1] row_vector[2] pars_log_Pb3;
  array[T-1] row_vector[2] pars_log_Pb4;
  array[T] row_vector[2] pars_log_In;
  
  // lognormal mean and sd parameters for Br and Nb summands
  array[T-1] row_vector[2] pars_log_Pb2Br;
  array[T-1] row_vector[2] pars_log_Pb3Br;
  array[T-1] row_vector[2] pars_log_Pb4Br;
  array[T-1] row_vector[2] pars_log_InBr;
  array[T-1] row_vector[2] pars_log_BrBr;
  array[T-1] row_vector[2] pars_log_NbBr;
  array[T-1] row_vector[2] pars_log_InNb;
  array[T-1] row_vector[2] pars_log_BrNb;
  array[T-1] row_vector[2] pars_log_NbNb;
  
  //// define the population model
  
  // immigrants at all times
  for (t in 1:T){
    // expected number per year
    lambda[t] = exp( mean_log_lambda + sd_log_lambda * z_log_lambda[t] );
    // moment-match lognormal to poisson
    pars_log_In[t] = get_lnorm_pars_pois( lambda[t] );         
    // log number immigrants per year
    log_In[t]  =  pars_log_In[t][1] +  pars_log_In[t][2] *  z_log_In[t];      
  }
  
  // initial log sizes
  log_Pb1[1] = log_Pb1_1;
  log_Pb2[1] = log_Pb2_1;
  log_Pb3[1] = log_Pb3_1;      
  log_Pb4[1] = log_Pb4_1;
  log_Br[1] = log_Br_1;
  log_Nb[1] = log_Nb_1;
  
  for (t in 1:(T-1)) {
    // compute mean and sd parameters
    //  mean and sd parameters for Pb1, ..., Pb4, In
    pars_log_Pb1[t] = get_lnorm_pars_binom( log_sum_exp( log_Br[t], log_In[t] ),
                                           0.5*sW[t]*s0[t] );
    pars_log_Pb2[t] = get_lnorm_pars_binom( log_Pb1[t], sN[t] );
    pars_log_Pb3[t] = get_lnorm_pars_binom( log_Pb2[t], (1 - f3[t]) * sN[t] );
    pars_log_Pb4[t] = get_lnorm_pars_binom( log_Pb3[t], (1 - f4[t]) * sN[t] );
    
    // mean and sd paramters for Br and Nb summands
    pars_log_Pb2Br[t] = get_lnorm_pars_binom( log_Pb2[t], f3[t] * sN[t] );
    pars_log_Pb3Br[t] = get_lnorm_pars_binom( log_Pb3[t], f4[t] * sN[t] );
    pars_log_Pb4Br[t] = get_lnorm_pars_binom( log_Pb4[t], sN[t] );
    pars_log_BrBr[t] = get_lnorm_pars_binom( log_Br[t], bb[t] * sB[t] );
    pars_log_InBr[t] = get_lnorm_pars_binom( log_In[t], bb[t] * sB[t] );
    pars_log_NbBr[t] = get_lnorm_pars_binom( log_Nb[t], nb[t] * sN[t] );
    pars_log_BrNb[t] = get_lnorm_pars_binom( log_Br[t], (1 - bb[t]) * sB[t] );
    pars_log_InNb[t] = get_lnorm_pars_binom( log_In[t], (1 - bb[t]) * sB[t] );
    pars_log_NbNb[t] = get_lnorm_pars_binom( log_Nb[t], (1 - nb[t]) * sN[t] );
    
    // compute log pop sizes at t+1
    // compute log sizes for Pb1, ..., Pb4 and In, at time t + 1
    log_Pb1[t+1] = pars_log_Pb1[t][1] + pars_log_Pb1[t][2] * z_log_Pb1[t];
    log_Pb2[t+1] = pars_log_Pb2[t][1] + pars_log_Pb2[t][2] * z_log_Pb2[t];
    log_Pb3[t+1] = pars_log_Pb3[t][1] + pars_log_Pb3[t][2] * z_log_Pb3[t];
    log_Pb4[t+1] = pars_log_Pb4[t][1] + pars_log_Pb4[t][2] * z_log_Pb4[t];
    
    // compute log size summands for Br, Nb at time t+1
    // note: the summands are indexed at time t (not t+1)
    log_Pb2Br[t] = pars_log_Pb2Br[t][1] + pars_log_Pb2Br[t][2] * z_log_Pb2Br[t];    
    log_Pb3Br[t] = pars_log_Pb3Br[t][1] + pars_log_Pb3Br[t][2] * z_log_Pb3Br[t];    
    log_Pb4Br[t] = pars_log_Pb4Br[t][1] + pars_log_Pb4Br[t][2] * z_log_Pb4Br[t];    
    log_BrBr[t] = pars_log_BrBr[t][1] + pars_log_BrBr[t][2] * z_log_BrBr[t];    
    log_InBr[t] = pars_log_InBr[t][1] + pars_log_InBr[t][2] * z_log_InBr[t];    
    log_NbBr[t] = pars_log_NbBr[t][1] + pars_log_NbBr[t][2] * z_log_NbBr[t];     
    log_BrNb[t] = pars_log_BrNb[t][1] + pars_log_BrNb[t][2] * z_log_BrNb[t];     
    log_InNb[t] = pars_log_InNb[t][1] + pars_log_InNb[t][2] * z_log_InNb[t];
    log_NbNb[t] = pars_log_NbNb[t][1] + pars_log_NbNb[t][2] * z_log_NbNb[t];    
    
    // compute log sizes for Br, Nb at time t + 1
    // note: the summands are indexed at time t (not t+1)
    log_Br[t+1] = log_sum_exp( [ log_Pb2Br[t], log_Pb3Br[t], log_Pb4Br[t],
                                 log_BrBr[t], log_InBr[t], log_NbBr[t] ] );
    log_Nb[t+1] = log_sum_exp( [ log_BrNb[t], log_InNb[t], log_NbNb[t] ] );
  }
}

model {
  // informative priors on vital rates 
  // based on logit-multinomial approximation to HMM posterior
  logit_vr ~ multi_normal_cholesky(mean_logit_vr, L);

  // informative prior on survival to weaning probability
  // see '10_setting_prior_pup_mortality.R'
  sW ~ beta(96*5,4*5);

  // priors on mean_log_lambda and sd_log_lambda 
  // see '11_setting_immigration_priors.R'
  mean_log_lambda ~ normal(log(50), 0.3);  
  sd_log_lambda ~ normal(0, 0.4);         
  
  // informative priors on initial class sizes 
  // based on scaled stable stage distribution and chosen omega value
  // see '15_set_prior_stage_sizes_1986.R'
  exp(log_Pb1_1) ~ normal(Pb1_1_pars[1],Pb1_1_pars[2])T[0,];    
  target += log_Pb1_1;
  exp(log_Pb2_1) ~ normal(Pb2_1_pars[1],Pb2_1_pars[2])T[0,];    
  target += log_Pb2_1;
  exp(log_Pb3_1) ~ normal(Pb3_1_pars[1],Pb3_1_pars[2])T[0,];    
  target += log_Pb3_1;
  exp(log_Pb4_1) ~ normal(Pb4_1_pars[1],Pb4_1_pars[2])T[0,];    
  target += log_Pb4_1;
  exp(log_Br_1) ~ normal(Br_1_pars[1],Br_1_pars[2])T[0,];
  target += log_Br_1;
  exp(log_Nb_1) ~ normal(Nb_1_pars[1],Nb_1_pars[2])T[0,];
  target += log_Nb_1;
  
  // vague prior on count error variance
  sig_c ~ normal(0, 20);
  
  //// standard normal priors on all 'raw' variables
  // for non-centered parametrisation
  // expected number of immigrants
  z_log_lambda ~ std_normal();  
  
  // log sizes for Pb1, ..., Pb4 and In
  z_log_Pb1 ~ std_normal();
  z_log_Pb2 ~ std_normal();
  z_log_Pb3 ~ std_normal();
  z_log_Pb4 ~ std_normal();
  z_log_In ~ std_normal();
  
  // log size summands for Br, Nb
  z_log_Pb2Br ~ std_normal();
  z_log_Pb3Br ~ std_normal();
  z_log_Pb4Br ~ std_normal();
  z_log_BrBr ~ std_normal();
  z_log_InBr ~ std_normal();
  z_log_NbBr ~ std_normal();
  z_log_BrNb ~ std_normal();
  z_log_InNb ~ std_normal();
  z_log_NbNb ~ std_normal();
  
  // count likelihood
  y_count ~ normal(exp(log_Br) + exp(log_In), sig_c); 
}

generated quantities {
  // pop stage sizes
  vector[T] Pb1 = exp(log_Pb1);
  vector[T] Pb2 = exp(log_Pb2);
  vector[T] Pb3 = exp(log_Pb3);
  vector[T] Pb4 = exp(log_Pb4);
  vector[T] Br = exp(log_Br);
  vector[T] In = exp(log_In);
  vector[T] Nb = exp(log_Nb);
  // total sizes
  vector[T] Pb_tot = Pb1 + Pb2 + Pb3 + Pb4;         // pre-br classes
  vector[T] Br_tot = Br + In;                      // br classes
  vector[T] Tot = Pb_tot + Nb + Br_tot;            // total population size
  // proportion of breeding, prebreeding and non-breeding individuals
  vector[T] pPb = Pb_tot ./ Tot;
  vector[T] pNb = Nb ./ Tot;               
  vector[T] pBr_tot = Br_tot ./ Tot;
}
