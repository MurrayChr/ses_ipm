// IPM_GP_RE with 'Gaussian Process' priors on vital rates and 
// 'Random Effects' priors on immigration

// Hidden states are labeled with 'Pb' for prebreeder, 'Br' for breeder and 
// 'Nb' for non-breeder

// 1 - P0, underyearling 1 tag
// 2 - P0, underyearling 2 tag
// 3 - P1, 1yr old, 1 tag
// 4 - P1, 1yr old, 2 tag
// 5 - P2, 2yr old, 1 tag
// 6 - P2, 2yr old, 2 tag
// 7 - P3, 3yr old pre-breeder, 1 tag
// 8 - P3, 3yr old pre-breeder, 2 tag
// 9 - P4, pre-breeder age 4, 1 tag
// 10 - P4, pre-breeder age 4, 2 tag
// 11 - Br, breeder, 1 tag
// 12 - Br, breeder, 2 tag
// 13 - Nb, non-breeder, 1 tag
// 14 - Nb, non-breeder, 2 tag
// 15 - dead

functions {
// get cholesky decomposition of gp covariance matrix
  matrix gp_cov_chol( data array[] real xn, data real sigma_intercept, 
                      data vector jitter, real ls, real  sigma) {
    int n = num_elements(xn);
    matrix[n, n] K = gp_exp_quad_cov(xn, sigma, ls) +  sigma_intercept;
    matrix[n, n] L = cholesky_decompose(add_diag(K, jitter));
    return L;
  }
  // get probs from logit-scale gp processes in cholesky form
  vector gp_cov_chol_to_prob( matrix L_mean, matrix L_sd, 
                              vector z_mean, vector z_sd, vector z_logit_p) {
    int n = num_elements(z_logit_p);
    vector[n] logit_p = L_mean * z_mean + exp(L_sd * z_sd) .* z_logit_p;
    vector[n] p = inv_logit(logit_p);
    return p;
  }
  // partial sum function to pass to reduce_sum
  real partial_sum (array[,] int y_slice,
    int start, int end,
    int T,
    array[] int fc,
    array[] int mult,
    array[] int tagloc,
    int K, int L,
    array[] matrix theta_in,
    array[] matrix theta_out,
    array[] matrix theta_newin,
    array[] matrix phi) {
    real ll_term = 0;
    for (n in start:end){
      int ny = n - start + 1;
      array[T] row_vector[K] gamma;
      // initialise forward algorithm at first capture
      gamma[fc[n]] = rep_row_vector(0,K);
      gamma[fc[n]][y_slice[ny,fc[n]]== 8 ? 2 : 1] = 1;
      // recursion
      if (tagloc[n]==1) {
        for (t in fc[n]:(T-1)){
          row_vector[K] gamma_theta;
          gamma_theta = gamma[t]*theta_in[t];  // theta_in used here
          gamma[t+1] =phi[t+1][:,y_slice[ny,t+1]]'.*gamma_theta; 
        }
      } else if (tagloc[n]==0) {
        for (t in fc[n]:(T-1)){
          row_vector[K] gamma_theta;
          gamma_theta = gamma[t]*theta_out[t]; // theta_out used here
          gamma[t+1] =phi[t+1][:,y_slice[ny,t+1]]'.*gamma_theta; 
        }
      } else if (tagloc[n]==2) {
        for (t in fc[n]:(T-1)){
          row_vector[K] gamma_theta;
          gamma_theta = gamma[t]*theta_newin[t]; // theta_newin used here
          gamma[t+1] =phi[t+1][:,y_slice[ny,t+1]]'.*gamma_theta; 
        }
      }
      ll_term = ll_term + mult[n]*log(sum(gamma[T]));
    }
    return ll_term;
  }
  
  // construct transition matrix
  matrix get_theta (int K, real s0, real sN, real sB, real f3, real f4, real nb, real bb, real tl){
    matrix[K,K] theta;
    // P0       P0   P1                    P2                    P3                                  P4                                  Br                          Nb                                              D
    theta[1] =  [0,0, s0*(1-tl),         0,         0,         0,                0,                0,                0,                0,            0,            0,                0,                0, 1-s0*(1-tl)];
    theta[2] =  [0,0,     s0*tl, s0*(1-tl),         0,         0,                0,                0,                0,                0,            0,            0,                0,                0,        1-s0];
    // P1
    theta[3] =  [0,0,         0,         0, sN*(1-tl),         0,                0,                0,                0,                0,            0,            0,                0,                0, 1-sN*(1-tl)];
    theta[4] =  [0,0,         0,         0,     sN*tl, sN*(1-tl),                0,                0,                0,                0,            0,            0,                0,                0,        1-sN];
    // P2
    theta[5] =  [0,0,         0,         0,         0,         0, sN*(1-f3)*(1-tl),                0,                0,                0, sN*f3*(1-tl),            0,                0,                0, 1-sN*(1-tl)];
    theta[6] =  [0,0,         0,         0,         0,         0,     sN*(1-f3)*tl, sN*(1-f3)*(1-tl),                0,                0,     sN*f3*tl, sN*f3*(1-tl),                0,                0,        1-sN];
    // P3
    theta[7] =  [0,0,         0,         0,         0,         0,                0,                0, sN*(1-f4)*(1-tl),                0, sN*f4*(1-tl),            0,                0,                0, 1-sN*(1-tl)];
    theta[8] =  [0,0,         0,         0,         0,         0,                0,                0,     sN*(1-f4)*tl, sN*(1-f4)*(1-tl),     sN*f4*tl, sN*f4*(1-tl),                0,                0,        1-sN];
     // P4 
    theta[9] =  [0,0,         0,         0,         0,         0,                0,                0,                0,                0,    sN*(1-tl),            0,                0,                0, 1-sN*(1-tl)];
    theta[10] = [0,0,         0,         0,         0,         0,                0,                0,                0,                0,        sN*tl,    sN*(1-tl),                0,                0,        1-sN];
    // Br
    theta[11] = [0,0,         0,         0,         0,         0,                0,                0,                0,                0, sB*bb*(1-tl),            0, sB*(1-bb)*(1-tl),                0, 1-sB*(1-tl)];
    theta[12] = [0,0,         0,         0,         0,         0,                0,                0,                0,                0,     sB*bb*tl, sB*bb*(1-tl),     sB*(1-bb)*tl, sB*(1-bb)*(1-tl),        1-sB];
    // Nb
    theta[13] = [0,0,         0,         0,         0,         0,                0,                0,                0,                0, sN*nb*(1-tl),            0, sN*(1-nb)*(1-tl),                0, 1-sN*(1-tl)];
    theta[14] = [0,0,         0,         0,         0,         0,                0,                0,                0,                0,     sN*nb*tl, sN*nb*(1-tl),     sN*(1-nb)*tl, sN*(1-nb)*(1-tl),        1-sN];
    // Dead
    theta[15] = [0,0,         0,         0,         0,         0,                0,                0,                0,                0,            0,            0,                0,                0,           1];
    return theta;
  }

  // contruct emission matrix
  matrix get_phi (int K, int L, real qN, real qB, real pBu, real pBe){
    matrix[K,L] phi;
    phi[1] =  [                  0,                  0,                       0,                       0,                       0,                       0,                           1,                           0,                       0,                       0,                          0,                          0,                          0,                           0,                            0];
    phi[2] =  [                  0,                  0,                       0,                       0,                       0,                       0,                           0,                           1,                       0,                       0,                          0,                          0,                          0,                           0,                            0];
    phi[3] =  [                  0,                  0,                       0,                       0,                       0,                       0,                          qN,                           0,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[4] =  [                  0,                  0,                       0,                       0,                       0,                       0,                           0,                          qN,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[5] =  [                  0,                  0,                       0,                       0,                       0,                       0,                          qN,                           0,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[6] =  [                  0,                  0,                       0,                       0,                       0,                       0,                           0,                          qN,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[7] =  [                  0,                  0,                       0,                       0,                       0,                       0,                          qN,                           0,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[8] =  [                  0,                  0,                       0,                       0,                       0,                       0,                           0,                          qN,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[9] =  [                  0,                  0,                       0,                       0,                       0,                       0,                          qN,                           0,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[10] = [                  0,                  0,                       0,                       0,                       0,                       0,                           0,                          qN,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[11] = [         qB*pBu*pBe,                  0,          qB*pBu*(1-pBe),                       0,          qB*(1-pBu)*pBe,                       0,          qB*(1-pBu)*(1-pBe),                           0,          (1-qB)*pBu*pBe,                       0,         (1-qB)*pBu*(1-pBe),                          0,         (1-qB)*(1-pBu)*pBe,                           0,       (1-qB)*(1-pBu)*(1-pBe)];
    phi[12] = [                  0,         qB*pBu*pBe,                       0,          qB*pBu*(1-pBe),                       0,          qB*(1-pBu)*pBe,                           0,          qB*(1-pBu)*(1-pBe),                       0,          (1-qB)*pBu*pBe,                          0,         (1-qB)*pBu*(1-pBe),                          0,          (1-qB)*(1-pBu)*pBe,       (1-qB)*(1-pBu)*(1-pBe)];
    phi[13] = [                  0,                  0,                       0,                       0,                       0,                       0,                          qN,                           0,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[14] = [                  0,                  0,                       0,                       0,                       0,                       0,                           0,                          qN,                       0,                       0,                          0,                          0,                          0,                           0,                         1-qN];
    phi[15] = [                  0,                  0,                       0,                       0,                       0,                       0,                           0,                           0,                       0,                       0,                          0,                          0,                          0,                           0,                            1];
    return phi;
  }
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
  // count data
  int<lower=1> Tco;                        // number of years with count data
  vector<lower=0>[Tco] y_count;            // counts
  // capture-mark-recapture data
  int<lower=1> Tmr;                         // number of years of cmr data      
  vector[Tmr-1] x;                         // time covariate
  int<lower=1> N;                          // no. unique capture histores in cmr data
  int<lower=1> K;                          // no. hidden states
  int<lower=1> L;                          // no. observable states
  array[N,Tmr] int<lower=0,upper=L> y_cmr; // unique capture histories
  array[N] int<lower=1,upper=Tmr-1> fc;    // first capture 
  array[N] int<lower=1> mult;              // capture history multiplicities
  array[N] int<lower=0,upper=2> tagloc;    // tag location factor
  
  // values for gp lengthscale and marginal standard deviation hyperparameters
  // log(std dev) process
  real<lower=0> ls_sd_s0;           // lengthscale for s0  log(std dev) process          
  real<lower=0> sigma_sd_s0;        // scale for s0 
  real<lower=0> ls_sd_sN;           // lengthscale for sN  log(std dev) process          
  real<lower=0> sigma_sd_sN;        // scale for sN
  real<lower=0> ls_sd_sB;           // lengthscale for sB  log(std dev) process          
  real<lower=0> sigma_sd_sB;        // scale for sB 
  real<lower=0> ls_sd_f3;           // lengthscale for f3  log(std dev) process          
  real<lower=0> sigma_sd_f3;        // scale for f3 
  real<lower=0> ls_sd_f4;           // lengthscale for f4  log(std dev) process          
  real<lower=0> sigma_sd_f4;        // scale for f4 
  real<lower=0> ls_sd_nb;           // lengthscale for nb  log(std dev) process          
  real<lower=0> sigma_sd_nb;        // scale for nb 
  real<lower=0> ls_sd_bb;           // lengthscale for bb  log(std dev) process          
  real<lower=0> sigma_sd_bb;        // scale for bb 
  // mean process
  real<lower=0> ls_mean_s0;           // lengthscale for s0 mean process          
  real<lower=0> sigma_mean_s0;        // scale for s0 
  real<lower=0> ls_mean_sN;           // lengthscale for sN mean process          
  real<lower=0> sigma_mean_sN;        // scale for sN
  real<lower=0> ls_mean_sB;           // lengthscale for sB mean process          
  real<lower=0> sigma_mean_sB;        // scale for sB 
  real<lower=0> ls_mean_f3;           // lengthscale for f3 mean process          
  real<lower=0> sigma_mean_f3;        // scale for f3 
  real<lower=0> ls_mean_f4;           // lengthscale for f4 mean process          
  real<lower=0> sigma_mean_f4;        // scale for f4 
  real<lower=0> ls_mean_nb;           // lengthscale for nb mean process          
  real<lower=0> sigma_mean_nb;        // scale for nb 
  real<lower=0> ls_mean_bb;           // lengthscale for bb mean process          
  real<lower=0> sigma_mean_bb;        // scale for bb 
  
  // mean and sd parameters for initial stage size priors
  array[2] real<lower=0> Pb1_1_pars;
  array[2] real<lower=0> Pb2_1_pars;
  array[2] real<lower=0> Pb3_1_pars;
  array[2] real<lower=0> Pb4_1_pars;
  array[2] real<lower=0> Br_1_pars;
  array[2] real<lower=0> Nb_1_pars;
}

transformed data {
  // normalise x values
  real xmean = mean(x);
  real xsd = sd(x);
  array[Tmr-1] real xn = to_array_1d((x - xmean)/xsd);    
  real sigma_intercept = 0.1;                   // ie value when sigma = 0
  vector[Tmr-1] jitter = rep_vector(1e-9, Tmr-1);
  
  // gaussian process covariance matrices in cholesky form
  matrix[Tmr-1, Tmr-1] L_mean_s0 = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_s0, sigma_mean_s0);
  matrix[Tmr-1, Tmr-1] L_sd_s0 = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_s0, sigma_sd_s0);
  matrix[Tmr-1, Tmr-1] L_mean_sN = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_sN, sigma_mean_sN);
  matrix[Tmr-1, Tmr-1] L_sd_sN = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_sN, sigma_sd_sN);
  matrix[Tmr-1, Tmr-1] L_mean_sB = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_sB, sigma_mean_sB);
  matrix[Tmr-1, Tmr-1] L_sd_sB = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_sB, sigma_sd_sB);
  matrix[Tmr-1, Tmr-1] L_mean_f3 = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_f3, sigma_mean_f3);
  matrix[Tmr-1, Tmr-1] L_sd_f3 = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_f3, sigma_sd_f3);
  matrix[Tmr-1, Tmr-1] L_mean_f4 = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_f4, sigma_mean_f4);
  matrix[Tmr-1, Tmr-1] L_sd_f4 = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_f4, sigma_sd_f4);
  matrix[Tmr-1, Tmr-1] L_mean_nb = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_nb, sigma_mean_nb);
  matrix[Tmr-1, Tmr-1] L_sd_nb = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_nb, sigma_sd_nb);
  matrix[Tmr-1, Tmr-1] L_mean_bb = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_bb, sigma_mean_bb);
  matrix[Tmr-1, Tmr-1] L_sd_bb = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_bb, sigma_sd_bb);
}

parameters {
  // s0
  vector[Tmr-1] z_mean_s0;              // non-centered variables for mean process
  vector[Tmr-1] z_sd_s0;                // non-centered variables for s0 log(std dev) process
  vector[Tmr-1] z_logit_s0;             // non-centered variables for logit(s0)
  // sN
  vector[Tmr-1] z_mean_sN;              // non-centered variables for mean process
  vector[Tmr-1] z_sd_sN;                // non-centered variables for sN log(std dev) process
  vector[Tmr-1] z_logit_sN;             // non-centered variables for logit(sN) 
  // sB
  vector[Tmr-1] z_mean_sB;              // non-centered variables for mean process
  vector[Tmr-1] z_sd_sB;                // non-centered variables for sB log(std dev) process
  vector[Tmr-1] z_logit_sB;             // non-centered variables for logit(sB) 
  // f3
  vector[Tmr-1] z_mean_f3;              // non-centered variables for mean process
  vector[Tmr-1] z_sd_f3;                // non-centered variables for f3 log(std dev) process
  vector[Tmr-1] z_logit_f3;             // non-centered variables for logit(f3) 
  // f4
  vector[Tmr-1] z_mean_f4;              // non-centered variables for mean process
  vector[Tmr-1] z_sd_f4;                // non-centered variables for f4 log(std dev) process
  vector[Tmr-1] z_logit_f4;             // non-centered variables for logit(f4) 
  // nb
  vector[Tmr-1] z_mean_nb;              // non-centered variables for mean process
  vector[Tmr-1] z_sd_nb;                // non-centered variables for nb log(std dev) process
  vector[Tmr-1] z_logit_nb;             // non-centered variables for logit(nb) 
  // bb
  vector[Tmr-1] z_mean_bb;              // non-centered variables for mean process
  vector[Tmr-1] z_sd_bb;                // non-centered variables for bb log(std dev) process
  vector[Tmr-1] z_logit_bb;             // non-centered variables for logit(bb) 
  
  // survival to weaning probability
  vector<lower=0, upper=1>[Tco-1] sW;
  
  // detection parameters
  vector<lower=0,upper=1>[Tmr] qN;       // nonbreeders in moult
  vector<lower=0,upper=1>[Tmr] qB;       // breeders in moult
  vector<lower=0,upper=1>[Tmr] pBu;      // breeders in uneven weeks of breeding season
  vector<lower=0,upper=1>[Tmr] pBe;      // breeders in even weeks of breeding season
  
  // location specific tag loss probabilities
  real<lower=0,upper=1> ti;        // inner interdigital webbing
  real<lower=0,upper=1> to;        // outer interdigital webbing
  real<lower=0,upper=1> ti_new;    // inner interdigital webbing, 2015 onwards
  
  // initial log pop stage sizes
  real log_Pb1_1;
  real log_Pb2_1;
  real log_Pb3_1;
  real log_Pb4_1;
  real log_Br_1;
  real log_Nb_1;
  
  // raw log pop stage sizes for Pb1, ..., Pb4, In
  // for non-centered parametrisation
  vector[Tco-1] z_log_Pb1;
  vector[Tco-1] z_log_Pb2;       
  vector[Tco-1] z_log_Pb3;      
  vector[Tco-1] z_log_Pb4;
  vector[Tco] z_log_In;
  
  // raw log summands for Br and Nb stage sizes
  // for non-centered parmetrisation
  vector[Tco-1] z_log_Pb2Br;    
  vector[Tco-1] z_log_Pb3Br;
  vector[Tco-1] z_log_Pb4Br;
  vector[Tco-1] z_log_BrBr;
  vector[Tco-1] z_log_InBr;
  vector[Tco-1] z_log_NbBr;
  vector[Tco-1] z_log_BrNb;
  vector[Tco-1] z_log_InNb;
  vector[Tco-1] z_log_NbNb;
  
  // 'poisson-lognormal' prior parameters on immigration
  vector[Tco] z_log_lambda;      // raw "poisson parameter" lambda on log scale
  real mean_log_lambda;          // mean for hyperprior on poisson parameter 
  real<lower=0> sd_log_lambda;   // std.dev for hyperprior on poisson parameter
}

transformed parameters {
  // Calculate  demogrpahic rates from gp parameters and non-centered versions
  vector[Tmr-1] s0 = gp_cov_chol_to_prob(L_mean_s0, L_sd_s0, z_mean_s0, z_sd_s0, z_logit_s0);
  vector[Tmr-1] sN = gp_cov_chol_to_prob(L_mean_sN, L_sd_sN, z_mean_sN, z_sd_sN, z_logit_sN);
  vector[Tmr-1] sB = gp_cov_chol_to_prob(L_mean_sB, L_sd_sB, z_mean_sB, z_sd_sB, z_logit_sB);
  vector[Tmr-1] f3 = gp_cov_chol_to_prob(L_mean_f3, L_sd_f3, z_mean_f3, z_sd_f3, z_logit_f3);
  vector[Tmr-1] f4 = gp_cov_chol_to_prob(L_mean_f4, L_sd_f4, z_mean_f4, z_sd_f4, z_logit_f4);
  vector[Tmr-1] bb = gp_cov_chol_to_prob(L_mean_bb, L_sd_bb, z_mean_bb, z_sd_bb, z_logit_bb);
  vector[Tmr-1] nb = gp_cov_chol_to_prob(L_mean_nb, L_sd_nb, z_mean_nb, z_sd_nb, z_logit_nb);
  
  //// declare variables for the population model
  
  // expected number of immigrants in each year
  vector<lower=0>[Tco] lambda;
  
  // log pop stage sizes 
  vector[Tco] log_Pb1;         // one-year-old prebreeders
  vector[Tco] log_Pb2;         // two-year-old prebreeders
  vector[Tco] log_Pb3;         // three-year-old prebreeders   
  vector[Tco] log_Pb4;         // four-year-old prebreeders
  vector[Tco] log_Br;          // all breeders (excluding new immigrants)
  vector[Tco] log_In;          // newly arrived breeding immigrants
  vector[Tco] log_Nb;          // all nonbreeders 
  
  // log summands for Br and Nb stage sizes
  vector[Tco-1] log_Pb2Br;    // three-year-old first-time breeders
  vector[Tco-1] log_Pb3Br;    // four-year-old first-time breeders
  vector[Tco-1] log_Pb4Br;    // five-year-old first-time breeders
  vector[Tco-1] log_BrBr;     // surviving breeders transition to breeders
  vector[Tco-1] log_InBr;     // surviving immigrants transition to breeders
  vector[Tco-1] log_NbBr;     // surviving nonbreeders transition to breeders
  vector[Tco-1] log_BrNb;     // surviving breeders transition to nonbreeders
  vector[Tco-1] log_InNb;     // surviving immigrants transition to nonbreeders
  vector[Tco-1] log_NbNb;     // surviving nonbreeders transition to nonbreeders
    
  // lognormal mean and sd parameters for Pb1, ..., Pb4, In
  array[Tco-1] row_vector[2] pars_log_Pb1;
  array[Tco-1] row_vector[2] pars_log_Pb2;
  array[Tco-1] row_vector[2] pars_log_Pb3;
  array[Tco-1] row_vector[2] pars_log_Pb4;
  array[Tco] row_vector[2] pars_log_In;
  
  // lognormal mean and sd parameters for Br and Nb summands
  array[Tco-1] row_vector[2] pars_log_Pb2Br;
  array[Tco-1] row_vector[2] pars_log_Pb3Br;
  array[Tco-1] row_vector[2] pars_log_Pb4Br;
  array[Tco-1] row_vector[2] pars_log_InBr;
  array[Tco-1] row_vector[2] pars_log_BrBr;
  array[Tco-1] row_vector[2] pars_log_NbBr;
  array[Tco-1] row_vector[2] pars_log_InNb;
  array[Tco-1] row_vector[2] pars_log_BrNb;
  array[Tco-1] row_vector[2] pars_log_NbNb;
  
  //// define the population model
  
  // immigrants at all times
  for (t in 1:Tco){
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
  
  for (t in 1:(Tco-1)) {  // t indexes time starting in first count year: 1986
    int tmr = t + 3;     // tmr indexes time starting in first cmr year: 1983
    // compute mean and sd parameters
    //  mean and sd parameters for Pb1, ..., Pb4, In
    pars_log_Pb1[t] = get_lnorm_pars_binom( log_sum_exp( log_Br[t], log_In[t] ),
                                           0.5*sW[t]*s0[tmr] );
    pars_log_Pb2[t] = get_lnorm_pars_binom( log_Pb1[t], sN[tmr] );
    pars_log_Pb3[t] = get_lnorm_pars_binom( log_Pb2[t], (1 - f3[tmr]) * sN[tmr] );
    pars_log_Pb4[t] = get_lnorm_pars_binom( log_Pb3[t], (1 - f4[tmr]) * sN[tmr] );
    
    // mean and sd paramters for Br and Nb summands
    pars_log_Pb2Br[t] = get_lnorm_pars_binom( log_Pb2[t], f3[tmr] * sN[tmr] );
    pars_log_Pb3Br[t] = get_lnorm_pars_binom( log_Pb3[t], f4[tmr] * sN[tmr] );
    pars_log_Pb4Br[t] = get_lnorm_pars_binom( log_Pb4[t], sN[tmr] );
    pars_log_BrBr[t] = get_lnorm_pars_binom( log_Br[t], bb[tmr] * sB[tmr] );
    pars_log_InBr[t] = get_lnorm_pars_binom( log_In[t], bb[tmr] * sB[tmr] );
    pars_log_NbBr[t] = get_lnorm_pars_binom( log_Nb[t], nb[tmr] * sN[tmr] );
    pars_log_BrNb[t] = get_lnorm_pars_binom( log_Br[t], (1 - bb[tmr]) * sB[tmr] );
    pars_log_InNb[t] = get_lnorm_pars_binom( log_In[t], (1 - bb[tmr]) * sB[tmr] );
    pars_log_NbNb[t] = get_lnorm_pars_binom( log_Nb[t], (1 - nb[tmr]) * sN[tmr] );
    
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
  // declare transition and emission matrices
  array[Tmr-1] matrix[K,K] theta_in;
  array[Tmr-1] matrix[K,K] theta_out;
  array[Tmr-1] matrix[K,K] theta_newin;
  array[Tmr] matrix[K,L] phi;
  
  // priors
  // s0
  z_mean_s0 ~ std_normal();
  z_sd_s0 ~ std_normal();
  z_logit_s0 ~ std_normal();
  // sN
  z_mean_sN ~ std_normal();
  z_sd_sN ~ std_normal();
  z_logit_sN ~ std_normal();
  // sB
  z_mean_sB ~ std_normal();
  z_sd_sB ~ std_normal();
  z_logit_sB ~ std_normal();
  // f3
  z_mean_f3 ~ std_normal();
  z_sd_f3 ~ std_normal();
  z_logit_f3 ~ std_normal();
  // f4
  z_mean_f4 ~ std_normal();
  z_sd_f4 ~ std_normal();
  z_logit_f4 ~ std_normal();
  // bb
  z_mean_bb ~ std_normal();
  z_sd_bb ~ std_normal();
  z_logit_bb ~ std_normal();
  // nb
  z_mean_nb ~ std_normal();
  z_sd_nb ~ std_normal();
  z_logit_nb ~ std_normal();
  
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
  y_count ~ normal(exp(log_Br) + exp(log_In), 10); 
  
  // define transition matrices
  // theta_in with parameter 'ti'
  for (t in 1:(Tmr-1)){
    theta_in[t] = get_theta(K, s0[t], sN[t], sB[t], f3[t], f4[t], nb[t], bb[t], ti);
  }
  // theta_out with parameter 'to'
  for (t in 1:(Tmr-1)){
    theta_out[t] = get_theta(K, s0[t], sN[t], sB[t], f3[t], f4[t], nb[t], bb[t], to);
  }
  // theta_newin with parameter 'ti_new'
  for (t in 1:(Tmr-1)){
    theta_newin[t] = get_theta(K, s0[t], sN[t], sB[t], f3[t], f4[t], nb[t], bb[t], ti_new);
  }
  // define emission matrices
  for (t in 1:Tmr){
    phi[t] = get_phi(K, L, qN[t], qB[t], pBu[t], pBe[t]);
  }
  // capture-recapture likelihood, reduce_sum implementation
  int grainsize = 1;
  target += reduce_sum(partial_sum,y_cmr,grainsize,
  Tmr,fc,mult,tagloc,K,L,theta_in,theta_out,theta_newin,phi);
}

generated quantities {
  // pop stage sizes
  vector[Tco] Pb1 = exp(log_Pb1);
  vector[Tco] Pb2 = exp(log_Pb2);
  vector[Tco] Pb3 = exp(log_Pb3);
  vector[Tco] Pb4 = exp(log_Pb4);
  vector[Tco] Br = exp(log_Br);
  vector[Tco] In = exp(log_In);
  vector[Tco] Nb = exp(log_Nb);
  // total sizes
  vector[Tco] Pb_tot = Pb1 + Pb2 + Pb3 + Pb4;         // pre-br classes
  vector[Tco] Br_tot = Br + In;                      // br classes
  vector[Tco] Tot = Pb_tot + Nb + Br_tot;            // total population size
  // proportion of breeding, prebreeding and non-breeding individuals
  vector[Tco] pPb = Pb_tot ./ Tot;
  vector[Tco] pNb = Nb ./ Tot;               
  vector[Tco] pBr_tot = Br_tot ./ Tot;
}
