// Hidden states are labeled with 'Pb' for prebreeder, 'Br' for breeder and 
// 'Nb' for non-breeder

// 1 - Pb0, underyearling 1 tag
// 2 - Pb0, underyearling 2 tag
// 3 - Pb1, 1yr old, 1 tag
// 4 - Pb1, 1yr old, 2 tag
// 5 - Pb2, 2yr old, 1 tag
// 6 - Pb2, 2yr old, 2 tag
// 7 - Pb3, 3yr old pre-breeder, 1 tag
// 8 - Pb3, 3yr old pre-breeder, 2 tag
// 9 - Pb4, pre-breeder age 4, 1 tag
// 10 - Pb4, pre-breeder age 4, 2 tag
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
}

data {
  int<lower=2> T;                   // number of years
  vector[T-1] x;                    // time covariate
  // capture-recapture data
  int<lower=1> N;                   // number of individuals
  int<lower=1> K;                   // number of hidden states
  int<lower=1> L;                   // number of observable states
  array[N,T] int<lower=0,upper=L> y;      // unique capture histories
  array[N] int<lower=1,upper=T-1> fc;     // first capture 
  array[N] int<lower=1> mult;             // capture history multiplicities
  array[N] int<lower=0,upper=2> tagloc;   // tag location factor
  
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
}

transformed data {
  // normalise x values
  real xmean = mean(x);
  real xsd = sd(x);
  array[T-1] real xn = to_array_1d((x - xmean)/xsd);    
  real sigma_intercept = 0.1;                   // ie value when sigma = 0
  vector[T-1] jitter = rep_vector(1e-9, T-1);
  
  // gaussian process covariance matrices in cholesky form
  matrix[T-1, T-1] L_mean_s0 = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_s0, sigma_mean_s0);
  matrix[T-1, T-1] L_sd_s0 = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_s0, sigma_sd_s0);
  matrix[T-1, T-1] L_mean_sN = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_sN, sigma_mean_sN);
  matrix[T-1, T-1] L_sd_sN = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_sN, sigma_sd_sN);
  matrix[T-1, T-1] L_mean_sB = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_sB, sigma_mean_sB);
  matrix[T-1, T-1] L_sd_sB = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_sB, sigma_sd_sB);
  matrix[T-1, T-1] L_mean_f3 = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_f3, sigma_mean_f3);
  matrix[T-1, T-1] L_sd_f3 = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_f3, sigma_sd_f3);
  matrix[T-1, T-1] L_mean_f4 = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_f4, sigma_mean_f4);
  matrix[T-1, T-1] L_sd_f4 = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_f4, sigma_sd_f4);
  matrix[T-1, T-1] L_mean_nb = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_nb, sigma_mean_nb);
  matrix[T-1, T-1] L_sd_nb = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_nb, sigma_sd_nb);
  matrix[T-1, T-1] L_mean_bb = gp_cov_chol(xn, sigma_intercept, jitter, ls_mean_bb, sigma_mean_bb);
  matrix[T-1, T-1] L_sd_bb = gp_cov_chol(xn, sigma_intercept, jitter, ls_sd_bb, sigma_sd_bb);
}

parameters {
  // s0
  vector[T-1] z_mean_s0;              // non-centered variables for mean process
  vector[T-1] z_sd_s0;                // non-centered variables for s0 log(std dev) process
  vector[T-1] z_logit_s0;             // non-centered variables for logit(s0)
  // sN
  vector[T-1] z_mean_sN;              // non-centered variables for mean process
  vector[T-1] z_sd_sN;                // non-centered variables for sN log(std dev) process
  vector[T-1] z_logit_sN;             // non-centered variables for logit(sN) 
  // sB
  vector[T-1] z_mean_sB;              // non-centered variables for mean process
  vector[T-1] z_sd_sB;                // non-centered variables for sB log(std dev) process
  vector[T-1] z_logit_sB;             // non-centered variables for logit(sB) 
  // f3
  vector[T-1] z_mean_f3;              // non-centered variables for mean process
  vector[T-1] z_sd_f3;                // non-centered variables for f3 log(std dev) process
  vector[T-1] z_logit_f3;             // non-centered variables for logit(f3) 
  // f4
  vector[T-1] z_mean_f4;              // non-centered variables for mean process
  vector[T-1] z_sd_f4;                // non-centered variables for f4 log(std dev) process
  vector[T-1] z_logit_f4;             // non-centered variables for logit(f4) 
  // nb
  vector[T-1] z_mean_nb;              // non-centered variables for mean process
  vector[T-1] z_sd_nb;                // non-centered variables for nb log(std dev) process
  vector[T-1] z_logit_nb;             // non-centered variables for logit(nb) 
  // bb
  vector[T-1] z_mean_bb;              // non-centered variables for mean process
  vector[T-1] z_sd_bb;                // non-centered variables for bb log(std dev) process
  vector[T-1] z_logit_bb;             // non-centered variables for logit(bb) 
  
  // detection parameters
  vector<lower=0,upper=1>[T] qN;       // nonbreeders in moult
  vector<lower=0,upper=1>[T] qB;       // breeders in moult
  vector<lower=0,upper=1>[T] pBu;      // breeders in uneven weeks of breeding season
  vector<lower=0,upper=1>[T] pBe;      // breeders in even weeks of breeding season
  
  // location specific tag loss probabilities
  real<lower=0,upper=1> ti;        // inner interdigital webbing
  real<lower=0,upper=1> to;        // outer interdigital webbing
  real<lower=0,upper=1> ti_new;    // inner interdigital webbing, 2015 onwards
}

transformed parameters {
  // Calculate  demogrpahic rates from gp parameters and non-centered versions
  vector[T-1] s0 = gp_cov_chol_to_prob(L_mean_s0, L_sd_s0, z_mean_s0, z_sd_s0, z_logit_s0);
  vector[T-1] sN = gp_cov_chol_to_prob(L_mean_sN, L_sd_sN, z_mean_sN, z_sd_sN, z_logit_sN);
  vector[T-1] sB = gp_cov_chol_to_prob(L_mean_sB, L_sd_sB, z_mean_sB, z_sd_sB, z_logit_sB);
  vector[T-1] f3 = gp_cov_chol_to_prob(L_mean_f3, L_sd_f3, z_mean_f3, z_sd_f3, z_logit_f3);
  vector[T-1] f4 = gp_cov_chol_to_prob(L_mean_f4, L_sd_f4, z_mean_f4, z_sd_f4, z_logit_f4);
  vector[T-1] bb = gp_cov_chol_to_prob(L_mean_bb, L_sd_bb, z_mean_bb, z_sd_bb, z_logit_bb);
  vector[T-1] nb = gp_cov_chol_to_prob(L_mean_nb, L_sd_nb, z_mean_nb, z_sd_nb, z_logit_nb);
}

model {
  // declare transition and emission matrices
  array[T-1] matrix[K,K] theta_in;
  array[T-1] matrix[K,K] theta_out;
  array[T-1] matrix[K,K] theta_newin;
  array[T] matrix[K,L] phi;
  
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

  // define transition matrices
  // theta_in with parameter 'ti'
  for (t in 1:(T-1)){
    theta_in[t] = get_theta(K, s0[t], sN[t], sB[t], f3[t], f4[t], nb[t], bb[t], ti);
  }
  // theta_out with parameter 'to'
  for (t in 1:(T-1)){
    theta_out[t] = get_theta(K, s0[t], sN[t], sB[t], f3[t], f4[t], nb[t], bb[t], to);
  }
  // theta_ with parameter 'ti_new'
  for (t in 1:(T-1)){
    theta_newin[t] = get_theta(K, s0[t], sN[t], sB[t], f3[t], f4[t], nb[t], bb[t], ti_new);
  }
  // define emission matrices
  for (t in 1:T){
    phi[t] = get_phi(K, L, qN[t], qB[t], pBu[t], pBe[t]);
  }
  // capture-recapture likelihood, reduce_sum implementation
  int grainsize = 1;
  target += reduce_sum(partial_sum,y,grainsize,
  T,fc,mult,tagloc,K,L,theta_in,theta_out,theta_newin,phi);
}
