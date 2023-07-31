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
  // capture-recapture data
  int<lower=1> N;                   // number of individuals
  int<lower=1> K;                   // number of hidden states
  int<lower=1> L;                   // number of observable states
  array[N,T] int<lower=0,upper=L> y;      // unique capture histories
  array[N] int<lower=1,upper=T-1> fc;     // first capture 
  array[N] int<lower=1> mult;             // capture history multiplicities
  array[N] int<lower=0,upper=2> tagloc;   // tag location factor
}

parameters {
  // detection parameters
  vector<lower=0,upper=1>[T] qN;       // non-breeders in moult
  vector<lower=0,upper=1>[T] qB;       // breeders in moult
  vector<lower=0,upper=1>[T] pBu;      // breeders in uneven weeks of breeding season
  vector<lower=0,upper=1>[T] pBe;      // breeders in even weeks of breeding season
  // location specific tag loss probabilities
  real<lower=0,upper=1> ti;        // inner interdigital webbing
  real<lower=0,upper=1> to;        // outer interdigital webbing
  real<lower=0,upper=1> ti_new;    // inner interdigital webbing, 2015 onwards
  // hyperprior parameters for logit-scale survival and conditional breeding probabilities
  real mu_logit_s0;
  real<lower=0> sig_logit_s0;
  real mu_logit_sN;
  real<lower=0> sig_logit_sN;
  real mu_logit_sB;
  real<lower=0> sig_logit_sB;
  real mu_logit_f3;
  real<lower=0> sig_logit_f3;
  real mu_logit_f4;
  real<lower=0> sig_logit_f4;
  real mu_logit_bb;
  real<lower=0> sig_logit_bb;
  real mu_logit_nb;
  real<lower=0> sig_logit_nb;
  // raw logit scale probabilities get std normal priors
  vector[T-1] z_logit_s0;
  vector[T-1] z_logit_sN;
  vector[T-1] z_logit_sB;
  vector[T-1] z_logit_f3;
  vector[T-1] z_logit_f4;
  vector[T-1] z_logit_bb;
  vector[T-1] z_logit_nb;
}

transformed parameters {
  // declare survival and conditional breeding probabilities
  vector<lower=0,upper=1>[T-1] s0;
  vector<lower=0,upper=1>[T-1] sN;
  vector<lower=0,upper=1>[T-1] sB;
  vector<lower=0,upper=1>[T-1] f3;    // prob first breeding age 3
  vector<lower=0,upper=1>[T-1] f4;    // prob first breeding age 4
  vector<lower=0,upper=1>[T-1] nb;    // non-breeder to breeder
  vector<lower=0,upper=1>[T-1] bb;    // breeder to breeder
  // define survival probabilities
  s0 = inv_logit(mu_logit_s0 + sig_logit_s0*z_logit_s0);
  sN = inv_logit(mu_logit_sN + sig_logit_sN*z_logit_sN);
  sB = inv_logit(mu_logit_sB + sig_logit_sB*z_logit_sB);
  // define state transition probabilities
  f3 = inv_logit(mu_logit_f3 + sig_logit_f3*z_logit_f3);
  f4 = inv_logit(mu_logit_f4 + sig_logit_f4*z_logit_f4);
  bb = inv_logit(mu_logit_bb + sig_logit_bb*z_logit_bb);
  nb = inv_logit(mu_logit_nb + sig_logit_nb*z_logit_nb);
}

model {
  // declare transition and emission matrices
  array[T-1] matrix[K,K] theta_in;
  array[T-1] matrix[K,K] theta_out;
  array[T-1] matrix[K,K] theta_newin;
  array[T] matrix[K,L] phi;
  
  // priors on hyperparameters 
  // imply approx uniform prior on parameters
  mu_logit_s0 ~ normal(0,1.25);
  sig_logit_s0 ~ exponential(1.5);
  mu_logit_sN ~ normal(0,1.25);
  sig_logit_sN ~ exponential(1.5);
  mu_logit_sB ~ normal(0,1.25);
  sig_logit_sB ~ exponential(1.5);
  mu_logit_f3 ~ normal(0,1.25);
  sig_logit_f3 ~ exponential(1.5);
  mu_logit_f4 ~ normal(0,1.25);
  sig_logit_f4 ~ exponential(1.5);
  mu_logit_bb ~ normal(0,1.25);
  sig_logit_bb ~ exponential(1.5);
  mu_logit_nb ~ normal(0,1.25);
  sig_logit_nb ~ exponential(1.5);
  
  // standard normal priors 
  z_logit_s0 ~ std_normal();
  z_logit_sN ~ std_normal();
  z_logit_sB ~ std_normal();
  z_logit_f3 ~ std_normal();
  z_logit_f4 ~ std_normal();
  z_logit_bb ~ std_normal();
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
