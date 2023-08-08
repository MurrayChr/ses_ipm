# Simulate count and capture-mark-recapture data under the IPM model
# with specified levels of temporal variation in immigration.

#' The count data from 1986 - 2021 are simulated at the population level, and 
#' mark-recapture data are simulated at the individual level for 1983-2021, for
#' cohorts tagged 1983 - 2020. We simulate the same number of capture histories
#' in each cohort as there are in the real data, and convert it to a reduced 
#' representation in terms of unique capture histories and their multiplicities.
#' 
#' Saves to 'data' subdirectory.

library(tidyverse)

sim_ipm_data <- function(
    s0, sN, sB, f3, f4, bb, nb, sW,                # vital rates
    qN, qB, pBu, pBe,                              # detection probabilities
    ti, to, ti_new,                                # tag loss probabilities
    Pb1_1, Pb2_1, Pb3_1, Pb4_1, Br_1, Nb_1,        # initial stage sizes
    mean_log_lambda, sd_log_lambda,                # lambda is expected number immigrants
    sim_id ) {              
  
  # checks size of inputs
  Tmr <- 39                                      # number of years mark_recapture, 1983:2021
  Tco <- 36                                      # number of years count data, 1986:2021   
  stopifnot( length(s0) == Tmr-1, length(sN) == Tmr-1, length(sB) == Tmr-1,
             length(f3) == Tmr-1, length(f4) == Tmr-1, length(bb) == Tmr-1, 
             length(nb) == Tmr-1, length(sW) == Tco -1 )
  stopifnot( length(qN) == Tmr, length(qB) == Tmr, length(pBu) == Tmr, 
             length(pBe) == Tmr )
  stopifnot( length(ti) == 1, length(to) == 1, length(ti_new) == 1, 
             length(Pb1_1) == 1, length(Pb2_1) == 1, length(Pb3_1) == 1, 
             length(Pb4_1) == 1, length(Br_1) == 1, length(Nb_1) == 1,
             length(mean_log_lambda) == 1, length(sd_log_lambda) == 1)
  
  # check that variation in immigration is consitent with sim_id name
  if ( str_detect( sim_id, "low") ) {
    stopifnot( str_detect( sim_id, "low") & (sd_log_lambda == 0)  )
  }
  if ( str_detect( sim_id, "med") ) {
    stopifnot( str_detect( sim_id, "med") & (sd_log_lambda == 0.2)  )
  }
  if ( str_detect( sim_id, "high") ) {
    stopifnot( str_detect( sim_id, "high") & (sd_log_lambda == 0.4)  )
  }
  
  # tensor product utility for defining transition and emission matrices 
  # to generate cmr data
  tensor_product <- function(A,B){
    m <- dim(A)[1]
    n <- dim(A)[2]
    M <- dim(B)[1]
    N <- dim(B)[2]
    C <- matrix(0,nrow=m*M,ncol=n*N)
    for (i in 1:m){
      for (j in 1:n){
        ifirst <- (i-1)*M + 1
        ilast <- i*M
        jfirst <- (j-1)*N + 1
        jlast <- j*N
        C[ifirst:ilast,jfirst:jlast] <- A[i,j]*B
      }
    }
    C
  }
  
  # Simulate population 1986 - 2021
  
  # draw number of immigrants in all years
  if ( sd_log_lambda != 0 ) {
    log_lambda <- rnorm(Tco, mean_log_lambda, sd_log_lambda)
    In <- rpois(Tco, exp( log_lambda ) )
  } else {
    In <- rep( floor( exp( mean_log_lambda ) ), Tco)
  }
  
  # vectors of population stage-class sizes
  Pb1 <- c(Pb1_1 ,rep(NA, Tco - 1))
  Pb2 <- c(Pb2_1 ,rep(NA, Tco - 1))
  Pb3 <- c(Pb3_1 ,rep(NA, Tco - 1))
  Pb4 <- c(Pb4_1 ,rep(NA, Tco - 1))
  Br <- c(Br_1 ,rep(NA, Tco - 1))
  Nb <- c(Nb_1 ,rep(NA, Tco - 1))
  
  # simulate the stage-structured population
  for ( t in 1:(Tco - 1) ) {
    # index for the vital rates s0, ..., nb but not sW
    tmr <- t + 3
    # draw survival and state transitions to t+1
    Pb1[t+1] <- rbinom( 1, Br[t] + In[t], 0.5*sW[t]*s0[tmr] )
    Pb2[t+1] <- rbinom( 1, Pb1[t], sN[tmr] )
    Pb3[t+1] <- rbinom( 1, Pb2[t], sN[tmr]*(1 - f3[tmr]) )
    Pb4[t+1] <- rbinom( 1, Pb3[t], sN[tmr]*(1 - f4[tmr]) )
    Br[t+1] <- rbinom( 1, Pb2[t], sN[tmr]*f3[tmr] ) + 
      rbinom( 1, Pb3[t], sN[tmr]*f4[tmr] ) + 
      rbinom( 1, Pb4[t], sN[tmr] ) +
      rbinom( 1 , Br[t] + In[t], sB[tmr]*bb[tmr] ) +
      rbinom( 1, Nb[t], sN[tmr]*nb[tmr] )
    Nb[t+1] <- rbinom( 1 , Br[t] + In[t], sB[tmr]*(1 - bb[tmr]) ) +
      rbinom( 1, Nb[t], sN[tmr]*(1 - nb[tmr]) )
  }
  
  # Simulate count data
  y_count <- rnorm(Tco, Br + In, 10)
  
  # Simulate cmr data 1983 - 2021 
  data_structure <- readRDS("data/cmr_dataset_structure.RDS")
  N <- data_structure$num_indiv
  fc <- data_structure$indiv_inits$fc
  init_no_tags <- data_structure$indiv_inits$init_no_tags
  tagloc <- data_structure$indiv_inits$tagloc
  
  # create transition matrices
  get_trans_mat <- function(s0, sN, sB, f3, f4, bb, nb, tau) {
    # transition between biological states
    trans <- matrix(c(0, s0,   0,         0,          0,      0,         0,   # Pb0
                      0,   0, sN,         0,          0,      0,         0,   # Pb1
                      0,   0,  0, sN*(1-f3),          0,  sN*f3,         0,   # Pb2
                      0,   0,  0,         0,  sN*(1-f4),  sN*f4,         0,   # Pb3
                      0,   0,  0,         0,          0,     sN,         0,   # Pb4
                      0,   0,  0,         0,          0,  sB*bb, sB*(1-bb),   # Br
                      0,   0,  0,         0,          0,  sN*nb, sN*(1-nb)),  # Nb
                    nrow=7, ncol=7, byrow=TRUE)
    # transitions between tag states
    tag <- matrix(c( 1- tau,   0,     # 1 tag
                     tau, 1- tau),    # 2 tags
                  nrow=2, ncol=2, byrow=TRUE)
    # transitions between biological-tag states
    trans <- tensor_product( trans, tag )
    # append last row for dead/emigrated/untagged state
    trans <- rbind( trans, rep(0,14) )
    # append last col, enforcing row-sum-to-one constraint
    trans <- cbind( trans, 1 - rowSums(trans) )
    # return transition matrix
    trans
  }
  
  # create three lists of transition matrices according to tag location
  trans_ti <- list()
  trans_to <- list()
  trans_ti_new <- list()
  for (t in 1:(Tmr - 1)) {
    trans_ti[[t]] <- get_trans_mat(s0[t], sN[t], sB[t], f3[t], f4[t], bb[t], nb[t], ti)
    trans_to[[t]] <- get_trans_mat(s0[t], sN[t], sB[t], f3[t], f4[t], bb[t], nb[t], to)
    trans_ti_new[[t]] <- get_trans_mat(s0[t], sN[t], sB[t], f3[t], f4[t], bb[t], nb[t], ti_new)
  }
  
  # matrix of true biological and tag states
  Z <- matrix(0, N, Tmr)
  for ( i in 1:N ) {
    t_fc <- fc[i]
    ntag <- init_no_tags[i]
    Z[i,t_fc] <- ifelse(ntag == 2, 2, 1) # initialise in Pb0 with given no of tags
    # set transition matrices according to tag location
    tag <- tagloc[i]   # 1 for 'ti', 0 for 'to', 2 for 'ti_new'
    trans <- switch(tag + 1, trans_to, trans_ti, trans_ti_new) # tag + 1 to shift into range 1:3
    for (t in t_fc:(Tmr - 1)) {
      Z[i, t + 1] <- sample(1:15, 1, prob = trans[[t]][ Z[i,t], ] )
    }
  } 
  
  # capture history matrix
  # define emission matrix
  get_emit_mat <- function(qN, qB, pBu, pBe) {
    emit <- matrix(c(          0,              0,              0,                  1,              0,                  0,                  0,   # Pb0
                               0,              0,              0,                 qN,              0,                  0,                  0,   # Pb1
                               0,              0,              0,                 qN,              0,                  0,                  0,   # Pb2
                               0,              0,              0,                 qN,              0,                  0,                  0,   # Pb3
                               0,              0,              0,                 qN,              0,                  0,                  0,   # Pb4
                               qB*pBu*pBe, qB*pBu*(1-pBe), qB*(1-pBu)*pBe, qB*(1-pBu)*(1-pBe), (1-qB)*pBu*pBe, (1-qB)*pBu*(1-pBe), (1-qB)*(1-pBu)*pBe,   # Br
                               0,              0,              0,                 qN,              0,                  0,                  0),  # Nb
                   nrow=7, ncol=7, byrow=TRUE)
    # observe number of tags without error
    tag <- matrix(c( 1, 0,     # 1 tag
                     0, 1),    # 2 tags
                  nrow=2, ncol=2, byrow=TRUE)
    # observation states with number of tags
    emit <- tensor_product( emit, tag )
    # append last row for never observed state
    emit <- rbind( emit, rep(0,14) )
    # append last col, enforcing row-sum-to-one constraint
    emit <- cbind( emit, 1 - rowSums(emit) )
    # return emission matrix
    emit
  }
  
  emit <- list()
  for ( t in 1:Tmr ) {
    emit[[t]] <- get_emit_mat( qN[t], qB[t], pBu[t], pBe[t] )
  }
  
  # capture history matrix
  Y <- matrix(0, nrow = N, ncol = Tmr)
  for ( i in 1:N ) {
    t_fc <- fc[i]
    for ( t in t_fc:Tmr ) {
      Y[i,t] <- sample(1:15, 1, prob = emit[[t]][ Z[i,t], ] )
    }
  }
  
  # create reduced data representation
  colnames(Y) <- str_c("yr", 1983:2021) 
  data <- as_tibble(Y) %>%
    add_column( fc = fc, tagloc = as.integer(tagloc) )
  
  reduced_data <- data %>%
    group_by_all() %>%           # groups individuals with same capture history (hence same fc, tagloc)
    summarise(count=n()) %>%     # counts the number of individuals with each unique capture history
    ungroup() %>%                # forgets the grouping structure 
    arrange(fc)
  
  # shuffle rows for efficient parallelisation of likelihood calculation
  Nred <- nrow(reduced_data)
  shuffle <- sample(1:Nred, Nred, replace = FALSE)
  shuffled_reduced_data <- reduced_data[shuffle, ]
  
  # save data and simulated values
  sim_args <- list( s0=s0, sN=sN, sB=sB, f3=f3, f4=f4, bb=bb, nb=nb, sW=sW, 
                    qN=qN, qB=qB, pBu=pBu, pBe=pBe, ti=ti, to=to, ti_new=ti_new,  
                    Pb1_1=Pb1_1, Pb2_1=Pb2_1, Pb3_1=Pb3_1, Pb4_1=Pb4_1, Br_1=Br_1, 
                    Nb_1=Nb_1, 
                    mean_log_lambda=mean_log_lambda, sd_log_lambda=sd_log_lambda, 
                    sim_id=sim_id )
  sim_values <- list( Pb1=Pb1, Pb2=Pb2, Pb3=Pb3, Pb4=Pb4, Br=Br, Nb=Nb, In=In,
                      Z = Z, Y = Y)
  sim_data <- list( cmr_data = shuffled_reduced_data, 
                    count_data = tibble( year = 1986:2021, count = floor(y_count) ) )
  output <- c( sim_args, sim_values, sim_data )
  saveRDS( output, str_c("data/sim_",sim_id,".RDS") )
}



