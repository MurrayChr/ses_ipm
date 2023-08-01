# Function to simulate replicate capture-mark-recapture data
# using as 'true' parameter values a posterior draw from an hmm fit
library(tidyverse)

# Simulating replicate cmr datasets --------------------------------------------

sim_rep_cmr_data <- function(fit,       # cmdstanr fit for hmm model
                             draw,      # posterior draw to use as 'true' parameter values
                             vars = c("s0", "sN", "sB", "f3", "f4", "bb", "nb",
                                      "ti", "to", "ti_new", "qN", "qB", "pBu", "pBe")
                             ) {

  # Unpack dataset structure attributes ----------------------------------------
  data_str <- readRDS("data/cmr_dataset_structure.RDS")
  N <- data_str$num_indiv  
  Tmr <- data_str$num_years
  fc <- data_str$indiv_inits$fc
  init_no_tags <- data_str$indiv_inits$init_no_tags
  tagloc <- data_str$indiv_inits$tagloc

  # Format posterior sample ----------------------------------------------------
  # utility functions to extract variable names and time index
  # e.g. "s0[7]" gives variable name "s0" and integer '7' 
  get_var_name <- function(str){str_split(str,"\\[")[[1]][1]}
  get_var_time <- function(str){
    s <- str_extract(str,"\\[[^()]+\\]")  
    s <- substring(s,2,str_length(s)-1)
    parse_integer(s)
  }
  
  # check draw is possible
  if (!all(draw > 0, draw <= dim(fit$draws(format = "df"))[1])){
    return(str_c("Error: ",draw," exceeds number of posterior samples."))
  }
  
  # format posterior sample
  posterior_sample <- fit$draws(variables = vars, format = "df") %>%
    as_tibble() %>%
    filter(.draw == draw) %>% 
    select(-c(".chain",".iteration",".draw")) %>%
    pivot_longer(everything()) %>%
    mutate(par = sapply(name, get_var_name), t = sapply(name,get_var_time), .after = name)
  
  # decant posterior sample values from tibble to vectors
  s0 <- filter(posterior_sample,par=="s0")$value
  sN <- filter(posterior_sample,par=="sN")$value
  sB <- filter(posterior_sample,par=="sB")$value
  f3 <- filter(posterior_sample,par=="f3")$value
  f4 <- filter(posterior_sample,par=="f4")$value
  bb <- filter(posterior_sample,par=="bb")$value
  nb <- filter(posterior_sample,par=="nb")$value
  ti <- filter(posterior_sample,par=="ti")$value
  to <- filter(posterior_sample,par=="to")$value
  ti_new <- filter(posterior_sample,par=="ti_new")$value
  qN <- filter(posterior_sample,par=="qN")$value
  qB <- filter(posterior_sample,par=="qB")$value
  pBu <- filter(posterior_sample,par=="pBu")$value
  pBe <- filter(posterior_sample,par=="pBe")$value
  
  # Define transition and emission matrices ------------------------------------
  # tensor product utility function
  # to construct transition and emission matrices we first construct preliminary versions that
  # contain only probabilities of transitioning between biological states, then take the
  # tensor product of these with matrices that contain transition probabilities between tag states
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
  
  # Simulate true states and observed capture histories ------------------------
  # matrix of true biological and tag states
  z <- matrix(0, N, Tmr)
  for ( i in 1:N ) {
    t_fc <- fc[i]
    ntag <- init_no_tags[i]
    z[i,t_fc] <- ifelse(ntag == 2, 2, 1) # initialise in Pb0 with given no of tags
    # set transition matrices according to tag location
    tag <- tagloc[i]   # 1 for 'ti', 0 for 'to', 2 for 'ti_new'
    trans <- switch(tag + 1, trans_to, trans_ti, trans_ti_new) # tag + 1 to shift into range 1:3
    for (t in t_fc:(Tmr - 1)) {
      z[i, t + 1] <- sample(1:15, 1, prob = trans[[t]][ z[i,t], ] )
    }
  } 
  
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
  y <- matrix(0, nrow = N, ncol = Tmr)
  for ( i in 1:N ) {
    t_fc <- fc[i]
    for ( t in t_fc:Tmr ) {
      y[i,t] <- sample(1:15, 1, prob = emit[[t]][ z[i,t], ] )
    }
  }
  colnames(y) <- 1983:2021
  
  return( list( draw = draw, values = posterior_sample, y = y,
                trans_ti = trans_ti, trans_to = trans_to,
                trans_ti_new = trans_ti_new, emit = emit,
                data_structure = data_str ) )
}

  