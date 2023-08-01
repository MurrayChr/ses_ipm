# Function to compute expected breeder counts
library(tidyverse)

get_expected_breeder_count <- function(rep_data) {
  # decant values for vital rates, initial population stage sizes and immigration 
  s0 <- filter(rep_data$vital_rates,par=="s0")$value
  sN <- filter(rep_data$vital_rates,par=="sN")$value
  sB <- filter(rep_data$vital_rates,par=="sB")$value
  f3 <- filter(rep_data$vital_rates,par=="f3")$value
  f4 <- filter(rep_data$vital_rates,par=="f4")$value
  bb <- filter(rep_data$vital_rates,par=="bb")$value
  nb <- filter(rep_data$vital_rates,par=="nb")$value
  sW <- filter(rep_data$vital_rates,par=="sW")$value
  
  Pb1_1 <- rep_data$init_pop_stage_sizes$`Pb1[1]`
  Pb2_1 <- rep_data$init_pop_stage_sizes$`Pb2[1]`
  Pb3_1 <- rep_data$init_pop_stage_sizes$`Pb3[1]`
  Pb4_1 <- rep_data$init_pop_stage_sizes$`Pb4[1]`
  Br_1 <- rep_data$init_pop_stage_sizes$`Br[1]`
  Nb_1 <- rep_data$init_pop_stage_sizes$`Nb[1]`
  
  In <- rep_data$In %>%
    pivot_longer( everything()) %>%
    pull(value) %>%
    floor()
  
  # Calculate expected values using projection matrix
  Tco <- length(In)
  pop_vec <- list()  # list of vectors of population stage sizes
  proj_mat <- list() # list of projection matrices
  for ( t in 1:(Tco - 1) ) {
    # index for the vital rates s0, ..., nb but not sW
    tmr <- t + 3
    proj_mat[[t]] <- matrix(NA, 6, 7)
    proj_mat[[t]][1,] <- c(      0,                   0,                   0,       0,   0.5*sW[t]*s0[tmr],                   0,   0.5*sW[t]*s0[tmr])
    proj_mat[[t]][2,] <- c(sN[tmr],                   0,                   0,       0,                   0,                   0,                   0)
    proj_mat[[t]][3,] <- c(      0, sN[tmr]*(1-f3[tmr]),                   0,       0,                   0,                   0,                   0)
    proj_mat[[t]][4,] <- c(      0,                   0, sN[tmr]*(1-f4[tmr]),       0,                   0,                   0,                   0)
    proj_mat[[t]][5,] <- c(      0,     sN[tmr]*f3[tmr],     sN[tmr]*f4[tmr], sN[tmr],     sB[tmr]*bb[tmr],     sN[tmr]*nb[tmr],     sB[tmr]*bb[tmr])
    proj_mat[[t]][6,] <- c(      0,                   0,                   0,       0, sB[tmr]*(1-bb[tmr]), sN[tmr]*(1-nb[tmr]), sB[tmr]*(1-bb[tmr]))
  }
  temp_pop_vec <- rep_data$init_pop_stage_sizes %>%
    pivot_longer( everything() ) %>%
    pull(value)
  pop_vec[[1]] <- matrix( c(temp_pop_vec, In[1]), ncol = 1 )
  for (t in 1:(Tco - 1)) {
    temp_pop_vec <- proj_mat[[t]] %*% pop_vec[[t]]
    pop_vec[[t+1]] <- matrix( c(temp_pop_vec, In[t+1]), ncol = 1 )
  }
  
  # expected count 
  #' equals the sum of expected number of Br and In, which are entries 
  #' 5 and 7 in the pop_vec vector
  expected_count <- rep(NA, Tco)
  for (t in 1:Tco) {
    expected_count[t] <- sum( pop_vec[[t]][c(5,7)] )
  }
  
  return( expected_count)
}

