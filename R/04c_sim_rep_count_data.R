# Simulate count data under our population model from 1986 to 2021
library(tidyverse)

sim_rep_count_data <- function(fit, draw) {
  vital_rates <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb", "sW")
  initial_pop_stage_sizes <- c("Pb1[1]","Pb2[1]","Pb3[1]", "Pb4[1]", "Br[1]", "Nb[1]")

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
  
  # format posterior sample of vital rates
  vr_posterior_sample <- fit$draws(variables = vital_rates, format = "df") %>%
    as_tibble() %>%
    filter(.draw == draw) %>% 
    select(-c(".chain",".iteration",".draw")) %>%
    pivot_longer(everything()) %>%
    mutate(par = sapply(name, get_var_name), t = sapply(name,get_var_time), .after = name)
  
  # decant posterior sample values from tibble to vectors
  s0 <- filter(vr_posterior_sample,par=="s0")$value
  sN <- filter(vr_posterior_sample,par=="sN")$value
  sB <- filter(vr_posterior_sample,par=="sB")$value
  f3 <- filter(vr_posterior_sample,par=="f3")$value
  f4 <- filter(vr_posterior_sample,par=="f4")$value
  bb <- filter(vr_posterior_sample,par=="bb")$value
  nb <- filter(vr_posterior_sample,par=="nb")$value
  sW <- filter(vr_posterior_sample,par=="sW")$value
  
  # format posterior sample of initial stage sizes
  ipss_posterior_sample <- fit$draws(variables = initial_pop_stage_sizes, format = "df") %>%
    as_tibble() %>%
    filter(.draw == draw) %>% 
    select(-c(".chain",".iteration",".draw")) 
  
  # decant initial population stage size samples
  Pb1_1 <- floor(ipss_posterior_sample$`Pb1[1]`)
  Pb2_1 <- floor(ipss_posterior_sample$`Pb2[1]`)
  Pb3_1 <- floor(ipss_posterior_sample$`Pb3[1]`)
  Pb4_1 <- floor(ipss_posterior_sample$`Pb4[1]`)
  Br_1 <- floor(ipss_posterior_sample$`Br[1]`)
  Nb_1 <- floor(ipss_posterior_sample$`Nb[1]`)
  
  # decant number of immigrants in all years
  In_posterior_sample <- fit$draws(variables = "In", format = "df") %>%
    as_tibble() %>%
    filter(.draw == draw) %>% 
    select(-c(".chain",".iteration",".draw")) 
  
  In <- In_posterior_sample %>%
    pivot_longer( everything()) %>%
    pull(value) %>%
    floor()
  
  # Simulate population -------------------------------------------------------
  Tco <- length(In)

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
  
  # Simulate count data --------------------------------------------------------
  count <- rnorm(Tco, Br + In, 10)
  
  # outputs
  return( list(draw = draw, vital_rates = vr_posterior_sample, 
               init_pop_stage_sizes = ipss_posterior_sample,
               In = In_posterior_sample, count = count) )
}
