# Sim ipm data and fit model
# One must set the simulation id, and specify the level of variation in
# immigration
library(tidyverse)
library(posterior)
library(cmdstanr)
source("R/11_sim_ipm_data_function.R")

# import parameter values
ipm_fit <- readRDS("outputs/ipm_01_fit.RDS")

# randomly draw posterior sample
post_sample <- sample(1:4000,1)

# decant posterior values
variables <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb", "sW", "qN", "qB", "pBu", "pBe", 
               "Pb1[1]", "Pb2[1]", "Pb3[1]", "Pb4[1]", "Br[1]","Nb[1]", 
               "ti", "to", "ti_new", 
               "mean_log_lambda", "sd_log_lambda")
post_vals <- ipm_fit$draws(variables, format = "df") %>%
  filter( .draw == post_sample ) %>%
  select( - starts_with(".") )

s0 <- select(post_vals, starts_with("s0") ) %>% as.numeric()
sN <- select(post_vals, starts_with("sN") ) %>% as.numeric()
sB <- select(post_vals, starts_with("sB") ) %>% as.numeric()
f3 <- select(post_vals, starts_with("f3") ) %>% as.numeric()
f4 <- select(post_vals, starts_with("f4") ) %>% as.numeric()
bb <- select(post_vals, starts_with("bb") ) %>% as.numeric()
nb <- select(post_vals, starts_with("nb", ignore.case = FALSE) ) %>% as.numeric()
sW <- select(post_vals, starts_with("sW") ) %>% as.numeric()
qN <- select(post_vals, starts_with("qN") ) %>% as.numeric()
qB <- select(post_vals, starts_with("qB") ) %>% as.numeric()
pBu <- select(post_vals, starts_with("pBu") ) %>% as.numeric()
pBe <- select(post_vals, starts_with("pBe") ) %>% as.numeric()
ti <- select(post_vals, "ti" ) %>% as.numeric()
to <- select(post_vals, "to" ) %>% as.numeric()
ti_new <- select(post_vals, "ti_new" ) %>% as.numeric()
Pb1_1 <- select(post_vals, "Pb1[1]" ) %>% as.numeric() %>% floor()
Pb2_1 <- select(post_vals, "Pb2[1]" ) %>% as.numeric() %>% floor()
Pb3_1 <- select(post_vals, "Pb3[1]" ) %>% as.numeric() %>% floor()
Pb4_1 <- select(post_vals, "Pb4[1]" ) %>% as.numeric() %>% floor()
Br_1 <- select(post_vals, "Br[1]" ) %>% as.numeric() %>% floor()
Nb_1 <- select(post_vals, "Nb[1]" ) %>% as.numeric() %>% floor()
mean_log_lambda <- select(post_vals, "mean_log_lambda" ) %>% as.numeric()

# select level of variation in immigration
sd_log_lambda <- 0

# simulation id
sim_id <- str_c("low_001")

# simulate ipm data and save output to 'data' sub-directory
sim_ipm_data(s0=s0, sN=sN, sB=sB, f3=f3, f4=f4, bb=bb, nb=nb, sW=sW, 
             qN=qN, qB=qB, pBu=pBu, pBe=pBe,
             ti=ti, to=to, ti_new=ti_new, 
             Pb1_1=Pb1_1, Pb2_1=Pb2_1, Pb3_1=Pb3_1, Pb4_1=Pb4_1, Br_1=Br_1, Nb_1=Nb_1,
             mean_log_lambda=mean_log_lambda, sd_log_lambda=sd_log_lambda, 
             sim_id=sim_id )
