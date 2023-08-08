# Fit population models to count data 
# using informative priors on vital rates given by a multivariate normal 
# approximation to the logit-transformed posterior of an HMM
library(tidyverse)
library(cmdstanr)

# read in selected hmm model fit -----------------------------------------------
hmm_fit <- readRDS("outputs/hmm_01_fixed_effects_fit.RDS")
# hmm_fit <- readRDS("outputs/hmm_02_random_effects_fit.RDS")
# hmm_fit <- readRDS("outputs/hmm_03_gaussian_process_fit.RDS")
# hmm_fit <- readRDS("outputs/hmm_04_cheap_gp_fit.RDS")

# match multivariate normal to logit-transformed posterior draws ---------------
# extract posterior draws of vital rates
hmm_posterior <- 
  hmm_fit$draws(variables = c("s0", "sN", "sB", "f3", "f4", "bb", "nb"),
                format = "df") %>%
  select( -starts_with("."))  # removes .chain, .iteration, .draws

# remove first three years from each parameter
# as count data is only available from 1986 onwards
hmm_posterior <- hmm_posterior %>%
  select( - contains( str_c("[",1:3,"]") ) ) 

# map to logit scale
logit <- function(p){
  eps <- 10^(-9)
  if (p == 1) {
    p <- 1 - eps
  }
  if (p == 0) {
    p <- eps
  }
  log(p) - log(1-p) 
}

logit_hmm_posterior <- apply(hmm_posterior, c(1,2), logit)

# compute mean and covariance
mean_logit_vr <- apply(logit_hmm_posterior,2,mean)
cov_logit_vr <- cov(logit_hmm_posterior)

# set indices for each variable s0, sN, ...
inds_s0 <- which( str_detect(names(mean_logit_vr), "s0") )
inds_sN <- which( str_detect(names(mean_logit_vr), "sN") )
inds_sB <- which( str_detect(names(mean_logit_vr), "sB") )
inds_f3 <- which( str_detect(names(mean_logit_vr), "f3") )
inds_f4 <- which( str_detect(names(mean_logit_vr), "f4") )
inds_nb <- which( str_detect(names(mean_logit_vr), "nb") )
inds_bb <- which( str_detect(names(mean_logit_vr), "bb") )
inds_list <- list(inds_s0 = inds_s0, inds_sN = inds_sN, inds_sB = inds_sB,
                  inds_f3 = inds_f3, inds_f4 = inds_f4, inds_bb = inds_bb,
                  inds_nb = inds_nb)

# Read in parameters for priors on intiial stage sizes -------------------------
prior_init <- readRDS("outputs/prior_stage_sizes_1986.RDS")

# choose omega value
omega_val <- 0.1

prior_init_pars <- list()
for ( p in c("Pb1", "Pb2", "Pb3", "Pb4", "Br", "Nb") ) {
  p_prior_inits <- prior_init %>% 
    filter( omega ==  omega_val) %>%
    filter(par_name == p)
  prior_init_pars[[ str_c(p,"_1_pars") ]] <- 
    c(p_prior_inits$mean, p_prior_inits$sd)
}

# Fit population model with cmdstanr -------------------------------------------
# read in 15 October count data
count_data <- readRDS("data/count_data_15oct.RDS")
count_data <- count_data$count

# select pop model
file <- "stan/pop_01_logit_mvn.stan"   # 'random effects' prior on immigration
# file <- "stan/pop_02_gp_in.stan"   # gaussian process prior on immigration

# compile 
mod <- cmdstan_model(file)

# input data
input_data <- c( list( T = length(count_data), 
                       y_count = count_data,
                       n_pars = length(mean_logit_vr),
                       mean_logit_vr = mean_logit_vr, 
                       cov_logit_vr = cov_logit_vr ), 
                 inds_list,
                 prior_init_pars )

# if using pop_02 with gaussian process prior on immigration, uncomment following line
# to add a required time covariate to the input data
# input_data <- c( input_data, list(x = 1:length(count_data)))

# fit
fit <- mod$sample(data = input_data, 
                  chains = 4, 
                  parallel_chains = 4,
                  adapt_delta = 0.9)

# Save model output ------------------------------------------------------------
# fit$save_object("outputs/pop_01_prior_hmm_02.RDS")
# fit$save_object("outputs/pop_01_prior_hmm_03.RDS")
# fit$save_object("outputs/pop_01_prior_hmm_04.RDS")
# fit$save_object("outputs/pop_02_prior_hmm_04.RDS")
