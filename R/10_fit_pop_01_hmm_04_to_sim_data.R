#' Fit ipm to simulated data with two-stage fit
#' First, fit hmm_04 ('cheap gp') to the cmr data, then fit pop_01 to the count 
#' data with priors on the vital rates given by a multi-normal approximation to
#' the posterior from the multi-event model  
library(tidyverse)
library(cmdstanr)

# Get simulated data -----------------------------------------------------------

# choose sim id
sim_id <- "low_001"

# read in simulated cmr and count data
sim_data <- readRDS( str_c("data/sim_",sim_id,".RDS") )
cmr_data <- sim_data$cmr_data
count_data <- sim_data$count_data

# prepare cmr data by separating data components
# i.e. first capture, tag location, unique capture histories and their counts
fc <- cmr_data$fc
tagloc <- cmr_data$tagloc
mult <- cmr_data$count
y_cmr <- cmr_data %>% 
  select(starts_with("yr")) %>%
  as.matrix()

# prepare count data
y_count <- count_data$count   # counts without the years

# Fit hmm_04 to the cmr data ---------------------------------------------------

# import gp lengthscale and marginal standard deviation values
full_gp_fit <- readRDS("outputs/hmm_03_gaussian_process_fit.RDS")
vr <- c("s0", "sN", "sB", "f3", "f4", "nb", "bb")    # vital rates
vr_hyperpars <- c( str_c("ls_mean_", vr), str_c("sigma_mean_", vr),
                   str_c("ls_sd_", vr), str_c("sigma_sd_", vr) )

gp_hyperpars <- full_gp_fit$summary(vr_hyperpars) %>%
  select(variable, median) %>%
  pivot_wider(names_from = variable, values_from = median) %>%
  as.list()

# list of all other data
data_list <- list(N = nrow(y_cmr),    # number of unique capture histories
                  T = ncol(y_cmr),    # number of years
                  x = 1:(ncol(y_cmr) - 1),       # time covariate for gp models
                  K = 15,           # number of hidden states
                  L = 15,           # number of observable states 
                  y = y_cmr,            # unique capture histories
                  mult = mult,     # multiplicity/count of each capture history
                  fc = fc,          # occasion of first capture
                  tagloc = tagloc)  # tag location covariate

# combined data list
data_list <- c( data_list, gp_hyperpars )

# path to .stan model file
file <- "stan/hmm_04_cheap_gp_parallel.stan"

# compile model with within-chain parallelisation enabled
mod <- cmdstan_model(file, cpp_options = list(stan_threads = TRUE))

# fit model (with 24 cores)
fit <- mod$sample(data = data_list,
                  chains = 4,
                  parallel_chains = 4,    
                  threads_per_chain = 6,   
                  refresh = 500)            

# save model output
fit$save_object( str_c("outputs/hmm_04_fit_",sim_id,".RDS" ) )

# Fit pop_01 model to count data with informative vital rate prior -------------

hmm_fit <- readRDS( str_c("outputs/hmm_04_fit_",sim_id,".RDS" ) )

# match multivariate normal to logit-transformed posterior draws 
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

# Read in parameters for priors on intiial stage sizes 
prior_init <- readRDS("outputs/prior_stage_sizes_1986.RDS")

# choose omega value
omega_val <- 0.1

prior_init_pars <- list()
for ( p in c("Pb1", "Pb2", "Pb3", "Pb4", "Br", "Nb") ) {
  p_prior_inits <- prior_init %>% 
    filter( omega ==  omega_val) %>%
    filter(par_name == p)
  prior_init_pars[[ str_c(p,"_1_pars") ]] <- 
    # c(p_prior_inits$mean, p_prior_inits$sd)
    c(p_prior_inits$mean, p_prior_inits$sd)
}

# fit pop_01 model
file <- "stan/pop_01_logit_mvn.stan"
mod <- cmdstan_model(file)
input_data <- c( list( T = length(y_count), 
                       y_count = y_count,
                       n_pars = length(mean_logit_vr),
                       mean_logit_vr = mean_logit_vr, 
                       cov_logit_vr = cov_logit_vr ), 
                 inds_list,
                 prior_init_pars )
fit <- mod$sample(data = input_data,
                  chains = 4,
                  parallel_chains = 4,
                  adapt_delta = 0.9,
                  refresh = 500)

# save model output 
# fit$save_object( str_c("outputs/pop_01_hmm_04_fit_",sim_id,".RDS" ) )