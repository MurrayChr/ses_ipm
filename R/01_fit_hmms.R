# Fit hidden markov model to capture-recapture data, with chosen priors
# on survival and state transition probabilities 
library(tidyverse)
library(cmdstanr)

# Read in data -----------------------------------------------------------------
data <- readRDS("data/cmr_data_shuffled_reduced.RDS")  
data <- data[[1]]

# Separate data components -----------------------------------------------------
# i.e. first capture, tag location, unique capture histories and their counts
fc <- data$fc
tagloc <- data$tagloc
count <- data$count
y <- data %>% 
  select(starts_with("yr")) %>%
  as.matrix()

# Create data list to pass to Stan ---------------------------------------------
data_list <- list(N = dim(y)[1],    # number of unique capture histories
                  T = dim(y)[2],    # number of years
                  K = 15,           # number of hidden states
                  L = 15,           # number of observable states 
                  y = y,            # unique capture histories
                  mult = count,     # multiplicity/count of each capture history
                  fc = fc,          # occasion of first capture
                  tagloc = tagloc)  # tag location covariate

# if fitting hmm_04 with 'cheap' gaussian process priors, uncomment following
# lines to add extra 'data' arguments with gaussian process hyperparameters

# hmm_fit <- readRDS("outputs/hmm_03_gaussian_process_fit.RDS")
# vr <- c("s0", "sN", "sB", "f3", "f4", "nb", "bb")    # vital rates
# vr_hyperpars <- c( str_c("ls_mean_", vr), str_c("sigma_mean_", vr),
#                    str_c("ls_sd_", vr), str_c("sigma_sd_", vr) )
# gp_hyperpars <- hmm_fit$summary(vr_hyperpars) %>%
#   select(variable, median) %>%
#   pivot_wider(names_from = variable, values_from = median) %>%
#   as.list()
# data_list <- c( data_list, gp_hyperpars )

# Fit model --------------------------------------------------------------------
# see https://mc-stan.org/cmdstanr/articles/cmdstanr.html

# path to chosen .stan model file
file <- "stan/hmm_01_fixed_effects_parallel.stan"
# file <- "stan/hmm_02_random_effects_parallel.stan"
# file <- "stan/hmm_03_gaussian_process_parallel.stan"
# file <- "stan/hmm_04_cheap_gp_parallel.stan"

# compile model with within-chain parallelisation enabled
mod <- cmdstan_model(file, cpp_options = list(stan_threads = TRUE))

# fit model (on 24 core machine)
fit <- mod$sample(data = data_list,
                  chains = 4,
                  parallel_chains = 4,    
                  threads_per_chain = 6,   
                  refresh = 1000)            # update progress every <refresh> iterations

# Save model output ------------------------------------------------------------
# fit$save_object(file = "outputs/hmm_01_fixed_effects_fit.RDS")
# fit$save_object(file = "outputs/hmm_02_random_effects_fit.RDS")
# fit$save_object(file = "outputs/hmm_03_gaussian_process_fit.RDS")
# fit$save_object(file = "outputs/hmm_04_cheap_gp_fit.RDS")

