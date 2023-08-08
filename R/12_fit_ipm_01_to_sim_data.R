# Fit ipm_01 to simulated data
# using "cheap" gaussian process priors on vital rates
library(tidyverse)
library(cmdstanr)

# Get simulated data -----------------------------------------------------------

# choose sim id
sim_id <- "low_1"

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

# Read in gp hyperparameters ---------------------------------------------
# import gp lengthscale and marginal standard deviation values
hmm_fit <- readRDS("outputs/hmm_03_gaussian_process_fit.RDS")
vr <- c("s0", "sN", "sB", "f3", "f4", "nb", "bb")    # vital rates
vr_hyperpars <- c( str_c("ls_mean_", vr), str_c("sigma_mean_", vr),
                   str_c("ls_sd_", vr), str_c("sigma_sd_", vr) )

gp_hyperpars <- hmm_fit$summary(vr_hyperpars) %>%
  select(variable, median) %>%
  pivot_wider(names_from = variable, values_from = median) %>%
  as.list()

# Read in paramaters for priors on initial pop stage sizes --------------------- 
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

# Create data list for Stan ----------------------------------------------------
data_list <- list(Tco = length(y_count),       # no. years count data
                  y_count = y_count,           # count data
                  Tmr = dim(y_cmr)[2],         # no. years cmr data
                  x = 1:(dim(y_cmr)[2] - 1),   # time covariate for gps
                  N = dim(y_cmr)[1],           # no. unique capture histories
                  K = 15,                      # no.hidden states
                  L = 15,                      # no. observable states 
                  y_cmr = y_cmr,               # unique capture histories
                  mult = mult,                 # multiplicity of each capture history
                  fc = fc,                     # occasion of first capture
                  tagloc = tagloc)             # tag location covariate

# combined data list
data_list <- c( data_list, gp_hyperpars, prior_init_pars )

# Fit model --------------------------------------------------------------------
# see https://mc-stan.org/cmdstanr/articles/cmdstanr.html

# path to .stan model file
file <- "stan/ipm_01.stan"

# compile model with within-chain parallelisation enabled
mod <- cmdstan_model(file, cpp_options = list(stan_threads = TRUE))

# fit model (using 24 cores)
fit <- mod$sample(data = data_list,
                  chains = 4,
                  parallel_chains = 4,
                  threads_per_chain = 6,
                  refresh = 500)            # update progress every <refresh> iterations

# Save model output ------------------------------------------------------------
# fit$save_object( str_c("outputs/ipm_01_fit_",sim_id,".RDS" ) )
