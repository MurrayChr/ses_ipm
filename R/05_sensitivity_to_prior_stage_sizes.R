# Examine the sensitivity of 1990 stage size estimates to the priors used for
# the initial stage sizes in 1986.
library(tidyverse)
library(cmdstanr)

# We use the iterative procedure to fit IPM_GP_RE repeatedly, with different
# priors on the initial stage sizes in 1986, defined by different assumed immigration
# rates, omega.

# Read in multievent model fit and format posterior for use as prior -----------
hmm_fit <- readRDS("outputs/hmm_04_cheap_gp_fit.RDS")

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

# Fit models with different priors on initial stage sizes ----------------------

# read in parameters for priors on initial stage sizes
prior_init <- readRDS("outputs/prior_stage_sizes_1986.RDS")

# loop over assumed immigration rate omega, fitting the model to count data
# initialise list to store fits
pop_01_hmm_04_fits <- list()

# values for omega
omega_vals <- unique(prior_init$omega)

# count data
count_data <- readRDS("data/count_data_15oct.RDS")
count_data <- count_data$count

# compile stan model
file <- "stan/pop_01_logit_mvn.stan"
mod <- cmdstan_model(file)

for ( i in 1:length(omega_vals) ) {

  # decant prior parameters for initial stage size distributions
  prior_init_pars <- list()
  for ( p in c("Pb1", "Pb2", "Pb3", "Pb4", "Br", "Nb") ) {
    p_prior_inits <- prior_init %>% 
      filter( omega ==  omega_vals[i]) %>%
      filter(par_name == p)
    prior_init_pars[[ str_c(p,"_1_pars") ]] <- 
      c(p_prior_inits$mean, p_prior_inits$sd)
  }
  
  # prepare input data
  input_data <- c( list( T = length(count_data), 
                         y_count = count_data,
                         n_pars = length(mean_logit_vr),
                         mean_logit_vr = mean_logit_vr, 
                         cov_logit_vr = cov_logit_vr ), 
                   inds_list,
                   prior_init_pars )
  
  # fit model and save to list
  pop_01_hmm_04_fits[[i]] <- mod$sample( data = input_data,
                    chains = 4, parallel_chains = 4,
                    adapt_delta = 0.9, refresh = 10 )
}

# save model fits
# saveRDS(pop_01_hmm_04_fits, "outputs/prior_stage_size_sensitivity_pop_01_hmm_04_fits.RDS")

# Look at posterior distributions of stage sizes in 1990 -----------------------

post_samples <- tibble()
var_names <- c("Pb1", "Pb2", "Pb3", "Pb4", "Nb", "Br", "In")
t_1990 <- which(1986:2021 == 1990)
vars <- str_c(var_names, "[", t_1990, "]")
for ( i in 1:length(omega_vals) ) {
  fit <- pop_01_hmm_04_fits[[i]]
  temp <- fit$draws( variables = vars, format = "df" ) %>%
    select( -starts_with(".") ) %>%
    pivot_longer( everything() ) %>%
    add_column( omega = omega_vals[i] )
  post_samples <- rbind(post_samples, temp)
}

# save 1990 posterior samples
# saveRDS(post_samples, "outputs/prior_stage_size_sensitivity_1990_posteriors.RDS")

# plot
facet_names = c("Pb1[5]" = "Pb1", "Pb2[5]" = "Pb2", "Pb3[5]" = "Pb3", "Pb4[5]" = "Pb4",
                "Br[5]" = "Br", "Nb[5]" = "Nb", "In[5]" = "In")
post_samples %>%
  ggplot() +
  facet_wrap( "name" , scales = "free", nrow = 2,
              labeller = as_labeller(facet_names) ) +
  stat_density(aes(x = value, colour = as.factor(omega) ),  # using stat density cf. geom_density to get lines in legend
               geom="line", position="identity") +
  labs( colour = "omega" ) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = c(0.9,0.2))







