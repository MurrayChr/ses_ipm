# Goodness of fit via posterior predictive checks
library(tidyverse)
source("R/04a_sim_rep_cmr_data.R")
source("R/04b_detection_functions.R")

# Utility function for naming and reading in files -----------------------------
# adds zeros to have, e.g. '001' not '1' in name so that files are listed neatly
file_number_helper <- function(n, num_digits = 3){
  str <- as.character(n)
  if (str_length(str) < num_digits){
    str <- str_c(str_dup("0",num_digits-str_length(str)), str)
  }
  str
}

# Simulate replicate datasets --------------------------------------------------
# Uncomment this section if you want to re-simulate replicate datasets
# The replicate datasets are in /ppc_rep_data

## select hmm model
# file <- "hmm_01_fixed_effects_fit_chpc.RDS"
# file <- "hmm_02_random_effects_fit_chpc.RDS"
# file <- "hmm_03_gaussian_process_fit_chpc.RDS"
# file <- "hmm_04_cheap_gp_fit_chpc.RDS"

## prefix for naming files
# prefix <- str_sub(file, 1, 6)

## simulate and save replicate datasets
# fit <- readRDS( str_c( "outputs/", file ) )
# n_reps <- 100
# draws <- sample(1:4000, n_reps, replace = FALSE)
# for (n in 1:n_reps) {
#   rep_data <- sim_rep_cmr_data(fit, draws[n])
#   out_file <- str_c( prefix, "_rep_", file_number_helper(n), ".RDS")
#   saveRDS( rep_data, str_c("ppc_rep_data/", out_file ) )
# }

# Compute expected and observed number of detections for simulated data --------

## select hmm model
# model_file <- "hmm_01_fixed_effects_fit.RDS"
# model_file <- "hmm_02_random_effects_fit.RDS"
# model_file <- "hmm_03_gaussian_process_fit.RDS"
model_file <- "hmm_04_cheap_gp_fit.RDS"

model_prefix <- str_sub( model_file, 1,6 )

# import simulated replicate datasets
rep_data <- list()
n_reps <- 100
file_prefix <- str_c("ppc_rep_data/",model_prefix,"_rep_")
for (n in 1:n_reps){
  file <- str_c(file_prefix, file_number_helper(n),".RDS")
  rep_data[[n]] <- readRDS(file)
}

# import real data and format for get_observed_detections()
data <- readRDS("data/cmr_data.RDS")
real_data <- list()
real_data$y <- data$data %>%
  select( starts_with("yr") ) %>%
  as.matrix()
real_data$data_structure <- readRDS("data/cmr_dataset_structure.RDS")

# create tibble of expected and observed number of detections of both types ("all", "breeder")
get_n_det <- function(model_prefix, rep, obs_exp, type, data, num_years, real = FALSE) {
  if ( real ) {
    count <- get_observed_detections(data, type)
  }
  else {
    if ( obs_exp == "observed" ) {
      count <- get_observed_detections(data[[rep]], type)
    }
    if ( obs_exp == "expected" ) {
      count <- get_expected_detections(data[[rep]], type)
    }
  }
  tibble( model = model_prefix,
          rep = rep(rep, num_years),
          type = !!(type),
          obs_exp = !!(obs_exp),
          year = 1983:2021,
          count = count )
}

num_years <- ncol(real_data$y)
num_detections <- tibble()
for (n in 1:n_reps) {
  for (t in c("all", "breeder")) {
    for (oe in c("observed", "expected")) {
      temp <- get_n_det(model_prefix, rep = n, obs_exp = oe, 
                                        type = t, rep_data, num_years)
      num_detections <- rbind(num_detections, temp)
    }
  }
}

# calculate observed detections for real data
num_detections_real <- tibble()
for (t in c("all", "breeder")) {
  temp <- get_n_det(model_prefix=NA, rep=NA, obs_exp="observed", 
                                    type=t, real_data, num_years, real = TRUE)
  num_detections_real <- rbind(num_detections_real, temp)
}

# save num_detections and num_detections_real for plotting later
# saveRDS(num_detections, str_c("outputs/ppc_",model_prefix ,"_rep_data_detections.RDS"))
# saveRDS(num_detections_real, str_c("outputs/ppc_real_data_detections.RDS"))

# Graphical ppc for all detections and breeder detections ----------------------
# type <- "all"
type <- "breeder"
num_detections %>%
  filter( obs_exp == "observed", type == !!(type) ) %>%
  ggplot( aes(x = year, y = count, group = rep) ) +
  geom_line( size = 0.2, alpha = 0.8, colour = "grey" ) +
  geom_line( data = filter(num_detections_real, type == !!(type)),
             aes(x = year, y = count), colour = "navyblue", size = 1.5) +
  labs( y = "number of seals detected") +
  theme_classic( ) 

# Compute Freeman-Tukey statistic for all detections and breeder detections ----

# function to compute Freeman-Tukey statistic
get_ft <- function(observed_counts, expected_counts){
  sum( ( sqrt(observed_counts) - sqrt(expected_counts) )^2 )
}

# compute tibble of pairs FT(y_rep[i], par[i]), FT(y, par[i])
ft_stats <- tibble()
for (type in c("all", "breeder")) {
  temp <- tibble(type = !!(type), y_rep = rep(NA, n_reps), y =  rep(NA, n_reps))
  for (n in 1:n_reps) {
    expected_n_det <- num_detections %>%
      filter(rep == n, type == !!(type), obs_exp == "expected") %>%
      pull(count)
    rep_obs_n_det <- num_detections %>%
      filter(rep == n, type == !!(type), obs_exp == "observed") %>%
      pull(count)
    real_obs_n_det <- num_detections_real %>%
      filter( type == !!(type) ) %>%
      pull(count)
    temp$y_rep[n] <- get_ft(expected_n_det, rep_obs_n_det)
    temp$y[n] <- get_ft(expected_n_det, real_obs_n_det)
  }
  ft_stats <- rbind(ft_stats, temp)
}

# save for plotting later
# saveRDS(ft_stats, str_c("outputs/ppc_",model_prefix ,"_freeman_tukey_detections.RDS"))

# Compute Bayesian p-values and plot -------------------------------------------
# type <- "all"
type <- "breeder"

# compute Bayesian p-value
pB <- ft_stats %>%
  filter(type == !!(type)) %>%
  transmute(bool = y_rep > y) %>%
  summarise(mean(bool))

# plot
ft_stats %>%
  filter(type == !!(type)) %>%
  ggplot(aes(x=y,y=y_rep)) +
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, size = 0.3, alpha = 0.5) +
  coord_fixed(ratio = 1, xlim = c(0,NA), ylim = c(0,NA)) +
  labs(x="real data discrepancy", y = "replicate data discrepancy") +
  annotate("text", x = 12, y = 15, size = 4,
           label = str_c("p[B] == ",pB), parse = TRUE) +
  labs( title = str_c("FT discrepancy for ", 
                      type, " detections."),
        subtitle = str_c( model_file ) ) +
  theme_classic() 




