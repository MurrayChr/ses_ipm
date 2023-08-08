# Posterior predictive checks for the breeder count data
library(tidyverse)
source("R/04c_sim_rep_count_data.R")
source("R/04d_expected_count_function.R")

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

## select ipm model
# file <- "ipm_01_fit.RDS"

## prefix for naming files
# prefix <- str_sub(file, 1, 6)

## simulate and save replicate datasets
# fit <- readRDS( str_c( "outputs/", file ) )
# n_reps <- 100
# draws <- sample(1:4000, n_reps, replace = FALSE)
# for (n in 1:n_reps) {
#   rep_data <- sim_rep_count_data(fit, draws[n])
#   out_file <- str_c( prefix, "_rep_", file_number_helper(n), ".RDS")
#   saveRDS( rep_data, str_c("ppc_rep_data/", out_file ) )
# }

# Read in replicate count data -------------------------------------------------

rep_count_data <- list()
n_reps <- 100
for (n in 1:n_reps){
  rep_count_data[[n]] <- readRDS( str_c("ppc_rep_data/ipm_01_rep_", file_number_helper(n),".RDS") )
}
  
# Plot the replicate and real count data ---------------------------------------
n_reps <- 100
counts <- tibble()
for (n in 1:n_reps) {
  temp <- tibble(rep_num = n, year = 1986:2021, 
                 count = rep_count_data[[n]]$count)
  counts <- rbind( counts, temp )
}
# save for plotting later
# saveRDS(counts, "outputs/ppc_ipm_01_rep_counts.RDS")

# read in real counts
real_count <- readRDS("data/count_data_15oct.RDS")

counts %>%
  ggplot(  ) +
  geom_line( aes(x = year, y = count, group = rep_num),
             size = 0.2, alpha = 0.8, colour = "grey" ) +
  geom_line( data = real_count, aes(x = year, y = count), colour = "navyblue", size = 1.5) +
  labs( y = "number of breeding seals") +
  theme_classic( ) 

# Compute a Bayesian p-value ---------------------------------------------------
# function to compute Freeman-Tukey statistic
get_ft <- function(observed_counts, expected_counts){
  sum( ( sqrt(observed_counts) - sqrt(expected_counts) )^2 )
}

# read in real counts
real_count <- readRDS("data/count_data_15oct.RDS")

# compute tibble of pairs FT(y_rep[i], par[i]), FT(y, par[i])
ft_stats <- tibble(y_rep = rep(NA, n_reps), y = rep(NA, n_reps))
for (n in 1:n_reps) {
  rep_data <- rep_count_data[[n]]
  expected_count <- get_expected_breeder_count(rep_data)
  observed_count_rep_data <- rep_data$count
  observed_count_real_data <- real_count$count
  ft_stats$y_rep[n] <- get_ft(observed_count_rep_data, expected_count)
  ft_stats$y[n] <- get_ft(observed_count_real_data, expected_count)
}
# save for plotting later
# saveRDS(ft_stats, "outputs/ppc_ipm_01_freeman_tukey_breeder_counts.RDS")

pB <- ft_stats %>%             # compute Bayesian p-value
  transmute(bool = y_rep > y) %>%
  summarise(mean(bool))

# find good plot limits
plot_lim <- ft_stats %>%
  summarise( min = min(y_rep, y),
             max = max(y_rep, y),
             diam = max(y_rep, y) -  min(y_rep, y))

ft_stats %>%
  ggplot(aes(x=y,y=y_rep)) +
  geom_point( size = 2, shape = 21) + 
  geom_abline( slope = 1, intercept = 0, size = 0.5, alpha = 0.2) +
  coord_fixed( xlim = with( plot_lim, c(0, max + 0.05*diam) ), 
               ylim = with( plot_lim, c(0, max + 0.05*diam) ) ) +
  annotate("text", x = 42, y = 28, size = 5,
           label = str_c("p[B] == ",pB), parse = TRUE) +
  labs( x = "real data discrepancy", y = "replicate data discrepancy") +
  theme_classic( ) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 12) )
