# Inspect model diagnostics and posterior estimates from hmm
library(tidyverse)
library(cmdstanr)
library(bayesplot)

# Select hmm model to inspect --------------------------------------------------
# file <- "hmm_01_fixed_effects_fit.RDS"
file <- "hmm_02_random_effects_fit.RDS"
# file <- "hmm_03_gaussian_process_fit.RDS"
# file <- "hmm_04_cheap_gp_fit.RDS"
fit <- readRDS(str_c("outputs/", file))

# MCMC diagnostics -------------------------------------------------------------
fit$diagnostic_summary()

# Effective sample size, Rhat etc ----------------------------------------------

# Choose a parameter e.g. "s0"
fit$summary(variables = "s0")

# Functions to plot posterior estimates ----------------------------------------

plot_vital_rate <- function(vr, fit) {
  vital_rates <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb")
  if (!(vr %in% vital_rates)) {
    return( str_c( "Error! vital rate must be one of: ", 
                   str_c(vital_rates, collapse = ", "),".") )
  }
  fit$summary(variables = vr) %>%
    add_column(year = 1983:2020) %>%
    ggplot( ) +
    geom_pointrange( aes(x = year, y = median, ymin = q5, ymax = q95)) +
    coord_cartesian( ylim = c(0,1)) +
    labs(y = "posterior median and 90% CrI",
         title = str_c("Estimate for ", vr, " from ", str_sub(file,1,6)))
}

plot_detection <- function(det, fit) {
  det_pars <- c("qN", "qB", "pBu", "pBe")
  if (!(det %in% det_pars)) {
    return( str_c( "Error! detection must be one of: ", 
                   str_c(det_pars, collapse = ", "),".") )
  }
  fit$summary(variables = det) %>%
    add_column(year = 1983:2021) %>%
    ggplot( ) +
    geom_pointrange( aes(x = year, y = median, ymin = q5, ymax = q95)) +
    coord_cartesian( ylim = c(0,1)) +
    labs(y = "posterior median and 90% CrI",
         title = str_c("Estimate for ", det, " from ", str_sub(file,1,6)))
}

# Plot posterior estimates of vital rates or detection parameters --------------

# vital rates
vr <- "s0"   # choose one of "s0", "sN", "sB", "f3", "f4", "bb", "nb"
plot_vital_rate(vr, fit)

# detection probabilities
det <- "qN"  # choose one of "qN", "qB", "pBu", "pBe"
plot_detection(det, fit)

# Plot estimates of tag loss probabilities -------------------------------------
fit$draws(variables = c("ti", "to", "ti_new")) %>%
  mcmc_hist() + 
  labs( title= "Posterior estimates for tag loss probabilities" )



