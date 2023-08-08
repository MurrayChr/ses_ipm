# Inspect estimates from ipms
library(tidyverse)
library(cmdstanr)
library(bayesplot)

# Select ipm model to inspect --------------------------------------------------
file <- "ipm_01_fit.RDS"     # IPM_GP_RE
# file <- "ipm_02_fit.RDS"     # IPM_RE_RE
# file <- "ipm_03_fit.RDS"     # IPM_GP_GP
fit <- readRDS(str_c("outputs/", file))

# MCMC diagnostics -------------------------------------------------------------
fit$diagnostic_summary()

# Effective sample size, Rhat etc ----------------------------------------------

# Choose a parameter e.g. "Br", "s0", etc
fit$summary(variables = "In")

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

plot_stage_size <- function(st, fit) {
  stage <- c("Pb1", "Pb2", "Pb3", "Pb4", "Br", "Nb", "In")
  if (!(st %in% stage)) {
    return( str_c( "Error! stage size must be one of: ", 
                   str_c(stage, collapse = ", "),".") )
  }
  fit$summary(variables = st) %>%
    add_column(year = 1986:2021) %>%
    ggplot( ) +
    geom_pointrange( aes(x = year, y = median, ymin = q5, ymax = q95)) +
    coord_cartesian( ylim = c(0,NA)) +
    labs(y = "posterior median and 90% CrI",
         title = str_c("Estimate for ", st, " from ", str_sub(file,1,6)))
}

# Plot posterior estimates of vital rates or detection parameters --------------

# vital rates
vr <- "s0"   # choose one of "s0", "sN", "sB", "f3", "f4", "bb", "nb"
plot_vital_rate(vr, fit)

# detection probabilities
det <- "qN"  # choose one of "qN", "qB", "pBu", "pBe"
plot_detection(det, fit)

# Plot population stage size ---------------------------------------------------

st <- "Br"  # choose one of "Pb1", "Pb2", "Pb3", "Pb4", "Br", "Nb", "In"
plot_stage_size(st, fit) 

# Immigration hyperparameter posteriors vs priors for ipm_GP_RE ----------------
#  priors are (see ipm_01.stan, lines 488, 489)
#  mean_log_lambda ~ normal(log(50), 0.3)
#  sd_log_lambda prior ~ halfnormal(0, 0.4)

prior_vs_posterior <- function(post_draws, prior_draws) {
  pr_tib <- as_tibble_col(prior_draws, column_name = "prior") %>%
    pivot_longer(everything())
  po_tib <- as_tibble_col(post_draws, column_name = "posterior") %>%
    pivot_longer(everything())
  rbind( po_tib, pr_tib )  %>%
    ggplot( aes( x = value, fill = name) ) +
    geom_density( alpha = 0.5, colour = NA) 
}

# mean_log_lambda
post_draws <- fit$draws(variables = "mean_log_lambda", format = "df") %>%
  select( -starts_with(".") ) %>%
  pull()

prior_draws <- rnorm(length(post_draws), log(50), 0.3)

prior_vs_posterior(post_draws, prior_draws) +
  labs( title = "mean_log_lambda",
        subtitle = "ipm_01")

# sd_log_lambda
post_draws <- fit$draws(variables = "sd_log_lambda", format = "df") %>%
  select( -starts_with(".") ) %>%
  pull()

# draw from half-normal by drawing from normal and discarding negative samples
prior_draws <- rnorm( 10*length(post_draws), 0, 0.4)
prior_draws <- prior_draws[prior_draws > 0][1:length(post_draws)]

prior_vs_posterior(post_draws, prior_draws) +
  labs( title = "sd_log_lambda",
        subtitle = "ipm_01")


