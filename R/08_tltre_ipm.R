# tLTRE analysis from IPM estimates, where the IPM has been fit using the 
# 'full' rather than iterative, two-stage procedure.
library(cmdstanr)
library(posterior)
library(tidyverse)

# extract posterior samples from model -----------------------------------------

# choose fitted model
file <- "ipm_01_fit.RDS"
# file <- "ipm_02_fit.RDS"
# file <- "ipm_02_fit.RDS"
fit <- readRDS(str_c("outputs/",file))

# parameters needed for tltre
vr_pars <- c("s0", "sN", "sB", "sW")       # survival probabilities 
pop_pars <- c("pNb", "pPb", "In", "Tot")   # pop stage-class proportions
ltre_pars <- c(vr_pars, pop_pars)

# posterior samples
posterior <- fit$draws( variables = ltre_pars, format = "df" ) %>%
  select( !starts_with(".") )
n_draws <- nrow(posterior)

# pop stage-class proportions available 1986-2021, sW from 1986-2020,
# but s0, sN, sB estimates available 1983-2020 
# hence we remove first three years of s0, sN, sB estimates
posterior <- posterior %>%
  select( - ( contains(str_c("[",1:3,"]")) & contains(c("s0", "sN", "sB")) ) )

# define timeframe for ltre analysis -------------------------------------------

# first and last years with estimates for all variables
# (with understanding that survival probabilities available to last_year - 1)
first_year <- 1986
last_year <- 2021    

# choose first and last years for tltre analysis
ltre_first_year <- 1990
ltre_last_year <- 2021
n_ltre_years <- length(ltre_first_year:ltre_last_year)

# first and last ltre indices
first_t <- which(first_year:last_year == ltre_first_year)
last_t <- which(first_year:last_year == ltre_last_year)

# create list of posterior samples for chosen ltre years by variable -----------

posterior_list <- list()
for (par in ltre_pars){
  par_cols_bool <- str_detect( colnames(posterior), par)
  par_cols <- colnames(posterior)[par_cols_bool]
  if (par %in% vr_pars) {
    par_cols <- par_cols[first_t:(last_t-1)]
  } else {
    par_cols <- par_cols[first_t:last_t]
  }
  posterior_list[[par]] <- posterior %>%
    select( all_of(par_cols) ) %>%
    as.matrix()
}

# calculate population growth rate and immigration rate ------------------------

# calculate omega - available in one less year than other variables
posterior_list$omega <- 
  posterior_list$In[ , 2:n_ltre_years] / posterior_list$Tot[ , 1:(n_ltre_years-1)]

# trim last time of all other variables
posterior_list <- lapply( posterior_list, function(A) { A[, 1:(n_ltre_years - 1)] } )

# expression for realised population growth rate (lambda)
expr_lambda <- expression( (0.5*sW*s0 + sB) * (1 - pNb - pPb) + 
                             sN * (pNb + pPb) + omega )

# realised population growth rate
posterior_list$lambda <- eval( expr_lambda, posterior_list )

# Calculate sensitivities evaluated at temporal means --------------------------
time_means <- lapply( posterior_list, rowMeans )
evaluate_dlambda_at_mean <- function(arg){
  eval( D( expr_lambda, arg ), time_means )
}
lambda_arg_names <- c(vr_pars, "omega", "pPb", "pNb")
sens <- list()
for (par in lambda_arg_names) {
  sens[[par]] <- evaluate_dlambda_at_mean(par)
  # dlambda/domega is identically one, so eval returns length one vector that
  # must be repeated to vector of length n_draws for calculating contributions
  if (par == "omega") { 
    sens[[par]] <- rep( evaluate_dlambda_at_mean(par), n_draws )
  }
}

# Calculate contributions to variance in population growth rate ----------------
contributions <- matrix(NA,n_draws,length(lambda_arg_names))
colnames(contributions) <- lambda_arg_names

for (j in 1:n_draws){
  # compute covariance between jth sample of all lambda arguments
  cov_mat <- posterior_list[lambda_arg_names] %>%
    lapply(function(A){A[j,]}) %>%     # selects jth draw 
    do.call(cbind, .) %>%              # creates matrix with col for each arg
    stats::var()                       # creates covariance matrix
  
  sens_vec <- sens %>%
    lapply(function(x){x[j]}) %>%      # selects jth draw
    do.call(rbind,.)                   # arranges values in column
  
  contributions[j,] <- t( sens_vec * ( cov_mat %*% sens_vec ) ) 
}

summary_contributions <- tibble(par = lambda_arg_names)
for (p in c(.05, .5, .95)){
  name <- str_c("p",p*100)
  summary_contributions[[name]] <- apply(contributions,2,function(x){quantile(x,p)})
}

# save contribution summary
# ipm <- str_split(file, "_fit.RDS")[[1]][1]
# saveRDS(summary_contributions, str_c("outputs/tltre_contributions_", ipm, ".RDS"))

# plot tLTRE contributions
ipm <- str_split(file, "_fit.RDS")[[1]][1]
summary_contributions %>%
  filter( par != "sW") %>%
  ggplot() +
  geom_bar(aes(x = factor(par, levels = unique(par)), y=p50),
           stat = "identity", alpha = 0.7, fill = "navyblue", width = 0.7) +
  geom_linerange(aes(x = par, ymin = p5, ymax = p95), alpha = 0.7) +
  labs(y = "contribution",
       x = "",
       title =
         str_c( "Contribution to variance in realised growth rate ",
                ltre_first_year, " - ", ltre_last_year ),
       subtitle = ipm) +
  theme_classic()

# plot relative difference between variance in realised growth rate
# and sum of tltre contributions

# calculate temporal variance in realised growth rate
temp_var_lambda <- apply( posterior_list$lambda , 1, stats::var )

# calculate sum of contributions
sum_contr <- rowSums( contributions )

# calculate percentage relative difference
rel_diff <- 100 * (temp_var_lambda - sum_contr) / temp_var_lambda

tibble( rel_diff  = rel_diff ) %>%
  ggplot( aes(x = rel_diff) ) +
  geom_density( fill = "navyblue", alpha = 0.6, colour = NA ) +
  theme_classic() +
  geom_vline( xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs( x = "Relative difference (%)" )
# ggsave("figs/si_tltre_linear_approx.pdf", height = 4, width = 6)
