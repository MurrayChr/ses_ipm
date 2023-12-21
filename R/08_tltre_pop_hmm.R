#' tLTRE analysis from ipm estimates fit with the iterative two-stage procedure
#' tLTRE contributions computed using get_tltre_contribution() function
library(cmdstanr)
library(posterior)
library(tidyverse)
source("R/08_get_tltre_contribution_function.R")

# Specify the file of an ipm fit with the two-stage procedure ------------------
# these are always of the form "pop_xx_hmm_xx_fit.RDS"

# file <- "pop_01_hmm_04_fit.RDS"   # IPM_GP_RE
# file <- "pop_01_hmm_02_fit.RDS"   # IPM_RE_RE
# file <- "pop_02_hmm_04_fit.RDS"   # IPM_GP_GP
# file <- "pop_03a_hmm_04_fit.RDS"   # IPM_GP_RE, sigma_c ~ informative prior
file <- "pop_03b_hmm_04_fit.RDS"   # IPM_GP_RE, sigma_c ~ vague prior


# extract posterior samples from model -----------------------------------------

# read in fitted model
fit <-readRDS( str_c( "outputs/", file ) )

# parameters needed for tltre
vr_pars <- c("s0", "sN", "sB", "sW")       # survival probabilities 
pop_pars <- c("pNb", "pPb", "In", "Tot")   # pop stage-class proportions
ltre_pars <- c(vr_pars, pop_pars)

# posterior samples
posterior <- fit$draws( variables = ltre_pars, format = "df" ) %>%
  select( !starts_with(".") )
n_draws <- nrow(posterior)

# unlike for ipms fit with the 'full' (not iterative procedure), s0, sN, sB are available 
# from 1986 as they were cropped before fitting the pop model so there is no need to crop them here

# create list of draws x time matrices for each variable
posterior_list <- list()
for (par in ltre_pars){
  par_cols_bool <- str_detect( colnames(posterior), par)
  par_cols <- colnames(posterior)[par_cols_bool]
  posterior_list[[par]] <- posterior %>%
    select( all_of(par_cols) ) %>%
    as.matrix()
}

# compute the tltre contribution for each posterior sample
posterior_contr <- list()
for ( i in 1:n_draws ) {
  posterior_contr[[i]] <-
    get_tltre_contribution(s0 = posterior_list$s0[i,], 
                           sN = posterior_list$sN[i,], 
                           sB = posterior_list$sB[i,], 
                           sW = posterior_list$sW[i,],
                           In = posterior_list$In[i,], 
                           Tot = posterior_list$Tot[i,], 
                           pNb = posterior_list$pNb[i,], 
                           pPb = posterior_list$pPb[i,])$contribution
}

# as a matrix
posterior_contr_mat <- posterior_contr %>%
  do.call(rbind,.)

summary_contributions <- tibble(par = colnames(posterior_contr_mat))
for (p in c(.05, .5, .95)){
  name <- str_c("p",p*100)
  summary_contributions[[name]] <- apply(posterior_contr_mat,2,function(x){quantile(x,p)})
}

# save
# pop_hmm <- str_split(file, "_fit.RDS")[[1]][1]
# out_file <- str_c("outputs/tltre_contributions_", pop_hmm, ".RDS")
# saveRDS( summary_contributions, out_file )
