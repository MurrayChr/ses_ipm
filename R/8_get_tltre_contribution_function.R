#' A function that takes values for all tltre variables from 1986 - 2021 (with
#' understanding that survival probs available 1986 - 2020 ) and computes
#' a tltre contribution
library(tidyverse)

get_tltre_contribution <- function(s0,  # underyearling survival 1986-2020
                                   sN,  # non- & pre-breeder survival 1986-2020
                                   sB,  # breeder survival 1986-2020
                                   sW,  # survival to weaning 1986-2020
                                   In,  # new immigrant breeders 1986 2021
                                   Tot, # total pop size 1986-2021
                                   pNb, # nonbr proportion of total pop 1986-2021
                                   pPb, # prebr proportion of total pop 1986-2021
                                   tltre_first_year = 1990, # tltre timeframe
                                   tltre_last_year = 2021) {
  # check that inputs have correct lengths
  stopifnot( length(s0)==length(1986:2020) )
  stopifnot( length(sN)==length(1986:2020) )
  stopifnot( length(sB)==length(1986:2020) )
  stopifnot( length(sW)==length(1986:2020) )
  stopifnot( length(In)==length(1986:2021) )
  stopifnot( length(Tot)==length(1986:2021) )
  stopifnot( length(pNb)==length(1986:2021) )
  stopifnot( length(pPb)==length(1986:2021) )
  
  # arrange values in list
  values_list <- list( s0 = s0, sN = sN, sB = sB, sW = sW, In = In, Tot = Tot, 
                       pNb = pNb, pPb = pPb )
  
  # get indices for chosen tltre years, relative to 1986 - 2021
  first_t <- which( 1986:2021 == tltre_first_year )
  last_t <- which( 1986:2021 == tltre_last_year )
  n_tltre_years <- length( tltre_first_year:tltre_last_year )
  
  # crop values list to tltre years
  for ( par in names(values_list) ) {
    if ( par %in% c("s0", "sN", "sB", "sW") ) { 
      values_list[[par]] <- values_list[[par]][first_t:(last_t - 1)]
    }
    else {
      values_list[[par]] <- values_list[[par]][first_t:last_t]
    }
  }
  
  # calculate immigration rate omega
  # available in one less year than other variables
  values_list$omega <- values_list$In[2:n_tltre_years] / values_list$Tot[1:(n_tltre_years-1)]
  
  # trim last time off all other variables
  values_list <- lapply( values_list, function(x) { x[1:(n_tltre_years - 1)] } )
  
  # calculate realised population growth rate lambda using model-based formula
  # (not Tot[t+1] / Tot[t] )
  expr_lambda <- expression( (0.5*sW*s0 + sB) * (1 - pNb - pPb) + 
                               sN * (pNb + pPb) + omega )
  
  # realised population growth rate
  values_list$lambda <- eval( expr_lambda, values_list )
  
  # calculate sensitivities evaluated at temporal means
  time_means <- lapply( values_list, mean )
  evaluate_dlambda_at_mean <- function(arg){
    eval( D( expr_lambda, arg ), time_means )
  }
  lambda_arg_names <- c("s0", "sN", "sB","sW", "omega", "pPb", "pNb")
  sens <- list()
  for ( par in lambda_arg_names ) {
    sens[[par]] <- evaluate_dlambda_at_mean(par)
  }
  
  # calculate contribution to variance in population growth rate
  contribution <- rep( NA, length(lambda_arg_names) )
  names(contribution) <- lambda_arg_names
  cov_mat <- values_list[lambda_arg_names] %>% 
    do.call(cbind, .) %>%              # creates matrix with col for each arg
    stats::var()                       # creates covariance matrix
  
  sens_vec <- sens[lambda_arg_names] %>%  # indexing lambda_arg_names ensures correct order
    do.call(c,.)                          # creates vector from list 
  
  contribution <- t( sens_vec * ( cov_mat %*% sens_vec ) )  # summation achieved in matrix multiplication
  
  # calculate variance in realised growth rate
  var_lambda <- var( values_list$lambda )
  
  # output
  output <- c( values_list, 
               list( cov_mat = cov_mat, sens_vec = sens_vec, 
               contribution = contribution, var_lambda = var_lambda,
               tltre_first_year = tltre_first_year,
               tltre_last_year = tltre_last_year ) )
  output
}


