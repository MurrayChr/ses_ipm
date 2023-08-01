# Functions to compute observed and expected number of detections 
library(tidyverse)

# get_expected_detections function ---------------------------------------------
# input is replicate data with all structural information necessary
get_expected_detections <- function(rep_data, type){  # type can be "all" or "breeder"
  if (!(type %in% c("all","breeder"))) {
    return("Please enter valid type argument, options: 'all' or 'breeder'")
  }
  # create tibble of cohort information ----------------------------------------
  # numbers in cohort and initial stage distribution (proportion with 1,2 tags)
  cohorts <- rep_data$data_structure$indiv_inits %>%
    group_by(fc, tagloc) %>%
    summarise(count = n()) %>%
    ungroup()
  
  cohorts$dbl_tag <- rep_data$data_structure$indiv_inits %>%
    group_by(fc,init_no_tags,tagloc) %>%
    summarise(count = n()) %>%
    filter(init_no_tags == 2) %>%
    pull(count)
  
  cohorts <- cohorts %>%
    mutate(prop_dbl_tag = dbl_tag/count)
  
  n_cohorts <- dim(cohorts)[1]
  
  # create list of initial distribution vectors 'delta'
  delta <- list()
  for (t in 1:n_cohorts){
    delta[[t]] <- matrix(rep(0,15),nrow=1)
    delta[[t]][1,2] <- cohorts$prop_dbl_tag[t]
    delta[[t]][1,1] <- 1 - cohorts$prop_dbl_tag[t]
  }
  
  # decant transition and emission matrices
  trans_ti <- rep_data$trans_ti
  trans_to <- rep_data$trans_to
  trans_ti_new <- rep_data$trans_ti_new
  emit <- rep_data$emit
  
  # matrix whose (i,t) entry is expected no. of individuals of cohort i
  # detected at time t
  num_years <- rep_data$data_structure$num_years
  expected_detections_by_cohort <- matrix(0,num_years,num_years)
  for (t in 2:num_years){               # zero detections of marked indiv at time 1
    for (i in 1:min((t-1),n_cohorts)){  # i indexes cohort, last cohort at n_cohorts
      # assign theta according to tag
      if (cohorts$tagloc[i] == 1) {           # '1' encodes inner tags
        theta <- trans_ti
      } else if (cohorts$tagloc[i] == 0) {
        theta <- trans_to
      } else {
        theta <- trans_ti_new
      }
      theta_prod <- theta[[i]]
      s <- i
      while(s < t-1){
        s <- s+1
        theta_prod <- theta_prod %*% theta[[s]]
      }
      # cols of emit with desired states
      det_cols <- switch(type, all = 1:14, breeder = c(1:6,9:14))   
      expected_detections_by_cohort[i,t] <- 
        cohorts$count[i]*rowSums(delta[[i]]%*%theta_prod%*%emit[[t]][, det_cols])
    }
  }
  expected_detections <- apply(expected_detections_by_cohort,2,sum)
  expected_detections
}

# get observed detections ------------------------------------------------------
get_observed_detections <- function(rep_data, type) { # type can be "all" or "breeder"
  if (!(type %in% c("all","breeder"))) {
    return("Please enter valid type argument, options: 'all' or 'breeder'")
  }
  # decant information from rep_data
  num_years <- rep_data$data_structure$num_years
  n_cohorts <- rep_data$data_structure$indiv_inits %>%
    pull(fc) %>%
    unique() %>%
    length()
  fc <- rep_data$data_structure$indiv_inits$fc
  y <- rep_data$y
  # observed_detections_by_cohort[i,t] = num observed detections from cohort i in year t
  observed_detections_by_cohort <- matrix(0,num_years,num_years)
  det_codes <- switch(type, all = 1:14, breeder = c(1:6,9:14)) 
  for (t in 2:num_years){
    for (i in 1:min(t-1,n_cohorts)){
      observed_detections_by_cohort[i,t] <- sum(y[fc==i,t] %in% det_codes)
    }
  }
  observed_detections <- apply(observed_detections_by_cohort,2,sum)
  observed_detections
}