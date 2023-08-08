# Setting priors for the initial population stage sizes in 1986
library(tidyverse)

# set parameter values sourced from literature (see Appendix for sources)
s0 <- 0.72
sN <- 0.83
sB <- 0.83
f3 <- 0.18
f4 <- 0.46
nb <- 0.7*0.69 + 0.3*0.38
bb <- 0.7*0.79 + 0.3*0.52
sW <- 1 - 0.038
omega_vec <- c(.01, .05, .1, .15)

# compute prior distribution parameters ----------------------------------------

# means of population stage sizes 
par_names <- c("Pb1", "Pb2", "Pb3", "Pb4", "Br", "Nb", "In")
prior_means <- tibble()

for ( omega in omega_vec ) {
  # construct the projection matrix
  P <- matrix( c(    0,         0,         0,     0, 0.5*sW*s0,         0, 0.5*sW*s0,
                    sN,         0,         0,     0,         0,         0,         0,
                     0, sN*(1-f3),         0,     0,         0,         0,         0,
                     0,         0, sN*(1-f4),     0,         0,         0,         0,
                     0,     sN*f3,     sN*f4,    sN,     sB*bb,     sN*nb,     sB*bb,
                     0,         0,         0,     0, sB*(1-bb), sN*(1-nb), sB*(1-bb),
                 omega,     omega,     omega, omega,     omega,     omega,     omega) ,
               nrow = 7, ncol = 7, byrow = TRUE)
  
  # calculate the stable stage distribution (ssd)
  eigen <- eigen(P)
  evec <- sapply(eigen$vectors[,1],Re) # leading eigenvalue
  eval <- Re(eigen$values[1])          # corresponding eigenvector
  ssd <- evec / sum(evec)              # normalise to ssd
  names(ssd) <- par_names
  
  # calculate the mean of each prior
  # based on the ssd scaled using 15Oct breeding count of 641 females in 1986
  Ntot <- 641 / ( ssd["Br"] + ssd["In"] )    # Ntot = total population size
  temp_means <- as.list(ssd*Ntot)
  temp_means <- as_tibble(temp_means) %>%
    pivot_longer( everything() ) %>%
    rename( par_name = name, mean = value ) %>%
    add_column(omega = omega)
  
  prior_means <- rbind(prior_means, temp_means)
}

# std deviations of prior population stage sizes 
# set sd relative to the mean to ensure positive size with high probability
prior_pars <- prior_means %>%
  mutate( sd = mean/abs(qnorm(.001)),  # ensures prob of negative size is 0.001
          .after = mean )

# save prior parameters
# saveRDS( prior_pars, "outputs/prior_stage_sizes_1986.RDS" )

# plot prior distributions of population stage sizes ---------------------------

# draw a large sample from all the priors
nsamples <- 10^4
prior_samples <- tibble()

for ( w in omega_vec ) {
  temp_samples <- tibble(.rows = nsamples)
  for ( par in par_names ) {
    mu <- filter( prior_pars, par_name == par, omega == w)$mean
    sig <- filter( prior_pars, par_name == par, omega == w)$sd    
    sample <- rnorm( nsamples, mu, sig )  
    temp_samples <- temp_samples %>%
      add_column(!!par := sample)       # https://stackoverflow.com/questions/68170333
  }
  temp_samples <- temp_samples %>%
    pivot_longer( everything() ) %>%
    add_column(omega = w)
  
  prior_samples <- rbind(prior_samples, temp_samples)
}

# plot 
prior_samples %>%
  filter(name != "In") %>%
  ggplot() +
  facet_wrap( "name" , scales = "free") +
  stat_density(aes(x = value, colour = as.factor(omega) ),  # using stat density cf. geom_density to get lines in legend
               geom="line", position="identity") +
  labs(#title = "Population stage size priors in 1986.",
       #subtitle = "The prior on In is not the shown here as it is based on the 'random effects' structure",
       colour = "omega") +
  theme_classic() +
  theme(axis.title.x = element_blank())




