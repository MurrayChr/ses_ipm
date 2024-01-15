# Stan code

Stan code to fit all multi-event models and IPMs. 

File prefixes indicate the type of model:

- `hmm`: a multi-event model
- `ipm`: an integrated population model (fit in a single step)
- `pop`: an integrated population model fit in two steps, using using the posterior from a multi-event model in the prior for a population model

The exact models being fit in each file are:

- `hmm_01_fixed_effects_parallel.stan`: the multi-event model with iid uniform priors on all vital rates
- `hmm_02_random_effects_parallel.stan`: the multi-event model with a year random effect on each vital rate
- `hmm_03_gaussian_process_parallel.stan`: the multi-event model with Gaussian Process prior structures on each vital rate
- `hmm_04_cheap_gp_parallel.stan`: the multi-event model with 'cheap' Gaussian Process priors on each vital rate, in which the values of the GP hyperparameters are fixed (rather than estimated).
   
   The 'parallel' suffix in these names refers to the fact that these models have been coded using Stan's `reduce_sum` function to enable within-chain parallelisation. 
   
- `ipm_01.stan`: *IPM<sub>GP,RE</sub>* with GP priors on vital rates, RE prior on immigration
- `ipm_02.stan`: *IPM<sub>RE,RE</sub>* with RE priors on vital rates, RE prior on immigration
- `ipm_03.stan`: *IPM<sub>GP,GP</sub>* with GP priors on vital rates, GP prior on immigration
- `pop_01_logit_mvn.stan`: fits the population model with RE prior on immigration to the count data, using the posterior from a multi-event model to create an informative logit-multivariate normal prior on the vital rates
- `pop_02_gp_in.stan`: as for `pop_01_logit_mvn.stan`, but with a GP prior on immigration
- `pop_03a_sig_c_inform.stan`: as for `pop_01_logit_mvn.stan`, but with the count error standard deviation estimated under an informative prior
- `pop_03a_sig_c_vague.stan`: as for `pop_01_logit_mvn.stan`, but with the count error standard deviation estimated under a vague prior

