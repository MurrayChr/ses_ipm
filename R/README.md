# R code 

R code to run all analyses, including model fitting, inspection of model outputs, goodness-of-fit, and
tLTRE analyses. Files sharing a numerical prefix all relate to a common task.

- `01_fit_hmms.R`: fits several versions of the multi-event model
- `02_fit_pop_hmm.R`: fits IPMs in two steps using the posterior from a multievent model in the prior for a population model
- `03_fit_ipms.R`: fits several versions of the IPM (in a single step)
- `04`: goodness-of-fit for multi-event and integrated population models using posterior predictive checks
  - `04a_sim_rep_cmr_data.R`: simulates replicate data from the posterior predictive distribution of a multi-event model
  - `04b_detection_functions.R`: functions to compute observed and expected number of detections of breeding seals and all seals
  - `04c_sim_rep_count_data.R`: simulates replicate count data from the posterior predictive distribution of an IPM
  - `04d_expected_count_function.R`: function to compute expected breeder counts
- `05_*.R`: prior sensitivity to initial population stage sizes
- `06_inspect_hmms.R` and `07_inspect_ipms.R`: examine MCMC diagnostics and plot estimates for fitted multi-event models resp. IPMs
- `08_*.R`: tLTRE analyses from fitted IPMs 
- `09_*.R`: simulate IPM data under different levels of variation in immigration
- `10_*.R`: fit IPMs to simulated data
- `11_make_figures.R`: makes final figures for publication


