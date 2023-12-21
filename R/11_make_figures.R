# Create Figures 2 - 5 in the paper, and most of the figures in the appendices
library(posterior)
library(tidyverse)
library(bayesplot)
library(RColorBrewer)
library(cowplot)

# Figure 2: tLTRE contributions, pop growth rate, annual estimates -------------
blues <- brewer.pal(n=9, name = "Blues")
greens <- brewer.pal(n=9, name = "Greens")
oranges <- brewer.pal(n=9, "Oranges")
pairs <- brewer.pal(n = 12, name = "Paired")

ipm_fit <- readRDS("outputs/ipm_01_fit.RDS")
get_var_name <- function(x){ str_split(x, "\\[")[[1]][1] }
ipm_summary <- ipm_fit$summary(variables = c("Br_tot", "Tot")) %>%
  add_column( year = rep(1986:2021,2), .after = "variable" ) %>%
  mutate( var_name = sapply( variable,  get_var_name), .after = "variable" )

# total population size, total breeders and count
count_data <- as_tibble(readRDS("data/count_data_15oct.RDS"))
colour_labels <- c(expression("Tot"), expression(Br["tot"]))
colour_values <- c("Tot" = oranges[7], "Br_tot" = pairs[7])

pop_plot <- ipm_summary %>%
  filter( year >= 1990 ) %>%
  ggplot( aes(x = year) ) +
  geom_pointrange( aes(y = median, ymin = q5, ymax = q95,
                       colour = var_name)) +
  geom_point( data = filter(count_data, year >= 1990),
              mapping = aes(y = count) ) +
  annotate("text", x = 2011.25, y = 1091, label = "Count", size = 5) +
  annotate("point", x = 2008.25, y = 1091, size = 2) +
  coord_cartesian( ylim = c(300,1500)) +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey94"),
         panel.grid.minor = element_line(colour = "grey96"),
         legend.title = element_blank(),
         legend.text = element_text(size = 15),
         legend.text.align = 0,
         legend.position = c(0.625,0.75),
         axis.text.x = element_text(size = 13),
         axis.text.y =  element_text(size = 12),
         axis.title = element_text(size = 13)
         ) +
  scale_color_manual(values = colour_values,
                     labels = colour_labels) +
  labs( y = "number of seals" )
pop_plot

# population growth rate
tot_posterior <- ipm_fit$draws(variables = "Tot", format = "df") %>%
  select( starts_with("Tot") )
n_years <- ncol(tot_posterior)
lambda_posterior <- tot_posterior[,2:n_years] / tot_posterior[,1:(n_years - 1)]
lambda_summary <- tibble(year = 1986:2020, median = NA, q5 = NA, q95 = NA)
lambda_summary$median <- apply(lambda_posterior,2,function(x){quantile(x,0.5)})
lambda_summary$q5 <- apply(lambda_posterior,2,function(x){quantile(x,0.05)})
lambda_summary$q95 <- apply(lambda_posterior,2,function(x){quantile(x,0.95)})

lambda_plot <- lambda_summary %>%
  filter(year >= 1990) %>%
  ggplot( aes(x = year) ) +
  geom_hline( yintercept = 1, colour = "grey", alpha = 0.4 ) +
  geom_pointrange( aes(y = median, ymin = q5, ymax = q95),
                   colour = oranges[7] ) +
  coord_cartesian( ylim = c(0.85,1.2) ) +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey94"),
         panel.grid.minor = element_line(colour = "grey96"),
         axis.text.x = element_text(size = 13),
         axis.text.y =  element_text(size = 12),
         axis.title = element_text(size = 13)) +
  labs( y = "population growth rate" )
lambda_plot

# plot of all three survival probabilities
colour_values <- c("sB" = blues[8], "sN" = blues[6], "s0" = blues[4])
colour_labels <- c(expression(s["B"]), expression(s["N"]), expression(s["0"]))
survival_plot <- ipm_fit$summary(variables = c("s0", "sN", "sB")) %>%
  add_column( year = rep(1983:2020,3), .after = "variable" ) %>%
  mutate( var_name = sapply( variable,  get_var_name), .after = "variable" ) %>%
  filter( year >= 1990 ) %>%
  ggplot( aes(x = year, colour = var_name) ) +
  geom_pointrange( aes(y = median, ymin =q5, ymax = q95)) +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey94"),
         panel.grid.minor = element_line(colour = "grey96"),
         legend.title = element_blank(),
         legend.text = element_text(size = 15),
         legend.position = c(0.65,0.85),
         axis.text.x = element_text(size = 13),
         axis.text.y =  element_text(size = 12),
         axis.title = element_text(size = 13)) +
  scale_color_manual(values = colour_values,
                     labels = colour_labels) +
  labs( y = "survival probability" ) +
  coord_cartesian(ylim = c(0.4,1))
survival_plot

# immigration
immigration_plot <- ipm_fit$summary("In") %>%
  add_column( year = 1986:2021 ) %>%
  filter( year >= 1990 ) %>%
  ggplot( aes(x = year) ) +
  geom_pointrange( aes(y = median, ymin =q5, ymax = q95),
                   colour = pairs[10] ) +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey94"),
         panel.grid.minor = element_line(colour = "grey96"),
         axis.text.x = element_text(size = 13),
         axis.text.y =  element_text(size = 12),
         axis.title = element_text(size = 13)) +
  coord_cartesian( ylim = c(0,NA)) +
  labs( y = "number of immigrants" )
immigration_plot

# pre-breeder and non-breeder proportion
pop_structure_plot <- ipm_fit$summary( variables = c("pPb", "pNb") ) %>%
  add_column( year = rep(1986:2021,2), .after = "variable" ) %>%
  mutate( var_name = sapply( variable,  get_var_name), .after = "variable" ) %>%
  filter( year >= 1990 ) %>%
  ggplot( aes(x = year, colour = var_name) ) +
  geom_pointrange( aes(y = median, ymin =q5, ymax = q95)) +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey94"),
         panel.grid.minor = element_line(colour = "grey96"),
         legend.title = element_blank(),
         legend.text = element_text(size = 15),
         legend.position = c(0.775,0.46),
         axis.text.x = element_text(size = 13),
         axis.text.y =  element_text(size = 12),
         axis.title = element_text(size = 13)) +
  scale_color_manual(values = c("pPb" = greens[5], "pNb" = greens[7]),
                     labels = c(expression(p["Pb"]), expression(p["Nb"]))) +
  labs( y = "proportion of total population" ) +
  coord_cartesian( ylim = c(0,0.45))
pop_structure_plot

# tLTRE contributions
summary_contr <- readRDS("outputs/tltre_contributions_ipm_01.RDS")
summary_contr <- summary_contr %>%
  filter(par != "sW") 

# plot contributions with variable names as x-axis labels
par_labels <- c(expression(s["0"]), expression(s["N"]), expression(s["B"]),
                expression(omega), expression(p["Pb"]), expression(p["Nb"]))
full_plot <- summary_contr %>%
  ggplot( aes(x = par, fill = par) ) +
  geom_bar( aes(y = p50), stat = "identity", width = 0.8, show.legend = FALSE) +
  scale_x_discrete( limits = summary_contr$par,                   # order columns as in tibble
                    labels = par_labels ) +
  scale_fill_manual(values = c("s0" = blues[4], "sB" = blues[8], "sN" = blues[6],
                               "omega" = pairs[10],"pPb" = greens[5], "pNb" = greens[7] ) ) +
  geom_linerange( aes(ymin = p5, ymax = p95), alpha = 0.9) +
  labs( y = "tLTRE contribution") +
  theme_classic( ) +
  theme( axis.title.x = element_blank(),            # remove x axis label "var"
         axis.text.x = element_text(size = 13),     # increase size of x-axis labels 
         axis.text.y =  element_text(size = 12),
         axis.title.y = element_text(size = 13))     
full_plot

# arrange six plots into grid
plot_grid(lambda_plot, pop_plot, full_plot,  
          immigration_plot, survival_plot, pop_structure_plot,
          nrow = 3, ncol = 2,
          align = "vh", axis = "tb",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
          label_size = 16,
          hjust = -5, vjust = 3)
# ggsave("figs/fig2_results.pdf", height = 14, width = 12)










# Fig 3: Prior sensitivity of immigration -----------------------------------
rm(list=ls())
# grab some colours 
blues <- brewer.pal(n=9, name = "Blues")
purples <- brewer.pal(n=9, name = "Purples")
greens <- brewer.pal(n=9, name = "Greens")

# tLTRE plot
# read in tltre contributions ('tc') for various models
tc_ipm_01 <- readRDS("outputs/tltre_contributions_ipm_01.RDS")
tc_ipm_02 <- readRDS("outputs/tltre_contributions_ipm_02.RDS")
tc_ipm_03 <- readRDS("outputs/tltre_contributions_ipm_03.RDS")

# order in 'models' controls order in contributions df which determines
# order in which the model bars are plotted
models <- list("ipm_02" = tc_ipm_02, "ipm_01" = tc_ipm_01,  "ipm_03" = tc_ipm_03)
n_models <- length(models)

# add model name to each summary
for (m in 1:n_models) {
  model_name <- names(models)[m]
  models[[model_name]] <- models[[model_name]] %>%
    add_column(model= model_name)
}

# combine contributions
contributions <- models[[1]]
for (m in 2:n_models) {
  contributions <- rbind(contributions, models[[m]])
}

# add correctly typeset labels
par_labs <- c(expression(s["0"]), expression(s["N"]), expression(s["B"]),
              expression(omega), expression(p["Pb"]), expression(p["Nb"]))

contributions <- contributions %>%
  filter( par != "sW" ) %>%
  add_column( par_labels = rep(par_labs,3) )

colour_values <- c("ipm_01" = blues[8],
                   "ipm_03" = greens[6],
                   "ipm_02"= purples[8])
colour_labels <- c("ipm_01" = expression(IPM["GP,RE"]),
                   "ipm_03" = expression(IPM["GP,GP"]),
                   "ipm_02"= expression(IPM["RE,RE"]))
tltre_plot <- contributions %>%
  ggplot( aes(x = factor(par, levels = unique(par)),
              fill = factor(model, levels = unique(model)), y = p50)  ) +
  geom_bar( stat = "identity", position = position_dodge(), alpha = 0.7 ) +
  scale_x_discrete( limits = unique(contributions$par),
                    labels = contributions$par_labels ) +
  geom_linerange( aes(ymin = p5, ymax = p95 ),
                  position = position_dodge(width = 0.9),
                  alpha = 0.7) +
  labs( y = "tLTRE contribution") +
  scale_fill_manual( values = colour_values,
                     labels = colour_labels) +
  theme_classic( ) +
  theme( axis.title.x = element_blank(), 
         axis.text.x = element_text(size = 13),    
         axis.text.y =  element_text(size = 12),
         axis.title.y = element_text(size = 13),
         legend.title = element_blank(),
         legend.text = element_text(size = 14),
         legend.position = c(0.8,0.8))
tltre_plot

# annual estimates of immigration (without IPM_RE_RE)
ipm_01 <- readRDS("outputs/ipm_01_fit.RDS")
ipm_03 <- readRDS("outputs/ipm_03_fit.RDS")

ipm_01_est <- ipm_01$summary(variables = "In") %>%
  select(variable, median, q5, q95) %>%
  add_column(model = "ipm_01", year = 1986:2021)
ipm_03_est <- ipm_03$summary(variables = "In") %>%
  select(variable, median, q5, q95) %>%
  add_column(model = "ipm_03", year = 1986:2021)

est <- rbind(ipm_01_est, ipm_03_est)
in_plot <- est %>%
  mutate( xshift = sapply(model, function(x){ifelse(x == "ipm_01", -0.15, 0.15)}) ) %>%
  filter( year >= 1990 ) %>%
  ggplot( ) +
  geom_pointrange( aes(x = year + xshift, y = median, ymin = q5, ymax = q95,
                       colour = model ), alpha = 0.9  ) +
  coord_cartesian( ylim = c(0,NA)) +
  labs( x = "year", y = "number of immigrants") +
  scale_colour_manual(values = c("ipm_01" = blues[8], "ipm_03" = greens[6]),
                      labels = c("ipm_01" = expression(IPM["GP,RE"]),
                                 "ipm_03" = expression(IPM["GP,GP"]))) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey94"),
        panel.grid.minor = element_line(colour = "grey96"),
        axis.text.x = element_text(size = 13),     # increase size of x-axis labels
        axis.text.y =  element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.65,0.875))
in_plot

# combine both plots
plot_grid(tltre_plot, in_plot,
          align = "h", axis = "tb",
          labels = c("(a)", "(b)"),
          label_size = 16,
          hjust = -0.1, vjust = 2)
# ggsave("figs/fig3_prior_sensitivity.pdf", height = 5, width = 12)










# Fig 4: PPC for breeder detections with sB estimates --------------------------
rm(list=ls())
ppc_plot <- function( model_prefix, type, 
                      inset_annotate_x, inset_annotate_y, inset_annotate_size,
                      inset_axis_text_size, inset_axis_title_size, 
                      inset_x, inset_y, inset_width, inset_height ) {
  num_detections_real <- readRDS("outputs/ppc_real_data_detections.RDS")
  num_detections <- readRDS(str_c("outputs/ppc_", model_prefix, "_rep_data_detections.RDS"))
  ft_stats <- readRDS(str_c("outputs/ppc_", model_prefix, "_freeman_tukey_detections.RDS"))
  
  # graphical ppc
  graph_plot <- num_detections %>%
    filter( obs_exp == "observed", type == !!(type) ) %>%
    ggplot( aes(x = year, y = count, group = rep) ) +
    geom_line( size = 0.2, alpha = 0.8, colour = "grey" ) +
    geom_line( data = filter(num_detections_real, type == !!(type)),
               aes(x = year, y = count), colour = "navyblue", size = 1) +
    theme_classic( ) +
    theme( axis.text.x = element_text(size = 13),
           axis.text.y =  element_text(size = 12),
           axis.title = element_text(size = 13))
  if (type == "breeder") {
    graph_plot <- graph_plot +
      labs( y = "number of seals detected in breeding season") +
      coord_cartesian(ylim = c(0, 350))
  }
  if (type == "all") {
    graph_plot <- graph_plot + 
      labs( y = "number of seals detected in the year" )
  }
  # bayesian p-value with freeman-tukey statistic
  pB <- ft_stats %>%             # compute Bayesian p-value
    filter(type == !!(type)) %>%
    transmute(bool = y_rep > y) %>%
    summarise(mean(bool))
  
  # find good plot limits
  plot_lim <- ft_stats %>%
    filter(type == !!(type)) %>%
    summarise( min = min(y_rep, y),
               max = max(y_rep, y),
               diam = max(y_rep, y) -  min(y_rep, y))
  
  pb_plot <- ft_stats %>%
    filter(type == !!(type)) %>%
    ggplot(aes(x=y,y=y_rep)) +
    geom_point( size = 2, shape = 21) + 
    geom_abline( slope = 1, intercept = 0, size = 0.5, alpha = 0.2) +
    coord_fixed( xlim = with( plot_lim, c(0, max + 0.05*diam) ), 
                 ylim = with( plot_lim, c(0, max + 0.05*diam) ) ) +
    annotate("text", x = inset_annotate_x, y = inset_annotate_y, 
             size = inset_annotate_size,
             label = str_c("p[B] == ",pB), parse = TRUE) +
    labs( x = "real data discrepancy", y = "replicate data discrepancy") +
    theme_classic( ) +
    theme( axis.text = element_text(size = inset_axis_text_size),
           axis.title = element_text(size = inset_axis_title_size) )
  
  # combine the two plots 
    # p-value as inset
  ppc_plot <- ggdraw(graph_plot) +
    draw_plot(pb_plot, x = inset_x, y = inset_y, width = inset_width, height = inset_height)
  
  ppc_plot
}

# make plot for breeder detections from hmm_02, with random effect priors
type <- "breeder"
model_prefix <- "hmm_02"
ppc_plot_1 <- ppc_plot(model_prefix, type, 
                       inset_annotate_x = 15, inset_annotate_y = 30, inset_annotate_size = 4,
                       inset_axis_text_size = 8, inset_axis_title_size = 9, 
                       inset_x = 0.62, inset_y = 0.13, inset_width = 0.37, inset_height = 0.37 )
ppc_plot_1

# make plot for breeder detections from hmm_04, with gaussian process priors
type <- "breeder"
model_prefix <- "hmm_04"
ppc_plot_2 <- ppc_plot(model_prefix, type, 
                       inset_annotate_x = 9, inset_annotate_y = 18, inset_annotate_size = 4,
                       inset_axis_text_size = 8, inset_axis_title_size = 9, 
                       inset_x = 0.62, inset_y = 0.13, inset_width = 0.37, inset_height = 0.37 )
ppc_plot_2

# add lower panels with breeder survival estimates
hmm_02_fit <- readRDS("outputs/hmm_02_random_effects_fit.RDS")
hmm_04_fit <- readRDS("outputs/hmm_04_cheap_gp_fit.RDS")

sb_hmm_02_plot <- hmm_02_fit$summary(variables = "sB") %>%
  add_column(year = 1983:2020) %>%
  ggplot( aes(x = year) ) +
  geom_pointrange( aes(y = median, ymin = q5, ymax = q95),
                   colour = "navyblue") +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey94"),
         panel.grid.minor = element_line(colour = "grey96"),
         legend.title = element_blank(),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(size = 13),
         axis.text.y =  element_text(size = 12),
         axis.title = element_text(size = 13)) +
  labs( y = "breeder survival probability" ) +
  coord_cartesian(ylim = c(0.2,1))
sb_hmm_02_plot

sb_hmm_04_plot <- hmm_04_fit$summary(variables = "sB") %>%
  add_column(year = 1983:2020) %>%
  ggplot( aes(x = year) ) +
  geom_pointrange( aes(y = median, ymin = q5, ymax = q95),
                   colour = "navyblue") +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey94"),
         panel.grid.minor = element_line(colour = "grey96"),
         legend.title = element_blank(),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(size = 13),
         axis.text.y =  element_text(size = 12),
         axis.title = element_text(size = 13)) +
  labs( y = "breeder survival probability" ) +
  coord_cartesian(ylim = c(0.2,1))
sb_hmm_04_plot

# combine four plots
br_det_w_sb_plot <- plot_grid(ppc_plot_1, ppc_plot_2,
                              sb_hmm_02_plot, sb_hmm_04_plot,
                              nrow = 2, ncol = 2,
                              align = "h", axis = "tb",
                              labels = c("(a)", "(b)", "(c)", "(d)"),
                              label_size = 16,
                              hjust = -3.5, vjust = 3)
br_det_w_sb_plot 
# ggsave(str_c("figs/fig4_ppc_breeder_detections_with_survival_estimates.pdf"),
#        plot = br_det_w_sb_plot, height = 10, width = 12)










# Fig 5: Simulation results ----------------------------------------------------
rm(list=ls())
# grab some colours
Dark2 <- brewer.pal(n=8, name = "Dark2")
RdYlBu <- brewer.pal(n=11, name = "RdYlBu")
colour_values <- c("low" = RdYlBu[1], "med" = RdYlBu[10], "high" = Dark2[6])

sim_contr <- readRDS("outputs/tltre_summary_sim_data.RDS")
# filter for immigration rate parameter 
sim_contr <- sim_contr %>%
  filter( par == "omega" )

# Plot of simulation posteriors by immigration level, and real data posterior
# real data posterior
tltre_full_ipm_01 <- readRDS("outputs/tltre_contributions_ipm_01.RDS")

sim_posteriors <- sim_contr %>%
  add_column(y_coord = (rep(1:100,3) + c(rep(0,100), rep(240,100), rep(120,100)))/10^5 ) %>%
  ggplot( ) +
  geom_pointrange( aes(y = y_coord, x = p50, xmin = p5, xmax = p95, colour = level), 
                   alpha = 0.5) +
  geom_pointrange(data = filter(tltre_full_ipm_01, par=="omega"),
                  aes(y = 50/10^5, x = p50, xmin = p5, xmax = p95),
                  size = 0.8) +
  geom_pointrange(data = filter(tltre_full_ipm_01, par=="omega"),
                  aes(y = 170/10^5, x = p50, xmin = p5, xmax = p95),
                  size = 0.8) +
  geom_pointrange(data = filter(tltre_full_ipm_01, par=="omega"),
                  aes(y = 290/10^5, x = p50, xmin = p5, xmax = p95),
                  size = 0.8) +
  scale_color_manual(values = colour_values) +
  theme_classic() +
  coord_fixed(ratio = 1) +
  theme(axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.75,0.67),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 11)) +
  labs( x = "estimated tLTRE contribution" )
sim_posteriors

# Plot true and estimated tLTRE contributions
# find good plot limits
plot_lim <- sim_contr %>%
  select( truth, p50 ) %>%
  summarise( min = min(truth, p50),
             max = max(truth, p50),
             diam = max(truth, p50) - min(truth, p50))

# choose colour values and plot
sim_truth_vs_est <- sim_contr %>%
  ggplot( aes(x = truth, y = p50, colour = level) ) +
  geom_point( size = 2, shape = 21  ) +
  geom_abline( slope = 1, intercept = 0, size = 0.5, alpha = 0.2) +
  coord_fixed( xlim = with( plot_lim, c(min - 0.05*diam, max + 0.05*diam) ),
               ylim = with( plot_lim, c(min - 0.05*diam, max + 0.05*diam) ) ) +
  theme_classic( ) +
  theme( panel.grid.major = element_line(colour = "grey95"),
         panel.grid.minor = element_line(colour = "grey97"),
         legend.title = element_blank(),
         legend.position = c(0.775,0.2),
         legend.text = element_text(size = 10),
         axis.text.x = element_text(size = 10),
         axis.text.y =  element_text(size = 10),
         axis.title = element_text(size = 11) ) +
  labs(x = "true tLTRE contribution", y = "estimated tLTRE contribution") +
  scale_color_manual(values = colour_values)
sim_truth_vs_est

# combine plots
plot_grid(sim_truth_vs_est, sim_posteriors,
          align = "h", axis = "tb",
          labels = c("(a)", "(b)"),
          label_size = 13,
          scale = c(1,1/1.14),
          vjust = 13, hjust = -2)
# ggsave("figs/fig5_sim_estimates.pdf", scale = 2)







# SI Fig: Vital rates estimates from three HMMs -------------------------------
rm(list=ls())
hmm_01 <- readRDS("outputs/hmm_01_fixed_effects_fit.RDS")
hmm_02 <- readRDS("outputs/hmm_02_random_effects_fit.RDS")
hmm_03 <- readRDS("outputs/hmm_03_gaussian_process_fit.RDS")

vr <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb")
get_var_name <- function(x){ str_split(x, "\\[")[[1]][1] }
get_summary <- function( hmm_fit, vr, model_name ) {
  vr_summary <- hmm_fit$summary( variables = vr ) 
  vr_summary <- vr_summary %>%
    add_column( .after = "variable", 
                var_name = sapply( vr_summary$variable,  get_var_name) ) %>%
    add_column( .after = "var_name", year = rep(1983:2020, length(vr) ) ) %>%
    select( "variable", "var_name", "year", "median", "q5", "q95" ) %>%
    add_column(model = model_name)
}

hmm_01_summary <- get_summary(hmm_01, vr, "fixed effects")
hmm_02_summary <- get_summary(hmm_02, vr, "random effects")
hmm_03_summary <- get_summary(hmm_03, vr, "gaussian process")
summary <- rbind(hmm_01_summary, hmm_02_summary, hmm_03_summary)

# prepare factor levels and labels for faceting
vr_levels <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb")
vr_labels <- c( expression(s[0]),
                expression(s[N]),
                expression(s[B]),
                expression(f[3]),
                expression(f[4]),
                expression(bb),
                expression(nb) )
mod_levels <- c("fixed effects", "gaussian process", "random effects")
mod_labels <- c(expression(fixed~effects), expression(gaussian~process), expression(random~effects))

# plot
summary %>%
  mutate(
    var_name = factor( var_name, levels = vr_levels, labels = vr_labels ),
    model = factor( model, levels = mod_levels, labels = mod_labels )
    ) %>%
  ggplot( aes(x = year) ) +
  geom_point( aes(y = median) ) +
  geom_linerange( aes(ymin = q5, ymax = q95)) +
  coord_cartesian( ylim = c(0,1) ) +
  facet_grid(  vars(var_name), vars(model), labeller = label_parsed ) +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey90"),
         strip.text.y = element_text(angle = 0),
         strip.text = element_text(size=12),
         axis.title.y = element_blank() ) 

# ggsave("figs/si_vital_rates_compared_fe_gp_re.pdf", height = 1.414*10, width = 10)










# SI Fig: Detection probability estimates from three HMMs ---------------------
rm(list=ls())
hmm_01 <- readRDS("outputs/hmm_01_fixed_effects_fit.RDS")
hmm_02 <- readRDS("outputs/hmm_02_random_effects_fit.RDS")
hmm_03 <- readRDS("outputs/hmm_03_gaussian_process_fit.RDS")

det <- c("pBu", "pBe", "qB", "qN")
get_var_name <- function(x){ str_split(x, "\\[")[[1]][1] }
get_summary <- function( hmm_fit, det, model_name ) {
  vr_summary <- hmm_fit$summary( variables = det ) 
  vr_summary <- vr_summary %>%
    add_column( .after = "variable", 
                var_name = sapply( vr_summary$variable,  get_var_name) ) %>%
    add_column( .after = "var_name", year = rep(1983:2021, length(det) ) ) %>%
    select( "variable", "var_name", "year", "median", "q5", "q95" ) %>%
    add_column(model = model_name)
}

hmm_01_summary <- get_summary(hmm_01, det, "fixed effects")
hmm_02_summary <- get_summary(hmm_02, det, "random effects")
hmm_03_summary <- get_summary(hmm_03, det, "gaussian process")
summary <- rbind(hmm_01_summary, hmm_02_summary, hmm_03_summary)

det_levels <- c("pBu", "pBe", "qB", "qN")
det_labels <- c( expression(p[Bu]), 
                 expression(p[Be]), 
                 expression(q[B]), 
                 expression(q[N]) )
mod_levels <- c("fixed effects", "gaussian process", "random effects")
mod_labels <- c(expression(fixed~effects), expression(gaussian~process), expression(random~effects))

summary %>%
  mutate(
    var_name = factor( var_name, levels = det_levels, labels = det_labels ),
    model = factor( model, levels = mod_levels, labels = mod_labels )
  ) %>%
  ggplot( aes(x = year) ) +
  geom_point( aes(y = median) ) +
  geom_linerange( aes(ymin = q5, ymax = q95)) +
  coord_cartesian( ylim = c(0,1) ) +
  facet_grid(  vars(var_name), vars(model), labeller = label_parsed ) +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey90"),
         strip.text.y = element_text(angle = 0),
         strip.text = element_text(size=12),
         axis.title.y = element_blank() ) 

# ggsave("figs/si_detection_probs_compared_fe_gp_re.pdf", height = 1.414*10*(4/7), width = 10)










# SI Fig: Vital rate comparison cheap vs. full GP -----------------------------
rm(list=ls())
hmm_gp <- readRDS("outputs/hmm_03_gaussian_process_fit.RDS")
hmm_cgp <- readRDS("outputs/hmm_04_cheap_gp_fit.RDS")

vr <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb")
get_var_name <- function(x){ str_split(x, "\\[")[[1]][1] }
get_summary <- function( hmm_fit, vr, model_name ) {
  vr_summary <- hmm_fit$summary( variables = vr ) 
  vr_summary <- vr_summary %>%
    add_column( .after = "variable", 
                var_name = sapply( vr_summary$variable,  get_var_name) ) %>%
    add_column( .after = "var_name", year = rep(1983:2020, length(vr) ) ) %>%
    select( "variable", "var_name", "year", "median", "q5", "q95" ) %>%
    add_column(model = model_name)
}

hmm_gp_summary <- get_summary(hmm_gp, vr, "full Gaussian Process")
hmm_cgp_summary <- get_summary(hmm_cgp, vr, "cheap Gaussian Process")
summary <- rbind(hmm_gp_summary, hmm_cgp_summary)

# prepare factor levels and labels for faceting
vr_levels <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb")
vr_labels <- c( expression(s[0]),
                expression(s[N]),
                expression(s[B]),
                expression(f[3]),
                expression(f[4]),
                expression(bb),
                expression(nb) )
mod_levels <- c("full Gaussian Process", "cheap Gaussian Process")
mod_labels <- c(expression(full~Gaussian~Process), 
                expression(cheap~Gaussian~Process))

summary %>%
  mutate(
    var_name = factor( var_name, levels = vr_levels, labels = vr_labels ),
    model = factor( model, levels = mod_levels, labels = mod_labels )
  ) %>%
  ggplot( aes(x = year) ) +
  geom_point( aes(y = median) ) +
  geom_linerange( aes(ymin = q5, ymax = q95)) +
  facet_grid(  vars(var_name), vars(model), labeller = label_parsed,
               scales = "free" ) +
  theme_classic() +
  theme( panel.grid.major = element_line(colour = "grey90"),
         strip.text.y = element_text(angle = 0),
         strip.text = element_text(size=12),
         axis.title.y = element_blank() ) 

# ggsave("figs/si_vital_rates_compared_fgp_cgp.pdf", height = 1.414*10, width = 10)










# SI Fig: Prior vs posterior plots for vital rate GP hyperparameters ----------
# There are four GP hyperparameters, a marginal standard deviation ('sigma') 
# and a lengthscale ('ls') for each of the mean process and the log(sd) process,
# for each of the vital rates!
rm(list=ls())
hmm_gp <- readRDS("outputs/hmm_03_gaussian_process_fit.RDS")
vr <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb")
# hyperparametrs for the mean process
sigma_mean <- str_c("sigma_mean_", vr)
ls_mean <- str_c("ls_mean_", vr)
# hyperparameters for the log(sd) process
sigma_sd <- str_c("sigma_sd_", vr)
ls_sd <- str_c("ls_sd_", vr)

# Prior vs posterior for the mean process
mean_proc_post <- hmm_gp$draws( variables = c(sigma_mean, ls_mean),
                            format = "df" ) %>%
  select( -starts_with(".") ) %>%
  add_column( draws = "posterior" )

mean_proc_prior <- list( draws = rep("prior", 4000) )
halfnormal_prior <- abs( rnorm( 4000, 0, 2) )
mean_proc_lognormal_prior <- rlnorm( 4000, log(0.3), 0.2)
for (p in sigma_mean) {
  mean_proc_prior[[p]] <- halfnormal_prior
}
for (p in ls_mean) {
  mean_proc_prior[[p]] <- mean_proc_lognormal_prior
}
mean_proc_prior <- as_tibble( mean_proc_prior )
mean_process_draws <- rbind( mean_proc_prior, mean_proc_post )

# some helper functions
get_vr_from_col <- function(col) {
  vr <- rep(NA, length(col))
  for (i in 1:length(col)) {
    vr[i] <- str_split(col[i], "_")[[1]][3]
  }
  vr
}
get_gp_par_from_col <- function(col) {
  gppar <- rep(NA, length(col))
  for (i in 1:length(col)) {
    gppar[i] <- str_split(col[i], "_")[[1]][1]
  }
  gppar
}

# prepare factor levels and labels for faceting
vr_levels <- c("s0", "sN", "sB", "f3", "f4", "bb", "nb")
vr_labels <- c( expression(s[0]),
                expression(s[N]),
                expression(s[B]),
                expression(f[3]),
                expression(f[4]),
                expression(bb),
                expression(nb) )
gpp_levels <- c("ls", "sigma")
gpp_labels <- c(expression(eta), expression(sigma))

# plot prior vs. posterior for mean process hyperparameters
mean_plot <- mean_process_draws %>%
  pivot_longer(-"draws") %>%
  mutate(
    draws = factor( draws, levels = c("prior", "posterior") ),
    vr = factor( get_vr_from_col(name), levels = vr_levels, labels = vr_labels ),
    gpp = factor( get_gp_par_from_col(name), levels = gpp_levels, labels = gpp_labels)
  ) %>%
  ggplot( aes(x = value, fill = draws ) ) +
  geom_density(
    colour = NA,
    alpha = 0.7,
    bounds = c(0,Inf)
  ) +
  scale_fill_manual(values = c("grey75", "navyblue")) +
  facet_grid( vars(vr), vars(gpp), scales = "free",
              labeller = label_parsed) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = c(0.85,0.35),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        strip.text.y = element_text(angle = 0))

# Priors vs posteriors for the sd process
sd_proc_post <- hmm_gp$draws( variables = c(sigma_sd, ls_sd),
                                format = "df" ) %>%
  select( -starts_with(".") ) %>%
  add_column( draws = "posterior" )

sd_proc_prior <- list( draws = rep("prior", 4000) )
sd_proc_lognormal_prior <- rlnorm( 4000, log(0.5), 0.2)
for (p in sigma_sd) {
  sd_proc_prior[[p]] <- halfnormal_prior
}
for (p in ls_sd) {
  sd_proc_prior[[p]] <- sd_proc_lognormal_prior
}
sd_proc_prior <- as_tibble( sd_proc_prior )
sd_process_draws <- rbind( sd_proc_prior, sd_proc_post )

# plot prior vs. posterior for sd process hyperparameters
sd_plot <- sd_process_draws %>%
  pivot_longer(-"draws") %>%
  mutate(
    draws = factor( draws, levels = c("prior", "posterior") ),
    vr = factor( get_vr_from_col(name), levels = vr_levels, labels = vr_labels ),
    gpp = factor( get_gp_par_from_col(name), levels = gpp_levels, labels = gpp_labels )
  ) %>%
  ggplot( aes(x = value, fill = draws ) ) +
  geom_density(
    colour = NA,
    alpha = 0.7,
    bounds = c(0,Inf)
  ) +
  scale_fill_manual(values = c("grey75", "navyblue")) +
  facet_grid( vars(vr), vars(gpp), scales = "free",
              labeller = label_parsed) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = c(0.85,0.45),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle = 0))

# combine plots
plot_grid( mean_plot, sd_plot,
          align = "h", axis = "tb",
          labels = c("(a)", "(b)"),
          label_size = 16,
          hjust = -0.15, vjust = 1.75 )

# ggsave("figs/si_vital_rates_gp_hyperpas_prior_v_posterior.pdf", height = 1.2*10, width = 10)









# SI Fig: Prior vs posterior plots for immigration RE hyperparameters --------
# Priors on immigration parameters are (see "ipm_01.stan" lines 488, 489)
# mean_log_lambda ~ normal(log(50), 0.3);  
# sd_log_lambda ~ normal(0, 0.4); (implies half-normal with the >0 variable constraint) 
rm(list=ls())
ipm_01 <- readRDS("outputs/ipm_01_fit.RDS")
post_draws <- ipm_01$draws(variables = c("mean_log_lambda", "sd_log_lambda"),
                           format = "df") %>%
  select(-starts_with(".")) %>%
  add_column(draws = "posterior")

prior_draws <- tibble(
  mean_log_lambda = rnorm(4000, log(50), 0.3),
  sd_log_lambda = abs( rnorm(4000, 0, 0.4) ),
  draws = "prior"
)

rbind(prior_draws, post_draws) %>%
  pivot_longer(cols = c("mean_log_lambda", "sd_log_lambda")) %>%
  mutate( name = factor(name, labels = c(expression(mu[log(Lambda)]), expression(sigma[log(Lambda)]))) ) %>%
  ggplot( aes(x = value, fill = factor(draws, levels = c("prior", "posterior"))) ) +
  geom_density( colour = NA,
                alpha = 0.7,
                bounds = c(0,Inf)) +
  facet_wrap( facets = vars(name), scales = "free", 
              labeller = label_parsed) +
  scale_fill_manual(values = c("grey75", "navyblue")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = c(0.85,0.8),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12))

# ggsave("figs/si_ipm_01_in_re_prior_hyperpars.pdf", width = 8, height = 4)








# SI Fig:  Prior vs posterior plots for immigration GP hyperparameters --------
# in IPM_GP_GP
# Priors on immigration parameters are (see "ipm_03.stan" lines 482-484)
# mean_log_In ~ normal( log(50), .4);   
# ls_log_In ~ lognormal(log(.2), .3);
# sigma_log_In ~ normal(0,2);  (implies half-normal with the >0 variable constraint) 
rm(list=ls())
ipm_03 <- readRDS("outputs/ipm_03_fit.RDS")
post_draws <- ipm_03$draws(variables = c("mean_log_In", "ls_log_In", "sigma_log_In"),
                           format = "df") %>%
  select(-starts_with(".")) %>%
  add_column(draws = "posterior")

prior_draws <- tibble(
  mean_log_In = rnorm( 4000, log(50), 0.4),
  ls_log_In = rlnorm( 4000, log(0.2), 0.3 ),
  sigma_log_In = abs( rnorm(4000, 0, 2) ),
  draws = "prior"
)

# prepare factor levels and labels for faceting
name_levels <- c("mean_log_In", "ls_log_In", "sigma_log_In")
name_labels <- c( expression(mu[log(In)]), expression(eta[log(In)]),
                  expression(sigma[log(In)]) )
rbind(prior_draws, post_draws) %>%
  pivot_longer(cols = c("mean_log_In", "ls_log_In", "sigma_log_In")) %>%
  mutate( name = factor(name, levels = name_levels, labels = name_labels) ) %>%
  ggplot( aes(x = value, fill = factor(draws, levels = c("prior", "posterior"))) ) +
  geom_density( colour = NA,
                alpha = 0.7,
                bounds = c(0,Inf)) +
  facet_wrap( facets = vars(name), scales = "free", 
              labeller = label_parsed) +
  scale_fill_manual(values = c("grey75", "navyblue")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = c(0.9,0.5),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12))

# ggsave("figs/si_ipm_03_in_gp_prior_hyperpars.pdf", width = 10, height = 4)










# SI Fig: Posterior predictive checks for multievent models -------------------
rm(list=ls())

ppc_plot <- function( model_prefix, type, 
                      inset_annotate_x, inset_annotate_y, inset_annotate_size,
                      inset_axis_text_size, inset_axis_title_size, 
                      inset_x, inset_y, inset_width, inset_height ) {
  num_detections_real <- readRDS("outputs/ppc_real_data_detections.RDS")
  num_detections <- readRDS(str_c("outputs/ppc_", model_prefix, "_rep_data_detections.RDS"))
  ft_stats <- readRDS(str_c("outputs/ppc_", model_prefix, "_freeman_tukey_detections.RDS"))
  
  # graphical ppc
  graph_plot <- num_detections %>%
    filter( obs_exp == "observed", type == !!(type) ) %>%
    ggplot( aes(x = year, y = count, group = rep) ) +
    geom_line( size = 0.2, alpha = 0.8, colour = "grey" ) +
    geom_line( data = filter(num_detections_real, type == !!(type)),
               aes(x = year, y = count), colour = "navyblue", size = 1) +
    theme_classic( ) +
    theme( axis.text.x = element_text(size = 13),
           axis.text.y =  element_text(size = 12),
           axis.title = element_text(size = 13))
  if (type == "breeder") {
    graph_plot <- graph_plot +
      labs( y = "number of seals detected in breeding season") +
      coord_cartesian(ylim = c(0, 350))
  }
  if (type == "all") {
    graph_plot <- graph_plot + 
      labs( y = "number of seals detected in the year" )
  }
  # bayesian p-value with freeman-tukey statistic
  pB <- ft_stats %>%             # compute Bayesian p-value
    filter(type == !!(type)) %>%
    transmute(bool = y_rep > y) %>%
    summarise(mean(bool))
  
  # find good plot limits
  plot_lim <- ft_stats %>%
    filter(type == !!(type)) %>%
    summarise( min = min(y_rep, y),
               max = max(y_rep, y),
               diam = max(y_rep, y) -  min(y_rep, y))
  
  pb_plot <- ft_stats %>%
    filter(type == !!(type)) %>%
    ggplot(aes(x=y,y=y_rep)) +
    geom_point( size = 2, shape = 21) + 
    geom_abline( slope = 1, intercept = 0, linewidth = 0.5, alpha = 0.2) +
    coord_fixed( xlim = with( plot_lim, c(0, max + 0.05*diam) ), 
                 ylim = with( plot_lim, c(0, max + 0.05*diam) ) ) +
    annotate("text", x = inset_annotate_x, y = inset_annotate_y, 
             size = inset_annotate_size,
             label = str_c("p[B] == ",pB), parse = TRUE) +
    labs( x = "real data discrepancy", y = "replicate data discrepancy") +
    theme_classic( ) +
    theme( axis.text = element_text(size = inset_axis_text_size),
           axis.title = element_text(size = inset_axis_title_size) )
  
  # combine the two plots 
    # p-value as inset
  ppc_plot <- ggdraw(graph_plot) +
    draw_plot(pb_plot, x = inset_x, y = inset_y, width = inset_width, height = inset_height)
  
  ppc_plot
}

# make plot for breeder detections from hmm_04, with gaussian process priors
type <- "breeder"
model_prefix <- "hmm_04"
ppc_hmm_04_breeder <- ppc_plot(model_prefix, type, 
                               inset_annotate_x = 9, inset_annotate_y = 18, inset_annotate_size = 4,
                               inset_axis_text_size = 8, inset_axis_title_size = 9, 
                               inset_x = 0.62, inset_y = 0.13, inset_width = 0.37, inset_height = 0.37 )
# ggsave(str_c("figs/ppc_", model_prefix, "_", type,"_detections.pdf"),
#        plot = ppc_hmm_04_breeder, height = 5, width = 7.5)

# make plot for all detections from hmm_04, with gaussian process priors
type <- "all"
model_prefix <- "hmm_04"
ppc_hmm_04_all <- ppc_plot(model_prefix, type, 
                           inset_annotate_x = 9, inset_annotate_y = 18, inset_annotate_size = 4,
                           inset_axis_text_size = 8, inset_axis_title_size = 9, 
                           inset_x = 0.62, inset_y = 0.13, inset_width = 0.37, inset_height = 0.37 )
# ggsave(str_c("figs/ppc_", model_prefix, "_", type,"_detections.pdf"),
#        plot = ppc_hmm_04_all, height = 5, width = 7.5)

# combine plots of breeder detections and all detections for hmm_04
hmm_04_ppc <- plot_grid(ppc_hmm_04_breeder, ppc_hmm_04_all,
                        align = "h", axis = "tb",
                        labels = c("(a)", "(b)"),
                        label_size = 16,
                        hjust = -3.5, vjust = 3)
hmm_04_ppc
# ggsave(str_c("figs/ppc_hmm_04_breeder_and_all_detections.pdf"),
#        plot = hmm_04_ppc, height = 5, width = 12)









# SI Fig: PPCs for population model based on the count data -------------------
rm(list=ls())

# read in replicate and real breeder counts
counts <- readRDS("outputs/ppc_ipm_01_rep_counts.RDS")
real_count <- readRDS("data/count_data_15oct.RDS")

# graphical ppc
graph_plot <- counts %>%
  ggplot( aes(x = year, y = count) ) +
  geom_line( aes(group = rep_num), size = 0.2, colour = "grey" ) +
  geom_line( data = real_count, colour = "navyblue", size = 1) +
  labs( y = "breeding count") +
  theme_classic( ) +
  theme( axis.text.x = element_text(size = 13), 
         axis.text.y =  element_text(size = 12),
         axis.title = element_text(size = 13) )
graph_plot

# bayesian p-value for freeman-tukey statistic
ft_stats <- readRDS("outputs/ppc_ipm_01_freeman_tukey_breeder_counts.RDS")

pB <- ft_stats %>%             # compute Bayesian p-value
  transmute(bool = y_rep > y) %>%
  summarise(mean(bool))

# find good plot limits
plot_lim <- ft_stats %>%
  summarise( min = min(y_rep, y),
             max = max(y_rep, y),
             diam = max(y_rep, y) -  min(y_rep, y))

pb_plot <- ft_stats %>%
  ggplot(aes(x=y,y=y_rep)) +
  geom_point( size = 2, shape = 21) + 
  geom_abline( slope = 1, intercept = 0, size = 0.5, alpha = 0.2) +
  coord_fixed( xlim = with( plot_lim, c(0, max + 0.05*diam) ), 
               ylim = with( plot_lim, c(0, max + 0.05*diam) ) ) +
  annotate("text", x = 45, y = 25, size = 4,
           label = str_c("p[B] == ",pB), parse = TRUE) +
  labs( x = "real data discrepancy", y = "replicate data discrepancy") +
  theme_classic( ) +
  theme( axis.text = element_text(size = 8),
         axis.title = element_text(size = 9) )

# combine plots with p-value as inset
ppc_plot <- ggdraw(graph_plot) +
  draw_plot(pb_plot, x = 0.24, y = 0.53, width = 0.37, height = 0.37)
# ppc_plot
# ggsave("figs/ppc_ipm_01_breeder_counts.pdf", ppc_plot, height = 5, width = 6)









# SI Fig: Estimating count error standard deviation sigma_c -------------------
# read in models
rm(list=ls())
ipm_01 <- readRDS("outputs/ipm_01_fit.RDS")  # sigma_c = 10
pop_03a_hmm_04 <- readRDS("outputs/pop_03a_hmm_04_fit.RDS")  # sigma_c ~ normal(10,2.5)
pop_03b_hmm_04 <- readRDS("outputs/pop_03b_hmm_04_fit.RDS")  # sigma_c ~ halfnormal(0,20)

# prior vs posterior plots for sigma_c

# informative prior 
inform_post <- pop_03a_hmm_04$draws(variables = "sig_c", format = "df") %>%
  select(-starts_with(".")) %>%
  add_column( draws = "posterior" )
inform_prior <- tibble( sig_c = rnorm(4000,10,2.5), draws = "prior" )
inform_plot <- rbind( inform_post, inform_prior ) %>%
  ggplot( aes(x = sig_c, fill = factor(draws, levels = c("prior", "posterior")) ) ) +
  geom_density( colour = NA, alpha = 0.7 ) +
  scale_fill_manual(values = c(posterior = "navyblue", prior = "grey75")) +
  theme_classic() +
  theme(
    legend.title = element_blank(), 
    legend.position = c(0.8,0.7),
    legend.text = element_text(size=12),
    axis.title.x = element_text(size=13)
  ) +
  labs( x = expression(sigma[c]) )
inform_plot

# vague prior
vague_post <- pop_03b_hmm_04$draws(variables = "sig_c", format = "df") %>%
  select(-starts_with(".")) %>%
  add_column( draws = "posterior" )
vague_prior <- tibble( sig_c = abs( rnorm(4000,0,20) ), draws = "prior" )
vague_plot <- rbind( vague_post, vague_prior ) %>%
  ggplot( aes(x = sig_c, fill = factor(draws, levels = c("prior", "posterior")) ) ) +
  geom_density( colour = NA, alpha = 0.7, bounds = c(0,Inf) ) +
  scale_fill_manual(values = c(posterior = "navyblue", prior = "grey75")) +
  theme_classic() +
  theme(
    legend.title = element_blank(), 
    legend.position = c(0.7,0.7),
    legend.text = element_text(size=12),
    axis.title.x = element_text(size=13)
  ) +
  labs( x = expression(sigma[c]) )
vague_plot

# combine both plots
plot_grid(inform_plot, vague_plot,
          align = "h", axis = "tb",
          labels = c("(a)", "(b)"),
          label_size = 16,
          hjust = -0.1, vjust = 2)
# ggsave("figs/si_fig_sigma_c_prior_v_posterior.pdf", height = 5, width = 12)


# comparison of immigration annual estimates
ipm_01_in <- ipm_01$summary(variables = "In") %>%
  add_column( model = "ipm_01", year = 1986:2021)
pop_03a_in <- pop_03a_hmm_04$summary(variables = "In") %>%
  add_column( model = "pop_03a", year = 1986:2021)
pop_03b_in <- pop_03b_hmm_04$summary(variables = "In") %>%
  add_column( model = "pop_03b", year = 1986:2021)
est <- rbind( ipm_01_in, pop_03a_in, pop_03b_in )
xshift <- function(x) {switch(x, ipm_01 = -1, pop_03a = 0, pop_03b = 1)}
blues <- brewer.pal(n=9, name = "Blues")
in_plot <- est %>%
  mutate( xshift = 0.2*sapply(model, xshift) ) %>%
  filter( year >= 1990 ) %>%
  ggplot( ) +
  geom_pointrange( aes(x = year + xshift, y = median, ymin = q5, ymax = q95,
                       colour = model ), alpha = 0.9  ) +
  coord_cartesian( ylim = c(0,NA)) +
  labs( x = "year", y = "number of immigrants") +
  scale_colour_manual(values = c("ipm_01" = blues[8], "pop_03a" = blues[6], 
                                 "pop_03b" = blues[4]),
                      labels = c("ipm_01" = "value fixed",
                                 "pop_03a" = "informative prior",
                                 "pop_03b" = "vague prior")) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey94"),
        panel.grid.minor = element_line(colour = "grey96"),
        axis.text.x = element_text(size = 13),     # increase size of x-axis labels
        axis.text.y =  element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = c(0.65,0.875))
in_plot

# tLTRE plot
# read in tltre contributions ('tc') for various models
tc_ipm_01 <- readRDS("outputs/tltre_contributions_ipm_01.RDS")
tc_pop_03a <- readRDS("outputs/tltre_contributions_pop_03a_hmm_04.RDS")
tc_pop_03b <- readRDS("outputs/tltre_contributions_pop_03b_hmm_04.RDS")

# order in 'models' controls order in contributions df which determines
# order in which the model bars are plotted
models <- list("pop_03b" = tc_pop_03b, "pop_03a" = tc_pop_03a,  "ipm_01" = tc_ipm_01)
n_models <- length(models)

# add model name to each summary
for (m in 1:n_models) {
  model_name <- names(models)[m]
  models[[model_name]] <- models[[model_name]] %>%
    add_column(model= model_name)
}

# combine contributions
contributions <- models[[1]]
for (m in 2:n_models) {
  contributions <- rbind(contributions, models[[m]])
}

# add correctly typeset labels
par_labs <- c(expression(s["0"]), expression(s["N"]), expression(s["B"]),
              expression(omega), expression(p["Pb"]), expression(p["Nb"]))

contributions <- contributions %>%
  filter( par != "sW" ) %>%
  add_column( par_labels = rep(par_labs,3) )

colour_values <- c("pop_03b" = blues[4], 
                   "pop_03a" = blues[6], 
                   "ipm_01" = blues[8] )                   
colour_labels <- c("pop_03b" = "vague prior",
                   "pop_03a" = "informative prior",
                   "ipm_01" = "value fixed")
tltre_plot <- contributions %>%
  ggplot( aes(x = factor(par, levels = unique(par)),
              fill = factor(model, levels = unique(model)), y = p50)  ) +
  geom_bar( stat = "identity", position = position_dodge(), alpha = 0.7 ) +
  scale_x_discrete( limits = unique(contributions$par),
                    labels = contributions$par_labels ) +
  geom_linerange( aes(ymin = p5, ymax = p95 ),
                  position = position_dodge(width = 0.9),
                  alpha = 0.7) +
  labs( y = "tLTRE contribution") +
  scale_fill_manual( values = colour_values,
                     labels = colour_labels) +
  theme_classic( ) +
  theme( axis.title.x = element_blank(), 
         axis.text.x = element_text(size = 13),    
         axis.text.y =  element_text(size = 12),
         axis.title.y = element_text(size = 13),
         legend.title = element_blank(),
         legend.text = element_text(size = 14),
         legend.position = c(0.2,0.8))
tltre_plot

# combine both plots
plot_grid(tltre_plot, in_plot,
          align = "h", axis = "tb",
          labels = c("(a)", "(b)"),
          label_size = 16,
          hjust = -0.1, vjust = 2)
# ggsave("figs/si_fig_in_sensitivity_to_sigma_c.pdf", height = 5, width = 12)









# Posterior tag loss probabilities --------------------------------------------
hmm_04 <- readRDS("outputs/hmm_04_cheap_gp_fit.RDS")

# prepare factor for plotting with facet_wrap
name_levels <- c("ti", "to", "ti_new")
name_labels <- c( expression(tau[inner]),
                  expression(tau[outer]),
                  expression(tau[new~inner]) )
hmm_04$draws(variables = c("ti", "to", "ti_new"), format = "df") %>%
  select(-starts_with(".")) %>%
  pivot_longer( everything() ) %>%
  mutate( name = factor(name, levels=name_levels, labels=name_labels ) ) %>%
  ggplot( aes(x = value) ) +
  geom_density( colour = NA, alpha = 0.7, fill = "navyblue" ) +
  facet_wrap( vars(name), labeller = label_parsed ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(),
        legend.text = element_text(size = 12),
        legend.position = c(0.9,0.5),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12))
# ggsave("figs/si_tag_loss_posterior.pdf", height = 2, width = 5.5)
