# Create Figures 2 - 5 in the paper
library(posterior)
library(tidyverse)
library(bayesplot)
library(RColorBrewer)

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
         axis.title = element_text(size = 13)) +
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
library(cowplot)
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
library(cowplot)
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
  library(cowplot)
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
library(cowplot)
plot_grid(sim_truth_vs_est, sim_posteriors,
          align = "h", axis = "tb",
          labels = c("(a)", "(b)"),
          label_size = 13,
          scale = c(1,1/1.14),
          vjust = 13, hjust = -2)
# ggsave("figs/fig5_sim_estimates.pdf", scale = 2)
