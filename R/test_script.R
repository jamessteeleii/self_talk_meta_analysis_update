library(metafor)
library(orchaRd)
library(tidyverse)
library(brms)
library(modelr)
library(tidybayes)
library(bayesplot)
library(rstan)
library(ggridges)
library(patchwork)

data <- read.csv("data/Final data.csv")

data$ri <- 0.7

data$pre_sd_pool <-
  sqrt(((data$n_st - 1) * data$pre_sd_st ^ 2 + (data$n_con - 1) * data$pre_sd_con ^
          2
  ) / (data$n_st + data$n_con - 2))

data$post_sd_pool <-
  sqrt(((data$n_st - 1) * data$post_sd_st ^ 2 + (data$n_con - 1) * data$post_sd_con ^
          2
  ) / (data$n_st + data$n_con - 2))

##### Effect sizes

### For effects where increase is good
data_increase <- subset(data, improvement == "Increase")

# Calculate pre-post control effect sizes
data_increase_ppc <-
  subset(data_increase, study_design == "between")

data_increase_ppc_st <- escalc(
  measure = "SMCR",
  m1i = post_m_st,
  m2i = pre_m_st,
  sd1i = pre_sd_pool,
  ni = n_st,
  ri = ri,
  data = data_increase_ppc
)
data_increase_ppc_con <- escalc(
  measure = "SMCR",
  m1i = post_m_con,
  m2i = pre_m_con,
  sd1i = pre_sd_pool,
  ni = n_con,
  ri = ri,
  data = data_increase_ppc
)

data_increase_ppc$yi <-
  (data_increase_ppc_st$yi - data_increase_ppc_con$yi)
data_increase_ppc$vi <-
  (data_increase_ppc_st$vi + data_increase_ppc_con$vi)

# Calculate post only effect sizes
data_increase_post <-
  subset(data_increase, study_design == "between-post")

data_increase_post$post_sd_st <-
  replmiss(data_increase_post$post_sd_st,
           with(data_increase_post, post_sd_con))

data_increase_post <- escalc(
  measure = "SMD",
  m1i = post_m_st,
  m2i = post_m_con,
  sd1i = post_sd_st,
  sd2i = post_sd_con,
  n1i = n_st,
  n2i = n_con,
  data = data_increase_post
)

# Calculate pre-post effect sizes
data_increase_pp <-
  subset(data_increase, study_design == "pre-post")

data_increase_pp <- escalc(
  measure = "SMCR",
  m1i = post_m_st,
  m2i = pre_m_st,
  sd1i = pre_sd_st,
  ni = n_st,
  ri = ri,
  data = data_increase_pp
)

### For effects where decrease is good
data_decrease <- subset(data, improvement == "Decrease")

# Calculate pre-post control effect sizes
data_decrease_ppc <-
  subset(data_decrease, study_design == "between")

data_decrease_ppc_st <- escalc(
  measure = "SMCR",
  m1i = pre_m_st,
  m2i = post_m_st,
  sd1i = pre_sd_pool,
  ni = n_st,
  ri = ri,
  data = data_decrease_ppc
)
data_decrease_ppc_con <- escalc(
  measure = "SMCR",
  m1i = pre_m_con,
  m2i = post_m_con,
  sd1i = pre_sd_pool,
  ni = n_con,
  ri = ri,
  data = data_decrease_ppc
)

data_decrease_ppc$yi <-
  (data_decrease_ppc_st$yi - data_decrease_ppc_con$yi)
data_decrease_ppc$vi <-
  (data_decrease_ppc_st$vi + data_decrease_ppc_con$vi)

# Calculate post only effect sizes
data_decrease_post <-
  subset(data_decrease, study_design == "between-post")

data_decrease_post$post_sd_st <-
  replmiss(data_decrease_post$post_sd_st,
           with(data_decrease_post, post_sd_con))

data_decrease_post <- escalc(
  measure = "SMD",
  m1i = post_m_con,
  m2i = post_m_st,
  sd1i = post_sd_st,
  sd2i = post_sd_con,
  n1i = n_con,
  n2i = n_st,
  data = data_decrease_post
)

# Calculate pre-post effect sizes
data_decrease_pp <-
  subset(data_decrease, study_design == "pre-post")

data_decrease_pp <- escalc(
  measure = "SMCR",
  m1i = pre_m_st,
  m2i = post_m_st,
  sd1i = pre_sd_st,
  ni = n_st,
  ri = ri,
  data = data_decrease_pp
)


### Combine all data
data <-
  rbind(
    data_decrease_post,
    data_decrease_pp,
    data_decrease_ppc,
    data_increase_post,
    data_increase_pp,
    data_increase_ppc
  )

      # ##### Fit initial model
      # MultiLevelModel_st <- rma.mv(
      #   yi,
      #   V = vi,
      #   data = data,
      #   slab = paste(label),
      #   random = list( ~ 1 |
      #                    study / group / effect),
      #   method = "REML",
      #   test = "t",
      #   control = list(optimizer = "optim", optmethod =
      #                    "Nelder-Mead")
      # )
      #
      # save(MultiLevelModel_st, file = "models/MultiLevelModel_st")
      #
      # ### Calculate I^2
      # I2_st <- i2_ml(MultiLevelModel_st)
      #
      # ### Calculate robust estimate from multi-level model
      # RobuEstMultiLevelModel_st <- robust(MultiLevelModel_st, data$study)
      #
      # save(RobuEstMultiLevelModel_st, file = "models/RobuEstMultiLevelModel_st")
      #
      # ### Caterpillar plot
      #
      # # Overall estimate
      # diamond_st <-
      #   data.frame(
      #     x = c(
      #       RobuEstMultiLevelModel_st$b[1] + (RobuEstMultiLevelModel_st$se * 1.96),
      #       RobuEstMultiLevelModel_st$b[1],
      #       RobuEstMultiLevelModel_st$b[1] - (RobuEstMultiLevelModel_st$se *
      #                                           1.96),
      #       RobuEstMultiLevelModel_st$b[1]
      #     ),
      #     y = c(-7.5, -12.5, -7.5, -2.5)
      #   )
      #
      # # Prediction interval
      # PI_st <- as.data.frame(predict(RobuEstMultiLevelModel_st))
      #
      # # I^2 labels
      # I2_st_lab <- data.frame(level = c("study", "arm", "es"),
      #                         I2 = I2_st[2:4]) %>%
      #   pivot_wider(names_from = "level", values_from = "I2")
      #
      # # Plot
      # forest_st <- data %>%
      #   mutate(
      #     number = seq(1:nrow(data)),
      #     es = factor(number, levels = number[order(yi)]),
      #     se = sqrt(vi)
      #   ) %>%
      #   ggplot(aes(x = yi, y = es)) +
      #   geom_vline(xintercept = 0, linetype = "dashed") +
      #   geom_point(size = 0.5, alpha = 0.8) +
      #   scale_x_continuous(limits = c(-2, 7.5),
      #                      breaks = c(-1, 0, 1, 2, 3, 4, 5, 6)) +
      #   geom_linerange(aes(xmin = yi - (se * 1.96), xmax = yi + (se * 1.96)), alpha =
      #                    0.8) +
      #   scale_y_discrete(limits = rev, expand = expansion(mult = c(0.075, 0))) +
      #   geom_text(
      #     data = mutate_if(PI_st,
      #                      is.numeric, round, 2),
      #     aes(
      #       label = glue::glue(
      #         "Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"
      #       ),
      #       x = 5,
      #       y = 45
      #     ),
      #     hjust = "centre",
      #     size = 3
      #   ) +
      #   geom_text(
      #     data = PI_st,
      #     aes(
      #       label = glue::glue(
      #         "[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"
      #       ),
      #       x = 5,
      #       y = 40
      #     ),
      #     hjust = "centre",
      #     size = 3
      #   ) +
      #   geom_text(
      #     data = I2_st_lab,
      #     aes(
      #       label = glue::glue(
      #         "I^2 [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"
      #       ),
      #       x = 5,
      #       y = 35
      #     ),
      #     hjust = "centre",
      #     size = 3
      #   ) +
      #   geom_segment(
      #     data = PI_st,
      #     aes(
      #       y = -7.5,
      #       yend = -7.5,
      #       x = pi.lb,
      #       xend = pi.ub
      #     ),
      #     arrow = arrow(
      #       length = unit(0.30, "cm"),
      #       angle = 90,
      #       ends = "both",
      #       type = "open"
      #     )
      #   ) +
      #   geom_polygon(data = diamond_st, aes(x = x, y = y)) +
      #   labs(y = "",
      #        x = "Standardised Mean Difference (Positive Values Favour Self-Talk)",
      #        title = "Overall Model") +
      #   theme_classic() +
      #   theme(
      #     axis.text.y = element_blank(),
      #     axis.ticks.y = element_blank(),
      #     axis.line.y = element_blank()
      #   )
      #
      # forest_st

# run rstan quicker - for bayesian analysis
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)


# Main model

get_prior(yi | se(1 / vi) ~ 1 + (1 | study / experiment / group / effect), data)

prior_meta <-
  prior("student_t(3, 0.48, 0.05)", class = "Intercept")

ggplot() +
  stat_function(
    data = data.frame(x = c(-2, 2)),
    aes(x),
    fun = dstudent_t,
    n = 101,
    args = list(
      df = 3,
      mu = 0.48,
      sigma = 0.05,
      log = FALSE
    ),
    color = "black"
  ) +
  geom_point(aes(y=0, x=0.48)) +
  geom_errorbarh(aes(y= 0, xmin=0.38, xmax=0.58)) +
  labs(x = "Prior for model intercept (Standardised Mean Difference)" ,
       y = "PDF",
       title = "Visualisation of prior distribution",
       subtitle = "Taken from main model intercept estimate from Hatzigeorgiadis et al., (2011)") +
  scale_x_continuous(limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA))


main_model <-
  brm(
    yi | se(sqrt(vi)) ~ 1 + (1 | study / experiment / group / effect),
    data = data,
    prior = prior_meta,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    control = list(adapt_delta = 0.99)
  )

# save(main_model, file = "models/main_model")

plot(main_model)

pp_check(main_model)

study_draws <- main_model %>%
  spread_draws(b_Intercept, r_study[study, ]) %>%
  mutate(b_Intercept = b_Intercept + r_study,
         study = as.factor(study)) %>%
  select(.chain, .iteration, .draw, b_Intercept, study) %>%
  merge(data_effect_sizes[, 1:2], by = "study")

study_summary <- group_by(study_draws, label) %>%
  mean_qi(b_Intercept) %>%
  mutate(b_Intercept = b_Intercept,
         .lower = .lower,
         .upper = .upper)

pooled_draws <- main_model %>%
  spread_draws(b_Intercept) %>%
  mutate(study = "Posterior Pooled Estimate",
         label = "Posterior Pooled Estimate")

pooled_summary <- group_by(pooled_draws, label) %>%
  mean_qi(b_Intercept) %>%
  mutate(b_Intercept = b_Intercept,
         .lower = .lower,
         .upper = .upper)

nd <- data.frame(study = "new", vi = 0)

pred_int_data <- posterior_predict(
  object = main_model,
  newdata = nd,
  re_formula = NULL,
  allow_new_levels = TRUE,
  sample_new_levels = "gaussian"
)

pred_int_data <- median_qi(pred_int_data) %>%
  mutate(label = as.character("Posterior Pooled Estimate"))

prior_plot <-
  data.frame(b_Intercept = rstudent_t(
    n = 1e6,
    df = 3,
    mu = 0.48,
    sigma = 0.05
  ),
  label = "Prior (Hatzigeorgiadis et al., 2011)")


forest_study <- ggplot(aes(x = b_Intercept,
                           y = reorder(label, b_Intercept)),
                       data = study_draws) +


  # # Add vertical lines and band for pooled effect and CI
  annotate(
    "rect",
    xmin = pred_int_data$ymin,
    xmax = pred_int_data$ymax,
    ymin = -Inf,
    ymax = Inf,
    alpha = .1,
    fill = "black"
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  # geom_vline(
  #   xintercept = fixef(main_model)[1, 1],
  #   color = "black",
  #   linewidth = 1
  # ) +
  # geom_vline(
  #   xintercept = fixef(main_model)[1, 3:4],
  #   color = "black",
  #   linewidth = 0.25
  # ) +
  scale_y_discrete() +
  # geom_segment(aes(y = y, yend = yend, x = x, xend = xend),
  #              size = 0.5,
  #              data.frame(y = 1, yend = 1, x = pred_int_data$ymin, xend = pred_int_data$ymax),
  #              arrow = arrow(angle = 90, length = unit(0.1, "inches"), ends = "both")) +

  # Add densities
  geom_density_ridges(
    fill = "darkgrey",
    rel_min_height = 0.01,
    col = NA,
    scale = 1,
    alpha = 0.8
  ) +
  stat_pointinterval(point_interval = mean_qi,
                     .width = .95,
                     size = 0.75) +

  # Add individual study data subset(Data_effects, !is.na(yi)), aes(x = yi, y = study_name)
  geom_point(
    data = subset(data_effect_sizes,!is.na(yi)),
    aes(x = yi, y = label),
    position = position_nudge(y = -0.1),
    shape = "|"
  ) +

  # Add text and labels
  geom_text(
    data = mutate_if(study_summary,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"),
      x = 4
    ),
    hjust = "inward",
    size = 3
  ) +
  labs(x = element_blank(), # summary measure
       y = element_blank(),
       title = "Study level estimates",
       subtitle = "Posterior distributions, mean and 95% quantile intervals, individual effects (ticks), and 95% prediction interval (grey band)") +
  scale_x_continuous(limits = c(-2, 4), breaks = c(-2, -1, 0, 1, 2, 3, 4)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        plot.subtitle = element_text(size=6))

prior_posterior <- rbind(pooled_draws[, c(4, 6)], prior_plot) %>%
  mutate(label = factor(label, levels = c("Prior (Hatzigeorgiadis et al., 2011)",
                                          "Posterior Pooled Estimate")))

posterior_update <- ggplot(data = prior_posterior,
                           aes(x = b_Intercept,
                               color = label, fill = label)) +

  geom_vline(xintercept = 0, linetype = 2) +
  geom_density(
    rel_min_height = 0.01,
    col = NA,
    scale = 1,
    alpha = 0.5
  ) +

  # Add text and labels
  geom_text(
    data = mutate_if(pooled_summary,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"),
      y = label,
      x = 1
    ),
    hjust = "inward",
    size = 3,
    color = "black"
  ) +

  scale_color_manual(values = c("#009E73", "#E69F00"), limits = rev) +
  scale_fill_manual(values = c("#009E73", "#E69F00"), limits = rev) +

  labs(x = "Standardised Mean Difference (Positive Values Favour Self-Talk)", # summary measure
       y = element_blank(),
       fill = "",
       title = "Updated Posterior Pooled Estimate",
       subtitle = "Note: x-axis rescaled from panel (A) for easier comparison of prior and posterior distributions") +
  scale_x_continuous(limits = c(-0.5, 1), breaks = c(-0.5, 0, 0.5, 1)) +
  theme_classic() +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA),
        plot.subtitle = element_text(size=6))

(forest_study / posterior_update) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  "plots/bayesian_forest_plot.png",
  width = 7.5,
  height = 10,
  device = "png",
  dpi = 300
)


##### subgroup and moderators models

get_prior(yi | se(1 / vi) ~ 1 + motor_demands + (1 | study / experiment / group / effect), data=data_effect_sizes)





motor_demands_model <-
  brm(
    yi | se(sqrt(vi)) ~ 1 + motor_demands + (1 | study / experiment / group / effect),
    data = data_effect_sizes,
    prior = prior_meta,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    control = list(adapt_delta = 0.99)
  )


plot(motor_demands_model)
