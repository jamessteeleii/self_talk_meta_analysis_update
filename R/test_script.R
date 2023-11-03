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

data <- data %>%
  mutate(
    # Convert p to t (Change scores)
    delta_t_value_st = replace_na(qt(delta_p_value_st/2, df=n_st-1, lower.tail=FALSE)),
    delta_t_value_con = replace_na(qt(delta_p_value_con/2, df=n_con-1, lower.tail=FALSE)),

    # Convert t to SE (Change scores)
    delta_se_st = replace_na(if_else(is.na(delta_m_st),
                                    (post_m_st - pre_m_st)/delta_t_value_st, delta_m_st/delta_t_value_st)),
    delta_se_con = replace_na(if_else(is.na(delta_m_con),
                                    (post_m_con - pre_m_con)/delta_t_value_con, delta_m_con/delta_t_value_con)),

    # Make positive
    delta_se_st = if_else(delta_se_st < 0, delta_se_st * -1, delta_se_st),
    delta_se_con = if_else(delta_se_con < 0, delta_se_con * -1, delta_se_con),

    # Convert SE to SD (Change scores)
    delta_sd_st = replace_na(delta_se_st * sqrt(n_st)),
    delta_sd_con = replace_na(delta_se_con * sqrt(n_con)),

    # Add missing deltas
    delta_m_st = replace_na(post_m_st - pre_m_st),
    delta_m_con = replace_na(post_m_con - pre_m_con),

    # Calculate pre-post correlation coefficient for those with pre, post, and delta SDs
    ri_st = replace_na((pre_sd_st^2 + post_sd_st^2 - delta_sd_st^2)/(2 * pre_sd_st * post_sd_st)),
    ri_con = replace_na((pre_sd_con^2 + post_sd_con^2 - delta_sd_con^2)/(2 * pre_sd_con * post_sd_con)),

    # Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
    ri_st = if_else(between(ri_st,-1,1) == FALSE, NA, ri_st),
    ri_con = if_else(between(ri_con,-1,1) == FALSE, NA, ri_con),

    # Add 0.7 assumed correlation where missing
    ri_st = replace_na(ri_st, 0.7),
    ri_con = replace_na(ri_con, 0.7)
  )


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


prior_main_model <-
  brm(
    yi | se(sqrt(vi)) ~ 1 + (1 | study / experiment / group / effect),
    data = data_effect_sizes,
    prior = prior_meta,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    control = list(adapt_delta = 0.99),
    sample_prior = "only",
  )


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



prior_preds <- predictions(prior_main_model, type = "response",
                           # newdata = nd,
                           re_formula = NA) %>%
  posterior_draws() %>%
  transform(type = "Response") %>%
  mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")


# Posterior distribution samples
posterior_preds <- predictions(main_model, type = "response",
                               # newdata = nd,
                               re_formula = NA) %>%
  posterior_draws() %>%
  transform(type = "Response") %>%
  mutate(label = "Posterior Pooled Estimate")

posterior_summary <- group_by(posterior_preds, label) %>%
  mean_qi(draw) %>%
  mutate(draw = draw,
         .lower = .lower,
         .upper = .upper)

# Combine prior and posterior samples for plot
prior_posterior <- rbind(prior_preds, posterior_preds) %>%
  mutate(label = factor(label, levels = c("Prior (Hatzigeorgiadis et al., 2011)",
                                          "Posterior Pooled Estimate")))


#
#
#
#
prior_draws <-
  prior_main_model %>%
  spread_draws(b_Intercept) %>%
    mutate(study = "Prior (Hatzigeorgiadis et al., 2011)",
           label = "Prior (Hatzigeorgiadis et al., 2011)")


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

posterior_draws <- main_model %>%
  spread_draws(b_Intercept) %>%
  mutate(study = "Posterior Pooled Estimate",
         label = "Posterior Pooled Estimate")

posterior_summary <- group_by(pooled_draws, label) %>%
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

prior_posterior <- rbind(posterior_draws[, c(4, 6)], prior_draws[, c(4, 6)]) %>%
  mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                          "Prior (Hatzigeorgiadis et al., 2011)"
                                          )))

posterior_update <- ggplot(data = prior_posterior,
                           aes(x = b_Intercept,
                               color = label, fill = label)) +

  geom_vline(xintercept = 0, linetype = 2) +

  stat_slab(alpha=0.6, linewidth=0) +

  # Add text and labels
  geom_text(
    data = mutate_if(posterior_summary,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("Posterior Pooled Estimate\n{b_Intercept} [{.lower}, {.upper}]"),
      y = 0.1,
      x = 1
    ),
    hjust = "inward",
    size = 3,
    color = "black"
  ) +

  scale_y_discrete(expand = c(0, 0)) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +

  labs(x = "Standardised Mean Difference (Positive Values Favour Self-Talk)", # summary measure
       y = element_blank(),
       fill = "",
       title = "Updated Posterior Pooled Estimate",
       subtitle = "Note: x-axis rescaled from panel (A) for easier comparison of prior and posterior distributions") +
  scale_x_continuous(limits = c(-0.5, 1), breaks = c(-0.5, 0, 0.5, 1)) +
  theme_classic() +
  guides(color = "none") +
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

get_prior(yi | se(1 / vi) ~ 0 + motor_demands + (1 | study / experiment / group / effect), data=data_effect_sizes)

motor_demands_prior <-
  c(
    prior("student_t(3, 0.67, 0.07)", class = "b", coef = "motor_demandsFine"),
    prior("student_t(3, 0.26, 0.04)", class = "b", coef = "motor_demandsGross")
  )



prior_motor_demands_model <-
  brm(
    yi | se(sqrt(vi)) ~ 0 + motor_demands + (1 | study / experiment / group / effect),
    data = data_effect_sizes,
    prior = motor_demands_prior,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    control = list(adapt_delta = 0.99),
    sample_prior = "only"
  )



plot(motor_demands_model)
pp_check(motor_demands_model)



# Prior distribution samples
nd <- datagrid(model = motor_demands_model,
         motor_demands = c("Fine", "Gross")
         )

prior_preds <- predictions(prior_motor_demands_model, type = "response",
            newdata = nd,
            re_formula = NA) %>%
  posterior_draws() %>%
  transform(type = "Response") %>%
  mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

# Posterior distribution samples
posterior_preds <- predictions(motor_demands_model, type = "response",
                           newdata = nd,
                           re_formula = NA) %>%
  posterior_draws() %>%
  transform(type = "Response") %>%
  mutate(label = "Posterior Pooled Estimate")

posterior_summary <- group_by(posterior_preds, motor_demands, label) %>%
  mean_qi(draw) %>%
  mutate(draw = draw,
         .lower = .lower,
         .upper = .upper)

# Combine prior and posterior samples for plot
prior_posterior <- rbind(prior_preds, posterior_preds) %>%
  mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                          "Prior (Hatzigeorgiadis et al., 2011)"
                                          )))

posterior_update <- ggplot(data = prior_posterior,
                           aes(x = draw, y = motor_demands,
                               color = label, fill = label)
                           ) +
  # Add reference line at zero
  geom_vline(xintercept = 0, linetype = 2) +

  stat_slab(alpha = 0.6, linewidth = 0) +

  # Add individual study data
  geom_point(
    data = data_effect_sizes %>%
      filter(!is.na(yi)) %>%
      mutate(label = NA),
    aes(x = yi, y = motor_demands),
    position = position_nudge(y = -0.05),
    shape = "|"
  ) +

  # Add text and labels
  geom_text(
    data = mutate_if(posterior_summary,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      y = motor_demands,
      x = 1.5
    ),
    hjust = "inward",
    size = 3,
    color = "black",
    position = position_nudge(y=0.1)
  ) +

  scale_color_manual(values = c("#009E73", "#E69F00")) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +

  labs(x = "Standardised Mean Difference", # summary measure
       y = element_blank(),
       fill = "",
       title = "Motor Demands (Fine versus Gross)",
       subtitle = "Prior and posterior distributions for pooled estimates, individual effects (ticks), and mean and 95% quantile interval for posterior (text label)"
  ) +
  scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
  guides(color = "none",
         shape = "none") +
  theme_classic() +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA),
        title = element_text(size=8),
        plot.subtitle = element_text(size=6))




















data_effect_sizes_study_design <- data_effect_sizes %>%
  mutate(study_design = if_else(study_design == "between", "Pre/post - experimental/control",
                                if_else(study_design == "between-post", "Post - experimental/control", "Pre/post - experimental" )))


get_prior(yi | se(1 / vi) ~ 0 + study_design + (1 | study / experiment / group / effect), data=data_effect_sizes_study_design)

study_design_prior <-
  c(
    prior("student_t(3, 0.37, 0.09)", class = "b", coef = "study_designPostMexperimentalDcontrol"),
    prior("student_t(3, 0.36, 0.1)", class = "b", coef = "study_designPreDpostMexperimental"),
    prior("student_t(3, 0.53, 0.06)", class = "b", coef = "study_designPreDpostMexperimentalDcontrol")

  )



prior_study_design_model <-
  brm(
    yi | se(sqrt(vi)) ~ 0 + study_design + (1 | study / experiment / group / effect),
    data = data_effect_sizes_study_design,
    prior = study_design_prior,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    control = list(adapt_delta = 0.99),
    sample_prior = "only"
  )

study_design_model <-
  brm(
    yi | se(sqrt(vi)) ~ 0 + study_design + (1 | study / experiment / group / effect),
    data = data_effect_sizes_study_design,
    prior = study_design_prior,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    control = list(adapt_delta = 0.99),
  )



plot(study_design_model)
pp_check(study_design_model)


# Prior distribution samples
nd <- datagrid(model = study_design_model,
               study_design = c("Pre/post - experimental/control", "Post - experimental/control", "Pre/post - experimental")
)

prior_preds <- predictions(prior_study_design_model, type = "response",
                           newdata = nd,
                           re_formula = NA) %>%
  posterior_draws() %>%
  transform(type = "Response") %>%
  mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

# Posterior distribution samples
posterior_preds <- predictions(study_design_model, type = "response",
                               newdata = nd,
                               re_formula = NA) %>%
  posterior_draws() %>%
  transform(type = "Response") %>%
  mutate(label = "Posterior Pooled Estimate")

posterior_summary <- group_by(posterior_preds, study_design, label) %>%
  mean_qi(draw) %>%
  mutate(draw = draw,
         .lower = .lower,
         .upper = .upper)

# Combine prior and posterior samples for plot
prior_posterior <- rbind(prior_preds, posterior_preds) %>%
  mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                          "Prior (Hatzigeorgiadis et al., 2011)"
  )))

posterior_update <- ggplot(data = prior_posterior,
                           aes(x = draw, y = study_design,
                               color = label, fill = label)
) +
  # Add reference line at zero
  geom_vline(xintercept = 0, linetype = 2) +

  stat_slab(alpha = 0.6, linewidth = 0) +

  # Add individual study data
  geom_point(
    data = data_effect_sizes_study_design %>%
      filter(!is.na(yi)) %>%
      mutate(label = NA),
    aes(x = yi, y = study_design),
    position = position_nudge(y = -0.05),
    shape = "|"
  ) +

  # Add text and labels
  geom_text(
    data = mutate_if(posterior_summary,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      y = study_design,
      x = 1.5
    ),
    hjust = "inward",
    size = 3,
    color = "black",
    position = position_nudge(y=0.1)
  ) +

  # scale_y_discrete(limits = c("Instructional",
  #                             "Motivational",
  #                             "Combined Instructional/Motivational",
  #                             "Rational")) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +

  labs(x = "Standardised Mean Difference", # summary measure
       y = element_blank(),
       fill = "",
       title = "Study Design (Pre/Post or Post Only With Experimental Only or Control)",
       subtitle = "Prior and posterior distributions for pooled estimates, individual effects (ticks), and mean and 95% quantile interval for posterior (text label)"
  ) +
  scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
  guides(color = "none",
         shape = "none") +
  theme_classic() +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA),
        title = element_text(size=8),
        plot.subtitle = element_text(size=6))







targets::tar_load(
  # motor_demands_model_plot,
  # participant_group_model_plot
  # selftalk_content_model_plot
  # matching_model_plot
  # task_novelty_model_plot
  # cue_selection_model_plot
  # training_model_plot
  study_design_model_plot
)


library(patchwork)

( (motor_demands_model_plot / participant_group_model_plot / selftalk_content_model_plot) |
  (matching_model_plot / task_novelty_model_plot / cue_selection_model_plot) |
  (training_model_plot / study_design_model_plot / plot_spacer()) ) +
  plot_annotation(tag_levels = "A",
                  title = "Updated Posterior Moderators Estimates",
                  subtitle = "Prior and posterior distributions for pooled estimates, individual effects (ticks), and mean and 95% quantile interval for posterior (text label)") +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

patchwork::plot_spacer()



# BF <- bayestestR::bayesfactor_parameters(main_model, null = 0)
#
# BF$
#
# si <- bayestestR::si(
#   posterior = main_model,
#   prior = NULL,
#   BF = ,
#   verbose = FALSE
# )
#
# plot(si)
#
# BF_gradient <- function(model, effects) {
#   bayestestR::bayesfactor_parameters(model, null = effects)
#
#   BF$log_BF
# }
#
#
#


install.packages("doSNOW")

install.packages("doParallel")

install.packages("doMPI")

install.packages("foreach")

#
library(foreach)
library(doParallel)

n.cores <- parallel::detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#ensure stopped on exit
on.exit(stopCluster(my.cluster))



system.time(

logBF_curve <- foreach (i = seq(0,1,length=100), .packages = c("bayestestR"), .combine = "rbind") %dopar% {
  data.frame(
    effect = i,
    BF = exp(bayesfactor_parameters(main_model, null = i)$log_BF)
  )
}

)


BF_curve_plot <- logBF_curve %>% ggplot(aes(x=effect, y=log10(BF))) +
  geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
  annotate("text", label = stringr::str_wrap("Negative", width = 10),
           x = 1, y=-0.25, size = 3) +
  annotate("text", label = stringr::str_wrap("Weak", width = 10),
           x = 1, y=0.25, size = 3) +
  annotate("text", label = stringr::str_wrap("Substantial", width = 10),
           x = 1, y=0.75, size = 3) +
  annotate("text", label = stringr::str_wrap("Strong", width = 10),
           x = 1, y=1.25, size = 3) +
  annotate("text", label = stringr::str_wrap("Very Strong", width = 10),
           x = 1, y=1.75, size = 3) +
  annotate("text", label = stringr::str_wrap("Decisive", width = 10),
           x = 1, y=2.25, size = 3) +
  geom_smooth(se=FALSE) +
  labs(x = "Point Effect Size BF Calculated Against",
       y = "log10(BF)",
       title = "Change in Evidence") +
  scale_y_continuous(limits = c(-0.5,5)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA))



  targets::tar_load(main_model_plot)

  library(patchwork)


  main_model_plot | BF_curve_plot







  plan(cluster, workers = 10)

  BF_curve <- future_map_dfr(.x = seq(0,1,length=10), .f = function(.x) {
    return(data.frame(effect = .x,
                      BF = bayesfactor_parameters(motor_demands_model, null = .x)))
  })

  plan(sequential)











