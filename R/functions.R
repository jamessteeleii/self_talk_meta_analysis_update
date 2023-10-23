# Read in and prepare data for effect size calculations and models
read_prepare_data <- function(file) {
  data <- read_csv(here("data", "Final data.csv")) %>%
    mutate_at(c(2:8, 30:32, 34:37), as.factor) %>%
    clean_names() %>%

    # add in pre-post correlation assumption
    mutate(ri = 0.7) %>%

    # add calculation of pooled baseline standard deviations
    mutate(pre_sd_pool = sqrt(((n_st - 1) * pre_sd_st ^ 2 +
                                 (n_con - 1) * pre_sd_con ^ 2
    ) /
      (n_st + n_con - 2)),
    post_sd_pool = sqrt(((n_st - 1) * post_sd_st ^ 2 +
                           (n_con - 1) * post_sd_con ^ 2
    ) /
      (n_st + n_con - 2)))
}

# Calculate all effect sizes and recombine data
calculate_effect_sizes <- function(data) {
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
}

# Misc functions
make_plot_tiff <- function(plot, width, height, path) {
  ggsave(
    path,
    plot,
    width = width,
    height = height,
    device = "tiff",
    dpi = 300
  )

}

get_tidy_model <- function(model) {
  tidy(model)
}

# Setup rstan to run quicker
rstan_setup <- function() {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)
}

# Priors
set_main_model_prior <- function() {

  # Taken form Hatsigeorgiadis et al. (2011) overall estimate
  main_model_prior <-
    prior("student_t(3, 0.48, 0.05)", class = "Intercept")
}

set_motor_demands_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  motor_demands_prior <-
    c(
      prior("student_t(3, 0.67, 0.07)", class = "b", coef = "motor_demandsFine"),
      prior("student_t(3, 0.26, 0.04)", class = "b", coef = "motor_demandsGross")
    )
}

prior_motor_demands_model <- function(data, prior) {
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
}



# Models
fit_main_model <- function(data, prior) {
  main_model <-
    brm(
      yi | se(sqrt(vi)) ~ 1 + (1 | study / experiment / group / effect),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99)
    )
}

plot_main_model <- function(data, model) {
  # Get study level draws and summarise
  study_draws <- model %>%
    spread_draws(b_Intercept, r_study[study,]) %>%
    mutate(b_Intercept = b_Intercept + r_study,
           study = as.factor(study)) %>%
    select(.chain, .iteration, .draw, b_Intercept, study) %>%
    merge(data[, 1:2], by = "study")

  study_summary <- group_by(study_draws, label) %>%
    mean_qi(b_Intercept) %>%
    mutate(b_Intercept = b_Intercept,
           .lower = .lower,
           .upper = .upper)

  # Get pooled estimate draws and summarise
  pooled_draws <- model %>%
    spread_draws(b_Intercept) %>%
    mutate(study = "Posterior Pooled Estimate",
           label = "Posterior Pooled Estimate")

  pooled_summary <- group_by(pooled_draws, label) %>%
    mean_qi(b_Intercept) %>%
    mutate(b_Intercept = b_Intercept,
           .lower = .lower,
           .upper = .upper)

  # Get prediction interval
  nd <- data.frame(study = "new", vi = 0)

  pred_int_data <- posterior_predict(
    object = model,
    newdata = nd,
    re_formula = NULL,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )

  pred_int_data <- median_qi(pred_int_data) %>%
    mutate(label = as.character("Posterior Pooled Estimate"))

  # Dataframe for the prior to plot
  prior_plot <-
    data.frame(b_Intercept = rstudent_t(
      n = 1e6,
      df = 3,
      mu = 0.48,
      sigma = 0.05
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Forest plot of study level estimates
  forest_study <- ggplot(aes(x = b_Intercept,
                             y = reorder(label, b_Intercept)),
                         data = study_draws) +

    # Add band for prediction interval
    annotate(
      "rect",
      xmin = pred_int_data$ymin,
      xmax = pred_int_data$ymax,
      ymin = -Inf,
      ymax = Inf,
      alpha = .1,
      fill = "black"
    ) +

    # Add ref line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    scale_y_discrete() +

    # Add densities and point intervals
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

    # Add individual study data
    geom_point(
      data = subset(data, !is.na(yi)),
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
    labs(
      x = element_blank(),
      # summary measure
      y = element_blank(),
      title = "Study level estimates",
      subtitle = "Posterior distributions, mean and 95% quantile intervals, individual effects (ticks), and 95% prediction interval (grey band)"
    ) +
    scale_x_continuous(limits = c(-2, 4),
                       breaks = c(-2,-1, 0, 1, 2, 3, 4)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

  # Combine prior dataframe with pooled draws
  prior_posterior <- rbind(pooled_draws[, c(4, 6)], prior_plot) %>%
    mutate(label = factor(
      label,
      levels = c(
        "Prior (Hatzigeorgiadis et al., 2011)",
        "Posterior Pooled Estimate"
      )
    ))

  # Plot showing updating of posterior estimate
  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = b_Intercept,
                                 color = label, fill = label)) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add densities
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

    labs(
      x = "Standardised Mean Difference (Positive Values Favour Self-Talk)",
      # summary measure
      y = element_blank(),
      fill = "",
      title = "Updated Posterior Pooled Estimate",
      subtitle = "Prior and posterior distributions for pooled estimates, and mean and 95% quantile interval for posterior (text label)\nNote: x-axis rescaled from panel (A) for easier comparison of prior and posterior distributions"
    ) +
    scale_x_continuous(limits = c(-0.5, 1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(fill = NA),
      plot.subtitle = element_text(size = 6)
    )

  # Combine plots
  (forest_study / posterior_update) +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(tag_levels = "A")

}

fit_motor_demands_model <- function(data, prior) {
  motor_demands_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + motor_demands + (1 | study / experiment / group / effect),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99)
    )
}

plot_motor_demands_model <- function(data, prior, prior_model, model) {

}


# Model checks
make_rhat_plot <- function(model) {
  mod_rhat <- enframe(brms::rhat(model)) %>%
    filter(!str_detect(name, "^r_id"))

  rhat_main_params <- mod_rhat$value

  mcmc_rhat(rhat_main_params) +
    scale_x_continuous(breaks = c(1, 1.01, 1.02, 1.03, 1.04, 1.05)) +
    geom_vline(xintercept = 1.01,
               linetype = "dashed",
               alpha = 0.25)
}

make_trace_plot <- function(model) {
  plot(model)
}

make_pp_check <- function(model) {
  pp_check(model)
}
