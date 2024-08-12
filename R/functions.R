# Read in and prepare data for effect size calculations and models
read_prepare_data <- function(file) {
  data <- read_csv(here("data", "Final data.csv")) |>
    mutate_at(c(2:9, 41:43, 45:48), as.factor) |>
    clean_names() |>

    # Convert SE to SD
    mutate(
      pre_sd_st = replace_na(if_else(is.na(pre_sd_st),
                                     pre_se_st * sqrt(n_st), pre_sd_st)),
      pre_sd_con = replace_na(if_else(is.na(pre_sd_con),
                                     pre_se_con * sqrt(n_con), pre_sd_con)),
      post_sd_st = replace_na(if_else(is.na(post_sd_st),
                                     post_se_st * sqrt(n_st), post_sd_st)),
      post_sd_con = replace_na(if_else(is.na(post_sd_con),
                                      post_se_con * sqrt(n_con), post_sd_con))
    ) |>

    # Add in pre-post correlations
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
      delta_sd_st = replace_na(if_else(is.na(delta_sd_st),
                                       delta_se_st * sqrt(n_st), delta_sd_st)),
      delta_sd_con = replace_na(if_else(is.na(delta_sd_con),
                                        delta_se_con * sqrt(n_con), delta_sd_con)),

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

    )  |>

    # add calculation of pooled baseline standard deviations
    mutate(pre_sd_pool = sqrt(((n_st - 1) * pre_sd_st ^ 2 +
                                 (n_con - 1) * pre_sd_con ^ 2
    ) /
      (n_st + n_con - 2)),
    post_sd_pool = sqrt(((n_st - 1) * post_sd_st ^ 2 +
                           (n_con - 1) * post_sd_con ^ 2
    ) /
      (n_st + n_con - 2))) |>

    # Formatting variables for moderator models
    unite(c("selftalk_content", "motor_demands"), col="matching", sep = "/",
          remove = FALSE) |>


    separate(novelty, into=c("selftalk_novelty1", "task_novelty1"), sep = '\\s*and\\s*',
             remove = FALSE) |>
    separate(novelty, into=c("selftalk_novelty2", "task_novelty2"), sep = '\\s*but\\s*',
             remove = FALSE) |>
    mutate(task_novelty = if_else(task_novelty2 == "with the sport", "Learned", "Novel", missing = "Novel")) |>
    select(-selftalk_novelty1, -selftalk_novelty2, -task_novelty1, -task_novelty2) |>

    mutate(cue_selection = if_else(cue_selection == "No", "Assigned", "Self-selected"),
           overtness_selection = if_else(overtness_selection == "No", "Assigned",
                                         if_else(overtness_selection == "Yes", "Self-selected", NA)),
           training = if_else(acute_chronic == "acute", "No training", "Training"),
           study_design = if_else(study_design == "between", "Pre/post - experimental/control",
                                  if_else(study_design == "between-post", "Post - experimental/control", "Pre/post - experimental" )))

  }

# Calculate all effect sizes and recombine data
calculate_effect_sizes <- function(data) {
  ### For effects where increase is good
  data_increase <- subset(data, improvement == "Increase")

  # Calculate pre-post control effect sizes
  data_increase_ppc <-
    subset(data_increase, study_design == "Pre/post - experimental/control")

  data_increase_ppc_st <- escalc(
    measure = "SMCR",
    m1i = post_m_st,
    m2i = pre_m_st,
    sd1i = pre_sd_pool,
    ni = n_st,
    ri = ri_st,
    data = data_increase_ppc
  )
  data_increase_ppc_con <- escalc(
    measure = "SMCR",
    m1i = post_m_con,
    m2i = pre_m_con,
    sd1i = pre_sd_pool,
    ni = n_con,
    ri = ri_con,
    data = data_increase_ppc
  )

  data_increase_ppc$yi <-
    (data_increase_ppc_st$yi - data_increase_ppc_con$yi)
  data_increase_ppc$vi <-
    (data_increase_ppc_st$vi + data_increase_ppc_con$vi)

  # Calculate post only effect sizes
  data_increase_post <-
    subset(data_increase, study_design == "Post - experimental/control")

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
    subset(data_increase, study_design == "Pre/post - experimental")

  data_increase_pp <- escalc(
    measure = "SMCR",
    m1i = post_m_st,
    m2i = pre_m_st,
    sd1i = pre_sd_st,
    ni = n_st,
    ri = ri_st,
    data = data_increase_pp
  )

  ### For effects where decrease is good
  data_decrease <- subset(data, improvement == "Decrease")

  # Calculate pre-post control effect sizes
  data_decrease_ppc <-
    subset(data_decrease, study_design == "Pre/post - experimental/control")

  data_decrease_ppc_st <- escalc(
    measure = "SMCR",
    m1i = pre_m_st,
    m2i = post_m_st,
    sd1i = pre_sd_pool,
    ni = n_st,
    ri = ri_st,
    data = data_decrease_ppc
  )
  data_decrease_ppc_con <- escalc(
    measure = "SMCR",
    m1i = pre_m_con,
    m2i = post_m_con,
    sd1i = pre_sd_pool,
    ni = n_con,
    ri = ri_con,
    data = data_decrease_ppc
  )

  data_decrease_ppc$yi <-
    (data_decrease_ppc_st$yi - data_decrease_ppc_con$yi)
  data_decrease_ppc$vi <-
    (data_decrease_ppc_st$vi + data_decrease_ppc_con$vi)

  # Calculate post only effect sizes
  data_decrease_post <-
    subset(data_decrease, study_design == "Post - experimental/control")

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
    subset(data_decrease, study_design == "Pre/post - experimental")

  data_decrease_pp <- escalc(
    measure = "SMCR",
    m1i = pre_m_st,
    m2i = post_m_st,
    sd1i = pre_sd_st,
    ni = n_st,
    ri = ri_st,
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

ggplot_2_grob_nback <- function(plot) {
  plot_grob <- ggplot2::ggplotGrob(plot)
  plot_new <- ggpubr::as_ggplot(plot_grob)

  remove(plot)
  remove(plot_grob)

  gc()

  return(plot_new)
}

get_logBF_curve <- function(model) {

  plan(multisession, workers = 5)

  BF_curve <- future_map_dfr(.x = seq(0,1,length=100), .f = function(.x) {
    return(data.frame(effect = .x,
                      BF = bayesfactor_parameters(model, null = .x)))
  })



  return(BF_curve)
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
set_main_prior <- function() {

  # Taken fromm Hatsigeorgiadis et al. (2011) overall estimate
  main_model_prior <-
    c(
      prior("student_t(60, 0.48, 0.05)", class = "Intercept")
    )
}

sample_prior_main_model <- function(data, prior) {



  prior_main_model <-
    brm(
      yi | se(sqrt(vi)) ~ 1 + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      future = TRUE,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only",
    )



  prior_main_model
}

set_motor_demands_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  motor_demands_prior <-
    c(
      prior("student_t(35, 0.67, 0.07)", class = "b", coef = "motor_demandsFine"),
      prior("student_t(25, 0.26, 0.04)", class = "b", coef = "motor_demandsGross")
    )
}

sample_prior_motor_demands_model <- function(data, prior) {
  prior_motor_demands_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + motor_demands + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

set_participant_group_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  participant_group_prior <-
    c(
      prior("student_t(30, 0.5, 0.08)", class = "b", coef = "participant_groupNonMathletes"),
      prior("student_t(20, 0.47, 0.08)", class = "b", coef = "participant_groupBeginnerathletes"),
      prior("student_t(6, 0.38, 0.14)", class = "b", coef = "participant_groupExperiencedathletes")
    )
}

sample_prior_participant_group_model <- function(data, prior) {
  prior_participant_group_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + participant_group + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

set_selftalk_content_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  # New groups are given the overall prior for the main model
  selftalk_content_prior <-
    c(
      prior("student_t(36, 0.55, 0.08)", class = "b", coef = "selftalk_contentInstructional"),
      prior("student_t(18, 0.37, 0.06)", class = "b", coef = "selftalk_contentMotivational"),
      prior("student_t(3, 0.48, 0.05)", class = "b", coef = "selftalk_contentCombinedInstructionalDMotivational"),
      prior("student_t(3, 0.48, 0.05)", class = "b", coef = "selftalk_contentRational")

    )
}

sample_prior_selftalk_content_model <- function(data, prior) {
  prior_selftalk_content_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + selftalk_content + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

set_matching_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  matching_prior <-
    c(
      prior("student_t(22, 0.83, 0.1)", class = "b", coef = "matchingInstructionalDFine"),
      prior("student_t(12, 0.22, 0.04)", class = "b", coef = "matchingInstructionalDGross"),
      prior("student_t(8, 0.41, 0.09)", class = "b", coef = "matchingMotivationalDFine"),
      prior("student_t(8, 0.33, 0.08)", class = "b", coef = "matchingMotivationalDGross")

    )
}

sample_prior_matching_model <- function(data, prior) {

  data <- data |>
    filter(matching == "Motivational/Fine" |
             matching == "Motivational/Gross" |
             matching == "Instructional/Fine" |
             matching == "Instructional/Gross")

  prior_matching_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + matching + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

set_task_novelty_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  task_novelty_prior <-
    c(
      prior("student_t(13, 0.41, 0.05)", class = "b", coef = "task_noveltyLearned"),
      prior("student_t(45, 0.73, 0.14)", class = "b", coef = "task_noveltyNovel")

    )
}

sample_prior_task_novelty_model <- function(data, prior) {

  prior_task_novelty_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + task_novelty + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

set_cue_selection_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  cue_selection_prior <-
    c(
      prior("student_t(44, 0.49, 0.07)", class = "b", coef = "cue_selectionAssigned"),
      prior("student_t(14, 0.44, 0.07)", class = "b", coef = "cue_selectionSelfMselected")

    )
}

sample_prior_cue_selection_model <- function(data, prior) {

  prior_cue_selection_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + cue_selection + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

set_overtness_selection_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  overtness_selection_prior <-
    c(
      prior("student_t(28, 0.49, 0.08)", class = "b", coef = "overtness_selectionAssigned"),
      prior("student_t(25, 0.48, 0.08)", class = "b", coef = "overtness_selectionSelfMselected")

    )
}

sample_prior_overtness_selection_model <- function(data, prior) {

  prior_overtness_selection_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + overtness_selection + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

set_training_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  training_prior <-
    c(
      prior("student_t(37, 0.37, 0.04)", class = "b", coef = "trainingNotraining"),
      prior("student_t(21, 0.8, 0.12)", class = "b", coef = "trainingTraining")

    )
}

sample_prior_training_model <- function(data, prior) {

  prior_training_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + training + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

set_study_design_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  study_design_prior <-
    c(
      prior("student_t(3, 0.37, 0.09)", class = "b", coef = "study_designPostMexperimentalDcontrol"),
      prior("student_t(14, 0.36, 0.1)", class = "b", coef = "study_designPreDpostMexperimental"),
      prior("student_t(33, 0.53, 0.06)", class = "b", coef = "study_designPreDpostMexperimentalDcontrol")

    )
}

sample_prior_study_design_model <- function(data, prior) {

  prior_study_design_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + study_design + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only"
    )
}

# Simulation of additional new study
additional_new_study_sim <- function(prior_model) {

  sim <- function(participant_n = as.double(),
                  effect_size = as.double(),
                  ... # helps the function work with pmap() below
  ) {

    # Taken fromm Hatsigeorgiadis et al. (2011) overall estimate
    main_model_prior <-
      c(
        brms::prior("student_t(60, 0.48, 0.05)", class = "Intercept")
      )

    effect_size_dat <- escalc(
      measure = "SMD1",
      m1i = effect_size,
      m2i = 0,
      sd2i = 1,
      n1i = participant_n/2,
      n2i = participant_n/2
    )

    effect_size_dat$study_code <- 1

    main_model <-
      brm(
        yi | se(sqrt(vi)) ~ 1 + (1 | study_code),
        data = effect_size_dat,
        prior = main_model_prior,
        chains = 4,
        cores = 4,
        seed = 1988,
        warmup = 2000,
        iter = 6000,
        control = list(adapt_delta = 0.99)
      )

    percentage_in_rope <- rope(main_model, ci = 1, range = c(0.38, 0.58))[1,5]

    # Sample draws from the prior distribution
    prior_draws <-
      prior_model |>
      spread_draws(b_Intercept) |>
      mutate(study = "Prior (Hatzigeorgiadis et al., 2011)",
             label = "Prior (Hatzigeorgiadis et al., 2011)")


    # Sample draws from the posterior distribution
    posterior_draws <- main_model |>
      spread_draws(b_Intercept) |>
      mutate(study = "Posterior Pooled Estimate",
             label = "Posterior Pooled Estimate")

    posterior_summary <- group_by(posterior_draws, label) |>
      mean_qi(b_Intercept) |>
      mutate(b_Intercept = b_Intercept,
             .lower = .lower,
             .upper = .upper)

    # Combine prior dataframe with pooled draws
    prior_posterior <- rbind(posterior_draws[, c(4, 6)], prior_draws[, c(4, 6)]) |>
      mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                              "Prior (Hatzigeorgiadis et al., 2011)"
      ))) |>
      mutate(percentage_in_rope = percentage_in_rope,
             yi = effect_size_dat$yi,
             vi = effect_size_dat$vi)

    remove(main_model)

    gc()

    return(prior_posterior)

  }


  plan(multisession, workers = 4)

  set.seed(1988)

  sims <- crossing(
    # rep = 1:1000, # number of replicates
    participant_n = c(10,20,40,80,160,320,640,1282,2560,5120), # range of participant N
    effect_size = seq(0,1, by = 0.2) # true effect size
  ) %>% # not sure why base pipe doesn't work here
    mutate(analysis = future_pmap(., sim)) %>% # not sure why base pipe doesn't work here
    unnest(analysis)

  plan(sequential)

  sims
}

plot_additional_new_study_sims <- function(sims) {
  sims_summary <- sims |>
    group_by(participant_n, effect_size) |>
    slice(1) |>
    mutate(effect_lab = "Effect Size")

  # Plot showing updating of posterior estimate from a single new study
  posterior_update_new_study <- sims |>
    mutate(effect_lab = "Effect Size") |>
    ggplot(aes(x = b_Intercept, y = factor(participant_n),
               color = label, fill = label)) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add densities
    stat_slab(alpha=0.6, linewidth=0) +

    # Add text and labels
    geom_text(
      data = mutate_if(sims_summary,
                       is.numeric, round, 3),
      aes(
        label = glue::glue("{scales::percent(percentage_in_rope)}"),
        x = 1.1
      ),
      hjust = "inward",
      size = 3,
      color = "black",
      position = position_nudge(y = 0.5)
    ) +
    ggh4x::facet_nested(. ~ effect_lab + effect_size) +

    scale_y_discrete(expand = c(0, 0)) +
    scale_color_manual(values = c("#009E73", "#E69F00")) +
    scale_fill_manual(values = c("#009E73", "#E69F00")) +

    labs(
      x = "Standardised Mean Difference (Positive Values Favour Self-Talk)",
      # summary measure
      y = "Sample Size (n)",
      fill = "",
      title = "Examining the effects of a single new trial varying the true effect size and sample size",
      subtitle = "Prior and posterior distributions for pooled estimates and percentage of posterior distribution within 95% quantile interval of prior distribution (text label)"
    ) +
    guides(color = "none") +
    scale_x_continuous(limits = c(-0.05, 1.1),
                       breaks = c(0, 0.5, 1)) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(fill = NA)
    )

  posterior_update_new_study
}

# Models
fit_main_model <- function(data, prior) {



  main_model <-
    brm(
      yi | se(sqrt(vi)) ~ 1 + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      future = TRUE,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )



  main_model
}

fit_motor_demands_model <- function(data, prior) {
  motor_demands_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + motor_demands + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_participant_group_model <- function(data, prior) {
  participant_group_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + participant_group + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_selftalk_content_model <- function(data, prior) {
  selftalk_content_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + selftalk_content + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_matching_model <- function(data, prior) {

  data <- data |>
    filter(matching == "Motivational/Fine" |
             matching == "Motivational/Gross" |
             matching == "Instructional/Fine" |
             matching == "Instructional/Gross")

  matching_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + matching + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_task_novelty_model <- function(data, prior) {

  task_novelty_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + task_novelty + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_cue_selection_model <- function(data, prior) {

  cue_selection_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + cue_selection + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_overtness_selection_model <- function(data, prior) {

  data <- data |>
    filter(overtness_selection == "Assigned" |
             overtness_selection == "Self-selected")

  overtness_selection_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + overtness_selection + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_training_model <- function(data, prior) {

  training_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + training + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_study_design_model <- function(data, prior) {

  study_design_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + study_design + (1 | study_code / experiment_code / group_code / effect_code),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_cumulative_main_model <- function(data) {

  data <- data |>
    mutate(year = as.numeric(as.character(year)))

  # Create empty data frame for draws
  cumulative_draws <- data.frame(year = as.numeric(),
                                  b_Intercept = as.numeric(),
                                 prior = as.character())

  # Initial posterior as prior
  posterior <- c(0.48, 0.05, 60)

  # Loop through each model using previous years posterior as prior
  for (i in unique(data$year)){

    data_year <- data |>
      filter(year == i)

    prior <- set_prior(paste("student_t(",posterior[3],",", posterior[1],",", posterior[2],")"),
                       class = "Intercept")

    main_model <-
      brm(
        yi | se(sqrt(vi)) ~ 1 + (1 | study_code / experiment_code / group_code / effect_code),
        data = data_year,
        prior = prior,
        chains = 4,
        cores = 4,
        seed = 1988,
        warmup = 2000,
        iter = 6000,
        control = list(adapt_delta = 0.99)
      )

    draws <- main_model |>
      spread_draws(b_Intercept)

    posterior <- MASS::fitdistr(draws$b_Intercept, "t")$estimate

    cumulative_draws <- rbind(cumulative_draws,
                               data.frame(year = i,
                                          b_Intercept = draws$b_Intercept,
                                          prior = paste("student_t(",posterior[3],",", posterior[1],",", posterior[2],")"))
    )

  }

  return(cumulative_draws)

}

# Small study/Publication bias effects
fit_pet_model <- function(data) {


  pet_main_model_rma.mv <- rma.mv(yi, vi,
                                  random = ~ 1 | study_code / experiment_code / group_code / effect_code,
                                  mods = ~ sqrt(vi),
                                  data = data,
                                  method = "REML",
                                  test = "t"
  )
}

fit_rma.mv_model <- function(data) {

  rma.mv_model <- rma.mv(yi, vi,
                         random = ~ 1 | study_code / experiment_code / group_code / effect_code,
                         data = data,
                         method = "REML",
                         test = "t")
}

fit_null_robma_model <- function(data) {

  RoBMA::RoBMA.options(max_cores = 5)

  data_RoBMA <- data |>
    filter(!is.na(vi))

  null_RoBMA <- RoBMA::RoBMA(d = data_RoBMA$yi, v = data_RoBMA$vi,
                      parallel = TRUE,
                      seed = 1988)
}

fit_prior_robma_model <- function(data) {

  RoBMA::RoBMA.options(max_cores = 5)

  data_RoBMA <- data |>
    filter(!is.na(vi))

  prior_RoBMA <- RoBMA::RoBMA(d = data_RoBMA$yi, v = data_RoBMA$vi,
                       priors_effect = RoBMA::prior(distribution = "t", parameters = list(location = 0.48, scale = 0.05, df = 60)),
                       parallel = TRUE,
                       seed = 1988,
  )
}

# P-hacking effects
fit_p_hack_model <- function(data) {
  p_hack_model <- phma(yi, vi, filter(data, !is.na(yi)))
}

fit_classic_model <- function(data) {
  classic_model <- cma(yi, vi, filter(data, !is.na(yi)))
}

fit_p_hack_model_prior <- function(data) {
  p_hack_model_prior <- phma(yi, vi, filter(data, !is.na(yi)),
                             prior = list(theta0_mean = 0.48, theta0_sd = 0.05))
}

fit_classic_model_prior <- function(data) {
  classic_model_prior <- cma(yi, vi, filter(data, !is.na(yi)),
                             prior = list(theta0_mean = 0.48, theta0_sd = 0.05))
}


# Plots
plot_main_model_forest <- function(data, model) {

  # Sample draws from the individual studies
  study_draws <- model |>
    spread_draws(b_Intercept, r_study_code[study_code, ]) |>
    mutate(b_Intercept = b_Intercept + r_study_code,
           study_code = as.factor(study_code)) |>
    select(.chain, .iteration, .draw, b_Intercept, study_code) |>
    merge(data[, c(1,3)], by = "study_code")

  study_summary <- group_by(study_draws, label) |>
    mean_qi(b_Intercept) |>
    mutate(b_Intercept = b_Intercept,
           .lower = .lower,
           .upper = .upper)

  # Prediction interval
  nd <- data.frame(study = "new", vi = 0)

  pred_int_data <- posterior_predict(
    object = model,
    newdata = nd,
    re_formula = NULL,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )

  pred_int_data <- median_qi(pred_int_data) |>
    mutate(label = as.character("Posterior Pooled Estimate"))

  remove(model)

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

}

plot_contour_funnel <- function(data, tidy_model, pet_model, tidy_pet_model) {

  data <- data |>
    mutate(
      wi = 1/sqrt(vi),
      size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))

  contours <- data.frame(se.seq = seq(0, max(sqrt(data$vi), na.rm = TRUE), 0.001),
                         ll95 = 0-(1.96*seq(0, max(sqrt(data$vi), na.rm = TRUE), 0.001)),
                         ul95 = 0+(1.96*seq(0, max(sqrt(data$vi), na.rm = TRUE), 0.001)),
                         ll99 = 0-(3.29*seq(0, max(sqrt(data$vi), na.rm = TRUE), 0.001)),
                         ul99 = 0+(3.29*seq(0, max(sqrt(data$vi), na.rm = TRUE), 0.001)))

  pet_preds <- as.data.frame(predict(pet_model, newmods=seq(0, sqrt(max(data$vi, na.rm = TRUE)), length.out = 500))) |>
    mutate(x = seq(0, sqrt(max(data$vi, na.rm = TRUE)), length.out = 500))

  ggplot(aes(y = yi, x = sqrt(vi)), data = data) +
    geom_point(aes(size = size), alpha = 0.25) +
    geom_pointrange(aes(x = -0.05, y = estimate, ymin = conf.low, ymax = conf.high),
                    data = tidy_model[1,]) +
    annotate("text", y = tidy_model$estimate[1], x = -0.1, label = "Main Model\nEstimate", size = 2) +
    geom_pointrange(aes(x = -0.05, y = estimate, ymin = conf.low, ymax = conf.high),
                    data = tidy_pet_model[1,]) +
    annotate("text", y = tidy_pet_model$estimate[1], x = -0.1, label = "Multilevel PET\nEstimate", size = 2) +
    geom_segment(aes(y = 0, yend = 0, x = 0, xend = max(sqrt(vi), na.rm=TRUE)),
                 linetype = 3) +
    geom_line(aes(y = ll95, x = se.seq), linetype = 'dotted', data = contours) +
    geom_line(aes(y = ul95, x = se.seq), linetype = 'dotted', data = contours) +
    geom_line(aes(y = ll99, x = se.seq), linetype = 'dashed', data = contours) +
    geom_line(aes(y = ul99, x = se.seq), linetype = 'dashed', data = contours) +
    annotate("text", y = 1.75, x = 0.65, label = "italic(p) < 0.05",
             parse = TRUE, size = 2.5) +
    annotate("text", y = 2.5, x = 0.65, label = "italic(p) < 0.01",
             parse = TRUE, size = 2.5) +
    annotate("text", y = -1.75, x = 0.65, label = "italic(p) < 0.05",
             parse = TRUE, size = 2.5) +
    annotate("text", y = -2.5, x = 0.65, label = "italic(p) < 0.01",
             parse = TRUE, size = 2.5) +
    geom_ribbon(aes(x = x, y = pred, ymin = ci.lb, ymax = ci.ub),
                alpha = 0.25,
                data = pet_preds) +
    geom_line(aes(x = x, y = pred),
              data = pet_preds) +
    labs(
      y = "Standardised Mean Difference",
      x = "Standard Error",
      title = "Contour Enhanced Funnel Plot",
      subtitle = "Estimates from frequentist multilevel main model and PET-PEESE model (regression line and ribbon)\nNote, models only include studies since 2011"
    ) +
    guides(size = "none") +
    coord_flip() +
    scale_x_reverse() +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_p_hack_models <- function(p_hack_model, classic_model) {
  samples_p_hack_model <- tibble(
    draws = rstan::extract(p_hack_model)$theta0,
    model = "p-hacking model"
  )

  posterior_summary_p_hack_model <- samples_p_hack_model |>
    mean_qi(draws) |>
    mutate(model = "p-hacking model")

  samples_classic_model <- tibble(
    draws = rstan::extract(classic_model)$theta0,
    model = "Classic model"
  )

  posterior_summary_classic_model <- samples_classic_model |>
    mean_qi(draws) |>
    mutate(model = "Classic model")

  samples_models <- rbind(
    samples_classic_model, samples_p_hack_model
  )

  plot_p_hack_classic <- samples_models |>
    ggplot(aes(x = draws,
               color = model, fill = model)
    ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary_classic_model,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("Classic model without prior\n{draws} [{.lower}, {.upper}]"),
        y = 0.1,
        x = 1.5
      ),
      hjust = "inward",
      size = 3,
      color = "black",
      position = position_nudge(y=0.1)
    ) +

    geom_text(
      data = mutate_if(posterior_summary_p_hack_model,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("p-hacking model without prior\n{draws} [{.lower}, {.upper}]"),
        y = 0.1,
        x = -0.75
      ),
      hjust = "inward",
      size = 3,
      color = "black",
      position = position_nudge(y=0.1)
    ) +

    scale_color_manual(values = c("#0072B2", "#D55E00")) +
    scale_fill_manual(values = c("#0072B2", "#D55E00")) +

    labs(x = "Standardised Mean Difference", # summary measure
         y = element_blank(),
         fill = "",
         title = "Models without using prior (Hatzigeorgiadis et al., 2011)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_p_hack_models_prior <- function(p_hack_model_prior, classic_model_prior) {
  samples_p_hack_model_prior <- tibble(
    draws = rstan::extract(p_hack_model_prior)$theta0,
    model = "p-hacking model"
  )

  posterior_summary_p_hack_model_prior <- samples_p_hack_model_prior |>
    mean_qi(draws) |>
    mutate(model = "p-hacking model")

  samples_classic_model_prior <- tibble(
    draws = rstan::extract(classic_model_prior)$theta0,
    model = "Classic model"
  )

  posterior_summary_classic_model_prior <- samples_classic_model_prior |>
    mean_qi(draws) |>
    mutate(model = "Classic model")

  samples_prior_models <- rbind(
    samples_classic_model_prior, samples_p_hack_model_prior
  )

  plot_p_hack_classic_prior <- samples_prior_models |>
    ggplot(aes(x = draws,
               color = model, fill = model)
    ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary_classic_model_prior,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("Classic model with prior\n{draws} [{.lower}, {.upper}]"),
        y = 0.1,
        x = 1.5
      ),
      hjust = "inward",
      size = 3,
      color = "black",
      position = position_nudge(y=0.1)
    ) +

    geom_text(
      data = mutate_if(posterior_summary_p_hack_model_prior,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("p-hacking model with prior\n{draws} [{.lower}, {.upper}]"),
        y = 0.1,
        x = -0.75
      ),
      hjust = "inward",
      size = 3,
      color = "black",
      position = position_nudge(y=0.1)
    ) +

    scale_color_manual(values = c("#0072B2", "#D55E00")) +
    scale_fill_manual(values = c("#0072B2", "#D55E00")) +

    labs(x = "Standardised Mean Difference", # summary measure
         y = element_blank(),
         fill = "",
         title = "Models using prior (Hatzigeorgiadis et al., 2011)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_main_model_update <- function(prior_model, model) {
  # Sample draws from the prior distribution
  prior_draws <-
    prior_model |>
    spread_draws(b_Intercept) |>
    mutate(study = "Prior (Hatzigeorgiadis et al., 2011)",
           label = "Prior (Hatzigeorgiadis et al., 2011)")


  # Sample draws from the posterior distribution
  posterior_draws <- model |>
    spread_draws(b_Intercept) |>
    mutate(study = "Posterior Pooled Estimate",
           label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_draws, label) |>
    mean_qi(b_Intercept) |>
    mutate(b_Intercept = b_Intercept,
           .lower = .lower,
           .upper = .upper)

  # Combine prior dataframe with pooled draws
  prior_posterior <- rbind(posterior_draws[, c(4, 6)], prior_draws[, c(4, 6)]) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  # Plot showing updating of posterior estimate
  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = b_Intercept,
                                 color = label, fill = label)) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add densities
    stat_slab(alpha=0.6, linewidth=0) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("Posterior Pooled Estimate\n{b_Intercept} [{.lower}, {.upper}]"),
        y = 0.1,
        x = 1.1
      ),
      hjust = "inward",
      size = 3,
      color = "black"
    ) +

    scale_y_discrete(expand = c(0, 0)) +
    scale_color_manual(values = c("#009E73", "#E69F00")) +
    scale_fill_manual(values = c("#009E73", "#E69F00")) +

    labs(
      x = "Standardised Mean Difference (Positive Values Favour Self-Talk)",
      # summary measure
      y = element_blank(),
      fill = "",
      title = "Updated Posterior Pooled Estimate",
      subtitle = "Prior and posterior distributions for pooled estimates, and mean and 95% quantile interval for posterior (text label)\nNote: x-axis rescaled from panel (A) for easier comparison of prior and posterior distributions"
    ) +
    guides(color = "none") +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(fill = NA),
      plot.subtitle = element_text(size = 6)
    )

}

plot_BF_curve_main_model <- function(BF_curve) {

  # Bayes Factor curve plot
  BF_curve_plot <- BF_curve |> ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          axis.title.y = element_text(vjust = -40),
          plot.subtitle = element_text(size = 6))
}

plot_motor_demands_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_motor_demandsFine, b_motor_demandsGross) |>
    mutate(motor_demands = case_when(
      .variable == "b_motor_demandsFine" ~ "Fine",
      .variable == "b_motor_demandsGross" ~ "Gross"

    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_motor_demandsFine, b_motor_demandsGross) |>
    mutate(motor_demands = case_when(
      .variable == "b_motor_demandsFine" ~ "Fine",
      .variable == "b_motor_demandsGross" ~ "Gross"

    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, motor_demands, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = motor_demands,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = motor_demands),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
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
         title = "Motor Demands",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_participant_group_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_participant_groupBeginnerathletes, b_participant_groupExperiencedathletes, b_participant_groupNonMathletes) |>
    mutate(participant_group = case_when(
      .variable == "b_participant_groupBeginnerathletes" ~ "Beginner athletes",
      .variable == "b_participant_groupExperiencedathletes" ~ "Experienced athletes",
      .variable == "b_participant_groupNonMathletes" ~ "Non-athletes"
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_participant_groupBeginnerathletes, b_participant_groupExperiencedathletes, b_participant_groupNonMathletes) |>
    mutate(participant_group = case_when(
      .variable == "b_participant_groupBeginnerathletes" ~ "Beginner athletes",
      .variable == "b_participant_groupExperiencedathletes" ~ "Experienced athletes",
      .variable == "b_participant_groupNonMathletes" ~ "Non-athletes"
    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, participant_group, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = participant_group,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = participant_group),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
        y = participant_group,
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
         title = "Participants",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_selftalk_content_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_selftalk_contentCombinedInstructionalDMotivational,
                 b_selftalk_contentInstructional,
                 b_selftalk_contentMotivational,
                 b_selftalk_contentRational) |>
    mutate(selftalk_content = case_when(
      .variable == "b_selftalk_contentCombinedInstructionalDMotivational" ~ "Combined Instructional/Motivational",
      .variable == "b_selftalk_contentInstructional" ~ "Instructional",
      .variable == "b_selftalk_contentMotivational" ~ "Motivational",
      .variable == "b_selftalk_contentRational" ~ "Rational"
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_selftalk_contentCombinedInstructionalDMotivational,
                 b_selftalk_contentInstructional,
                 b_selftalk_contentMotivational,
                 b_selftalk_contentRational) |>
    mutate(selftalk_content = case_when(
      .variable == "b_selftalk_contentCombinedInstructionalDMotivational" ~ "Combined Instructional/Motivational",
      .variable == "b_selftalk_contentInstructional" ~ "Instructional",
      .variable == "b_selftalk_contentMotivational" ~ "Motivational",
      .variable == "b_selftalk_contentRational" ~ "Rational"
    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, selftalk_content, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = selftalk_content,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = selftalk_content),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
        y = selftalk_content,
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
         title = "Self-Talk Content",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_matching_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_matchingInstructionalDFine,
                 b_matchingInstructionalDGross,
                 b_matchingMotivationalDFine,
                 b_matchingMotivationalDGross) |>
    mutate(matching = case_when(
      .variable == "b_matchingInstructionalDFine" ~ "Instructional/Fine",
      .variable == "b_matchingInstructionalDGross" ~ "Instructional/Gross",
      .variable == "b_matchingMotivationalDFine" ~ "Motivational/Fine",
      .variable == "b_matchingMotivationalDGross" ~ "Motivational/Gross"
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_matchingInstructionalDFine,
                 b_matchingInstructionalDGross,
                 b_matchingMotivationalDFine,
                 b_matchingMotivationalDGross) |>
    mutate(matching = case_when(
      .variable == "b_matchingInstructionalDFine" ~ "Instructional/Fine",
      .variable == "b_matchingInstructionalDGross" ~ "Instructional/Gross",
      .variable == "b_matchingMotivationalDFine" ~ "Motivational/Fine",
      .variable == "b_matchingMotivationalDGross" ~ "Motivational/Gross"
    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, matching, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = matching,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = matching),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
        y = matching,
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
         title = "Matching Hypothesis",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_task_novelty_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_task_noveltyLearned,
                 b_task_noveltyNovel) |>
    mutate(task_novelty = case_when(
      .variable == "b_task_noveltyLearned" ~ "Learned",
      .variable == "b_task_noveltyNovel" ~ "Novel"
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_task_noveltyLearned,
                 b_task_noveltyNovel) |>
    mutate(task_novelty = case_when(
      .variable == "b_task_noveltyLearned" ~ "Learned",
      .variable == "b_task_noveltyNovel" ~ "Novel"
    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, task_novelty, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = task_novelty,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = task_novelty),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
        y = task_novelty,
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
         title = "Task Novelty",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_cue_selection_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_cue_selectionAssigned,
                 b_cue_selectionSelfMselected) |>
    mutate(cue_selection = case_when(
      .variable == "b_cue_selectionAssigned" ~ "Assigned",
      .variable == "b_cue_selectionSelfMselected" ~ "Self-selected"
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_cue_selectionAssigned,
                 b_cue_selectionSelfMselected) |>
    mutate(cue_selection = case_when(
      .variable == "b_cue_selectionAssigned" ~ "Assigned",
      .variable == "b_cue_selectionSelfMselected" ~ "Self-selected"
    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, cue_selection, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = cue_selection,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = cue_selection),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
        y = cue_selection,
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
         title = "Cue Selection",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_overtness_selection_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_overtness_selectionAssigned,
                 b_overtness_selectionSelfMselected) |>
    mutate(overtness_selection = case_when(
      .variable == "b_overtness_selectionAssigned" ~ "Assigned",
      .variable == "b_overtness_selectionSelfMselected" ~ "Self-selected"
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_overtness_selectionAssigned,
                 b_overtness_selectionSelfMselected) |>
    mutate(overtness_selection = case_when(
      .variable == "b_overtness_selectionAssigned" ~ "Assigned",
      .variable == "b_overtness_selectionSelfMselected" ~ "Self-selected"
    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, overtness_selection, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = overtness_selection,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = overtness_selection),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
        y = overtness_selection,
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
         title = "Overtness Selection",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_training_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_trainingNotraining,
                 b_trainingTraining) |>
    mutate(training = case_when(
      .variable == "b_trainingNotraining" ~ "No training",
      .variable == "b_trainingTraining" ~ "Training"
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_trainingNotraining,
                 b_trainingTraining) |>
    mutate(training = case_when(
      .variable == "b_trainingNotraining" ~ "No training",
      .variable == "b_trainingTraining" ~ "Training"
    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, training, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = training,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = training),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
        y = training,
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
         title = "Training Intervention",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_study_design_model <- function(data, prior_model, model) {

  prior_preds <- prior_model |>
    gather_draws(b_study_designPostMexperimentalDcontrol,
                 b_study_designPreDpostMexperimental,
                 b_study_designPreDpostMexperimentalDcontrol) |>
    mutate(study_design = case_when(
      .variable == "b_study_designPostMexperimentalDcontrol" ~ "Post - experimental/control",
      .variable == "b_study_designPreDpostMexperimental" ~ "Pre/post - experimental",
      .variable == "b_study_designPreDpostMexperimentalDcontrol" ~ "Pre/post - experimental/control"
    ),
    label = "Prior (Hatzigeorgiadis et al., 2011)")

  posterior_preds <- model |>
    gather_draws(b_study_designPostMexperimentalDcontrol,
                 b_study_designPreDpostMexperimental,
                 b_study_designPreDpostMexperimentalDcontrol) |>
    mutate(study_design = case_when(
      .variable == "b_study_designPostMexperimentalDcontrol" ~ "Post - experimental/control",
      .variable == "b_study_designPreDpostMexperimental" ~ "Pre/post - experimental",
      .variable == "b_study_designPreDpostMexperimentalDcontrol" ~ "Pre/post - experimental/control"
    ),
    label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, study_design, label) |>
    mean_qi(.value)

  # Combine prior and posterior samples for plot
  prior_posterior <- rbind(prior_preds, posterior_preds) |>
    mutate(label = factor(label, levels = c("Posterior Pooled Estimate",
                                            "Prior (Hatzigeorgiadis et al., 2011)"
    )))

  posterior_update <- ggplot(data = prior_posterior,
                             aes(x = .value, y = study_design,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = model$data |>
        mutate(label = NA),
      aes(x = yi, y = study_design),
      fill = NA,
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{.value} [{.lower}, {.upper}]"),
        y = study_design,
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
         title = "Study Design",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8))

}

plot_panel_moderators <- function(plot1, plot2, plot3,
                                  plot4, plot5, plot6,
                                  plot7, plot8, plot9) {

  ( (plot1 / plot2 / plot3) |
      (plot4 / plot5 / plot6) |
      (plot7 / plot8 / plot9) ) +
    plot_annotation(tag_levels = "A",
                    title = "Updated Posterior Moderators Estimates",
                    subtitle = "Prior and posterior distributions for pooled estimates, individual effects (ticks), and mean and 95% quantile interval for posterior (text label)") +
    plot_layout(guides = "collect",
                axes = "collect",
                axis_titles = "collect") & theme(legend.position = 'bottom')

}

plot_panel_main_model <- function(plot1, plot2, plot3) {
  # Combine plots
  (plot1 / plot2 / plot3) +
    plot_layout(heights = c(2, 1,1)) +
    plot_annotation(tag_levels = "A")
}

plot_panel_p_hack <- function(plot1, plot2) {
  (plot1 / plot2) +
    plot_annotation(title = "Adjusted estimates assuming possible presence of p-hacking",
                    subtitle = "Estimates from Bayesian random effects main model and mixture model for p-hacking both with, and without, using prior\nNote, models ignore multilevel structure",
                    theme = theme(plot.subtitle = element_text(size = 8))) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

# Supplemental plots (moderator change in evidence and cummulative plot)
plot_BF_curve_motor_demands <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_motor_demandsFine'='Fine',
                               'b_motor_demandsGross'='Gross')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Motor Demands)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_participant_group <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_participant_groupNonMathletes'='Non-athletes',
                               'b_participant_groupBeginnerathletes'='Beginner athletes',
                               'b_participant_groupExperiencedathletes'='Experienced athletes')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Participants)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_selftalk_content <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_selftalk_contentInstructional'='Instructional',
                               'b_selftalk_contentMotivational'='Motivational',
                               'b_selftalk_contentCombinedInstructionalDMotivational'='Combined Instructional/Motivational',
                               'b_selftalk_contentRational'='Rational')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Self-Talk Content)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_matching <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_matchingInstructionalDFine'='Instructional/Fine',
                               'b_matchingInstructionalDGross'='Instructional/Gross',
                               'b_matchingMotivationalDFine'='Motivational/Fine',
                               'b_matchingMotivationalDGross'='Motivational/Gross')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Matching Hypothesis)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_task_novelty <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_task_noveltyLearned'='Learned',
                               'b_task_noveltyNovel'='Novel')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Task Novelty)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_cue_selection <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_cue_selectionAssigned'='Assigned',
                               'b_cue_selectionSelfMselected'='Self-selected')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Cue Selection)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_overtness_selection <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_overtness_selectionAssigned'='Assigned',
                               'b_overtness_selectionSelfMselected'='Self-selected')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Overtness Selection)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_training <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_trainingNotraining'='No-training',
                               'b_trainingTraining'='Training')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Motor Demands)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_study_design <- function(BF_curve) {

  BF_curve |>

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_study_designPostMexperimentalDcontrol'='Post - experimental/control',
                               'b_study_designPreDpostMexperimental'='Pre/post - experimental',
                               'b_study_designPreDpostMexperimentalDcontrol'='Pre/post - experimental/control')) |>

    ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add bands for Jeffreys thresholds for log10BF
    geom_hline(yintercept = c(0,0.5,1,1.5,2), linetype = "dashed") +
    annotate("text", label = stringr::str_wrap("Negative", width = 15),
             x = 1, y=-0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Weak", width = 15),
             x = 1, y=0.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Substantial", width = 15),
             x = 1, y=0.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Strong", width = 15),
             x = 1, y=1.25, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Very Strong", width = 15),
             x = 1, y=1.75, size = 1.75) +
    annotate("text", label = stringr::str_wrap("Decisive", width = 15),
             x = 1, y=2.25, size = 1.75) +
    geom_smooth(se=FALSE, color="black") +
    labs(x = expression(paste("Standardised Mean Difference BF Calculated ", italic("Against"))),
         y = "log10(BF)",
         title = "Change in Evidence (Motor Demands)",
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_cumulative_main_model <- function(data, prior_model, cumulative_draws) {

  prior_draws <- prior_model |>
    spread_draws(b_Intercept) |>
    mutate(year = "Hatzigeorgiadis et al., (2011)",
           prior = NA) |>
    select(year, b_Intercept, prior)

  cumulative_draws <- rbind(prior_draws, cumulative_draws) |>
    mutate(year = factor(year, levels = c(
      "Hatzigeorgiadis et al., (2011)",
      "2011",
      "2012",
      "2013",
      "2014",
      "2015",
      "2016",
      "2017",
      "2018",
      "2019",
      "2020",
      "2021",
      "2022",
      "2023"
    )))

  posterior_summary <- group_by(cumulative_draws, year) |>
    mean_qi(b_Intercept)

  cumulative_draws |>
    ggplot(aes(y = year, x = b_Intercept)) +

    # Add ref line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    # Add densities and point intervals
    geom_density_ridges(
      fill = "darkgrey",
      rel_min_height = 0.01,
      scale = 1,
      alpha = 0.8
    ) +
    stat_pointinterval(point_interval = mean_qi,
                       .width = .95,
                       size = 1) +

    # Add individual study data subset(Data_effects, !is.na(yi)), aes(x = yi, y = study_name)
    geom_point(
      data = subset(data,!is.na(yi)),
      aes(x = yi, y = as.factor(year)),
      position = position_nudge(y = -0.1),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"),
        x = 4
      ),
      hjust = "inward",
      size = 3
    ) +

    # Add labs and theme
    labs(y = "Year",
         x = "Standardised Mean Difference (Positive Values Favour Self-Talk)",
         title = "Cumulative Updating of Posterior Pooled Estimate",
         subtitle = "Each year uses the previous years posterior pooled estimate as its prior\nPosterior distributions, mean and 95% quantile intervals, and individual effects (ticks)") +
    # scale_y_continuous(limits = c(-0.5, 1), breaks = c(-0.5, 0, 0.5, 1)) +
    scale_y_discrete(limits = rev,
                     labels = function(x) str_wrap(x, width = 15)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size=6))
}

# Model checks
make_rhat_plot <- function(model) {
  mod_rhat <- enframe(brms::rhat(model)) |>
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
