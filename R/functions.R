# Read in and prepare data for effect size calculations and models
read_prepare_data <- function(file) {
  data <- read_csv(here("data", "Final data.csv")) %>%
    mutate_at(c(2:9, 41:43, 45:48), as.factor) %>%
    clean_names() %>%

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

    )  %>%

    # add calculation of pooled baseline standard deviations
    mutate(pre_sd_pool = sqrt(((n_st - 1) * pre_sd_st ^ 2 +
                                 (n_con - 1) * pre_sd_con ^ 2
    ) /
      (n_st + n_con - 2)),
    post_sd_pool = sqrt(((n_st - 1) * post_sd_st ^ 2 +
                           (n_con - 1) * post_sd_con ^ 2
    ) /
      (n_st + n_con - 2))) %>%

    # Formatting variables for moderator models
    unite(c("selftalk_content", "motor_demands"), col="matching", sep = "/",
          remove = FALSE) %>%


    separate(novelty, into=c("selftalk_novelty1", "task_novelty1"), sep = '\\s*and\\s*',
             remove = FALSE) %>%
    separate(novelty, into=c("selftalk_novelty2", "task_novelty2"), sep = '\\s*but\\s*',
             remove = FALSE) %>%
    mutate(task_novelty = if_else(task_novelty2 == "with the sport", "Learned", "Novel", missing = "Novel")) %>%
    select(-selftalk_novelty1, -selftalk_novelty2, -task_novelty1, -task_novelty2) %>%

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

get_logBF_curve <- function(model) {

  plan(cluster, workers = 10)

  BF_curve <- future_map_dfr(.x = seq(0,1,length=100), .f = function(.x) {
    return(data.frame(effect = .x,
                      BF = bayesfactor_parameters(model, null = .x)))
  })

  plan(sequential)

  return(BF_curve)
}

# Setup rstan to run quicker
rstan_setup <- function() {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)
}

# Priors
set_main_prior <- function() {

  # Taken form Hatsigeorgiadis et al. (2011) overall estimate
  main_model_prior <-
    prior("student_t(3, 0.48, 0.05)", class = "Intercept")
}

sample_prior_main_model <- function(data, prior) {
  prior_main_model <-
    brm(
      yi | se(sqrt(vi)) ~ 1 + (1 | study / experiment / group / effect),
      data = data,
      prior = prior,
      chains = 4,
      cores = 4,
      seed = 1988,
      warmup = 2000,
      iter = 8000,
      control = list(adapt_delta = 0.99),
      sample_prior = "only",
    )
}

set_motor_demands_prior <- function() {

  # Note, this model omits the intercept in order to include the priors on the groups directly
  motor_demands_prior <-
    c(
      prior("student_t(3, 0.67, 0.07)", class = "b", coef = "motor_demandsFine"),
      prior("student_t(3, 0.26, 0.04)", class = "b", coef = "motor_demandsGross")
    )
}

sample_prior_motor_demands_model <- function(data, prior) {
  prior_motor_demands_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + motor_demands + (1 | study / experiment / group / effect),
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
      prior("student_t(3, 0.5, 0.08)", class = "b", coef = "participant_groupNonMathletes"),
      prior("student_t(3, 0.47, 0.08)", class = "b", coef = "participant_groupBeginnerathletes"),
      prior("student_t(3, 0.38, 0.14)", class = "b", coef = "participant_groupExperiencedathletes")
    )
}

sample_prior_participant_group_model <- function(data, prior) {
  prior_participant_group_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + participant_group + (1 | study / experiment / group / effect),
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
      prior("student_t(3, 0.55, 0.08)", class = "b", coef = "selftalk_contentInstructional"),
      prior("student_t(3, 0.37, 0.06)", class = "b", coef = "selftalk_contentMotivational"),
      prior("student_t(3, 0.48, 0.05)", class = "b", coef = "selftalk_contentCombinedInstructionalDMotivational"),
      prior("student_t(3, 0.48, 0.05)", class = "b", coef = "selftalk_contentRational")

    )
}

sample_prior_selftalk_content_model <- function(data, prior) {
  prior_selftalk_content_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + selftalk_content + (1 | study / experiment / group / effect),
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
      prior("student_t(3, 0.83, 0.1)", class = "b", coef = "matchingInstructionalDFine"),
      prior("student_t(3, 0.22, 0.04)", class = "b", coef = "matchingInstructionalDGross"),
      prior("student_t(3, 0.41, 0.09)", class = "b", coef = "matchingMotivationalDFine"),
      prior("student_t(3, 0.33, 0.08)", class = "b", coef = "matchingMotivationalDGross")

    )
}

sample_prior_matching_model <- function(data, prior) {

  data <- data %>%
    filter(matching == "Motivational/Fine" |
             matching == "Motivational/Gross" |
             matching == "Instructional/Fine" |
             matching == "Instructional/Gross")

  prior_matching_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + matching + (1 | study / experiment / group / effect),
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
      prior("student_t(3, 0.41, 0.05)", class = "b", coef = "task_noveltyLearned"),
      prior("student_t(3, 0.73, 0.14)", class = "b", coef = "task_noveltyNovel")

    )
}

sample_prior_task_novelty_model <- function(data, prior) {

  prior_task_novelty_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + task_novelty + (1 | study / experiment / group / effect),
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
      prior("student_t(3, 0.49, 0.07)", class = "b", coef = "cue_selectionAssigned"),
      prior("student_t(3, 0.44, 0.07)", class = "b", coef = "cue_selectionSelfMselected")

    )
}

sample_prior_cue_selection_model <- function(data, prior) {

  prior_cue_selection_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + cue_selection + (1 | study / experiment / group / effect),
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
      prior("student_t(3, 0.49, 0.08)", class = "b", coef = "overtness_selectionAssigned"),
      prior("student_t(3, 0.48, 0.08)", class = "b", coef = "overtness_selectionSelfMselected")

    )
}

sample_prior_overtness_selection_model <- function(data, prior) {

  prior_overtness_selection_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + overtness_selection + (1 | study / experiment / group / effect),
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
      prior("student_t(3, 0.37, 0.04)", class = "b", coef = "trainingNotraining"),
      prior("student_t(3, 0.8, 0.12)", class = "b", coef = "trainingTraining")

    )
}

sample_prior_training_model <- function(data, prior) {

  prior_training_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + training + (1 | study / experiment / group / effect),
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
      prior("student_t(3, 0.36, 0.1)", class = "b", coef = "study_designPreDpostMexperimental"),
      prior("student_t(3, 0.53, 0.06)", class = "b", coef = "study_designPreDpostMexperimentalDcontrol")

    )
}

sample_prior_study_design_model <- function(data, prior) {

  prior_study_design_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + study_design + (1 | study / experiment / group / effect),
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
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
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
      warmup = 4000,
      iter = 40000,
      control = list(adapt_delta = 0.99)
    )
}

fit_participant_group_model <- function(data, prior) {
  participant_group_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + participant_group + (1 | study / experiment / group / effect),
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
      yi | se(sqrt(vi)) ~ 0 + selftalk_content + (1 | study / experiment / group / effect),
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

  data <- data %>%
    filter(matching == "Motivational/Fine" |
             matching == "Motivational/Gross" |
             matching == "Instructional/Fine" |
             matching == "Instructional/Gross")

  matching_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + matching + (1 | study / experiment / group / effect),
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
      yi | se(sqrt(vi)) ~ 0 + task_novelty + (1 | study / experiment / group / effect),
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
      yi | se(sqrt(vi)) ~ 0 + cue_selection + (1 | study / experiment / group / effect),
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

  data <- data %>%
    filter(overtness_selection == "Assigned" |
             overtness_selection == "Self-selected")

  overtness_selection_model <-
    brm(
      yi | se(sqrt(vi)) ~ 0 + overtness_selection + (1 | study / experiment / group / effect),
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
      yi | se(sqrt(vi)) ~ 0 + training + (1 | study / experiment / group / effect),
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
      yi | se(sqrt(vi)) ~ 0 + study_design + (1 | study / experiment / group / effect),
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

  data <- data %>%
    mutate(year = as.numeric(as.character(year)))

  # Create empty data frame for draws
  cumulative_draws <- data.frame(year = as.numeric(),
                                  b_Intercept = as.numeric,
                                 prior = as.character())

  # Loop through each model using previous years posterior as prior
  for (i in unique(data$year)){

    data_year <- data %>%
      filter(year == i)

    posterior <- c(0.48, 0.05, 3)

    prior <- set_prior(paste("student_t(",posterior[3],",", posterior[1],",", posterior[2],")"),
                       class = "Intercept")

    main_model <-
      brm(
        yi | se(sqrt(vi)) ~ 1 + (1 | study / experiment / group / effect),
        data = data_year,
        prior = prior,
        chains = 4,
        cores = 4,
        seed = 1988,
        warmup = 2000,
        iter = 6000,
        control = list(adapt_delta = 0.99)
      )

    draws <- main_model %>%
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

# Plots
plot_main_model <- function(data, prior_model, model, BF_curve) {

  # Sample draws from the prior distribution
  prior_draws <-
    prior_model %>%
    spread_draws(b_Intercept) %>%
    mutate(study = "Prior (Hatzigeorgiadis et al., 2011)",
           label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Sample draws from the individual studies
  study_draws <- model %>%
    spread_draws(b_Intercept, r_study[study, ]) %>%
    mutate(b_Intercept = b_Intercept + r_study,
           study = as.factor(study)) %>%
    select(.chain, .iteration, .draw, b_Intercept, study) %>%
    merge(data[, c(1,3)], by = "study")

  study_summary <- group_by(study_draws, label) %>%
    mean_qi(b_Intercept) %>%
    mutate(b_Intercept = b_Intercept,
           .lower = .lower,
           .upper = .upper)

  # Sample draws from the posterior distribution
  posterior_draws <- model %>%
    spread_draws(b_Intercept) %>%
    mutate(study = "Posterior Pooled Estimate",
           label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_draws, label) %>%
    mean_qi(b_Intercept) %>%
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

  pred_int_data <- median_qi(pred_int_data) %>%
    mutate(label = as.character("Posterior Pooled Estimate"))

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
  prior_posterior <- rbind(posterior_draws[, c(4, 6)], prior_draws[, c(4, 6)]) %>%
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

  # Bayes Factor curve plot
  BF_curve_plot <- BF_curve %>% ggplot(aes(x=effect, y=log10(exp(BF.log_BF)))) +

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          axis.title.y = element_text(vjust = -40),
          plot.subtitle = element_text(size = 6))

  # Combine plots
  (forest_study / posterior_update / BF_curve_plot) +
    plot_layout(heights = c(2, 1,1)) +
    plot_annotation(tag_levels = "A")

}

plot_motor_demands_model <- function(data, prior_model, model) {

  # Prior distribution samples
  nd <- datagrid(model = model,
                 motor_demands = c("Fine", "Gross")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
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
      data = data %>%
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
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
}

plot_participant_group_model <- function(data, prior_model, model) {

  # Prior distribution samples
  nd <- datagrid(model = model,
                 participant_group = c("Non-athletes",
                                       "Experienced athletes",
                                       "Beginner athletes")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
                                 newdata = nd,
                                 re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, participant_group, label) %>%
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
                             aes(x = draw, y = participant_group,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = data %>%
        filter(!is.na(yi)) %>%
        mutate(label = NA),
      aes(x = yi, y = participant_group),
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{draw} [{.lower}, {.upper}]"),
        y = participant_group,
        x = 1.5
      ),
      hjust = "inward",
      size = 3,
      color = "black",
      position = position_nudge(y=0.1)
    ) +

    scale_y_discrete(limits = c("Non-athletes",
                                "Beginner athletes",
                                "Experienced athletes")) +
    scale_color_manual(values = c("#009E73", "#E69F00")) +
    scale_fill_manual(values = c("#009E73", "#E69F00")) +

    labs(x = "Standardised Mean Difference", # summary measure
         y = element_blank(),
         fill = "",
         title = "Participants (Non-athletes versus Beginner Athletes versus Experienced Athletes)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
}

plot_selftalk_content_model <- function(data, prior_model, model) {

  # Prior distribution samples
  nd <- datagrid(model = model,
                 selftalk_content = c("Instructional",
                                      "Motivational",
                                      "Combined Instructional/Motivational",
                                      "Rational")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
                                 newdata = nd,
                                 re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, selftalk_content, label) %>%
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
                             aes(x = draw, y = selftalk_content,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = data %>%
        filter(!is.na(yi)) %>%
        mutate(label = NA),
      aes(x = yi, y = selftalk_content),
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{draw} [{.lower}, {.upper}]"),
        y = selftalk_content,
        x = 1.5
      ),
      hjust = "inward",
      size = 3,
      color = "black",
      position = position_nudge(y=0.1)
    ) +

    scale_y_discrete(limits = c("Rational",
                                "Combined Instructional/Motivational",
                                "Instructional",
                                "Motivational"
                                )) +
    scale_color_manual(values = c("#009E73", "#E69F00")) +
    scale_fill_manual(values = c("#009E73", "#E69F00")) +

    labs(x = "Standardised Mean Difference", # summary measure
         y = element_blank(),
         fill = "",
         title = "Self-Talk Content (Instructional versus Motivational versus Combined versus Rational)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
}

plot_matching_model <- function(data, prior_model, model) {

  # Filter for effects in model
  data <- data %>%
    filter(matching == "Motivational/Fine" |
             matching == "Motivational/Gross" |
             matching == "Instructional/Fine" |
             matching == "Instructional/Gross")

  # Prior distribution samples
  nd <- datagrid(model = model,
                 matching = c("Instructional/Fine",
                              "Instructional/Gross",
                              "Motivational/Fine",
                              "Motivational/Gross")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
                                 newdata = nd,
                                 re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, matching, label) %>%
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
                             aes(x = draw, y = matching,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = data %>%
        filter(!is.na(yi)) %>%
        mutate(label = NA),
      aes(x = yi, y = matching),
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{draw} [{.lower}, {.upper}]"),
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
         title = "Matching Hypothesis (Instructional/Motivational and Gross/Fine)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
}

plot_task_novelty_model <- function(data, prior_model, model) {

  # Prior distribution samples
  nd <- datagrid(model = model,
                 task_novelty = c("Learned",
                                  "Novel")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
                                 newdata = nd,
                                 re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, task_novelty, label) %>%
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
                             aes(x = draw, y = task_novelty,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = data %>%
        filter(!is.na(yi)) %>%
        mutate(label = NA),
      aes(x = yi, y = task_novelty),
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{draw} [{.lower}, {.upper}]"),
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
         title = "Task Novelty (Novel versus Learned)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
}

plot_cue_selection_model <- function(data, prior_model, model) {

  # Prior distribution samples
  nd <- datagrid(model = model,
                 cue_selection = c("Assigned",
                                  "Self-selected")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
                                 newdata = nd,
                                 re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, cue_selection, label) %>%
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
                             aes(x = draw, y = cue_selection,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = data %>%
        filter(!is.na(yi)) %>%
        mutate(label = NA),
      aes(x = yi, y = cue_selection),
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{draw} [{.lower}, {.upper}]"),
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
         title = "Cue Selection (Assigned versus Self-selected)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
}

plot_overtness_selection_model <- function(data, prior_model, model) {

  data <- data %>%
    filter(overtness_selection == "Assigned" |
             overtness_selection == "Self-selected")

  # Prior distribution samples
  nd <- datagrid(model = model,
                 overtness_selection = c("Assigned",
                                   "Self-selected")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
                                 newdata = nd,
                                 re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, overtness_selection, label) %>%
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
                             aes(x = draw, y = overtness_selection,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = data %>%
        filter(!is.na(yi)) %>%
        mutate(label = NA),
      aes(x = yi, y = overtness_selection),
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{draw} [{.lower}, {.upper}]"),
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
         title = "Overtness Selection (Assigned versus Self-selected)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
}

plot_training_model <- function(data, prior_model, model) {

  # Prior distribution samples
  nd <- datagrid(model = model,
                 training = c("Training",
                              "No training")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
                                 newdata = nd,
                                 re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Posterior Pooled Estimate")

  posterior_summary <- group_by(posterior_preds, training, label) %>%
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
                             aes(x = draw, y = training,
                                 color = label, fill = label)
  ) +
    # Add reference line at zero
    geom_vline(xintercept = 0, linetype = 2) +

    stat_slab(alpha = 0.6, linewidth = 0) +

    # Add individual study data
    geom_point(
      data = data %>%
        filter(!is.na(yi)) %>%
        mutate(label = NA),
      aes(x = yi, y = training),
      position = position_nudge(y = -0.05),
      shape = "|"
    ) +

    # Add text and labels
    geom_text(
      data = mutate_if(posterior_summary,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{draw} [{.lower}, {.upper}]"),
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
         title = "Training Intervention? (No-training versus Training)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
}

plot_study_design_model <- function(data, prior_model, model) {

  # Prior distribution samples
  nd <- datagrid(model = model,
                 study_design = c("Pre/post - experimental/control", "Post - experimental/control", "Pre/post - experimental")
  )

  prior_preds <- predictions(prior_model, type = "response",
                             newdata = nd,
                             re_formula = NA) %>%
    posterior_draws() %>%
    transform(type = "Response") %>%
    mutate(label = "Prior (Hatzigeorgiadis et al., 2011)")

  # Posterior distribution samples
  posterior_preds <- predictions(model, type = "response",
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
      data = data %>%
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

    scale_color_manual(values = c("#009E73", "#E69F00")) +
    scale_fill_manual(values = c("#009E73", "#E69F00")) +

    labs(x = "Standardised Mean Difference", # summary measure
         y = element_blank(),
         fill = "",
         title = "Study Design (Pre/Post or Post Only With Experimental Only or Control)",
    ) +
    scale_x_continuous(limits = c(-1, 1.5), breaks = c(-1,-0.5, 0, 0.5, 1)) +
    guides(color = "none",
           shape = "none") +
    theme_classic() +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA),
          title = element_text(size=8),
          plot.subtitle = element_text(size=6))
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
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')

}

# Supplemental plots (moderator change in evidence and cummulative plot)
plot_BF_curve_motor_demands <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_motor_demandsFine'='Fine',
                               'b_motor_demandsGross'='Gross')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_participant_group <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_participant_groupNonMathletes'='Non-athletes',
                               'b_participant_groupBeginnerathletes'='Beginner athletes',
                               'b_participant_groupExperiencedathletes'='Experienced athletes')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_selftalk_content <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_selftalk_contentInstructional'='Instructional',
                               'b_selftalk_contentMotivational'='Motivational',
                               'b_selftalk_contentCombinedInstructionalDMotivational'='Combined Instructional/Motivational',
                               'b_selftalk_contentRational'='Rational')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_matching <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_matchingInstructionalDFine'='Instructional/Fine',
                               'b_matchingInstructionalDGross'='Instructional/Gross',
                               'b_matchingMotivationalDFine'='Motivational/Fine',
                               'b_matchingMotivationalDGross'='Motivational/Gross')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_task_novelty <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_task_noveltyLearned'='Learned',
                               'b_task_noveltyNovel'='Novel')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_cue_selection <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_cue_selectionAssigned'='Assigned',
                               'b_cue_selectionSelfMselected'='Self-selected')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_overtness_selection <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_overtness_selectionAssigned'='Assigned',
                               'b_overtness_selectionSelfMselected'='Self-selected')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_training <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_trainingNotraining'='No-training',
                               'b_trainingTraining'='Training')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_BF_curve_study_design <- function(BF_curve) {

  BF_curve %>%

    # Recode levels
    mutate(BF.Parameter=recode(BF.Parameter,
                               'b_study_designPostMexperimentalDcontrol'='Post - experimental/control',
                               'b_study_designPreDpostMexperimental'='Pre/post - experimental',
                               'b_study_designPreDpostMexperimentalDcontrol'='Pre/post - experimental/control')) %>%

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
         subtitle = "Positive values indicate greater evidence against standardised mean difference after updating prior\nNote, thresholds for evidence are indicated for Jeffreys (1961) scale"
    ) +
    scale_x_continuous(limits = c(-0.55, 1.1),
                       breaks = c(-0.5, 0, 0.5, 1)) +
    facet_wrap(~BF.Parameter) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          plot.subtitle = element_text(size = 6))

}

plot_cumulative_main_model <- function(data, prior_model, cumulative_draws) {

  prior_draws <- prior_model %>%
    spread_draws(b_Intercept) %>%
    mutate(year = "Hatzigeorgiadis et al., (2011)",
           prior = NA) %>%
    select(year, b_Intercept, prior)

  cumulative_draws <- rbind(prior_draws, cumulative_draws) %>%
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
      "2022"
    )))

  posterior_summary <- group_by(cumulative_draws, year) %>%
    mean_qi(b_Intercept)

  cumulative_draws %>%
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
