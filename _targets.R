# _targets.R file
library(targets)
library(tarchetypes)
library(future)
library(future.callr)
plan(callr)
source("R/functions.R")
tar_option_set(
  packages = c(
    "here",
    # "readxl",
    "metafor",
    "brms",
    "modelr",
    "tidybayes",
    "bayesplot",
    "bayestestR",
    "rstan",
    "ggridges",
    "janitor",
    "tidyverse",
    "base",
    "furrr",
    "patchwork",
    # "scales",
    # "ggtext",
    # "zoo",
    # "performance",
    # "see",
    # "lme4",
    "marginaleffects"
    # "broom.mixed",
    # "kableExtra",
    # "knitr",
    # "quarto",
    # "officer",
    # "officedown"
  )
)

list(
  # Load in and prepare data
  tar_target(file, here("data","Final data.csv"), format = "file"),
  tar_target(data, read_prepare_data(file)),
  tar_target(data_effect_sizes, calculate_effect_sizes(data)),

  # Setup rstan to run chains in parallel
  tar_target(rstan, rstan_setup()),

  # Fit, check, and plot main model
  tar_target(main_prior, set_main_prior()),
  tar_target(prior_main_model, sample_prior_main_model(data_effect_sizes, main_prior)),
  tar_target(main_model, fit_main_model(data_effect_sizes, main_prior)),
  tar_target(rhat_main_model, make_rhat_plot(main_model)),
  tar_target(trace_plot_main_model, make_trace_plot(main_model)),
  tar_target(pp_check_main_model, make_pp_check(main_model)),
  tar_target(main_model_logBF_curve, get_logBF_curve(main_model)),
  tar_target(main_model_plot, plot_main_model(data_effect_sizes, prior_main_model, main_model, main_model_logBF_curve)),

  # Fit, check, and plot motor demands model
  tar_target(motor_demands_prior, set_motor_demands_prior()),
  tar_target(prior_motor_demands_model, sample_prior_motor_demands_model(data_effect_sizes, motor_demands_prior)),
  tar_target(motor_demands_model, fit_motor_demands_model(data_effect_sizes, motor_demands_prior)),
  tar_target(rhat_motor_demands_model, make_rhat_plot(motor_demands_model)),
  tar_target(trace_plot_motor_demands_model, make_trace_plot(motor_demands_model)),
  tar_target(pp_check_motor_demands_model, make_pp_check(motor_demands_model)),
  tar_target(motor_demands_model_plot, plot_motor_demands_model(data_effect_sizes, prior_motor_demands_model, motor_demands_model)),
  tar_target(motor_demands_model_logBF_curve, get_logBF_curve(motor_demands_model)),

  # Fit, check, and plot participant group model
  tar_target(participant_group_prior, set_participant_group_prior()),
  tar_target(prior_participant_group_model, sample_prior_participant_group_model(data_effect_sizes, participant_group_prior)),
  tar_target(participant_group_model, fit_participant_group_model(data_effect_sizes, participant_group_prior)),
  tar_target(rhat_participant_group_model, make_rhat_plot(participant_group_model)),
  tar_target(trace_plot_participant_group_model, make_trace_plot(participant_group_model)),
  tar_target(pp_check_participant_group_model, make_pp_check(participant_group_model)),
  tar_target(participant_group_model_plot, plot_participant_group_model(data_effect_sizes, prior_participant_group_model, participant_group_model)),
  tar_target(participant_group_model_logBF_curve, get_logBF_curve(participant_group_model)),

  # Fit, check, and plot self-talk content model
  tar_target(selftalk_content_prior, set_selftalk_content_prior()),
  tar_target(prior_selftalk_content_model, sample_prior_selftalk_content_model(data_effect_sizes, selftalk_content_prior)),
  tar_target(selftalk_content_model, fit_selftalk_content_model(data_effect_sizes, selftalk_content_prior)),
  tar_target(rhat_selftalk_content_model, make_rhat_plot(selftalk_content_model)),
  tar_target(trace_plot_selftalk_content_model, make_trace_plot(selftalk_content_model)),
  tar_target(pp_check_selftalk_content_model, make_pp_check(selftalk_content_model)),
  tar_target(selftalk_content_model_plot, plot_selftalk_content_model(data_effect_sizes, prior_selftalk_content_model, selftalk_content_model)),
  tar_target(selftalk_content_model_logBF_curve, get_logBF_curve(selftalk_content_model)),

  # Fit, check, and plot matching hypothesis model
  tar_target(matching_prior, set_matching_prior()),
  tar_target(prior_matching_model, sample_prior_matching_model(data_effect_sizes, matching_prior)),
  tar_target(matching_model, fit_matching_model(data_effect_sizes, matching_prior)),
  tar_target(rhat_matching_model, make_rhat_plot(matching_model)),
  tar_target(trace_plot_matching_model, make_trace_plot(matching_model)),
  tar_target(pp_check_matching_model, make_pp_check(matching_model)),
  tar_target(matching_model_plot, plot_matching_model(data_effect_sizes, prior_matching_model, matching_model)),
  tar_target(matching_model_logBF_curve, get_logBF_curve(matching_model)),

  # Fit, check, and plot task_novelty model
  tar_target(task_novelty_prior, set_task_novelty_prior()),
  tar_target(prior_task_novelty_model, sample_prior_task_novelty_model(data_effect_sizes, task_novelty_prior)),
  tar_target(task_novelty_model, fit_task_novelty_model(data_effect_sizes, task_novelty_prior)),
  tar_target(rhat_task_novelty_model, make_rhat_plot(task_novelty_model)),
  tar_target(trace_plot_task_novelty_model, make_trace_plot(task_novelty_model)),
  tar_target(pp_check_task_novelty_model, make_pp_check(task_novelty_model)),
  tar_target(task_novelty_model_plot, plot_task_novelty_model(data_effect_sizes, prior_task_novelty_model, task_novelty_model)),
  tar_target(task_novelty_model_logBF_curve, get_logBF_curve(task_novelty_model)),

  # Fit, check, and plot cue_selection model
  tar_target(cue_selection_prior, set_cue_selection_prior()),
  tar_target(prior_cue_selection_model, sample_prior_cue_selection_model(data_effect_sizes, cue_selection_prior)),
  tar_target(cue_selection_model, fit_cue_selection_model(data_effect_sizes, cue_selection_prior)),
  tar_target(rhat_cue_selection_model, make_rhat_plot(cue_selection_model)),
  tar_target(trace_plot_cue_selection_model, make_trace_plot(cue_selection_model)),
  tar_target(pp_check_cue_selection_model, make_pp_check(cue_selection_model)),
  tar_target(cue_selection_model_plot, plot_cue_selection_model(data_effect_sizes, prior_cue_selection_model, cue_selection_model)),
  tar_target(cue_selection_model_logBF_curve, get_logBF_curve(cue_selection_model)),

  # Fit, check, and plot overtness_selection model
  tar_target(overtness_selection_prior, set_overtness_selection_prior()),
  tar_target(prior_overtness_selection_model, sample_prior_overtness_selection_model(data_effect_sizes, overtness_selection_prior)),
  tar_target(overtness_selection_model, fit_overtness_selection_model(data_effect_sizes, overtness_selection_prior)),
  tar_target(rhat_overtness_selection_model, make_rhat_plot(overtness_selection_model)),
  tar_target(trace_plot_overtness_selection_model, make_trace_plot(overtness_selection_model)),
  tar_target(pp_check_overtness_selection_model, make_pp_check(overtness_selection_model)),
  tar_target(overtness_selection_model_plot, plot_overtness_selection_model(data_effect_sizes, prior_overtness_selection_model, overtness_selection_model)),
  tar_target(overtness_selection_model_logBF_curve, get_logBF_curve(overtness_selection_model)),

  # Fit, check, and plot training model
  tar_target(training_prior, set_training_prior()),
  tar_target(prior_training_model, sample_prior_training_model(data_effect_sizes, training_prior)),
  tar_target(training_model, fit_training_model(data_effect_sizes, training_prior)),
  tar_target(rhat_training_model, make_rhat_plot(training_model)),
  tar_target(trace_plot_training_model, make_trace_plot(training_model)),
  tar_target(pp_check_training_model, make_pp_check(training_model)),
  tar_target(training_model_plot, plot_training_model(data_effect_sizes, prior_training_model, training_model)),
  tar_target(training_model_logBF_curve, get_logBF_curve(training_model)),

  # Fit, check, and plot study_design model
  tar_target(study_design_prior, set_study_design_prior()),
  tar_target(prior_study_design_model, sample_prior_study_design_model(data_effect_sizes, study_design_prior)),
  tar_target(study_design_model, fit_study_design_model(data_effect_sizes, study_design_prior)),
  tar_target(rhat_study_design_model, make_rhat_plot(study_design_model)),
  tar_target(trace_plot_study_design_model, make_trace_plot(study_design_model)),
  tar_target(pp_check_study_design_model, make_pp_check(study_design_model)),
  tar_target(study_design_model_plot, plot_study_design_model(data_effect_sizes, prior_study_design_model, study_design_model)),
  tar_target(study_design_model_logBF_curve, get_logBF_curve(study_design_model)),

  # Make panel plot of moderators
  tar_target(moderators_panel_plot, plot_panel_moderators(
    motor_demands_model_plot, participant_group_model_plot, selftalk_content_model_plot,
      matching_model_plot, task_novelty_model_plot, cue_selection_model_plot,
    overtness_selection_model_plot, training_model_plot, study_design_model_plot
  )),

  # Make supplemental plots of evidence change for moderators
  tar_target(BF_curve_motor_demands_plot, plot_BF_curve_motor_demands(motor_demands_model_logBF_curve)),
  tar_target(BF_curve_participant_group_plot, plot_BF_curve_participant_group(participant_group_model_logBF_curve)),
  tar_target(BF_curve_selftalk_content_plot, plot_BF_curve_selftalk_content(selftalk_content_model_logBF_curve)),
  tar_target(BF_curve_matching_plot, plot_BF_curve_matching(matching_model_logBF_curve)),
  tar_target(BF_curve_task_novelty_plot, plot_BF_curve_task_novelty(task_novelty_model_logBF_curve)),
  tar_target(BF_curve_cue_selection_plot, plot_BF_curve_cue_selection(cue_selection_model_logBF_curve)),
  tar_target(BF_curve_overtness_selection_plot, plot_BF_curve_overtness_selection(overtness_selection_model_logBF_curve)),
  tar_target(BF_curve_training_plot, plot_BF_curve_training(training_model_logBF_curve)),
  tar_target(BF_curve_study_design_plot, plot_BF_curve_study_design(study_design_model_logBF_curve)),


  # Make plots tiffs
  tar_target(main_model_plot_tiff, make_plot_tiff(main_model_plot, 7.5, 13.33, "plots/main_model_plot.tiff")),
  tar_target(moderators_panel_plot_tiff, make_plot_tiff(moderators_panel_plot, 21, 9, "plots/moderators_panel_plot.tiff")),
  tar_target(BF_curve_motor_demands_plot_tiff, make_plot_tiff(BF_curve_motor_demands_plot, 7.5, 5, "plots/BF_curve_motor_demands_plot_plot.tiff")),
  tar_target(BF_curve_participant_group_plot_tiff, make_plot_tiff(BF_curve_participant_group_plot, 7.5, 5, "plots/BF_curve_participant_group_plot_plot.tiff")),
  tar_target(BF_curve_selftalk_content_plot_tiff, make_plot_tiff(BF_curve_selftalk_content_plot, 7.5, 5, "plots/BF_curve_selftalk_content_plot_plot.tiff")),
  tar_target(BF_curve_matching_plot_tiff, make_plot_tiff(BF_curve_matching_plot, 7.5, 5, "plots/BF_curve_matching_plot_plot.tiff")),
  tar_target(BF_curve_task_novelty_plot_tiff, make_plot_tiff(BF_curve_task_novelty_plot, 7.5, 5, "plots/BF_curve_task_novelty_plot_plot.tiff")),
  tar_target(BF_curve_cue_selection_plot_tiff, make_plot_tiff(BF_curve_cue_selection_plot, 7.5, 5, "plots/BF_curve_cue_selection_plot_plot.tiff")),
  tar_target(BF_curve_overtness_selection_plot_tiff, make_plot_tiff(BF_curve_overtness_selection_plot, 7.5, 5, "plots/BF_curve_overtness_selection_plot_plot.tiff")),
  tar_target(BF_curve_training_plot_tiff, make_plot_tiff(BF_curve_training_plot, 7.5, 5, "plots/BF_curve_training_plot_plot.tiff")),
  tar_target(BF_curve_study_design_plot_tiff, make_plot_tiff(BF_curve_study_design_plot, 7.5, 5, "plots/BF_curve_study_design_plot_plot.tiff"))


#
#   # Model checks
#   tar_target(model_checks, make_model_checks_tiff(model)),
#
#   # Make and save plots
#   tar_target(individual_data_plot, plot_individual_data(data)),
#   tar_target(individual_data_plot_tiff, make_individual_data_plot_tiff(individual_data_plot)),
#
#   tar_target(model_plot, plot_model(data, model)),
#   tar_target(model_plot_tiff, make_model_plot_tiff(model_plot)),
#
#   tar_target(individual_preds_plot, plot_individual_preds(data, model)),
#   tar_target(individual_preds_plot_tiff, make_individual_preds_plot_tiff(individual_preds_plot)),
#
#   tar_target(main_plot, combine_plots(individual_data_plot, individual_preds_plot, model_plot)),
#
#   # Add baseline adjustment to data
#   tar_target(data_adj, add_adj_data(data)),
#
#   # Fit new model adjusting for each participants baseline
#   tar_target(model_adj, fit_model_adj(data_adj)),
#   tar_target(tidy_model_adj, get_tidy_model(model_adj)),
#
#   # Model checks
#   tar_target(model_adj_checks, make_model_adj_checks_tiff(model_adj)),
#
#   # Make and save plots
#
#   tar_target(model_adj_plot, plot_model_adj(data_adj, model_adj)),
#   tar_target(model_adj_plot_tiff, make_model_adj_plot_tiff(model_adj_plot)),
#
#   tar_target(individual_adj_preds_plot, plot_individual_adj_preds(data_adj, model_adj)),
#   tar_target(individual_adj_preds_plot_tiff, make_individual_adj_preds_plot_tiff(individual_adj_preds_plot)),
#
#   tar_target(main_adj_plot, combine_adj_plots(individual_data_plot, individual_adj_preds_plot, model_adj_plot)),
#
#   # Compare models
#   tar_target(model_comparison_2logBF, compare_models(model, model_adj)),
#
#   # Render the report
#   tar_quarto(report, "report.qmd")
#
#
)
