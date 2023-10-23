# _targets.R file
library(targets)
library(tarchetypes)
source("R/functions.R")
tar_option_set(
  packages = c(
    "here",
    # "readxl",
    "metafor",
    "orchaRd",
    "brms",
    "modelr",
    "tidybayes",
    "bayesplot",
    "rstan",
    "ggridges",
    "janitor",
    "tidyverse",
    "base",
    # "purrr",
    # "scales",
    # "ggtext",
    # "zoo",
    # "performance",
    # "see",
    # "lme4",
    # "marginaleffects",
    # "broom.mixed",
    "patchwork"
    # "kableExtra",
    # "knitr",
    # "bayestestR",
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

  # Fit, check, and plot main model
  # tar_target(rstan_setup, rstan_setup()), # to run chains in parallel
  tar_target(main_model_prior, set_main_model_prior()),
  tar_target(main_model, fit_main_model(data_effect_sizes, main_model_prior)),
  tar_target(rhat_main_model, make_rhat_plot(main_model)),
  tar_target(trace_plot_main_model, make_trace_plot(main_model)),
  tar_target(pp_check_main_model, make_pp_check(main_model)),
  tar_target(main_model_plot, plot_main_model(data_effect_sizes, main_model)),

  # Fit, check, and plot motor demands model
  tar_target(motor_demands_prior, set_motor_demands_prior()),
  tar_target(motor_demands_model, fit_motor_demands_model(data_effect_sizes, motor_demands_prior)),
  tar_target(rhat_motor_demands_model, make_rhat_plot(motor_demands_model)),
  tar_target(trace_plot_motor_demands_model, make_trace_plot(motor_demands_model)),
  tar_target(pp_check_motor_demands_model, make_pp_check(motor_demands_model)),
  # tar_target(main_model_plot, plot_main_model(data_effect_sizes, main_model)),

  # Make plots tiffs
  tar_target(main_model_plot_tiff, make_plot_tiff(main_model_plot, 7.5, 10, "plots/main_model_plot.tiff"))


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
