#!/usr/bin/env Rscript

# Simulation Script: Intersample Meta Analysis Summary
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Summarize intersample simulation estimates against truth across patient-count settings.
# - Compare PANORAMIC against naive baselines for mean and heterogeneity parameters.
# - Export performance tables and comparison figures for manuscript reporting.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

source("analysis/simulation/sim_utils.R")

project_root <- sim_find_project_root()
paths <- sim_load_paths(project_root)

in_data_dir <- file.path(paths$data, "interim", "simulation", "intersample")
out_tbl_dir <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intersample"))
out_fig_dir <- sim_mkdir(file.path(paths$output, "figures", "simulation", "intersample"))

in_rds <- file.path(in_data_dir, "intersample_meta_estimates.rds")
if (!file.exists(in_rds)) stop("Missing input file: ", in_rds)

df <- readRDS(in_rds)
need <- c(
  "n_patients", "replicate", "feature_id",
  "mu_hat", "se_mu", "tau2_patient", "tau2_sample",
  "true_mu", "true_tau2_patient", "true_tau2_sample",
  "mu_naive", "tau2_patient_naive", "tau2_sample_naive"
)
miss <- setdiff(need, colnames(df))
if (length(miss) > 0L) stop("Input is missing required columns: ", paste(miss, collapse = ", "))

df <- df %>%
  mutate(
    mu_error = mu_hat - true_mu,
    tau2_patient_error = tau2_patient - true_tau2_patient,
    tau2_sample_error = tau2_sample - true_tau2_sample,
    mu_covered_95 = is.finite(se_mu) & abs(mu_hat - true_mu) <= 1.96 * se_mu
  )

by_feature <- df %>%
  group_by(n_patients, feature_id) %>%
  summarise(
    mu_bias = mean(mu_error, na.rm = TRUE),
    mu_rmse = sqrt(mean(mu_error^2, na.rm = TRUE)),
    mu_coverage_95 = mean(mu_covered_95, na.rm = TRUE),
    tau2_patient_bias = mean(tau2_patient_error, na.rm = TRUE),
    tau2_patient_rmse = sqrt(mean(tau2_patient_error^2, na.rm = TRUE)),
    tau2_sample_bias = mean(tau2_sample_error, na.rm = TRUE),
    tau2_sample_rmse = sqrt(mean(tau2_sample_error^2, na.rm = TRUE)),
    n_reps = n(),
    .groups = "drop"
  )

summary_tbl <- by_feature %>%
  group_by(n_patients) %>%
  summarise(
    mu_bias = mean(mu_bias, na.rm = TRUE),
    mu_rmse = mean(mu_rmse, na.rm = TRUE),
    mu_coverage_95 = mean(mu_coverage_95, na.rm = TRUE),
    tau2_patient_bias = mean(tau2_patient_bias, na.rm = TRUE),
    tau2_patient_rmse = mean(tau2_patient_rmse, na.rm = TRUE),
    tau2_sample_bias = mean(tau2_sample_bias, na.rm = TRUE),
    tau2_sample_rmse = mean(tau2_sample_rmse, na.rm = TRUE),
    n_features = n(),
    .groups = "drop"
  )

summary_tbl <- summary_tbl %>% arrange(n_patients)

summary_csv <- file.path(out_tbl_dir, "intersample_meta_summary.csv")
feature_csv <- file.path(out_tbl_dir, "intersample_meta_by_feature.csv")

sim_safe_write_csv(summary_tbl, summary_csv)
sim_safe_write_csv(by_feature, feature_csv)

# Per user request, no standalone RMSE/coverage/bias-vs-n_patients figures are
# generated here. Figure outputs are limited to PANORAMIC-vs-naive comparisons.

# PANORAMIC vs naive comparison for mean and heterogeneity parameters
cmp_rows <- list(
  data.frame(
    n_patients = df$n_patients,
    replicate = df$replicate,
    feature_id = df$feature_id,
    parameter = "mu",
    method = "panoramic",
    estimate = df$mu_hat,
    truth = df$true_mu,
    stringsAsFactors = FALSE
  ),
  data.frame(
    n_patients = df$n_patients,
    replicate = df$replicate,
    feature_id = df$feature_id,
    parameter = "mu",
    method = "naive",
    estimate = df$mu_naive,
    truth = df$true_mu,
    stringsAsFactors = FALSE
  ),
  data.frame(
    n_patients = df$n_patients,
    replicate = df$replicate,
    feature_id = df$feature_id,
    parameter = "tau2_patient",
    method = "panoramic",
    estimate = df$tau2_patient,
    truth = df$true_tau2_patient,
    stringsAsFactors = FALSE
  ),
  data.frame(
    n_patients = df$n_patients,
    replicate = df$replicate,
    feature_id = df$feature_id,
    parameter = "tau2_patient",
    method = "naive",
    estimate = df$tau2_patient_naive,
    truth = df$true_tau2_patient,
    stringsAsFactors = FALSE
  ),
  data.frame(
    n_patients = df$n_patients,
    replicate = df$replicate,
    feature_id = df$feature_id,
    parameter = "tau2_sample",
    method = "panoramic",
    estimate = df$tau2_sample,
    truth = df$true_tau2_sample,
    stringsAsFactors = FALSE
  ),
  data.frame(
    n_patients = df$n_patients,
    replicate = df$replicate,
    feature_id = df$feature_id,
    parameter = "tau2_sample",
    method = "naive",
    estimate = df$tau2_sample_naive,
    truth = df$true_tau2_sample,
    stringsAsFactors = FALSE
  )
)

cmp_long <- dplyr::bind_rows(cmp_rows) %>%
  dplyr::mutate(error = estimate - truth)

  cmp_by_feature <- cmp_long %>%
    dplyr::group_by(n_patients, parameter, method, feature_id) %>%
    dplyr::summarise(
      bias = mean(error, na.rm = TRUE),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      .groups = "drop"
    )

  cmp_summary <- cmp_by_feature %>%
    dplyr::group_by(n_patients, parameter, method) %>%
    dplyr::summarise(
      bias = mean(bias, na.rm = TRUE),
      rmse = mean(rmse, na.rm = TRUE),
      n_features = dplyr::n(),
      .groups = "drop"
    )

  cmp_summary <- cmp_summary %>%
    dplyr::arrange(parameter, n_patients, method)

  cmp_csv <- file.path(out_tbl_dir, "intersample_naive_comparison_summary.csv")
  cmp_feature_csv <- file.path(out_tbl_dir, "intersample_naive_comparison_by_feature.csv")
  sim_safe_write_csv(cmp_summary, cmp_csv)
  sim_safe_write_csv(cmp_by_feature, cmp_feature_csv)

  p_cmp_rmse <- ggplot(cmp_by_feature, aes(x = n_patients, y = rmse, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_grid(parameter ~ feature_id, scales = "free_y") +
    labs(
      title = "PANORAMIC vs naive: RMSE by feature",
      x = "Number of patients",
      y = "RMSE",
      color = "Method"
    ) +
    scale_color_manual(values = c(panoramic = "#56B4E9", naive = "#E69F00")) +
    theme_minimal()

  p_cmp_bias <- ggplot(cmp_by_feature, aes(x = n_patients, y = bias, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_grid(parameter ~ feature_id, scales = "free_y") +
    labs(
      title = "PANORAMIC vs naive: bias by feature",
      x = "Number of patients",
      y = "Bias",
      color = "Method"
    ) +
    scale_color_manual(values = c(panoramic = "#56B4E9", naive = "#E69F00")) +
    theme_minimal()

  cmp_rmse_png <- file.path(out_fig_dir, "intersample_naive_comparison_rmse.png")
  cmp_rmse_pdf <- file.path(out_fig_dir, "intersample_naive_comparison_rmse.pdf")
  cmp_bias_png <- file.path(out_fig_dir, "intersample_naive_comparison_bias.png")
  cmp_bias_pdf <- file.path(out_fig_dir, "intersample_naive_comparison_bias.pdf")

  ggsave(cmp_rmse_png, p_cmp_rmse, width = 13.0, height = 8.5, dpi = 300)
  ggsave(cmp_rmse_pdf, p_cmp_rmse, width = 13.0, height = 8.5)
  ggsave(cmp_bias_png, p_cmp_bias, width = 13.0, height = 8.5, dpi = 300)
  ggsave(cmp_bias_pdf, p_cmp_bias, width = 13.0, height = 8.5)

  # Failure-regime diagnostics: explicitly report where naive beats PANORAMIC.
  cmp_pair <- cmp_long %>%
    dplyr::filter(parameter %in% c("mu", "tau2_patient")) %>%
    dplyr::mutate(abs_error = abs(error)) %>%
    dplyr::select(n_patients, replicate, feature_id, parameter, method, error, abs_error) %>%
    tidyr::pivot_wider(
      names_from = method,
      values_from = c(error, abs_error),
      names_sep = "__"
    ) %>%
    dplyr::mutate(
      naive_better = abs_error__naive < abs_error__panoramic,
      panoramic_better = abs_error__panoramic < abs_error__naive,
      tied = abs_error__naive == abs_error__panoramic,
      abs_error_margin = abs_error__panoramic - abs_error__naive
    )

  fail_summary <- cmp_pair %>%
    dplyr::group_by(parameter, feature_id, n_patients) %>%
    dplyr::summarise(
      n_reps = dplyr::n(),
      p_naive_better = mean(naive_better, na.rm = TRUE),
      p_panoramic_better = mean(panoramic_better, na.rm = TRUE),
      p_tied = mean(tied, na.rm = TRUE),
      mean_abs_error_margin = mean(abs_error_margin, na.rm = TRUE),
      median_abs_error_margin = stats::median(abs_error_margin, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(parameter, feature_id, n_patients)

  fail_csv <- file.path(out_tbl_dir, "intersample_failure_regimes_by_feature.csv")
  sim_safe_write_csv(fail_summary, fail_csv)

  p_fail_heatmap <- ggplot(
    fail_summary,
    aes(x = factor(n_patients), y = feature_id, fill = p_naive_better)
  ) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_wrap(~ parameter, ncol = 1, scales = "free_y") +
    scale_fill_gradient(low = "#56B4E9", high = "#D55E00", limits = c(0, 1)) +
    labs(
      title = "Failure regimes: probability naive has lower absolute error",
      x = "Number of patients",
      y = "Feature",
      fill = "P(naive\nbetter)"
    ) +
    theme_minimal()

  fail_heatmap_png <- file.path(out_fig_dir, "intersample_failure_regimes_heatmap.png")
  fail_heatmap_pdf <- file.path(out_fig_dir, "intersample_failure_regimes_heatmap.pdf")
  ggsave(fail_heatmap_png, p_fail_heatmap, width = 11.5, height = 8.5, dpi = 300)
  ggsave(fail_heatmap_pdf, p_fail_heatmap, width = 11.5, height = 8.5)

  # Dedicated A|B diagnostics.
  ab_candidates <- unique(df$feature_id[grepl("^A\\|B\\|", df$feature_id) &
    grepl("\\|local_comp_enrichment$", df$feature_id)])
  ab_feature_id <- if ("A|B|25|local_comp_enrichment" %in% ab_candidates) {
    "A|B|25|local_comp_enrichment"
  } else if (length(ab_candidates) > 0L) {
    ab_candidates[[1]]
  } else {
    NA_character_
  }

  if (!is.na(ab_feature_id)) {
    ab_cmp <- cmp_pair %>%
      dplyr::filter(feature_id == ab_feature_id) %>%
      dplyr::arrange(parameter, n_patients, replicate)

    ab_summary <- ab_cmp %>%
      dplyr::group_by(parameter, n_patients) %>%
      dplyr::summarise(
        n_reps = dplyr::n(),
        p_naive_better = mean(naive_better, na.rm = TRUE),
        mean_abs_error_margin = mean(abs_error_margin, na.rm = TRUE),
        median_abs_error_margin = stats::median(abs_error_margin, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(parameter, n_patients)

    ab_cmp_csv <- file.path(out_tbl_dir, "intersample_ab_failure_cases.csv")
    ab_summary_csv <- file.path(out_tbl_dir, "intersample_ab_failure_summary.csv")
    sim_safe_write_csv(ab_cmp, ab_cmp_csv)
    sim_safe_write_csv(ab_summary, ab_summary_csv)

    ab_long <- cmp_long %>%
      dplyr::filter(feature_id == ab_feature_id, parameter %in% c("mu", "tau2_patient")) %>%
      dplyr::mutate(abs_error = abs(error))

    p_ab_abs_error <- ggplot(ab_long, aes(x = n_patients, y = abs_error, color = method)) +
      geom_point(
        alpha = 0.55,
        size = 2.1,
        position = position_jitter(width = 0.8, height = 0)
      ) +
      stat_summary(fun = mean, geom = "line", linewidth = 0.95, aes(group = method)) +
      stat_summary(fun = mean, geom = "point", size = 2.8) +
      facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
      labs(
        title = paste0("A|B diagnostics (", ab_feature_id, "): absolute error by method"),
        x = "Number of patients",
        y = "Absolute error",
        color = "Method"
      ) +
      scale_color_manual(values = c(panoramic = "#56B4E9", naive = "#E69F00")) +
      theme_minimal()

    ab_abs_png <- file.path(out_fig_dir, "intersample_ab_abs_error.png")
    ab_abs_pdf <- file.path(out_fig_dir, "intersample_ab_abs_error.pdf")
    ggsave(ab_abs_png, p_ab_abs_error, width = 10.0, height = 8.0, dpi = 300)
    ggsave(ab_abs_pdf, p_ab_abs_error, width = 10.0, height = 8.0)

    if (all(c("mean_vi", "mean_n_cells_total", "n_samples_obs") %in% colnames(df))) {
      ab_precision <- df %>%
        dplyr::filter(feature_id == ab_feature_id) %>%
        dplyr::select(
          n_patients, replicate, mean_vi, mean_n_cells_total, n_samples_obs,
          mu_hat, mu_naive, true_mu, tau2_patient, tau2_patient_naive, true_tau2_patient
        ) %>%
        dplyr::mutate(
          vi_quartile = dplyr::ntile(mean_vi, 4L),
          size_quartile = dplyr::ntile(mean_n_cells_total, 4L)
        )

      ab_precision_long <- dplyr::bind_rows(
        ab_precision %>%
          dplyr::transmute(
            n_patients, replicate, vi_quartile, size_quartile, n_samples_obs,
            parameter = "mu", method = "panoramic",
            abs_error = abs(mu_hat - true_mu)
          ),
        ab_precision %>%
          dplyr::transmute(
            n_patients, replicate, vi_quartile, size_quartile, n_samples_obs,
            parameter = "mu", method = "naive",
            abs_error = abs(mu_naive - true_mu)
          ),
        ab_precision %>%
          dplyr::transmute(
            n_patients, replicate, vi_quartile, size_quartile, n_samples_obs,
            parameter = "tau2_patient", method = "panoramic",
            abs_error = abs(tau2_patient - true_tau2_patient)
          ),
        ab_precision %>%
          dplyr::transmute(
            n_patients, replicate, vi_quartile, size_quartile, n_samples_obs,
            parameter = "tau2_patient", method = "naive",
            abs_error = abs(tau2_patient_naive - true_tau2_patient)
          )
      )

      ab_precision_summary <- ab_precision_long %>%
        dplyr::group_by(parameter, vi_quartile, method) %>%
        dplyr::summarise(
          mean_abs_error = mean(abs_error, na.rm = TRUE),
          median_abs_error = stats::median(abs_error, na.rm = TRUE),
          n = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::arrange(parameter, vi_quartile, method)

      ab_precision_csv <- file.path(out_tbl_dir, "intersample_ab_precision_strata.csv")
      sim_safe_write_csv(ab_precision_summary, ab_precision_csv)
    }
  }
}

message("Saved intersample analysis outputs:")
message("  summary: ", summary_csv)
message("  by-feature: ", feature_csv)
if (has_naive) {
  message("  naive comparison: ", file.path(out_tbl_dir, "intersample_naive_comparison_summary.csv"))
  message("  naive comparison figures: ",
          file.path(out_fig_dir, "intersample_naive_comparison_rmse.png"), ", ",
          file.path(out_fig_dir, "intersample_naive_comparison_bias.png"))
  message("  failure regime table: ", file.path(out_tbl_dir, "intersample_failure_regimes_by_feature.csv"))
  message("  failure regime heatmap: ", file.path(out_fig_dir, "intersample_failure_regimes_heatmap.png"))
  ab_diag <- file.path(out_tbl_dir, "intersample_ab_failure_summary.csv")
  if (file.exists(ab_diag)) message("  A|B diagnostics: ", ab_diag)
}
