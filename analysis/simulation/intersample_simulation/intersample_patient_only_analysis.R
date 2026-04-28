#!/usr/bin/env Rscript

# Simulation Script: Patient-Only Meta Analysis Summary
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Summarize patient-only simulation outputs for PANORAMIC and naive estimators.
# - Compute bias/RMSE/MAE and method win-rate diagnostics across patient counts.
# - Export comparison tables and faceted performance figures.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

source("analysis/simulation/sim_utils.R")

project_root <- sim_find_project_root()
paths <- sim_load_paths(project_root)

in_data_dir <- file.path(paths$data, "interim", "simulation", "intersample_patient_only")
out_tbl_dir <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intersample_patient_only"))
out_fig_dir <- sim_mkdir(file.path(paths$output, "figures", "simulation", "intersample_patient_only"))

in_rds <- file.path(in_data_dir, "intersample_patient_only_estimates_intersample_patient_only.rds")

message("Reading estimates: ", in_rds)
df <- readRDS(in_rds)

need <- c("n_patients", "replicate", "feature_id", "mu_hat", "tau2_patient", "mu_naive", "tau2_patient_naive", "true_mu", "true_tau2_patient")
miss <- setdiff(need, colnames(df))
if (length(miss) > 0L) stop("Input missing columns: ", paste(miss, collapse = ", "))

cmp_long <- dplyr::bind_rows(
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
  )
) %>%
  dplyr::mutate(error = estimate - truth, abs_error = abs(error))

cmp_by_feature <- cmp_long %>%
  dplyr::group_by(n_patients, parameter, method, feature_id) %>%
  dplyr::summarise(
    bias = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mae = mean(abs_error, na.rm = TRUE),
    .groups = "drop"
  )

cmp_summary <- cmp_by_feature %>%
  dplyr::group_by(n_patients, parameter, method) %>%
  dplyr::summarise(
    bias = mean(bias, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
    mae = mean(mae, na.rm = TRUE),
    n_features = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(parameter, n_patients, method)

pair_df <- cmp_long %>%
  dplyr::filter(parameter %in% c("mu", "tau2_patient")) %>%
  dplyr::select(n_patients, replicate, feature_id, parameter, method, abs_error) %>%
  tidyr::pivot_wider(names_from = method, values_from = abs_error, names_sep = "__") %>%
  dplyr::mutate(
    naive_better = abs_error__naive < abs_error__panoramic,
    margin_abs_error = abs_error__panoramic - abs_error__naive
  )

failure_summary <- pair_df %>%
  dplyr::group_by(parameter, feature_id, n_patients) %>%
  dplyr::summarise(
    p_naive_better = mean(naive_better, na.rm = TRUE),
    mean_margin_abs_error = mean(margin_abs_error, na.rm = TRUE),
    median_margin_abs_error = stats::median(margin_abs_error, na.rm = TRUE),
    .groups = "drop"
  )

ab_feature <- "A|B|25|local_comp_enrichment"
if (!ab_feature %in% unique(cmp_long$feature_id)) stop("A|B|25|local_comp_enrichment is missing from estimates.")

ab_long <- cmp_long %>%
  dplyr::filter(feature_id == ab_feature, parameter %in% c("mu", "tau2_patient"))

p_rmse <- ggplot(cmp_by_feature, aes(x = n_patients, y = rmse, color = method)) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.1) +
  facet_grid(parameter ~ feature_id, scales = "free_y") +
  scale_color_manual(values = c(panoramic = "#56B4E9", naive = "#E69F00")) +
  labs(
    title = "Controlled patient-only simulation: RMSE by feature",
    x = "Number of patients",
    y = "RMSE",
    color = "Method"
  ) +
  theme_minimal()

p_bias <- ggplot(cmp_by_feature, aes(x = n_patients, y = bias, color = method)) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  facet_grid(parameter ~ feature_id, scales = "free_y") +
  scale_color_manual(values = c(panoramic = "#56B4E9", naive = "#E69F00")) +
  labs(
    title = "Controlled patient-only simulation: Bias by feature",
    x = "Number of patients",
    y = "Bias",
    color = "Method"
  ) +
  theme_minimal()

p_fail <- ggplot(failure_summary, aes(x = factor(n_patients), y = feature_id, fill = p_naive_better)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ parameter, ncol = 1, scales = "free_y") +
  scale_fill_gradient(low = "#56B4E9", high = "#D55E00", limits = c(0, 1)) +
  labs(
    title = "Failure map: probability naive has lower absolute error",
    x = "Number of patients",
    y = "Feature",
    fill = "P(naive\nbetter)"
  ) +
  theme_minimal()

p_ab <- ggplot(ab_long, aes(x = n_patients, y = abs_error, color = method)) +
  geom_point(alpha = 0.5, size = 1.8, position = position_jitter(width = 0.8, height = 0)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.95, aes(group = method)) +
  stat_summary(fun = mean, geom = "point", size = 2.8) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c(panoramic = "#56B4E9", naive = "#E69F00")) +
  labs(
    title = paste0("A|B controlled diagnostic (", ab_feature, "): absolute error"),
    x = "Number of patients",
    y = "Absolute error",
    color = "Method"
  ) +
  theme_minimal()

summary_csv <- file.path(out_tbl_dir, "intersample_patient_only_naive_comparison_summary.csv")
feature_csv <- file.path(out_tbl_dir, "intersample_patient_only_naive_comparison_by_feature.csv")
failure_csv <- file.path(out_tbl_dir, "intersample_patient_only_failure_summary.csv")

rmse_png <- file.path(out_fig_dir, "intersample_patient_only_naive_comparison_rmse.png")
bias_png <- file.path(out_fig_dir, "intersample_patient_only_naive_comparison_bias.png")
fail_png <- file.path(out_fig_dir, "intersample_patient_only_failure_heatmap.png")
ab_png <- file.path(out_fig_dir, "intersample_patient_only_ab_abs_error.png")

sim_safe_write_csv(cmp_summary, summary_csv)
sim_safe_write_csv(cmp_by_feature, feature_csv)
sim_safe_write_csv(failure_summary, failure_csv)

ggsave(rmse_png, p_rmse, width = 13, height = 8.5, dpi = 300)
ggsave(bias_png, p_bias, width = 13, height = 8.5, dpi = 300)
ggsave(fail_png, p_fail, width = 11.5, height = 8.5, dpi = 300)
ggsave(ab_png, p_ab, width = 9, height = 8, dpi = 300)

message("Saved controlled patient-only analysis outputs:")
message("  summary: ", summary_csv)
message("  by-feature: ", feature_csv)
message("  failure summary: ", failure_csv)
message("  figures: ", rmse_png, ", ", bias_png, ", ", fail_png, ", ", ab_png)
