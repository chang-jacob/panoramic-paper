#!/usr/bin/env Rscript

# Simulation Script: Intrasample Variance Analysis
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Analyze intrasample simulation outputs for variance calibration and coverage.
# - Build feature-level and replicate-level calibration summaries by pattern.
# - Export manuscript-ready summary tables and figures.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source("analysis/simulation/sim_utils.R")

project_root <- sim_find_project_root()
paths <- sim_load_paths(project_root)

in_data_dir <- file.path(paths$data, "interim", "simulation", "intrasample")
out_tbl_root <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intrasample"))
out_fig_root <- sim_mkdir(file.path(paths$output, "figures", "simulation", "intrasample"))

cfg <- list(pattern_arg = "all")

pattern_env <- Sys.getenv("PANORAMIC_INTRASAMPLE_PATTERN", unset = "")
if (nzchar(pattern_env)) cfg$pattern_arg <- pattern_env

all_patterns <- c("uniform", "opposing_gradient", "clustered")
pattern_arg <- tolower(trimws(as.character(cfg$pattern_arg)))
if (identical(pattern_arg, "all")) {
  patterns <- all_patterns
} else {
  patterns <- unique(trimws(strsplit(pattern_arg, ",", fixed = TRUE)[[1]]))
  patterns <- patterns[patterns %in% all_patterns]
}
if (length(patterns) == 0L) {
  stop(
    "Invalid pattern_arg: ", pattern_arg,
    ". Use one of: all, ", paste(all_patterns, collapse = ", ")
  )
}

pattern_tag <- if (length(patterns) == length(all_patterns)) "combined" else paste(patterns, collapse = "_")
out_tbl_dir <- sim_mkdir(file.path(out_tbl_root, pattern_tag))
out_fig_dir <- sim_mkdir(file.path(out_fig_root, pattern_tag))

load_pattern_input <- function(in_data_dir, patterns, all_patterns) {
  canonical_rds <- file.path(in_data_dir, "intrasample_spatialstats_long.rds")
  df <- readRDS(canonical_rds)
  if (length(patterns) < length(all_patterns)) {
    df <- df %>% filter(pattern %in% patterns)
  }
  df
}

df_long <- load_pattern_input(in_data_dir, patterns, all_patterns)

needed <- c("pattern", "replicate", "feature_id", "ct1", "ct2", "radius_um", "stat", "yi", "vi")
miss <- setdiff(needed, colnames(df_long))
if (length(miss) > 0L) stop("Input is missing required columns: ", paste(miss, collapse = ", "))

if (nrow(df_long) == 0L) {
  stop("No rows available for selected pattern(s): ", paste(patterns, collapse = ", "))
}

truth_df <- sim_build_feature_truth(df_long)
joined <- sim_join_truth(df_long, truth_df)

metrics <- bind_rows(lapply(split(joined, joined$pattern), function(x) {
  out <- sim_variance_calibration_metrics(x)
  out$pattern <- unique(x$pattern)[1]
  out
}))

rep_corr <- sim_per_replicate_variance_correlation(df_long, truth_df)

feature_scatter <- truth_df %>%
  filter(is.finite(empirical_variance), is.finite(mean_bootstrap_variance), empirical_variance >= 0, mean_bootstrap_variance >= 0) %>%
  mutate(
    log_empirical_variance = log1p(empirical_variance),
    log_bootstrap_variance = log1p(mean_bootstrap_variance)
  )

metrics_csv <- file.path(out_tbl_dir, paste0("intrasample_calibration_metrics_", pattern_tag, ".csv"))
truth_csv <- file.path(out_tbl_dir, paste0("intrasample_feature_truth_", pattern_tag, ".csv"))
rep_corr_csv <- file.path(out_tbl_dir, paste0("intrasample_replicate_correlations_", pattern_tag, ".csv"))
cover_feature_csv <- file.path(out_tbl_dir, paste0("intrasample_feature_coverage95_", pattern_tag, ".csv"))

sim_safe_write_csv(metrics, metrics_csv)
sim_safe_write_csv(truth_df, truth_csv)
sim_safe_write_csv(rep_corr, rep_corr_csv)

cover_feature <- joined %>%
  mutate(
    se_hat = sqrt(pmax(vi, 0)),
    lower_95 = yi - 1.96 * se_hat,
    upper_95 = yi + 1.96 * se_hat,
    covered_95 = is.finite(empirical_mean) & is.finite(lower_95) & is.finite(upper_95) &
      empirical_mean >= lower_95 & empirical_mean <= upper_95
  ) %>%
  group_by(pattern, feature_id, ct1, ct2, radius_um, stat) %>%
  summarise(
    coverage_95 = mean(covered_95, na.rm = TRUE),
    n_reps = dplyr::n(),
    .groups = "drop"
  )
sim_safe_write_csv(cover_feature, cover_feature_csv)

rep_corr_medians <- rep_corr %>%
  group_by(pattern) %>%
  summarise(
    median_corr = median(pearson_corr, na.rm = TRUE),
    y_max = max(pearson_corr, na.rm = TRUE),
    y_min = min(pearson_corr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_pos = y_max + pmax(0.03, 0.08 * (y_max - y_min)),
    label = sprintf("median = %.2f", median_corr)
  )

p_corr <- ggplot(rep_corr, aes(x = pattern, y = pearson_corr, fill = pattern)) +
  geom_violin(trim = FALSE, alpha = 0.55, color = "black") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_jitter(
    aes(color = pattern),
    width = 0.12,
    height = 0,
    alpha = 0.5,
    size = 1.1,
    show.legend = FALSE
  ) +
  geom_text(
    data = rep_corr_medians,
    aes(x = pattern, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 3.3,
    vjust = 0
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Bootstrap vs empirical variance correlation by simulation replicate",
    subtitle = "Points are individual simulation replicates",
    x = "Pattern",
    y = "Pearson correlation"
  ) +
  theme_classic() +
  theme(legend.position = "none")

p_scatter <- ggplot(feature_scatter, aes(x = log_bootstrap_variance, y = log_empirical_variance)) +
  geom_point(alpha = 0.65, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  facet_wrap(~ pattern, scales = "free") +
  labs(
    title = "Feature-level variance calibration",
    subtitle = "Each point is one feature; axes are log1p-transformed",
    x = "log1p(mean bootstrap variance)",
    y = "log1p(empirical variance)"
  ) +
  theme_bw()

corr_summary <- bind_rows(
  metrics %>%
    transmute(pattern, metric = "Pearson", value = pearson_feature_correlation),
  metrics %>%
    transmute(pattern, metric = "Spearman", value = spearman_feature_correlation)
)

p_corr_summary <- ggplot(corr_summary, aes(x = metric, y = value, fill = metric)) +
  geom_col(width = 0.62, alpha = 0.7, color = "black") +
  geom_point(size = 2.0, color = "black") +
  facet_wrap(~ pattern, scales = "free_x") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Feature-level correlation summaries",
    subtitle = "Correlation between empirical variance and mean bootstrap variance",
    x = NULL,
    y = "Correlation"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p_rmse <- ggplot(metrics, aes(x = pattern, y = log_rmse_variance, fill = pattern)) +
  geom_col(width = 0.62, alpha = 0.75, color = "black") +
  geom_text(aes(label = sprintf("%.3f", log_rmse_variance)), vjust = -0.5, size = 3.2) +
  labs(
    title = "Log-scale RMSE of variance calibration",
    subtitle = "RMSE between log1p(empirical variance) and log1p(mean bootstrap variance)",
    x = "Pattern",
    y = "log-RMSE"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p_slope <- ggplot(metrics, aes(x = pattern, y = variance_calibration_slope, color = pattern)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.0) +
  geom_text(aes(label = sprintf("%.2f", variance_calibration_slope)), nudge_y = 0.04, size = 3.2) +
  labs(
    title = "Variance calibration slope (through origin)",
    subtitle = "Dashed line = ideal slope of 1",
    x = "Pattern",
    y = "Slope"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p_cov_feature <- ggplot(cover_feature, aes(x = pattern, y = coverage_95, fill = pattern)) +
  geom_violin(trim = FALSE, alpha = 0.55, color = "black") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_jitter(aes(color = pattern), width = 0.12, height = 0, alpha = 0.45, size = 1.0, show.legend = FALSE) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  geom_point(
    data = metrics,
    aes(x = pattern, y = mean_interval_coverage_95),
    inherit.aes = FALSE,
    shape = 23,
    size = 3.0,
    stroke = 0.9,
    fill = "white",
    color = "black"
  ) +
  geom_text(
    data = metrics,
    aes(x = pattern, y = mean_interval_coverage_95, label = sprintf("mean=%.2f", mean_interval_coverage_95)),
    inherit.aes = FALSE,
    nudge_y = 0.03,
    size = 3.1
  ) +
  labs(
    title = "95% interval coverage",
    subtitle = "Violin/box/jitter: feature-level coverage; dashed line = nominal 0.95; diamond = mean coverage",
    x = "Pattern",
    y = "Coverage"
  ) +
  theme_bw() +
  theme(legend.position = "none")

corr_png <- file.path(out_fig_dir, paste0("intrasample_replicate_corr_violin_", pattern_tag, ".png"))
corr_pdf <- file.path(out_fig_dir, paste0("intrasample_replicate_corr_violin_", pattern_tag, ".pdf"))
scatter_png <- file.path(out_fig_dir, paste0("intrasample_feature_calibration_scatter_", pattern_tag, ".png"))
scatter_pdf <- file.path(out_fig_dir, paste0("intrasample_feature_calibration_scatter_", pattern_tag, ".pdf"))
corr_summary_png <- file.path(out_fig_dir, paste0("intrasample_correlation_summary_", pattern_tag, ".png"))
corr_summary_pdf <- file.path(out_fig_dir, paste0("intrasample_correlation_summary_", pattern_tag, ".pdf"))
rmse_png <- file.path(out_fig_dir, paste0("intrasample_logrmse_summary_", pattern_tag, ".png"))
rmse_pdf <- file.path(out_fig_dir, paste0("intrasample_logrmse_summary_", pattern_tag, ".pdf"))
slope_png <- file.path(out_fig_dir, paste0("intrasample_slope_summary_", pattern_tag, ".png"))
slope_pdf <- file.path(out_fig_dir, paste0("intrasample_slope_summary_", pattern_tag, ".pdf"))
cov_png <- file.path(out_fig_dir, paste0("intrasample_coverage95_summary_", pattern_tag, ".png"))
cov_pdf <- file.path(out_fig_dir, paste0("intrasample_coverage95_summary_", pattern_tag, ".pdf"))

ggsave(corr_png, p_corr, width = 9, height = 5.4, dpi = 300)
ggsave(corr_pdf, p_corr, width = 9, height = 5.4)
ggsave(scatter_png, p_scatter, width = 10.5, height = 6.8, dpi = 300)
ggsave(scatter_pdf, p_scatter, width = 10.5, height = 6.8)
ggsave(corr_summary_png, p_corr_summary, width = 9.2, height = 5.8, dpi = 300)
ggsave(corr_summary_pdf, p_corr_summary, width = 9.2, height = 5.8)
ggsave(rmse_png, p_rmse, width = 8.6, height = 5.4, dpi = 300)
ggsave(rmse_pdf, p_rmse, width = 8.6, height = 5.4)
ggsave(slope_png, p_slope, width = 8.6, height = 5.4, dpi = 300)
ggsave(slope_pdf, p_slope, width = 8.6, height = 5.4)
ggsave(cov_png, p_cov_feature, width = 9.4, height = 5.8, dpi = 300)
ggsave(cov_pdf, p_cov_feature, width = 9.4, height = 5.8)

if (identical(pattern_tag, "combined")) {
  sim_safe_write_csv(metrics, file.path(out_tbl_root, "intrasample_calibration_metrics.csv"))
  sim_safe_write_csv(truth_df, file.path(out_tbl_root, "intrasample_feature_truth.csv"))
  sim_safe_write_csv(rep_corr, file.path(out_tbl_root, "intrasample_replicate_correlations.csv"))
  sim_safe_write_csv(cover_feature, file.path(out_tbl_root, "intrasample_feature_coverage95.csv"))
  ggsave(file.path(out_fig_root, "intrasample_replicate_corr_violin.png"), p_corr, width = 9, height = 5.4, dpi = 300)
  ggsave(file.path(out_fig_root, "intrasample_replicate_corr_violin.pdf"), p_corr, width = 9, height = 5.4)
  ggsave(file.path(out_fig_root, "intrasample_feature_calibration_scatter.png"), p_scatter, width = 10.5, height = 6.8, dpi = 300)
  ggsave(file.path(out_fig_root, "intrasample_feature_calibration_scatter.pdf"), p_scatter, width = 10.5, height = 6.8)
  ggsave(file.path(out_fig_root, "intrasample_correlation_summary.png"), p_corr_summary, width = 9.2, height = 5.8, dpi = 300)
  ggsave(file.path(out_fig_root, "intrasample_correlation_summary.pdf"), p_corr_summary, width = 9.2, height = 5.8)
  ggsave(file.path(out_fig_root, "intrasample_logrmse_summary.png"), p_rmse, width = 8.6, height = 5.4, dpi = 300)
  ggsave(file.path(out_fig_root, "intrasample_logrmse_summary.pdf"), p_rmse, width = 8.6, height = 5.4)
  ggsave(file.path(out_fig_root, "intrasample_slope_summary.png"), p_slope, width = 8.6, height = 5.4, dpi = 300)
  ggsave(file.path(out_fig_root, "intrasample_slope_summary.pdf"), p_slope, width = 8.6, height = 5.4)
  ggsave(file.path(out_fig_root, "intrasample_coverage95_summary.png"), p_cov_feature, width = 9.4, height = 5.8, dpi = 300)
  ggsave(file.path(out_fig_root, "intrasample_coverage95_summary.pdf"), p_cov_feature, width = 9.4, height = 5.8)
}

message("Saved intrasample analysis outputs for: ", pattern_tag)
message("  metrics: ", metrics_csv)
message("  feature truth: ", truth_csv)
message("  replicate correlations: ", rep_corr_csv)
message("  feature coverage: ", cover_feature_csv)
message("  figures: ", corr_png, ", ", scatter_png, ", ", corr_summary_png, ", ", rmse_png, ", ", slope_png, ", ", cov_png)
