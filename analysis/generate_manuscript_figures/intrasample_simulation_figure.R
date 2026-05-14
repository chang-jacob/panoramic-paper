#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

find_project_root <- function(start = getwd(), max_depth = 10L) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_depth)) {
    if (file.exists(file.path(cur, "config", "default.yml"))) return(cur)
    nxt <- dirname(cur)
    if (identical(nxt, cur)) break
    cur <- nxt
  }
  stop("Could not find project root containing config/default.yml from: ", start)
}

pick_existing <- function(paths, label = "file") {
  for (p in paths) {
    if (file.exists(p)) return(p)
  }
  stop("Could not find ", label, ". Checked:\n- ", paste(paths, collapse = "\n- "))
}

save_plot_formats <- function(plot_obj, stem, width, height, dpi = 300) {
  ggsave(
    filename = paste0(stem, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  ggsave(
    filename = paste0(stem, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in"
  )
}

project_root <- find_project_root()
source(file.path(project_root, "analysis", "simulation", "sim_utils.R"))
paths <- sim_load_paths(project_root)

in_tbl_root <- file.path(paths$output, "tables", "simulation", "intrasample")
in_data_root <- file.path(paths$data, "interim", "simulation", "intrasample")
out_tbl_dir <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intrasample", "publication"))
out_fig_dir <- sim_mkdir(file.path(paths$output, "figures", "simulation", "intrasample", "publication"))

patterns <- c("uniform", "opposing_gradient", "clustered")
all_patterns <- c("uniform", "opposing_gradient", "clustered")
n_examples_per_pattern <- 3L

pattern_labels <- c(
  uniform = "Uniform",
  opposing_gradient = "Opposing Gradient",
  clustered = "Clustered"
)
pattern_palette <- c(
  "Uniform" = "gray",
  "Opposing Gradient" = "gray",
  "Clustered" = "gray"
)
cell_type_palette <- c(
  "A" = "#1F77B4",
  "B" = "#FF7F0E",
  "C" = "#B3B3B3",
  "Rare 1" = "#D62728",
  "Rare 2" = "#2CA02C"
)

params_path <- pick_existing(
  c(
    file.path(in_data_root, "intrasample_params.rds"),
    file.path(in_data_root, "intrasample_params_combined.rds"),
    file.path(in_data_root, "intrasample_params_cluster.rds"),
    file.path(in_data_root, "intrasample_params_diverging.rds")
  ),
  label = "intrasample params RDS"
)
run_params <- readRDS(params_path)

seed_base <- as.integer(if (!is.null(run_params$seed_base)) run_params$seed_base else 123L)
square_size <- as.numeric(if (!is.null(run_params$square_size)) run_params$square_size else 600)
density_scale <- as.numeric(if (!is.null(run_params$density_scale)) run_params$density_scale else 1)

cell_types <- c("A", "B", "C", "Rare 1", "Rare 2")
lambda_baseline <- c(
  "A" = 0.2 / 1000,
  "B" = 0.15 / 1000,
  "C" = 0.1 / 1000,
  "Rare 1" = 0.01 / 1000,
  "Rare 2" = 0.005 / 1000
)
lambda_scaled <- lambda_baseline * density_scale

set.seed(seed_base + 8001L)
gradient_rate <- stats::runif(1, min = 1e-10, max = 0.005)

set.seed(seed_base + 9001L)
cluster_parent_lambda <- (0.00001 + stats::runif(1, 0, 0.0001)) * density_scale
cluster_offspring_mean <- as.numeric(3 + stats::rnbinom(1, size = 10, mu = 2))
cluster_offspring_sd <- 175

pattern_params <- list(
  uniform = list(lambda = lambda_scaled),
  opposing_gradient = list(
    gradient_rate = gradient_rate,
    target_density_a = as.numeric(lambda_scaled["A"]),
    target_density_b = as.numeric(lambda_scaled["B"]),
    lambda_c = as.numeric(lambda_scaled["C"]),
    lambda_rare1 = as.numeric(lambda_scaled["Rare 1"]),
    lambda_rare2 = as.numeric(lambda_scaled["Rare 2"])
  ),
  clustered = list(
    parent_lambda = cluster_parent_lambda,
    offspring_mean = cluster_offspring_mean,
    offspring_sd = cluster_offspring_sd,
    lambda_c = as.numeric(lambda_scaled["C"]),
    lambda_rare1 = as.numeric(lambda_scaled["Rare 1"]),
    lambda_rare2 = as.numeric(lambda_scaled["Rare 2"])
  )
)

simulate_example_points <- function(pattern, rep_id) {
  pat_idx <- match(pattern, all_patterns)
  sim_seed <- seed_base + as.integer(rep_id) + pat_idx * 100000L
  pp <- sim_simulate_pattern_ppp(
    pattern = pattern,
    square_size = square_size,
    cell_types = cell_types,
    params = pattern_params[[pattern]],
    seed = sim_seed
  )
  df <- as.data.frame(pp)
  if ("marks" %in% colnames(df)) names(df)[names(df) == "marks"] <- "cell_type"
  if (is.list(df$cell_type)) df$cell_type <- unlist(df$cell_type, use.names = FALSE)

  df |>
    mutate(
      pattern = pattern,
      replicate = as.integer(rep_id),
      pattern_display = factor(pattern_labels[pattern], levels = pattern_labels[patterns]),
      cell_type = factor(as.character(cell_type), levels = names(cell_type_palette))
    )
}

preview_points <- bind_rows(lapply(patterns, function(pat) {
  bind_rows(lapply(seq_len(n_examples_per_pattern), function(rep_id) {
    simulate_example_points(pat, rep_id)
  }))
}))

load_pattern_long <- function(pattern) {
  csv_path <- pick_existing(
    c(
      file.path(in_tbl_root, pattern, "intrasample_spatialstats_long.csv"),
      file.path(in_tbl_root, pattern, paste0("intrasample_spatialstats_long_", pattern, ".csv"))
    ),
    label = paste0("intrasample_spatialstats_long CSV for pattern=", pattern)
  )
  df <- readr::read_csv(csv_path, show_col_types = FALSE)
  needed <- c("pattern", "replicate", "feature_id", "ct1", "ct2", "radius_um", "stat", "yi", "vi")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0L) {
    stop("Long CSV missing required columns (", pattern, "): ", paste(miss, collapse = ", "))
  }
  df
}

df_long <- bind_rows(lapply(patterns, load_pattern_long))
truth_df <- sim_build_feature_truth(df_long)
rep_corr <- sim_per_replicate_variance_correlation(df_long, truth_df)
rep_corr <- rep_corr |>
  mutate(pattern_display = factor(pattern_labels[pattern], levels = pattern_labels[patterns]))

rep_corr_summary <- rep_corr |>
  group_by(pattern_display) |>
  summarise(
    n_reps = dplyr::n(),
    median_corr = median(pearson_corr, na.rm = TRUE),
    mean_corr = mean(pearson_corr, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(rep_corr, file.path(out_tbl_dir, "intrasample_publication_replicate_correlations.csv"))
readr::write_csv(rep_corr_summary, file.path(out_tbl_dir, "intrasample_publication_replicate_correlations_summary.csv"))

make_pattern_plot <- function(df_pattern) {
  ggplot(df_pattern, aes(x = x, y = y, color = cell_type)) +
    geom_point(size = 0.42, alpha = 0.88) +
    scale_color_manual(values = cell_type_palette, drop = FALSE) +
    coord_equal(
      xlim = c(0, square_size),
      ylim = c(0, square_size),
      expand = FALSE
    ) +
    theme_void(base_size = 8) +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

pattern_plots <- list()
for (pat in patterns) {
  for (rep_id in seq_len(n_examples_per_pattern)) {
    p_pat <- make_pattern_plot(
      df_pattern = preview_points |> filter(pattern == pat, replicate == rep_id)
    )
    pattern_plots[[length(pattern_plots) + 1L]] <- p_pat
    save_plot_formats(
      plot_obj = p_pat,
      stem = file.path(out_fig_dir, paste0("intrasample_pattern_example_", pat, "_rep", rep_id)),
      width = 1,
      height = 1
    )
  }
}

pattern_examples_panel <- cowplot::plot_grid(
  plotlist = pattern_plots,
  ncol = n_examples_per_pattern
)

rep_corr_medians <- rep_corr |>
  group_by(pattern_display) |>
  summarise(
    median_corr = median(pearson_corr, na.rm = TRUE),
    y_max = max(pearson_corr, na.rm = TRUE),
    y_min = min(pearson_corr, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    y_pos = y_max + pmax(0.03, 0.08 * (y_max - y_min)),
    label = sprintf("median = %.2f", median_corr)
  )

p_corr <- ggplot(rep_corr, aes(x = pattern_display, y = pearson_corr, fill = pattern_display)) +
  geom_violin(trim = FALSE, alpha = 0.55, color = "black") +
  geom_jitter(
    aes(color = pattern_display),
    width = 0.12,
    height = 0,
    alpha = 0.45,
    size = 0.65,
    show.legend = FALSE
  ) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_text(
    data = rep_corr_medians,
    aes(x = pattern_display, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    vjust = 0
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  scale_fill_manual(values = pattern_palette, guide = "none") +
  scale_color_manual(values = pattern_palette, guide = "none") +
  labs(
    x = "Pattern",
    y = "Pearson correlation"
  ) +
  theme_classic(base_size = 11)

combined_panel <- cowplot::plot_grid(
  pattern_examples_panel,
  p_corr,
  ncol = 1,
  labels = c("A", "B"),
  rel_heights = c(2.1, 1)
)

print(pattern_examples_panel)
print(p_corr)
print(combined_panel)

save_plot_formats(
  plot_obj = pattern_examples_panel,
  stem = file.path(out_fig_dir, "intrasample_pattern_examples_9panel"),
  width = 3,
  height = 3
)
save_plot_formats(
  plot_obj = p_corr,
  stem = file.path(out_fig_dir, "intrasample_bootstrap_vs_groundtruth_corr_violin"),
  width = 3.25,
  height = 4
)
save_plot_formats(
  plot_obj = combined_panel,
  stem = file.path(out_fig_dir, "intrasample_simulation_publication_panel"),
  width = 10.5,
  height = 13.0
)

message("Saved intrasample publication figures and summary tables.")
message("  Figure dir: ", out_fig_dir)
message("  Table dir: ", out_tbl_dir)
