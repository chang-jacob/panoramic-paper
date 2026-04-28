#!/usr/bin/env Rscript

# Simulation Script: Intrasample Variance Simulation
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Simulate spatial point patterns under uniform/gradient/clustered regimes.
# - Run PANORAMIC local composition enrichment with block bootstrap per replicate.
# - Export long-form simulation outputs, previews, and run parameter artifacts.

suppressPackageStartupMessages({
  library(panoramic)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
})

source("analysis/simulation/sim_utils.R")

project_root <- sim_find_project_root()
paths <- sim_load_paths(project_root)

out_data_root <- sim_mkdir(file.path(paths$data, "interim", "simulation", "intrasample"))
out_tbl_root <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intrasample"))
out_fig_root <- sim_mkdir(file.path(paths$output, "figures", "simulation", "intrasample"))

cfg <- list(
  n_sims = 300L,
  boot_nsim = 100L,
  seed_base = 123L,
  pattern_arg = "all",
  n_workers = 10L,
  square_size = 600,
  radii_um = 25,
  boot_mode = "block",
  tile_size = 62.5,
  crc_spe_key = "one_per_patient"
)

pattern_env <- Sys.getenv("PANORAMIC_INTRASAMPLE_PATTERN", unset = "")
if (nzchar(pattern_env)) {
  cfg$pattern_arg <- pattern_env
}

cfg$n_sims <- as.integer(cfg$n_sims)
cfg$boot_nsim <- as.integer(cfg$boot_nsim)
cfg$seed_base <- as.integer(cfg$seed_base)
cfg$n_workers <- as.integer(cfg$n_workers)
cfg$pattern_arg <- tolower(trimws(as.character(cfg$pattern_arg)))
cfg$tile_size <- 2.5 * cfg$radii_um

outer_bpparam <- if (cfg$n_workers > 1L) {
  BiocParallel::SnowParam(
    workers = cfg$n_workers,
    progressbar = TRUE,
    type = "SOCK",
    tasks = cfg$n_sims
  )
} else {
  BiocParallel::SerialParam(progressbar = TRUE)
}

sim_save_preview_panels <- function(points_plot, lollipop_plot, png_file, pdf_file,
                                    width = 11.6, height = 6.2,
                                    panel_widths = c(1.6, 1.0)) {
  draw_panels <- function() {
    grid::grid.newpage()
    lay <- grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(panel_widths, "null"))
    grid::pushViewport(grid::viewport(layout = lay))
    print(points_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(lollipop_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::popViewport()
  }

  grDevices::png(filename = png_file, width = width, height = height, units = "in", res = 300)
  draw_panels()
  grDevices::dev.off()

  grDevices::pdf(file = pdf_file, width = width, height = height, onefile = FALSE)
  draw_panels()
  grDevices::dev.off()
}

cell_types <- c("A", "B", "C", "Rare 1", "Rare 2")
lambda_baseline <- c(
  "A" = 0.2 / 1000,
  "B" = 0.15 / 1000,
  "C" = 0.1 / 1000,
  "Rare 1" = 0.01 / 1000,
  "Rare 2" = 0.005 / 1000
)

crc_spe_path <- sim_crc_tma_spe_list_path(paths, key = as.character(cfg$crc_spe_key))
spe_list_crc <- readRDS(crc_spe_path)
avg_crc_cells <- mean(vapply(spe_list_crc, function(spe) as.numeric(ncol(spe)), numeric(1)), na.rm = TRUE)
target_total_density <- avg_crc_cells / (cfg$square_size^2)
density_scale <- target_total_density / sum(lambda_baseline)
lambda_scaled <- lambda_baseline * density_scale

set.seed(cfg$seed_base + 8001L)
gradient_rate <- stats::runif(1, min = 1e-10, max = 0.005)

set.seed(cfg$seed_base + 9001L)
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

all_patterns <- c("uniform", "opposing_gradient", "clustered")
if (identical(cfg$pattern_arg, "all")) {
  patterns <- all_patterns
} else {
  patterns <- unique(trimws(strsplit(cfg$pattern_arg, ",", fixed = TRUE)[[1]]))
  patterns <- patterns[patterns %in% all_patterns]
}
if (length(patterns) == 0L) {
  stop(
    "Invalid pattern_arg: ", cfg$pattern_arg,
    ". Use one of: all, ", paste(all_patterns, collapse = ", ")
  )
}
pattern_tag <- if (length(patterns) == 1L) patterns else "combined"
out_data_dir <- sim_mkdir(file.path(out_data_root, pattern_tag))
out_tbl_dir <- sim_mkdir(file.path(out_tbl_root, pattern_tag))
out_fig_dir <- sim_mkdir(file.path(out_fig_root, pattern_tag))

run_params <- list(
  n_sims = cfg$n_sims,
  boot_nsim = cfg$boot_nsim,
  seed_base = cfg$seed_base,
  pattern_arg = cfg$pattern_arg,
  radii_um = cfg$radii_um,
  boot_mode = cfg$boot_mode,
  tile_size = cfg$tile_size,
  square_size = cfg$square_size,
  avg_crc_cells = avg_crc_cells,
  density_scale = density_scale,
  n_workers = cfg$n_workers,
  crc_spe_path = crc_spe_path
)

print(run_params[c("n_sims", "boot_nsim", "radii_um", "boot_mode", "tile_size", "square_size", "avg_crc_cells", "density_scale", "n_workers", "pattern_arg")])

sim_safe_save_rds(run_params, file.path(out_data_dir, paste0("intrasample_params_", pattern_tag, ".rds")))
sim_safe_write_csv(
  data.frame(
    square_size_um = cfg$square_size,
    avg_crc_cells = avg_crc_cells,
    target_total_density = target_total_density,
    baseline_total_density = sum(lambda_baseline),
    density_scale = density_scale,
    expected_cells_in_sim = target_total_density * cfg$square_size^2,
    stringsAsFactors = FALSE
  ),
  file.path(out_tbl_dir, paste0("density_match_summary_", pattern_tag, ".csv"))
)

# Also write simplified names within the run-specific subdirectory.
sim_safe_save_rds(run_params, file.path(out_data_dir, "intrasample_params.rds"))
sim_safe_write_csv(
  data.frame(
    square_size_um = cfg$square_size,
    avg_crc_cells = avg_crc_cells,
    target_total_density = target_total_density,
    baseline_total_density = sum(lambda_baseline),
    density_scale = density_scale,
    expected_cells_in_sim = target_total_density * cfg$square_size^2,
    stringsAsFactors = FALSE
  ),
  file.path(out_tbl_dir, "density_match_summary.csv")
)

all_results <- list()

for (pat in patterns) {
  pat_idx <- match(pat, all_patterns)
  message("\nPattern: ", pat)
  pat_data_dir <- sim_mkdir(file.path(out_data_root, pat))
  pat_tbl_dir <- sim_mkdir(file.path(out_tbl_root, pat))
  pat_fig_dir <- sim_mkdir(file.path(out_fig_root, pat))

  preview_seed <- cfg$seed_base + 1 + pat_idx * 100000L
  preview_pp <- sim_simulate_pattern_ppp(
    pattern = pat,
    square_size = cfg$square_size,
    cell_types = cell_types,
    params = pattern_params[[pat]],
    seed = preview_seed
  )
  preview_df <- as.data.frame(preview_pp)
  names(preview_df)[names(preview_df) == "marks"] <- "cell_type"
  if (is.list(preview_df$cell_type)) {
    preview_df$cell_type <- unlist(preview_df$cell_type, use.names = FALSE)
  }
  preview_df$cell_type <- as.character(preview_df$cell_type)
  preview_df$pattern <- pat
  preview_df$replicate <- 1L

  preview_points_plot <- ggplot(preview_df, aes(x = x, y = y, color = cell_type)) +
    geom_point(size = 0.55, alpha = 0.85) +
    scale_x_continuous(limits = c(0, cfg$square_size), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, cfg$square_size), expand = expansion(mult = 0)) +
    coord_equal() +
    theme_bw() +
    labs(
      title = paste0("Simulated pattern preview: ", pat),
      subtitle = "Replicate 1",
      x = "x (um)",
      y = "y (um)",
      color = "Cell type"
    )

  preview_tab <- as.data.frame(table(factor(preview_df$cell_type, levels = cell_types)), stringsAsFactors = FALSE)
  names(preview_tab) <- c("cell_type", "cell_count")
  preview_counts <- preview_tab %>% mutate(cell_count = as.integer(cell_count))

  total_cells <- sum(preview_counts$cell_count)
  denom_cells <- if (total_cells > 0) total_cells else 1L
  preview_counts <- preview_counts %>%
    mutate(
      proportion = cell_count / denom_cells,
      label = sprintf("%.1f%% (%d)", 100 * proportion, cell_count),
      cell_type = factor(cell_type, levels = rev(cell_types))
    )

  lollipop_xmax <- max(preview_counts$proportion, na.rm = TRUE)
  if (!is.finite(lollipop_xmax) || lollipop_xmax <= 0) lollipop_xmax <- 1
  lollipop_xmax <- min(1, lollipop_xmax * 1.3 + 0.03)

  preview_lollipop_plot <- ggplot(preview_counts, aes(x = proportion, y = cell_type, color = cell_type)) +
    geom_segment(aes(x = 0, xend = proportion, yend = cell_type), linewidth = 0.9, alpha = 0.8, show.legend = FALSE) +
    geom_point(size = 2.6, show.legend = FALSE) +
    geom_text(aes(label = label), nudge_x = lollipop_xmax * 0.03, hjust = 0, size = 3.35, color = "black", show.legend = FALSE) +
    scale_x_continuous(
      limits = c(0, lollipop_xmax),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    theme_bw() +
    labs(
      title = "Cell type composition",
      subtitle = sprintf("Total cells: %d", total_cells),
      x = "Proportion",
      y = NULL
    ) +
    theme(legend.position = "none")

  preview_png <- file.path(pat_fig_dir, paste0("preview_", pat, "_rep1.png"))
  preview_pdf <- file.path(pat_fig_dir, paste0("preview_", pat, "_rep1.pdf"))
  preview_csv <- file.path(pat_tbl_dir, paste0("preview_", pat, "_rep1_points.csv"))
  preview_comp_csv <- file.path(pat_tbl_dir, paste0("preview_", pat, "_rep1_composition.csv"))

  sim_save_preview_panels(
    points_plot = preview_points_plot,
    lollipop_plot = preview_lollipop_plot,
    png_file = preview_png,
    pdf_file = preview_pdf
  )
  sim_safe_write_csv(preview_df, preview_csv)
  sim_safe_write_csv(preview_counts, preview_comp_csv)
  message("  Preview saved: ", preview_png)

  sim_fun <- function(i) {
    sim_seed <- cfg$seed_base + i + pat_idx * 100000L
    sample_id <- sprintf("%s_sim_%04d", pat, i)
    pp <- sim_simulate_pattern_ppp(
      pattern = pat,
      square_size = cfg$square_size,
      cell_types = cell_types,
      params = pattern_params[[pat]],
      seed = sim_seed
    )
    spe <- sim_ppp_to_spe(pp, sample_id = sample_id)
    se_stats <- sim_run_spatialstats_one(
      spe = spe,
      sample_id = sample_id,
      cell_type_col = "cell_type",
      radii_um = cfg$radii_um,
      stat = "local_comp_enrichment",
      nsim = cfg$boot_nsim,
      correction = "translate",
      min_cells = 1L,
      window = "rect",
      seed = sim_seed,
      boot = cfg$boot_mode,
      tile_size = cfg$tile_size,
      BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
    )

    out <- sim_flatten_spatialstats(se_stats)
    out$pattern <- pat
    out$replicate <- i
    out$sample_id <- sample_id
    out$square_size <- cfg$square_size
    out$density_scale <- density_scale
    out
  }

  one_pattern <- BiocParallel::bplapply(
    X = seq_len(cfg$n_sims),
    FUN = sim_fun,
    BPPARAM = outer_bpparam
  )

  all_results[[pat]] <- bind_rows(one_pattern)

  # Save per-pattern simulation result files into pattern-specific directories.
  sim_safe_save_rds(all_results[[pat]], file.path(pat_data_dir, "intrasample_spatialstats_long.rds"))
  sim_safe_write_csv(all_results[[pat]], file.path(pat_tbl_dir, "intrasample_spatialstats_long.csv"))
  message(sprintf("  Completed replicates: %d / %d", cfg$n_sims, cfg$n_sims))
}

df_long <- bind_rows(all_results)

long_rds <- file.path(out_data_dir, paste0("intrasample_spatialstats_long_", pattern_tag, ".rds"))
long_csv <- file.path(out_tbl_dir, paste0("intrasample_spatialstats_long_", pattern_tag, ".csv"))

sim_safe_save_rds(df_long, long_rds)
sim_safe_write_csv(df_long, long_csv)
sim_safe_save_rds(df_long, file.path(out_data_dir, "intrasample_spatialstats_long.rds"))
sim_safe_write_csv(df_long, file.path(out_tbl_dir, "intrasample_spatialstats_long.csv"))

canonical_rds <- file.path(out_data_root, "intrasample_spatialstats_long.rds")
canonical_csv <- file.path(out_tbl_root, "intrasample_spatialstats_long.csv")
canonical_params <- file.path(out_data_root, "intrasample_params.rds")
canonical_density <- file.path(out_tbl_root, "density_match_summary.csv")

sim_safe_save_rds(df_long, canonical_rds)
sim_safe_write_csv(df_long, canonical_csv)
sim_safe_save_rds(run_params, canonical_params)
sim_safe_write_csv(
  data.frame(
    square_size_um = cfg$square_size,
    avg_crc_cells = avg_crc_cells,
    target_total_density = target_total_density,
    baseline_total_density = sum(lambda_baseline),
    density_scale = density_scale,
    expected_cells_in_sim = target_total_density * cfg$square_size^2,
    stringsAsFactors = FALSE
  ),
  canonical_density
)

message("\nSaved intrasample simulation outputs:")
message("  run directory: ", out_data_dir)
message("  params: ", file.path(out_data_dir, paste0("intrasample_params_", pattern_tag, ".rds")))
message("  long rds: ", long_rds)
message("  long csv: ", long_csv)
