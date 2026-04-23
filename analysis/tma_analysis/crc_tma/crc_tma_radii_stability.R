#!/usr/bin/env Rscript

# CRC TMA Script: Radii Stability
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Re-run CRC TMA analysis across radii with per-radius block size.
# - Apply within-radius FDR across cell-type pairs and build stability summaries.
# - Export radii stability tables and discrete/continuous heatmaps.

suppressPackageStartupMessages({
  library(config)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(BiocParallel)
  library(panoramic)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

parse_radii <- function(x) {
  if (length(x) == 1L && is.character(x)) {
    y <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
  } else {
    y <- as.numeric(x)
  }
  y <- y[is.finite(y) & y > 0]
  sort(unique(y))
}

clean_pair_component <- function(x) {
  x <- as.character(x)
  x <- gsub("[(){}]", "", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

min_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

mean_abs_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(abs(x))
}

build_radii_stability_tables <- function(
    contrast_tbl,
    radii_expected,
    fdr_mode = "within_radius",
    alpha_pair_include = 0.10,
    sig_pairs_only = TRUE,
    max_pairs = NULL
) {
  if (!identical(fdr_mode, "within_radius")) {
    stop("fdr_mode must be 'within_radius' for this analysis script.")
  }
  req <- c("ct1", "ct2", "radius_um", "p_diff", "fdr_diff", "beta_diff", "z_diff", "coloc_direction")
  missing_req <- setdiff(req, colnames(contrast_tbl))
  if (length(missing_req) > 0L) {
    stop("contrast table missing required columns: ", paste(missing_req, collapse = ", "))
  }

  df <- contrast_tbl |>
    dplyr::transmute(
      ct1 = as.character(.data$ct1),
      ct2 = as.character(.data$ct2),
      radius_um = as.numeric(.data$radius_um),
      p_diff = as.numeric(.data$p_diff),
      fdr_diff_global = as.numeric(.data$fdr_diff),
      beta_diff = as.numeric(.data$beta_diff),
      z_diff = as.numeric(.data$z_diff),
      coloc_direction = as.character(.data$coloc_direction)
    ) |>
    dplyr::filter(is.finite(.data$radius_um))

  if (nrow(df) == 0L) {
    stop("No finite radius values in contrast table.")
  }

  df$pair <- clean_pair_component(df$coloc_direction)

  # Collapse to one row per pair/radius before multiple-testing correction.
  df <- df |>
    dplyr::arrange(.data$pair, .data$radius_um, .data$p_diff) |>
    dplyr::group_by(.data$pair, .data$radius_um) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup()

  # FDR correction among pairs within each radius.
  df$fdr_diff_within_radius <- NA_real_
  radius_vals <- sort(unique(df$radius_um))
  for (r in radius_vals) {
    idx <- which(df$radius_um == r & is.finite(df$p_diff))
    if (length(idx) > 0L) {
      df$fdr_diff_within_radius[idx] <- p.adjust(df$p_diff[idx], method = "BH")
    }
  }

  df$fdr_for_heatmap <- df$fdr_diff_within_radius

  pair_summary <- df |>
    dplyr::group_by(.data$pair) |>
    dplyr::summarise(
      n_radii_total = dplyr::n_distinct(.data$radius_um),
      n_radii_sig = dplyr::n_distinct(.data$radius_um[is.finite(.data$fdr_for_heatmap) & .data$fdr_for_heatmap <= alpha_pair_include]),
      radii_sig = paste(sort(unique(.data$radius_um[is.finite(.data$fdr_for_heatmap) & .data$fdr_for_heatmap <= alpha_pair_include])), collapse = ", "),
      min_fdr_for_heatmap = min_or_na(.data$fdr_for_heatmap),
      min_fdr_global = min_or_na(.data$fdr_diff_global),
      min_fdr_within_radius = min_or_na(.data$fdr_diff_within_radius),
      min_p = min_or_na(.data$p_diff),
      mean_abs_z = mean_abs_or_na(.data$z_diff),
      mean_abs_beta = mean_abs_or_na(.data$beta_diff),
      .groups = "drop"
    ) |>
    dplyr::mutate(radii_sig = ifelse(nzchar(.data$radii_sig), .data$radii_sig, NA_character_)) |>
    dplyr::arrange(dplyr::desc(.data$n_radii_sig), .data$min_fdr_for_heatmap, .data$min_p, .data$pair)

  pair_summary_plot <- pair_summary
  if (isTRUE(sig_pairs_only) && any(pair_summary_plot$n_radii_sig > 0L, na.rm = TRUE)) {
    pair_summary_plot <- pair_summary_plot |>
      dplyr::filter(.data$n_radii_sig > 0L)
  }
  if (!is.null(max_pairs)) {
    max_pairs_num <- as.integer(max_pairs)
    if (length(max_pairs_num) != 1L || is.na(max_pairs_num) || !is.finite(max_pairs_num) || max_pairs_num <= 0L) {
      stop("max_pairs must be NULL or a positive integer.")
    }
    if (nrow(pair_summary_plot) > max_pairs_num) {
      pair_summary_plot <- dplyr::slice_head(pair_summary_plot, n = max_pairs_num)
    }
  }
  if (nrow(pair_summary_plot) == 0L) {
    stop("No pairs available for plotting after filtering.")
  }

  pair_levels <- pair_summary_plot$pair
  radii_levels <- sort(unique(as.numeric(radii_expected)))
  if (length(radii_levels) == 0L) {
    stop("radii_expected must contain at least one numeric radius.")
  }

  plot_data <- df |>
    dplyr::filter(.data$pair %in% pair_levels) |>
    dplyr::select(
      .data$pair, .data$radius_um, .data$p_diff,
      .data$fdr_diff_global, .data$fdr_diff_within_radius, .data$fdr_for_heatmap,
      .data$beta_diff, .data$z_diff
    ) |>
    tidyr::complete(
      pair = pair_levels,
      radius_um = radii_levels,
      fill = list(
        p_diff = NA_real_,
        fdr_diff_global = NA_real_,
        fdr_diff_within_radius = NA_real_,
        fdr_for_heatmap = NA_real_,
        beta_diff = NA_real_,
        z_diff = NA_real_
      )
    )

  alpha_label <- paste0("FDR < ", format(alpha_pair_include, trim = TRUE))
  plot_data <- plot_data |>
    dplyr::mutate(
      sig_level = dplyr::case_when(
        !is.finite(.data$fdr_for_heatmap) ~ "Missing",
        .data$fdr_for_heatmap < 0.001 ~ "FDR < 0.001",
        .data$fdr_for_heatmap < 0.01 ~ "FDR < 0.01",
        .data$fdr_for_heatmap < 0.05 ~ "FDR < 0.05",
        alpha_pair_include > 0.05 & .data$fdr_for_heatmap < alpha_pair_include ~ alpha_label,
        TRUE ~ "Not sig."
      ),
      p_level = dplyr::case_when(
        !is.finite(.data$p_diff) ~ "Missing",
        .data$p_diff < 0.001 ~ "p < 0.001",
        .data$p_diff < 0.01 ~ "p < 0.01",
        .data$p_diff <= 0.05 ~ "p <= 0.05",
        TRUE ~ "p >= 0.05"
      ),
      fdr_sig_symbol = dplyr::case_when(
        !is.finite(.data$fdr_for_heatmap) ~ "",
        .data$fdr_for_heatmap < 0.001 ~ "***",
        .data$fdr_for_heatmap < 0.01 ~ "**",
        .data$fdr_for_heatmap < 0.05 ~ "*",
        alpha_pair_include > 0.05 & .data$fdr_for_heatmap < alpha_pair_include ~ ".",
        TRUE ~ ""
      ),
      fdr_sig_label = ifelse(nzchar(.data$fdr_sig_symbol), .data$fdr_sig_symbol, NA_character_),
      p_plot = ifelse(is.finite(.data$p_diff), pmax(.data$p_diff, .Machine$double.xmin), NA_real_),
      p_plot_sig = ifelse(is.finite(.data$p_diff) & .data$p_diff <= 0.05, pmax(.data$p_diff, 1e-4), NA_real_),
      p_is_nonsig = is.finite(.data$p_diff) & .data$p_diff > 0.05,
      fdr_capped = ifelse(is.finite(.data$fdr_for_heatmap), pmin(.data$fdr_for_heatmap, alpha_pair_include), NA_real_),
      pair = factor(.data$pair, levels = rev(pair_levels)),
      radius_um = factor(.data$radius_um, levels = radii_levels)
    )

  sig_levels <- c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05")
  if (alpha_pair_include > 0.05) sig_levels <- c(sig_levels, alpha_label)
  sig_levels <- c(sig_levels, "Not sig.", "Missing")
  plot_data$sig_level <- factor(plot_data$sig_level, levels = unique(sig_levels))
  plot_data$p_level <- factor(plot_data$p_level, levels = c("p < 0.001", "p < 0.01", "p <= 0.05", "p >= 0.05", "Missing"))

  list(
    processed = df,
    pair_summary = pair_summary,
    pair_summary_plot = pair_summary_plot,
    plot_data = plot_data
  )
}

plot_radii_stability_discrete <- function(
    plot_data,
    alpha_pair_include = 0.10,
    pair_axis_text_size = 7
) {
  fill_values <- c(
    "p < 0.001" = "#b2182b",
    "p < 0.01" = "#ef8a62",
    "p <= 0.05" = "#fddbc7",
    "p >= 0.05" = "#bdbdbd",
    "Missing" = "#d9d9d9"
  )

  ggplot(plot_data, aes(x = .data$radius_um, y = .data$pair, fill = .data$p_level)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = .data$fdr_sig_label), size = 2.8, na.rm = TRUE) +
    scale_fill_manual(values = fill_values, drop = FALSE, name = "Raw p-value") +
    labs(
      x = "Radius (um)",
      y = "Cell type pair",
      title = "CRC TMA colocalization stability across radii",
      subtitle = "Cell symbols denote within-radius FDR tiers: *** < 0.001, ** < 0.01, * < 0.05"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = pair_axis_text_size),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

plot_radii_stability_continuous <- function(
    plot_data,
    alpha_pair_include = 0.10,
    pair_axis_text_size = 7
) {
  breaks <- c(1e-4, 1e-3, 1e-2, 0.05)

  ggplot(plot_data, aes(x = .data$radius_um, y = .data$pair)) +
    geom_tile(
      data = dplyr::filter(plot_data, .data$p_is_nonsig),
      fill = "#bdbdbd",
      color = "white",
      linewidth = 0.35
    ) +
    geom_tile(
      data = dplyr::filter(plot_data, !.data$p_is_nonsig),
      aes(fill = .data$p_plot_sig),
      color = "white",
      linewidth = 0.35
    ) +
    geom_text(aes(label = .data$fdr_sig_label), size = 2.8, na.rm = TRUE) +
    scale_fill_gradientn(
      colors = c("#b2182b", "#ef8a62", "#fddbc7"),
      values = scales::rescale(c(1e-4, 1e-3, 0.05), from = c(1e-4, 0.05)),
      trans = "log10",
      limits = c(1e-4, 0.05),
      breaks = breaks,
      labels = formatC(breaks, format = "f", digits = 3),
      oob = scales::squish,
      na.value = "#d9d9d9",
      name = "Raw p-value"
    ) +
    labs(
      x = "Radius (um)",
      y = "Cell type pair",
      title = "CRC TMA colocalization stability (continuous raw p-value)",
      subtitle = "Cell symbols denote within-radius FDR tiers: *** < 0.001, ** < 0.01, * < 0.05"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = pair_axis_text_size),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

# ---- Config ----
project_root <- getwd()
config_file <- file.path(project_root, "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)
source(file.path(project_root, "analysis", "tma_analysis", "tma_shared_helpers.R"))

# ---- Analysis Parameters (edit in script) ----
analysis_params <- list(
  spe_list_key = "one_per_patient",
  spe_list_path = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_one_per_patient.rds"),
  columns = list(
    cell_type = "cell_type",
    patient = "patient",
    group = "group_label",
    region = "spot"
  ),
  panoramic = list(
    radii_um = c(25, 50, 75, 100, 125, 150),
    stat = "local_comp_enrichment",
    nsim = 100L,
    correction = "translate",
    min_cells = 5L,
    concavity = 50,
    window = "convex",
    seed = 123L,
    boot = "block",
    tile_size = NULL, # NULL -> 2.5 * each radius_um (per-radius loop)
    workers = 10L,
    progressbar = TRUE
  ),
  meta = list(
    method_mv = "REML",
    group_tau2 = "none",
    vi_floor = "group_median",
    tau_structure = "patient"
  ),
  stability = list(
    fdr_mode = "within_radius", # within_radius | global
    alpha_pair_include = 0.10,
    sig_pairs_only = TRUE,
    max_pairs = NULL
  ),
  plotting = list(
    show_plots = FALSE,
    heatmap_width = 9.0,
    heatmap_height_min = 6.5,
    heatmap_height_per_pair = 0.20,
    pair_axis_text_size = 7
  )
)

# ---- Parameters ----
spe_list_key <- as.character(analysis_params$spe_list_key)
spe_list_path <- as.character(analysis_params$spe_list_path)
cell_type_col <- as.character(analysis_params$columns$cell_type)
patient_col <- as.character(analysis_params$columns$patient)
group_col <- as.character(analysis_params$columns$group)
region_col <- as.character(analysis_params$columns$region)

radii_um <- parse_radii(analysis_params$panoramic$radii_um)
if (length(radii_um) == 0L) stop("analysis_params$panoramic$radii_um must contain positive numeric radii.")
stat <- as.character(analysis_params$panoramic$stat)
nsim <- as.integer(analysis_params$panoramic$nsim)
if (length(nsim) != 1L || is.na(nsim) || !is.finite(nsim) || nsim < 1L) stop("analysis_params$panoramic$nsim must be a positive integer.")
correction <- as.character(analysis_params$panoramic$correction)
min_cells <- as.integer(analysis_params$panoramic$min_cells)
if (length(min_cells) != 1L || is.na(min_cells) || !is.finite(min_cells) || min_cells < 1L) stop("analysis_params$panoramic$min_cells must be a positive integer.")
concavity <- as.numeric(analysis_params$panoramic$concavity)
if (length(concavity) != 1L || is.na(concavity) || !is.finite(concavity) || concavity <= 0) stop("analysis_params$panoramic$concavity must be > 0.")
window <- as.character(analysis_params$panoramic$window)
seed <- as.integer(analysis_params$panoramic$seed)
if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) stop("analysis_params$panoramic$seed must be an integer.")
boot <- as.character(analysis_params$panoramic$boot)
tile_size_raw <- analysis_params$panoramic$tile_size
if (is.null(tile_size_raw) || length(tile_size_raw) == 0L) {
  tile_size_by_radius <- 2.5 * radii_um
} else {
  tile_size_input <- as.numeric(tile_size_raw)
  if (length(tile_size_input) != 1L || is.na(tile_size_input) || !is.finite(tile_size_input) || tile_size_input <= 0) {
    stop("analysis_params$panoramic$tile_size must be NULL or one positive numeric value.")
  }
  tile_size_by_radius <- rep(tile_size_input, length(radii_um))
}
workers <- as.integer(analysis_params$panoramic$workers)
if (length(workers) != 1L || is.na(workers) || !is.finite(workers) || workers < 1L) stop("analysis_params$panoramic$workers must be a positive integer.")
progressbar <- isTRUE(analysis_params$panoramic$progressbar)
BPPARAM <- BiocParallel::SnowParam(workers = workers, progressbar = progressbar)

method_mv <- as.character(analysis_params$meta$method_mv)
group_tau2 <- as.character(analysis_params$meta$group_tau2)
vi_floor <- as.character(analysis_params$meta$vi_floor)
tau_structure <- as.character(analysis_params$meta$tau_structure)

fdr_mode <- as.character(analysis_params$stability$fdr_mode)
if (!identical(fdr_mode, "within_radius")) stop("analysis_params$stability$fdr_mode must be 'within_radius'.")
alpha_pair_include <- as.numeric(analysis_params$stability$alpha_pair_include)
if (length(alpha_pair_include) != 1L || is.na(alpha_pair_include) || !is.finite(alpha_pair_include) || alpha_pair_include <= 0 || alpha_pair_include > 1) {
  stop("analysis_params$stability$alpha_pair_include must be in (0,1].")
}
sig_pairs_only <- isTRUE(analysis_params$stability$sig_pairs_only)
max_pairs <- analysis_params$stability$max_pairs

show_plots <- isTRUE(analysis_params$plotting$show_plots)
heatmap_width <- as.numeric(analysis_params$plotting$heatmap_width)
if (length(heatmap_width) != 1L || is.na(heatmap_width) || !is.finite(heatmap_width) || heatmap_width <= 0) stop("analysis_params$plotting$heatmap_width must be > 0.")
heatmap_height_min <- as.numeric(analysis_params$plotting$heatmap_height_min)
if (length(heatmap_height_min) != 1L || is.na(heatmap_height_min) || !is.finite(heatmap_height_min) || heatmap_height_min <= 0) stop("analysis_params$plotting$heatmap_height_min must be > 0.")
heatmap_height_per_pair <- as.numeric(analysis_params$plotting$heatmap_height_per_pair)
if (length(heatmap_height_per_pair) != 1L || is.na(heatmap_height_per_pair) || !is.finite(heatmap_height_per_pair) || heatmap_height_per_pair <= 0) {
  stop("analysis_params$plotting$heatmap_height_per_pair must be > 0.")
}
pair_axis_text_size <- as.numeric(analysis_params$plotting$pair_axis_text_size)
if (length(pair_axis_text_size) != 1L || is.na(pair_axis_text_size) || !is.finite(pair_axis_text_size) || pair_axis_text_size <= 0) stop("analysis_params$plotting$pair_axis_text_size must be > 0.")

# ---- Outputs ----
interim_dir <- file.path(paths$data, "interim", "crc_tma", "radii_stability")
tab_dir <- file.path(paths$output, "tables", "crc_tma", spe_list_key, "radii_stability")
fig_dir <- file.path(paths$output, "figures", "crc_tma", spe_list_key, "radii_stability")
dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stat_tag <- gsub("[^a-z0-9]+", "_", tolower(stat))
radii_tag <- paste(radii_um, collapse = "_")
tile_tag <- paste(format(signif(tile_size_by_radius, 6), scientific = FALSE, trim = TRUE), collapse = "_")
tile_slug <- gsub("[^a-z0-9]+", "", tolower(gsub("\\.", "p", tile_tag)))
run_prefix <- paste0(
  "crc_tma_radii_stability_", spe_list_key,
  "_", stat_tag,
  "_r", radii_tag,
  "_nsim", nsim,
  "_", boot,
  "_tile", tile_slug,
  "_mv"
)

prep_path <- file.path(interim_dir, paste0(run_prefix, "_prep.rds"))
stats_path <- file.path(interim_dir, paste0(run_prefix, "_se_stats.rds"))
meta_path <- file.path(interim_dir, paste0(run_prefix, "_se_meta.rds"))
runtime_csv <- file.path(tab_dir, paste0(run_prefix, "_runtime_metrics.csv"))

spatial_csv <- file.path(tab_dir, paste0(run_prefix, "_spatialstats.csv"))
meta_csv <- file.path(tab_dir, paste0(run_prefix, "_meta.csv"))
contrast_csv <- file.path(tab_dir, paste0(run_prefix, "_contrast.csv"))
results_processed_csv <- file.path(tab_dir, paste0(run_prefix, "_radii_stability_processed.csv"))
pair_summary_csv <- file.path(tab_dir, paste0(run_prefix, "_radii_stability_pair_summary.csv"))
heatmap_data_csv <- file.path(tab_dir, paste0(run_prefix, "_radii_stability_heatmap_data.csv"))

heatmap_discrete_png <- file.path(fig_dir, paste0(run_prefix, "_radii_stability_heatmap_discrete.png"))
heatmap_discrete_pdf <- file.path(fig_dir, paste0(run_prefix, "_radii_stability_heatmap_discrete.pdf"))
heatmap_cont_png <- file.path(fig_dir, paste0(run_prefix, "_radii_stability_heatmap_continuous.png"))
heatmap_cont_pdf <- file.path(fig_dir, paste0(run_prefix, "_radii_stability_heatmap_continuous.pdf"))

# ---- Runtime tracking ----
runtime_tracker <- create_runtime_tracker(run_prefix)
append_runtime_metric <- runtime_tracker$append_runtime_metric
write_runtime_metrics <- runtime_tracker$write_runtime_metrics

# ---- Load preprocessed data ----
if (!file.exists(spe_list_path)) {
  stop("Missing preprocessed input: ", spe_list_path)
}
message("Reading preprocessed data: ", spe_list_path)
spe_list <- readRDS(spe_list_path)

sample_meta <- lapply(spe_list, function(spe) {
  data.frame(
    patient = get_unique(spe, patient_col),
    group_label = get_unique(spe, group_col),
    region = get_unique(spe, region_col),
    stringsAsFactors = FALSE
  )
})
sample_meta <- do.call(rbind, sample_meta)
sample_meta$sample <- sample_meta$region
names(spe_list) <- sample_meta$sample

design <- data.frame(
  sample = sample_meta$sample,
  group = sample_meta$group_label,
  stringsAsFactors = FALSE
)

# ---- Step 1: Prepare ----
t_prepare_start <- Sys.time()
prep <- panoramic_prepare(
  spe_list = spe_list,
  design = design,
  cell_type = cell_type_col,
  min_cells = min_cells,
  concavity = concavity,
  window = window,
  BPPARAM = BPPARAM
)
saveRDS(prep, prep_path)
t_prepare_end <- Sys.time()
append_runtime_metric(
  stage = "prepare",
  status = "computed",
  started_at = t_prepare_start,
  ended_at = t_prepare_end,
  cache_path = prep_path,
  params = list(min_cells = min_cells, concavity = concavity, window = window, workers = workers),
  obj = prep
)

# ---- Step 2: Spatial statistics + bootstrap ----
t_stats_start <- Sys.time()
se_stats_by_radius <- vector("list", length(radii_um))
for (i in seq_along(radii_um)) {
  radius_i <- radii_um[[i]]
  tile_size_i <- tile_size_by_radius[[i]]
  message("Running panoramic_spatialstats for radius ", radius_i, " (tile_size=", tile_size_i, ")")
  se_stats_by_radius[[i]] <- panoramic_spatialstats(
    prep = prep,
    pairs = "auto",
    radii_um = radius_i,
    stat = stat,
    nsim = nsim,
    correction = correction,
    seed = seed,
    boot = boot,
    tile_size = tile_size_i,
    BPPARAM = BPPARAM,
    verbose = FALSE
  )
}
se_stats <- do.call(SummarizedExperiment::rbind, se_stats_by_radius)
se_stats <- attach_sample_metadata(
  se = se_stats,
  sample_meta_df = sample_meta,
  patient_col = "patient",
  group_col = "group_label",
  region_col = "region",
  sample_meta_patient_col = "patient",
  sample_meta_group_col = "group_label"
)
saveRDS(se_stats, stats_path)
t_stats_end <- Sys.time()
append_runtime_metric(
  stage = "spatialstats",
  status = "computed",
  started_at = t_stats_start,
  ended_at = t_stats_end,
  cache_path = stats_path,
  params = list(
    stat = stat,
    radii_um = paste(radii_um, collapse = ","),
    nsim = nsim,
    correction = correction,
    seed = seed,
    boot = boot,
    tile_size_by_radius = paste0(radii_um, ":", tile_size_by_radius, collapse = ","),
    workers = workers
  ),
  obj = se_stats
)

# ---- Step 3: Meta-analysis ----
t_meta_start <- Sys.time()
se_meta <- panoramic_meta_mv(
  se = se_stats,
  patient_col = patient_col,
  group_col = group_col,
  sample_col = "sample",
  tau_structure = tau_structure,
  method = method_mv,
  group_tau2 = group_tau2,
  vi_floor = vi_floor,
  BPPARAM = BPPARAM
)
saveRDS(se_meta, meta_path)
t_meta_end <- Sys.time()
append_runtime_metric(
  stage = "meta_analysis",
  status = "computed",
  started_at = t_meta_start,
  ended_at = t_meta_end,
  cache_path = meta_path,
  params = list(
    method = method_mv,
    tau_structure = tau_structure,
    group_tau2 = group_tau2,
    vi_floor = vi_floor,
    workers = workers
  ),
  obj = se_meta
)

# ---- Extract tables ----
spatial_tbl <- extract_spatialstats_table(se_stats)
meta_tbl <- extract_meta_table(se_meta) |>
  dplyr::arrange(.data$p_diff)
contrast_tbl <- extract_contrast_table(se_meta)
if (nrow(contrast_tbl) == 0L) {
  stop("No contrast table available from se_meta; cannot build radii stability heatmap.")
}
contrast_tbl <- contrast_tbl |>
  dplyr::arrange(.data$radius_um, .data$p_diff)

readr::write_csv(spatial_tbl, spatial_csv)
readr::write_csv(meta_tbl, meta_csv)
readr::write_csv(contrast_tbl, contrast_csv)

# ---- Radii stability summaries ----
stability_out <- build_radii_stability_tables(
  contrast_tbl = contrast_tbl,
  radii_expected = radii_um,
  fdr_mode = fdr_mode,
  alpha_pair_include = alpha_pair_include,
  sig_pairs_only = sig_pairs_only,
  max_pairs = max_pairs
)

readr::write_csv(stability_out$processed, results_processed_csv)
readr::write_csv(stability_out$pair_summary, pair_summary_csv)
readr::write_csv(stability_out$plot_data, heatmap_data_csv)

# ---- Radii stability heatmaps ----
n_pairs_plot <- nrow(stability_out$pair_summary_plot)
heatmap_height <- max(heatmap_height_min, heatmap_height_per_pair * max(n_pairs_plot, 12L))

p_heatmap_discrete <- plot_radii_stability_discrete(
  plot_data = stability_out$plot_data,
  alpha_pair_include = alpha_pair_include,
  pair_axis_text_size = pair_axis_text_size
)
p_heatmap_cont <- plot_radii_stability_continuous(
  plot_data = stability_out$plot_data,
  alpha_pair_include = alpha_pair_include,
  pair_axis_text_size = pair_axis_text_size
)

ggsave(heatmap_discrete_png, p_heatmap_discrete, width = heatmap_width, height = heatmap_height, dpi = 300)
ggsave(heatmap_discrete_pdf, p_heatmap_discrete, width = heatmap_width, height = heatmap_height)
ggsave(heatmap_cont_png, p_heatmap_cont, width = heatmap_width, height = heatmap_height, dpi = 300)
ggsave(heatmap_cont_pdf, p_heatmap_cont, width = heatmap_width, height = heatmap_height)

if (isTRUE(show_plots)) {
  print(p_heatmap_discrete)
  print(p_heatmap_cont)
}

write_runtime_metrics(runtime_csv)

message("Finished CRC TMA radii stability analysis.")
message("Run prefix: ", run_prefix)
message("Radii tested: ", paste(radii_um, collapse = ", "))
message("FDR mode for heatmap significance tiers: within_radius")
message("Core tables: ", spatial_csv, " | ", meta_csv, " | ", contrast_csv)
message("Stability tables: ", results_processed_csv, " | ", pair_summary_csv, " | ", heatmap_data_csv)
message("Heatmap (discrete): ", heatmap_discrete_png, " | ", heatmap_discrete_pdf)
message("Heatmap (continuous): ", heatmap_cont_png, " | ", heatmap_cont_pdf)
message("Runtime metrics: ", runtime_csv)
