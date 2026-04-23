#!/usr/bin/env Rscript

# CRC TMA Script: Network Distance Analysis
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Compute nearest-neighbor distance and compactness metrics from CRC SPE samples.
# - Run CLR vs DII t-tests for pairwise distances and compactness outcomes.
# - Export network/representative sample plots and composite figure panels.

suppressPackageStartupMessages({
  library(config)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
  library(panoramic)
  library(RANN)
})

get_unique <- function(spe, col) {
  vals <- unique(as.character(SummarizedExperiment::colData(spe)[[col]]))
  vals <- vals[!is.na(vals)]
  if (length(vals) != 1L) {
    stop("Expected one unique value for colData[['", col, "']] per sample, got: ",
         paste(vals, collapse = ", "))
  }
  vals
}

extract_xy <- function(spe) {
  coords <- as.matrix(SpatialExperiment::spatialCoords(spe))
  if (ncol(coords) < 2L) {
    stop("spatialCoords(spe) must have at least two coordinate columns.")
  }
  out <- coords[, 1:2, drop = FALSE]
  colnames(out) <- c("x", "y")
  out
}

compute_nearest_dist <- function(source_coords, target_coords, same_type = FALSE) {
  n_source <- nrow(source_coords)
  n_target <- nrow(target_coords)
  if (n_source == 0L) return(numeric(0))
  if (n_target == 0L) return(rep(Inf, n_source))

  if (same_type && n_target < 2L) {
    return(rep(Inf, n_source))
  }

  if (same_type) {
    nn <- RANN::nn2(data = target_coords, query = source_coords, k = 2)
    return(as.numeric(nn$nn.dists[, 2]))
  }
  nn <- RANN::nn2(data = target_coords, query = source_coords, k = 1)
  as.numeric(nn$nn.dists[, 1])
}

compute_nearest_any_excluding_self <- function(all_coords, source_idx) {
  n_all <- nrow(all_coords)
  if (n_all <= 1L || length(source_idx) == 0L) {
    return(rep(Inf, length(source_idx)))
  }

  query <- all_coords[source_idx, , drop = FALSE]

  nn <- RANN::nn2(data = all_coords, query = query, k = 2)
  as.numeric(nn$nn.dists[, 2])
}

run_one_sided_ttests <- function(sample_means) {
  pairs <- unique(sample_means[, c("source_type", "target_type")])
  out <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    src <- as.character(pairs$source_type[i])
    tgt <- as.character(pairs$target_type[i])
    sub <- sample_means |>
      dplyr::filter(source_type == src, target_type == tgt) |>
      dplyr::filter(group %in% c("CLR", "DII"), is.finite(mean_distance_um))

    x <- sub$mean_distance_um[sub$group == "CLR"]
    y <- sub$mean_distance_um[sub$group == "DII"]

    mean_clr <- if (length(x) > 0L) mean(x, na.rm = TRUE) else NA_real_
    mean_dii <- if (length(y) > 0L) mean(y, na.rm = TRUE) else NA_real_
    delta <- mean_clr - mean_dii

    row <- data.frame(
      source_type = src,
      target_type = tgt,
      n_clr = length(x),
      n_dii = length(y),
      mean_clr = mean_clr,
      mean_dii = mean_dii,
      delta_clr_minus_dii = delta,
      alternative_observed = NA_character_,
      p_one_sided_observed = NA_real_,
      p_clr_greater_dii = NA_real_,
      p_clr_less_dii = NA_real_,
      statistic_t = NA_real_,
      df = NA_real_,
      stringsAsFactors = FALSE
    )

    if (length(x) >= 2L && length(y) >= 2L) {
      alt_obs <- if (is.finite(delta) && delta >= 0) "greater" else "less"
      tt_obs <- stats::t.test(x, y, alternative = alt_obs)
      tt_gt <- stats::t.test(x, y, alternative = "greater")
      tt_lt <- stats::t.test(x, y, alternative = "less")

      row$alternative_observed <- alt_obs
      row$p_one_sided_observed <- tt_obs$p.value
      row$p_clr_greater_dii <- tt_gt$p.value
      row$p_clr_less_dii <- tt_lt$p.value
      row$statistic_t <- unname(tt_obs$statistic)
      row$df <- unname(tt_obs$parameter)
    }

    out[[i]] <- row
  }
  dplyr::bind_rows(out)
}

run_compactness_ttest <- function(sample_compactness, measure_col = "compactness_mean_distance_um") {
  stopifnot(measure_col %in% colnames(sample_compactness))
  sub <- sample_compactness |>
    dplyr::filter(group %in% c("CLR", "DII"), is.finite(.data[[measure_col]]))

  x <- sub[[measure_col]][sub$group == "CLR"]
  y <- sub[[measure_col]][sub$group == "DII"]

  out <- data.frame(
    measure = measure_col,
    n_clr = length(x),
    n_dii = length(y),
    mean_clr = if (length(x) > 0L) mean(x, na.rm = TRUE) else NA_real_,
    mean_dii = if (length(y) > 0L) mean(y, na.rm = TRUE) else NA_real_,
    delta_clr_minus_dii = NA_real_,
    p_two_sided = NA_real_,
    p_clr_greater_dii = NA_real_,
    p_clr_less_dii = NA_real_,
    statistic_t = NA_real_,
    df = NA_real_,
    stringsAsFactors = FALSE
  )
  out$delta_clr_minus_dii <- out$mean_clr - out$mean_dii

  if (length(x) >= 2L && length(y) >= 2L) {
    tt2 <- stats::t.test(x, y, alternative = "two.sided")
    ttg <- stats::t.test(x, y, alternative = "greater")
    ttl <- stats::t.test(x, y, alternative = "less")
    out$p_two_sided <- tt2$p.value
    out$p_clr_greater_dii <- ttg$p.value
    out$p_clr_less_dii <- ttl$p.value
    out$statistic_t <- unname(tt2$statistic)
    out$df <- unname(tt2$parameter)
  }
  out
}

# ---- Config ----
config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)

# ---- Required columns/inputs ----
cell_type_col <- "cell_type"
group_col <- "group_label"
patient_col <- "patient"
sample_col <- "spot"

network_run_prefix <- "crc_tma_one_per_patient_local_comp_enrichment_r25_nsim100_block_tile625_mv"
se_meta_path <- file.path(paths$data, "interim", paste0(network_run_prefix, "_se_meta.rds"))
if (!file.exists(se_meta_path)) {
  stop("Missing se_meta RDS: ", se_meta_path)
}

# ---- Analysis parameters ----
# Cell types in compactness analysis.
target_cell_types <- c(
  "b_cells",
  "cd8_t_cells",
  "cd4_t_cells_cd45ro",
  "stroma"
)
max_dist_um <- 25
groups_keep <- c("CLR", "DII")
clr_rank_for_composite <- 2L
dii_rank_for_composite <- 3L
n_representative_samples <- max(clr_rank_for_composite, dii_rank_for_composite)

# ---- Input ----
spe_list_path <- file.path(
  paths$data, "processed", "crc_tma", "crc_tma_spe_list_one_per_patient.rds"
)
if (!file.exists(spe_list_path)) {
  stop("Missing SPE list input: ", spe_list_path)
}
spe_list <- readRDS(spe_list_path)
if (length(spe_list) == 0L) stop("SPE list is empty: ", spe_list_path)
required_cols <- c(cell_type_col, group_col, patient_col, sample_col)
missing_cols <- setdiff(required_cols, colnames(as.data.frame(SummarizedExperiment::colData(spe_list[[1]]))))
if (length(missing_cols) > 0L) {
  stop("Missing required colData columns in SPE list: ", paste(missing_cols, collapse = ", "))
}

# ---- Output ----
out_dir <- file.path(
  paths$output, "tables", "crc_tma", "one_sample_per_patient", "network_distance_analysis"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Distance analysis ----
rows <- list()
i_out <- 0L
sample_compactness_rows <- list()
i_comp <- 0L

for (i in seq_along(spe_list)) {
  spe <- spe_list[[i]]
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  xy <- extract_xy(spe)

  cell_type <- as.character(cd[[cell_type_col]])
  group <- get_unique(spe, group_col)
  if (!group %in% groups_keep) next
  patient <- get_unique(spe, patient_col)
  sample_id <- get_unique(spe, sample_col)

  keep <- which(
    !is.na(cell_type) &
      is.finite(xy[, 1]) &
      is.finite(xy[, 2])
  )
  if (length(keep) == 0L) next

  xy <- xy[keep, , drop = FALSE]
  cell_type <- cell_type[keep]

  idx_by_type <- lapply(target_cell_types, function(ct) which(cell_type == ct))
  names(idx_by_type) <- target_cell_types

  sample_cell_values <- numeric(0)
  sample_cell_values_relative <- numeric(0)

  for (src in target_cell_types) {
    idx_src <- idx_by_type[[src]]
    if (length(idx_src) == 0L) next

    dist_mat <- matrix(NA_real_, nrow = length(idx_src), ncol = length(target_cell_types))
    colnames(dist_mat) <- target_cell_types
    j <- 0L

    for (tgt in target_cell_types) {
      j <- j + 1L
      idx_tgt <- idx_by_type[[tgt]]
      d_raw <- compute_nearest_dist(
        source_coords = xy[idx_src, , drop = FALSE],
        target_coords = xy[idx_tgt, , drop = FALSE],
        same_type = identical(src, tgt)
      )
      d_cap <- pmin(d_raw, max_dist_um)
      dist_mat[, j] <- d_cap

      i_out <- i_out + 1L
      rows[[i_out]] <- data.frame(
        sample = sample_id,
        patient = patient,
        group = group,
        source_type = src,
        target_type = tgt,
        n_source_cells = length(d_cap),
        mean_distance_um = mean(d_cap, na.rm = TRUE),
        median_distance_um = stats::median(d_cap, na.rm = TRUE),
        sd_distance_um = stats::sd(d_cap, na.rm = TRUE),
        prop_capped_at_max = mean(!is.finite(d_raw) | d_raw >= max_dist_um),
        stringsAsFactors = FALSE
      )
    }

    # Per-cell compactness for this source type in this sample:
    # mean nearest distance across all target types for each source cell.
    cell_compactness_src <- rowMeans(dist_mat, na.rm = TRUE)
    sample_cell_values <- c(sample_cell_values, cell_compactness_src)

    # Density-adjusted compactness: relative to nearest neighbor distance
    # to any cell in the sample.
    d_any <- compute_nearest_any_excluding_self(all_coords = xy, source_idx = idx_src)
    cell_compactness_src_relative <- cell_compactness_src / pmax(d_any, 1e-8)
    cell_compactness_src_relative[!is.finite(cell_compactness_src_relative)] <- NA_real_
    sample_cell_values_relative <- c(
      sample_cell_values_relative,
      cell_compactness_src_relative[is.finite(cell_compactness_src_relative)]
    )
  }

  if (length(sample_cell_values) > 0L) {
    i_comp <- i_comp + 1L
    sample_compactness_rows[[i_comp]] <- data.frame(
      sample = sample_id,
      patient = patient,
      group = group,
      n_cells_used = length(sample_cell_values),
      compactness_mean_distance_um = mean(sample_cell_values, na.rm = TRUE),
      compactness_median_distance_um = stats::median(sample_cell_values, na.rm = TRUE),
      n_cells_used_relative = length(sample_cell_values_relative),
      compactness_relative_mean = if (length(sample_cell_values_relative) > 0L) mean(sample_cell_values_relative, na.rm = TRUE) else NA_real_,
      compactness_relative_median = if (length(sample_cell_values_relative) > 0L) stats::median(sample_cell_values_relative, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

sample_means <- dplyr::bind_rows(rows)
if (nrow(sample_means) == 0L) {
  stop("No distance summaries were produced. Check input SPE list and cell-type labels.")
}

group_summary <- sample_means |>
  dplyr::group_by(source_type, target_type, group) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    mean_of_sample_means = mean(mean_distance_um, na.rm = TRUE),
    median_of_sample_means = median(mean_distance_um, na.rm = TRUE),
    sd_of_sample_means = sd(mean_distance_um, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(source_type, target_type, group)

ttest_results <- run_one_sided_ttests(sample_means) |>
  dplyr::arrange(source_type, target_type)

sample_compactness <- dplyr::bind_rows(sample_compactness_rows)
if (nrow(sample_compactness) == 0L) {
  stop("No sample compactness values were produced.")
}
sample_compactness <- sample_compactness |>
  dplyr::mutate(
    group = factor(group, levels = c("CLR", "DII"))
  ) |>
  dplyr::arrange(group, sample)

compactness_ttest <- run_compactness_ttest(
  sample_compactness,
  measure_col = "compactness_relative_mean"
)

analysis_meta <- data.frame(
  max_distance_um = max_dist_um,
  groups_included = paste(groups_keep, collapse = ","),
  primary_compactness_measure = "compactness_relative_mean",
  source_spe_list = spe_list_path,
  target_cell_types = paste(target_cell_types, collapse = ","),
  stringsAsFactors = FALSE
)

readr::write_csv(sample_means, file.path(out_dir, "nearest_distance_sample_means.csv"))
readr::write_csv(group_summary, file.path(out_dir, "nearest_distance_group_summary.csv"))
readr::write_csv(ttest_results, file.path(out_dir, "nearest_distance_one_sided_ttests.csv"))
readr::write_csv(sample_compactness, file.path(out_dir, "sample_compactness.csv"))
readr::write_csv(compactness_ttest, file.path(out_dir, "sample_compactness_ttest.csv"))
readr::write_csv(analysis_meta, file.path(out_dir, "nearest_distance_analysis_metadata.csv"))

# ---- Representative samples nearest group means ----
group_means <- sample_compactness |>
  dplyr::filter(group %in% c("CLR", "DII"), is.finite(compactness_relative_mean)) |>
  dplyr::group_by(group) |>
  dplyr::summarise(group_mean = mean(compactness_relative_mean, na.rm = TRUE), .groups = "drop")

representative_samples <- sample_compactness |>
  dplyr::filter(group %in% c("CLR", "DII"), is.finite(compactness_relative_mean)) |>
  dplyr::left_join(group_means, by = "group") |>
  dplyr::mutate(abs_diff_from_group_mean = abs(compactness_relative_mean - group_mean)) |>
  dplyr::group_by(group) |>
  dplyr::arrange(abs_diff_from_group_mean, sample, .by_group = TRUE) |>
  dplyr::mutate(rank_within_group = dplyr::row_number()) |>
  dplyr::slice_head(n = n_representative_samples) |>
  dplyr::ungroup()

readr::write_csv(
  representative_samples,
  file.path(out_dir, "representative_samples_nearest_group_mean.csv")
)

# ---- Plot B: compactness violin + one-sided p-value annotation ----
compactness_plot <- ggplot(
  sample_compactness,
  aes(x = group, y = compactness_relative_mean, fill = group, color = group)
) +
  geom_violin(trim = FALSE, alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.35, color = "black") +
  geom_jitter(width = 0.08, height = 0, size = 2.2, alpha = 0.9) +
  scale_fill_manual(values = c("CLR" = "#88c6ae", "DII" = "#88c6ae")) +
  scale_color_manual(values = c("CLR" = "#88c6ae", "DII" = "#88c6ae")) +
  labs(
    x = NULL,
    y = "Sample Compactness (relative mean distance)",
    title = "CRC one-per-patient compactness by group"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

compactness_p_clr_less <- compactness_ttest$p_clr_less_dii[1]
p_label <- if (is.finite(compactness_p_clr_less)) {
  paste0("One-sided t-test (CLR < DII): p = ", format(compactness_p_clr_less, digits = 3))
} else {
  "One-sided t-test (CLR < DII): p = NA"
}

y_vals <- sample_compactness$compactness_relative_mean[
  is.finite(sample_compactness$compactness_relative_mean)
]
y_max <- if (length(y_vals) > 0L) max(y_vals) else 1
y_min <- if (length(y_vals) > 0L) min(y_vals) else 0
y_pad <- max(0.05 * (y_max - y_min), 0.02)

compactness_plot_annotated <- compactness_plot +
  annotate("text", x = 1.5, y = y_max + y_pad, label = p_label, size = 3.4) +
  coord_cartesian(ylim = c(y_min, y_max + 2 * y_pad), clip = "off")

# ---- Plot A: network ----
se_meta <- readRDS(se_meta_path)

network_plot <- panoramic::plot_spatial_network(
  se_diff = se_meta,
  fdr_threshold = 0.05,
  leiden_resolution = 1.2,
  z_sign = "negative",
  include_nonsig = TRUE,
  nonsig_max_fdr = 1.0,
  directed = FALSE,
  layout = "fr"
)

format_network_plot_legend <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  p +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.key.height = grid::unit(0.3, "cm"),
      legend.key.width = grid::unit(0.6, "cm"),
      legend.spacing.x = grid::unit(0.08, "cm")
    ) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      fill = guide_legend(nrow = 2, byrow = TRUE),
      shape = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 2, byrow = TRUE),
      size = guide_legend(nrow = 2, byrow = TRUE),
      alpha = guide_legend(nrow = 2, byrow = TRUE)
    )
}

network_plot <- format_network_plot_legend(network_plot)

sample_id_lookup <- vapply(
  spe_list,
  function(spe) get_unique(spe, sample_col),
  FUN.VALUE = character(1)
)

# ---- Plots C/D: representative CLR and DII sample maps ----
plot_representative_sample <- function(sample_id, panel_title) {
  idx <- which(sample_id_lookup == sample_id)
  if (length(idx) == 0L) {
    return(ggplot() + theme_void() + labs(title = paste0(panel_title, ": sample not found")))
  }

  spe <- spe_list[[idx[1]]]
  xy <- extract_xy(spe)
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  ct <- as.character(cd[[cell_type_col]])
  keep <- which(is.finite(xy[, 1]) & is.finite(xy[, 2]))
  xy <- xy[keep, , drop = FALSE]
  ct <- ct[keep]

  plot_type <- ifelse(
    ct == "b_cells", "B Cell",
    ifelse(
      ct == "cd8_t_cells", "CD8+ T Cell",
      ifelse(
        ct == "cd4_t_cells_cd45ro", "CD4+ CD45RO+ T Cell",
        ifelse(ct == "stroma", "Stroma", "Other")
      )
    )
  )

  plot_df <- data.frame(
    x = xy[, 1],
    y = xy[, 2],
    plot_type = factor(
      plot_type,
      levels = c("B Cell", "CD8+ T Cell", "CD4+ CD45RO+ T Cell", "Stroma", "Other")
    ),
    stringsAsFactors = FALSE
  )

  compactness_val <- sample_compactness |>
    dplyr::filter(sample == sample_id) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(compactness_relative_mean)

  subtitle <- paste0(
    sample_id,
    " | compactness_relative_mean = ",
    if (length(compactness_val) == 1L && is.finite(compactness_val)) {
      format(round(compactness_val, 4), nsmall = 4)
    } else {
      "NA"
    }
  )

  ggplot(plot_df, aes(x = x, y = y, color = plot_type)) +
    geom_point(size = 0.5, alpha = 0.9) +
    scale_color_manual(
      values = c(
        "B Cell" = "#1F77B4",
        "CD8+ T Cell" = "#D62728",
        "CD4+ CD45RO+ T Cell" = "#2CA02C",
        "Stroma" = "#9467BD",
        "Other" = "grey80"
      ),
      drop = FALSE
    ) +
    coord_equal() +
    labs(
      title = panel_title,
      subtitle = subtitle,
      x = "x",
      y = "y",
      color = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom")
}

pick_composite_sample <- function(group_name, rank_value) {
  hit <- representative_samples |>
    dplyr::filter(group == group_name, rank_within_group == rank_value)
  if (nrow(hit) == 0L) {
    stop(
      "Requested ", group_name, " representative sample at rank ", rank_value,
      " not available. Increase n_representative_samples or check available samples."
    )
  }
  as.character(hit$sample[1])
}

clr_rep_sample <- pick_composite_sample("CLR", clr_rank_for_composite)
dii_rep_sample <- pick_composite_sample("DII", dii_rank_for_composite)

clr_rep_plot <- plot_representative_sample(
  sample_id = clr_rep_sample,
  panel_title = "CLR representative sample"
)
dii_rep_plot <- plot_representative_sample(
  sample_id = dii_rep_sample,
  panel_title = "DII representative sample"
)

# Print/save top 3 nearest-to-mean samples for each group.
top_rep_plot_dir <- file.path(
  out_dir,
  paste0("top", n_representative_samples, "_nearest_group_mean_sample_maps")
)
dir.create(top_rep_plot_dir, recursive = TRUE, showWarnings = FALSE)

message("Top nearest-to-mean samples by group:")
print(
  representative_samples |>
    dplyr::select(group, rank_within_group, sample, compactness_relative_mean, group_mean, abs_diff_from_group_mean)
)

for (i in seq_len(nrow(representative_samples))) {
  sid <- as.character(representative_samples$sample[i])
  grp <- as.character(representative_samples$group[i])
  rk <- representative_samples$rank_within_group[i]
  p_top <- plot_representative_sample(
    sample_id = sid,
    panel_title = paste0(grp, " sample nearest group mean #", rk)
  )
  print(p_top)
  stem <- file.path(top_rep_plot_dir, paste0("group_", grp, "_rank_", rk, "_", sid))
  ggsave(paste0(stem, ".png"), p_top, width = 4.5, height = 4.2, units = "in", dpi = 300)
  ggsave(paste0(stem, ".pdf"), p_top, width = 4.5, height = 4.2, units = "in")
}

# Print the four panel plots individually.
print(network_plot)
print(compactness_plot_annotated)
print(clr_rep_plot)
print(dii_rep_plot)

# Save the four panel plots individually.
ggsave(
  filename = file.path(out_dir, "plot_A_network.png"),
  plot = network_plot,
  width = 4.5,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "plot_A_network.pdf"),
  plot = network_plot,
  width = 4.5,
  height = 4,
  units = "in"
)
ggsave(
  filename = file.path(out_dir, "plot_B_compactness_violin.png"),
  plot = compactness_plot_annotated,
  width = 4.5,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "plot_B_compactness_violin.pdf"),
  plot = compactness_plot_annotated,
  width = 4.5,
  height = 4,
  units = "in"
)
ggsave(
  filename = file.path(out_dir, "plot_C_clr_representative.png"),
  plot = clr_rep_plot,
  width = 4.5,
  height = 4.2,
  units = "in",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "plot_C_clr_representative.pdf"),
  plot = clr_rep_plot,
  width = 4.5,
  height = 4.2,
  units = "in"
)
ggsave(
  filename = file.path(out_dir, "plot_D_dii_representative.png"),
  plot = dii_rep_plot,
  width = 4.5,
  height = 4.2,
  units = "in",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "plot_D_dii_representative.pdf"),
  plot = dii_rep_plot,
  width = 4.5,
  height = 4.2,
  units = "in"
)

# ---- Four-panel composite ----
composite_4panel <- cowplot::plot_grid(
  network_plot,
  compactness_plot_annotated,
  clr_rep_plot,
  dii_rep_plot,
  ncol = 2,
  labels = c("A", "B", "C (i)", "C (ii)")
)

print(composite_4panel)

ggsave(
  filename = file.path(out_dir, "compactness_network_composite_4panel.png"),
  plot = composite_4panel,
  width = 6.5,
  height = 5,
  units = "in",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "compactness_network_composite_4panel.pdf"),
  plot = composite_4panel,
  width = 6.5,
  height = 5,
  units = "in"
)

# ---- Alternative composite: one-row layout ----
# A: network (legend at bottom), B: violin, C column: C(i) and C(ii) stacked.
network_plot_one_row <- network_plot

c_column <- cowplot::plot_grid(
  clr_rep_plot + theme(legend.position = "none"),
  dii_rep_plot + theme(legend.position = "none"),
  ncol = 1,
  labels = c("C (i)", "C (ii)")
)

composite_1row <- cowplot::plot_grid(
  network_plot_one_row,
  compactness_plot_annotated,
  c_column,
  nrow = 1,
  labels = c("A", "B", ""),
  rel_widths = c(1.05, 1, 1)
)

print(composite_1row)

ggsave(
  filename = file.path(out_dir, "compactness_network_composite_1row.png"),
  plot = composite_1row,
  width = 9.5,
  height = 3.8,
  units = "in",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "compactness_network_composite_1row.pdf"),
  plot = composite_1row,
  width = 9.5,
  height = 3.8,
  units = "in"
)

message("Finished network-distance analysis and four-panel plotting.")
message("Output directory: ", out_dir)
