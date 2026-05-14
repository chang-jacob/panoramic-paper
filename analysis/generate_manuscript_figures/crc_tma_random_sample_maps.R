#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(config)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(SpatialExperiment)
  library(SummarizedExperiment)
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

save_plot_formats <- function(plot_obj, stem, width, height, dpi = 300) {
  ggsave(
    filename = paste0(stem, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    filename = paste0(stem, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

add_missing_palette_colors <- function(base_palette, observed_levels) {
  observed_levels <- unique(as.character(observed_levels))
  missing_levels <- setdiff(observed_levels, names(base_palette))
  if (length(missing_levels) == 0L) {
    return(base_palette)
  }
  extra_cols <- scales::hue_pal()(length(missing_levels))
  names(extra_cols) <- missing_levels
  c(base_palette, extra_cols)
}

get_unique_value <- function(x, label) {
  vals <- unique(as.character(stats::na.omit(x)))
  if (length(vals) != 1L) {
    stop("Expected exactly one ", label, " value, found: ", paste(vals, collapse = ", "))
  }
  vals[[1]]
}

replace_cell_type_labels <- function(x, label_map) {
  out <- as.character(x)
  idx <- match(out, names(label_map))
  hit <- !is.na(idx)
  out[hit] <- unname(label_map[idx[hit]])
  out
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

extract_sample_metadata <- function(spe_list) {
  sample_names <- names(spe_list)
  if (is.null(sample_names)) sample_names <- rep(NA_character_, length(spe_list))

  dplyr::bind_rows(lapply(seq_along(spe_list), function(i) {
    spe <- spe_list[[i]]
    cd <- as.data.frame(SummarizedExperiment::colData(spe))
    sample_id <- if ("spot" %in% colnames(cd)) {
      get_unique_value(cd$spot, "spot")
    } else if (!is.na(sample_names[[i]]) && nzchar(sample_names[[i]])) {
      sample_names[[i]]
    } else {
      paste0("sample_", i)
    }

    data.frame(
      sample = sample_id,
      patient = get_unique_value(cd$patient, "patient"),
      group_label = get_unique_value(cd$group_label, "group_label"),
      n_cells = ncol(spe),
      stringsAsFactors = FALSE
    )
  }))
}

make_square_limits <- function(x, y) {
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  xmid <- mean(xr)
  ymid <- mean(yr)
  half_span <- max(diff(xr), diff(yr)) / 2
  if (!is.finite(half_span) || half_span <= 0) half_span <- 1
  list(
    xlim = c(xmid - half_span, xmid + half_span),
    ylim = c(ymid - half_span, ymid + half_span)
  )
}

plot_sample_map <- function(spe, palette_raw, label_map, point_size = 0.7, point_alpha = 0.95) {
  xy <- SpatialExperiment::spatialCoords(spe)
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  keep <- is.finite(xy[, 1]) & is.finite(xy[, 2])
  xy <- xy[keep, , drop = FALSE]
  cell_type_raw <- as.character(cd$cell_type[keep])

  palette_raw <- add_missing_palette_colors(palette_raw, sort(unique(cell_type_raw)))
  palette_display <- stats::setNames(
    unname(palette_raw),
    replace_cell_type_labels(names(palette_raw), label_map)
  )
  display_levels <- unique(replace_cell_type_labels(names(palette_raw), label_map))

  plot_df <- data.frame(
    x = xy[, 1],
    y = xy[, 2],
    cell_type = factor(
      replace_cell_type_labels(cell_type_raw, label_map),
      levels = display_levels
    ),
    stringsAsFactors = FALSE
  )

  lims <- make_square_limits(plot_df$x, plot_df$y)

  ggplot(plot_df, aes(x = x, y = y, color = cell_type)) +
    geom_point(size = point_size, alpha = point_alpha, shape = 16, stroke = 0) +
    scale_color_manual(values = palette_display, drop = FALSE) +
    coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

make_legend_plot <- function(palette_raw, label_map) {
  palette_display <- stats::setNames(
    unname(palette_raw),
    replace_cell_type_labels(names(palette_raw), label_map)
  )
  display_levels <- names(palette_display)

  legend_df <- data.frame(
    x = seq_along(display_levels),
    y = 1,
    cell_type = factor(display_levels, levels = display_levels),
    stringsAsFactors = FALSE
  )

  base_plot <- ggplot(legend_df, aes(x = x, y = y, color = cell_type)) +
    geom_point(size = 2.4, alpha = 1) +
    scale_color_manual(values = palette_display, drop = FALSE) +
    labs(color = NULL) +
    theme_void(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3, alpha = 1)))

  cowplot::ggdraw(cowplot::get_legend(base_plot))
}

project_root <- find_project_root()
config_file <- file.path(project_root, "config", "default.yml")
paths <- config::get("paths", file = config_file)

analysis_params <- list(
  selection_seed = 20260428L,
  n_samples_per_group = 3L,
  point_size = 0.5,
  point_alpha = 0.95,
  panel_width_in = 1,
  panel_height_in = 1,
  legend_width_in = 4.8,
  legend_height_in = 2.5
)

spe_list_path <- file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_one_per_patient.rds")
if (!file.exists(spe_list_path)) {
  stop("Missing one-per-patient CRC TMA SPE list: ", spe_list_path)
}

fig_dir <- file.path(paths$output, "figures", "crc_tma", "one_sample_per_patient", "random_sample_maps")
tab_dir <- file.path(paths$output, "tables", "crc_tma", "one_sample_per_patient", "random_sample_maps")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

celltype_palette <- c(
  "b_cells" = "#CF6B97",
  "cd4_t_cells_cd45ro" = "#56BAE9",
  "cd68_cd163_macrophages" = "#1171B8",
  "cd68_macrophages" = "#2596D3",
  "cd8_t_cells" = "#12B283",
  "granulocytes" = "#38A1B6",
  "plasma_cells" = "#6174A4",
  "tregs" = "#89B64A",
  "smooth_muscle" = "#F6C414",
  "stroma" = "#F3AE16",
  "vasculature" = "#F09C1C",
  "tumor_cells" = "#D44F1C",
  "undefined" = "#D3D3D3"
)

cell_type_label_map <- c(
  tumor_cells = "Tumor",
  smooth_muscle = "Smooth Muscle",
  granulocytes = "Granulocyte",
  cd8_t_cells = "CD8+ T Cell",
  stroma = "Stroma",
  plasma_cells = "Plasma Cell",
  cd68_cd163_macrophages = "CD68+ CD163+ Macrophage",
  cd4_t_cells_cd45ro = "CD4+ CD45RO+ T Cell",
  b_cells = "B Cell",
  tregs = "Treg",
  cd68_macrophages = "CD68+ Macrophage",
  vaculature = "Vasculature",
  vasculature = "Vasculature",
  undefined = "Undefined"
)

spe_list <- readRDS(spe_list_path)
sample_meta <- extract_sample_metadata(spe_list)

if (!all(c("CLR", "DII") %in% unique(sample_meta$group_label))) {
  stop("Expected CLR and DII groups in one-per-patient SPE list.")
}

group_counts <- sample_meta %>%
  dplyr::count(group_label, name = "n_samples_available")
missing_groups <- group_counts %>%
  dplyr::filter(group_label %in% c("CLR", "DII"), n_samples_available < analysis_params$n_samples_per_group)
if (nrow(missing_groups) > 0L) {
  stop(
    "Not enough one-per-patient samples available for requested random draw:\n",
    paste(
      paste0(missing_groups$group_label, ": ", missing_groups$n_samples_available),
      collapse = "\n"
    )
  )
}

set.seed(as.integer(analysis_params$selection_seed))
selected_samples <- sample_meta %>%
  dplyr::filter(group_label %in% c("CLR", "DII")) %>%
  dplyr::mutate(group_label = factor(group_label, levels = c("CLR", "DII"))) %>%
  dplyr::group_by(group_label) %>%
  dplyr::slice_sample(n = as.integer(analysis_params$n_samples_per_group)) %>%
  dplyr::mutate(group_rank = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group_label, group_rank)

sample_lookup <- selected_samples$sample
old_spe_names <- names(spe_list)
names(spe_list) <- vapply(seq_along(spe_list), function(i) {
  cd <- as.data.frame(SummarizedExperiment::colData(spe_list[[i]]))
  if ("spot" %in% colnames(cd)) {
    get_unique_value(cd$spot, "spot")
  } else {
    if (!is.null(old_spe_names) && nzchar(old_spe_names[i])) old_spe_names[i] else paste0("sample_", i)
  }
}, character(1))

palette_raw_full <- add_missing_palette_colors(
  celltype_palette,
  sort(unique(unlist(lapply(spe_list[sample_lookup], function(spe) {
    as.character(SummarizedExperiment::colData(spe)$cell_type)
  }))))
)
palette_key <- data.frame(
  cell_type_raw = names(palette_raw_full),
  cell_type_label = replace_cell_type_labels(names(palette_raw_full), cell_type_label_map),
  hex = unname(palette_raw_full),
  stringsAsFactors = FALSE
)
utils::write.csv(palette_key, file = file.path(tab_dir, "random_sample_map_palette_key.csv"), row.names = FALSE)

plot_list <- vector("list", nrow(selected_samples))
manifest_rows <- vector("list", nrow(selected_samples))

for (i in seq_len(nrow(selected_samples))) {
  sample_id <- as.character(selected_samples$sample[i])
  group_label <- as.character(selected_samples$group_label[i])
  group_rank <- as.integer(selected_samples$group_rank[i])
  patient_id <- as.character(selected_samples$patient[i])
  spe <- spe_list[[sample_id]]

  p <- plot_sample_map(
    spe = spe,
    palette_raw = palette_raw_full,
    label_map = cell_type_label_map,
    point_size = analysis_params$point_size,
    point_alpha = analysis_params$point_alpha
  )
  plot_list[[i]] <- p

  file_stub <- paste0(
    tolower(group_label),
    "_random_",
    sprintf("%02d", group_rank),
    "_",
    sanitize_filename(sample_id)
  )
  save_plot_formats(
    plot_obj = p,
    stem = file.path(fig_dir, file_stub),
    width = analysis_params$panel_width_in,
    height = analysis_params$panel_height_in
  )

  manifest_rows[[i]] <- data.frame(
    selection_seed = as.integer(analysis_params$selection_seed),
    group_label = group_label,
    group_rank = group_rank,
    patient = patient_id,
    sample = sample_id,
    n_cells = as.integer(selected_samples$n_cells[i]),
    file_stub = file_stub,
    stringsAsFactors = FALSE
  )
}

selection_manifest <- dplyr::bind_rows(manifest_rows)
utils::write.csv(
  selection_manifest,
  file = file.path(tab_dir, "random_sample_map_selection_manifest.csv"),
  row.names = FALSE
)

legend_plot <- make_legend_plot(
  palette_raw = palette_raw_full,
  label_map = cell_type_label_map
)
save_plot_formats(
  plot_obj = legend_plot,
  stem = file.path(fig_dir, "cell_type_legend"),
  width = analysis_params$legend_width_in,
  height = analysis_params$legend_height_in
)

combined_grid <- cowplot::plot_grid(
  plotlist = plot_list,
  ncol = as.integer(analysis_params$n_samples_per_group),
  align = "hv",
  axis = "tblr"
)
save_plot_formats(
  plot_obj = combined_grid,
  stem = file.path(fig_dir, "random_sample_maps_grid"),
  width = analysis_params$panel_width_in * analysis_params$n_samples_per_group,
  height = analysis_params$panel_height_in * 2
)

message("Saved random CRC TMA sample maps to: ", fig_dir)
message("Saved manifest/palette tables to: ", tab_dir)
print(selection_manifest)
