#!/usr/bin/env Rscript

# CRC TMA Script: Network Results Publication Figures
# Author: Jacob Chang
# Date: 2026-04-28
# Summary:
# - Rebuild CRC TMA network-distance plots from saved analysis outputs.
# - Export publication-sized individual panels and a compact composite.
# - Write PNG/PDF assets into the publication figures directory.

suppressPackageStartupMessages({
  library(config)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(panoramic)
  library(png)
  library(readr)
  library(S4Vectors)
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

require_file <- function(path, label = "file") {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", path)
  }
  path
}

require_first_existing <- function(paths, label = "file") {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) {
    stop(
      "Missing ", label, ". Checked:\n",
      paste0("  - ", paths, collapse = "\n")
    )
  }
  hit[[1]]
}

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

replace_cell_type_labels <- function(x, label_map) {
  if (!is.character(x)) return(x)
  out <- x
  for (k in names(label_map)) {
    pat <- paste0("(?<![A-Za-z0-9_])", k, "(?![A-Za-z0-9_])")
    out <- gsub(pat, label_map[[k]], out, perl = TRUE)
  }
  out
}

relabel_df_chars <- function(df, label_map) {
  if (!is.data.frame(df) || nrow(df) == 0L) return(df)
  char_cols <- names(df)[vapply(df, is.character, logical(1))]
  for (nm in char_cols) {
    df[[nm]] <- replace_cell_type_labels(df[[nm]], label_map)
  }
  df
}

relabel_ggplot <- function(p, label_map) {
  if (!inherits(p, "ggplot")) return(p)

  if (is.data.frame(p$data)) {
    p$data <- relabel_df_chars(p$data, label_map)
  }

  if (!is.null(p$labels)) {
    for (nm in names(p$labels)) {
      if (is.character(p$labels[[nm]])) {
        p$labels[[nm]] <- replace_cell_type_labels(p$labels[[nm]], label_map)
      }
    }
  }

  for (i in seq_along(p$layers)) {
    layer_data <- p$layers[[i]]$data
    if (is.data.frame(layer_data)) {
      p$layers[[i]]$data <- relabel_df_chars(layer_data, label_map)
    }
    if (!is.null(p$layers[[i]]$aes_params$label) &&
        is.character(p$layers[[i]]$aes_params$label)) {
      p$layers[[i]]$aes_params$label <- replace_cell_type_labels(
        p$layers[[i]]$aes_params$label,
        label_map
      )
    }
  }
  p
}

attach_sample_metadata <- function(se, sample_meta_df) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if (!"sample" %in% colnames(cd)) {
    cd$sample <- colnames(se)
  }
  cd$region <- cd$sample
  cd$patient <- sample_meta_df$patient[match(cd$sample, sample_meta_df$sample)]
  cd$group_label <- sample_meta_df$group_label[match(cd$sample, sample_meta_df$sample)]
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd)
  se
}

save_plot_formats <- function(plot_obj, out_stem, width_in, height_in, dpi = 300,
                              vector_mode = c("raster_wrapped", "native")) {
  vector_mode <- match.arg(vector_mode)

  png_path <- paste0(out_stem, ".png")
  pdf_path <- paste0(out_stem, ".pdf")

  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width_in,
    height = height_in,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  invisible(gc())

  if (identical(vector_mode, "native")) {
    ggsave(
      filename = pdf_path,
      plot = plot_obj,
      width = width_in,
      height = height_in,
      units = "in",
      bg = "white"
    )
    invisible(gc())
    return(invisible(NULL))
  }

  img <- png::readPNG(png_path)

  if (capabilities("cairo")) {
    grDevices::cairo_pdf(filename = pdf_path, width = width_in, height = height_in, bg = "white")
  } else {
    grDevices::pdf(file = pdf_path, width = width_in, height = height_in, bg = "white")
  }
  grid::grid.newpage()
  grid::grid.raster(
    img,
    x = 0.5,
    y = 0.5,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    interpolate = TRUE
  )
  grDevices::dev.off()

  rm(img)
  invisible(gc())
}

set_fixed_network_node_size <- function(p, node_size = 8) {
  if (!inherits(p, "ggplot")) return(p)
  node_size <- as.numeric(node_size)
  if (length(node_size) != 1L || is.na(node_size) || !is.finite(node_size) || node_size <= 0) {
    return(p)
  }

  for (i in seq_along(p$layers)) {
    geom_classes <- class(p$layers[[i]]$geom)
    if (any(geom_classes %in% c("GeomNodePoint", "GeomPoint"))) {
      if (!is.null(p$layers[[i]]$mapping$size)) {
        p$layers[[i]]$mapping$size <- NULL
      }
      p$layers[[i]]$aes_params$size <- node_size
    }
  }

  p
}

plot_representative_sample <- function(sample_id, panel_title, spe_list, sample_id_lookup,
                                       sample_compactness, cell_type_col = "cell_type",
                                       show_title = FALSE) {
  idx <- which(sample_id_lookup == sample_id)
  if (length(idx) == 0L) {
    return(ggplot() + theme_void())
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

  p <- ggplot(plot_df, aes(x = x, y = y, color = plot_type)) +
    geom_point(size = 0.2, alpha = 0.85) +
    scale_color_manual(
      values = c(
        "B Cell" = "#1F77B4",
        "CD8+ T Cell" = "#D62728",
        "CD4+ CD45RO+ T Cell" = "#2CA02C",
        "Stroma" = "#9467BD",
        "Other" = "grey85"
      ),
      drop = FALSE
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    coord_equal() +
    labs(
      title = if (show_title) panel_title else NULL,
      subtitle = sample_id,
      x = NULL,
      y = NULL,
      color = NULL
    ) +
    theme_void(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7, hjust = 0.5),
      plot.margin = margin(4, 4, 4, 4)
    )

  p
}

pick_ranked_sample <- function(representative_samples, group_name, rank_value) {
  hit <- representative_samples |>
    dplyr::filter(group == group_name, rank_within_group == rank_value)
  if (nrow(hit) == 0L) {
    stop(
      "Requested ", group_name, " representative sample at rank ", rank_value,
      " not available."
    )
  }
  as.character(hit$sample[1])
}

config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)

analysis_params <- list(
  network_run_prefix = "crc_tma_one_per_patient_local_comp_enrichment_r25_nsim100_block_tile62p5_mv",
  clr_rank = 2L,
  dii_rank = 3L,
  show_titles = FALSE,
  network_width_in = 4.25,
  network_height_in = 3.2,
  network_with_legend_width_in = 8.5,
  network_with_legend_height_in = 4.5,
  network_legend_width_in = 8.5,
  network_legend_height_in = 2.5,
  network_node_size = 8,
  compactness_width_in = 3.25,
  compactness_height_in = 3.0,
  representative_width_in = 2.6,
  representative_height_in = 2.85,
  representative_with_legend_width_in = 4.25,
  representative_with_legend_height_in = 3.6,
  representative_legend_width_in = 3.0,
  representative_legend_height_in = 0.45,
  composite_width_in = 9.0,
  composite_height_in = 3.4
)

cell_type_col <- "cell_type"
group_col <- "group_label"
patient_col <- "patient"
sample_col <- "spot"

network_run_prefix <- as.character(analysis_params$network_run_prefix)
fig_dir <- file.path(paths$output, "figures", "crc_tma", "one_sample_per_patient")
tab_dir <- file.path(paths$output, "tables", "crc_tma", "one_sample_per_patient")
out_dir <- file.path(fig_dir, "publication_panels", "network_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

network_tab_dir <- file.path(tab_dir, "network_distance_analysis")
sample_compactness_path <- require_file(
  file.path(network_tab_dir, "sample_compactness.csv"),
  label = "sample compactness CSV"
)
compactness_ttest_path <- require_file(
  file.path(network_tab_dir, "sample_compactness_ttest.csv"),
  label = "sample compactness t-test CSV"
)
representative_samples_path <- require_file(
  file.path(network_tab_dir, "representative_samples_nearest_group_mean.csv"),
  label = "representative sample ranking CSV"
)

interim_dir <- file.path(paths$data, "interim")
network_run_prefix_candidates <- unique(c(
  network_run_prefix,
  sub("block_tile625_mv$", "block_tile62p5_mv", network_run_prefix)
))

se_meta_path <- require_first_existing(
  file.path(interim_dir, paste0(network_run_prefix_candidates, "_se_meta.rds")),
  label = "se_meta RDS"
)
spe_list_path <- require_file(
  file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_one_per_patient.rds"),
  label = "CRC one-per-patient SPE list"
)

sample_compactness <- readr::read_csv(sample_compactness_path, show_col_types = FALSE)
compactness_ttest <- readr::read_csv(compactness_ttest_path, show_col_types = FALSE)
representative_samples <- readr::read_csv(representative_samples_path, show_col_types = FALSE)
se_meta <- readRDS(se_meta_path)
spe_list <- readRDS(spe_list_path)

if (length(spe_list) == 0L) {
  stop("SPE list is empty: ", spe_list_path)
}

required_cols <- c(cell_type_col, group_col, patient_col, sample_col)
missing_cols <- setdiff(required_cols, colnames(as.data.frame(SummarizedExperiment::colData(spe_list[[1]]))))
if (length(missing_cols) > 0L) {
  stop("Missing required colData columns in SPE list: ", paste(missing_cols, collapse = ", "))
}

sample_meta <- lapply(spe_list, function(spe) {
  data.frame(
    patient = get_unique(spe, patient_col),
    group_label = get_unique(spe, group_col),
    sample = get_unique(spe, sample_col),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

names(spe_list) <- sample_meta$sample
sample_id_lookup <- vapply(
  spe_list,
  function(spe) get_unique(spe, sample_col),
  FUN.VALUE = character(1)
)
se_meta <- attach_sample_metadata(se_meta, sample_meta)

compactness_plot <- ggplot(
  sample_compactness,
  aes(x = group, y = compactness_relative_mean, fill = group, color = group)
) +
  geom_violin(trim = FALSE, alpha = 0.2, color = NA) +
  geom_jitter(width = 0.08, height = 0, size = 0.9, alpha = 0.85) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.35, color = "black") +
  scale_fill_manual(values = c("CLR" = "#88c6ae", "DII" = "#88c6ae")) +
  scale_color_manual(values = c("CLR" = "#88c6ae", "DII" = "#88c6ae")) +
  labs(
    x = NULL,
    y = "Sample Compactness"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.margin = margin(4, 8, 4, 4)
  )

compactness_p_clr_less <- compactness_ttest$p_clr_less_dii[1]
p_label <- if (length(compactness_p_clr_less) == 1L && is.finite(compactness_p_clr_less)) {
  paste0("CLR < DII: p = ", format(compactness_p_clr_less, digits = 3))
} else {
  "CLR < DII: p = NA"
}

y_vals <- sample_compactness$compactness_relative_mean[
  is.finite(sample_compactness$compactness_relative_mean)
]
y_max <- if (length(y_vals) > 0L) max(y_vals) else 1
y_min <- if (length(y_vals) > 0L) min(y_vals) else 0
y_pad <- max(0.06 * (y_max - y_min), 0.03)

compactness_plot <- compactness_plot +
  annotate("text", x = 1.5, y = y_max + y_pad, label = p_label, size = 2.6) +
  coord_cartesian(ylim = c(0, y_max + 2 * y_pad), clip = "off")

network_plot <- panoramic::plot_spatial_network(
  se_diff = se_meta,
  fdr_threshold = 0.05,
  leiden_resolution = 1.2,
  z_sign = "negative",
  include_nonsig = TRUE,
  nonsig_max_fdr = 1.0,
  directed = FALSE,
  layout = "fr",
  sig_operator = "gt"
)

network_plot <- relabel_ggplot(network_plot, c(
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
  vasculature = "Vasculature"
) )
network_plot <- set_fixed_network_node_size(
  network_plot,
  node_size = analysis_params$network_node_size
)

network_plot_legend_source <- network_plot +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(4, 8, 4, 8),
    legend.margin = margin(2, 2, 2, 2),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

network_legend_plot <- ggdraw() +
  draw_grob(
    cowplot::get_legend(network_plot_legend_source),
    x = 0.5, y = 0.5,
    width = 0.98, height = 0.98
  )

network_plot <- network_plot +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

clr_rep_sample <- pick_ranked_sample(representative_samples, "CLR", analysis_params$clr_rank)
dii_rep_sample <- pick_ranked_sample(representative_samples, "DII", analysis_params$dii_rank)

clr_rep_plot <- plot_representative_sample(
  sample_id = clr_rep_sample,
  panel_title = "CLR representative sample",
  spe_list = spe_list,
  sample_id_lookup = sample_id_lookup,
  sample_compactness = sample_compactness,
  cell_type_col = cell_type_col,
  show_title = isTRUE(analysis_params$show_titles)
)

dii_rep_plot <- plot_representative_sample(
  sample_id = dii_rep_sample,
  panel_title = "DII representative sample",
  spe_list = spe_list,
  sample_id_lookup = sample_id_lookup,
  sample_compactness = sample_compactness,
  cell_type_col = cell_type_col,
  show_title = isTRUE(analysis_params$show_titles)
)

legend_only <- cowplot::get_legend(clr_rep_plot)
legend_plot <- ggdraw(legend_only)
clr_rep_plot_nolegend <- clr_rep_plot + theme(legend.position = "none")
dii_rep_plot_nolegend <- dii_rep_plot + theme(legend.position = "none")

rep_column <- plot_grid(
  clr_rep_plot_nolegend,
  dii_rep_plot_nolegend,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1)
)

composite_plot <- plot_grid(
  network_plot,
  compactness_plot,
  rep_column,
  nrow = 1,
  rel_widths = c(1.6, 1, 1.02),
  align = "h"
)

save_plot_formats(
  network_plot,
  file.path(out_dir, "plot_A_network_publication"),
  width_in = analysis_params$network_width_in,
  height_in = analysis_params$network_height_in,
  vector_mode = "native"
)

save_plot_formats(
  network_plot_legend_source,
  file.path(out_dir, "plot_A_network_with_legend_publication"),
  width_in = analysis_params$network_with_legend_width_in,
  height_in = analysis_params$network_with_legend_height_in,
  vector_mode = "native"
)

save_plot_formats(
  network_legend_plot,
  file.path(out_dir, "plot_A_network_legend_publication"),
  width_in = analysis_params$network_legend_width_in,
  height_in = analysis_params$network_legend_height_in,
  vector_mode = "native"
)

save_plot_formats(
  compactness_plot,
  file.path(out_dir, "plot_B_compactness_violin_publication"),
  width_in = analysis_params$compactness_width_in,
  height_in = analysis_params$compactness_height_in,
  vector_mode = "native"
)

save_plot_formats(
  clr_rep_plot_nolegend,
  file.path(out_dir, "plot_C_clr_representative_publication"),
  width_in = analysis_params$representative_width_in,
  height_in = analysis_params$representative_height_in,
  vector_mode = "raster_wrapped"
)

save_plot_formats(
  clr_rep_plot,
  file.path(out_dir, "plot_C_clr_representative_with_legend_publication"),
  width_in = analysis_params$representative_with_legend_width_in,
  height_in = analysis_params$representative_with_legend_height_in,
  vector_mode = "raster_wrapped"
)

save_plot_formats(
  dii_rep_plot_nolegend,
  file.path(out_dir, "plot_D_dii_representative_publication"),
  width_in = analysis_params$representative_width_in,
  height_in = analysis_params$representative_height_in,
  vector_mode = "raster_wrapped"
)

save_plot_formats(
  dii_rep_plot,
  file.path(out_dir, "plot_D_dii_representative_with_legend_publication"),
  width_in = analysis_params$representative_with_legend_width_in,
  height_in = analysis_params$representative_with_legend_height_in,
  vector_mode = "raster_wrapped"
)

save_plot_formats(
  legend_plot,
  file.path(out_dir, "plot_CD_representative_legend_publication"),
  width_in = analysis_params$representative_legend_width_in,
  height_in = analysis_params$representative_legend_height_in,
  vector_mode = "native"
)

save_plot_formats(
  composite_plot,
  file.path(out_dir, "compactness_network_composite_publication"),
  width_in = analysis_params$composite_width_in,
  height_in = analysis_params$composite_height_in,
  vector_mode = "raster_wrapped"
)

message("Saved publication network figures to: ", out_dir)
