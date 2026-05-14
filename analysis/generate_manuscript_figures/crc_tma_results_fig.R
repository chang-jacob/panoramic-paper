#!/usr/bin/env Rscript

# CRC TMA Script: Publication Figure Assembly
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Assemble CRC TMA manuscript panels from explicit PANORAMIC output files.
# - Re-export clean-label volcano, representative, forest, and network panels.
# - Write publication-ready combined and individual figure assets.

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

call_with_supported_args <- function(fun, args) {
  keep <- intersect(names(args), names(formals(fun)))
  do.call(fun, args[keep])
}

pick_forest_from_index <- function(index_csv, ct1_value, ct2_value, variant = "full") {
  idx <- readr::read_csv(index_csv, show_col_types = FALSE)
  req <- c("ct1", "ct2", "variant", "status", "png")
  missing_req <- setdiff(req, colnames(idx))
  if (length(missing_req) > 0L) {
    stop("Forest index missing required columns: ", paste(missing_req, collapse = ", "))
  }
  keep <- idx[idx$ct1 == ct1_value & idx$ct2 == ct2_value & idx$variant == variant & idx$status == "ok", , drop = FALSE]
  if (nrow(keep) == 0L) {
    stop("No forest row found for ct1=", ct1_value, ", ct2=", ct2_value, ", variant=", variant, ".")
  }
  require_file(as.character(keep$png[1]), label = "forest PNG from index")
}

make_panel_from_png <- function(path, target_ar, fit = c("cover", "contain")) {
  fit <- match.arg(fit)
  img <- png::readPNG(path)
  if (length(dim(img)) < 3L) {
    img <- array(rep(img, 3L), dim = c(dim(img), 3L))
  }

  # Trim near-white border from source panel to reduce excess margins.
  rgb <- img[, , 1:3, drop = FALSE]
  alpha <- if (dim(img)[3] >= 4L) img[, , 4] else matrix(1, nrow = dim(img)[1], ncol = dim(img)[2])
  is_white <- alpha <= 0.01 | (rgb[, , 1] > 0.985 & rgb[, , 2] > 0.985 & rgb[, , 3] > 0.985)
  content_mask <- !is_white
  rows <- which(rowSums(content_mask) > 0L)
  cols <- which(colSums(content_mask) > 0L)
  if (length(rows) > 0L && length(cols) > 0L) {
    pad <- 2L
    r1 <- max(1L, min(rows) - pad)
    r2 <- min(dim(img)[1], max(rows) + pad)
    c1 <- max(1L, min(cols) - pad)
    c2 <- min(dim(img)[2], max(cols) + pad)
    img <- img[r1:r2, c1:c2, , drop = FALSE]
  }

  h <- dim(img)[1]
  w <- dim(img)[2]
  img_ar <- h / w

  if (fit == "cover") {
    # Fill panel completely (may crop).
    if (img_ar > target_ar) { # image is taller than target
      draw_w <- 1
      draw_h <- img_ar / target_ar
      x0 <- 0
      y0 <- (1 - draw_h) / 2
    } else {
      draw_h <- 1
      draw_w <- target_ar / img_ar
      x0 <- (1 - draw_w) / 2
      y0 <- 0
    }
  } else {
    # Keep image fully inside panel (no overflow).
    if (img_ar > target_ar) { # image is taller than target
      draw_h <- 1
      draw_w <- target_ar / img_ar
      x0 <- (1 - draw_w) / 2
      y0 <- 0
    } else {
      draw_w <- 1
      draw_h <- img_ar / target_ar
      x0 <- 0
      y0 <- (1 - draw_h) / 2
    }
  }

  ggdraw() +
    draw_grob(
      rasterGrob(img, interpolate = TRUE),
      x = x0, y = y0, width = draw_w, height = draw_h
    )
}

config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)

analysis_params <- list(
  run_prefix = "crc_tma_one_per_patient_local_comp_enrichment_r25_nsim100_block_tile625_mv",
  forest_variant = "full"
)

fig_dir <- file.path(paths$output, "figures", "crc_tma", "one_sample_per_patient")
tab_dir <- file.path(paths$output, "tables", "crc_tma", "one_sample_per_patient")
out_dir <- file.path(fig_dir, "publication_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_prefix <- as.character(analysis_params$run_prefix)
if (!nzchar(run_prefix)) {
  stop("analysis_params$run_prefix must be a non-empty string.")
}
forest_variant <- as.character(analysis_params$forest_variant)
if (!forest_variant %in% c("full", "minimal")) {
  stop("analysis_params$forest_variant must be 'full' or 'minimal'.")
}
forest_dir <- file.path(fig_dir, paste0(run_prefix, "_forest"))
forest_index_csv <- file.path(forest_dir, "forest_plot_index.csv")
require_file(forest_index_csv, label = "forest plot index CSV")

volcano_path <- require_file(
  file.path(fig_dir, paste0(run_prefix, "_volcano_fdr_diff.png")),
  label = "volcano panel PNG"
)

rep_b_cd8_path <- require_file(
  file.path(fig_dir, paste0(run_prefix, "_representative_samples_fdr_b_cells_to_cd8_t_cells_r25um.png")),
  label = "representative panel (b_cells -> cd8_t_cells)"
)

rep_mac_cd8_path <- require_file(
  file.path(fig_dir, paste0(run_prefix, "_representative_samples_fdr_cd68_cd163_macrophages_to_cd8_t_cells_r25um.png")),
  label = "representative panel (cd68_cd163_macrophages -> cd8_t_cells)"
)

forest_b_cd8_path <- pick_forest_from_index(
  index_csv = forest_index_csv,
  ct1_value = "cd8_t_cells",
  ct2_value = "b_cells",
  variant = forest_variant
)

forest_mac_cd8_path <- pick_forest_from_index(
  index_csv = forest_index_csv,
  ct1_value = "cd8_t_cells",
  ct2_value = "cd68_cd163_macrophages",
  variant = forest_variant
)

network_path <- require_file(
  file.path(fig_dir, paste0(run_prefix, "_network_negative_with_nonsig.png")),
  label = "network panel (CLR up-localization)"
)

# Figure geometry for 9 x 6.5 in output.
fig_width <- 9
fig_height <- 6.5
bottom_h <- 1.625
top_h <- (fig_height - bottom_h) / 2
mid_h <- top_h
row_heights <- c(top = top_h, middle = mid_h, bottom = bottom_h)
half_w <- fig_width / 2

# Make top and middle panels share the same aspect ratio so widths match visually.
top_ar <- top_h / half_w
rep_ar <- (top_h / 2) / half_w
mid_ar <- mid_h / half_w
bottom_panel_w <- 2
net_ar <- bottom_h / bottom_panel_w

volcano_panel <- make_panel_from_png(volcano_path, target_ar = top_ar)
rep_b_cd8_panel <- make_panel_from_png(rep_b_cd8_path, target_ar = rep_ar)
rep_mac_cd8_panel <- make_panel_from_png(rep_mac_cd8_path, target_ar = rep_ar)
forest_b_cd8_panel <- make_panel_from_png(forest_b_cd8_path, target_ar = mid_ar)
forest_mac_cd8_panel <- make_panel_from_png(forest_mac_cd8_path, target_ar = mid_ar)
network_panel <- make_panel_from_png(network_path, target_ar = net_ar, fit = "contain")

make_dummy_panel <- function(label) {
  ggdraw() +
    draw_label(label, size = 10, color = "grey45", fontface = "plain")
}

dummy_1 <- make_dummy_panel("Dummy 1")
dummy_2 <- make_dummy_panel("Dummy 2")
dummy_3 <- make_dummy_panel("Dummy 3")

top_right <- plot_grid(rep_b_cd8_panel, rep_mac_cd8_panel, ncol = 1, rel_heights = c(1, 1))
top_row <- plot_grid(volcano_panel, top_right, ncol = 2, rel_widths = c(1, 1))
middle_row <- plot_grid(forest_b_cd8_panel, forest_mac_cd8_panel, ncol = 2, rel_widths = c(1, 1))
bottom_row <- plot_grid(
  ggdraw(),
  network_panel,
  dummy_1,
  dummy_2,
  dummy_3,
  ggdraw(),
  ncol = 6,
  rel_widths = c(0.5, 2, 2, 2, 2, 0.5)
)

final_panel <- plot_grid(
  top_row,
  middle_row,
  bottom_row,
  ncol = 1,
  rel_heights = unname(row_heights)
)

out_prefix <- file.path(out_dir, paste0(run_prefix, "_crc_tma_results_panel_9x6p5in"))
message("Run prefix: ", run_prefix)
message("Panel output path stem (not written): ", out_prefix)
message("Selected inputs:")
message("  volcano: ", volcano_path)
message("  rep b->cd8: ", rep_b_cd8_path)
message("  rep mac->cd8: ", rep_mac_cd8_path)
message("  forest b->cd8: ", forest_b_cd8_path)
message("  forest mac->cd8: ", forest_mac_cd8_path)
message("  network: ", network_path)

# Reference only: keep the cowplot composition code above for layout planning.
# print(final_panel)

# ---- Cell type label mapping for individual exports ----
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
  vasculature = "Vasculature"
)

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

set_text_layer_size_pt <- function(p, text_pt = 10) {
  if (!inherits(p, "ggplot")) return(p)
  text_pt <- as.numeric(text_pt)
  if (length(text_pt) != 1L || is.na(text_pt) || !is.finite(text_pt) || text_pt <= 0) return(p)
  text_mm <- text_pt / ggplot2::.pt

  for (i in seq_along(p$layers)) {
    geom_classes <- class(p$layers[[i]]$geom)
    if (any(geom_classes %in% c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel"))) {
      p$layers[[i]]$aes_params$size <- text_mm
    }
  }
  p
}

drop_text_layers <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  keep <- vapply(
    p$layers,
    function(layer) {
      geom_classes <- class(layer$geom)
      !any(geom_classes %in% c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel"))
    },
    logical(1)
  )
  p$layers <- p$layers[keep]
  p
}

apply_uniform_title_size <- function(p, title_pt = 8, subtitle_pt = 7) {
  if (!inherits(p, "ggplot")) return(p)
  p + theme(
    plot.title = element_text(size = title_pt, face = "bold"),
    plot.subtitle = element_text(size = subtitle_pt)
  )
}

get_unique <- function(spe, col) {
  vals <- unique(as.character(SummarizedExperiment::colData(spe)[[col]]))
  vals <- vals[!is.na(vals)]
  if (length(vals) != 1L) {
    stop("Expected one unique value for colData[[", col, "]], got: ", paste(vals, collapse = ", "))
  }
  vals
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

choose_representative_sample <- function(values, target) {
  ok <- is.finite(values)
  if (!any(ok)) stop("No finite sample-level values available for representative selection.")
  if (!is.finite(target)) stop("Representative target value must be finite.")
  idx <- which(ok)
  idx[which.min(abs(values[idx] - target))]
}

make_feature_key <- function(df) {
  paste(df$ct1, df$ct2, df$radius_um, df$stat, sep = "|")
}

plot_representative_single <- function(se_stats, se_meta, spe_list,
                                       ct1, ct2, radius_um = 25,
                                       stat = "local_comp_enrichment",
                                       group_col = "group_label",
                                       cell_type_col = "cell_type",
                                       label_map = cell_type_label_map) {
  rd_stats <- as.data.frame(SummarizedExperiment::rowData(se_stats))
  rd_meta <- as.data.frame(SummarizedExperiment::rowData(se_meta))
  cd_stats <- as.data.frame(SummarizedExperiment::colData(se_stats))
  yi <- SummarizedExperiment::assay(se_stats, "yi")

  key_stats <- setNames(seq_len(nrow(rd_stats)), make_feature_key(rd_stats))
  key_meta <- setNames(seq_len(nrow(rd_meta)), make_feature_key(rd_meta))
  key <- paste(ct1, ct2, radius_um, stat, sep = "|")
  if (!key %in% names(key_stats) || !key %in% names(key_meta)) {
    stop("Representative feature not found for key: ", key)
  }

  row_stat <- key_stats[[key]]
  row_meta <- key_meta[[key]]
  y_row <- yi[row_stat, ]

  groups <- unique(as.character(cd_stats[[group_col]]))
  groups <- groups[!is.na(groups)]
  groups <- c(intersect(c("CLR", "DII"), groups), setdiff(groups, c("CLR", "DII")))
  groups <- unique(groups)
  if (length(groups) < 2L) {
    stop("Expected at least two groups in se_stats colData[['", group_col, "']].")
  }

  rep_list <- list()
  out_i <- 0L
  for (g in groups[1:2]) {
    idx_group <- which(as.character(cd_stats[[group_col]]) == g)
    if (length(idx_group) == 0L) next
    mu_col <- paste0(make.names(g), "_mu_hat")
    if (!mu_col %in% colnames(rd_meta)) {
      stop("Missing expected group mean column in se_meta rowData: ", mu_col)
    }
    mu_target <- rd_meta[row_meta, mu_col]
    pick_local <- choose_representative_sample(y_row[idx_group], mu_target)
    pick_col <- idx_group[pick_local]

    sample_id <- if ("sample" %in% colnames(cd_stats)) {
      as.character(cd_stats$sample[pick_col])
    } else {
      as.character(colnames(yi)[pick_col])
    }
    if (!sample_id %in% names(spe_list)) {
      stop("Representative sample id not found in SPE list: ", sample_id)
    }
    spe <- spe_list[[sample_id]]

    coords <- as.data.frame(SpatialExperiment::spatialCoords(spe))
    if (ncol(coords) < 2L) {
      stop("spatialCoords(spe) must contain at least two columns.")
    }
    colnames(coords)[1:2] <- c("x", "y")
    ct_vals <- as.character(SummarizedExperiment::colData(spe)[[cell_type_col]])

    role <- if (identical(ct1, ct2)) {
      ifelse(ct_vals == ct2, ct2, "Other")
    } else {
      ifelse(ct_vals == ct2, ct2, ifelse(ct_vals == ct1, ct1, "Other"))
    }

    out_i <- out_i + 1L
    rep_list[[out_i]] <- data.frame(
      x = coords$x,
      y = coords$y,
      role = role,
      group = g,
      sample = sample_id,
      panel = paste0(g, " | ", sample_id),
      stringsAsFactors = FALSE
    )
  }

  if (length(rep_list) == 0L) {
    stop("No representative samples were selected for ", ct2, " -> ", ct1)
  }

  rep_df <- dplyr::bind_rows(rep_list)
  rep_df$panel <- factor(rep_df$panel, levels = unique(rep_df$panel))

  source_label <- replace_cell_type_labels(ct2, label_map)
  target_label <- replace_cell_type_labels(ct1, label_map)
  rep_df$role <- replace_cell_type_labels(rep_df$role, label_map)

  color_levels <- if (identical(ct1, ct2)) c("Other", source_label) else c("Other", source_label, target_label)
  color_map <- if (identical(ct1, ct2)) {
    stats::setNames(c("grey80", "#1F77B4"), color_levels)
  } else {
    stats::setNames(c("grey80", "#1F77B4", "#FF7F0E"), color_levels)
  }

  ggplot(rep_df, aes(x = x, y = y, color = role)) +
    geom_point(size = 0.55, alpha = 0.9) +
    scale_color_manual(values = color_map, breaks = color_levels) +
    scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    scale_y_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    coord_cartesian(clip = "off") +
    facet_wrap(~panel, nrow = 1, scales = "free") +
    labs(
      title = paste0(source_label, " -> ", target_label),
      subtitle = paste0("Representative samples at r=", radius_um, " um"),
      x = "x",
      y = "y",
      color = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      strip.text = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8),
      aspect.ratio = 1,
      panel.spacing = grid::unit(0.15, "lines"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_plot_formats <- function(plot_obj, out_stem, width_in, height_in, dpi = 300,
                              vector_mode = c("raster_wrapped", "native")) {
  vector_mode <- match.arg(vector_mode)

  png_path <- paste0(out_stem, ".png")
  pdf_path <- paste0(out_stem, ".pdf")
  svg_path <- paste0(out_stem, ".svg")

  # PNG is always generated directly from ggplot.
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
    # Native vector export can be large/unstable for point-heavy plots.
    ggsave(
      filename = pdf_path,
      plot = plot_obj,
      width = width_in,
      height = height_in,
      units = "in",
      bg = "white"
    )
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggsave(
        filename = svg_path,
        plot = plot_obj,
        width = width_in,
        height = height_in,
        units = "in",
        bg = "white",
        device = svglite::svglite
      )
    } else {
      grDevices::svg(filename = svg_path, width = width_in, height = height_in, bg = "white")
      print(plot_obj)
      grDevices::dev.off()
    }
    invisible(gc())
    return(invisible(NULL))
  }

  # Stable mode: wrap rasterized PNG into PDF/SVG containers.
  img <- png::readPNG(png_path)

  if (capabilities("cairo")) {
    grDevices::cairo_pdf(filename = pdf_path, width = width_in, height = height_in, bg = "white")
  } else {
    grDevices::pdf(file = pdf_path, width = width_in, height = height_in, bg = "white")
  }
  grid::grid.newpage()
  grid::grid.raster(img, x = 0.5, y = 0.5, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)
  grDevices::dev.off()

  grDevices::svg(filename = svg_path, width = width_in, height = height_in, bg = "white")
  grid::grid.newpage()
  grid::grid.raster(img, x = 0.5, y = 0.5, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)
  grDevices::dev.off()

  rm(img)
  invisible(gc())
}

# ---- Build source objects needed for clean-label re-export ----
interim_dir <- file.path(paths$data, "interim")
se_stats_path <- file.path(interim_dir, paste0(run_prefix, "_se_stats.rds"))
se_meta_path <- file.path(interim_dir, paste0(run_prefix, "_se_meta.rds"))
spe_list_path <- file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_one_per_patient.rds")

if (!file.exists(se_stats_path) || !file.exists(se_meta_path)) {
  stop("Missing interim se_stats/se_meta RDS for run prefix: ", run_prefix)
}
if (!file.exists(spe_list_path)) {
  stop("Missing SPE list: ", spe_list_path)
}

se_stats <- readRDS(se_stats_path)
se_meta <- readRDS(se_meta_path)
spe_list <- readRDS(spe_list_path)

patient_col <- "patient"
group_col_spe <- "group_label"
region_col_spe <- "spot"
cell_type_col <- "cell_type"
required_spe_cols <- c(patient_col, group_col_spe, region_col_spe, cell_type_col)
missing_spe_cols <- setdiff(required_spe_cols, colnames(as.data.frame(SummarizedExperiment::colData(spe_list[[1]]))))
if (length(missing_spe_cols) > 0L) {
  stop("Missing required colData columns in one_per_patient SPE list: ", paste(missing_spe_cols, collapse = ", "))
}

sample_meta <- lapply(spe_list, function(spe) {
  data.frame(
    patient = get_unique(spe, patient_col),
    group_label = get_unique(spe, group_col_spe),
    region = get_unique(spe, region_col_spe),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()
sample_meta$sample <- sample_meta$region
names(spe_list) <- sample_meta$sample

se_stats <- attach_sample_metadata(se_stats, sample_meta)
group_col_meta <- "group_label"
if (!group_col_meta %in% colnames(as.data.frame(SummarizedExperiment::colData(se_meta)))) {
  stop("Missing required group column in se_meta colData: ", group_col_meta)
}

exports_dir <- file.path(out_dir, "individual_exports")
dir.create(exports_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Generate plots with clean labels ----
volcano_plot <- call_with_supported_args(
  panoramic::plot_volcano,
  list(
    se_diff = se_meta,
    sig_col = "fdr_diff",
    alpha = 0.05,
    x_scale = "log2fc",
    label_top = 12,
    label_text_pt = 6,
    title = "CRC TMA colocalization volcano"
  )
)
volcano_plot <- relabel_ggplot(volcano_plot, cell_type_label_map)
volcano_plot <- set_text_layer_size_pt(volcano_plot, text_pt = 6)
volcano_plot <- volcano_plot +
  labs(color = NULL) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.text = element_text(size = 7)
  )
volcano_plot_no_legend <- volcano_plot + theme(legend.position = "none")
volcano_plot_no_annotations <- drop_text_layers(volcano_plot)
volcano_plot_no_annotations_no_legend <- volcano_plot_no_annotations + theme(legend.position = "none")

rep_b_cd8_plot <- plot_representative_single(
  se_stats = se_stats,
  se_meta = se_meta,
  spe_list = spe_list,
  ct1 = "cd8_t_cells",
  ct2 = "b_cells",
  radius_um = 25,
  group_col = "group_label",
  cell_type_col = cell_type_col,
  label_map = cell_type_label_map
)

rep_mac_cd8_plot <- plot_representative_single(
  se_stats = se_stats,
  se_meta = se_meta,
  spe_list = spe_list,
  ct1 = "cd8_t_cells",
  ct2 = "cd68_cd163_macrophages",
  radius_um = 25,
  group_col = "group_label",
  cell_type_col = cell_type_col,
  label_map = cell_type_label_map
)

forest_b_cd8_plot <- call_with_supported_args(
  panoramic::plot_forest,
  list(
    se_meta = se_meta,
    ct1 = "cd8_t_cells",
    ct2 = "b_cells",
    radius_um = 25,
    group_col = group_col_meta,
    show_est_se = TRUE,
    show_ci = TRUE
  )
)
forest_b_cd8_plot <- relabel_ggplot(forest_b_cd8_plot, cell_type_label_map)
forest_b_cd8_plot_no_stats <- call_with_supported_args(
  panoramic::plot_forest,
  list(
    se_meta = se_meta,
    ct1 = "cd8_t_cells",
    ct2 = "b_cells",
    radius_um = 25,
    group_col = group_col_meta,
    show_est_se = FALSE,
    show_ci = FALSE
  )
)
forest_b_cd8_plot_no_stats <- relabel_ggplot(forest_b_cd8_plot_no_stats, cell_type_label_map)

forest_mac_cd8_plot <- call_with_supported_args(
  panoramic::plot_forest,
  list(
    se_meta = se_meta,
    ct1 = "cd8_t_cells",
    ct2 = "cd68_cd163_macrophages",
    radius_um = 25,
    group_col = group_col_meta,
    show_est_se = TRUE,
    show_ci = TRUE
  )
)
forest_mac_cd8_plot <- relabel_ggplot(forest_mac_cd8_plot, cell_type_label_map)
forest_mac_cd8_plot_no_stats <- call_with_supported_args(
  panoramic::plot_forest,
  list(
    se_meta = se_meta,
    ct1 = "cd8_t_cells",
    ct2 = "cd68_cd163_macrophages",
    radius_um = 25,
    group_col = group_col_meta,
    show_est_se = FALSE,
    show_ci = FALSE
  )
)
forest_mac_cd8_plot_no_stats <- relabel_ggplot(forest_mac_cd8_plot_no_stats, cell_type_label_map)

network_plot <- call_with_supported_args(
  panoramic::plot_spatial_network,
  list(
    se_diff = se_meta,
    fdr_threshold = 0.05,
    z_sign = "negative",
    include_nonsig = TRUE,
    nonsig_max_fdr = 1.0,
    directed = FALSE,
    layout = "fr",
    sig_operator = "gt"
  )
)
network_plot <- relabel_ggplot(network_plot, cell_type_label_map)
network_plot <- network_plot + theme(
  legend.position = "none",
  plot.margin = margin(2, 2, 2, 2)
)

# ---- Enforce uniform title sizing across all exported plots ----
uniform_title_pt <- 8
uniform_subtitle_pt <- 7

volcano_plot <- apply_uniform_title_size(volcano_plot, uniform_title_pt, uniform_subtitle_pt)
volcano_plot_no_legend <- apply_uniform_title_size(volcano_plot_no_legend, uniform_title_pt, uniform_subtitle_pt)
volcano_plot_no_annotations <- apply_uniform_title_size(volcano_plot_no_annotations, uniform_title_pt, uniform_subtitle_pt)
volcano_plot_no_annotations_no_legend <- apply_uniform_title_size(volcano_plot_no_annotations_no_legend, uniform_title_pt, uniform_subtitle_pt)
rep_b_cd8_plot <- apply_uniform_title_size(rep_b_cd8_plot, uniform_title_pt, uniform_subtitle_pt)
rep_mac_cd8_plot <- apply_uniform_title_size(rep_mac_cd8_plot, uniform_title_pt, uniform_subtitle_pt)
forest_b_cd8_plot <- apply_uniform_title_size(forest_b_cd8_plot, uniform_title_pt, uniform_subtitle_pt)
forest_b_cd8_plot_no_stats <- apply_uniform_title_size(forest_b_cd8_plot_no_stats, uniform_title_pt, uniform_subtitle_pt)
forest_mac_cd8_plot <- apply_uniform_title_size(forest_mac_cd8_plot, uniform_title_pt, uniform_subtitle_pt)
forest_mac_cd8_plot_no_stats <- apply_uniform_title_size(forest_mac_cd8_plot_no_stats, uniform_title_pt, uniform_subtitle_pt)
network_plot <- apply_uniform_title_size(network_plot, uniform_title_pt, uniform_subtitle_pt)

# ---- Export individual panels (dimensions are height x width in spec) ----
save_plot_formats(
  volcano_plot,
  file.path(exports_dir, paste0(run_prefix, "_volcano_fdr_clean_labels_annotated_with_legend")),
  width_in = 3,
  height_in = 4
)
save_plot_formats(
  volcano_plot_no_legend,
  file.path(exports_dir, paste0(run_prefix, "_volcano_fdr_clean_labels_annotated_no_legend")),
  width_in = 3.25,
  height_in = 4
)
save_plot_formats(
  volcano_plot_no_annotations,
  file.path(exports_dir, paste0(run_prefix, "_volcano_fdr_clean_labels_no_annotations_with_legend")),
  width_in = 3.25,
  height_in = 4
)
save_plot_formats(
  volcano_plot_no_annotations_no_legend,
  file.path(exports_dir, paste0(run_prefix, "_volcano_fdr_clean_labels_no_annotations_no_legend")),
  width_in = 3.25,
  height_in = 4
)
save_plot_formats(
  rep_b_cd8_plot,
  file.path(exports_dir, paste0(run_prefix, "_representative_b_cell_to_cd8_t_cell_clean_labels")),
  width_in = 3.25,
  height_in = 2
)
save_plot_formats(
  rep_mac_cd8_plot,
  file.path(exports_dir, paste0(run_prefix, "_representative_cd68_cd163_macrophage_to_cd8_t_cell_clean_labels")),
  width_in = 3.25,
  height_in = 2
)
save_plot_formats(
  forest_b_cd8_plot,
  file.path(exports_dir, paste0(run_prefix, "_forest_b_cell_to_cd8_t_cell_clean_labels")),
  width_in = 3.25,
  height_in = 3
)
save_plot_formats(
  forest_b_cd8_plot_no_stats,
  file.path(exports_dir, paste0(run_prefix, "_forest_b_cell_to_cd8_t_cell_clean_labels_no_stats_text")),
  width_in = 3.25,
  height_in = 3
)
save_plot_formats(
  forest_mac_cd8_plot,
  file.path(exports_dir, paste0(run_prefix, "_forest_cd68_cd163_macrophage_to_cd8_t_cell_clean_labels")),
  width_in = 3.25,
  height_in = 3
)
save_plot_formats(
  forest_mac_cd8_plot_no_stats,
  file.path(exports_dir, paste0(run_prefix, "_forest_cd68_cd163_macrophage_to_cd8_t_cell_clean_labels_no_stats_text")),
  width_in = 3.25,
  height_in = 3
)
save_plot_formats(
  network_plot,
  file.path(exports_dir, paste0(run_prefix, "_network_clr_up_clean_labels")),
  width_in = 3,
  height_in = 1.5
)
save_plot_formats(
  network_plot,
  file.path(exports_dir, paste0(run_prefix, "_network_clr_up_clean_labels_larger")),
  width_in = 4,
  height_in = 2.25
)

message("Individual exports written to: ", exports_dir)
