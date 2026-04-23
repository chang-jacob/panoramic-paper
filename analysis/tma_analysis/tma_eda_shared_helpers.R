# TMA Analysis Helper Script: Shared EDA Plot/Table Utilities
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Shared helper functions for TMA exploratory analyses.
# - Provides reusable plotting, palette, composition, and PCA utility routines.

eda_save_plot <- function(p, file_stub, width, height, dpi = 300) {
  png_path <- paste0(file_stub, ".png")
  pdf_path <- paste0(file_stub, ".pdf")
  ggplot2::ggsave(filename = png_path, plot = p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(filename = pdf_path, plot = p, width = width, height = height)
  invisible(list(png = png_path, pdf = pdf_path))
}

eda_add_missing_palette_colors <- function(base_palette, observed_levels) {
  observed_levels <- unique(as.character(observed_levels))
  missing_levels <- setdiff(observed_levels, names(base_palette))
  if (length(missing_levels) == 0L) {
    return(base_palette)
  }

  extra_cols <- scales::hue_pal()(length(missing_levels))
  names(extra_cols) <- missing_levels
  c(base_palette, extra_cols)
}

eda_downsample_for_maps <- function(x, sample_col, max_per_sample, seed = 123L) {
  set.seed(seed)
  x |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sample_col))) |>
    dplyr::group_modify(~ {
      if (nrow(.x) > max_per_sample) {
        dplyr::slice_sample(.x, n = max_per_sample)
      } else {
        .x
      }
    }) |>
    dplyr::ungroup()
}

eda_plot_spatial_overview <- function(
    df_sub,
    sample_col,
    x_col,
    y_col,
    cell_type_col,
    palette,
    title,
    file_stub,
    max_cells_per_sample = 75000L,
    seed = 123L,
    map_dot_size = 0.5,
    map_dot_alpha = 1,
    legend_position = "right",
    facet_ncol = NULL
) {
  if (nrow(df_sub) == 0L) {
    message("Skipping spatial plot (no rows): ", title)
    return(invisible(NULL))
  }

  plot_df <- eda_downsample_for_maps(
    x = df_sub,
    sample_col = sample_col,
    max_per_sample = max_cells_per_sample,
    seed = seed
  )
  n_samples <- dplyr::n_distinct(plot_df[[sample_col]])
  if (n_samples == 0L) {
    message("Skipping spatial plot (no samples): ", title)
    return(invisible(NULL))
  }

  if (is.null(facet_ncol)) {
    ncol_facets <- if (n_samples <= 24L) 6L else if (n_samples <= 48L) 8L else 10L
  } else {
    ncol_facets <- suppressWarnings(as.integer(facet_ncol))
    if (!is.finite(ncol_facets) || ncol_facets < 1L) {
      ncol_facets <- if (n_samples <= 24L) 6L else if (n_samples <= 48L) 8L else 10L
    }
  }
  nrow_facets <- ceiling(n_samples / ncol_facets)
  plot_width <- min(24, max(10, ncol_facets * 1.9))
  plot_height <- min(30, max(8, nrow_facets * 1.6))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[cell_type_col]]
    )
  ) +
    ggplot2::geom_point(size = map_dot_size, alpha = map_dot_alpha, shape = 16, stroke = 0) +
    ggplot2::coord_equal() +
    ggplot2::facet_wrap(stats::as.formula(paste0("~", sample_col)), ncol = ncol_facets) +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::labs(title = title, color = "Cell type") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.text = ggplot2::element_text(size = 7, face = "bold"),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.key.width = grid::unit(0.35, "cm")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2.8, alpha = 1)))

  eda_save_plot(p, file_stub = file_stub, width = plot_width, height = plot_height)
  invisible(p)
}

eda_compute_composition_tables <- function(
    df,
    group_col,
    strata_col,
    sample_id_col,
    patient_col = "patient",
    cell_type_col = "cell_type"
) {
  df_std <- data.frame(
    group_label = as.character(df[[group_col]]),
    strata = as.character(df[[strata_col]]),
    patient = as.character(df[[patient_col]]),
    sample_id = as.character(df[[sample_id_col]]),
    cell_type = as.character(df[[cell_type_col]]),
    stringsAsFactors = FALSE
  )

  celltype_overall <- df_std |>
    dplyr::count(cell_type, name = "n_cells") |>
    dplyr::mutate(prop = n_cells / sum(n_cells)) |>
    dplyr::arrange(dplyr::desc(n_cells))

  celltype_by_group <- df_std |>
    dplyr::count(group_label, cell_type, name = "n_cells") |>
    dplyr::group_by(group_label) |>
    dplyr::mutate(prop = n_cells / sum(n_cells)) |>
    dplyr::ungroup() |>
    dplyr::arrange(group_label, dplyr::desc(n_cells))

  celltype_by_sample <- df_std |>
    dplyr::count(group_label, strata, patient, sample_id, cell_type, name = "n_cells") |>
    dplyr::group_by(group_label, strata, patient, sample_id) |>
    dplyr::mutate(sample_total_cells = sum(n_cells), prop = n_cells / sample_total_cells) |>
    dplyr::ungroup()

  sample_summary <- celltype_by_sample |>
    dplyr::group_by(group_label, strata, patient, sample_id) |>
    dplyr::summarise(
      n_cells = dplyr::first(sample_total_cells),
      n_cell_types = sum(n_cells > 0),
      shannon = {
        p <- prop[prop > 0]
        -sum(p * log(p))
      },
      .groups = "drop"
    ) |>
    dplyr::arrange(group_label, strata, dplyr::desc(n_cells))

  celltype_prevalence_by_group <- celltype_by_sample |>
    dplyr::mutate(present = n_cells > 0) |>
    dplyr::group_by(group_label, cell_type) |>
    dplyr::summarise(
      n_samples = dplyr::n_distinct(sample_id),
      n_samples_present = sum(present),
      prevalence = n_samples_present / n_samples,
      .groups = "drop"
    ) |>
    dplyr::arrange(group_label, dplyr::desc(prevalence), cell_type)

  list(
    df_std = df_std,
    celltype_overall = celltype_overall,
    celltype_by_group = celltype_by_group,
    celltype_by_sample = celltype_by_sample,
    sample_summary = sample_summary,
    celltype_prevalence_by_group = celltype_prevalence_by_group
  )
}

eda_run_composition_pca <- function(celltype_by_sample, max_components = 5L) {
  mat_df <- celltype_by_sample |>
    dplyr::select(sample_id, group_label, strata, patient, cell_type, prop) |>
    dplyr::group_by(sample_id, group_label, strata, patient, cell_type) |>
    dplyr::summarise(prop = sum(prop), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = cell_type,
      values_from = prop,
      values_fill = 0
    )

  feature_cols <- setdiff(colnames(mat_df), c("sample_id", "group_label", "strata", "patient"))
  if (length(feature_cols) < 2L) return(NULL)

  x <- as.matrix(mat_df[, feature_cols, drop = FALSE])
  rownames(x) <- mat_df$sample_id
  keep_features <- apply(x, 2, stats::sd) > 0
  if (sum(keep_features) < 2L) return(NULL)

  x <- x[, keep_features, drop = FALSE]
  pca <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  n_pc <- min(max_components, ncol(pca$x))

  scores <- as.data.frame(pca$x[, seq_len(n_pc), drop = FALSE], stringsAsFactors = FALSE)
  scores$sample_id <- rownames(scores)
  scores <- dplyr::left_join(scores, mat_df[, c("sample_id", "group_label", "strata", "patient")], by = "sample_id")

  loadings <- as.data.frame(pca$rotation[, seq_len(n_pc), drop = FALSE], stringsAsFactors = FALSE)
  loadings$cell_type <- rownames(loadings)
  loadings <- dplyr::relocate(loadings, cell_type)

  var_explained <- data.frame(
    pc = paste0("PC", seq_along(pca$sdev)),
    variance_explained = (pca$sdev^2) / sum(pca$sdev^2),
    stringsAsFactors = FALSE
  )

  list(scores = scores, loadings = loadings, variance = var_explained)
}

eda_plot_pca_scores <- function(scores, variance_tbl, title, group_palette) {
  vx <- variance_tbl$variance_explained[variance_tbl$pc == "PC1"][1]
  vy <- variance_tbl$variance_explained[variance_tbl$pc == "PC2"][1]
  x_lab <- paste0("PC1 (", scales::percent(vx, accuracy = 0.1), ")")
  y_lab <- paste0("PC2 (", scales::percent(vy, accuracy = 0.1), ")")

  ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = group_label, shape = strata)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::scale_color_manual(values = group_palette, drop = FALSE) +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = "Group",
      shape = "Strata"
    ) +
    ggplot2::theme_classic(base_size = 11)
}
