#!/usr/bin/env Rscript

# CRC TMA Script: Exploratory Data Analysis
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Generate descriptive cohort summaries and exploratory CRC TMA visualizations.
# - Produce manuscript/supplement-ready EDA tables and figures.
# - Save outputs under the CRC TMA EDA figure/table directories.

suppressPackageStartupMessages({
  library(config)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(forcats)
})

# ---- Load paths ----
config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)
source(file.path(getwd(), "analysis", "tma_analysis", "tma_eda_shared_helpers.R"))

input_file <- file.path(paths$data, "processed", "crc_tma", "processed_crc_tma.csv")
fig_dir <- file.path(paths$output, "figures", "crc_tma", "eda")
tab_dir <- file.path(paths$output, "tables", "crc_tma", "eda")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Missing processed input file: ", input_file)
}

# ---- Parameters (edit in script) ----
analysis_params <- list(
  max_cells_per_sample_map = 75000L,
  seed = 123L,
  map_dot_size = 0.5,
  map_dot_alpha = 1.0,
  map_legend_position = "bottom",
  pca_components = 5L
)

# ---- Load data ----
df <- readr::read_csv(input_file, show_col_types = FALSE)

# ---- Validation ----
required_cols <- c(
  "cell_id", "cell_type", "x", "y",
  "patient", "group", "group_label",
  "file_name", "spot", "tma_ab"
)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Input data missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- Cell type palette ----
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

group_palette <- c(
  "CLR" = "#2F5597",
  "DII" = "#C00000"
)

# ---- Standardize fields used by helper workflows ----
df <- df |>
  mutate(
    cell_type = as.character(cell_type),
    group_label = as.character(group_label),
    tma_ab = toupper(as.character(tma_ab)),
    sample_id = as.character(spot),
    strata = as.character(tma_ab)
  )

celltype_palette <- eda_add_missing_palette_colors(
  base_palette = celltype_palette,
  observed_levels = sort(unique(df$cell_type))
)
group_palette <- eda_add_missing_palette_colors(
  base_palette = group_palette,
  observed_levels = sort(unique(df$group_label))
)

group_levels <- c("CLR", "DII", setdiff(sort(unique(df$group_label)), c("CLR", "DII")))
df <- df |>
  mutate(
    group_label = factor(group_label, levels = unique(group_levels)),
    strata = factor(strata, levels = c("A", "B", setdiff(sort(unique(strata)), c("A", "B")))),
    cell_type = factor(cell_type, levels = names(celltype_palette))
  )

# ---- Core summary tables ----
tbls <- eda_compute_composition_tables(
  df = df,
  group_col = "group_label",
  strata_col = "strata",
  sample_id_col = "sample_id",
  patient_col = "patient",
  cell_type_col = "cell_type"
)

sample_balance <- tbls$sample_summary |>
  distinct(group_label, strata, sample_id) |>
  count(group_label, strata, name = "n_samples") |>
  arrange(group_label, strata)

patient_balance <- tbls$sample_summary |>
  distinct(group_label, patient) |>
  count(group_label, name = "n_patients") |>
  arrange(group_label)

regions_per_patient <- tbls$sample_summary |>
  distinct(group_label, patient, sample_id) |>
  count(group_label, patient, name = "n_regions") |>
  arrange(group_label, desc(n_regions), patient)

readr::write_csv(tbls$celltype_overall, file.path(tab_dir, "eda_celltype_overall.csv"))
readr::write_csv(tbls$celltype_by_group, file.path(tab_dir, "eda_celltype_by_group.csv"))
readr::write_csv(tbls$celltype_by_sample, file.path(tab_dir, "eda_celltype_by_sample.csv"))
readr::write_csv(tbls$sample_summary, file.path(tab_dir, "eda_sample_summary.csv"))
readr::write_csv(tbls$celltype_prevalence_by_group, file.path(tab_dir, "eda_celltype_prevalence_by_group.csv"))
readr::write_csv(sample_balance, file.path(tab_dir, "eda_sample_balance_by_group_tma.csv"))
readr::write_csv(patient_balance, file.path(tab_dir, "eda_patient_balance_by_group.csv"))
readr::write_csv(regions_per_patient, file.path(tab_dir, "eda_regions_per_patient.csv"))

# ---- Composition plots ----
p_overall <- ggplot(
  tbls$celltype_overall,
  aes(x = fct_reorder(as.character(cell_type), n_cells), y = prop, fill = cell_type)
) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = celltype_palette, drop = FALSE) +
  labs(
    title = "CRC TMA cell type composition (overall)",
    x = "Cell type",
    y = "Proportion of all cells"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

eda_save_plot(p_overall, file.path(fig_dir, "eda_celltype_overall"), width = 8.5, height = 5.5)

p_by_group <- ggplot(
  tbls$celltype_by_group,
  aes(x = fct_reorder(as.character(cell_type), prop, .fun = max), y = prop, fill = group_label)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  labs(
    title = "Cell type composition by group",
    x = "Cell type",
    y = "Proportion within group",
    fill = "Group"
  ) +
  theme_classic(base_size = 11)

eda_save_plot(p_by_group, file.path(fig_dir, "eda_celltype_by_group"), width = 9, height = 6.2)

# ---- Sample-level distribution plots ----
metrics_long <- tbls$sample_summary |>
  pivot_longer(
    cols = c(n_cells, n_cell_types, shannon),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(metric = factor(metric, levels = c("n_cells", "n_cell_types", "shannon")))

p_metrics <- ggplot(metrics_long, aes(x = group_label, y = value, fill = group_label)) +
  geom_violin(trim = FALSE, alpha = 0.25, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, height = 0, size = 1, alpha = 0.6) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  scale_y_continuous(labels = function(x) scales::comma(x, accuracy = 1)) +
  labs(
    title = "Sample-level metric distributions by group",
    x = "Group",
    y = "Value",
    fill = "Group"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

eda_save_plot(p_metrics, file.path(fig_dir, "eda_sample_metric_distributions_by_group"), width = 10.5, height = 4.5)

p_celltype_violin <- ggplot(
  tbls$celltype_by_sample,
  aes(x = group_label, y = prop, fill = group_label)
) +
  geom_violin(trim = FALSE, alpha = 0.25, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.08, height = 0, size = 0.5, alpha = 0.45) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Sample-level cell type proportions by group",
    x = "Group",
    y = "Proportion within sample",
    fill = "Group"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(size = 8, face = "bold"))

eda_save_plot(p_celltype_violin, file.path(fig_dir, "eda_sample_celltype_proportion_violin_by_group"), width = 12.5, height = 11.5)

p_sample_counts <- ggplot(tbls$sample_summary, aes(x = strata, y = n_cells, fill = group_label)) +
  geom_violin(trim = FALSE, alpha = 0.25, color = NA, position = position_dodge(width = 0.85)) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.9, position = position_dodge(width = 0.85)) +
  geom_jitter(aes(color = group_label), width = 0.1, height = 0, size = 1, alpha = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  scale_color_manual(values = group_palette, drop = FALSE) +
  scale_y_continuous(labels = function(x) scales::comma(x, accuracy = 1)) +
  labs(
    title = "Number of cells per sample by TMA and group",
    x = "TMA",
    y = "Cells per sample",
    fill = "Group"
  ) +
  theme_classic(base_size = 11)

eda_save_plot(p_sample_counts, file.path(fig_dir, "eda_sample_cell_counts_by_tma_group"), width = 8.2, height = 5.2)

# ---- Additional supplementary EDA ----
p_prevalence <- ggplot(
  tbls$celltype_prevalence_by_group,
  aes(
    x = prevalence,
    y = fct_reorder(as.character(cell_type), prevalence, .fun = max),
    color = group_label,
    group = cell_type
  )
) +
  geom_line(color = "grey75", linewidth = 0.35) +
  geom_point(size = 2.2) +
  scale_color_manual(values = group_palette, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Cell type prevalence across samples by group",
    x = "Fraction of samples containing cell type",
    y = "Cell type",
    color = "Group"
  ) +
  theme_classic(base_size = 11)

eda_save_plot(p_prevalence, file.path(fig_dir, "eda_celltype_prevalence_by_group"), width = 9.8, height = 6.8)

p_balance <- ggplot(sample_balance, aes(x = strata, y = n_samples, fill = group_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  labs(
    title = "Sample balance by TMA and group",
    x = "TMA",
    y = "Number of samples",
    fill = "Group"
  ) +
  theme_classic(base_size = 11)

eda_save_plot(p_balance, file.path(fig_dir, "eda_sample_balance_by_group_tma"), width = 8.2, height = 5.0)

p_regions_per_patient <- ggplot(regions_per_patient, aes(x = group_label, y = n_regions, fill = group_label)) +
  geom_violin(trim = FALSE, alpha = 0.25, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, alpha = 0.7) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  labs(
    title = "Number of regions per patient by group",
    x = "Group",
    y = "Regions per patient",
    fill = "Group"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

eda_save_plot(p_regions_per_patient, file.path(fig_dir, "eda_regions_per_patient_by_group"), width = 7.0, height = 4.8)

pca_out <- eda_run_composition_pca(
  celltype_by_sample = tbls$celltype_by_sample,
  max_components = analysis_params$pca_components
)
if (!is.null(pca_out)) {
  readr::write_csv(pca_out$scores, file.path(tab_dir, "eda_sample_composition_pca_scores.csv"))
  readr::write_csv(pca_out$loadings, file.path(tab_dir, "eda_sample_composition_pca_loadings.csv"))
  readr::write_csv(pca_out$variance, file.path(tab_dir, "eda_sample_composition_pca_variance.csv"))

  p_pca <- eda_plot_pca_scores(
    scores = pca_out$scores,
    variance_tbl = pca_out$variance,
    title = "CRC sample composition PCA",
    group_palette = group_palette
  )
  eda_save_plot(p_pca, file.path(fig_dir, "eda_sample_composition_pca"), width = 8.3, height = 6.0)
}

# ---- Spatial overview panels (all + by TMA + by group) ----
eda_plot_spatial_overview(
  df_sub = df,
  sample_col = "sample_id",
  x_col = "x",
  y_col = "y",
  cell_type_col = "cell_type",
  palette = celltype_palette,
  title = "CRC TMA: all sample spatial architectures",
  file_stub = file.path(fig_dir, "eda_spatial_samples_all"),
  max_cells_per_sample = analysis_params$max_cells_per_sample_map,
  seed = analysis_params$seed,
  map_dot_size = analysis_params$map_dot_size,
  map_dot_alpha = analysis_params$map_dot_alpha,
  legend_position = analysis_params$map_legend_position
)

eda_plot_spatial_overview(
  df_sub = df |> filter(strata == "A"),
  sample_col = "sample_id",
  x_col = "x",
  y_col = "y",
  cell_type_col = "cell_type",
  palette = celltype_palette,
  title = "CRC TMA A: sample spatial architectures",
  file_stub = file.path(fig_dir, "eda_spatial_samples_tma_A"),
  max_cells_per_sample = analysis_params$max_cells_per_sample_map,
  seed = analysis_params$seed,
  map_dot_size = analysis_params$map_dot_size,
  map_dot_alpha = analysis_params$map_dot_alpha,
  legend_position = analysis_params$map_legend_position
)

eda_plot_spatial_overview(
  df_sub = df |> filter(strata == "B"),
  sample_col = "sample_id",
  x_col = "x",
  y_col = "y",
  cell_type_col = "cell_type",
  palette = celltype_palette,
  title = "CRC TMA B: sample spatial architectures",
  file_stub = file.path(fig_dir, "eda_spatial_samples_tma_B"),
  max_cells_per_sample = analysis_params$max_cells_per_sample_map,
  seed = analysis_params$seed,
  map_dot_size = analysis_params$map_dot_size,
  map_dot_alpha = analysis_params$map_dot_alpha,
  legend_position = analysis_params$map_legend_position
)

eda_plot_spatial_overview(
  df_sub = df |> filter(group_label == "CLR"),
  sample_col = "sample_id",
  x_col = "x",
  y_col = "y",
  cell_type_col = "cell_type",
  palette = celltype_palette,
  title = "CLR samples: spatial architectures",
  file_stub = file.path(fig_dir, "eda_spatial_samples_group_CLR"),
  max_cells_per_sample = analysis_params$max_cells_per_sample_map,
  seed = analysis_params$seed,
  map_dot_size = analysis_params$map_dot_size,
  map_dot_alpha = analysis_params$map_dot_alpha,
  legend_position = analysis_params$map_legend_position
)

eda_plot_spatial_overview(
  df_sub = df |> filter(group_label == "DII"),
  sample_col = "sample_id",
  x_col = "x",
  y_col = "y",
  cell_type_col = "cell_type",
  palette = celltype_palette,
  title = "DII samples: spatial architectures",
  file_stub = file.path(fig_dir, "eda_spatial_samples_group_DII"),
  max_cells_per_sample = analysis_params$max_cells_per_sample_map,
  seed = analysis_params$seed,
  map_dot_size = analysis_params$map_dot_size,
  map_dot_alpha = analysis_params$map_dot_alpha,
  legend_position = analysis_params$map_legend_position
)

message("CRC TMA EDA complete.")
message("Figures written to: ", fig_dir)
message("Tables written to: ", tab_dir)
