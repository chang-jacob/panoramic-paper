#!/usr/bin/env Rscript

# HNSCC TMA Script: Exploratory Data Analysis
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Generate descriptive cohort summaries and exploratory HNSCC TMA visualizations.
# - Focus on LN_benign vs LN_met_adjacent comparisons for EDA outputs.
# - Produce manuscript/supplement-ready EDA tables and figures.

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

input_file <- file.path(paths$data, "processed", "hnscc_tma", "processed_hnscc_tma.csv")
fig_dir <- file.path(paths$output, "figures", "hnscc_tma", "eda", "ln_benign_vs_ln_met_adjacent")
tab_dir <- file.path(paths$output, "tables", "hnscc_tma", "eda", "ln_benign_vs_ln_met_adjacent")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Missing processed input file: ", input_file)
}

# ---- Parameters (edit in script) ----
analysis_params <- list(
  comparison_core_categories = c("LN_benign", "LN_met_adjacent"),
  max_cells_per_sample_map = 75000L,
  seed = 123L,
  map_dot_size = 0.5,
  map_dot_alpha = 1.0,
  map_legend_position = "bottom",
  pca_components = 5L,
  map_facet_ncol = 10L
)

# ---- Load data ----
df <- readr::read_csv(input_file, show_col_types = FALSE)

required_cols <- c("cell_id", "cell_type", "x", "y", "patient", "region", "core_category")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Input data missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- Palette definitions ----
celltype_palette <- c(
  "B" = "#6A4C93",
  "CAF" = "#D95F02",
  "CD4T" = "#1B9E77",
  "CD8T" = "#2CA25F",
  "DC" = "#4C78A8",
  "LEC" = "#B279A2",
  "NK" = "#5E3C99",
  "TFH" = "#998EC3",
  "Treg" = "#66A61E",
  "VSMC_myofibroblast" = "#E6AB02",
  "gdT" = "#1F78B4",
  "granulocyte" = "#33A02C",
  "macrophages" = "#A6CEE3",
  "mast" = "#FB9A99",
  "nerve" = "#B15928",
  "plasma" = "#CAB2D6",
  "stromal" = "#FDBF6F",
  "tumor" = "#E31A1C",
  "tumor_CA9" = "#FF7F00",
  "vasculature" = "#A6761D"
)

comparison_palette <- c(
  "LN_benign" = "#59A14F",
  "LN_met_adjacent" = "#E15759"
)

comparison_levels <- analysis_params$comparison_core_categories

# ---- Standardize and filter to comparison ----
df <- df |>
  mutate(
    cell_type = as.character(cell_type),
    core_category = as.character(core_category),
    sample_id = as.character(region),
    patient = as.character(patient)
  ) |>
  filter(core_category %in% comparison_levels)

if (nrow(df) == 0L) {
  stop(
    "No rows found for requested comparison categories: ",
    paste(comparison_levels, collapse = ", ")
  )
}

df <- df |>
  mutate(
    comparison_group = factor(core_category, levels = comparison_levels),
    cell_type = factor(cell_type, levels = names(eda_add_missing_palette_colors(celltype_palette, unique(cell_type))))
  )

celltype_palette <- eda_add_missing_palette_colors(
  base_palette = celltype_palette,
  observed_levels = sort(unique(as.character(df$cell_type)))
)
comparison_palette <- eda_add_missing_palette_colors(
  base_palette = comparison_palette,
  observed_levels = sort(unique(as.character(df$comparison_group)))
)

# ---- Core comparison tables ----
tbls <- eda_compute_composition_tables(
  df = df,
  group_col = "comparison_group",
  strata_col = "comparison_group",
  sample_id_col = "sample_id",
  patient_col = "patient",
  cell_type_col = "cell_type"
)

sample_balance <- tbls$sample_summary |>
  distinct(group_label, sample_id) |>
  count(group_label, name = "n_samples") |>
  arrange(group_label)

patient_balance <- tbls$sample_summary |>
  distinct(group_label, patient) |>
  count(group_label, name = "n_patients") |>
  arrange(group_label)

regions_per_patient <- tbls$sample_summary |>
  distinct(group_label, patient, sample_id) |>
  count(group_label, patient, name = "n_regions") |>
  arrange(group_label, desc(n_regions), patient)

readr::write_csv(tbls$celltype_overall, file.path(tab_dir, "eda_cmp_celltype_overall.csv"))
readr::write_csv(tbls$celltype_by_group, file.path(tab_dir, "eda_cmp_celltype_by_group.csv"))
readr::write_csv(tbls$celltype_by_sample, file.path(tab_dir, "eda_cmp_celltype_by_sample.csv"))
readr::write_csv(tbls$sample_summary, file.path(tab_dir, "eda_cmp_sample_summary.csv"))
readr::write_csv(tbls$celltype_prevalence_by_group, file.path(tab_dir, "eda_cmp_celltype_prevalence_by_group.csv"))
readr::write_csv(sample_balance, file.path(tab_dir, "eda_cmp_sample_balance.csv"))
readr::write_csv(patient_balance, file.path(tab_dir, "eda_cmp_patient_balance.csv"))
readr::write_csv(regions_per_patient, file.path(tab_dir, "eda_cmp_regions_per_patient.csv"))

# ---- Composition plots ----
comparison_title <- "HNSCC LN_benign vs LN_met_adjacent"

p_overall <- ggplot(
  tbls$celltype_overall,
  aes(x = fct_reorder(as.character(cell_type), n_cells), y = prop, fill = cell_type)
) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = celltype_palette, drop = FALSE) +
  labs(
    title = paste0(comparison_title, ": cell type composition (overall)"),
    x = "Cell type",
    y = "Proportion of all cells"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

eda_save_plot(p_overall, file.path(fig_dir, "eda_cmp_celltype_overall"), width = 8.8, height = 6.0)

p_by_group <- ggplot(
  tbls$celltype_by_group,
  aes(x = fct_reorder(as.character(cell_type), prop, .fun = max), y = prop, fill = group_label)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = comparison_palette, drop = FALSE) +
  labs(
    title = paste0(comparison_title, ": cell type composition by comparison group"),
    x = "Cell type",
    y = "Proportion within group",
    fill = "Comparison group"
  ) +
  theme_classic(base_size = 11)

eda_save_plot(p_by_group, file.path(fig_dir, "eda_cmp_celltype_by_group"), width = 9.4, height = 6.5)

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
  scale_fill_manual(values = comparison_palette, drop = FALSE) +
  scale_y_continuous(labels = function(x) scales::comma(x, accuracy = 1)) +
  labs(
    title = paste0(comparison_title, ": sample-level metrics"),
    x = "Comparison group",
    y = "Value",
    fill = "Comparison group"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

eda_save_plot(p_metrics, file.path(fig_dir, "eda_cmp_sample_metric_distributions"), width = 10.5, height = 4.5)

p_celltype_violin <- ggplot(
  tbls$celltype_by_sample,
  aes(x = group_label, y = prop, fill = group_label)
) +
  geom_violin(trim = FALSE, alpha = 0.25, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.08, height = 0, size = 0.5, alpha = 0.45) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = comparison_palette, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = paste0(comparison_title, ": sample-level cell type proportions"),
    x = "Comparison group",
    y = "Proportion within sample",
    fill = "Comparison group"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(size = 8, face = "bold"))

eda_save_plot(p_celltype_violin, file.path(fig_dir, "eda_cmp_sample_celltype_proportion_violin"), width = 12.5, height = 11.5)

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
  scale_color_manual(values = comparison_palette, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = paste0(comparison_title, ": cell type prevalence by group"),
    x = "Fraction of samples containing cell type",
    y = "Cell type",
    color = "Comparison group"
  ) +
  theme_classic(base_size = 11)

eda_save_plot(p_prevalence, file.path(fig_dir, "eda_cmp_celltype_prevalence_by_group"), width = 9.8, height = 6.8)

p_sample_balance <- ggplot(sample_balance, aes(x = group_label, y = n_samples, fill = group_label)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = comparison_palette, drop = FALSE) +
  labs(
    title = paste0(comparison_title, ": sample balance"),
    x = "Comparison group",
    y = "Number of samples"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

eda_save_plot(p_sample_balance, file.path(fig_dir, "eda_cmp_sample_balance"), width = 7.0, height = 4.8)

p_patient_balance <- ggplot(patient_balance, aes(x = group_label, y = n_patients, fill = group_label)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = comparison_palette, drop = FALSE) +
  labs(
    title = paste0(comparison_title, ": patient balance"),
    x = "Comparison group",
    y = "Number of patients"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

eda_save_plot(p_patient_balance, file.path(fig_dir, "eda_cmp_patient_balance"), width = 7.0, height = 4.8)

p_regions_per_patient <- ggplot(regions_per_patient, aes(x = group_label, y = n_regions, fill = group_label)) +
  geom_violin(trim = FALSE, alpha = 0.25, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, alpha = 0.7) +
  scale_fill_manual(values = comparison_palette, drop = FALSE) +
  labs(
    title = paste0(comparison_title, ": regions per patient"),
    x = "Comparison group",
    y = "Regions per patient",
    fill = "Comparison group"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

eda_save_plot(p_regions_per_patient, file.path(fig_dir, "eda_cmp_regions_per_patient"), width = 7.0, height = 4.8)

pca_out <- eda_run_composition_pca(
  celltype_by_sample = tbls$celltype_by_sample,
  max_components = analysis_params$pca_components
)
if (!is.null(pca_out)) {
  readr::write_csv(pca_out$scores, file.path(tab_dir, "eda_cmp_sample_composition_pca_scores.csv"))
  readr::write_csv(pca_out$loadings, file.path(tab_dir, "eda_cmp_sample_composition_pca_loadings.csv"))
  readr::write_csv(pca_out$variance, file.path(tab_dir, "eda_cmp_sample_composition_pca_variance.csv"))

  p_pca <- eda_plot_pca_scores(
    scores = pca_out$scores,
    variance_tbl = pca_out$variance,
    title = paste0(comparison_title, ": sample composition PCA"),
    group_palette = comparison_palette
  )
  eda_save_plot(p_pca, file.path(fig_dir, "eda_cmp_sample_composition_pca"), width = 8.3, height = 6.0)
}

# ---- Spatial overview panels for comparison ----
eda_plot_spatial_overview(
  df_sub = df,
  sample_col = "sample_id",
  x_col = "x",
  y_col = "y",
  cell_type_col = "cell_type",
  palette = celltype_palette,
  title = "HNSCC LN_benign vs LN_met_adjacent: all samples",
  file_stub = file.path(fig_dir, "eda_cmp_spatial_samples_all"),
  max_cells_per_sample = analysis_params$max_cells_per_sample_map,
  seed = analysis_params$seed,
  map_dot_size = analysis_params$map_dot_size,
  map_dot_alpha = analysis_params$map_dot_alpha,
  legend_position = analysis_params$map_legend_position,
  facet_ncol = analysis_params$map_facet_ncol
)

eda_plot_spatial_overview(
  df_sub = df |> filter(comparison_group == "LN_benign"),
  sample_col = "sample_id",
  x_col = "x",
  y_col = "y",
  cell_type_col = "cell_type",
  palette = celltype_palette,
  title = "HNSCC LN_benign samples",
  file_stub = file.path(fig_dir, "eda_cmp_spatial_samples_ln_benign"),
  max_cells_per_sample = analysis_params$max_cells_per_sample_map,
  seed = analysis_params$seed,
  map_dot_size = analysis_params$map_dot_size,
  map_dot_alpha = analysis_params$map_dot_alpha,
  legend_position = analysis_params$map_legend_position,
  facet_ncol = analysis_params$map_facet_ncol
)

eda_plot_spatial_overview(
  df_sub = df |> filter(comparison_group == "LN_met_adjacent"),
  sample_col = "sample_id",
  x_col = "x",
  y_col = "y",
  cell_type_col = "cell_type",
  palette = celltype_palette,
  title = "HNSCC LN_met_adjacent samples",
  file_stub = file.path(fig_dir, "eda_cmp_spatial_samples_ln_met_adjacent"),
  max_cells_per_sample = analysis_params$max_cells_per_sample_map,
  seed = analysis_params$seed,
  map_dot_size = analysis_params$map_dot_size,
  map_dot_alpha = analysis_params$map_dot_alpha,
  legend_position = analysis_params$map_legend_position,
  facet_ncol = analysis_params$map_facet_ncol
)

message("HNSCC TMA comparison EDA complete.")
message("Figures written to: ", fig_dir)
message("Tables written to: ", tab_dir)
