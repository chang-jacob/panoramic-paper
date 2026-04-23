# CRC TMA Script: Preprocess
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Load raw CRC TMA cell-level and patient annotation tables.
# - Harmonize labels/columns and convert spatial coordinates to microns.
# - Filter excluded cell types and write processed CRC TMA CSV for downstream scripts.

suppressPackageStartupMessages({
  library(config)
  library(dplyr)
  library(janitor)
  library(readr)
  library(stringr)
})

# ---- Load paths ----
config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)

raw_dir <- file.path(paths$data, "crc_tma")
processed_dir <- file.path(paths$data, "processed", "crc_tma")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Input files ----
raw_cells_file <- file.path(raw_dir, "CRC_clusters_neighborhoods_markers.csv")
raw_anno_file  <- file.path(raw_dir, "Patient_data_TMA_annotations.csv")

if (!file.exists(raw_cells_file)) {
  stop("Missing raw cells file: ", raw_cells_file)
}
if (!file.exists(raw_anno_file)) {
  stop("Missing raw annotations file: ", raw_anno_file)
}

# ---- Load ----
# Note: spicyR paper drops the first column (empty) and uses clean_names()
raw_cells <- read_csv(raw_cells_file, show_col_types = FALSE) |>
  select(-1) |>
  clean_names()

raw_anno <- read_csv(raw_anno_file, show_col_types = FALSE, locale = locale(encoding = "UTF-8")) |>
  clean_names()

# ---- Standardize ----
# Match spicyR naming: cluster_name -> cell_type, x_x/y_y -> x/y
cells <- raw_cells |>
  rename(
    cell_id = cell_id,
    cell_type = cluster_name,
    file_name = file_name,
    patient = patients,
    group = groups,
    spot = spots,
    x = x_x,
    y = y_y
  ) |>
  mutate(
    cell_type = as.character(cell_type)
  )

anno <- raw_anno |>
  rename(
    patient = patient,
    group = group,
    tma_spot = tma_spot,
    full_histology = full_histology,
    simple_tumor_location = simple_tumor_location
  )

# ---- Merge annotations (by patient) ----
# The spicyR paper maps image -> patient via codex data, and then joins to clinical data.
# Here we join by patient directly to attach annotations.
if (!"patient" %in% colnames(cells)) {
  stop("cells is missing required 'patient' column.")
}
if (!"patient" %in% colnames(anno)) {
  stop("anno is missing required 'patient' column.")
}
df <- cells |>
  left_join(anno, by = "patient", suffix = c("", "_anno")) |>
  mutate(group = dplyr::coalesce(group_anno, group)) |>
  select(-group_anno)

# ---- Match spicyR naming + filtering ----
# spicyR: cellType <- factor(cluster_name, levels = lev, labels = make_clean_names(lev))
lev <- unique(df$cell_type)
df <- df |>
  mutate(
    cell_type = factor(cell_type, levels = lev, labels = make_clean_names(lev)),
    cell_type = as.character(cell_type)
  )

# ---- Rescale coordinates to microns ----
# Source coords are in pixels; convert using dataset-specific micron-per-pixel.
MICRONS_PER_PIXEL <- 0.37742 # This conversion is documented in Schürch et al. manuscript
df <- df |>
  mutate(
    x = x * MICRONS_PER_PIXEL,
    y = y * MICRONS_PER_PIXEL
  )

# spicyR: drop cell types with >21% zero images; remove "dirt"
cell_prop <- table(df$file_name, df$cell_type)
drop_cell_types <- names(which(colMeans(cell_prop == 0) > 0.21))
drop_cell_types <- unique(c(drop_cell_types, "dirt"))

df <- df |>
  filter(!cell_type %in% drop_cell_types)

# Rename groups
df <- df |>
  mutate(
    group = as.character(group),
    group_label = case_when(
      group %in% c("1", "CLR") ~ "CLR",
      group %in% c("2", "DII") ~ "DII",
      TRUE ~ group
    )
  )

df <- df %>% 
  dplyr::filter(cell_type != "undefined")

# ---- Basic validation ----
required_cols <- c("cell_id", "cell_type", "x", "y", "patient", "group", "file_name")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Processed data missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- Save ----
processed_file <- file.path(processed_dir, "processed_crc_tma.csv")
write_csv(df, processed_file)

message("Wrote processed data to: ", processed_file)
