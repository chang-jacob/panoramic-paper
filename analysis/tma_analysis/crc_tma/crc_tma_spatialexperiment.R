# CRC TMA Script: Build SpatialExperiment Lists
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Read processed CRC TMA table and core-level annotations.
# - Build per-spot SpatialExperiment objects for each analysis subset.
# - Save all/cohort-filtered/one-per-patient SPE list RDS files.

suppressPackageStartupMessages({
  library(config)
  library(dplyr)
  library(janitor)
  library(readr)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
})

# ---- Load paths ----
config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)

processed_dir <- file.path(paths$data, "processed", "crc_tma")
raw_dir <- file.path(paths$data, "crc_tma")
input_file <- file.path(processed_dir, "processed_crc_tma.csv")
output_file <- file.path(processed_dir, "crc_tma_spe_list.rds")
output_file_one_per_patient <- file.path(processed_dir, "crc_tma_spe_list_one_per_patient.rds")
output_file_core_filtered <- file.path(processed_dir, "crc_tma_spe_list_core_filtered.rds")
output_file_diffuse_only <- file.path(processed_dir, "crc_tma_spe_list_diffuse_only.rds")

if (!file.exists(input_file)) {
  stop("Missing processed input file: ", input_file)
}

# ---- Load processed data ----
df <- read_csv(input_file, show_col_types = FALSE)

# ---- Basic validation ----
required_cols <- c("cell_id", "cell_type", "x", "y", "patient", "spot", "file_name", "group_label")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Input data missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- Load core annotations (file_name -> core_annotation) ----
core_anno_file <- file.path(raw_dir, "crc_tma_core_annotation.csv")
if (!file.exists(core_anno_file)) {
  stop("Missing core annotations file: ", core_anno_file)
}

core_anno <- read_csv(core_anno_file, show_col_types = FALSE) |>
  clean_names() %>% 
  dplyr::rename(core_annotation = core_type) %>% 
  dplyr::mutate(file_name = paste0("reg", file_name))

required_core_cols <- c("file_name", "core_annotation")
missing_core_cols <- setdiff(required_core_cols, colnames(core_anno))
if (length(missing_core_cols) > 0) {
  stop(
    "Core annotations file missing required columns: ",
    paste(missing_core_cols, collapse = ", ")
  )
}

core_anno <- core_anno |>
  select(file_name, core_annotation) |>
  distinct()

conflicting_core <- core_anno |>
  dplyr::count(file_name) |>
  dplyr::filter(n > 1)
if (nrow(conflicting_core) > 0) {
  stop("Core annotations have multiple labels for some file_name values.")
}

df <- df |>
  left_join(core_anno, by = "file_name")

if (any(is.na(df$core_annotation))) {
  stop("Some rows are missing core_annotation after joining by file_name.")
}

# ---- Build SpatialExperiment list ----
make_spe_list <- function(input_df) {
  lapply(split(input_df, input_df$spot), function(subdf) {
    coords <- as.matrix(subdf[, c("x", "y")])
    rownames(coords) <- subdf$cell_id

    col_data <- S4Vectors::DataFrame(subdf)
    rownames(col_data) <- subdf$cell_id

    # Placeholder assay (PANORAMIC doesn't need expression values here)
    counts <- matrix(0, nrow = 1, ncol = nrow(subdf),
                     dimnames = list("placeholder", subdf$cell_id))

    SpatialExperiment(
      assays = list(counts = counts),
      colData = col_data,
      spatialCoords = coords
    )
  })
}

spe_list_all <- make_spe_list(df)

core_filtered_df <- df |>
  filter(
    (group_label == "CLR" & core_annotation == "LA") |
      (group_label == "DII" & core_annotation == "Diffuse")
  )

spe_list_core_filtered <- make_spe_list(core_filtered_df)

diffuse_only_df <- df |>
  filter(core_annotation == "Diffuse")

spe_list_diffuse_only <- make_spe_list(diffuse_only_df)

set.seed(1)
sampled_spots <- core_filtered_df |>
  distinct(patient, spot, group_label, core_annotation) |>
  group_by(patient) |>
  slice_sample(n = 1) |>
  ungroup()

one_per_patient_df <- core_filtered_df |>
  semi_join(sampled_spots, by = c("patient", "spot"))

spe_list_one_per_patient <- make_spe_list(one_per_patient_df)

# ---- Save ----
write_rds(spe_list_all, output_file)
message("Wrote SpatialExperiment list to: ", output_file)

write_rds(spe_list_one_per_patient, output_file_one_per_patient)
message("Wrote SpatialExperiment list to: ", output_file_one_per_patient)

write_rds(spe_list_core_filtered, output_file_core_filtered)
message("Wrote SpatialExperiment list to: ", output_file_core_filtered)

write_rds(spe_list_diffuse_only, output_file_diffuse_only)
message("Wrote SpatialExperiment list to: ", output_file_diffuse_only)
