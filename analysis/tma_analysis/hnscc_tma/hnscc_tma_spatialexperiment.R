#!/usr/bin/env Rscript

# HNSCC TMA Script: Build SpatialExperiment Lists
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Read processed HNSCC TMA data and validate required spatial/annotation columns.
# - Build region-level SpatialExperiment objects for PANORAMIC workflows.
# - Save HNSCC SpatialExperiment list RDS files for downstream analysis scripts.

suppressPackageStartupMessages({
  library(config)
  library(dplyr)
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

processed_dir <- file.path(paths$data, "processed", "hnscc_tma")
required_cols <- c("cell_id", "cell_type", "x", "y", "patient", "region")

# ---- Build helpers ----
build_spe_list <- function(df) {
  lapply(split(df, df$region), function(subdf) {
    coords <- as.matrix(subdf[, c("x", "y")])
    rownames(coords) <- subdf$cell_id

    col_data <- S4Vectors::DataFrame(subdf)
    rownames(col_data) <- subdf$cell_id

    counts <- matrix(
      0,
      nrow = 1,
      ncol = nrow(subdf),
      dimnames = list("placeholder", subdf$cell_id)
    )

    SpatialExperiment(
      assays = list(counts = counts),
      colData = col_data,
      spatialCoords = coords
    )
  })
}

write_spe_list_from_processed <- function(input_basename, output_basename, alias_basenames = character(0)) {
  input_file <- file.path(processed_dir, input_basename)
  output_file <- file.path(processed_dir, output_basename)

  if (!file.exists(input_file)) {
    stop("Missing processed input file: ", input_file)
  }

  df <- read_csv(input_file, show_col_types = FALSE)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Processed data missing required columns in ", input_basename, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  spe_list <- build_spe_list(df)
  write_rds(spe_list, output_file)
  message("Wrote SpatialExperiment list to: ", output_file)

  if (length(alias_basenames) > 0L) {
    for (alias_basename in alias_basenames) {
      alias_file <- file.path(processed_dir, alias_basename)
      if (!identical(alias_file, output_file)) {
        write_rds(spe_list, alias_file)
        message("Wrote alias SpatialExperiment list to: ", alias_file)
      }
    }
  }
}

write_pt_one_per_patient_spe_list <- function(input_basename, output_basename, alias_basenames = character(0)) {
  input_file <- file.path(processed_dir, input_basename)
  output_file <- file.path(processed_dir, output_basename)

  if (!file.exists(input_file)) {
    stop("Missing processed input file: ", input_file)
  }

  df <- read_csv(input_file, show_col_types = FALSE)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Processed data missing required columns in ", input_basename, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Keep one region per patient: choose region with most cells; break ties by region name.
  chosen_regions <- df |>
    dplyr::count(patient, region, name = "n_cells") |>
    dplyr::arrange(patient, dplyr::desc(n_cells), region) |>
    dplyr::group_by(patient) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::select(patient, region)

  df_one_per_patient <- df |>
    dplyr::semi_join(chosen_regions, by = c("patient", "region"))

  n_patients <- dplyr::n_distinct(df$patient)
  n_regions_total <- dplyr::n_distinct(df$region)
  n_regions_kept <- dplyr::n_distinct(df_one_per_patient$region)

  if (n_regions_kept != n_patients) {
    stop(
      "Expected one retained region per patient, but got ", n_regions_kept,
      " regions for ", n_patients, " patients."
    )
  }

  message(
    "PT one-per-patient selection: kept ", n_regions_kept, " regions across ",
    n_patients, " patients (from ", n_regions_total, " total regions)."
  )

  spe_list <- build_spe_list(df_one_per_patient)
  write_rds(spe_list, output_file)
  message("Wrote one-per-patient SpatialExperiment list to: ", output_file)

  if (length(alias_basenames) > 0L) {
    for (alias_basename in alias_basenames) {
      alias_file <- file.path(processed_dir, alias_basename)
      if (!identical(alias_file, output_file)) {
        write_rds(spe_list, alias_file)
        message("Wrote alias one-per-patient SpatialExperiment list to: ", alias_file)
      }
    }
  }
}

write_lnmix_nodal_matched_spe_list <- function(input_basename, output_basename, alias_basenames = character(0)) {
  input_file <- file.path(processed_dir, input_basename)
  output_file <- file.path(processed_dir, output_basename)

  if (!file.exists(input_file)) {
    stop("Missing processed input file: ", input_file)
  }

  df <- read_csv(input_file, show_col_types = FALSE)
  required_cols_local <- c(required_cols, "group", "core_category")
  missing_cols <- setdiff(required_cols_local, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Processed data missing required columns in ", input_basename, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df <- df |>
    dplyr::mutate(
      group = trimws(as.character(group)),
      core_category = as.character(core_category)
    )

  keep_idx <- (df$core_category == "LN_benign" & df$group %in% c("N-", "Nneg")) |
    (df$core_category == "LN_met_adjacent" & df$group %in% c("N+", "Npos"))

  df_matched <- df[keep_idx, , drop = FALSE]
  if (nrow(df_matched) == 0L) {
    stop("No rows retained for LN_benign/N- and LN_met_adjacent/N+ matched subset.")
  }

  n_regions <- dplyr::n_distinct(df_matched$region)
  n_patients <- dplyr::n_distinct(df_matched$patient)
  message(
    "Matched nodal LN subset: ", nrow(df_matched), " cells across ",
    n_regions, " regions and ", n_patients, " patients."
  )

  spe_list <- build_spe_list(df_matched)
  write_rds(spe_list, output_file)
  message("Wrote matched nodal LN SpatialExperiment list to: ", output_file)

  if (length(alias_basenames) > 0L) {
    for (alias_basename in alias_basenames) {
      alias_file <- file.path(processed_dir, alias_basename)
      if (!identical(alias_file, output_file)) {
        write_rds(spe_list, alias_file)
        message("Wrote alias matched nodal LN SpatialExperiment list to: ", alias_file)
      }
    }
  }
}

# ---- Build requested lists ----
write_spe_list_from_processed(
  input_basename = "processed_hnscc_tma_LNbenign_vs_LNmetadj.csv",
  output_basename = "hnscc_tma_LNbenign_vs_LNmetadj_spe_list.rds",
  alias_basenames = c("hnscc_tma_spe_list.rds", "hnscc_tma_spe_list_main_comparison.rds")
)

write_spe_list_from_processed(
  input_basename = "processed_hnscc_tma_PT.csv",
  output_basename = "hnscc_tma_PT_spe_list.rds",
  alias_basenames = c(
    "hnscc_tma_spe_list_pt.rds",
    "hnscc_tma_primary_tumors_spe_list.rds"
  )
)

write_pt_one_per_patient_spe_list(
  input_basename = "processed_hnscc_tma_PT.csv",
  output_basename = "hnscc_tma_PT_one_per_patient_spe_list.rds",
  alias_basenames = c("hnscc_tma_primary_tumors_one_per_patient_spe_list.rds")
)

write_lnmix_nodal_matched_spe_list(
  input_basename = "processed_hnscc_tma_LNbenign_vs_LNmetadj.csv",
  output_basename = "hnscc_tma_LNbenign_Nneg_and_LNmetadj_Npos_spe_list.rds"
)
