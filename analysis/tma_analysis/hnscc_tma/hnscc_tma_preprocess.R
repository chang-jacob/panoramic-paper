#!/usr/bin/env Rscript

# HNSCC TMA Script: Preprocess
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Load raw HNSCC TMA cell-level and clinical metadata tables.
# - Standardize columns/labels and define comparison subsets for downstream analyses.
# - Write processed HNSCC TMA CSV outputs used by analysis and EDA scripts.

suppressPackageStartupMessages({
  library(config)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ---- Load paths ----
config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)

raw_dir <- file.path(paths$data, "hnscc_tma")
processed_dir <- file.path(paths$data, "processed", "hnscc_tma")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

codex_file <- file.path(raw_dir, "HNSCCdiscovery_CODEX.csv")
clin_file <- file.path(raw_dir, "HNSCCdiscovery_clinmetadata.csv")
processed_file <- file.path(processed_dir, "processed_hnscc_tma.csv")
processed_main_comparison_file <- file.path(
  processed_dir,
  "processed_hnscc_tma_LNbenign_vs_LNmetadj.csv"
)

if (!file.exists(codex_file)) stop("Missing CODEX file: ", codex_file)
if (!file.exists(clin_file)) stop("Missing clinical file: ", clin_file)

# ---- Load ----
fast_read_selected_csv <- function(path, select_cols) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    # fread() is substantially faster on large CSVs and supports column projection.
    dt <- data.table::fread(
      input = path,
      select = select_cols,
      data.table = FALSE,
      showProgress = TRUE,
      check.names = FALSE
    )
    return(tibble::as_tibble(dt))
  }

  message("Package 'data.table' not found; falling back to readr::read_csv().")
  readr::read_csv(
    file = path,
    col_select = all_of(select_cols),
    show_col_types = FALSE
  )
}

id_cols <- c(
  "UID",
  "unique_region",
  "x",
  "y",
  "Xcorr",
  "Ycorr",
  "Pat_ID",
  "Tissue_type_broad",
  "cell_types_lowres",
  "cell_types_mediumres",
  "neighb_name_new",
  "neighb_name_short"
)

codex_meta <- fast_read_selected_csv(
  path = codex_file,
  select_cols = id_cols
)

# codex <- read_csv(codex_file, show_col_types = FALSE)
clin <- read_csv(clin_file, show_col_types = FALSE)

# ---- Basic validation ----
required_codex <- c("UID", "unique_region", "x", "y", "Pat_ID",
                    "cell_types_lowres", "cell_types_mediumres")
missing_codex <- setdiff(required_codex, colnames(codex_meta))
if (length(missing_codex) > 0) {
  stop("CODEX file missing required columns: ", paste(missing_codex, collapse = ", "))
}

required_clin <- c("CaseID", "unique_region")
missing_clin <- setdiff(required_clin, colnames(clin))
if (length(missing_clin) > 0) {
  stop("Clinical file missing required columns: ", paste(missing_clin, collapse = ", "))
}

# ---- Join clinical metadata ----
df <- codex_meta %>%
  left_join(clin, by = "unique_region", suffix = c("", "_clin"))

# ---- Standardize columns ----
df <- df %>%
  mutate(
    cell_id = UID,
    region = unique_region,
    patient = Pat_ID,
    cell_type = as.character(cell_types_mediumres),
    cell_type_lowres = as.character(cell_types_lowres),
    group = if ("Nodal_disease" %in% colnames(df)) as.character(Nodal_disease) else NA_character_, 
    core_category = Tissue_type_broad
  )

MICRONS_PER_PIXEL <- 0.37742
df <- df |>
  mutate(
    x = x * MICRONS_PER_PIXEL,
    y = y * MICRONS_PER_PIXEL
  ) 
# ---- Basic validation (processed) ----
required_processed <- c("cell_id", "cell_type", "x", "y", "patient", "region")
missing_processed <- setdiff(required_processed, colnames(df))
if (length(missing_processed) > 0) {
  stop("Processed data missing required columns: ", paste(missing_processed, collapse = ", "))
}

# > unique(df$Tissue_type_broad)
# [1] "PT"              "LN_met"          "LN_benign"       "LN_met_adjacent"

pt_file <- file.path(processed_dir, "processed_hnscc_tma_PT.csv")
ln_met_file <- file.path(processed_dir, "processed_hnscc_tma_LNmet.csv")
ln_benign_file <- file.path(processed_dir, "processed_hnscc_tma_LNbenign.csv")
ln_metadj_file <- file.path(processed_dir, "processed_hnscc_tma_LNmetadj.csv")

write_csv(df %>% filter(core_category == "PT"), pt_file)
write_csv(df %>% filter(core_category == "LN_met"), ln_met_file)
write_csv(df %>% filter(core_category == "LN_benign"), ln_benign_file)
write_csv(df %>% filter(core_category == "LN_met_adjacent"), ln_metadj_file)

# Main comparison subset: LN_benign vs LN_met_adjacent
df_main_comparison <- df %>%
  filter(core_category %in% c("LN_benign", "LN_met_adjacent"))
write_csv(df_main_comparison, processed_main_comparison_file)
message("Wrote main comparison subset to: ", processed_main_comparison_file)

pt_comparison_file <- pt_file
df_pt_comparison <- df %>%
  filter(core_category %in% c("PT"))
write_csv(df_pt_comparison, pt_comparison_file)
message("Wrote PT subset to: ", pt_comparison_file)

# ---- Save ----
write_csv(df, processed_file)
message("Wrote processed data to: ", processed_file)


foo <- df %>% 
  group_by(region) %>% 
  slice(1)
table(foo$group, foo$core_category)  

foo2 <- df %>% 
  group_by(patient) %>% 
  slice(1)
table(foo2$group, foo2$core_category)  
  

df |>
  distinct(Pat_ID, group, core_category) |>
  count(group, core_category, name = "n_unique_patients") |>
  tidyr::pivot_wider(
    names_from = core_category,
    values_from = n_unique_patients,
    values_fill = 0
  )


# TODO: Need to find the conversion of these pixel coordinates to micron!
