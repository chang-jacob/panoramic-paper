#!/usr/bin/env Rscript

# CRC TMA Script: Runtime Tables
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Collect runtime metric CSVs from CRC TMA analysis runs.
# - Render per-run runtime LaTeX tables and an aggregate runtime summary table.
# - Write summary CSV/TeX artifacts for manuscript reporting.

suppressPackageStartupMessages({
  library(config)
  library(dplyr)
  library(readr)
})

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#$])", "\\\\\\1", x, perl = TRUE)
  x
}

format_run_scope <- function(path) {
  if (grepl("/radii_stability/", path, fixed = FALSE)) {
    return("Radii stability")
  }
  if (grepl("/one_per_patient/", path, fixed = FALSE)) {
    return("Primary CRC (one per patient)")
  }
  if (grepl("/one_sample_per_patient/", path, fixed = FALSE)) {
    return("Primary CRC (legacy directory)")
  }
  "CRC runtime run"
}

format_stage_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x <- gsub("\\bmeta analysis\\b", "meta-analysis", x, ignore.case = TRUE)
  tools::toTitleCase(x)
}

render_runtime_stage_latex <- function(tbl, caption, label) {
  header <- paste(c("Stage", "Status", "Elapsed (s)", "Elapsed (min)", "Share (\\%)", "Features", "Samples"), collapse = " & ")
  align <- ">{\\raggedright\\arraybackslash}p{2.8cm}llrrrr"

  fmt_num <- function(x, digits = 3L) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
  fmt_int <- function(x) ifelse(is.finite(x), as.character(as.integer(round(x))), "NA")

  if (nrow(tbl) == 0L) {
    row_lines <- "No runtime stages found & NA & NA & NA & NA & NA & NA\\\\"
  } else {
    row_lines <- vapply(seq_len(nrow(tbl)), function(i) {
      vals <- c(
        latex_escape(tbl$Stage[[i]]),
        latex_escape(tbl$Status[[i]]),
        fmt_num(tbl$Elapsed_sec[[i]], digits = 3L),
        fmt_num(tbl$Elapsed_min[[i]], digits = 2L),
        fmt_num(tbl$Share_pct[[i]], digits = 1L),
        fmt_int(tbl$Features[[i]]),
        fmt_int(tbl$Samples[[i]])
      )
      paste0(paste(vals, collapse = " & "), "\\\\")
    }, FUN.VALUE = character(1))
  }

  lines <- c(
    paste0("\\begin{longtable}[t]{", align, "}"),
    paste0("\\caption{", latex_escape(caption), "}"),
    paste0("\\label{", label, "}\\\\"),
    "\\toprule",
    paste0(header, "\\\\"),
    "\\midrule",
    "\\endfirsthead",
    "",
    "\\toprule",
    paste0(header, "\\\\"),
    "\\midrule",
    "\\endhead",
    "",
    "\\midrule",
    "\\multicolumn{7}{r}{Continued on next page}\\\\",
    "\\endfoot",
    "",
    "\\bottomrule",
    "\\endlastfoot",
    "",
    row_lines,
    "\\end{longtable}"
  )
  paste(lines, collapse = "\n")
}

render_runtime_summary_latex <- function(tbl, caption, label) {
  header <- paste(c("Run", "Scope", "Total (min)", "Prepare (min)", "Spatialstats (min)", "Meta-analysis (min)", "Features", "Samples"), collapse = " & ")
  align <- ">{\\raggedright\\arraybackslash}p{4.8cm}>{\\raggedright\\arraybackslash}p{2.8cm}rrrrrr"

  fmt_num <- function(x, digits = 2L) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
  fmt_int <- function(x) ifelse(is.finite(x), as.character(as.integer(round(x))), "NA")

  row_lines <- vapply(seq_len(nrow(tbl)), function(i) {
    vals <- c(
      latex_escape(tbl$run_name[[i]]),
      latex_escape(tbl$scope[[i]]),
      fmt_num(tbl$total_min[[i]]),
      fmt_num(tbl$prepare_min[[i]]),
      fmt_num(tbl$spatialstats_min[[i]]),
      fmt_num(tbl$meta_analysis_min[[i]]),
      fmt_int(tbl$n_features[[i]]),
      fmt_int(tbl$n_samples[[i]])
    )
    paste0(paste(vals, collapse = " & "), "\\\\")
  }, FUN.VALUE = character(1))

  lines <- c(
    paste0("\\begin{longtable}[t]{", align, "}"),
    paste0("\\caption{", latex_escape(caption), "}"),
    paste0("\\label{", label, "}\\\\"),
    "\\toprule",
    paste0(header, "\\\\"),
    "\\midrule",
    "\\endfirsthead",
    "",
    "\\toprule",
    paste0(header, "\\\\"),
    "\\midrule",
    "\\endhead",
    "",
    "\\bottomrule",
    "\\endfoot",
    row_lines,
    "\\end{longtable}"
  )
  paste(lines, collapse = "\n")
}

build_stage_table <- function(df) {
  total_sec <- sum(df$elapsed_sec, na.rm = TRUE)
  safe_max <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    max(x)
  }
  share_pct <- if (is.finite(total_sec) && total_sec > 0) {
    100 * as.numeric(df$elapsed_sec) / total_sec
  } else {
    rep(NA_real_, nrow(df))
  }

  stage_df <- df |>
    transmute(
      Stage = format_stage_label(.data$stage),
      Status = as.character(.data$status),
      Elapsed_sec = as.numeric(.data$elapsed_sec),
      Elapsed_min = as.numeric(.data$elapsed_sec) / 60,
      Share_pct = share_pct,
      Features = as.numeric(.data$n_features),
      Samples = as.numeric(.data$n_samples)
    )

  total_row <- data.frame(
    Stage = "Total",
    Status = "computed",
    Elapsed_sec = total_sec,
    Elapsed_min = total_sec / 60,
    Share_pct = 100,
    Features = safe_max(df$n_features),
    Samples = safe_max(df$n_samples),
    stringsAsFactors = FALSE
  )

  stage_df <- bind_rows(stage_df, total_row)
  stage_df$Features[!is.finite(stage_df$Features)] <- NA_real_
  stage_df$Samples[!is.finite(stage_df$Samples)] <- NA_real_
  stage_df
}

config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)

crc_tab_root <- file.path(paths$output, "tables", "crc_tma")
runtime_files <- list.files(
  path = crc_tab_root,
  pattern = "runtime_metrics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(runtime_files) == 0L) {
  stop("No CRC TMA runtime metrics CSV files found under: ", crc_tab_root)
}

runtime_summary_rows <- vector("list", length(runtime_files))

for (i in seq_along(runtime_files)) {
  f <- runtime_files[[i]]
  run_df <- readr::read_csv(f, show_col_types = FALSE)
  if (nrow(run_df) == 0L) next

  run_name <- sub("_runtime_metrics\\.csv$", "", basename(f))
  scope <- format_run_scope(f)
  stage_tbl <- build_stage_table(run_df)

  out_tex <- file.path(dirname(f), paste0(run_name, "_runtime_table.tex"))
  caption <- paste0("Runtime breakdown for ", scope, " (", run_name, ").")
  label <- paste0("tab:s_", gsub("[^a-z0-9]+", "_", tolower(run_name)), "_runtime")
  latex_txt <- render_runtime_stage_latex(stage_tbl, caption = caption, label = label)
  writeLines(latex_txt, con = out_tex)

  prep_sec <- sum(run_df$elapsed_sec[run_df$stage == "prepare"], na.rm = TRUE)
  spatial_sec <- sum(run_df$elapsed_sec[run_df$stage == "spatialstats"], na.rm = TRUE)
  meta_sec <- sum(run_df$elapsed_sec[run_df$stage == "meta_analysis"], na.rm = TRUE)
  total_sec <- sum(run_df$elapsed_sec, na.rm = TRUE)
  n_features_vals <- as.numeric(run_df$n_features)
  n_features_vals <- n_features_vals[is.finite(n_features_vals)]
  n_samples_vals <- as.numeric(run_df$n_samples)
  n_samples_vals <- n_samples_vals[is.finite(n_samples_vals)]

  runtime_summary_rows[[i]] <- data.frame(
    run_name = run_name,
    scope = scope,
    total_min = total_sec / 60,
    prepare_min = prep_sec / 60,
    spatialstats_min = spatial_sec / 60,
    meta_analysis_min = meta_sec / 60,
    n_features = if (length(n_features_vals) > 0L) max(n_features_vals) else NA_real_,
    n_samples = if (length(n_samples_vals) > 0L) max(n_samples_vals) else NA_real_,
    runtime_csv = f,
    runtime_tex = out_tex,
    stringsAsFactors = FALSE
  )

  message("Runtime LaTeX table: ", out_tex)
}

summary_tbl <- bind_rows(runtime_summary_rows) |>
  mutate(
    n_features = ifelse(is.finite(.data$n_features), .data$n_features, NA_real_),
    n_samples = ifelse(is.finite(.data$n_samples), .data$n_samples, NA_real_)
  ) |>
  arrange(scope, run_name)

summary_dir <- file.path(crc_tab_root, "runtime_tables")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(summary_dir, "crc_tma_runtime_summary.csv")
summary_tex <- file.path(summary_dir, "crc_tma_runtime_summary_table.tex")
readr::write_csv(summary_tbl, summary_csv)

summary_latex <- render_runtime_summary_latex(
  summary_tbl,
  caption = "CRC TMA runtime summary across available runs.",
  label = "tab:s_crc_tma_runtime_summary"
)
writeLines(summary_latex, con = summary_tex)

message("Runtime summary CSV: ", summary_csv)
message("Runtime summary LaTeX table: ", summary_tex)
