#!/usr/bin/env Rscript

# CRC TMA Script: Benchmark (PANORAMIC vs t-test vs spicyR)
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Compare PANORAMIC contrast outputs against per-feature Welch t-tests and spicyR.
# - Build directional/pair-level benchmark merges and correlation summaries.
# - Export benchmark figures and PANORAMIC-style LaTeX tables for supplement use.

suppressPackageStartupMessages({
  library(config)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(BiocParallel)
  library(spicyR)
})

# ---- Analysis Parameters (edit in script) ----
analysis_params <- list(
  spe_list_key = "one_per_patient", # all | one_per_patient | core_filtered | diffuse_only
  columns = list(
    cell_type = "cell_type",
    patient = "patient",
    group = "group_label",
    region = "spot"
  ),
  # Explicit comparison order for effect_ttest = mean(group2) - mean(group1).
  group_order = c("CLR", "DII"),
  panoramic = list(
    radii_um = c(25),
    stat = "local_comp_enrichment",
    tile_size = 62.5,
    meta_engine = "mv" # mv | re
  ),
  ttest = list(
    min_samples_per_group = 2L
  ),
  spicy = list(
    workers = 10L,
    rs_step_um = 0.5,
    include_zero_cells = FALSE,
    subject_col = "patient", # set to NULL explicitly if no subject term is desired
    top_n = 10000L,
    # Match legacy spicyR paper-style filtering for CRC.
    drop_zero_fraction_gt = 0.21,
    drop_cell_types = c("undefined", "dirt")
  ),
  reporting = list(
    write_latex_results = TRUE,
    latex_sig_only = TRUE,
    latex_alpha = 0.05,
    latex_print_console = TRUE
  )
)

# ---- Helpers ----
parse_radii <- function(x) {
  if (length(x) == 1L && is.character(x)) {
    y <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
  } else {
    y <- as.numeric(x)
  }
  y <- y[is.finite(y) & y > 0]
  sort(unique(y))
}

pair_key_dir <- function(ct1, ct2) {
  paste(as.character(ct1), as.character(ct2), sep = "||")
}

pair_key_undir <- function(ct1, ct2) {
  a <- as.character(ct1)
  b <- as.character(ct2)
  lo <- ifelse(a <= b, a, b)
  hi <- ifelse(a <= b, b, a)
  paste(lo, hi, sep = "||")
}

get_unique <- function(spe, col) {
  vals <- unique(as.character(SummarizedExperiment::colData(spe)[[col]]))
  vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
  if (length(vals) != 1L) {
    stop("Expected exactly one unique value in colData[['", col, "']] per sample.")
  }
  vals[[1]]
}

get_sample_group_value <- function(spe, col, sample_id = NA_character_) {
  vals_raw <- as.character(SummarizedExperiment::colData(spe)[[col]])
  vals <- vals_raw[!is.na(vals_raw) & nzchar(trimws(vals_raw))]
  vals <- unique(vals)
  if (length(vals) != 1L) {
    stop("Expected exactly one unique group in colData[['", col, "']] for sample ", sample_id, ".")
  }
  vals[[1]]
}

compute_ttest_table <- function(spatial_tbl, group_col, group_order, min_samples_per_group = 2L) {
  dat <- spatial_tbl |>
    dplyr::mutate(
      group_cmp = as.character(.data[[group_col]]),
      yi = as.numeric(.data$yi)
    ) |>
    dplyr::filter(group_cmp %in% group_order, is.finite(yi))

  split_idx <- interaction(dat$ct1, dat$ct2, dat$radius_um, dat$stat, drop = TRUE)
  pieces <- split(dat, split_idx)

  out <- vector("list", length(pieces))
  i <- 0L
  for (df in pieces) {
    i <- i + 1L
    x <- df$yi[df$group_cmp == group_order[[1]]]
    y <- df$yi[df$group_cmp == group_order[[2]]]
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]

    effect <- if (length(x) > 0L && length(y) > 0L) mean(y) - mean(x) else NA_real_
    p <- if (length(x) >= min_samples_per_group && length(y) >= min_samples_per_group) {
      stats::t.test(x = x, y = y, alternative = "two.sided", var.equal = FALSE)$p.value
    } else {
      NA_real_
    }

    out[[i]] <- data.frame(
      ct1 = as.character(df$ct1[[1]]),
      ct2 = as.character(df$ct2[[1]]),
      radius_um = as.numeric(df$radius_um[[1]]),
      stat = as.character(df$stat[[1]]),
      n_control = length(x),
      n_case = length(y),
      mean_control = if (length(x) > 0L) mean(x) else NA_real_,
      mean_case = if (length(y) > 0L) mean(y) else NA_real_,
      effect_ttest = effect,
      p_ttest = p,
      status = "ok",
      stringsAsFactors = FALSE
    )
  }

  ttest_tbl <- dplyr::bind_rows(out)
  ttest_tbl$fdr_ttest <- p.adjust(ttest_tbl$p_ttest, method = "fdr")
  ttest_tbl$pair_key_dir <- pair_key_dir(ttest_tbl$ct1, ttest_tbl$ct2)
  ttest_tbl$pair_key_undir <- pair_key_undir(ttest_tbl$ct1, ttest_tbl$ct2)
  ttest_tbl
}

build_sample_meta <- function(spe_list, patient_col, group_col, region_col) {
  rows <- lapply(seq_along(spe_list), function(i) {
    spe <- spe_list[[i]]
    region_val <- get_unique(spe, region_col)
    data.frame(
      sample = as.character(region_val),
      patient = get_unique(spe, patient_col),
      group = get_sample_group_value(spe, group_col, sample_id = as.character(region_val)),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

build_pair_level <- function(tbl, p_col, fdr_col, effect_col) {
  tbl <- tbl |>
    dplyr::mutate(
      p_tmp = as.numeric(.data[[p_col]]),
      fdr_tmp = as.numeric(.data[[fdr_col]]),
      pair_key_undir = pair_key_undir(.data$ct1, .data$ct2)
    )
  tbl |>
    dplyr::group_by(pair_key_undir) |>
    dplyr::arrange(p_tmp, .by_group = TRUE) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      pair_key_undir = .data$pair_key_undir,
      ct1_rep = .data$ct1,
      ct2_rep = .data$ct2,
      !!p_col := .data$p_tmp,
      !!fdr_col := .data$fdr_tmp,
      !!effect_col := as.numeric(.data[[effect_col]])
    )
}

build_method_latex_meta <- function(
    tbl,
    group_order,
    effect_col,
    p_col,
    fdr_col,
    group1_col = NULL,
    group2_col = NULL
) {
  group1_name <- paste0(group_order[[1]], "_mu_hat")
  group2_name <- paste0(group_order[[2]], "_mu_hat")

  group1_vals <- if (is.null(group1_col)) {
    rep(NA_real_, nrow(tbl))
  } else {
    as.numeric(tbl[[group1_col]])
  }
  group2_vals <- if (is.null(group2_col)) {
    rep(NA_real_, nrow(tbl))
  } else {
    as.numeric(tbl[[group2_col]])
  }

  out <- data.frame(
    coloc_direction = paste0(as.character(tbl$ct2), " -> ", as.character(tbl$ct1)),
    beta_diff = as.numeric(tbl[[effect_col]]),
    p_diff = as.numeric(tbl[[p_col]]),
    fdr_diff = as.numeric(tbl[[fdr_col]]),
    stringsAsFactors = FALSE
  )
  out[[group1_name]] <- group1_vals
  out[[group2_name]] <- group2_vals
  out[, c("coloc_direction", group1_name, group2_name, "beta_diff", "p_diff", "fdr_diff"), drop = FALSE]
}

write_scatter <- function(df, x_col, y_col, x_lab, y_lab, title, out_stub) {
  x <- as.numeric(df[[x_col]])
  y <- as.numeric(df[[y_col]])
  keep <- is.finite(x) & is.finite(y)
  plot_df <- data.frame(x = x[keep], y = y[keep])
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(alpha = 0.7, size = 1.6, color = "#2C7FB8") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray45") +
    theme_bw(base_size = 11) +
    labs(title = title, x = x_lab, y = y_lab)
  ggsave(paste0(out_stub, ".png"), p, width = 6.5, height = 5, dpi = 300)
  ggsave(paste0(out_stub, ".pdf"), p, width = 6.5, height = 5)
}

write_pvalue_hist <- function(df, p_col, title, out_stub, bins = 20L, alpha_line = 0.05) {
  p_vals <- as.numeric(df[[p_col]])
  keep <- is.finite(p_vals) & p_vals >= 0 & p_vals <= 1
  hist_df <- data.frame(p_value = p_vals[keep])
  p <- ggplot(hist_df, aes(x = p_value)) +
    geom_histogram(bins = as.integer(bins), fill = "#2C7FB8", color = "white") +
    geom_vline(xintercept = alpha_line, color = "red", linetype = "dashed") +
    labs(title = title, x = "p-value", y = "Count") +
    theme_bw(base_size = 11)

  out_png <- paste0(out_stub, ".png")
  out_pdf <- paste0(out_stub, ".pdf")
  ggsave(out_png, p, width = 7, height = 5, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 5)
  invisible(list(png = out_png, pdf = out_pdf))
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#$])", "\\\\\\1", x, perl = TRUE)
  x
}

format_cell_pair_latex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s*->\\s*", " XXARROWTOKENXX ", x)
  x <- latex_escape(x)
  x <- gsub("XXARROWTOKENXX", "$\\\\to$", x, fixed = TRUE)
  x
}

build_pan_vs_method_results_table <- function(
    pair_merge,
    method = c("ttest", "spicyr"),
    sig_only = TRUE,
    alpha = 0.05
) {
  method <- match.arg(method)
  method_cols <- switch(
    method,
    ttest = c("effect_ttest", "p_ttest", "fdr_ttest"),
    spicyr = c("effect_spicyr", "p_spicyr", "fdr_spicyr")
  )
  method_label <- if (identical(method, "ttest")) "t-test" else "spicyR"

  tbl <- pair_merge
  pair_label <- ifelse(
    is.na(tbl$ct1_rep) | is.na(tbl$ct2_rep),
    gsub("\\|\\|", " -> ", as.character(tbl$pair_key_undir)),
    paste0(tbl$ct1_rep, " -> ", tbl$ct2_rep)
  )

  out <- data.frame(
    "Cell type pair" = pair_label,
    "PANORAMIC effect size" = as.numeric(tbl$effect_panoramic),
    "PANORAMIC p-value" = as.numeric(tbl$p_panoramic),
    "PANORAMIC adjusted p-value" = as.numeric(tbl$fdr_panoramic),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out[[paste0(method_label, " effect size")]] <- as.numeric(tbl[[method_cols[[1]]]])
  out[[paste0(method_label, " p-value")]] <- as.numeric(tbl[[method_cols[[2]]]])
  out[[paste0(method_label, " adjusted p-value")]] <- as.numeric(tbl[[method_cols[[3]]]])

  ord <- order(
    out[["PANORAMIC adjusted p-value"]],
    out[[paste0(method_label, " adjusted p-value")]],
    out[["PANORAMIC p-value"]],
    out[[paste0(method_label, " p-value")]],
    na.last = TRUE
  )
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  if (isTRUE(sig_only)) {
    keep_pan <- is.finite(out[["PANORAMIC adjusted p-value"]]) &
      out[["PANORAMIC adjusted p-value"]] <= alpha
    keep_method <- is.finite(out[[paste0(method_label, " adjusted p-value")]]) &
      out[[paste0(method_label, " adjusted p-value")]] <= alpha
    out <- out[keep_pan | keep_method, , drop = FALSE]
    rownames(out) <- NULL
  }

  out
}

render_pan_vs_method_results_latex <- function(
    tbl,
    caption,
    label,
    first_col_width = "3.2cm",
    continued_text = "Continued on next page"
) {
  col_align <- paste0(">{\\raggedright\\arraybackslash}p{", first_col_width, "}rrrrrr")
  header <- paste(latex_escape(names(tbl)), collapse = " & ")

  fmt_est <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 4), "NA")
  fmt_p <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 3), "NA")

  if (nrow(tbl) == 0L) {
    row_lines <- "No rows meeting selection criteria & NA & NA & NA & NA & NA & NA\\\\"
  } else {
    row_lines <- vapply(seq_len(nrow(tbl)), function(i) {
      vals <- c(
        format_cell_pair_latex(tbl[[1]][i]),
        fmt_est(tbl[[2]][i]),
        fmt_p(tbl[[3]][i]),
        fmt_p(tbl[[4]][i]),
        fmt_est(tbl[[5]][i]),
        fmt_p(tbl[[6]][i]),
        fmt_p(tbl[[7]][i])
      )
      paste0(paste(vals, collapse = " & "), "\\\\")
    }, FUN.VALUE = character(1))
  }

  lines <- c(
    paste0("\\begin{longtable}[t]{", col_align, "}"),
    paste0("\\caption{", latex_escape(caption), "}"),
    if (!is.null(label) && nzchar(as.character(label))) paste0("\\label{", as.character(label), "}\\\\") else "\\\\",
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
    paste0("\\multicolumn{7}{r}{", latex_escape(continued_text), "}\\\\"),
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

write_pan_vs_method_latex_results <- function(
    pair_merge,
    method = c("ttest", "spicyr"),
    out_path,
    sig_only = TRUE,
    alpha = 0.05,
    caption,
    label,
    first_col_width = "3.2cm",
    print_console = FALSE
) {
  method <- match.arg(method)
  tbl <- build_pan_vs_method_results_table(
    pair_merge = pair_merge,
    method = method,
    sig_only = sig_only,
    alpha = alpha
  )
  latex_txt <- render_pan_vs_method_results_latex(
    tbl = tbl,
    caption = caption,
    label = label,
    first_col_width = first_col_width
  )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(latex_txt, con = out_path)

  if (isTRUE(print_console)) {
    cat(latex_txt, sep = "\n")
    cat("\n")
  }

  invisible(list(path = out_path, n_rows = nrow(tbl), latex = latex_txt, table = tbl))
}

render_benchmark_summary_latex <- function(
    summary_tbl,
    caption = "CRC TMA benchmark summary metrics.",
    label = "tab:s_crc_tma_benchmark_summary",
    metric_col_width = "5.0cm",
    value_col_width = "8.0cm"
) {
  tbl <- summary_tbl |>
    dplyr::transmute(
      Metric = gsub("_", " ", as.character(.data$metric)),
      Value = as.character(.data$value)
    )

  header <- paste(latex_escape(names(tbl)), collapse = " & ")
  align <- paste0(
    ">{\\raggedright\\arraybackslash}p{", metric_col_width, "}",
    ">{\\raggedright\\arraybackslash}p{", value_col_width, "}"
  )

  row_lines <- vapply(seq_len(nrow(tbl)), function(i) {
    paste0(
      latex_escape(tbl$Metric[[i]]), " & ",
      latex_escape(tbl$Value[[i]]), "\\\\"
    )
  }, FUN.VALUE = character(1))

  lines <- c(
    paste0("\\begin{longtable}[t]{", align, "}"),
    paste0("\\caption{", latex_escape(caption), "}"),
    if (!is.null(label) && nzchar(as.character(label))) paste0("\\label{", as.character(label), "}\\\\") else "\\\\",
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

write_benchmark_summary_latex <- function(
    summary_tbl,
    out_path,
    caption = "CRC TMA benchmark summary metrics.",
    label = "tab:s_crc_tma_benchmark_summary",
    print_console = FALSE
) {
  latex_txt <- render_benchmark_summary_latex(
    summary_tbl = summary_tbl,
    caption = caption,
    label = label
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(latex_txt, con = out_path)

  if (isTRUE(print_console)) {
    cat(latex_txt, sep = "\n")
    cat("\n")
  }

  invisible(list(path = out_path, n_rows = nrow(summary_tbl), latex = latex_txt))
}

# ---- Run ----
project_root <- getwd()
config_file <- file.path(project_root, "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)

reporting_params <- analysis_params$reporting
write_latex_results <- isTRUE(reporting_params$write_latex_results)
latex_sig_only <- isTRUE(reporting_params$latex_sig_only)
latex_alpha <- as.numeric(reporting_params$latex_alpha)
latex_print_console <- isTRUE(reporting_params$latex_print_console)

shared_helpers <- new.env(parent = baseenv())
sys.source(file.path(project_root, "analysis", "tma_analysis", "tma_shared_helpers.R"), envir = shared_helpers)

radii_um <- parse_radii(analysis_params$panoramic$radii_um)
stat <- as.character(analysis_params$panoramic$stat)
meta_engine <- as.character(analysis_params$panoramic$meta_engine)

tile_size_raw <- analysis_params$panoramic$tile_size
tile_size <- as.numeric(tile_size_raw)
if (length(tile_size) != 1L || is.na(tile_size) || !is.finite(tile_size) || tile_size <= 0) {
  stop("analysis_params$panoramic$tile_size must be one positive numeric value.")
}

spe_list_key <- as.character(analysis_params$spe_list_key)
spe_list_paths <- list(
  all = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list.rds"),
  one_per_patient = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_one_per_patient.rds"),
  core_filtered = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_core_filtered.rds"),
  diffuse_only = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_diffuse_only.rds")
)
list_output_dirs <- c(
  one_per_patient = "one_sample_per_patient",
  core_filtered = "la_vs_diffuse",
  diffuse_only = "all_diffuse",
  all = "all"
)
spe_list_path <- spe_list_paths[[spe_list_key]]
list_dir <- list_output_dirs[[spe_list_key]]
if (is.null(spe_list_path) || is.null(list_dir)) {
  stop("Invalid analysis_params$spe_list_key: ", spe_list_key)
}
if (!file.exists(spe_list_path)) {
  stop("Missing SPE list: ", spe_list_path)
}

if (length(radii_um) == 0L) {
  stop("analysis_params$panoramic$radii_um must contain at least one positive radius.")
}
if (!identical(meta_engine, "mv")) {
  stop("analysis_params$panoramic$meta_engine must be 'mv'.")
}

group_order <- as.character(analysis_params$group_order[seq_len(2L)])
if (length(group_order) != 2L || any(!nzchar(group_order))) {
  stop("analysis_params$group_order must contain exactly two non-empty group labels.")
}

tab_dir <- file.path(paths$output, "tables", "crc_tma", list_dir)
fig_dir <- file.path(paths$output, "figures", "crc_tma", list_dir)
benchmark_tab_dir <- file.path(tab_dir, "benchmark")
benchmark_fig_dir <- file.path(fig_dir, "benchmark")
dir.create(benchmark_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(benchmark_fig_dir, recursive = TRUE, showWarnings = FALSE)

stat_tag <- gsub("[^a-z0-9]+", "_", tolower(stat))
radii_tag <- paste(radii_um, collapse = "_")
tile_tag <- format(signif(tile_size, 6), scientific = FALSE, trim = TRUE)
meta_tag <- "multilevel_mv"
run_prefix <- paste0(
  "crc_tma_", spe_list_key,
  "_", meta_tag,
  "_stat", stat_tag,
  "_r", radii_tag,
  "_t", tile_tag
)

contrast_csv <- file.path(tab_dir, paste0(run_prefix, "_contrast.csv"))
spatial_csv <- file.path(tab_dir, paste0(run_prefix, "_spatialstats.csv"))
meta_csv <- file.path(tab_dir, paste0(run_prefix, "_meta.csv"))
if (!file.exists(contrast_csv)) stop("Missing PANORAMIC contrast CSV: ", contrast_csv)
if (!file.exists(spatial_csv)) stop("Missing PANORAMIC spatialstats CSV: ", spatial_csv)
if (isTRUE(write_latex_results) && !file.exists(meta_csv)) stop("Missing PANORAMIC meta CSV: ", meta_csv)

message("Reading PANORAMIC contrast: ", contrast_csv)
message("Reading PANORAMIC spatial stats: ", spatial_csv)
contrast_tbl <- readr::read_csv(contrast_csv, show_col_types = FALSE)
spatial_tbl <- readr::read_csv(spatial_csv, show_col_types = FALSE)
meta_tbl <- if (isTRUE(write_latex_results)) readr::read_csv(meta_csv, show_col_types = FALSE) else NULL

pan_tbl <- contrast_tbl |>
  dplyr::transmute(
    ct1 = as.character(.data$ct1),
    ct2 = as.character(.data$ct2),
    radius_um = as.numeric(.data$radius_um),
    stat = as.character(.data$stat),
    beta_panoramic = as.numeric(.data$beta_diff),
    p_panoramic = as.numeric(.data$p_diff),
    fdr_panoramic = as.numeric(.data$fdr_diff)
  )
pan_tbl$pair_key_dir <- pair_key_dir(pan_tbl$ct1, pan_tbl$ct2)
pan_tbl$pair_key_undir <- pair_key_undir(pan_tbl$ct1, pan_tbl$ct2)

ttest_min_n <- as.integer(analysis_params$ttest$min_samples_per_group)
ttest_tbl <- compute_ttest_table(
  spatial_tbl = spatial_tbl,
  group_col = analysis_params$columns$group,
  group_order = group_order,
  min_samples_per_group = ttest_min_n
)

spe_list <- readRDS(spe_list_path)
sample_meta <- build_sample_meta(
  spe_list = spe_list,
  patient_col = analysis_params$columns$patient,
  group_col = analysis_params$columns$group,
  region_col = analysis_params$columns$region
)
sample_meta <- sample_meta |>
  dplyr::filter(group %in% group_order)
keep_samples <- intersect(sample_meta$sample, unique(as.character(spatial_tbl$sample)))
sample_meta <- sample_meta |>
  dplyr::filter(sample %in% keep_samples)

spe_by_sample <- spe_list
names(spe_by_sample) <- vapply(spe_by_sample, function(spe) {
  get_unique(spe, analysis_params$columns$region)
}, FUN.VALUE = character(1))
spe_by_sample <- spe_by_sample[sample_meta$sample]

cell_rows <- vector("list", length(spe_by_sample))
i <- 0L
for (sample_id in names(spe_by_sample)) {
  i <- i + 1L
  spe <- spe_by_sample[[sample_id]]
  coords <- SpatialExperiment::spatialCoords(spe)
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  ct_col <- analysis_params$columns$cell_type
  cell_rows[[i]] <- data.frame(
    x = as.numeric(coords[, 1]),
    y = as.numeric(coords[, 2]),
    cellType = as.character(cd[[ct_col]]),
    imageID = as.character(sample_id),
    stringsAsFactors = FALSE
  )
}
cell_df <- dplyr::bind_rows(cell_rows) |>
  dplyr::filter(
    is.finite(x), is.finite(y),
    !is.na(cellType), nzchar(trimws(cellType))
  )

drop_cell_types <- as.character(analysis_params$spicy$drop_cell_types)
drop_cell_types <- drop_cell_types[nzchar(drop_cell_types)]
cell_df <- cell_df |>
  dplyr::filter(!(.data$cellType %in% drop_cell_types))

zero_thr <- as.numeric(analysis_params$spicy$drop_zero_fraction_gt)
prop_tab <- table(cell_df$imageID, cell_df$cellType)
drop_sparse <- names(which(colMeans(prop_tab == 0) > zero_thr))
cell_df <- cell_df |>
  dplyr::filter(!(.data$cellType %in% drop_sparse))

cell_exp <- spicyR::SegmentedCells(
  cellData = cell_df,
  spatialCoords = c("x", "y"),
  imageIDString = "imageID"
)

pheno <- sample_meta |>
  dplyr::transmute(
    imageID = as.character(sample),
    patient = as.factor(patient),
    group = factor(as.character(group), levels = group_order)
  ) |>
  as.data.frame()
spicyR::imagePheno(cell_exp) <- pheno

spicy_subject_col <- analysis_params$spicy$subject_col
spicy_subject_col <- as.character(spicy_subject_col)
if (length(spicy_subject_col) != 1L || !nzchar(spicy_subject_col)) {
  stop("analysis_params$spicy$subject_col must be one non-empty string.")
}

rs_step <- as.numeric(analysis_params$spicy$rs_step_um)
Rs <- sort(unique(c(0, seq(0, max(radii_um), by = rs_step))))

workers <- as.integer(analysis_params$spicy$workers)
BPPARAM <- if (workers > 1L) {
  BiocParallel::SnowParam(workers = workers, type = "SOCK", progressbar = TRUE)
} else {
  BiocParallel::SerialParam(progressbar = TRUE)
}

spicy_fit <- spicyR::spicy(
  cell_exp,
  condition = "group",
  subject = spicy_subject_col,
  BPPARAM = BPPARAM,
  Rs = Rs,
  includeZeroCells = isTRUE(analysis_params$spicy$include_zero_cells)
)

spicy_prefix <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark"))
saveRDS(spicy_fit, paste0(spicy_prefix, "_spicyR_object.rds"))

top_n <- as.integer(analysis_params$spicy$top_n)
spicy_raw <- spicyR::topPairs(spicy_fit, n = top_n)
readr::write_csv(as.data.frame(spicy_raw), paste0(spicy_prefix, "_spicyR_topPairs_raw.csv"))

spicy_tbl <- data.frame(
  ct1 = as.character(spicy_raw$from),
  ct2 = as.character(spicy_raw$to),
  effect_spicyr = as.numeric(spicy_raw$coefficient),
  p_spicyr = as.numeric(spicy_raw$p.value),
  stringsAsFactors = FALSE
)
spicy_tbl <- spicy_tbl |>
  dplyr::filter(!is.na(ct1), !is.na(ct2), nzchar(ct1), nzchar(ct2))
spicy_tbl$fdr_spicyr <- p.adjust(spicy_tbl$p_spicyr, method = "fdr")
spicy_tbl$pair_key_dir <- pair_key_dir(spicy_tbl$ct1, spicy_tbl$ct2)
spicy_tbl$pair_key_undir <- pair_key_undir(spicy_tbl$ct1, spicy_tbl$ct2)

pan_dir_merge <- dplyr::left_join(
  pan_tbl,
  ttest_tbl |>
    dplyr::select(ct1, ct2, radius_um, stat, effect_ttest, p_ttest, fdr_ttest, n_control, n_case, status),
  by = c("ct1", "ct2", "radius_um", "stat")
)

pan_pair <- build_pair_level(pan_tbl, p_col = "p_panoramic", fdr_col = "fdr_panoramic", effect_col = "beta_panoramic")
ttest_pair <- build_pair_level(ttest_tbl, p_col = "p_ttest", fdr_col = "fdr_ttest", effect_col = "effect_ttest")
spicy_pair <- build_pair_level(spicy_tbl, p_col = "p_spicyr", fdr_col = "fdr_spicyr", effect_col = "effect_spicyr")

pair_merge <- dplyr::full_join(
  dplyr::full_join(
    pan_pair,
    ttest_pair |>
      dplyr::select(pair_key_undir, p_ttest, fdr_ttest, effect_ttest),
    by = "pair_key_undir"
  ),
  spicy_pair |>
    dplyr::select(pair_key_undir, p_spicyr, fdr_spicyr, effect_spicyr),
  by = "pair_key_undir"
)

alpha <- 0.05
set_pan <- unique(pair_merge$pair_key_undir[is.finite(pair_merge$fdr_panoramic) & pair_merge$fdr_panoramic <= alpha])
set_ttest <- unique(pair_merge$pair_key_undir[is.finite(pair_merge$fdr_ttest) & pair_merge$fdr_ttest <= alpha])
set_spicy <- unique(pair_merge$pair_key_undir[is.finite(pair_merge$fdr_spicyr) & pair_merge$fdr_spicyr <= alpha])

corr_pan_t <- cor(-log10(pair_merge$p_panoramic), -log10(pair_merge$p_ttest), use = "pairwise.complete.obs", method = "spearman")
corr_pan_s <- cor(-log10(pair_merge$p_panoramic), -log10(pair_merge$p_spicyr), use = "pairwise.complete.obs", method = "spearman")
corr_eff_pan_t <- cor(pan_dir_merge$beta_panoramic, pan_dir_merge$effect_ttest, use = "pairwise.complete.obs", method = "pearson")

summary_tbl <- data.frame(
  metric = c(
    "cohort",
    "run_prefix",
    "group_control",
    "group_case",
    "n_features_panoramic_directional",
    "n_features_ttest_directional",
    "n_pairs_spicyr_directional",
    "n_pairs_compared_pair_level",
    "n_sig_panoramic_fdr_0_05",
    "n_sig_ttest_fdr_0_05",
    "n_sig_spicyr_fdr_0_05",
    "n_overlap_sig_pan_ttest",
    "n_overlap_sig_pan_spicyr",
    "n_overlap_sig_ttest_spicyr",
    "spearman_logp_pan_vs_ttest",
    "spearman_logp_pan_vs_spicyr",
    "pearson_effect_pan_vs_ttest"
  ),
  value = c(
    "crc_tma",
    run_prefix,
    group_order[[1]],
    group_order[[2]],
    nrow(pan_tbl),
    nrow(ttest_tbl),
    nrow(spicy_tbl),
    nrow(pair_merge),
    sum(is.finite(pair_merge$fdr_panoramic) & pair_merge$fdr_panoramic <= alpha),
    sum(is.finite(pair_merge$fdr_ttest) & pair_merge$fdr_ttest <= alpha),
    sum(is.finite(pair_merge$fdr_spicyr) & pair_merge$fdr_spicyr <= alpha),
    length(intersect(set_pan, set_ttest)),
    length(intersect(set_pan, set_spicy)),
    length(intersect(set_ttest, set_spicy)),
    corr_pan_t,
    corr_pan_s,
    corr_eff_pan_t
  ),
  stringsAsFactors = FALSE
)

pan_csv <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_panoramic.csv"))
ttest_csv <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_ttest.csv"))
spicy_csv <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_spicyr.csv"))
merge_dir_csv <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_merge_directional.csv"))
merge_pair_csv <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_merge_pair_level.csv"))
summary_csv <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_summary.csv"))

readr::write_csv(pan_tbl, pan_csv)
readr::write_csv(ttest_tbl, ttest_csv)
readr::write_csv(spicy_tbl, spicy_csv)
readr::write_csv(pan_dir_merge, merge_dir_csv)
readr::write_csv(pair_merge, merge_pair_csv)
readr::write_csv(summary_tbl, summary_csv)

if (isTRUE(write_latex_results)) {
  alpha_slug <- gsub("\\.", "p", format(latex_alpha, scientific = FALSE, trim = TRUE))

  if (is.data.frame(meta_tbl) && nrow(meta_tbl) > 0L) {
    pan_latex_full <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_panoramic_results_table.tex"))
    pan_latex_selected <- if (isTRUE(latex_sig_only)) {
      file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_panoramic_results_table_sig", alpha_slug, ".tex"))
    } else {
      pan_latex_full
    }

    shared_helpers$write_panoramic_latex_results(
      meta_tbl = meta_tbl,
      out_path = pan_latex_full,
      sig_only = FALSE,
      alpha = latex_alpha,
      caption = "CRC TMA PANORAMIC results (benchmark context).",
      label = "tab:s_crc_tma_benchmark_panoramic",
      print_console = FALSE
    )
    pan_out <- shared_helpers$write_panoramic_latex_results(
      meta_tbl = meta_tbl,
      out_path = pan_latex_selected,
      sig_only = latex_sig_only,
      alpha = latex_alpha,
      caption = if (isTRUE(latex_sig_only)) {
        paste0("CRC TMA PANORAMIC results (FDR <= ", format(latex_alpha, trim = TRUE), ").")
      } else {
        "CRC TMA PANORAMIC results (benchmark context)."
      },
      label = "tab:s_crc_tma_benchmark_panoramic",
      print_console = latex_print_console
    )
    message("LaTeX PANORAMIC results table: ", pan_out$path, " (rows=", pan_out$n_rows, ")")
  }

  ttest_meta_tbl <- build_method_latex_meta(
    tbl = ttest_tbl,
    group_order = group_order,
    effect_col = "effect_ttest",
    p_col = "p_ttest",
    fdr_col = "fdr_ttest",
    group1_col = "mean_control",
    group2_col = "mean_case"
  )
  ttest_latex_full <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_ttest_results_table.tex"))
  ttest_latex_selected <- if (isTRUE(latex_sig_only)) {
    file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_ttest_results_table_sig", alpha_slug, ".tex"))
  } else {
    ttest_latex_full
  }
  shared_helpers$write_panoramic_latex_results(
    meta_tbl = ttest_meta_tbl,
    out_path = ttest_latex_full,
    sig_only = FALSE,
    alpha = latex_alpha,
    caption = "CRC TMA t-test results (benchmark context).",
    label = "tab:s_crc_tma_benchmark_ttest",
    print_console = FALSE
  )
  ttest_out <- shared_helpers$write_panoramic_latex_results(
    meta_tbl = ttest_meta_tbl,
    out_path = ttest_latex_selected,
    sig_only = latex_sig_only,
    alpha = latex_alpha,
    caption = if (isTRUE(latex_sig_only)) {
      paste0("CRC TMA t-test results (FDR <= ", format(latex_alpha, trim = TRUE), ").")
    } else {
      "CRC TMA t-test results (benchmark context)."
    },
    label = "tab:s_crc_tma_benchmark_ttest",
    print_console = latex_print_console
  )
  message("LaTeX t-test results table: ", ttest_out$path, " (rows=", ttest_out$n_rows, ")")

  spicyr_meta_tbl <- build_method_latex_meta(
    tbl = spicy_tbl,
    group_order = group_order,
    effect_col = "effect_spicyr",
    p_col = "p_spicyr",
    fdr_col = "fdr_spicyr"
  )
  spicyr_latex_full <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_spicyr_results_table.tex"))
  spicyr_latex_selected <- if (isTRUE(latex_sig_only)) {
    file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_spicyr_results_table_sig", alpha_slug, ".tex"))
  } else {
    spicyr_latex_full
  }
  shared_helpers$write_panoramic_latex_results(
    meta_tbl = spicyr_meta_tbl,
    out_path = spicyr_latex_full,
    sig_only = FALSE,
    alpha = latex_alpha,
    caption = "CRC TMA spicyR results (benchmark context).",
    label = "tab:s_crc_tma_benchmark_spicyr",
    print_console = FALSE
  )
  spicyr_out <- shared_helpers$write_panoramic_latex_results(
    meta_tbl = spicyr_meta_tbl,
    out_path = spicyr_latex_selected,
    sig_only = latex_sig_only,
    alpha = latex_alpha,
    caption = if (isTRUE(latex_sig_only)) {
      paste0("CRC TMA spicyR results (FDR <= ", format(latex_alpha, trim = TRUE), ").")
    } else {
      "CRC TMA spicyR results (benchmark context)."
    },
    label = "tab:s_crc_tma_benchmark_spicyr",
    print_console = latex_print_console
  )
  message("LaTeX spicyR results table: ", spicyr_out$path, " (rows=", spicyr_out$n_rows, ")")

  pan_ttest_full <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_pan_vs_ttest_results_table.tex"))
  pan_ttest_selected <- if (isTRUE(latex_sig_only)) {
    file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_pan_vs_ttest_results_table_sig", alpha_slug, ".tex"))
  } else {
    pan_ttest_full
  }
  pan_ttest_out <- write_pan_vs_method_latex_results(
    pair_merge = pair_merge,
    method = "ttest",
    out_path = pan_ttest_selected,
    sig_only = latex_sig_only,
    alpha = latex_alpha,
    caption = if (isTRUE(latex_sig_only)) {
      paste0("CRC TMA PANORAMIC vs t-test benchmark results (FDR <= ", format(latex_alpha, trim = TRUE), " in PANORAMIC or t-test).")
    } else {
      "CRC TMA PANORAMIC vs t-test benchmark results."
    },
    label = "tab:s_crc_tma_pan_vs_ttest",
    print_console = latex_print_console
  )
  if (!identical(pan_ttest_selected, pan_ttest_full)) {
    write_pan_vs_method_latex_results(
      pair_merge = pair_merge,
      method = "ttest",
      out_path = pan_ttest_full,
      sig_only = FALSE,
      alpha = latex_alpha,
      caption = "CRC TMA PANORAMIC vs t-test benchmark results.",
      label = "tab:s_crc_tma_pan_vs_ttest",
      print_console = FALSE
    )
  }
  message("LaTeX PANORAMIC vs t-test table: ", pan_ttest_out$path, " (rows=", pan_ttest_out$n_rows, ")")

  pan_spicy_full <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_pan_vs_spicyr_results_table.tex"))
  pan_spicy_selected <- if (isTRUE(latex_sig_only)) {
    file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_pan_vs_spicyr_results_table_sig", alpha_slug, ".tex"))
  } else {
    pan_spicy_full
  }
  pan_spicy_out <- write_pan_vs_method_latex_results(
    pair_merge = pair_merge,
    method = "spicyr",
    out_path = pan_spicy_selected,
    sig_only = latex_sig_only,
    alpha = latex_alpha,
    caption = if (isTRUE(latex_sig_only)) {
      paste0("CRC TMA PANORAMIC vs spicyR benchmark results (FDR <= ", format(latex_alpha, trim = TRUE), " in PANORAMIC or spicyR).")
    } else {
      "CRC TMA PANORAMIC vs spicyR benchmark results."
    },
    label = "tab:s_crc_tma_pan_vs_spicyr",
    print_console = latex_print_console
  )
  if (!identical(pan_spicy_selected, pan_spicy_full)) {
    write_pan_vs_method_latex_results(
      pair_merge = pair_merge,
      method = "spicyr",
      out_path = pan_spicy_full,
      sig_only = FALSE,
      alpha = latex_alpha,
      caption = "CRC TMA PANORAMIC vs spicyR benchmark results.",
      label = "tab:s_crc_tma_pan_vs_spicyr",
      print_console = FALSE
    )
  }
  message("LaTeX PANORAMIC vs spicyR table: ", pan_spicy_out$path, " (rows=", pan_spicy_out$n_rows, ")")

  summary_tex <- file.path(benchmark_tab_dir, paste0(run_prefix, "_benchmark_summary_table.tex"))
  summary_out <- write_benchmark_summary_latex(
    summary_tbl = summary_tbl,
    out_path = summary_tex,
    caption = "CRC TMA benchmark summary metrics.",
    label = "tab:s_crc_tma_benchmark_summary",
    print_console = latex_print_console
  )
  message("LaTeX benchmark summary table: ", summary_out$path, " (rows=", summary_out$n_rows, ")")
}

write_scatter(
  df = pan_dir_merge,
  x_col = "p_panoramic",
  y_col = "p_ttest",
  x_lab = "-log10(p) PANORAMIC",
  y_lab = "-log10(p) t-test",
  title = "CRC TMA: PANORAMIC vs t-test",
  out_stub = file.path(benchmark_fig_dir, paste0(run_prefix, "_benchmark_pan_vs_ttest_logp"))
)
write_scatter(
  df = pair_merge,
  x_col = "p_panoramic",
  y_col = "p_spicyr",
  x_lab = "-log10(p) PANORAMIC",
  y_lab = "-log10(p) spicyR",
  title = "CRC TMA: PANORAMIC vs spicyR",
  out_stub = file.path(benchmark_fig_dir, paste0(run_prefix, "_benchmark_pan_vs_spicyr_logp"))
)

write_pvalue_hist(
  df = ttest_tbl,
  p_col = "p_ttest",
  title = "CRC TMA t-test p-value histogram",
  out_stub = file.path(benchmark_fig_dir, paste0(run_prefix, "_benchmark_ttest_pvalue_histogram")),
  bins = 20L
)
write_pvalue_hist(
  df = spicy_tbl,
  p_col = "p_spicyr",
  title = "CRC TMA spicyR p-value histogram",
  out_stub = file.path(benchmark_fig_dir, paste0(run_prefix, "_benchmark_spicyr_pvalue_histogram")),
  bins = 20L
)

message("CRC benchmark complete.")
message("  PANORAMIC table: ", pan_csv)
message("  t-test table: ", ttest_csv)
message("  spicyR table: ", spicy_csv)
message("  merged pair table: ", merge_pair_csv)
message("  summary: ", summary_csv)
