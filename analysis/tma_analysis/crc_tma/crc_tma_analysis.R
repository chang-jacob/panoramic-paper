#!/usr/bin/env Rscript

# CRC TMA Script: Primary PANORAMIC Analysis
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Run prepare -> spatial statistics/bootstrap -> multilevel meta-analysis.
# - Export spatial/meta/contrast tables and runtime metrics.
# - Generate volcano, representative sample, network, forest, and LaTeX outputs.

suppressPackageStartupMessages({
  library(config)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(BiocParallel)
  library(panoramic)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ---- Config ----
config_file <- file.path(getwd(), "config", "default.yml")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}
paths <- config::get("paths", file = config_file)
source(file.path(getwd(), "analysis", "tma_analysis", "tma_shared_helpers.R"))

# ---- Analysis Parameters (edit in script) ----
analysis_params <- list(
  # Dataset selection for this analysis script.
  spe_list_key = "one_per_patient",
  spe_list_path = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_one_per_patient.rds"),

  # Column names in colData.
  columns = list(
    cell_type = "cell_type",
    patient = "patient",
    group = "group_label",
    region = "spot"
  ),

  # PANORAMIC spatial stats.
  panoramic = list(
    radii_um = c(25),
    stat = "local_comp_enrichment", # local_comp_enrichment | Lcross | Kcross | Lest | Kest
    nsim = 100L,
    correction = "translate",
    min_cells = 5L,
    concavity = 50,
    window = "convex", # concave | convex | rect
    seed = 123L,
    boot = "block", # approx | block
    tile_size = 62.5,
    workers = 10L,
    progressbar = TRUE
  ),

  # Meta-analysis.
  meta = list(
    method_mv = "REML",
    group_tau2 = "none", # none | separate
    vi_floor = "group_median", # group_median | median | none
    tau_structure = "patient" # patient | patient_sample
    
  ),

  # Plot/table outputs.
  plotting = list(
    volcano_point_size = 0.8,
    volcano_x_scale = "log2fc", # beta_diff | log2fc
    write_plot_debug = FALSE,
    network_leiden_res = 1.2,
    network_z_sign = "both", # both | positive | negative
    network_sig_operator = "gt", # lt | gt
    network_include_nonsig = TRUE,
    forest_top_n = 10L,
    forest_alpha = 0.05,
    forest_sig_col = "fdr_diff", # fdr_diff | p_diff
    forest_variant = "both", # both | full | minimal
    forest_show_est_se = TRUE,
    forest_show_ci = TRUE,
    show_plots = FALSE,
    show_plots_max = 6L,
    save_p_hist = FALSE
  ),

  reporting = list(
    write_latex_results = TRUE,
    latex_sig_only = TRUE,
    latex_alpha = 0.05,
    latex_print_console = TRUE
  ),

  output = list(
    analysis_subdir = NULL
  )
)

merge_nested_lists <- function(base, override) {
  if (is.null(override)) {
    return(base)
  }

  override_names <- names(override)
  if (is.null(override_names)) {
    stop("analysis_params_override must be a named list.")
  }

  for (nm in override_names) {
    if (nm %in% names(base) && is.list(base[[nm]]) && is.list(override[[nm]])) {
      base[[nm]] <- merge_nested_lists(base[[nm]], override[[nm]])
    } else {
      base[[nm]] <- override[[nm]]
    }
  }

  base
}

if (exists("analysis_params_override", inherits = TRUE)) {
  analysis_params <- merge_nested_lists(
    analysis_params,
    get("analysis_params_override", inherits = TRUE)
  )
}

# ---- Run controls ----
spe_list_key <- analysis_params$spe_list_key
spe_list_path <- analysis_params$spe_list_path
list_dir <- spe_list_key

# ---- Column names in colData(spe) ----
cell_type_col <- analysis_params$columns$cell_type
patient_col <- analysis_params$columns$patient
group_col <- analysis_params$columns$group
region_col <- analysis_params$columns$region

# ---- PANORAMIC parameters ----
radii_um <- analysis_params$panoramic$radii_um
stat <- analysis_params$panoramic$stat
nsim <- analysis_params$panoramic$nsim
correction <- analysis_params$panoramic$correction
min_cells <- analysis_params$panoramic$min_cells
concavity <- analysis_params$panoramic$concavity
window <- analysis_params$panoramic$window
seed <- analysis_params$panoramic$seed
boot <- analysis_params$panoramic$boot
tile_size <- analysis_params$panoramic$tile_size
workers <- analysis_params$panoramic$workers
progressbar <- analysis_params$panoramic$progressbar
BPPARAM <- BiocParallel::SnowParam(workers = workers, progressbar = progressbar)

# ---- Meta-analysis parameters ----
method_mv <- analysis_params$meta$method_mv
group_tau2 <- analysis_params$meta$group_tau2
vi_floor <- analysis_params$meta$vi_floor
tau_structure <- analysis_params$meta$tau_structure

# ---- Plot parameters ----
volcano_point_size <- analysis_params$plotting$volcano_point_size
volcano_x_scale <- analysis_params$plotting$volcano_x_scale
write_plot_debug <- analysis_params$plotting$write_plot_debug
network_leiden_res <- analysis_params$plotting$network_leiden_res
network_z_sign <- analysis_params$plotting$network_z_sign
network_sig_operator <- analysis_params$plotting$network_sig_operator
network_include_nonsig <- analysis_params$plotting$network_include_nonsig
forest_top_n <- analysis_params$plotting$forest_top_n
forest_alpha <- analysis_params$plotting$forest_alpha
forest_sig_col <- analysis_params$plotting$forest_sig_col
forest_variant <- analysis_params$plotting$forest_variant
forest_show_est_se <- analysis_params$plotting$forest_show_est_se
forest_show_ci <- analysis_params$plotting$forest_show_ci
show_plots <- analysis_params$plotting$show_plots
show_plots_max <- analysis_params$plotting$show_plots_max
save_p_hist <- analysis_params$plotting$save_p_hist

# ---- Reporting parameters ----
write_latex_results <- analysis_params$reporting$write_latex_results
latex_sig_only <- analysis_params$reporting$latex_sig_only
latex_alpha <- analysis_params$reporting$latex_alpha
latex_print_console <- analysis_params$reporting$latex_print_console

# ---- Output parameters ----
analysis_subdir <- analysis_params$output$analysis_subdir
if (is.null(analysis_subdir)) {
  use_analysis_subdir <- FALSE
} else {
  analysis_subdir <- as.character(analysis_subdir[[1]])
  use_analysis_subdir <- nzchar(analysis_subdir)
}

# ---- Outputs ----
interim_dir <- if (use_analysis_subdir) {
  file.path(paths$data, "interim", "crc_tma", list_dir, analysis_subdir)
} else {
  file.path(paths$data, "interim")
}
tab_dir <- if (use_analysis_subdir) {
  file.path(paths$output, "tables", "crc_tma", list_dir, analysis_subdir)
} else {
  file.path(paths$output, "tables", "crc_tma", list_dir)
}
fig_dir <- if (use_analysis_subdir) {
  file.path(paths$output, "figures", "crc_tma", list_dir, analysis_subdir)
} else {
  file.path(paths$output, "figures", "crc_tma", list_dir)
}
dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stat_tag <- gsub("[^a-z0-9]+", "_", tolower(stat))
radii_tag <- paste(radii_um, collapse = "_")
tile_tag <- format(signif(tile_size, 6), scientific = FALSE, trim = TRUE)
tile_slug <- gsub("[^a-z0-9]+", "", tolower(gsub("\\.", "p", tile_tag)))
run_prefix <- paste0(
  "crc_tma_", spe_list_key,
  "_", stat_tag,
  "_r", radii_tag,
  "_nsim", nsim,
  "_", boot,
  "_tile", tile_slug,
  "_mv"
)
prep_path <- file.path(interim_dir, paste0(run_prefix, "_prep.rds"))
stats_path <- file.path(interim_dir, paste0(run_prefix, "_se_stats.rds"))
meta_path <- file.path(interim_dir, paste0(run_prefix, "_se_meta.rds"))
runtime_csv <- file.path(tab_dir, paste0(run_prefix, "_runtime_metrics.csv"))

# ---- Helpers ----
runtime_tracker <- create_runtime_tracker(run_prefix)
append_runtime_metric <- runtime_tracker$append_runtime_metric
write_runtime_metrics <- runtime_tracker$write_runtime_metrics

build_volcano_df <- function(contrast_df, control_group, case_group,
                             sig_col = c("fdr_diff", "p_diff"),
                             alpha = 0.05) {
  sig_col <- match.arg(sig_col)
  df <- contrast_df
  # beta_diff is the modeled group contrast; expose it as log2fc for plotting.
  df$log2fc <- df$beta_diff
  df$feature_label <- as.character(df$coloc_direction)
  p_for_axis <- pmax(df$p_diff, .Machine$double.xmin)
  df$neglog10_p <- -log10(p_for_axis)
  sig_val <- df[[sig_col]]
  df$is_sig <- is.finite(sig_val) & (sig_val <= alpha)
  df$direction_group <- ifelse(df$log2fc >= 0, case_group, control_group)
  df$color_class <- ifelse(df$is_sig, df$direction_group, "Not significant")
  df
}

sanitize_slug <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "NA")
}

force_volcano_geom_sizes <- function(p, point_size = NULL, label_text_pt = NULL) {
  if (!inherits(p, "ggplot")) return(p)

  point_size <- as.numeric(point_size)
  if (length(point_size) != 1L || is.na(point_size) || !is.finite(point_size) || point_size <= 0) point_size <- NA_real_

  label_text_pt <- as.numeric(label_text_pt)
  if (length(label_text_pt) != 1L || is.na(label_text_pt) || !is.finite(label_text_pt) || label_text_pt <= 0) {
    label_size <- NA_real_
  } else {
    label_size <- label_text_pt / ggplot2::.pt
  }

  for (i in seq_along(p$layers)) {
    geom_i <- p$layers[[i]]$geom
    geom_classes <- class(geom_i)

    if (is.finite(point_size) && any(geom_classes %in% c("GeomPoint"))) {
      p$layers[[i]]$aes_params$size <- point_size
    }
    if (is.finite(label_size) && any(geom_classes %in% c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel"))) {
      p$layers[[i]]$aes_params$size <- label_size
    }
  }

  p
}

rank_forest_hits <- function(contrast_tbl, sig_col = "fdr_diff", alpha = 0.05) {
  if (!is.data.frame(contrast_tbl)) {
    stop("contrast_tbl must be a data.frame.")
  }
  if (nrow(contrast_tbl) == 0L) {
    stop("contrast_tbl has 0 rows.")
  }
  required_cols <- c("ct1", "ct2", "radius_um", "p_diff")
  if (!all(required_cols %in% colnames(contrast_tbl))) {
    missing_cols <- setdiff(required_cols, colnames(contrast_tbl))
    stop("contrast_tbl missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!sig_col %in% colnames(contrast_tbl)) {
    stop("sig_col='", sig_col, "' missing from contrast_tbl.")
  }

  hits <- if (identical(sig_col, "fdr_diff")) {
    contrast_tbl |>
      dplyr::filter(is.finite(fdr_diff), fdr_diff <= alpha) |>
      dplyr::arrange(fdr_diff, p_diff)
  } else {
    contrast_tbl |>
      dplyr::filter(is.finite(p_diff), p_diff <= alpha) |>
      dplyr::arrange(p_diff)
  }

  if (nrow(hits) == 0L) {
    stop("No rows passed ", sig_col, " <= ", alpha, " in rank_forest_hits().")
  }

  hits
}

select_forest_hits <- function(contrast_tbl, sig_col = "fdr_diff", alpha = 0.05, top_n = 10L,
                               beta_sign = c("negative", "nonnegative")) {
  beta_sign <- match.arg(beta_sign)
  hits <- rank_forest_hits(contrast_tbl, sig_col = sig_col, alpha = alpha)

  if (!"beta_diff" %in% colnames(hits)) {
    stop("beta_diff column is required in contrast_tbl.")
  }
  hits <- if (identical(beta_sign, "negative")) {
    dplyr::filter(hits, is.finite(beta_diff), beta_diff < 0)
  } else {
    dplyr::filter(hits, is.finite(beta_diff), beta_diff >= 0)
  }
  if (nrow(hits) == 0L) {
    stop("No rows passed beta_sign='", beta_sign, "' after significance filtering.")
  }

  hits |>
    dplyr::distinct(ct1, ct2, radius_um, .keep_all = TRUE) |>
    dplyr::slice_head(n = top_n)
}

save_network_outputs <- function(se_meta, fig_dir, run_prefix) {
  net_resolution <- network_leiden_res
  net_z_sign_local <- network_z_sign
  z_sign_values <- if (identical(net_z_sign_local, "both")) c("positive", "negative") else net_z_sign_local
  net_sig_operator <- network_sig_operator
  net_include_nonsig <- network_include_nonsig

  png_files <- character(0)
  pdf_files <- character(0)
  rds_files <- character(0)
  summary_files <- character(0)

  for (z_sign_value in z_sign_values) {
    sign_suffix <- paste0("_", z_sign_value)
    nonsig_suffix <- if (isTRUE(net_include_nonsig)) "_with_nonsig" else ""
    out_prefix <- file.path(fig_dir, paste0(run_prefix, "_network", sign_suffix, nonsig_suffix))

    net_out <- panoramic::plot_spatial_network(
      se_diff = se_meta,
      leiden_resolution = net_resolution,
      z_sign = z_sign_value,
      include_nonsig = net_include_nonsig,
      sig_operator = net_sig_operator,
      return_net = TRUE
    )
    net <- net_out$net
    p_net <- net_out$plot + ggplot2::labs(title = NULL, subtitle = NULL)

    out_png <- paste0(out_prefix, ".png")
    out_pdf <- paste0(out_prefix, ".pdf")
    out_rds <- paste0(out_prefix, ".rds")
    out_summary <- paste0(out_prefix, "_summary.csv")

    ggsave(out_png, p_net, width = 8.5, height = 7, dpi = 300)
    ggsave(out_pdf, p_net, width = 8.5, height = 7)
    saveRDS(net, out_rds)
    readr::write_csv(
      data.frame(
        z_sign = z_sign_value,
        directed = FALSE,
        sig_operator = net_sig_operator,
        include_nonsig = net_include_nonsig,
        nonsig_max_fdr = 1.0,
        n_clusters = net$n_clusters,
        modularity = net$modularity,
        n_nodes = igraph::vcount(net$graph),
        n_edges = igraph::ecount(net$graph),
        stringsAsFactors = FALSE
      ),
      out_summary
    )

    png_files <- c(png_files, out_png)
    pdf_files <- c(pdf_files, out_pdf)
    rds_files <- c(rds_files, out_rds)
    summary_files <- c(summary_files, out_summary)
  }

  message("Network PNG files: ", paste(png_files, collapse = ", "))
  message("Network PDF files: ", paste(pdf_files, collapse = ", "))
  message("Network summary CSV files: ", paste(summary_files, collapse = ", "))

  list(
    png = png_files,
    pdf = pdf_files,
    rds = rds_files,
    summary = summary_files,
    z_sign = z_sign_values
  )
}

save_forest_outputs <- function(se_meta, contrast_tbl, fig_dir, run_prefix, group_col) {
  top_n <- forest_top_n
  alpha <- forest_alpha
  sig_col <- forest_sig_col

  variant_specs <- if (identical(forest_variant, "both")) {
    list(
      list(variant = "full", show_est_se = TRUE, show_ci = TRUE),
      list(variant = "minimal", show_est_se = FALSE, show_ci = FALSE)
    )
  } else {
    if (is.null(forest_show_est_se) || is.null(forest_show_ci)) {
      stop("For forest_variant != 'both', forest_show_est_se and forest_show_ci must be explicitly set.")
    }
    list(
      list(
        variant = forest_variant,
        show_est_se = forest_show_est_se,
        show_ci = forest_show_ci
      )
    )
  }

  contrast_meta <- S4Vectors::metadata(se_meta)$panoramic$contrast
  if (is.null(contrast_meta$control) || is.null(contrast_meta$case)) {
    stop("contrast metadata must include explicit control/case labels.")
  }
  control_group <- as.character(contrast_meta$control)
  case_group <- as.character(contrast_meta$case)

  maybe_select_forest_hits <- function(beta_sign, direction_group) {
    hits <- tryCatch(
      select_forest_hits(
        contrast_tbl,
        sig_col = sig_col,
        alpha = alpha,
        top_n = top_n,
        beta_sign = beta_sign
      ),
      error = function(e) NULL
    )
    if (is.null(hits) || nrow(hits) == 0L) {
      return(NULL)
    }
    hits$direction_group <- direction_group
    hits$direction_rank <- seq_len(nrow(hits))
    hits
  }

  hits_ctrl <- maybe_select_forest_hits(beta_sign = "negative", direction_group = control_group)
  hits_case <- maybe_select_forest_hits(beta_sign = "nonnegative", direction_group = case_group)
  hit_tables <- Filter(function(x) !is.null(x) && nrow(x) > 0L, list(hits_ctrl, hits_case))
  if (length(hit_tables) == 0L) {
    message("No significant forest hits at alpha=", alpha, " using ", sig_col, "; skipping forest outputs.")
    return(list(index = data.frame(), index_csv = NA_character_, dir = NA_character_, plots = list()))
  }
  hits <- dplyr::bind_rows(hit_tables)

  forest_dir <- file.path(fig_dir, paste0(run_prefix, "_forest"))
  dir.create(forest_dir, recursive = TRUE, showWarnings = FALSE)

  index_rows <- list()
  plots <- list()
  idx_out <- 0L
  shown_count <- 0L

  for (i in seq_len(nrow(hits))) {
    hit <- hits[i, , drop = FALSE]
    ct1 <- as.character(hit$ct1)
    ct2 <- as.character(hit$ct2)
    direction_group <- as.character(hit$direction_group[1])
    direction_rank <- as.integer(hit$direction_rank[1])
    direction_slug <- paste0("higher_in_", sanitize_slug(direction_group))
    radius_um <- as.numeric(hit$radius_um)
    radius_slug <- gsub("\\.", "p", format(radius_um, trim = TRUE, scientific = FALSE))
    file_stub <- file.path(
      forest_dir,
      paste0(
        sprintf("%03d_", i),
        direction_slug, "_",
        sanitize_slug(ct1), "_to_", sanitize_slug(ct2),
        "_r", radius_slug, "um"
      )
    )
    for (spec in variant_specs) {
      variant <- spec$variant
      show_est_se <- isTRUE(spec$show_est_se)
      show_ci <- isTRUE(spec$show_ci)
      out_png <- paste0(file_stub, "_", variant, ".png")
      out_pdf <- paste0(file_stub, "_", variant, ".pdf")
      idx_out <- idx_out + 1L

      p <- NULL
      forest_args <- list(
        se_meta = se_meta,
        ct1 = ct1,
        ct2 = ct2,
        radius_um = radius_um,
        group_col = group_col
      )
      if (!show_est_se) forest_args$show_est_se <- FALSE
      if (!show_ci) forest_args$show_ci <- FALSE
      p <- do.call(panoramic::plot_forest, forest_args) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 5)
        )
      plots[[idx_out]] <- p
      if (isTRUE(show_plots) && shown_count < show_plots_max) {
        print(p)
        shown_count <- shown_count + 1L
      }
      ggsave(out_png, p, width = 12, height = 4, dpi = 300)
      ggsave(out_pdf, p, width = 12, height = 4)

      sig_value <- hit[[sig_col]][1]
      index_rows[[idx_out]] <- data.frame(
        rank = idx_out,
        hit_rank = i,
        variant = variant,
        show_est_se = show_est_se,
        show_ci = show_ci,
        direction_group = direction_group,
        direction_rank = direction_rank,
        ct1 = ct1,
        ct2 = ct2,
        radius_um = radius_um,
        p_diff = hit$p_diff[1],
        fdr_diff = hit$fdr_diff[1],
        sig_col = sig_col,
        sig_value = sig_value,
        status = "ok",
        error = NA_character_,
        png = out_png,
        pdf = out_pdf,
        stringsAsFactors = FALSE
      )
    }
  }

  index_df <- dplyr::bind_rows(index_rows)
  index_csv <- file.path(forest_dir, "forest_plot_index.csv")
  readr::write_csv(index_df, index_csv)
  message("Forest plots attempted: ", nrow(index_df))
  message("Forest hits (higher in ", control_group, "): ", if (is.null(hits_ctrl)) 0L else nrow(hits_ctrl))
  message("Forest hits (higher in ", case_group, "): ", if (is.null(hits_case)) 0L else nrow(hits_case))
  if (isTRUE(show_plots) && nrow(index_df) > shown_count) {
    message("Displayed first ", shown_count, " forest plots in-session (adjust analysis_params$plotting$show_plots_max to change).")
  }
  message("Forest index CSV: ", index_csv)
  list(index = index_df, index_csv = index_csv, dir = forest_dir, plots = plots)
}

# ---- Read preprocessed data ----
message("Reading preprocessed data: ", spe_list_path)
spe_list <- readRDS(spe_list_path)

sample_meta <- lapply(spe_list, function(spe) {
  data.frame(
    patient = get_unique(spe, patient_col),
    group_label = get_unique(spe, group_col),
    region = get_unique(spe, region_col),
    stringsAsFactors = FALSE
  )
})
sample_meta <- do.call(rbind, sample_meta)
sample_meta$sample <- sample_meta$region
names(spe_list) <- sample_meta$sample

design <- data.frame(
  sample = sample_meta$sample,
  group = sample_meta$group_label,
  stringsAsFactors = FALSE
)

# ---- Step 1: Prepare ----
t_prepare_start <- Sys.time()
prep <- panoramic_prepare(
  spe_list = spe_list,
  design = design,
  cell_type = cell_type_col,
  min_cells = min_cells,
  concavity = concavity,
  window = window,
  BPPARAM = BPPARAM
)
saveRDS(prep, prep_path)
t_prepare_end <- Sys.time()
append_runtime_metric(
  stage = "prepare",
  status = "computed",
  started_at = t_prepare_start,
  ended_at = t_prepare_end,
  cache_path = prep_path,
  params = list(
    min_cells = min_cells,
    concavity = concavity,
    window = window,
    workers = workers
  ),
  obj = prep
)

# ---- Step 2: Spatial statistics + bootstrap ----
t_stats_start <- Sys.time()
se_stats <- panoramic_spatialstats(
  prep = prep,
  pairs = "auto",
  radii_um = radii_um,
  stat = stat,
  nsim = nsim,
  correction = correction,
  seed = seed,
  boot = boot,
  tile_size = tile_size,
  BPPARAM = BPPARAM,
  verbose = FALSE
)
se_stats <- attach_sample_metadata(
  se = se_stats,
  sample_meta_df = sample_meta,
  patient_col = "patient",
  group_col = "group_label",
  region_col = "region",
  sample_meta_patient_col = "patient",
  sample_meta_group_col = "group_label"
)
saveRDS(se_stats, stats_path)
t_stats_end <- Sys.time()
append_runtime_metric(
  stage = "spatialstats",
  status = "computed",
  started_at = t_stats_start,
  ended_at = t_stats_end,
  cache_path = stats_path,
  params = list(
    stat = stat,
    radii_um = paste(radii_um, collapse = ","),
    nsim = nsim,
    correction = correction,
    seed = seed,
    boot = boot,
    tile_size = tile_size,
    workers = workers
  ),
  obj = se_stats
)

# ---- Step 3: Meta-analysis ----
t_meta_start <- Sys.time()
se_meta <- panoramic_meta_mv(
  se = se_stats,
  patient_col = patient_col,
  group_col = group_col,
  sample_col = "sample",
  tau_structure = tau_structure,
  method = method_mv,
  group_tau2 = group_tau2,
  vi_floor = vi_floor,
  BPPARAM = BPPARAM
)
saveRDS(se_meta, meta_path)
t_meta_end <- Sys.time()
append_runtime_metric(
  stage = "meta_analysis",
  status = "computed",
  started_at = t_meta_start,
  ended_at = t_meta_end,
  cache_path = meta_path,
  params = list(
    method = method_mv,
    tau_structure = tau_structure,
    group_tau2 = group_tau2,
    vi_floor = vi_floor,
    workers = workers
  ),
  obj = se_meta
)



# ---- Extract tables ----
spatial_tbl <- extract_spatialstats_table(se_stats)
meta_tbl <- extract_meta_table(se_meta) |> arrange(p_diff)
contrast_tbl <- extract_contrast_table(se_meta)
contrast_tbl <- contrast_tbl |> arrange(p_diff)

spatial_csv <- file.path(tab_dir, paste0(run_prefix, "_spatialstats.csv"))
meta_csv <- file.path(tab_dir, paste0(run_prefix, "_meta.csv"))
write_csv(spatial_tbl, spatial_csv)
write_csv(meta_tbl, meta_csv)

if (isTRUE(write_latex_results)) {
  alpha_slug <- gsub("\\.", "p", format(latex_alpha, scientific = FALSE, trim = TRUE))
  latex_full_path <- file.path(tab_dir, paste0(run_prefix, "_results_table.tex"))
  latex_selected_path <- if (isTRUE(latex_sig_only)) {
    file.path(tab_dir, paste0(run_prefix, "_results_table_sig", alpha_slug, ".tex"))
  } else {
    latex_full_path
  }

  write_panoramic_latex_results(
    meta_tbl = meta_tbl,
    out_path = latex_full_path,
    sig_only = FALSE,
    alpha = latex_alpha,
    caption = "CRC TMA PANORAMIC results.",
    print_console = FALSE
  )

  latex_out <- write_panoramic_latex_results(
    meta_tbl = meta_tbl,
    out_path = latex_selected_path,
    sig_only = latex_sig_only,
    alpha = latex_alpha,
    caption = if (isTRUE(latex_sig_only)) {
      paste0("CRC TMA PANORAMIC results (FDR <= ", format(latex_alpha, trim = TRUE), ").")
    } else {
      "CRC TMA PANORAMIC results."
    },
    print_console = latex_print_console
  )
  message("LaTeX results table: ", latex_out$path, " (rows=", latex_out$n_rows, ")")
}

contrast_csv <- file.path(tab_dir, paste0(run_prefix, "_contrast.csv"))
write_csv(contrast_tbl, contrast_csv)
message("Contrast table written: ", contrast_csv)

# ---- Figure outputs for paper panels ----
contrast_meta <- S4Vectors::metadata(se_meta)$panoramic$contrast
group_levels <- unique(as.character(as.data.frame(SummarizedExperiment::colData(se_stats))[[group_col]]))
group_levels <- group_levels[!is.na(group_levels)]
if (length(group_levels) != 2L) {
  stop("Expected exactly two groups in se_stats colData[['", group_col, "']].")
}
if (is.null(contrast_meta$control) || is.null(contrast_meta$case)) {
  stop("contrast metadata must include explicit control/case labels.")
}
control_group <- as.character(contrast_meta$control)
case_group <- as.character(contrast_meta$case)

volcano_prefix <- file.path(fig_dir, paste0(run_prefix, "_volcano"))
vol_fdr_plot <- panoramic::plot_volcano(
  se_diff = se_meta,
  x_scale = volcano_x_scale
)
vol_raw_plot <- panoramic::plot_volcano(
  se_diff = se_meta,
  sig_col = "p_diff",
  x_scale = volcano_x_scale
)
vol_fdr_plot <- force_volcano_geom_sizes(
  vol_fdr_plot,
  point_size = volcano_point_size
)
vol_raw_plot <- force_volcano_geom_sizes(
  vol_raw_plot,
  point_size = volcano_point_size
)
vol_fdr_png <- paste0(volcano_prefix, "_fdr_diff.png")
vol_fdr_pdf <- paste0(volcano_prefix, "_fdr_diff.pdf")
vol_raw_png <- paste0(volcano_prefix, "_p_diff.png")
vol_raw_pdf <- paste0(volcano_prefix, "_p_diff.pdf")
ggsave(vol_fdr_png, vol_fdr_plot, width = 9.2, height = 6.5, dpi = 300)
ggsave(vol_fdr_pdf, vol_fdr_plot, width = 9.2, height = 6.5)
ggsave(vol_raw_png, vol_raw_plot, width = 9.2, height = 6.5, dpi = 300)
ggsave(vol_raw_pdf, vol_raw_plot, width = 9.2, height = 6.5)
message("Volcano (FDR) PNG: ", vol_fdr_png)
message("Volcano (raw p) PNG: ", vol_raw_png)

# Publication-ready volcano (FDR): title 6 pt, subtitle/legend 5 pt, axes 4 pt; 2.25 in x 3 in.
vol_pub_plot <- vol_fdr_plot +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 6),
    plot.subtitle = ggplot2::element_text(size = 5),
    axis.title = ggplot2::element_text(size = 4),
    axis.text = ggplot2::element_text(size = 4),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 5),
    legend.position = "bottom"
  ) +
  ggplot2::labs(color = NULL)
vol_pub_png <- paste0(volcano_prefix, "_publication_fdr.png")
vol_pub_pdf <- paste0(volcano_prefix, "_publication_fdr.pdf")
ggsave(vol_pub_png, vol_pub_plot, width = 2.25, height = 3, units = "in", dpi = 300)
ggsave(vol_pub_pdf, vol_pub_plot, width = 2.25, height = 3, units = "in")
message("Volcano (publication FDR) PNG: ", vol_pub_png)

if (write_plot_debug) {
  vol_fdr <- build_volcano_df(contrast_tbl, control_group, case_group, sig_col = "fdr_diff")
  vol_raw <- build_volcano_df(contrast_tbl, control_group, case_group, sig_col = "p_diff")
  fdr_debug_csv <- paste0(volcano_prefix, "_fdr_debug.csv")
  raw_debug_csv <- paste0(volcano_prefix, "_raw_debug.csv")
  write_csv(vol_fdr, fdr_debug_csv)
  write_csv(vol_raw, raw_debug_csv)
}

rep_prefix <- file.path(fig_dir, paste0(run_prefix, "_representative_fdr"))
rep_out <- panoramic::plot_representative_samples(
  se_stats = se_stats,
  se_meta = se_meta,
  spe_list = spe_list,
  group_col = group_col,
  cell_type_col = cell_type_col,
  out_prefix = rep_prefix
)
if (!is.data.frame(rep_out$index) || nrow(rep_out$index) == 0L) {
  message("Representative panel generation returned no valid index rows; skipping representative sample index export.")
} else {
  rep_index_csv <- paste0(rep_prefix, "_index.csv")
  write_csv(rep_out$index, rep_index_csv)
  message("Representative hit panels written: ", nrow(rep_out$index))
  message("Representative panel index CSV: ", rep_index_csv)
}

save_network_outputs(se_meta = se_meta, fig_dir = fig_dir, run_prefix = run_prefix)

save_forest_outputs(
  se_meta = se_meta,
  contrast_tbl = contrast_tbl,
  fig_dir = fig_dir,
  run_prefix = run_prefix,
  group_col = group_col
)

# Optional supplement plot.
if (isTRUE(save_p_hist)) {
  p_hist_df <- contrast_tbl |>
    dplyr::transmute(p_value = as.numeric(.data$p_diff)) |>
    dplyr::filter(is.finite(.data$p_value), .data$p_value >= 0, .data$p_value <= 1)

  if (nrow(p_hist_df) == 0L) {
    stop("No finite p_diff values in [0,1] for p-value histogram.")
  }
  p_hist_plot <- ggplot(p_hist_df, aes(x = p_value)) +
    geom_histogram(bins = 15, fill = "#2C7FB8", color = "white") +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
    labs(
      title = "PANORAMIC p-value histogram",
      x = "p-value",
      y = "Count"
    ) +
    theme_bw(base_size = 11)

  p_hist_png <- file.path(fig_dir, paste0(run_prefix, "_pvalue_histogram.png"))
  p_hist_pdf <- file.path(fig_dir, paste0(run_prefix, "_pvalue_histogram.pdf"))
  ggsave(p_hist_png, p_hist_plot, width = 7, height = 5, dpi = 300)
  ggsave(p_hist_pdf, p_hist_plot, width = 7, height = 5)
  message("P-value histogram PNG: ", p_hist_png)
  message("P-value histogram PDF: ", p_hist_pdf)
}

write_runtime_metrics(runtime_csv)

message("Finished CRC TMA analysis.")
message("prep: ", prep_path)
message("se_stats: ", stats_path)
message("se_meta: ", meta_path)
message("spatial table: ", spatial_csv)
message("meta table: ", meta_csv)
message("runtime metrics: ", runtime_csv)
print(utils::head(meta_tbl, n = 10))
