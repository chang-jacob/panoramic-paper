#!/usr/bin/env Rscript

# HNSCC TMA Script: Primary PANORAMIC Analysis
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Run prepare -> spatial statistics/bootstrap -> multilevel meta-analysis.
# - Execute LN_benign vs LN_met_adjacent and PT Nneg vs Npos analysis branches.
# - Export spatial/meta/contrast tables, diagnostics, and manuscript-ready figures.

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
  # Column names in colData.
  columns = list(
    cell_type = "cell_type",
    patient = "patient",
    group = "core_category",
    region = "region"
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
    backend = "snow", # auto | multicore | snow | serial
    progressbar = TRUE
  ),

  # Meta-analysis.
  meta = list(
    engine = "mv",
    method_mv = "ML", # ML | REML
    group_tau2 = "none", # none | separate
    vi_floor = "group_median", # group_median | median | none
    tau_structure_mode = "auto", # auto | patient | patient_sample
    workers = 10L,
    backend = "snow", # auto | multicore | snow | serial
    progressbar = TRUE,
    warn_sigma2 = FALSE
  ),

  # Plot/table outputs.
  plotting = list(
    volcano_alpha = 0.05,
    volcano_label_top = 12L,
    volcano_label_pt = 4,
    rep_alpha = 0.05,
    rep_top_n = 10L,
    volcano_x_scale = "beta_diff", # beta_diff | log2fc
    network_fdr = 1.0,
    network_leiden_res = 1.2,
    network_directed = FALSE,
    network_layout = "fr",
    forest_top_n = 10L,
    forest_alpha = 0.05,
    forest_sig_col = "fdr_diff", # fdr_diff | p_diff
    save_p_hist = FALSE
  ),

  reporting = list(
    write_latex_results = TRUE,
    latex_sig_only = TRUE,
    latex_alpha = 0.05,
    latex_print_console = TRUE
  )
)

# ---- Inputs ----
dataset_tag <- "ln_benign_nneg_vs_ln_metadj_npos"
spe_list_path <- file.path(
  paths$data, "processed", "hnscc_tma", "hnscc_tma_LNbenign_Nneg_and_LNmetadj_Npos_spe_list.rds"
)
list_dir <- dataset_tag

# ---- Column names in colData(spe) ----
cell_type_col <- analysis_params$columns$cell_type
patient_col <- analysis_params$columns$patient
group_col <- analysis_params$columns$group
region_col <- analysis_params$columns$region
group_order <- c("LN_benign", "LN_met_adjacent")

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

# ---- Meta-analysis parameters ----
method_mv <- analysis_params$meta$method_mv
group_tau2 <- analysis_params$meta$group_tau2
vi_floor <- analysis_params$meta$vi_floor
tau_structure_mode <- analysis_params$meta$tau_structure_mode
meta_warn_sigma2 <- analysis_params$meta$warn_sigma2

workers <- analysis_params$panoramic$workers
panoramic_progressbar <- analysis_params$panoramic$progressbar
BPPARAM <- BiocParallel::SnowParam(workers = workers, progressbar = panoramic_progressbar)

meta_workers <- analysis_params$meta$workers
meta_progressbar <- analysis_params$meta$progressbar
BPPARAM_meta <- BiocParallel::SnowParam(workers = meta_workers, progressbar = meta_progressbar)

# ---- Plot parameters ----
volcano_alpha <- analysis_params$plotting$volcano_alpha
volcano_label_top <- analysis_params$plotting$volcano_label_top
volcano_label_pt <- analysis_params$plotting$volcano_label_pt
rep_alpha <- analysis_params$plotting$rep_alpha
rep_top_n <- analysis_params$plotting$rep_top_n
volcano_x_scale <- analysis_params$plotting$volcano_x_scale
network_fdr <- analysis_params$plotting$network_fdr
network_leiden_res <- analysis_params$plotting$network_leiden_res
network_directed <- analysis_params$plotting$network_directed
network_layout <- analysis_params$plotting$network_layout
forest_top_n <- analysis_params$plotting$forest_top_n
forest_alpha <- analysis_params$plotting$forest_alpha
forest_sig_col <- analysis_params$plotting$forest_sig_col
save_p_hist <- analysis_params$plotting$save_p_hist

# ---- Reporting parameters ----
write_latex_results <- analysis_params$reporting$write_latex_results
latex_sig_only <- analysis_params$reporting$latex_sig_only
latex_alpha <- analysis_params$reporting$latex_alpha
latex_print_console <- analysis_params$reporting$latex_print_console

# ---- Outputs ----
interim_dir <- file.path(paths$data, "interim")
tab_dir <- file.path(paths$output, "tables", "hnscc_tma", list_dir)
fig_dir <- file.path(paths$output, "figures", "hnscc_tma", list_dir)
dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stat_tag <- gsub("[^a-z0-9]+", "_", tolower(stat))
radii_tag <- paste(radii_um, collapse = "_")
tile_tag <- format(signif(tile_size, 6), scientific = FALSE, trim = TRUE)
meta_tag <- paste0("multilevel_", tolower(method_mv))
group_col_tag <- gsub("[^a-z0-9]+", "_", tolower(group_col))
encode_group_token <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("\\+", "pos", x)
  x <- gsub("-", "neg", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x
}
group_order_tag <- paste(vapply(group_order, encode_group_token, FUN.VALUE = character(1)), collapse = "_vs_")
run_prefix <- paste0(
  "hnscc_tma_", dataset_tag,
  "_", meta_tag,
  "_gcol", group_col_tag,
  "_gcmp", group_order_tag,
  "_stat", stat_tag,
  "_r", radii_tag,
  "_t", tile_tag
)
prep_path <- file.path(interim_dir, paste0(run_prefix, "_prep.rds"))
stats_path <- file.path(interim_dir, paste0(run_prefix, "_se_stats.rds"))
meta_path <- file.path(interim_dir, paste0(run_prefix, "_se_meta.rds"))
runtime_csv <- file.path(tab_dir, paste0(run_prefix, "_runtime_metrics.csv"))

# ---- Helpers ----
runtime_tracker <- create_runtime_tracker(run_prefix)
append_runtime_metric <- runtime_tracker$append_runtime_metric
write_runtime_metrics <- runtime_tracker$write_runtime_metrics

choose_tau_structure <- function(patient_vec) {
  per_patient <- table(patient_vec)
  if (any(per_patient > 1L)) "patient_sample" else "patient"
}

sanitize_slug <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "NA")
}

# ---- Read preprocessed data ----
message("Reading preprocessed data: ", spe_list_path)
spe_list <- readRDS(spe_list_path)

sample_meta <- lapply(spe_list, function(spe) {
  region_val <- get_unique(spe, region_col)
  data.frame(
    patient = get_unique(spe, patient_col),
    group = get_unique(spe, group_col),
    region = region_val,
    stringsAsFactors = FALSE
  )
})
sample_meta <- do.call(rbind, sample_meta)
sample_meta$sample <- sample_meta$region
sample_meta$group <- factor(as.character(sample_meta$group), levels = group_order)
names(spe_list) <- sample_meta$sample

message("Applied comparison groups/order: ", paste(group_order, collapse = " -> "))

spe_list <- spe_list[sample_meta$sample]
names(spe_list) <- sample_meta$sample

design <- data.frame(
  sample = sample_meta$sample,
  group = sample_meta$group,
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
  params = list(min_cells = min_cells, concavity = concavity, window = window, workers = workers),
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
  patient_col = patient_col,
  group_col = group_col,
  region_col = region_col
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
patient_sample_counts <- table(sample_meta$patient)
n_patients <- length(patient_sample_counts)
n_multi_sample_patients <- sum(patient_sample_counts > 1L)
tau_structure_auto <- choose_tau_structure(sample_meta$patient)
tau_structure_used <- if (identical(tau_structure_mode, "auto")) tau_structure_auto else tau_structure_mode
message(
  "HNSCC nested structure: ",
  n_multi_sample_patients, "/", n_patients,
  " patients have >1 sample; tau_structure='", tau_structure_used, "'."
)

t_meta_start <- Sys.time()
se_meta <- panoramic::panoramic_meta_mv(
  se = se_stats,
  patient_col = patient_col,
  group_col = group_col,
  sample_col = region_col,
  tau_structure = tau_structure_used,
  method = method_mv,
  group_tau2 = group_tau2,
  vi_floor = vi_floor,
  warn_sigma2 = meta_warn_sigma2,
  BPPARAM = BPPARAM_meta
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
    tau_structure = tau_structure_used,
    group_tau2 = group_tau2,
    vi_floor = vi_floor,
    warn_sigma2 = meta_warn_sigma2,
    n_patients = n_patients,
    n_multi_sample_patients = n_multi_sample_patients,
    workers = meta_workers
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
contrast_csv <- file.path(tab_dir, paste0(run_prefix, "_contrast.csv"))
write_csv(spatial_tbl, spatial_csv)
write_csv(meta_tbl, meta_csv)
write_csv(contrast_tbl, contrast_csv)
message("Contrast table written: ", contrast_csv)

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
    caption = "HNSCC TMA PANORAMIC results.",
    print_console = FALSE
  )

  latex_out <- write_panoramic_latex_results(
    meta_tbl = meta_tbl,
    out_path = latex_selected_path,
    sig_only = latex_sig_only,
    alpha = latex_alpha,
    caption = if (isTRUE(latex_sig_only)) {
      paste0("HNSCC TMA PANORAMIC results (FDR <= ", format(latex_alpha, trim = TRUE), ").")
    } else {
      "HNSCC TMA PANORAMIC results."
    },
    print_console = latex_print_console
  )
  message("LaTeX results table: ", latex_out$path, " (rows=", latex_out$n_rows, ")")
}

# ---- Figure outputs for paper panels ----
volcano_fdr <- panoramic::plot_volcano(
  se_diff = se_meta,
  sig_col = "fdr_diff",
  alpha = volcano_alpha,
  x_scale = volcano_x_scale,
  label_top = volcano_label_top,
  label_text_pt = volcano_label_pt,
  title = "HNSCC TMA colocalization volcano"
)
volcano_raw <- panoramic::plot_volcano(
  se_diff = se_meta,
  sig_col = "p_diff",
  alpha = volcano_alpha,
  x_scale = volcano_x_scale,
  label_top = volcano_label_top,
  label_text_pt = volcano_label_pt,
  title = "HNSCC TMA colocalization volcano"
)

volcano_fdr_png <- file.path(fig_dir, paste0(run_prefix, "_volcano_fdr_diff.png"))
volcano_fdr_pdf <- file.path(fig_dir, paste0(run_prefix, "_volcano_fdr_diff.pdf"))
volcano_raw_png <- file.path(fig_dir, paste0(run_prefix, "_volcano_p_diff.png"))
volcano_raw_pdf <- file.path(fig_dir, paste0(run_prefix, "_volcano_p_diff.pdf"))

ggsave(volcano_fdr_png, volcano_fdr, width = 9.2, height = 6.5, dpi = 300)
ggsave(volcano_fdr_pdf, volcano_fdr, width = 9.2, height = 6.5)
ggsave(volcano_raw_png, volcano_raw, width = 9.2, height = 6.5, dpi = 300)
ggsave(volcano_raw_pdf, volcano_raw, width = 9.2, height = 6.5)

message("Volcano (FDR) PNG: ", volcano_fdr_png)
message("Volcano (raw p) PNG: ", volcano_raw_png)

rep_prefix <- file.path(fig_dir, paste0(run_prefix, "_representative_samples_fdr"))
rep_out <- panoramic::plot_representative_samples(
  se_stats = se_stats,
  se_meta = se_meta,
  spe_list = spe_list,
  sig_col = "fdr_diff",
  alpha = rep_alpha,
  top_n = rep_top_n,
  group_col = group_col,
  cell_type_col = cell_type_col,
  sample_col = "sample",
  out_prefix = rep_prefix
)
rep_index_csv <- paste0(rep_prefix, "_index.csv")
write_csv(rep_out$index, rep_index_csv)
message("Representative hit panels written: ", nrow(rep_out$index))
message("Representative panel index CSV: ", rep_index_csv)

network_prefix <- file.path(fig_dir, paste0(run_prefix, "_colocalization_network"))
network_result <- panoramic::create_spatial_network(
  se_diff = se_meta,
  fdr_threshold = network_fdr,
  directed = network_directed,
  leiden_resolution = network_leiden_res
)
network_plot <- panoramic::plot_spatial_network(
  net_result = network_result,
  layout = network_layout
) +
  ggplot2::labs(title = NULL, subtitle = NULL) +
  ggplot2::theme(legend.position = "none")
network_png <- paste0(network_prefix, ".png")
network_pdf <- paste0(network_prefix, ".pdf")
network_rds <- paste0(network_prefix, ".rds")
network_summary_csv <- paste0(network_prefix, "_summary.csv")
ggsave(network_png, network_plot, width = 8.5, height = 7, dpi = 300)
ggsave(network_pdf, network_plot, width = 8.5, height = 7)
saveRDS(network_result, network_rds)
readr::write_csv(
  data.frame(
    n_clusters = network_result$n_clusters,
    modularity = network_result$modularity,
    n_nodes = igraph::vcount(network_result$graph),
    n_edges = igraph::ecount(network_result$graph),
    stringsAsFactors = FALSE
  ),
  network_summary_csv
)

forest_hits <- contrast_tbl |>
  dplyr::arrange(.data[[forest_sig_col]], p_diff) |>
  dplyr::distinct(ct1, ct2, radius_um, .keep_all = TRUE) |>
  dplyr::slice_head(n = forest_top_n)

forest_dir <- file.path(fig_dir, paste0(run_prefix, "_forest_plots"))
dir.create(forest_dir, recursive = TRUE, showWarnings = FALSE)
forest_index_rows <- vector("list", nrow(forest_hits))
for (i in seq_len(nrow(forest_hits))) {
  hit <- forest_hits[i, , drop = FALSE]
  ct1 <- as.character(hit$ct1)
  ct2 <- as.character(hit$ct2)
  radius_um <- as.numeric(hit$radius_um)
  radius_slug <- gsub("\\.", "p", format(radius_um, trim = TRUE, scientific = FALSE))
  forest_stub <- file.path(
    forest_dir,
    paste0(
      sprintf("%03d_", i),
      sanitize_slug(ct1), "_to_", sanitize_slug(ct2),
      "_r", radius_slug, "um"
    )
  )
  forest_png <- paste0(forest_stub, ".png")
  forest_pdf <- paste0(forest_stub, ".pdf")
  forest_plot <- panoramic::plot_forest(
    se_meta = se_meta,
    ct1 = ct1,
    ct2 = ct2,
    radius_um = radius_um,
    group_col = group_col
  ) +
    ggplot2::guides(color = ggplot2::guide_legend(title = NULL, override.aes = list(shape = 16, size = 2))) +
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 5))
  ggsave(forest_png, forest_plot, width = 12, height = 4, dpi = 300)
  ggsave(forest_pdf, forest_plot, width = 12, height = 4)
  forest_index_rows[[i]] <- data.frame(
    rank = i,
    ct1 = ct1,
    ct2 = ct2,
    radius_um = radius_um,
    p_diff = hit$p_diff[1],
    fdr_diff = hit$fdr_diff[1],
    sig_col = forest_sig_col,
    sig_value = hit[[forest_sig_col]][1],
    png = forest_png,
    pdf = forest_pdf,
    stringsAsFactors = FALSE
  )
}
forest_index_df <- dplyr::bind_rows(forest_index_rows)
forest_index_csv <- file.path(forest_dir, "forest_plot_index.csv")
readr::write_csv(forest_index_df, forest_index_csv)

# Optional supplement plot.
if (isTRUE(save_p_hist) && !is.null(contrast_tbl) && "p_diff" %in% colnames(contrast_tbl)) {
  p_hist_path <- file.path(fig_dir, paste0(run_prefix, "_pvalue_histogram.png"))
  png(p_hist_path, width = 1400, height = 1000, res = 200)
  hist(contrast_tbl$p_diff, main = "HNSCC contrast p-value histogram", xlab = "p-value")
  dev.off()
  message("P-value histogram PNG: ", p_hist_path)
}

write_runtime_metrics(runtime_csv)

message("Finished HNSCC TMA analysis.")
message("prep: ", prep_path)
message("se_stats: ", stats_path)
message("se_meta: ", meta_path)
message("spatial table: ", spatial_csv)
message("meta table: ", meta_csv)
message("runtime metrics: ", runtime_csv)
print(utils::head(meta_tbl, n = 10))

hist(contrast_tbl$p_diff)
