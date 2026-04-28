#!/usr/bin/env Rscript

# Simulation Script: Shared Utilities
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Provide shared path, I/O, parallel, simulation, and extraction helpers used by simulation workflows.
# - Build SpatialExperiment/SummarizedExperiment objects for intrasample and intersample runs.
# - Compute shared calibration summaries and locate CRC TMA reference inputs.

suppressPackageStartupMessages({
  library(config)
  library(panoramic)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(spatstat.geom)
  library(spatstat.random)
  library(BiocParallel)
  library(dplyr)
})

sim_find_project_root <- function(start = getwd(), max_depth = 8L) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_depth)) {
    cfg <- file.path(cur, "config", "default.yml")
    if (file.exists(cfg)) return(cur)
    nxt <- dirname(cur)
    if (identical(nxt, cur)) break
    cur <- nxt
  }
  stop("Could not locate project root containing config/default.yml from: ", start)
}

sim_load_paths <- function(project_root) {
  config::get("paths", file = file.path(project_root, "config", "default.yml"))
}

sim_mkdir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

sim_safe_save_rds <- function(object, file) {
  sim_mkdir(dirname(file))
  saveRDS(object, file = file)
  invisible(file)
}

sim_safe_write_csv <- function(df, file) {
  sim_mkdir(dirname(file))
  utils::write.csv(df, file = file, row.names = FALSE)
  invisible(file)
}

sim_clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))

sim_make_bpparam <- function(n_workers, backend = NULL) {
  n_workers <- as.integer(n_workers)
  if (!is.finite(n_workers) || n_workers < 1L) {
    stop("n_workers must be >= 1")
  }

  backend <- if (is.null(backend)) "" else tolower(trimws(as.character(backend)))
  if (!nzchar(backend)) {
    in_rstudio <- nzchar(Sys.getenv("RSTUDIO", unset = ""))
    backend <- if ((.Platform$OS.type == "unix") && !in_rstudio) "multicore" else "snow"
  }

  if (!backend %in% c("multicore", "snow", "serial")) {
    stop("backend must be one of: multicore, snow, serial")
  }

  bp <- if (n_workers <= 1L || identical(backend, "serial")) {
    BiocParallel::SerialParam(progressbar = TRUE)
  } else if (identical(backend, "multicore")) {
    BiocParallel::MulticoreParam(workers = n_workers, progressbar = TRUE)
  } else {
    BiocParallel::SnowParam(workers = n_workers, type = "SOCK", progressbar = TRUE)
  }

  BiocParallel::bpstopOnError(bp) <- FALSE
  BiocParallel::bpexportglobals(bp) <- TRUE
  BiocParallel::bpexportvariables(bp) <- TRUE
  bp
}

sim_draw_n_samples_per_patient <- function(cfg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- if (identical(cfg$n_samples_dist, "nbinom")) {
    stats::rnbinom(1, size = cfg$n_samples_nbinom_size, mu = cfg$n_samples_per_patient_mean)
  } else {
    stats::rpois(1, lambda = cfg$n_samples_per_patient_mean)
  }
  n <- max(cfg$n_samples_per_patient_min, as.integer(n))
  n <- min(cfg$n_samples_per_patient_max, as.integer(n))
  as.integer(n)
}

sim_ppp_to_spe <- function(ppp, sample_id = "sample_1", cell_type_col = "cell_type") {
  raw_df <- as.data.frame(ppp)
  coords <- as.matrix(raw_df[, c("x", "y"), drop = FALSE])
  cell_ids <- paste0(sample_id, "_cell_", seq_len(nrow(coords)))
  rownames(coords) <- cell_ids

  col_data <- S4Vectors::DataFrame(
    cell_id = cell_ids,
    sample = sample_id,
    stringsAsFactors = FALSE
  )
  col_data[[cell_type_col]] <- as.character(raw_df$marks)
  rownames(col_data) <- cell_ids

  counts <- matrix(
    0L,
    nrow = 1L,
    ncol = nrow(coords),
    dimnames = list("placeholder", cell_ids)
  )

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = col_data,
    spatialCoords = coords
  )
}

sim_simulate_pattern_ppp <- function(
    pattern = c("uniform", "opposing_gradient", "clustered"),
    square_size = 2000,
    cell_types = c("A", "B", "C", "Rare 1", "Rare 2"),
    params = list(),
    seed = NULL
) {
  pattern <- match.arg(pattern)
  if (!is.null(seed)) set.seed(seed)

  win <- spatstat.geom::square(square_size)
  stopifnot(length(cell_types) >= 2L)
  if (pattern %in% c("opposing_gradient", "clustered") && length(cell_types) < 5L) {
    stop("Patterns 'opposing_gradient' and 'clustered' require at least 5 cell types.")
  }

  if (identical(pattern, "uniform")) {
    lambda <- params$lambda
    if (is.null(lambda)) {
      lambda <- c(A = 2e-4, B = 1.5e-4, C = 1e-4, Rare1 = 1e-5, Rare2 = 5e-6)
    }
    lambda <- lambda[seq_len(length(cell_types))]
    names(lambda) <- cell_types
    pp <- spatstat.random::rmpoispp(lambda = as.numeric(lambda), win = win, types = cell_types)
    return(pp)
  }

  if (identical(pattern, "opposing_gradient")) {
    rate <- params$gradient_rate
    if (is.null(rate)) rate <- 0.0015
    lambda_c <- if (is.null(params$lambda_c)) 1e-4 else params$lambda_c
    lambda_r1 <- if (is.null(params$lambda_rare1)) 1e-5 else params$lambda_rare1
    lambda_r2 <- if (is.null(params$lambda_rare2)) 5e-6 else params$lambda_rare2
    target_a <- if (is.null(params$target_density_a)) NA_real_ else as.numeric(params$target_density_a)
    target_b <- if (is.null(params$target_density_b)) NA_real_ else as.numeric(params$target_density_b)

    mean_dexp <- if (rate > 0) {
      (1 - exp(-rate * square_size)) / (rate * square_size)
    } else {
      1
    }
    scale_a <- if (is.finite(target_a)) target_a / max(mean_dexp, 1e-10) else 1
    scale_b <- if (is.finite(target_b)) target_b / max(mean_dexp, 1e-10) else 1

    # Opposing-gradient profiles are proportional to exp(-rate * x), then scaled
    # so spatially averaged density matches target_a/target_b.
    la <- function(x, y) scale_a * exp(-rate * x)
    lb <- function(x, y) scale_b * exp(-rate * (square_size - x))

    pp_a <- spatstat.random::rpoispp(lambda = la, win = win)
    pp_b <- spatstat.random::rpoispp(lambda = lb, win = win)
    pp_c <- spatstat.random::rpoispp(lambda = lambda_c, win = win)
    pp_r1 <- spatstat.random::rpoispp(lambda = lambda_r1, win = win)
    pp_r2 <- spatstat.random::rpoispp(lambda = lambda_r2, win = win)

    pp_a$marks <- factor(rep(cell_types[1], pp_a$n), levels = cell_types)
    pp_b$marks <- factor(rep(cell_types[2], pp_b$n), levels = cell_types)
    pp_c$marks <- factor(rep(cell_types[3], pp_c$n), levels = cell_types)
    pp_r1$marks <- factor(rep(cell_types[4], pp_r1$n), levels = cell_types)
    pp_r2$marks <- factor(rep(cell_types[5], pp_r2$n), levels = cell_types)

    pp <- spatstat.geom::superimpose(pp_a, pp_b, pp_c, pp_r1, pp_r2, W = win)
    spatstat.geom::marks(pp) <- factor(spatstat.geom::marks(pp), levels = cell_types)
    return(pp)
  }

  parent_lambda <- if (is.null(params$parent_lambda)) 4e-5 else params$parent_lambda
  offspring_mean <- if (is.null(params$offspring_mean)) 7 else params$offspring_mean
  offspring_sd <- if (is.null(params$offspring_sd)) 175 else params$offspring_sd
  lambda_c <- if (is.null(params$lambda_c)) 1e-4 else params$lambda_c
  lambda_r1 <- if (is.null(params$lambda_rare1)) 1e-5 else params$lambda_rare1
  lambda_r2 <- if (is.null(params$lambda_rare2)) 5e-6 else params$lambda_rare2

  parents <- spatstat.random::rpoispp(lambda = parent_lambda, win = win)
  off_x <- numeric()
  off_y <- numeric()

  if (parents$n > 0) {
    n_off <- stats::rpois(parents$n, offspring_mean) + 1L
    tmp_x <- stats::rnorm(sum(n_off), mean = rep(parents$x, n_off), sd = offspring_sd / 2)
    tmp_y <- stats::rnorm(sum(n_off), mean = rep(parents$y, n_off), sd = offspring_sd / 2)
    keep <- tmp_x >= 0 & tmp_x <= square_size & tmp_y >= 0 & tmp_y <= square_size
    off_x <- tmp_x[keep]
    off_y <- tmp_y[keep]
  }

  pp_a <- spatstat.geom::ppp(
    x = parents$x,
    y = parents$y,
    window = win,
    marks = factor(rep(cell_types[1], parents$n), levels = cell_types)
  )
  pp_b <- spatstat.geom::ppp(
    x = off_x,
    y = off_y,
    window = win,
    marks = factor(rep(cell_types[2], length(off_x)), levels = cell_types)
  )
  pp_c <- spatstat.random::rpoispp(lambda = lambda_c, win = win)
  pp_r1 <- spatstat.random::rpoispp(lambda = lambda_r1, win = win)
  pp_r2 <- spatstat.random::rpoispp(lambda = lambda_r2, win = win)
  pp_c$marks <- factor(rep(cell_types[3], pp_c$n), levels = cell_types)
  pp_r1$marks <- factor(rep(cell_types[4], pp_r1$n), levels = cell_types)
  pp_r2$marks <- factor(rep(cell_types[5], pp_r2$n), levels = cell_types)

  pp <- spatstat.geom::superimpose(pp_a, pp_b, pp_c, pp_r1, pp_r2, W = win)
  spatstat.geom::marks(pp) <- factor(spatstat.geom::marks(pp), levels = cell_types)
  pp
}

sim_simulate_ab_cluster_ppp <- function(
    square_size = 600,
    parent_lambda = 3e-4,
    offspring_mean = 14,
    offspring_sd = 55,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  win <- spatstat.geom::square(square_size)
  parent_lambda <- as.numeric(parent_lambda)
  offspring_mean <- as.numeric(offspring_mean)
  offspring_sd <- as.numeric(offspring_sd)

  pp_parents <- spatstat.random::rpoispp(lambda = parent_lambda, win = win)

  off_x <- numeric()
  off_y <- numeric()
  if (pp_parents$n > 0L) {
    n_off <- stats::rpois(pp_parents$n, lambda = pmax(offspring_mean, 0))
    if (sum(n_off) > 0L) {
      tmp_x <- stats::rnorm(sum(n_off), mean = rep(pp_parents$x, n_off), sd = offspring_sd)
      tmp_y <- stats::rnorm(sum(n_off), mean = rep(pp_parents$y, n_off), sd = offspring_sd)
      keep <- tmp_x >= 0 & tmp_x <= square_size & tmp_y >= 0 & tmp_y <= square_size
      off_x <- tmp_x[keep]
      off_y <- tmp_y[keep]
    }
  }

  pp_a <- spatstat.geom::ppp(
    x = pp_parents$x,
    y = pp_parents$y,
    window = win,
    marks = factor(rep("A", pp_parents$n), levels = c("A", "B"))
  )
  pp_b <- spatstat.geom::ppp(
    x = off_x,
    y = off_y,
    window = win,
    marks = factor(rep("B", length(off_x)), levels = c("A", "B"))
  )

  pp <- spatstat.geom::superimpose(pp_a, pp_b, W = win)
  spatstat.geom::marks(pp) <- factor(spatstat.geom::marks(pp), levels = c("A", "B"))
  pp
}

sim_run_spatialstats_one <- function(
    spe,
    sample_id = "sample_1",
    cell_type_col = "cell_type",
    radii_um = 25,
    stat = "local_comp_enrichment",
    nsim = 100,
    correction = "translate",
    min_cells = 1L,
    concavity = 50,
    window = "rect",
    seed = 123,
    boot = "block",
    tile_size = 2.5 * max(radii_um),
    BPPARAM = BiocParallel::SerialParam()
) {
  design <- data.frame(sample = sample_id, group = "sim", stringsAsFactors = FALSE)
  prep <- panoramic::panoramic_prepare(
    spe_list = setNames(list(spe), sample_id),
    design = design,
    cell_type = cell_type_col,
    min_cells = min_cells,
    concavity = concavity,
    window = window,
    BPPARAM = BPPARAM
  )

  panoramic::panoramic_spatialstats(
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
}

sim_flatten_spatialstats <- function(se_stats) {
  out <- panoramic::panoramic_extract_contrast(se_stats, what = "spatialstats")
  out$feature_id <- paste(out$ct1, out$ct2, out$radius_um, out$stat, sep = "|")
  out
}

sim_build_feature_truth <- function(df_long) {
  df_long %>%
    group_by(pattern, feature_id, ct1, ct2, radius_um, stat) %>%
    summarise(
      empirical_mean = mean(yi, na.rm = TRUE),
      empirical_variance = stats::var(yi, na.rm = TRUE),
      mean_bootstrap_variance = mean(vi, na.rm = TRUE),
      n_reps = dplyr::n(),
      .groups = "drop"
    )
}

sim_join_truth <- function(df_long, truth_df) {
  dplyr::left_join(
    df_long,
    truth_df %>% select(pattern, feature_id, empirical_mean, empirical_variance, mean_bootstrap_variance),
    by = c("pattern", "feature_id")
  )
}

sim_variance_calibration_metrics <- function(df_joined) {
  feat <- df_joined %>%
    distinct(pattern, feature_id, empirical_variance, mean_bootstrap_variance) %>%
    filter(is.finite(empirical_variance), is.finite(mean_bootstrap_variance), empirical_variance >= 0, mean_bootstrap_variance >= 0)

  pear <- if (nrow(feat) > 2L) stats::cor(feat$empirical_variance, feat$mean_bootstrap_variance, method = "pearson") else NA_real_
  spear <- if (nrow(feat) > 2L) stats::cor(feat$empirical_variance, feat$mean_bootstrap_variance, method = "spearman") else NA_real_
  lrmse <- if (nrow(feat) > 0L) sqrt(mean((log1p(feat$mean_bootstrap_variance) - log1p(feat$empirical_variance))^2, na.rm = TRUE)) else NA_real_

  slope <- NA_real_
  if (nrow(feat) > 2L && any(feat$mean_bootstrap_variance > 0, na.rm = TRUE)) {
    fit <- stats::lm(empirical_variance ~ 0 + mean_bootstrap_variance, data = feat)
    slope <- as.numeric(stats::coef(fit)[1])
  }

  cover_df <- df_joined %>%
    mutate(
      se_hat = sqrt(pmax(vi, 0)),
      lower_95 = yi - 1.96 * se_hat,
      upper_95 = yi + 1.96 * se_hat,
      covered_95 = is.finite(empirical_mean) & is.finite(lower_95) & is.finite(upper_95) &
        empirical_mean >= lower_95 & empirical_mean <= upper_95
    )
  coverage <- mean(cover_df$covered_95, na.rm = TRUE)

  data.frame(
    pearson_feature_correlation = pear,
    spearman_feature_correlation = spear,
    log_rmse_variance = lrmse,
    variance_calibration_slope = slope,
    mean_interval_coverage_95 = coverage,
    n_feature_rows = nrow(feat),
    n_total_rows = nrow(df_joined),
    stringsAsFactors = FALSE
  )
}

sim_per_replicate_variance_correlation <- function(df_long, truth_df) {
  df_long %>%
    left_join(truth_df %>% select(pattern, feature_id, empirical_variance), by = c("pattern", "feature_id")) %>%
    group_by(pattern, replicate) %>%
    summarise(
      pearson_corr = stats::cor(vi, empirical_variance, use = "complete.obs", method = "pearson"),
      spearman_corr = stats::cor(vi, empirical_variance, use = "complete.obs", method = "spearman"),
      .groups = "drop"
    )
}

sim_effects_to_se <- function(df_effects) {
  feat <- df_effects %>%
    distinct(feature_id, ct1, ct2, radius_um, stat) %>%
    arrange(feature_id)
  col_keep <- c("sample", "patient", "group", "tissue_size_um", "tissue_area_um2")
  col_keep <- col_keep[col_keep %in% colnames(df_effects)]
  samp <- df_effects %>%
    distinct(dplyr::across(all_of(col_keep))) %>%
    arrange(sample)

  yi <- matrix(NA_real_, nrow = nrow(feat), ncol = nrow(samp), dimnames = list(feat$feature_id, samp$sample))
  vi <- matrix(NA_real_, nrow = nrow(feat), ncol = nrow(samp), dimnames = list(feat$feature_id, samp$sample))

  rix <- match(df_effects$feature_id, feat$feature_id)
  cix <- match(df_effects$sample, samp$sample)
  yi[cbind(rix, cix)] <- df_effects$yi
  vi[cbind(rix, cix)] <- df_effects$vi

  SummarizedExperiment::SummarizedExperiment(
    assays = list(yi = yi, vi = vi),
    rowData = S4Vectors::DataFrame(
      feature_id = feat$feature_id,
      ct1 = feat$ct1,
      ct2 = feat$ct2,
      radius_um = feat$radius_um,
      stat = feat$stat
    ),
    colData = S4Vectors::DataFrame(
      sample = samp$sample,
      patient = samp$patient,
      group = samp$group,
      tissue_size_um = samp$tissue_size_um,
      tissue_area_um2 = samp$tissue_area_um2
    )
  )
}

sim_extract_meta_estimates <- function(se_meta) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se_meta))
  keep <- intersect(
    c("feature_id", "ct1", "ct2", "radius_um", "stat", "mu_hat", "se_mu", "p_mu", "k", "tau2_patient", "tau2_sample"),
    colnames(rd)
  )
  rd[, keep, drop = FALSE]
}

sim_crc_tma_spe_list_path <- function(paths, key = "one_per_patient") {
  spe_list_paths <- list(
    all = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list.rds"),
    one_per_patient = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_one_per_patient.rds"),
    core_filtered = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_core_filtered.rds"),
    diffuse_only = file.path(paths$data, "processed", "crc_tma", "crc_tma_spe_list_diffuse_only.rds")
  )
  spe_list_paths[[key]]
}
