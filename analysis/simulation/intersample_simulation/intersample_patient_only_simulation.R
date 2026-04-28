#!/usr/bin/env Rscript

# Simulation Script: Intersample Patient-Only Simulation
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Simulate one-sample-per-patient A/B spatial effects under controlled heterogeneity.
# - Estimate parameters with PANORAMIC and naive baselines over progressive patient subsets.
# - Export truth, observed data, estimates, and patient-order artifacts.

suppressPackageStartupMessages({
  library(panoramic)
  library(BiocParallel)
  library(dplyr)
})

source("analysis/simulation/sim_utils.R")

project_root <- sim_find_project_root()
paths <- sim_load_paths(project_root)

out_data_dir <- sim_mkdir(file.path(paths$data, "interim", "simulation", "intersample_patient_only"))
out_tbl_dir <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intersample_patient_only"))

cfg <- list(
  cache_tag = "intersample_patient_only",
  n_reps = 20L,
  n_patients_max = 50L,
  n_patients_grid = c(10L, 20L, 30L, 40L, 50L),
  n_truth_realizations_per_patient = 80L,
  square_size_um = 600,
  target_cells = 1900,
  radius_um = 25,
  stat = "local_comp_enrichment",
  boot_nsim = 60L,
  boot_mode = "block",
  tile_size = 62.5,
  parent_lambda_min = 1e-5,
  parent_lambda_span = 1e-4,
  offspring_nb_size = 10,
  offspring_nb_mu = 2,
  offspring_count_offset = 3,
  offspring_sd_center_at_2000 = 200,
  offspring_sd_sd_at_2000 = 100,
  patient_lambda_jitter_sd = 0.25,
  sample_param_jitter_sd = 0.25,
  focus_feature_id = "A|B|25|local_comp_enrichment",
  seed_base = 4242L,
  n_workers = 10L
)

cfg$n_reps <- as.integer(cfg$n_reps)
cfg$n_patients_max <- as.integer(cfg$n_patients_max)
cfg$n_patients_grid <- sort(unique(as.integer(cfg$n_patients_grid)))
cfg$n_truth_realizations_per_patient <- as.integer(cfg$n_truth_realizations_per_patient)
cfg$square_size_um <- as.numeric(cfg$square_size_um)
cfg$target_cells <- as.numeric(cfg$target_cells)
cfg$radius_um <- as.numeric(cfg$radius_um)
cfg$stat <- as.character(cfg$stat)
cfg$boot_nsim <- as.integer(cfg$boot_nsim)
cfg$boot_mode <- as.character(cfg$boot_mode)
cfg$tile_size <- as.numeric(cfg$tile_size)
cfg$parent_lambda_min <- as.numeric(cfg$parent_lambda_min)
cfg$parent_lambda_span <- as.numeric(cfg$parent_lambda_span)
cfg$offspring_nb_size <- as.numeric(cfg$offspring_nb_size)
cfg$offspring_nb_mu <- as.numeric(cfg$offspring_nb_mu)
cfg$offspring_count_offset <- as.numeric(cfg$offspring_count_offset)
cfg$offspring_sd_center_at_2000 <- as.numeric(cfg$offspring_sd_center_at_2000)
cfg$offspring_sd_sd_at_2000 <- as.numeric(cfg$offspring_sd_sd_at_2000)
cfg$patient_lambda_jitter_sd <- as.numeric(cfg$patient_lambda_jitter_sd)
cfg$sample_param_jitter_sd <- as.numeric(cfg$sample_param_jitter_sd)
cfg$focus_feature_id <- as.character(cfg$focus_feature_id)
cfg$seed_base <- as.integer(cfg$seed_base)
cfg$n_workers <- as.integer(cfg$n_workers)
cfg$cache_tag <- as.character(cfg$cache_tag)
cfg$target_cell_density <- cfg$target_cells / (cfg$square_size_um^2)

if (cfg$n_workers < 1L) stop("n_workers must be >= 1")
if (cfg$n_patients_max < 2L) stop("n_patients_max must be >= 2")
if (min(cfg$n_patients_grid) < 2L) stop("n_patients_grid must be >= 2")
if (max(cfg$n_patients_grid) > cfg$n_patients_max) {
  stop("max(n_patients_grid) cannot exceed n_patients_max")
}
if (!identical(cfg$stat, "local_comp_enrichment")) {
  stop("This controlled fallback simulation is fixed to stat = 'local_comp_enrichment'.")
}

sim_draw_patient_latents <- function(n_patients, cfg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  patient_ids <- sprintf("patient_%03d", seq_len(n_patients))

  old_parent_mean <- cfg$parent_lambda_min + cfg$parent_lambda_span / 2
  old_offspring_mean <- cfg$offspring_count_offset + cfg$offspring_nb_mu
  target_parent_mean <- cfg$target_cell_density / (1 + old_offspring_mean)
  parent_scale <- target_parent_mean / old_parent_mean
  spatial_scale <- cfg$square_size_um / 2000

  parent_lambda_patient <- (
    cfg$parent_lambda_min + stats::runif(n_patients, 0, cfg$parent_lambda_span)
  ) * parent_scale * exp(stats::rnorm(n_patients, mean = 0, sd = cfg$patient_lambda_jitter_sd))

  offspring_mean_patient <- pmax(
    1,
    cfg$offspring_count_offset +
      stats::rnbinom(n_patients, size = cfg$offspring_nb_size, mu = cfg$offspring_nb_mu)
  )

  offspring_sd_patient <- pmax(
    5,
    abs(
      cfg$offspring_sd_center_at_2000 * spatial_scale +
        stats::rnorm(
          n_patients,
          mean = 0,
          sd = cfg$offspring_sd_sd_at_2000 * spatial_scale
        )
    )
  )

  data.frame(
    patient = patient_ids,
    parent_lambda_patient = parent_lambda_patient,
    offspring_mean_patient = offspring_mean_patient,
    offspring_sd_patient = offspring_sd_patient,
    stringsAsFactors = FALSE
  )
}

sim_eval_ab_sample <- function(patient_row, sample_id, cfg, seed, use_bootstrap = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  parent_lambda <- pmax(
    cfg$target_cell_density / 300,
    patient_row$parent_lambda_patient * exp(stats::rnorm(1, 0, cfg$sample_param_jitter_sd))
  )
  offspring_mean <- pmax(
    1,
    patient_row$offspring_mean_patient * exp(stats::rnorm(1, 0, cfg$sample_param_jitter_sd))
  )
  offspring_sd <- pmax(
    5,
    patient_row$offspring_sd_patient * exp(stats::rnorm(1, 0, cfg$sample_param_jitter_sd / 2))
  )

  pp <- sim_simulate_ab_cluster_ppp(
    square_size = cfg$square_size_um,
    parent_lambda = parent_lambda,
    offspring_mean = offspring_mean,
    offspring_sd = offspring_sd
  )

  if (isTRUE(use_bootstrap)) {
    spe <- sim_ppp_to_spe(pp, sample_id = sample_id, cell_type_col = "cell_type")
    se_stats <- sim_run_spatialstats_one(
      spe = spe,
      sample_id = sample_id,
      cell_type_col = "cell_type",
      radii_um = cfg$radius_um,
      stat = cfg$stat,
      nsim = cfg$boot_nsim,
      seed = seed,
      boot = cfg$boot_mode,
      tile_size = cfg$tile_size,
      BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
    )
    fx <- sim_flatten_spatialstats(se_stats)
    fx <- fx %>%
      dplyr::filter(
        ct1 %in% c("A", "B"),
        ct2 %in% c("A", "B"),
        radius_um == cfg$radius_um,
        stat == cfg$stat
      )
  } else {
    precompute_local_comp_cache <- getFromNamespace(".precompute_local_comp_cache", "panoramic")
    local_comp_score_fn <- getFromNamespace(".local_comp_score", "panoramic")
    cache <- precompute_local_comp_cache(pp, cfg$radius_um)
    ridx <- match(cfg$radius_um, cache$radii_um)
    if (is.na(ridx)) stop("Could not match radius in local_comp cache.")

    pairs <- expand.grid(ct1 = c("A", "B"), ct2 = c("A", "B"), stringsAsFactors = FALSE)
    calc_yi <- function(ct1, ct2) {
      tab <- table(as.character(spatstat.geom::marks(pp)))
      n1 <- unname(tab[ct1]); n2 <- unname(tab[ct2])
      if (is.na(n1) || n1 < 2L || is.na(n2) || n2 < 1L) return(NA_real_)

      anchor_idx <- cache$indices_by_type[[ct1]]
      to_idx <- match(ct2, cache$cell_types)
      if (length(anchor_idx) < 1L || is.na(to_idx)) return(NA_real_)

      frac <- cache$local_fraction[[ridx]][anchor_idx, to_idx]
      n_neighbors <- cache$local_total_neighbors[[ridx]][anchor_idx]
      overlap_area <- cache$local_overlap_area[anchor_idx, ridx]
      global_density <- cache$global_density[to_idx]
      if (identical(ct1, ct2) &&
          is.finite(cache$area_window) &&
          cache$area_window > 0) {
        n_to <- length(cache$indices_by_type[[ct2]])
        global_density <- max((n_to - 1) / cache$area_window, 0)
      }
      score <- local_comp_score_fn(
        frac = frac,
        n_neighbors = n_neighbors,
        overlap_area = overlap_area,
        global_density = global_density
      )
      yy <- mean(score, na.rm = TRUE)
      if (is.nan(yy)) NA_real_ else as.numeric(yy)
    }

    fx <- pairs
    fx$radius_um <- cfg$radius_um
    fx$stat <- cfg$stat
    fx$yi <- vapply(seq_len(nrow(fx)), function(i) calc_yi(fx$ct1[[i]], fx$ct2[[i]]), numeric(1))
    fx$vi <- NA_real_
    fx$feature_id <- paste(fx$ct1, fx$ct2, fx$radius_um, fx$stat, sep = "|")
  }

  marks_chr <- as.character(spatstat.geom::marks(pp))
  fx$sample <- sample_id
  fx$tissue_size_um <- cfg$square_size_um
  fx$tissue_area_um2 <- cfg$square_size_um^2
  fx$n_cells_A <- sum(marks_chr == "A")
  fx$n_cells_B <- sum(marks_chr == "B")
  fx$n_cells_total <- fx$n_cells_A + fx$n_cells_B
  fx$parent_lambda <- parent_lambda
  fx$offspring_mean <- offspring_mean
  fx$offspring_sd <- offspring_sd
  fx
}

sim_interleave_extremes <- function(x) {
  x <- as.character(x)
  if (length(x) <= 2L) return(x)
  left <- x[seq(1L, length(x), by = 2L)]
  right <- rev(x[seq(2L, length(x), by = 2L)])
  unique(c(left, right))
}

run_params <- list(
  cache_tag = cfg$cache_tag,
  n_reps = cfg$n_reps,
  n_patients_max = cfg$n_patients_max,
  n_patients_grid = cfg$n_patients_grid,
  n_truth_realizations_per_patient = cfg$n_truth_realizations_per_patient,
  square_size_um = cfg$square_size_um,
  target_cells = cfg$target_cells,
  radius_um = cfg$radius_um,
  stat = cfg$stat,
  boot_nsim = cfg$boot_nsim,
  boot_mode = cfg$boot_mode,
  tile_size = cfg$tile_size,
  sample_param_jitter_sd = cfg$sample_param_jitter_sd,
  seed_base = cfg$seed_base,
  n_workers = cfg$n_workers
)

print(run_params)

bp <- sim_make_bpparam(cfg$n_workers)
bp_started <- FALSE
if (!BiocParallel::bpisup(bp)) {
  bp <- BiocParallel::bpstart(bp)
  bp_started <- TRUE
}

patient_latent <- sim_draw_patient_latents(
  n_patients = cfg$n_patients_max,
  cfg = cfg,
  seed = cfg$seed_base + 13L
)

message("Generating truth realizations ...")
truth_plan <- do.call(rbind, lapply(seq_len(cfg$n_patients_max), function(p_idx) {
  data.frame(
    patient_idx = p_idx,
    truth_idx = seq_len(cfg$n_truth_realizations_per_patient),
    stringsAsFactors = FALSE
  )
}))
truth_plan$seed <- cfg$seed_base + 1000000L + truth_plan$patient_idx * 1000L + truth_plan$truth_idx
truth_plan$patient <- patient_latent$patient[truth_plan$patient_idx]
truth_plan$sample <- sprintf("%s_truth_%03d", truth_plan$patient, truth_plan$truth_idx)

truth_rows <- BiocParallel::bplapply(
  seq_len(nrow(truth_plan)),
  function(j, plan_arg, lat_arg, cfg_arg) {
    row <- plan_arg[j, , drop = FALSE]
    p_idx <- as.integer(row$patient_idx[[1]])
    p_row <- lat_arg[p_idx, , drop = FALSE]
    out <- sim_eval_ab_sample(
      patient_row = p_row,
      sample_id = row$sample[[1]],
      cfg = cfg_arg,
      seed = as.integer(row$seed[[1]]),
      use_bootstrap = FALSE
    )
    out$patient <- row$patient[[1]]
    out$replicate <- NA_integer_
    out
  },
  plan_arg = truth_plan,
  lat_arg = patient_latent,
  cfg_arg = cfg,
  BPPARAM = bp
)
df_truth <- dplyr::bind_rows(truth_rows)

truth_patient <- df_truth %>%
  dplyr::group_by(feature_id, ct1, ct2, radius_um, stat, patient) %>%
  dplyr::summarise(
    true_patient_mean = mean(yi, na.rm = TRUE),
    true_patient_var = stats::var(yi, na.rm = TRUE),
    .groups = "drop"
  )
truth_feature <- truth_patient %>%
  dplyr::group_by(feature_id, ct1, ct2, radius_um, stat) %>%
  dplyr::summarise(
    true_mu = mean(true_patient_mean, na.rm = TRUE),
    true_tau2_patient = stats::var(true_patient_mean, na.rm = TRUE),
    true_tau2_sample = mean(true_patient_var, na.rm = TRUE),
    .groups = "drop"
  )

ord_src <- truth_patient %>%
  dplyr::filter(feature_id == cfg$focus_feature_id) %>%
  dplyr::arrange(dplyr::desc(true_patient_var)) %>%
  dplyr::pull(patient)
if (!length(ord_src)) stop("focus_feature_id not found in truth: ", cfg$focus_feature_id)
patient_order <- sim_interleave_extremes(ord_src)

message("Generating observed one-sample-per-patient data ...")
obs_plan <- do.call(rbind, lapply(seq_len(cfg$n_reps), function(rep_idx) {
  data.frame(
    replicate = rep_idx,
    patient_idx = seq_len(cfg$n_patients_max),
    stringsAsFactors = FALSE
  )
}))
obs_plan$seed <- cfg$seed_base + 2000000L + obs_plan$replicate * 100000L + obs_plan$patient_idx * 1000L
obs_plan$patient <- patient_latent$patient[obs_plan$patient_idx]
obs_plan$sample <- sprintf("%s_rep_%03d", obs_plan$patient, obs_plan$replicate)

obs_rows <- BiocParallel::bplapply(
  seq_len(nrow(obs_plan)),
  function(j, plan_arg, lat_arg, cfg_arg) {
    row <- plan_arg[j, , drop = FALSE]
    p_idx <- as.integer(row$patient_idx[[1]])
    p_row <- lat_arg[p_idx, , drop = FALSE]
    out <- sim_eval_ab_sample(
      patient_row = p_row,
      sample_id = row$sample[[1]],
      cfg = cfg_arg,
      seed = as.integer(row$seed[[1]]),
      use_bootstrap = TRUE
    )
    out$patient <- row$patient[[1]]
    out$replicate <- as.integer(row$replicate[[1]])
    out
  },
  plan_arg = obs_plan,
  lat_arg = patient_latent,
  cfg_arg = cfg,
  BPPARAM = bp
)
df_obs <- dplyr::bind_rows(obs_rows)

message("Fitting meta-analysis over progressive patient subsets ...")
est_rows <- list()
counter <- 1L
total_steps <- length(cfg$n_patients_grid) * cfg$n_reps
pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)

for (rep_idx in seq_len(cfg$n_reps)) {
  rep_obs <- df_obs %>% dplyr::filter(replicate == rep_idx)
  for (n_pat in cfg$n_patients_grid) {
    keep_patients <- patient_order[seq_len(n_pat)]
    df_sub <- rep_obs %>% dplyr::filter(patient %in% keep_patients)

    se_stats <- sim_effects_to_se(df_sub %>%
      dplyr::select(feature_id, ct1, ct2, radius_um, stat, sample, patient, tissue_size_um, tissue_area_um2, yi, vi) %>%
      dplyr::mutate(group = "all"))

    se_meta <- panoramic::panoramic_meta_mv(
      se = se_stats,
      patient_col = "patient",
      group_col = NULL,
      sample_col = "sample",
      tau_structure = "patient",
      method = "REML",
      group_tau2 = "none",
      vi_floor = "median",
      BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
    )

    truth_sub <- truth_patient %>%
      dplyr::filter(patient %in% keep_patients) %>%
      dplyr::group_by(feature_id, ct1, ct2, radius_um, stat) %>%
      dplyr::summarise(
        true_mu = mean(true_patient_mean, na.rm = TRUE),
        true_tau2_patient = stats::var(true_patient_mean, na.rm = TRUE),
        .groups = "drop"
      )

    naive_sub <- df_sub %>%
      dplyr::group_by(feature_id, ct1, ct2, radius_um, stat) %>%
      dplyr::summarise(
        mu_naive = mean(yi, na.rm = TRUE),
        tau2_patient_naive = stats::var(yi, na.rm = TRUE),
        .groups = "drop"
      )

    est <- sim_extract_meta_estimates(se_meta) %>%
      dplyr::left_join(naive_sub %>% dplyr::select(feature_id, mu_naive, tau2_patient_naive), by = "feature_id") %>%
      dplyr::left_join(truth_sub %>% dplyr::select(feature_id, true_mu, true_tau2_patient), by = "feature_id") %>%
      dplyr::mutate(
        n_patients = n_pat,
        replicate = rep_idx
      )

    est_rows[[length(est_rows) + 1L]] <- est
    utils::setTxtProgressBar(pb, counter)
    counter <- counter + 1L
  }
}
close(pb)

if (isTRUE(bp_started) && BiocParallel::bpisup(bp)) {
  try(BiocParallel::bpstop(bp), silent = TRUE)
}

df_est <- dplyr::bind_rows(est_rows)
if (!nrow(df_est)) stop("No estimate rows produced in controlled patient-only simulation.")

params_rds <- file.path(out_data_dir, paste0("intersample_patient_only_params_", cfg$cache_tag, ".rds"))
truth_rds <- file.path(out_data_dir, paste0("intersample_patient_only_truth_", cfg$cache_tag, ".rds"))
obs_rds <- file.path(out_data_dir, paste0("intersample_patient_only_observed_", cfg$cache_tag, ".rds"))
est_rds <- file.path(out_data_dir, paste0("intersample_patient_only_estimates_", cfg$cache_tag, ".rds"))
est_csv <- file.path(out_tbl_dir, paste0("intersample_patient_only_estimates_", cfg$cache_tag, ".csv"))
order_csv <- file.path(out_tbl_dir, paste0("intersample_patient_only_patient_order_", cfg$cache_tag, ".csv"))

sim_safe_save_rds(run_params, params_rds)
sim_safe_save_rds(list(truth_long = df_truth, truth_patient = truth_patient, truth_feature = truth_feature), truth_rds)
sim_safe_save_rds(df_obs, obs_rds)
sim_safe_save_rds(df_est, est_rds)
sim_safe_write_csv(df_est, est_csv)
sim_safe_write_csv(data.frame(order_rank = seq_along(patient_order), patient = patient_order), order_csv)

message("\nSaved controlled patient-only intersample outputs:")
message("  params: ", params_rds)
message("  truth: ", truth_rds)
message("  observed: ", obs_rds)
message("  estimates rds: ", est_rds)
message("  estimates csv: ", est_csv)
message("  patient order: ", order_csv)
