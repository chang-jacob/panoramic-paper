#!/usr/bin/env Rscript

# Simulation Script: Intersample Meta Simulation
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Simulate patient/sample-level heterogeneity with realistic A/B spatial structure.
# - Generate truth and observed effects, then fit PANORAMIC multilevel meta-analysis.
# - Export estimate tables, completion summaries, and simulation metadata.

suppressPackageStartupMessages({
  library(panoramic)
  library(BiocParallel)
  library(dplyr)
  library(spatstat.geom)
})

source("analysis/simulation/sim_utils.R")

project_root <- sim_find_project_root()
paths <- sim_load_paths(project_root)

out_data_dir <- sim_mkdir(file.path(paths$data, "interim", "simulation", "intersample"))
out_tbl_dir <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intersample"))

cfg <- list(
  cache_tag = "intersample_meta",
  n_reps = 20L,
  n_patients_grid = c(10L, 20L, 30L, 40L, 50L),
  n_truth_realizations_per_patient = 80L,
  n_samples_per_patient_mean = 2.5,
  n_samples_per_patient_min = 1L,
  n_samples_per_patient_max = 5L,
  n_samples_dist = "nbinom",
  n_samples_nbinom_size = 1.2,
  square_size_um = 600,
  target_cells = 1900,
  radius_um = 25,
  stat = "local_comp_enrichment",
  boot_nsim = 60L,
  boot_mode = "block",
  tile_size = 62.5,
  patient_log_sd = 0.45,
  sample_log_sd = 0.30,
  base_offspring_mean = 12,
  base_offspring_sd = 55,
  thinning_enabled = TRUE,
  thinning_keep_mean = 0.82,
  thinning_keep_sd = 0.10,
  thinning_min_keep = 0.35,
  thinning_max_keep = 0.98,
  thinning_min_cells_total = 140L,
  thinning_min_cells_per_type = 30L,
  seed_base = 321L,
  n_workers = 10L
)

cfg$cache_tag <- as.character(cfg$cache_tag)
cfg$n_reps <- as.integer(cfg$n_reps)
cfg$n_patients_grid <- sort(unique(as.integer(cfg$n_patients_grid)))
cfg$n_truth_realizations_per_patient <- as.integer(cfg$n_truth_realizations_per_patient)
cfg$n_samples_per_patient_mean <- as.numeric(cfg$n_samples_per_patient_mean)
cfg$n_samples_per_patient_min <- as.integer(cfg$n_samples_per_patient_min)
cfg$n_samples_per_patient_max <- as.integer(cfg$n_samples_per_patient_max)
cfg$n_samples_dist <- tolower(as.character(cfg$n_samples_dist))
cfg$n_samples_nbinom_size <- as.numeric(cfg$n_samples_nbinom_size)
cfg$square_size_um <- as.numeric(cfg$square_size_um)
cfg$target_cells <- as.numeric(cfg$target_cells)
cfg$radius_um <- as.numeric(cfg$radius_um)
cfg$stat <- as.character(cfg$stat)
cfg$boot_nsim <- as.integer(cfg$boot_nsim)
cfg$boot_mode <- as.character(cfg$boot_mode)
cfg$tile_size <- as.numeric(cfg$tile_size)
cfg$patient_log_sd <- as.numeric(cfg$patient_log_sd)
cfg$sample_log_sd <- as.numeric(cfg$sample_log_sd)
cfg$base_offspring_mean <- as.numeric(cfg$base_offspring_mean)
cfg$base_offspring_sd <- as.numeric(cfg$base_offspring_sd)
cfg$thinning_enabled <- isTRUE(cfg$thinning_enabled)
cfg$thinning_keep_mean <- as.numeric(cfg$thinning_keep_mean)
cfg$thinning_keep_sd <- as.numeric(cfg$thinning_keep_sd)
cfg$thinning_min_keep <- as.numeric(cfg$thinning_min_keep)
cfg$thinning_max_keep <- as.numeric(cfg$thinning_max_keep)
cfg$thinning_min_cells_total <- as.integer(cfg$thinning_min_cells_total)
cfg$thinning_min_cells_per_type <- as.integer(cfg$thinning_min_cells_per_type)
cfg$seed_base <- as.integer(cfg$seed_base)
cfg$n_workers <- as.integer(cfg$n_workers)

if (!identical(cfg$stat, "local_comp_enrichment")) {
  stop("This simulation is designed for stat = 'local_comp_enrichment'.")
}
if (cfg$n_reps < 1L) stop("n_reps must be >= 1")
if (length(cfg$n_patients_grid) < 1L || min(cfg$n_patients_grid) < 2L) {
  stop("n_patients_grid must contain integers >= 2")
}
if (cfg$n_truth_realizations_per_patient < 2L) {
  stop("n_truth_realizations_per_patient must be >= 2")
}
if (cfg$n_samples_per_patient_min < 1L) stop("n_samples_per_patient_min must be >= 1")
if (cfg$n_samples_per_patient_max < cfg$n_samples_per_patient_min) stop("n_samples_per_patient_max must be >= n_samples_per_patient_min")
if (cfg$n_samples_per_patient_mean < cfg$n_samples_per_patient_min) stop("n_samples_per_patient_mean must be >= n_samples_per_patient_min")
if (!cfg$n_samples_dist %in% c("poisson", "nbinom")) stop("n_samples_dist must be 'poisson' or 'nbinom'")
if (!is.finite(cfg$n_samples_nbinom_size) || cfg$n_samples_nbinom_size <= 0) stop("n_samples_nbinom_size must be > 0")
if (cfg$n_workers < 1L) stop("n_workers must be >= 1")
if (!is.finite(cfg$thinning_min_keep) || cfg$thinning_min_keep <= 0 || cfg$thinning_min_keep >= 1) stop("thinning_min_keep must be in (0, 1)")
if (!is.finite(cfg$thinning_max_keep) || cfg$thinning_max_keep <= cfg$thinning_min_keep || cfg$thinning_max_keep > 1) {
  stop("thinning_max_keep must be in (thinning_min_keep, 1]")
}

cfg$tissue_area_um2 <- cfg$square_size_um^2
cfg$target_cell_density <- cfg$target_cells / cfg$tissue_area_um2
cfg$base_parent_lambda <- cfg$target_cell_density / (1 + cfg$base_offspring_mean)

sim_draw_patient_latents <- function(n_patients, cfg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  patient_ids <- sprintf("patient_%03d", seq_len(n_patients))

  strength <- stats::rnorm(n_patients, mean = 0, sd = cfg$patient_log_sd)

  parent_lambda_patient <- pmax(
    cfg$target_cell_density / 400,
    cfg$base_parent_lambda * exp(strength)
  )
  offspring_mean_patient <- pmax(
    1,
    cfg$base_offspring_mean * exp(0.7 * strength)
  )
  offspring_sd_patient <- pmax(
    8,
    cfg$base_offspring_sd * exp(-0.25 * strength)
  )

  data.frame(
    patient = patient_ids,
    strength = strength,
    parent_lambda_patient = parent_lambda_patient,
    offspring_mean_patient = offspring_mean_patient,
    offspring_sd_patient = offspring_sd_patient,
    stringsAsFactors = FALSE
  )
}

sim_apply_global_thinning <- function(pp, cfg) {
  if (!isTRUE(cfg$thinning_enabled) || pp$n < 1L) return(pp)

  df_pp <- as.data.frame(pp)
  marks_chr <- as.character(df_pp$marks)
  if (length(marks_chr) < 1L) return(pp)

  keep_prob_base <- sim_clamp(
    stats::rnorm(1, mean = cfg$thinning_keep_mean, sd = cfg$thinning_keep_sd),
    cfg$thinning_min_keep,
    cfg$thinning_max_keep
  )

  keep <- rep(TRUE, length(marks_chr))
  for (attempt in seq_len(4L)) {
    keep_prob <- sim_clamp(
      keep_prob_base + 0.08 * (attempt - 1L),
      cfg$thinning_min_keep,
      cfg$thinning_max_keep
    )
    keep <- stats::runif(length(marks_chr)) < keep_prob

    n_total <- sum(keep)
    n_a <- sum(keep & marks_chr == "A")
    n_b <- sum(keep & marks_chr == "B")
    if (n_total >= cfg$thinning_min_cells_total &&
        n_a >= cfg$thinning_min_cells_per_type &&
        n_b >= cfg$thinning_min_cells_per_type) {
      break
    }
  }

  if (!any(keep) || sum(keep) < cfg$thinning_min_cells_total) return(pp)

  spatstat.geom::ppp(
    x = df_pp$x[keep],
    y = df_pp$y[keep],
    window = spatstat.geom::square(cfg$square_size_um),
    marks = factor(marks_chr[keep], levels = c("A", "B"))
  )
}

sim_eval_ab_sample <- function(patient_row, patient_id, sample_id, cfg, seed, use_bootstrap = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  parent_lambda <- pmax(
    cfg$target_cell_density / 400,
    patient_row$parent_lambda_patient * exp(stats::rnorm(1, 0, cfg$sample_log_sd))
  )
  offspring_mean <- pmax(
    1,
    patient_row$offspring_mean_patient * exp(stats::rnorm(1, 0, cfg$sample_log_sd))
  )
  offspring_sd <- pmax(
    8,
    patient_row$offspring_sd_patient * exp(stats::rnorm(1, 0, cfg$sample_log_sd / 2))
  )

  pp <- sim_simulate_ab_cluster_ppp(
    square_size = cfg$square_size_um,
    parent_lambda = parent_lambda,
    offspring_mean = offspring_mean,
    offspring_sd = offspring_sd
  )
  pp <- sim_apply_global_thinning(pp, cfg)

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
    fx <- sim_flatten_spatialstats(se_stats) %>%
      dplyr::filter(
        ct1 %in% c("A", "B"),
        ct2 %in% c("A", "B"),
        radius_um == cfg$radius_um,
        stat == cfg$stat
      )
    if (!nrow(fx)) stop("No A/B rows returned by panoramic_spatialstats.")
  } else {
    precompute_local_comp_cache <- getFromNamespace(".precompute_local_comp_cache", "panoramic")
    local_comp_score_fn <- getFromNamespace(".local_comp_score", "panoramic")

    cache <- precompute_local_comp_cache(pp, cfg$radius_um)
    ridx <- match(cfg$radius_um, cache$radii_um)
    if (is.na(ridx)) {
      tol <- sqrt(.Machine$double.eps) * max(1, abs(cfg$radius_um))
      idx <- which(abs(cache$radii_um - cfg$radius_um) <= tol)
      ridx <- if (length(idx) > 0L) idx[1] else NA_integer_
    }
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
      if (identical(ct1, ct2) && is.finite(cache$area_window) && cache$area_window > 0) {
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
  n_a <- sum(marks_chr == "A")
  n_b <- sum(marks_chr == "B")

  fx$sample <- sample_id
  fx$patient <- patient_id
  fx$group <- "all"
  fx$tissue_size_um <- cfg$square_size_um
  fx$tissue_area_um2 <- cfg$tissue_area_um2
  fx$n_cells_A <- n_a
  fx$n_cells_B <- n_b
  fx$n_cells_total <- n_a + n_b
  fx$parent_lambda <- parent_lambda
  fx$offspring_mean <- offspring_mean
  fx$offspring_sd <- offspring_sd
  fx
}

run_params <- list(
  cache_tag = cfg$cache_tag,
  n_reps = cfg$n_reps,
  n_patients_grid = cfg$n_patients_grid,
  n_truth_realizations_per_patient = cfg$n_truth_realizations_per_patient,
  n_samples_per_patient_mean = cfg$n_samples_per_patient_mean,
  n_samples_per_patient_min = cfg$n_samples_per_patient_min,
  n_samples_per_patient_max = cfg$n_samples_per_patient_max,
  n_samples_dist = cfg$n_samples_dist,
  n_samples_nbinom_size = cfg$n_samples_nbinom_size,
  square_size_um = cfg$square_size_um,
  tissue_area_um2 = cfg$tissue_area_um2,
  target_cells = cfg$target_cells,
  target_cell_density = cfg$target_cell_density,
  radius_um = cfg$radius_um,
  stat = cfg$stat,
  boot_nsim = cfg$boot_nsim,
  boot_mode = cfg$boot_mode,
  tile_size = cfg$tile_size,
  patient_log_sd = cfg$patient_log_sd,
  sample_log_sd = cfg$sample_log_sd,
  thinning_enabled = cfg$thinning_enabled,
  thinning_keep_mean = cfg$thinning_keep_mean,
  thinning_keep_sd = cfg$thinning_keep_sd,
  seed_base = cfg$seed_base,
  n_workers = cfg$n_workers
)

print(run_params)

max_n_patients <- max(cfg$n_patients_grid)
patient_latent <- sim_draw_patient_latents(max_n_patients, cfg, seed = cfg$seed_base + 13L)

bp <- sim_make_bpparam(cfg$n_workers)
bp_started <- FALSE
if (!BiocParallel::bpisup(bp)) {
  bp <- BiocParallel::bpstart(bp)
  bp_started <- TRUE
}

message("Generating truth references for each patient ...")
truth_plan <- do.call(rbind, lapply(seq_len(max_n_patients), function(p_idx) {
  data.frame(
    patient_idx = p_idx,
    truth_idx = seq_len(cfg$n_truth_realizations_per_patient),
    stringsAsFactors = FALSE
  )
}))
truth_plan$seed <- cfg$seed_base + 1000000L + truth_plan$patient_idx * 1000L + truth_plan$truth_idx
truth_plan$patient <- patient_latent$patient[truth_plan$patient_idx]
truth_plan$sample <- sprintf("%s_truth_%03d", truth_plan$patient, truth_plan$truth_idx)

truth_res <- BiocParallel::bplapply(
  seq_len(nrow(truth_plan)),
  function(j, plan_arg, lat_arg, cfg_arg) {
    row <- plan_arg[j, , drop = FALSE]
    p_idx <- as.integer(row$patient_idx[[1]])
    p_row <- lat_arg[p_idx, , drop = FALSE]
    out <- tryCatch({
      sim_eval_ab_sample(
        patient_row = p_row,
        patient_id = row$patient[[1]],
        sample_id = row$sample[[1]],
        cfg = cfg_arg,
        seed = as.integer(row$seed[[1]]),
        use_bootstrap = FALSE
      )
    }, error = function(e) e)

    if (inherits(out, "error")) {
      return(list(ok = FALSE, error = conditionMessage(out), patient = row$patient[[1]], sample = row$sample[[1]]))
    }
    list(ok = TRUE, data = out)
  },
  plan_arg = truth_plan,
  lat_arg = patient_latent,
  cfg_arg = cfg,
  BPPARAM = bp
)

truth_ok <- vapply(truth_res, function(x) is.list(x) && isTRUE(x$ok), logical(1))
if (!any(truth_ok)) {
  stop("No truth rows were generated. Inspect simulation settings and cell-count thresholds.")
}

df_truth <- dplyr::bind_rows(lapply(truth_res[truth_ok], `[[`, "data"))
truth_fail <- dplyr::bind_rows(lapply(truth_res[!truth_ok], function(x) {
  data.frame(
    stage = "truth",
    patient = if (!is.null(x$patient)) x$patient else NA_character_,
    sample = if (!is.null(x$sample)) x$sample else NA_character_,
    error = if (!is.null(x$error)) x$error else "unknown truth error",
    stringsAsFactors = FALSE
  )
}))

truth_patient <- df_truth %>%
  dplyr::group_by(feature_id, ct1, ct2, radius_um, stat, patient) %>%
  dplyr::summarise(
    true_patient_mean = mean(yi, na.rm = TRUE),
    true_patient_var = stats::var(yi, na.rm = TRUE),
    .groups = "drop"
  )

message("Generating observed samples and fitting meta-analysis ...")

est_rows <- list()
fail_rows <- list()
if (nrow(truth_fail) > 0L) fail_rows[[length(fail_rows) + 1L]] <- truth_fail
sample_count_rows <- list()

total_steps <- cfg$n_reps * length(cfg$n_patients_grid)
progress_counter <- 0L
pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)

for (rep_idx in seq_len(cfg$n_reps)) {
  n_samples_per_patient <- vapply(
    seq_len(max_n_patients),
    function(p_idx) {
      sim_draw_n_samples_per_patient(
        cfg,
        seed = cfg$seed_base + rep_idx * 100000L + p_idx * 1000L + 17L
      )
    },
    integer(1)
  )

  obs_plan <- do.call(rbind, lapply(seq_len(max_n_patients), function(p_idx) {
    data.frame(
      patient_idx = p_idx,
      sample_ord = seq_len(n_samples_per_patient[[p_idx]]),
      stringsAsFactors = FALSE
    )
  }))
  obs_plan$seed <- cfg$seed_base + 2000000L + rep_idx * 100000L + obs_plan$patient_idx * 1000L + obs_plan$sample_ord
  obs_plan$patient <- patient_latent$patient[obs_plan$patient_idx]
  obs_plan$sample <- sprintf("%s_rep_%03d_s%02d", obs_plan$patient, rep_idx, obs_plan$sample_ord)

  obs_res <- BiocParallel::bplapply(
    seq_len(nrow(obs_plan)),
    function(j, plan_arg, lat_arg, cfg_arg) {
      row <- plan_arg[j, , drop = FALSE]
      p_idx <- as.integer(row$patient_idx[[1]])
      p_row <- lat_arg[p_idx, , drop = FALSE]
      out <- tryCatch({
        sim_eval_ab_sample(
          patient_row = p_row,
          patient_id = row$patient[[1]],
          sample_id = row$sample[[1]],
          cfg = cfg_arg,
          seed = as.integer(row$seed[[1]]),
          use_bootstrap = TRUE
        )
      }, error = function(e) e)

      if (inherits(out, "error")) {
        return(list(ok = FALSE, error = conditionMessage(out), patient = row$patient[[1]], sample = row$sample[[1]]))
      }
      list(ok = TRUE, data = out)
    },
    plan_arg = obs_plan,
    lat_arg = patient_latent,
    cfg_arg = cfg,
    BPPARAM = bp
  )

  obs_ok <- vapply(obs_res, function(x) is.list(x) && isTRUE(x$ok), logical(1))
  if (!any(obs_ok)) {
    fail_rows[[length(fail_rows) + 1L]] <- data.frame(
      stage = "observed",
      replicate = rep_idx,
      n_patients = NA_integer_,
      patient = NA_character_,
      sample = NA_character_,
      error = "All observed sample jobs failed for this replicate",
      stringsAsFactors = FALSE
    )
    for (n_pat in cfg$n_patients_grid) {
      progress_counter <- progress_counter + 1L
      utils::setTxtProgressBar(pb, progress_counter)
    }
    next
  }

  df_obs_rep <- dplyr::bind_rows(lapply(obs_res[obs_ok], `[[`, "data"))
  obs_fail <- dplyr::bind_rows(lapply(obs_res[!obs_ok], function(x) {
    data.frame(
      stage = "observed",
      replicate = rep_idx,
      n_patients = NA_integer_,
      patient = if (!is.null(x$patient)) x$patient else NA_character_,
      sample = if (!is.null(x$sample)) x$sample else NA_character_,
      error = if (!is.null(x$error)) x$error else "unknown observed error",
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(obs_fail) > 0L) fail_rows[[length(fail_rows) + 1L]] <- obs_fail

  for (n_pat in cfg$n_patients_grid) {
    keep_patients <- patient_latent$patient[seq_len(n_pat)]

    df_sub <- df_obs_rep %>%
      dplyr::filter(patient %in% keep_patients)

    have_patients <- sort(unique(df_sub$patient))
    if (length(have_patients) < n_pat) {
      miss_pat <- setdiff(keep_patients, have_patients)
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        stage = "subset",
        replicate = rep_idx,
        n_patients = n_pat,
        patient = paste(miss_pat, collapse = ","),
        sample = NA_character_,
        error = "Missing observed rows for one or more required patients",
        stringsAsFactors = FALSE
      )
      progress_counter <- progress_counter + 1L
      utils::setTxtProgressBar(pb, progress_counter)
      next
    }

    sample_counts <- df_sub %>%
      dplyr::group_by(patient) %>%
      dplyr::summarise(n_samples = dplyr::n_distinct(sample), .groups = "drop") %>%
      dplyr::mutate(n_patients = n_pat, replicate = rep_idx)
    sample_count_rows[[length(sample_count_rows) + 1L]] <- sample_counts

    truth_sub <- truth_patient %>%
      dplyr::filter(patient %in% keep_patients) %>%
      dplyr::group_by(feature_id, ct1, ct2, radius_um, stat) %>%
      dplyr::summarise(
        true_mu = mean(true_patient_mean, na.rm = TRUE),
        true_tau2_patient = stats::var(true_patient_mean, na.rm = TRUE),
        true_tau2_sample = mean(true_patient_var, na.rm = TRUE),
        .groups = "drop"
      )

    naive_sub <- df_sub %>%
      dplyr::group_by(feature_id, ct1, ct2, radius_um, stat) %>%
      dplyr::summarise(
        mu_naive = {
          pm <- tapply(yi, patient, function(z) mean(z, na.rm = TRUE))
          mean(as.numeric(pm), na.rm = TRUE)
        },
        tau2_overall_naive = stats::var(yi, na.rm = TRUE),
        tau2_patient_naive = {
          pm <- tapply(yi, patient, function(z) mean(z, na.rm = TRUE))
          stats::var(as.numeric(pm), na.rm = TRUE)
        },
        tau2_sample_naive = {
          pv <- tapply(yi, patient, function(z) stats::var(z, na.rm = TRUE))
          mean(as.numeric(pv), na.rm = TRUE)
        },
        .groups = "drop"
      )

    precision_sub <- df_sub %>%
      dplyr::group_by(feature_id) %>%
      dplyr::summarise(
        n_samples_obs = dplyr::n_distinct(sample),
        n_patients_obs = dplyr::n_distinct(patient),
        mean_vi = mean(vi, na.rm = TRUE),
        mean_n_cells_total = mean(n_cells_total, na.rm = TRUE),
        .groups = "drop"
      )

    se_stats <- sim_effects_to_se(df_sub %>%
      dplyr::select(
        feature_id, ct1, ct2, radius_um, stat,
        sample, patient, group, tissue_size_um, tissue_area_um2, yi, vi
      ))

    se_meta <- tryCatch({
      panoramic::panoramic_meta_mv(
        se = se_stats,
        patient_col = "patient",
        group_col = NULL,
        sample_col = "sample",
        tau_structure = "patient_sample",
        method = "REML",
        group_tau2 = "none",
        vi_floor = "median",
        BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
      )
    }, error = function(e) e)

    if (inherits(se_meta, "error")) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        stage = "meta",
        replicate = rep_idx,
        n_patients = n_pat,
        patient = NA_character_,
        sample = NA_character_,
        error = conditionMessage(se_meta),
        stringsAsFactors = FALSE
      )
      progress_counter <- progress_counter + 1L
      utils::setTxtProgressBar(pb, progress_counter)
      next
    }

    est <- sim_extract_meta_estimates(se_meta) %>%
      dplyr::left_join(
        naive_sub %>% dplyr::select(feature_id, mu_naive, tau2_overall_naive, tau2_patient_naive, tau2_sample_naive),
        by = "feature_id"
      ) %>%
      dplyr::left_join(
        precision_sub,
        by = "feature_id"
      ) %>%
      dplyr::left_join(
        truth_sub %>% dplyr::select(feature_id, true_mu, true_tau2_patient, true_tau2_sample),
        by = "feature_id"
      ) %>%
      dplyr::mutate(
        true_tau2 = true_tau2_patient,
        n_patients = n_pat,
        replicate = rep_idx
      )

    if (!nrow(est)) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        stage = "meta_extract",
        replicate = rep_idx,
        n_patients = n_pat,
        patient = NA_character_,
        sample = NA_character_,
        error = "No rows returned by sim_extract_meta_estimates",
        stringsAsFactors = FALSE
      )
    } else {
      est_rows[[length(est_rows) + 1L]] <- est
    }

    progress_counter <- progress_counter + 1L
    utils::setTxtProgressBar(pb, progress_counter)
  }
}
close(pb)

if (isTRUE(bp_started) && BiocParallel::bpisup(bp)) {
  try(BiocParallel::bpstop(bp), silent = TRUE)
}

df_est <- dplyr::bind_rows(est_rows)
df_fail <- dplyr::bind_rows(fail_rows)
df_sample_counts <- dplyr::bind_rows(sample_count_rows)

if (!nrow(df_est)) {
  fail_csv <- file.path(out_tbl_dir, "intersample_meta_failures.csv")
  if (nrow(df_fail) > 0L) sim_safe_write_csv(df_fail, fail_csv)
  stop("No estimate rows were produced. See failures at: ", fail_csv)
}

completion_tbl <- data.frame(
  n_patients = cfg$n_patients_grid,
  n_reps_done = vapply(
    cfg$n_patients_grid,
    function(np) length(unique(df_est$replicate[df_est$n_patients == np])),
    integer(1)
  ),
  n_reps_target = cfg$n_reps,
  stringsAsFactors = FALSE
)

params_rds <- file.path(out_data_dir, "intersample_meta_params.rds")
truth_rds <- file.path(out_data_dir, "intersample_meta_truth.rds")
est_rds <- file.path(out_data_dir, "intersample_meta_estimates.rds")
est_csv <- file.path(out_tbl_dir, "intersample_meta_estimates.csv")
fail_csv <- file.path(out_tbl_dir, "intersample_meta_failures.csv")
sample_count_csv <- file.path(out_tbl_dir, "intersample_samples_per_patient.csv")
completion_csv <- file.path(out_tbl_dir, "intersample_completion_status.csv")

sim_safe_save_rds(run_params, params_rds)
sim_safe_save_rds(
  list(truth_long = df_truth, truth_patient = truth_patient),
  truth_rds
)
sim_safe_save_rds(df_est, est_rds)
sim_safe_write_csv(df_est, est_csv)
if (nrow(df_fail) > 0L) sim_safe_write_csv(df_fail, fail_csv)
if (nrow(df_sample_counts) > 0L) sim_safe_write_csv(df_sample_counts, sample_count_csv)
sim_safe_write_csv(completion_tbl, completion_csv)

message("\nSaved intersample simulation outputs:")
message("  params: ", params_rds)
message("  truth: ", truth_rds)
message("  estimates rds: ", est_rds)
message("  estimates csv: ", est_csv)
if (nrow(df_sample_counts) > 0L) message("  sample-count table: ", sample_count_csv)
message("  completion status: ", completion_csv)
if (nrow(df_fail) > 0L) message("  failures: ", fail_csv)
if (any(completion_tbl$n_reps_done < completion_tbl$n_reps_target)) {
  message("  WARNING: not all n_patients settings reached target replicates.")
}
