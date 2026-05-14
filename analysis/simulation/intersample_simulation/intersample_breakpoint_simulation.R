#!/usr/bin/env Rscript

# Simulation Script: Intersample Breakpoint Simulation
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Simulate contamination/quality-breakpoint scenarios across patient-count regimes.
# - Fit PANORAMIC meta-analysis and quantify regime-dependent performance degradation.
# - Export breakpoint diagnostics, feature/pair summaries, and manuscript figures.

suppressPackageStartupMessages({
  library(panoramic)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(spatstat.geom)
})

source("analysis/simulation/sim_utils.R")

project_root <- sim_find_project_root()
paths <- sim_load_paths(project_root)

out_data_dir <- sim_mkdir(file.path(paths$data, "interim", "simulation", "intersample_breakpoint"))
out_tbl_dir <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intersample_breakpoint"))
out_fig_dir <- sim_mkdir(file.path(paths$output, "figures", "simulation", "intersample_breakpoint"))

sim_stamp <- function() format(Sys.time(), "%H:%M:%S")
sim_elapsed <- function(t0) sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
sim_log <- function(...) message(sprintf("[%s] %s", sim_stamp(), paste0(..., collapse = "")))

cfg <- list(
  cache_tag = "intersample_breakpoint",
  # Manuscript-quality defaults for the sample degradation study.
  n_reps = 5L,
  n_patients_grid = c(10L, 20L, 30L, 40L, 50L),
  n_truth_realizations_per_patient = 100L,
  n_samples_per_patient_mean = 2.2,
  n_samples_per_patient_min = 1L,
  n_samples_per_patient_max = 4L,
  n_samples_dist = "nbinom",
  n_samples_nbinom_size = 1.1,
  square_size_um = 600,
  target_cells = 1900,
  radius_um = 25,
  stat = "local_comp_enrichment",
  boot_nsim = 100L,
  boot_mode = "block",
  tile_size = 62.5,
  patient_log_sd = 0.35,
  sample_log_sd_base = 0.20,
  base_offspring_mean = 12,
  base_offspring_sd = 55,
  min_cells_total = 130L,
  min_cells_per_type = 28L,
  min_keep = 0.10,
  max_keep = 0.99,
  max_thinning_attempts = 5L,
  max_sample_attempts = 4L,
  parallel_backend = "multicore",
  seed_base = 811L,
  n_workers = 10L
)

if (exists("sim_cfg", inherits = TRUE)) {
  cfg <- utils::modifyList(cfg, get("sim_cfg", inherits = TRUE))
}

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
cfg$sample_log_sd_base <- as.numeric(cfg$sample_log_sd_base)
cfg$base_offspring_mean <- as.numeric(cfg$base_offspring_mean)
cfg$base_offspring_sd <- as.numeric(cfg$base_offspring_sd)
cfg$min_cells_total <- as.integer(cfg$min_cells_total)
cfg$min_cells_per_type <- as.integer(cfg$min_cells_per_type)
cfg$min_keep <- as.numeric(cfg$min_keep)
cfg$max_keep <- as.numeric(cfg$max_keep)
cfg$max_thinning_attempts <- as.integer(cfg$max_thinning_attempts)
cfg$max_sample_attempts <- as.integer(cfg$max_sample_attempts)
cfg$parallel_backend <- tolower(as.character(cfg$parallel_backend))
cfg$seed_base <- as.integer(cfg$seed_base)
cfg$n_workers <- as.integer(cfg$n_workers)

if (!identical(cfg$stat, "local_comp_enrichment")) {
  stop("This breakpoint simulation is designed for stat='local_comp_enrichment'.")
}
if (cfg$n_workers < 1L) stop("n_workers must be >= 1")
if (cfg$n_reps < 1L) stop("n_reps must be >= 1")
if (length(cfg$n_patients_grid) < 1L || min(cfg$n_patients_grid) < 2L) stop("n_patients_grid must be >= 2")
if (cfg$n_truth_realizations_per_patient < 3L) stop("n_truth_realizations_per_patient must be >= 3")
if (!cfg$n_samples_dist %in% c("poisson", "nbinom")) stop("n_samples_dist must be 'poisson' or 'nbinom'")
if (!is.finite(cfg$n_samples_nbinom_size) || cfg$n_samples_nbinom_size <= 0) stop("n_samples_nbinom_size must be > 0")
if (cfg$n_samples_per_patient_min < 1L) stop("n_samples_per_patient_min must be >= 1")
if (cfg$n_samples_per_patient_max < cfg$n_samples_per_patient_min) stop("n_samples_per_patient_max must be >= n_samples_per_patient_min")
if (cfg$n_samples_per_patient_mean < cfg$n_samples_per_patient_min) stop("n_samples_per_patient_mean must be >= n_samples_per_patient_min")
if (!cfg$parallel_backend %in% c("multicore", "snow", "serial")) {
  stop("parallel_backend must be one of: multicore, snow, serial")
}

cfg$tissue_area_um2 <- cfg$square_size_um^2
cfg$target_cell_density <- cfg$target_cells / cfg$tissue_area_um2
cfg$base_parent_lambda <- cfg$target_cell_density / (1 + cfg$base_offspring_mean)

scenario_tbl <- data.frame(
  scenario = c("baseline", "mild_break", "moderate_break", "severe_break", "extreme_break"),
  contamination_frac = c(0.00, 0.20, 0.35, 0.50, 0.65),
  quality_sd_clean = c(0.04, 0.04, 0.04, 0.04, 0.04),
  quality_sd_bad = c(0.05, 0.20, 0.30, 0.45, 0.60),
  quality_coupling = c(0.00, 0.50, 0.70, 0.90, 1.10),
  quality_to_anchor_keep = c(0.00, 0.25, 0.40, 0.58, 0.75),
  quality_to_nonanchor_keep = c(0.00, 0.06, 0.12, 0.20, 0.30),
  quality_to_tile = c(0.00, 0.18, 0.28, 0.42, 0.60),
  directional_bias_target = c("A", "A", "A", "A", "A"),
  directional_bias_strength = c(0.00, 0.18, 0.32, 0.48, 0.65),
  sample_log_sd_extra_bad = c(0.00, 0.10, 0.18, 0.28, 0.40),
  base_keep_mean = c(0.90, 0.87, 0.85, 0.83, 0.80),
  base_keep_sd = c(0.02, 0.05, 0.07, 0.10, 0.14),
  stringsAsFactors = FALSE
)

need_scn <- c(
  "scenario", "contamination_frac", "quality_sd_clean", "quality_sd_bad", "quality_coupling",
  "quality_to_anchor_keep", "quality_to_nonanchor_keep", "quality_to_tile",
  "directional_bias_target", "directional_bias_strength",
  "sample_log_sd_extra_bad",
  "base_keep_mean", "base_keep_sd"
)
miss_scn <- setdiff(need_scn, colnames(scenario_tbl))
if (length(miss_scn) > 0L) stop("scenario_tbl missing required columns: ", paste(miss_scn, collapse = ", "))

scenario_tbl <- scenario_tbl %>%
  mutate(
    scenario = as.character(scenario),
    contamination_frac = as.numeric(contamination_frac),
    quality_sd_clean = as.numeric(quality_sd_clean),
    quality_sd_bad = as.numeric(quality_sd_bad),
    quality_coupling = as.numeric(quality_coupling),
    quality_to_anchor_keep = as.numeric(quality_to_anchor_keep),
    quality_to_nonanchor_keep = as.numeric(quality_to_nonanchor_keep),
    quality_to_tile = as.numeric(quality_to_tile),
    directional_bias_target = toupper(as.character(directional_bias_target)),
    directional_bias_strength = as.numeric(directional_bias_strength),
    sample_log_sd_extra_bad = as.numeric(sample_log_sd_extra_bad),
    base_keep_mean = as.numeric(base_keep_mean),
    base_keep_sd = as.numeric(base_keep_sd)
  )
scenario_tbl$directional_bias_target[!scenario_tbl$directional_bias_target %in% c("A", "B")] <- "A"
scenario_tbl$directional_bias_strength[!is.finite(scenario_tbl$directional_bias_strength)] <- 0
scenario_tbl$directional_bias_strength <- pmax(0, scenario_tbl$directional_bias_strength)

sim_draw_patient_latents <- function(n_patients, cfg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  strength <- stats::rnorm(n_patients, mean = 0, sd = cfg$patient_log_sd)
  patient_ids <- sprintf("patient_%03d", seq_len(n_patients))
  data.frame(
    patient = patient_ids,
    strength = strength,
    parent_lambda_patient = pmax(cfg$target_cell_density / 450, cfg$base_parent_lambda * exp(strength)),
    offspring_mean_patient = pmax(1, cfg$base_offspring_mean * exp(0.7 * strength)),
    offspring_sd_patient = pmax(8, cfg$base_offspring_sd * exp(-0.25 * strength)),
    stringsAsFactors = FALSE
  )
}

sim_apply_quality_degradation <- function(pp, cfg, scenario_row, quality_abs) {
  if (pp$n < 1L) return(pp)
  df_pp <- as.data.frame(pp)
  marks_chr <- as.character(df_pp$marks)

  keep_mean <- sim_clamp(
    scenario_row$base_keep_mean - scenario_row$quality_to_anchor_keep * 0.10 * quality_abs,
    cfg$min_keep,
    cfg$max_keep
  )
  keep_base <- sim_clamp(
    stats::rnorm(1, mean = keep_mean, sd = scenario_row$base_keep_sd),
    cfg$min_keep,
    cfg$max_keep
  )

  keep_a <- sim_clamp(
    keep_base - scenario_row$quality_to_anchor_keep * quality_abs,
    cfg$min_keep,
    cfg$max_keep
  )
  keep_b <- sim_clamp(
    keep_base - scenario_row$quality_to_nonanchor_keep * quality_abs,
    cfg$min_keep,
    cfg$max_keep
  )
  bias_target <- toupper(as.character(scenario_row$directional_bias_target))
  if (!bias_target %in% c("A", "B")) bias_target <- "A"
  bias_strength <- as.numeric(scenario_row$directional_bias_strength)
  if (!is.finite(bias_strength) || bias_strength < 0) bias_strength <- 0
  bias_delta <- sim_clamp(bias_strength * quality_abs, 0, 0.55)
  if (bias_delta > 0) {
    if (identical(bias_target, "A")) {
      keep_a <- sim_clamp(keep_a + 0.55 * bias_delta, cfg$min_keep, cfg$max_keep)
      keep_b <- sim_clamp(keep_b - 1.20 * bias_delta, cfg$min_keep, cfg$max_keep)
    } else {
      keep_b <- sim_clamp(keep_b + 0.55 * bias_delta, cfg$min_keep, cfg$max_keep)
      keep_a <- sim_clamp(keep_a - 1.20 * bias_delta, cfg$min_keep, cfg$max_keep)
    }
  }

  keep <- rep(TRUE, length(marks_chr))
  for (attempt in seq_len(cfg$max_thinning_attempts)) {
    keep_prob_a <- sim_clamp(keep_a + 0.08 * (attempt - 1L), cfg$min_keep, cfg$max_keep)
    keep_prob_b <- sim_clamp(keep_b + 0.08 * (attempt - 1L), cfg$min_keep, cfg$max_keep)
    keep_prob <- ifelse(marks_chr == "A", keep_prob_a, keep_prob_b)
    keep <- stats::runif(length(marks_chr)) < keep_prob

    n_total <- sum(keep)
    n_a <- sum(keep & marks_chr == "A")
    n_b <- sum(keep & marks_chr == "B")
    if (n_total >= cfg$min_cells_total && n_a >= cfg$min_cells_per_type && n_b >= cfg$min_cells_per_type) break
  }

  if (!any(keep)) return(pp)
  marks_keep <- marks_chr[keep]
  if (bias_delta > 0) {
    flip_frac <- sim_clamp(0.90 * bias_delta, 0, 0.65)
    source_type <- if (identical(bias_target, "A")) "B" else "A"
    source_idx <- which(marks_keep == source_type)
    max_flip <- max(0L, length(source_idx) - cfg$min_cells_per_type)
    n_flip <- min(as.integer(floor(length(source_idx) * flip_frac)), max_flip)
    if (n_flip > 0L) {
      flip_idx <- sample(source_idx, size = n_flip, replace = FALSE)
      marks_keep[flip_idx] <- bias_target
    }
  }

  out <- spatstat.geom::ppp(
    x = df_pp$x[keep],
    y = df_pp$y[keep],
    window = spatstat.geom::square(cfg$square_size_um),
    marks = factor(marks_keep, levels = c("A", "B"))
  )
  out
}

sim_calc_truth_yi_ab <- function(pp, radius_um, stat) {
  if (!identical(stat, "local_comp_enrichment")) {
    stop("sim_calc_truth_yi_ab only supports local_comp_enrichment")
  }
  precompute_local_comp_cache <- getFromNamespace(".precompute_local_comp_cache", "panoramic")
  local_comp_score_fn <- getFromNamespace(".local_comp_score", "panoramic")
  cache <- precompute_local_comp_cache(pp, radius_um)

  ridx <- match(radius_um, cache$radii_um)
  if (is.na(ridx)) {
    tol <- sqrt(.Machine$double.eps) * max(1, abs(radius_um))
    idx <- which(abs(cache$radii_um - radius_um) <= tol)
    ridx <- if (length(idx) > 0L) idx[1] else NA_integer_
  }
  if (is.na(ridx)) stop("Could not match radius in local_comp cache")

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

  out <- pairs
  out$radius_um <- radius_um
  out$stat <- stat
  out$yi <- vapply(seq_len(nrow(out)), function(i) calc_yi(out$ct1[[i]], out$ct2[[i]]), numeric(1))
  out$vi <- NA_real_
  out$feature_id <- paste(out$ct1, out$ct2, out$radius_um, out$stat, sep = "|")
  out
}

sim_eval_one_sample <- function(patient_row, sample_id, cfg, seed, scenario_row = NULL, contaminated = FALSE, use_bootstrap = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(scenario_row)) {
    scenario_row <- as.list(data.frame(
      quality_sd_clean = 0,
      quality_sd_bad = 0,
      quality_coupling = 0,
      quality_to_anchor_keep = 0,
      quality_to_nonanchor_keep = 0,
      quality_to_tile = 0,
      directional_bias_target = "A",
      directional_bias_strength = 0,
      sample_log_sd_extra_bad = 0,
      base_keep_mean = 0.95,
      base_keep_sd = 0.02
    ))
  }

  quality_sd <- if (isTRUE(contaminated)) scenario_row$quality_sd_bad else scenario_row$quality_sd_clean
  quality_z <- scenario_row$quality_coupling * patient_row$strength + stats::rnorm(1, 0, quality_sd)
  quality_abs <- abs(quality_z)

  sample_log_sd <- cfg$sample_log_sd_base + if (isTRUE(contaminated)) scenario_row$sample_log_sd_extra_bad else 0

  parent_lambda <- pmax(
    cfg$target_cell_density / 450,
    patient_row$parent_lambda_patient * exp(stats::rnorm(1, 0, sample_log_sd))
  )
  offspring_mean <- pmax(
    1,
    patient_row$offspring_mean_patient * exp(stats::rnorm(1, 0, sample_log_sd))
  )
  offspring_sd <- pmax(
    8,
    patient_row$offspring_sd_patient * exp(stats::rnorm(1, 0, sample_log_sd / 2))
  )

  pp <- sim_simulate_ab_cluster_ppp(
    square_size = cfg$square_size_um,
    parent_lambda = parent_lambda,
    offspring_mean = offspring_mean,
    offspring_sd = offspring_sd
  )

  if (isTRUE(contaminated)) {
    pp <- sim_apply_quality_degradation(pp = pp, cfg = cfg, scenario_row = scenario_row, quality_abs = quality_abs)
  }

  marks_chr <- as.character(spatstat.geom::marks(pp))
  n_a <- sum(marks_chr == "A")
  n_b <- sum(marks_chr == "B")
  n_total <- length(marks_chr)

  if (n_total < cfg$min_cells_total || n_a < cfg$min_cells_per_type || n_b < cfg$min_cells_per_type) {
    stop("Insufficient cells after quality degradation")
  }

  if (isTRUE(use_bootstrap)) {
    tile_size_use <- max(
      1.5 * cfg$radius_um,
      cfg$tile_size * exp(scenario_row$quality_to_tile * quality_abs)
    )
    spe <- sim_ppp_to_spe(pp, sample_id = sample_id, cell_type_col = "cell_type")
    se_stats <- sim_run_spatialstats_one(
      spe = spe,
      sample_id = sample_id,
      cell_type_col = "cell_type",
      radii_um = cfg$radius_um,
      stat = cfg$stat,
      nsim = cfg$boot_nsim,
      correction = "translate",
      seed = seed,
      boot = cfg$boot_mode,
      tile_size = tile_size_use,
      BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
    )
    fx <- sim_flatten_spatialstats(se_stats) %>%
      dplyr::filter(
        ct1 %in% c("A", "B"),
        ct2 %in% c("A", "B"),
        radius_um == cfg$radius_um,
        stat == cfg$stat
      )
  } else {
    fx <- sim_calc_truth_yi_ab(pp = pp, radius_um = cfg$radius_um, stat = cfg$stat)
  }

  if (!nrow(fx)) stop("No A/B local_comp_enrichment rows returned")

  fx$sample <- sample_id
  fx$patient <- patient_row$patient
  fx$group <- "all"
  fx$tissue_size_um <- cfg$square_size_um
  fx$tissue_area_um2 <- cfg$tissue_area_um2
  fx$n_cells_A <- n_a
  fx$n_cells_B <- n_b
  fx$n_cells_total <- n_total
  fx$parent_lambda <- parent_lambda
  fx$offspring_mean <- offspring_mean
  fx$offspring_sd <- offspring_sd
  fx$quality_z <- quality_z
  fx$quality_abs <- quality_abs
  fx$contaminated <- contaminated
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
  target_cells = cfg$target_cells,
  radius_um = cfg$radius_um,
  stat = cfg$stat,
  boot_nsim = cfg$boot_nsim,
  boot_mode = cfg$boot_mode,
  tile_size = cfg$tile_size,
  patient_log_sd = cfg$patient_log_sd,
  sample_log_sd_base = cfg$sample_log_sd_base,
  parallel_backend = cfg$parallel_backend,
  seed_base = cfg$seed_base,
  n_workers = cfg$n_workers,
  scenarios = scenario_tbl
)

print(run_params[c(
  "cache_tag", "n_reps", "n_patients_grid", "n_truth_realizations_per_patient",
  "n_samples_per_patient_mean", "n_samples_dist", "square_size_um", "radius_um",
  "boot_nsim", "boot_mode", "tile_size", "patient_log_sd", "sample_log_sd_base",
  "parallel_backend", "n_workers"
)])
message("Scenarios:")
print(scenario_tbl)

max_n_patients <- max(cfg$n_patients_grid)
patient_latent <- sim_draw_patient_latents(max_n_patients, cfg, seed = cfg$seed_base + 17L)

bp <- sim_make_bpparam(cfg$n_workers, backend = cfg$parallel_backend)
bp_started <- FALSE
if (!BiocParallel::bpisup(bp)) {
  bp <- BiocParallel::bpstart(bp)
  bp_started <- TRUE
}
sim_log(
  "Parallel backend ready: ", class(bp)[1],
  " (workers=", BiocParallel::bpworkers(bp), ", requested=", cfg$parallel_backend, ")"
)
utils_file <- "analysis/simulation/sim_utils.R"
if (inherits(bp, "SnowParam")) {
  sim_log("SnowParam detected; simulation helpers will be sourced lazily on workers: ", utils_file)
}

truth_stage_t0 <- Sys.time()
sim_log("Precomputing truth (clean, no quality degradation) ...")
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
sim_log(
  "Truth plan: ", nrow(truth_plan), " jobs (",
  max_n_patients, " patients x ",
  cfg$n_truth_realizations_per_patient, " realizations/patient)"
)

truth_res <- BiocParallel::bplapply(
  seq_len(nrow(truth_plan)),
  function(j, plan_arg, lat_arg, cfg_arg, utils_file_arg,
           sim_clamp_fn, sim_apply_quality_degradation_fn, sim_calc_truth_yi_ab_fn, sim_eval_one_sample_fn) {
    if (!exists(".sim_breakpoint_worker_ready", envir = .GlobalEnv, inherits = FALSE)) {
      source(utils_file_arg, local = .GlobalEnv)
      assign("sim_clamp", sim_clamp_fn, envir = .GlobalEnv)
      assign("sim_apply_quality_degradation", sim_apply_quality_degradation_fn, envir = .GlobalEnv)
      assign("sim_calc_truth_yi_ab", sim_calc_truth_yi_ab_fn, envir = .GlobalEnv)
      assign("sim_eval_one_sample", sim_eval_one_sample_fn, envir = .GlobalEnv)
      assign(".sim_breakpoint_worker_ready", TRUE, envir = .GlobalEnv)
    }

    row <- plan_arg[j, , drop = FALSE]
    p_row <- lat_arg[row$patient_idx[[1]], , drop = FALSE]
    out <- tryCatch({
      sim_eval_one_sample(
        patient_row = p_row,
        sample_id = row$sample[[1]],
        cfg = cfg_arg,
        seed = as.integer(row$seed[[1]]),
        scenario_row = NULL,
        contaminated = FALSE,
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
  utils_file_arg = utils_file,
  sim_clamp_fn = sim_clamp,
  sim_apply_quality_degradation_fn = sim_apply_quality_degradation,
  sim_calc_truth_yi_ab_fn = sim_calc_truth_yi_ab,
  sim_eval_one_sample_fn = sim_eval_one_sample,
  BPPARAM = bp
)

truth_ok <- vapply(truth_res, function(x) is.list(x) && isTRUE(x$ok), logical(1))
truth_fail <- dplyr::bind_rows(lapply(truth_res[!truth_ok], function(x) {
  data.frame(
    stage = "truth",
    scenario = NA_character_,
    replicate = NA_integer_,
    n_patients = NA_integer_,
    patient = if (!is.null(x$patient)) x$patient else NA_character_,
    sample = if (!is.null(x$sample)) x$sample else NA_character_,
    error = if (!is.null(x$error)) x$error else "unknown truth error",
    stringsAsFactors = FALSE
  )
}))
if (!any(truth_ok)) {
  truth_fail_csv <- file.path(out_tbl_dir, "intersample_breakpoint_truth_failures.csv")
  if (nrow(truth_fail) > 0L) sim_safe_write_csv(truth_fail, truth_fail_csv)

  top_truth_msg <- "no captured worker errors"
  if (nrow(truth_fail) > 0L) {
    top_truth <- truth_fail %>%
      dplyr::count(error, sort = TRUE) %>%
      dplyr::slice_head(n = 5)
    top_truth_msg <- paste(sprintf("[%d] %s", top_truth$n, top_truth$error), collapse = " | ")
  }

  stop(
    "No truth rows generated. Top truth errors: ", top_truth_msg,
    if (nrow(truth_fail) > 0L) paste0(". Failure log: ", truth_fail_csv) else "",
    ". If running interactively, clear stale overrides with rm(sim_cfg)."
  )
}
df_truth <- dplyr::bind_rows(lapply(truth_res[truth_ok], `[[`, "data"))
sim_log(
  "Truth stage complete in ", sim_elapsed(truth_stage_t0), ": ",
  sum(truth_ok), "/", length(truth_ok), " jobs succeeded",
  if (any(!truth_ok)) paste0(" (", sum(!truth_ok), " failed)") else ""
)

truth_patient <- df_truth %>%
  dplyr::group_by(feature_id, ct1, ct2, radius_um, stat, patient) %>%
  dplyr::summarise(
    true_patient_mean = mean(yi, na.rm = TRUE),
    true_patient_var = stats::var(yi, na.rm = TRUE),
    .groups = "drop"
  )

obs_stage_t0 <- Sys.time()
sim_log("Generating observed data and fitting models ...")

est_rows <- list()
fail_rows <- list()
sample_meta_rows <- list()
if (nrow(truth_fail) > 0L) fail_rows[[length(fail_rows) + 1L]] <- truth_fail

total_steps <- nrow(scenario_tbl) * cfg$n_reps * length(cfg$n_patients_grid)
step_counter <- 0L
pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
sim_log(
  "Observed/meta plan: ", nrow(scenario_tbl), " scenarios x ",
  cfg$n_reps, " replicates x ", length(cfg$n_patients_grid),
  " patient-grid settings = ", total_steps, " meta fit steps"
)

for (s_idx in seq_len(nrow(scenario_tbl))) {
  scenario_t0 <- Sys.time()
  scn <- scenario_tbl[s_idx, , drop = FALSE]
  scenario_label <- scn$scenario[[1]]
  sim_log("Scenario ", s_idx, "/", nrow(scenario_tbl), ": ", scenario_label)

  for (rep_idx in seq_len(cfg$n_reps)) {
    rep_t0 <- Sys.time()
    sim_log("  Replicate ", rep_idx, "/", cfg$n_reps, " (", scenario_label, "): drawing observed samples")
    contaminated_flag <- stats::runif(max_n_patients) < scn$contamination_frac[[1]]

    n_samples_vec <- vapply(
      seq_len(max_n_patients),
      function(p_idx) {
        sim_draw_n_samples_per_patient(cfg, seed = cfg$seed_base + s_idx * 1000000L + rep_idx * 10000L + p_idx * 101L)
      },
      integer(1)
    )

    obs_plan <- do.call(rbind, lapply(seq_len(max_n_patients), function(p_idx) {
      data.frame(
        patient_idx = p_idx,
        sample_ord = seq_len(n_samples_vec[[p_idx]]),
        stringsAsFactors = FALSE
      )
    }))
    obs_plan$patient <- patient_latent$patient[obs_plan$patient_idx]
    obs_plan$sample <- sprintf("%s_%s_rep%03d_s%02d", obs_plan$patient, scenario_label, rep_idx, obs_plan$sample_ord)
    obs_plan$seed <- cfg$seed_base + 2000000L + s_idx * 1000000L + rep_idx * 10000L + obs_plan$patient_idx * 101L + obs_plan$sample_ord
    obs_plan$contaminated <- contaminated_flag[obs_plan$patient_idx]
    sim_log(
      "    Observed sample jobs: ", nrow(obs_plan),
      " [samples/patient min/med/max=",
      min(n_samples_vec), "/", as.integer(stats::median(n_samples_vec)), "/", max(n_samples_vec), "]"
    )
    sim_log(
      "    Contaminated observed jobs: ", sum(obs_plan$contaminated), "/", nrow(obs_plan),
      " (target_frac=", signif(scn$contamination_frac[[1]], 3), ")"
    )

    obs_res <- BiocParallel::bplapply(
      seq_len(nrow(obs_plan)),
      function(j, plan_arg, lat_arg, cfg_arg, scn_arg, rep_arg, scn_label, utils_file_arg,
               sim_clamp_fn, sim_apply_quality_degradation_fn, sim_calc_truth_yi_ab_fn, sim_eval_one_sample_fn) {
        if (!exists(".sim_breakpoint_worker_ready", envir = .GlobalEnv, inherits = FALSE)) {
          source(utils_file_arg, local = .GlobalEnv)
          assign("sim_clamp", sim_clamp_fn, envir = .GlobalEnv)
          assign("sim_apply_quality_degradation", sim_apply_quality_degradation_fn, envir = .GlobalEnv)
          assign("sim_calc_truth_yi_ab", sim_calc_truth_yi_ab_fn, envir = .GlobalEnv)
          assign("sim_eval_one_sample", sim_eval_one_sample_fn, envir = .GlobalEnv)
          assign(".sim_breakpoint_worker_ready", TRUE, envir = .GlobalEnv)
        }

        row <- plan_arg[j, , drop = FALSE]
        p_row <- lat_arg[row$patient_idx[[1]], , drop = FALSE]

        out <- NULL
        err_msg <- NA_character_
        for (attempt in seq_len(cfg_arg$max_sample_attempts)) {
          out_try <- tryCatch({
            sim_eval_one_sample(
              patient_row = p_row,
              sample_id = row$sample[[1]],
              cfg = cfg_arg,
              seed = as.integer(row$seed[[1]] + attempt - 1L),
              scenario_row = scn_arg,
              contaminated = isTRUE(row$contaminated[[1]]),
              use_bootstrap = TRUE
            )
          }, error = function(e) e)
          if (!inherits(out_try, "error")) {
            out <- out_try
            break
          }
          err_msg <- conditionMessage(out_try)
        }

        if (is.null(out)) {
          return(list(
            ok = FALSE,
            error = err_msg,
            stage = "observed",
            scenario = scn_label,
            replicate = rep_arg,
            n_patients = NA_integer_,
            patient = row$patient[[1]],
            sample = row$sample[[1]]
          ))
        }

        out$scenario <- scn_label
        out$replicate <- rep_arg

        sample_meta <- out %>%
          dplyr::distinct(
            scenario, replicate, patient, sample, contaminated,
            n_cells_A, n_cells_B, n_cells_total,
            quality_z, quality_abs, parent_lambda, offspring_mean, offspring_sd
          )

        list(ok = TRUE, data = out, sample_meta = sample_meta)
      },
      plan_arg = obs_plan,
      lat_arg = patient_latent,
      cfg_arg = cfg,
      scn_arg = scn,
      rep_arg = rep_idx,
      scn_label = scenario_label,
      utils_file_arg = utils_file,
      sim_clamp_fn = sim_clamp,
      sim_apply_quality_degradation_fn = sim_apply_quality_degradation,
      sim_calc_truth_yi_ab_fn = sim_calc_truth_yi_ab,
      sim_eval_one_sample_fn = sim_eval_one_sample,
      BPPARAM = bp
    )

    obs_ok <- vapply(obs_res, function(x) is.list(x) && isTRUE(x$ok), logical(1))
    obs_fail <- dplyr::bind_rows(lapply(obs_res[!obs_ok], function(x) {
      data.frame(
        stage = if (!is.null(x$stage)) x$stage else "observed",
        scenario = if (!is.null(x$scenario)) x$scenario else scenario_label,
        replicate = if (!is.null(x$replicate)) x$replicate else rep_idx,
        n_patients = if (!is.null(x$n_patients)) x$n_patients else NA_integer_,
        patient = if (!is.null(x$patient)) x$patient else NA_character_,
        sample = if (!is.null(x$sample)) x$sample else NA_character_,
        error = if (!is.null(x$error)) x$error else "unknown observed error",
        stringsAsFactors = FALSE
      )
    }))
    if (nrow(obs_fail) > 0L) fail_rows[[length(fail_rows) + 1L]] <- obs_fail
    sim_log(
      "    Observed generation complete in ", sim_elapsed(rep_t0), ": ",
      sum(obs_ok), "/", length(obs_ok), " sample jobs succeeded",
      if (nrow(obs_fail) > 0L) paste0(" (", nrow(obs_fail), " failed)") else ""
    )

    if (!any(obs_ok)) {
      for (dummy in cfg$n_patients_grid) {
        step_counter <- step_counter + 1L
        utils::setTxtProgressBar(pb, step_counter)
      }
      next
    }

    df_obs_rep <- dplyr::bind_rows(lapply(obs_res[obs_ok], `[[`, "data"))
    df_sample_meta_rep <- dplyr::bind_rows(lapply(obs_res[obs_ok], `[[`, "sample_meta"))

    for (n_pat in cfg$n_patients_grid) {
      sim_log("      Fitting meta for n_patients=", n_pat)
      keep_patients <- patient_latent$patient[seq_len(n_pat)]
      df_sub <- df_obs_rep %>% dplyr::filter(patient %in% keep_patients)

      if (!nrow(df_sub)) {
        fail_rows[[length(fail_rows) + 1L]] <- data.frame(
          stage = "subset_empty",
          scenario = scenario_label,
          replicate = rep_idx,
          n_patients = n_pat,
          patient = NA_character_,
          sample = NA_character_,
          error = "No observed rows after patient subset filtering",
          stringsAsFactors = FALSE
        )
        step_counter <- step_counter + 1L
        utils::setTxtProgressBar(pb, step_counter)
        next
      }

      sub_sample_meta <- df_sample_meta_rep %>%
        dplyr::filter(patient %in% keep_patients) %>%
        dplyr::mutate(n_patients = n_pat)
      sample_meta_rows[[length(sample_meta_rows) + 1L]] <- sub_sample_meta

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
          median_vi = stats::median(vi, na.rm = TRUE),
          mean_n_cells_A = mean(n_cells_A, na.rm = TRUE),
          mean_n_cells_B = mean(n_cells_B, na.rm = TRUE),
          mean_n_cells_total = mean(n_cells_total, na.rm = TRUE),
          contamination_frac_obs = mean(contaminated, na.rm = TRUE),
          mean_quality_abs = mean(quality_abs, na.rm = TRUE),
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
          scenario = scenario_label,
          replicate = rep_idx,
          n_patients = n_pat,
          patient = NA_character_,
          sample = NA_character_,
          error = conditionMessage(se_meta),
          stringsAsFactors = FALSE
        )
        step_counter <- step_counter + 1L
        utils::setTxtProgressBar(pb, step_counter)
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
          scenario = scenario_label,
          n_patients = n_pat,
          replicate = rep_idx
        )

      if (!nrow(est)) {
        fail_rows[[length(fail_rows) + 1L]] <- data.frame(
          stage = "meta_extract",
          scenario = scenario_label,
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

      step_counter <- step_counter + 1L
      utils::setTxtProgressBar(pb, step_counter)
    }
    sim_log("  Replicate ", rep_idx, "/", cfg$n_reps, " complete in ", sim_elapsed(rep_t0))
  }
  sim_log("Scenario ", scenario_label, " complete in ", sim_elapsed(scenario_t0))
}
close(pb)
sim_log("Observed/meta stage complete in ", sim_elapsed(obs_stage_t0))

if (isTRUE(bp_started) && BiocParallel::bpisup(bp)) {
  try(BiocParallel::bpstop(bp), silent = TRUE)
}

df_est <- dplyr::bind_rows(est_rows)
df_fail <- dplyr::bind_rows(fail_rows)
df_sample_meta <- dplyr::bind_rows(sample_meta_rows)

if (!nrow(df_est)) {
  fail_csv <- file.path(out_tbl_dir, "intersample_breakpoint_failures.csv")
  if (nrow(df_fail) > 0L) sim_safe_write_csv(df_fail, fail_csv)
  stop("No estimate rows were produced. See failures: ", fail_csv)
}

cmp_long <- dplyr::bind_rows(
  df_est %>%
    dplyr::transmute(scenario, n_patients, replicate, feature_id, parameter = "mu", method = "panoramic", estimate = mu_hat, truth = true_mu),
  df_est %>%
    dplyr::transmute(scenario, n_patients, replicate, feature_id, parameter = "mu", method = "naive", estimate = mu_naive, truth = true_mu),
  df_est %>%
    dplyr::transmute(scenario, n_patients, replicate, feature_id, parameter = "tau2_patient", method = "panoramic", estimate = tau2_patient, truth = true_tau2_patient),
  df_est %>%
    dplyr::transmute(scenario, n_patients, replicate, feature_id, parameter = "tau2_patient", method = "naive", estimate = tau2_patient_naive, truth = true_tau2_patient),
  df_est %>%
    dplyr::transmute(scenario, n_patients, replicate, feature_id, parameter = "tau2_sample", method = "panoramic", estimate = tau2_sample, truth = true_tau2_sample),
  df_est %>%
    dplyr::transmute(scenario, n_patients, replicate, feature_id, parameter = "tau2_sample", method = "naive", estimate = tau2_sample_naive, truth = true_tau2_sample)
) %>%
  dplyr::mutate(error = estimate - truth, abs_error = abs(error))

perf_by_replicate <- cmp_long %>%
  dplyr::group_by(scenario, n_patients, parameter, method, replicate) %>%
  dplyr::summarise(
    bias = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mae = mean(abs_error, na.rm = TRUE),
    n_features = dplyr::n_distinct(feature_id),
    .groups = "drop"
  )

perf_by_feature <- cmp_long %>%
  dplyr::group_by(scenario, n_patients, parameter, method, feature_id) %>%
  dplyr::summarise(
    bias = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mae = mean(abs_error, na.rm = TRUE),
    .groups = "drop"
  )

feature_meta <- df_est %>%
  dplyr::distinct(feature_id, ct1, ct2, radius_um, stat) %>%
  dplyr::mutate(celltype_pair = paste(ct1, ct2, sep = "->"))

perf_by_pair <- perf_by_feature %>%
  dplyr::left_join(feature_meta, by = "feature_id") %>%
  dplyr::select(
    scenario, n_patients, parameter, method,
    feature_id, ct1, ct2, celltype_pair, radius_um, stat,
    bias, rmse, mae
  )

perf_summary <- perf_by_replicate %>%
  dplyr::group_by(scenario, n_patients, parameter, method) %>%
  dplyr::summarise(
    bias = mean(bias, na.rm = TRUE),
    rmse_mean = mean(rmse, na.rm = TRUE),
    mae = mean(mae, na.rm = TRUE),
    n_features = max(n_features, na.rm = TRUE),
    n_reps = sum(is.finite(rmse)),
    rmse_sd = stats::sd(rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rmse_se = dplyr::if_else(n_reps > 0L, rmse_sd / sqrt(n_reps), NA_real_),
    rmse_ci_mult = dplyr::if_else(n_reps > 1L, stats::qt(0.975, df = n_reps - 1L), NA_real_),
    rmse_ci_lower = dplyr::if_else(
      n_reps > 1L,
      pmax(0, rmse_mean - rmse_ci_mult * rmse_se),
      rmse_mean
    ),
    rmse_ci_upper = dplyr::if_else(
      n_reps > 1L,
      rmse_mean + rmse_ci_mult * rmse_se,
      rmse_mean
    ),
    rmse = rmse_mean
  ) %>%
  dplyr::select(-rmse_mean, -rmse_ci_mult)

pair_cmp <- cmp_long %>%
  dplyr::select(scenario, n_patients, replicate, feature_id, parameter, method, abs_error) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = abs_error,
    names_prefix = "abs_error__"
  ) %>%
  dplyr::mutate(
    meta_better = abs_error__panoramic < abs_error__naive,
    naive_better = abs_error__naive < abs_error__panoramic,
    abs_error_margin = abs_error__naive - abs_error__panoramic
  )

breakpoint_signal <- pair_cmp %>%
  dplyr::group_by(scenario, n_patients, parameter) %>%
  dplyr::summarise(
    p_meta_better = mean(meta_better, na.rm = TRUE),
    mean_abs_error_margin = mean(abs_error_margin, na.rm = TRUE),
    median_abs_error_margin = stats::median(abs_error_margin, na.rm = TRUE),
    .groups = "drop"
  )

baseline_rmse <- perf_summary %>%
  dplyr::filter(scenario == "baseline") %>%
  dplyr::select(n_patients, parameter, method, baseline_rmse = rmse)

breakpoint_tbl <- perf_summary %>%
  dplyr::left_join(baseline_rmse, by = c("n_patients", "parameter", "method")) %>%
  dplyr::left_join(
    perf_summary %>%
      dplyr::select(scenario, n_patients, parameter, method, rmse) %>%
      tidyr::pivot_wider(names_from = method, values_from = rmse, names_prefix = "rmse_"),
    by = c("scenario", "n_patients", "parameter")
  ) %>%
  dplyr::left_join(breakpoint_signal, by = c("scenario", "n_patients", "parameter")) %>%
  dplyr::mutate(
    rmse_inflation_vs_baseline = rmse / baseline_rmse,
    rmse_ratio_naive_over_meta = rmse_naive / rmse_panoramic,
    break_flag = (!is.na(rmse_ratio_naive_over_meta) & rmse_ratio_naive_over_meta >= 1.25) |
      (!is.na(p_meta_better) & p_meta_better >= 0.70)
  )

breakpoint_first <- breakpoint_tbl %>%
  dplyr::filter(method == "naive", scenario != "baseline") %>%
  dplyr::arrange(scenario, parameter, n_patients) %>%
  dplyr::group_by(scenario, parameter) %>%
  dplyr::summarise(
    first_break_n_patients = {
      idx <- which(break_flag)
      if (length(idx) > 0L) n_patients[min(idx)] else NA_integer_
    },
    max_rmse_ratio_naive_over_meta = max(rmse_ratio_naive_over_meta, na.rm = TRUE),
    max_p_meta_better = max(p_meta_better, na.rm = TRUE),
    .groups = "drop"
  )

p_win <- ggplot(
  breakpoint_signal,
  aes(x = factor(n_patients), y = scenario, fill = p_meta_better)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ parameter, ncol = 1) +
  scale_fill_gradient(low = "#E69F00", high = "#0072B2", limits = c(0, 1)) +
  labs(
    title = "Probability PANORAMIC has lower absolute error",
    x = "Number of patients",
    y = "Scenario",
    fill = "P(meta\nbetter)"
  ) +
  theme_minimal()

p_rmse <- ggplot(
  perf_summary,
  aes(x = n_patients, y = rmse, color = method, group = method)
) +
  geom_ribbon(
    aes(ymin = rmse_ci_lower, ymax = rmse_ci_upper, fill = method),
    alpha = 0.18,
    linewidth = 0,
    color = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  facet_grid(parameter ~ scenario, scales = "free_y") +
  scale_color_manual(values = c(panoramic = "#0072B2", naive = "#D55E00")) +
  scale_fill_manual(values = c(panoramic = "#0072B2", naive = "#D55E00"), guide = "none") +
  labs(
    title = "RMSE comparison across breakpoint scenarios",
    x = "Number of patients",
    y = "RMSE",
    color = "Method"
  ) +
  theme_minimal()

p_rmse_pair <- ggplot(
  perf_by_pair,
  aes(x = n_patients, y = rmse, color = celltype_pair, linetype = method, group = interaction(method, celltype_pair))
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(parameter ~ scenario, scales = "free_y") +
  scale_color_manual(values = c("A->A" = "#1B9E77", "A->B" = "#D95F02", "B->A" = "#7570B3", "B->B" = "#E7298A")) +
  labs(
    title = "RMSE broken down by cell-type pair",
    x = "Number of patients",
    y = "RMSE",
    color = "Cell-type pair",
    linetype = "Method"
  ) +
  theme_minimal()

tau2_levels <- c("tau2_patient", "tau2_sample")

p_rmse_tau2 <- ggplot(
  perf_summary %>% dplyr::filter(parameter %in% tau2_levels),
  aes(x = n_patients, y = rmse, color = method, group = method)
) +
  geom_ribbon(
    aes(ymin = rmse_ci_lower, ymax = rmse_ci_upper, fill = method),
    alpha = 0.18,
    linewidth = 0,
    color = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  facet_grid(parameter ~ scenario, scales = "free_y") +
  scale_color_manual(values = c(panoramic = "#0072B2", naive = "#D55E00")) +
  scale_fill_manual(values = c(panoramic = "#0072B2", naive = "#D55E00"), guide = "none") +
  labs(
    title = "RMSE comparison across breakpoint scenarios for tau2",
    x = "Number of patients",
    y = "RMSE",
    color = "Method"
  ) +
  theme_minimal()

p_rmse_pair_tau2 <- ggplot(
  perf_by_pair %>% dplyr::filter(parameter %in% tau2_levels),
  aes(x = n_patients, y = rmse, color = celltype_pair, linetype = method, group = interaction(method, celltype_pair))
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(parameter ~ scenario, scales = "free_y") +
  scale_color_manual(values = c("A->A" = "#1B9E77", "A->B" = "#D95F02", "B->A" = "#7570B3", "B->B" = "#E7298A")) +
  labs(
    title = "RMSE by cell-type pair across breakpoint scenarios for tau2",
    x = "Number of patients",
    y = "RMSE",
    color = "Cell-type pair",
    linetype = "Method"
  ) +
  theme_minimal()

params_rds <- file.path(out_data_dir, "intersample_breakpoint_params.rds")
truth_rds <- file.path(out_data_dir, "intersample_breakpoint_truth.rds")
est_rds <- file.path(out_data_dir, "intersample_breakpoint_estimates.rds")
meta_rds <- file.path(out_data_dir, "intersample_breakpoint_sample_meta.rds")

est_csv <- file.path(out_tbl_dir, "intersample_breakpoint_estimates.csv")
fail_csv <- file.path(out_tbl_dir, "intersample_breakpoint_failures.csv")
sample_meta_csv <- file.path(out_tbl_dir, "intersample_breakpoint_sample_meta.csv")
perf_replicate_csv <- file.path(out_tbl_dir, "intersample_breakpoint_perf_by_replicate.csv")
perf_feature_csv <- file.path(out_tbl_dir, "intersample_breakpoint_perf_by_feature.csv")
perf_pair_csv <- file.path(out_tbl_dir, "intersample_breakpoint_perf_by_celltype_pair.csv")
perf_summary_csv <- file.path(out_tbl_dir, "intersample_breakpoint_perf_summary.csv")
breakpoint_signal_csv <- file.path(out_tbl_dir, "intersample_breakpoint_signal.csv")
breakpoint_csv <- file.path(out_tbl_dir, "intersample_breakpoint_table.csv")
breakpoint_first_csv <- file.path(out_tbl_dir, "intersample_breakpoint_first_break.csv")

win_png <- file.path(out_fig_dir, "intersample_breakpoint_meta_win_heatmap.png")
rmse_png <- file.path(out_fig_dir, "intersample_breakpoint_rmse_lines.png")
rmse_pair_png <- file.path(out_fig_dir, "intersample_breakpoint_rmse_by_celltype_pair.png")
rmse_tau2_png <- file.path(out_fig_dir, "intersample_breakpoint_rmse_tau2_lines.png")
rmse_pair_tau2_png <- file.path(out_fig_dir, "intersample_breakpoint_rmse_tau2_by_celltype_pair.png")

sim_safe_save_rds(run_params, params_rds)
sim_safe_save_rds(list(truth_long = df_truth, truth_patient = truth_patient), truth_rds)
sim_safe_save_rds(df_est, est_rds)
sim_safe_save_rds(df_sample_meta, meta_rds)

sim_safe_write_csv(df_est, est_csv)
if (nrow(df_fail) > 0L) sim_safe_write_csv(df_fail, fail_csv)
if (nrow(df_sample_meta) > 0L) sim_safe_write_csv(df_sample_meta, sample_meta_csv)
sim_safe_write_csv(perf_by_replicate, perf_replicate_csv)
sim_safe_write_csv(perf_by_feature, perf_feature_csv)
sim_safe_write_csv(perf_by_pair, perf_pair_csv)
sim_safe_write_csv(perf_summary, perf_summary_csv)
sim_safe_write_csv(breakpoint_signal, breakpoint_signal_csv)
sim_safe_write_csv(breakpoint_tbl, breakpoint_csv)
sim_safe_write_csv(breakpoint_first, breakpoint_first_csv)

ggsave(win_png, p_win, width = 10.5, height = 7.5, dpi = 300)
ggsave(rmse_png, p_rmse, width = 13.5, height = 8.5, dpi = 300)
ggsave(rmse_pair_png, p_rmse_pair, width = 14.5, height = 9.5, dpi = 300)
ggsave(rmse_tau2_png, p_rmse_tau2, width = 13.5, height = 8.5, dpi = 300)
ggsave(rmse_pair_tau2_png, p_rmse_pair_tau2, width = 14.5, height = 9.5, dpi = 300)

message("\nSaved intersample breakpoint outputs:")
message("  params: ", params_rds)
message("  truth: ", truth_rds)
message("  estimates: ", est_rds)
message("  sample meta: ", meta_rds)
message("  performance summary: ", perf_summary_csv)
message("  performance by cell-type pair: ", perf_pair_csv)
message("  breakpoint summary: ", breakpoint_first_csv)
message("  figures: ", win_png, " ; ", rmse_png, " ; ", rmse_pair_png, " ; ", rmse_tau2_png, " ; ", rmse_pair_tau2_png)
if (nrow(df_fail) > 0L) message("  failures: ", fail_csv)

if (interactive()) {
  assign("sim_breakpoint", list(
    config = cfg,
    scenarios = scenario_tbl,
    estimates = df_est,
    sample_meta = df_sample_meta,
    failures = df_fail,
    performance_by_feature = perf_by_feature,
    performance_by_pair = perf_by_pair,
    performance_summary = perf_summary,
    breakpoint_signal = breakpoint_signal,
    breakpoint_table = breakpoint_tbl,
    breakpoint_first = breakpoint_first,
    plots = list(
      meta_win_heatmap = p_win,
      rmse_lines = p_rmse,
      rmse_by_celltype_pair = p_rmse_pair,
      rmse_tau2_lines = p_rmse_tau2,
      rmse_tau2_by_celltype_pair = p_rmse_pair_tau2
    )
  ), envir = .GlobalEnv)

  message("\nInteractive object created: sim_breakpoint")
  print(head(breakpoint_first, 20))
  print(p_win)
  print(p_rmse)
  print(p_rmse_pair)
  print(p_rmse_tau2)
  print(p_rmse_pair_tau2)
}
