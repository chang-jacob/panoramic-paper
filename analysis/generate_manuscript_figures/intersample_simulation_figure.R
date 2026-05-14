#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

find_project_root <- function(start = getwd(), max_depth = 10L) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_depth)) {
    if (file.exists(file.path(cur, "config", "default.yml"))) return(cur)
    nxt <- dirname(cur)
    if (identical(nxt, cur)) break
    cur <- nxt
  }
  stop("Could not find project root containing config/default.yml from: ", start)
}

pick_existing <- function(paths, label = "file") {
  for (p in paths) {
    if (file.exists(p)) return(p)
  }
  stop("Could not find ", label, ". Checked:\n- ", paste(paths, collapse = "\n- "))
}

save_plot_formats <- function(plot_obj, stem, width, height, dpi = 300) {
  ggsave(
    filename = paste0(stem, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  ggsave(
    filename = paste0(stem, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in"
  )
}

sim_clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))

sim_draw_patient_latents_preview <- function(n_patients, cfg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  strength <- stats::rnorm(n_patients, mean = 0, sd = cfg$patient_log_sd)
  data.frame(
    patient = sprintf("patient_%03d", seq_len(n_patients)),
    strength = strength,
    parent_lambda_patient = pmax(cfg$target_cell_density / 450, cfg$base_parent_lambda * exp(strength)),
    offspring_mean_patient = pmax(1, cfg$base_offspring_mean * exp(0.7 * strength)),
    offspring_sd_patient = pmax(8, cfg$base_offspring_sd * exp(-0.25 * strength)),
    stringsAsFactors = FALSE
  )
}

sim_apply_quality_degradation_preview <- function(pp, cfg, scenario_row, quality_abs) {
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

  spatstat.geom::ppp(
    x = df_pp$x[keep],
    y = df_pp$y[keep],
    window = spatstat.geom::square(cfg$square_size_um),
    marks = factor(marks_keep, levels = c("A", "B"))
  )
}

sim_eval_preview_sample <- function(patient_row, scenario_row, cfg, seed, contaminated = FALSE) {
  set.seed(seed)

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
    pp <- sim_apply_quality_degradation_preview(
      pp = pp,
      cfg = cfg,
      scenario_row = scenario_row,
      quality_abs = quality_abs
    )
  }

  marks_chr <- as.character(spatstat.geom::marks(pp))
  n_a <- sum(marks_chr == "A")
  n_b <- sum(marks_chr == "B")
  n_total <- length(marks_chr)
  if (n_total < cfg$min_cells_total || n_a < cfg$min_cells_per_type || n_b < cfg$min_cells_per_type) {
    stop("Insufficient cells after quality degradation.")
  }

  list(
    pp = pp,
    contaminated = isTRUE(contaminated),
    quality_z = quality_z,
    quality_abs = quality_abs,
    parent_lambda = parent_lambda,
    offspring_mean = offspring_mean,
    offspring_sd = offspring_sd,
    n_cells_A = n_a,
    n_cells_B = n_b,
    n_cells_total = n_total
  )
}

sim_eval_preview_with_retry <- function(patient_row, scenario_row, cfg, seed, contaminated = FALSE) {
  last_err <- "unknown error"
  for (attempt in seq_len(cfg$max_sample_attempts)) {
    out <- tryCatch(
      sim_eval_preview_sample(
        patient_row = patient_row,
        scenario_row = scenario_row,
        cfg = cfg,
        seed = as.integer(seed + attempt - 1L),
        contaminated = contaminated
      ),
      error = function(e) e
    )
    if (!inherits(out, "error")) return(out)
    last_err <- conditionMessage(out)
  }
  stop(last_err)
}

project_root <- find_project_root()
source(file.path(project_root, "analysis", "simulation", "sim_utils.R"))
paths <- sim_load_paths(project_root)

out_tbl_dir <- sim_mkdir(file.path(paths$output, "tables", "simulation", "intersample", "publication"))
out_fig_dir <- sim_mkdir(file.path(paths$output, "figures", "simulation", "intersample", "publication"))

perf_summary_path <- pick_existing(
  c(
    file.path(paths$output, "tables", "simulation", "intersample_breakpoint", "intersample_breakpoint_perf_summary.csv"),
    file.path(project_root, "output", "tables", "simulation", "intersample_breakpoint", "intersample_breakpoint_perf_summary.csv")
  ),
  label = "intersample breakpoint performance summary CSV"
)

params_path <- pick_existing(
  c(
    file.path(paths$data, "interim", "simulation", "intersample_breakpoint", "intersample_breakpoint_params.rds"),
    file.path(project_root, "data", "interim", "simulation", "intersample_breakpoint", "intersample_breakpoint_params.rds")
  ),
  label = "intersample breakpoint params RDS"
)

estimates_path <- pick_existing(
  c(
    file.path(paths$data, "interim", "simulation", "intersample_breakpoint", "intersample_breakpoint_estimates.rds"),
    file.path(project_root, "data", "interim", "simulation", "intersample_breakpoint", "intersample_breakpoint_estimates.rds")
  ),
  label = "intersample breakpoint estimates RDS"
)

perf_summary <- readr::read_csv(perf_summary_path, show_col_types = FALSE)
run_params <- readRDS(params_path)
breakpoint_estimates <- readRDS(estimates_path)

need_perf <- c("scenario", "n_patients", "parameter", "method", "rmse", "bias", "mae", "n_features")
miss_perf <- setdiff(need_perf, colnames(perf_summary))
if (length(miss_perf) > 0L) {
  stop("Performance summary is missing required columns: ", paste(miss_perf, collapse = ", "))
}

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

scenario_pretty <- c(
  baseline = "Baseline",
  mild_break = "Mild",
  moderate_break = "Moderate",
  severe_break = "Severe",
  extreme_break = "Extreme"
)

perf_scenarios <- sort(unique(as.character(perf_summary$scenario)))
scenario_levels <- scenario_tbl$scenario[scenario_tbl$scenario %in% perf_scenarios]
if (length(scenario_levels) < 1L) {
  stop(
    "No overlap between breakpoint scenario ladder and performance summary scenarios.\n",
    "Expected at least one of: ", paste(scenario_tbl$scenario, collapse = ", "), "\n",
    "Observed: ", paste(perf_scenarios, collapse = ", ")
  )
}
unknown_perf_scenarios <- setdiff(perf_scenarios, scenario_tbl$scenario)
if (length(unknown_perf_scenarios) > 0L) {
  message(
    "Ignoring non-ladder scenarios in performance summary: ",
    paste(unknown_perf_scenarios, collapse = ", ")
  )
}

scenario_display_map <- setNames(
  paste0(
    ifelse(
      scenario_tbl$scenario %in% names(scenario_pretty),
      scenario_pretty[scenario_tbl$scenario],
      scenario_tbl$scenario
    ),
    " (",
    sprintf("%d%%", as.integer(round(100 * scenario_tbl$contamination_frac))),
    " contaminated)"
  ),
  scenario_tbl$scenario
)

scenario_display_levels <- unname(scenario_display_map[scenario_levels])
scenario_use <- scenario_tbl %>%
  dplyr::filter(scenario %in% scenario_levels)

cfg <- list(
  square_size_um = as.numeric(if (!is.null(run_params$square_size_um)) run_params$square_size_um else 600),
  target_cells = as.numeric(if (!is.null(run_params$target_cells)) run_params$target_cells else 1900),
  patient_log_sd = as.numeric(if (!is.null(run_params$patient_log_sd)) run_params$patient_log_sd else 0.35),
  sample_log_sd_base = as.numeric(if (!is.null(run_params$sample_log_sd_base)) run_params$sample_log_sd_base else 0.20),
  base_offspring_mean = as.numeric(if (!is.null(run_params$base_offspring_mean)) run_params$base_offspring_mean else 12),
  base_offspring_sd = as.numeric(if (!is.null(run_params$base_offspring_sd)) run_params$base_offspring_sd else 55),
  min_cells_total = as.integer(if (!is.null(run_params$min_cells_total)) run_params$min_cells_total else 130L),
  min_cells_per_type = as.integer(if (!is.null(run_params$min_cells_per_type)) run_params$min_cells_per_type else 28L),
  min_keep = as.numeric(if (!is.null(run_params$min_keep)) run_params$min_keep else 0.10),
  max_keep = as.numeric(if (!is.null(run_params$max_keep)) run_params$max_keep else 0.99),
  max_thinning_attempts = as.integer(if (!is.null(run_params$max_thinning_attempts)) run_params$max_thinning_attempts else 5L),
  max_sample_attempts = as.integer(if (!is.null(run_params$max_sample_attempts)) run_params$max_sample_attempts else 4L),
  seed_base = as.integer(if (!is.null(run_params$seed_base)) run_params$seed_base else 811L)
)
cfg$tissue_area_um2 <- cfg$square_size_um^2
cfg$target_cell_density <- cfg$target_cells / cfg$tissue_area_um2
cfg$base_parent_lambda <- cfg$target_cell_density / (1 + cfg$base_offspring_mean)

analysis_params <- list(
  multi_patient_n_patients = 3L,
  multi_patient_n_samples_per_patient = c(4L, 3L, 5L),
  degradation_patient_index = 1L,
  degradation_sample_index = 1L,
  multi_point_size = 0.08,
  progressive_point_size = 0.025,
  multi_plot_keep_frac = 0.25,
  progressive_plot_keep_frac = 0.25,
  point_alpha = 0.78
)

thin_plot_points <- function(df, keep_frac, seed) {
  keep_frac <- as.numeric(keep_frac)
  if (!is.finite(keep_frac) || keep_frac <= 0 || keep_frac >= 1) {
    return(df)
  }
  n_keep <- max(1L, floor(nrow(df) * keep_frac))
  if (n_keep >= nrow(df)) return(df)
  set.seed(as.integer(seed))
  dplyr::slice_sample(df, n = n_keep)
}

n_preview_patients <- as.integer(analysis_params$multi_patient_n_patients)
samples_per_patient <- as.integer(analysis_params$multi_patient_n_samples_per_patient)
if (length(samples_per_patient) == 1L) {
  samples_per_patient <- rep(samples_per_patient, n_preview_patients)
}
if (length(samples_per_patient) != n_preview_patients) {
  stop(
    "analysis_params$multi_patient_n_samples_per_patient must have length 1 or ",
    "match multi_patient_n_patients."
  )
}
if (any(!is.finite(samples_per_patient)) || any(samples_per_patient < 1L)) {
  stop("All multi-patient sample counts must be integers >= 1.")
}
max_samples_per_patient <- max(samples_per_patient)

patient_latent <- sim_draw_patient_latents_preview(
  n_patients = n_preview_patients,
  cfg = cfg,
  seed = cfg$seed_base + 17L
)

baseline_scenario <- scenario_tbl %>% dplyr::filter(scenario == "baseline")
if (nrow(baseline_scenario) != 1L) {
  stop("Expected exactly one baseline row in scenario_tbl.")
}

sample_levels <- paste0("Sample ", seq_len(max_samples_per_patient))
patient_levels <- paste0("Patient ", sprintf("%02d", seq_len(n_preview_patients)))

multi_patient_points <- list()
multi_patient_meta <- list()

for (p_idx in seq_len(n_preview_patients)) {
  p_row <- patient_latent[p_idx, , drop = FALSE]
  n_samples_this_patient <- samples_per_patient[[p_idx]]
  for (samp_idx in seq_len(n_samples_this_patient)) {
    seed_use <- as.integer(cfg$seed_base + 700000L + p_idx * 100L + samp_idx)
    sim_out <- sim_eval_preview_with_retry(
      patient_row = p_row,
      scenario_row = baseline_scenario,
      cfg = cfg,
      seed = seed_use,
      contaminated = FALSE
    )

    point_df <- as.data.frame(sim_out$pp)
    if ("marks" %in% colnames(point_df)) names(point_df)[names(point_df) == "marks"] <- "cell_type"
    if (is.list(point_df$cell_type)) point_df$cell_type <- unlist(point_df$cell_type, use.names = FALSE)

    patient_display <- paste0("Patient ", sprintf("%02d", p_idx))
    sample_display <- paste0("Sample ", samp_idx)
    sample_id <- sprintf("%s_baseline_s%02d", p_row$patient[[1]], samp_idx)

    multi_patient_points[[length(multi_patient_points) + 1L]] <- point_df %>%
      dplyr::mutate(
        patient = p_row$patient[[1]],
        patient_display = patient_display,
        sample = sample_id,
        sample_display = sample_display,
        cell_type = factor(as.character(cell_type), levels = c("A", "B"))
      )

    multi_patient_meta[[length(multi_patient_meta) + 1L]] <- data.frame(
      patient = p_row$patient[[1]],
      patient_display = patient_display,
      sample = sample_id,
      sample_display = sample_display,
      quality_z = sim_out$quality_z,
      quality_abs = sim_out$quality_abs,
      parent_lambda = sim_out$parent_lambda,
      offspring_mean = sim_out$offspring_mean,
      offspring_sd = sim_out$offspring_sd,
      n_cells_A = sim_out$n_cells_A,
      n_cells_B = sim_out$n_cells_B,
      n_cells_total = sim_out$n_cells_total,
      prop_A = sim_out$n_cells_A / sim_out$n_cells_total,
      prop_B = sim_out$n_cells_B / sim_out$n_cells_total,
      stringsAsFactors = FALSE
    )
  }
}

if (length(multi_patient_points) < 1L) {
  stop("No multi-sample patient preview point patterns were generated.")
}

multi_patient_points <- dplyr::bind_rows(multi_patient_points) %>%
  dplyr::mutate(
    patient_display = factor(patient_display, levels = patient_levels),
    sample_display = factor(sample_display, levels = sample_levels)
  )

multi_patient_points_plot <- multi_patient_points %>%
  dplyr::group_by(patient_display, sample_display) %>%
  dplyr::group_modify(~ thin_plot_points(
    df = .x,
    keep_frac = analysis_params$multi_plot_keep_frac,
    seed = cfg$seed_base +
      match(as.character(.y$patient_display[[1]]), patient_levels) * 1000L +
      match(as.character(.y$sample_display[[1]]), sample_levels) * 10L
  )) %>%
  dplyr::ungroup()

multi_patient_meta <- dplyr::bind_rows(multi_patient_meta) %>%
  dplyr::mutate(
    patient_display = factor(patient_display, levels = patient_levels),
    sample_display = factor(sample_display, levels = sample_levels)
  ) %>%
  dplyr::arrange(patient_display, sample_display)

multi_patient_plot <- ggplot(multi_patient_points_plot, aes(x = x, y = y, color = cell_type)) +
  geom_point(size = analysis_params$multi_point_size, alpha = analysis_params$point_alpha) +
  scale_color_manual(values = c("A" = "#1F77B4", "B" = "#FF7F0E"), drop = FALSE) +
  coord_equal(
    xlim = c(0, cfg$square_size_um),
    ylim = c(0, cfg$square_size_um),
    expand = FALSE
  ) +
  facet_grid(rows = vars(patient_display), cols = vars(sample_display), drop = FALSE) +
  theme_void(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "#F3F3F3", color = "#C7C7C7"),
    strip.text = element_text(size = 8),
    legend.position = "none",
    panel.spacing = grid::unit(0.08, "in"),
    plot.margin = margin(6, 10, 6, 10)
  )

progressive_points <- list()
progressive_meta <- list()

degradation_patient_index <- min(max(as.integer(analysis_params$degradation_patient_index), 1L), nrow(patient_latent))
degradation_sample_index <- max(as.integer(analysis_params$degradation_sample_index), 1L)
patient_progressive <- patient_latent[degradation_patient_index, , drop = FALSE]

for (s_idx in seq_len(nrow(scenario_use))) {
  scn <- scenario_use[s_idx, , drop = FALSE]
  scn_name <- scn$scenario[[1]]
  scn_label <- if (!is.na(scenario_display_map[[scn_name]])) scenario_display_map[[scn_name]] else scn_name
  seed_use <- as.integer(cfg$seed_base + 900000L + degradation_patient_index * 100L + degradation_sample_index)
  sim_out <- sim_eval_preview_with_retry(
    patient_row = patient_progressive,
    scenario_row = scn,
    cfg = cfg,
    seed = seed_use,
    contaminated = !identical(scn_name, "baseline")
  )

  point_df <- as.data.frame(sim_out$pp)
  if ("marks" %in% colnames(point_df)) names(point_df)[names(point_df) == "marks"] <- "cell_type"
  if (is.list(point_df$cell_type)) point_df$cell_type <- unlist(point_df$cell_type, use.names = FALSE)

  sample_id <- sprintf("%s_progressive_s%02d", patient_progressive$patient[[1]], degradation_sample_index)
  progressive_points[[length(progressive_points) + 1L]] <- point_df %>%
    dplyr::mutate(
      scenario = scn_name,
      scenario_display = scn_label,
      patient = patient_progressive$patient[[1]],
      sample = sample_id,
      cell_type = factor(as.character(cell_type), levels = c("A", "B"))
    )

  progressive_meta[[length(progressive_meta) + 1L]] <- data.frame(
    scenario = scn_name,
    scenario_display = scn_label,
    patient = patient_progressive$patient[[1]],
    sample = sample_id,
    contaminated = sim_out$contaminated,
    quality_z = sim_out$quality_z,
    quality_abs = sim_out$quality_abs,
    parent_lambda = sim_out$parent_lambda,
    offspring_mean = sim_out$offspring_mean,
    offspring_sd = sim_out$offspring_sd,
    n_cells_A = sim_out$n_cells_A,
    n_cells_B = sim_out$n_cells_B,
    n_cells_total = sim_out$n_cells_total,
    prop_A = sim_out$n_cells_A / sim_out$n_cells_total,
    prop_B = sim_out$n_cells_B / sim_out$n_cells_total,
    stringsAsFactors = FALSE
  )
}

if (length(progressive_points) < 1L) {
  stop("No progressive degradation point-pattern examples were generated.")
}

progressive_points <- dplyr::bind_rows(progressive_points) %>%
  dplyr::mutate(
    scenario_display = factor(scenario_display, levels = scenario_display_levels)
  )

progressive_points_plot <- progressive_points %>%
  dplyr::group_by(scenario_display) %>%
  dplyr::group_modify(~ thin_plot_points(
    df = .x,
    keep_frac = analysis_params$progressive_plot_keep_frac,
    seed = cfg$seed_base + match(as.character(.y$scenario_display[[1]]), scenario_display_levels) * 1000L
  )) %>%
  dplyr::ungroup()

progressive_meta <- dplyr::bind_rows(progressive_meta) %>%
  dplyr::mutate(
    scenario_display = factor(scenario_display, levels = scenario_display_levels)
  ) %>%
  dplyr::arrange(scenario_display)

progressive_plot_list <- lapply(scenario_display_levels, function(scn_label) {
  df_sub <- progressive_points_plot %>%
    dplyr::filter(as.character(scenario_display) == scn_label)

  ggplot(df_sub, aes(x = x, y = y, color = cell_type)) +
    geom_point(size = analysis_params$progressive_point_size, alpha = analysis_params$point_alpha) +
    scale_color_manual(values = c("A" = "#1F77B4", "B" = "#FF7F0E"), drop = FALSE) +
    coord_equal(
      xlim = c(0, cfg$square_size_um),
      ylim = c(0, cfg$square_size_um),
      expand = FALSE
    ) +
    labs(title = scn_label) +
    theme_void(base_size = 8) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.margin = margin(2, 2, 2, 2),
      panel.background = element_rect(fill = "white", color = "#D0D0D0", linewidth = 0.2),
      plot.background = element_rect(fill = "white", color = NA)
    )
})

progressive_plot <- cowplot::plot_grid(
  plotlist = progressive_plot_list,
  ncol = 1,
  align = "v"
)

perf_use <- perf_summary %>%
  dplyr::filter(
    parameter %in% c("mu", "tau2_patient"),
    method %in% c("panoramic", "naive"),
    scenario %in% scenario_levels
  ) %>%
  dplyr::mutate(
    scenario_display = factor(
      ifelse(
        scenario %in% names(scenario_pretty),
        scenario_pretty[scenario],
        scenario
      ),
      levels = unname(scenario_pretty[scenario_levels])
    ),
    parameter_display = factor(
      ifelse(parameter == "mu", "mu", "tau2"),
      levels = c("mu", "tau2")
    ),
    method = factor(method, levels = c("panoramic", "naive"))
  ) %>%
  dplyr::arrange(parameter_display, scenario_display, n_patients, method)

if (nrow(perf_use) < 1L) {
  stop("No rows available for RMSE panel after filtering to mu/tau2 and PANORAMIC/naive.")
}

perf_use <- perf_use %>%
  dplyr::mutate(
    method_display = factor(
      ifelse(method == "panoramic", "PANORAMIC", "Naive"),
      levels = c("PANORAMIC", "Naive")
    )
  )

has_rmse_ci <- all(c("rmse_ci_lower", "rmse_ci_upper") %in% colnames(perf_use)) &&
  any(is.finite(perf_use$rmse_ci_lower)) &&
  any(is.finite(perf_use$rmse_ci_upper))

rmse_plot <- ggplot(
  perf_use,
  aes(x = n_patients, y = rmse, color = method_display, group = method_display)
) 

if (has_rmse_ci) {
  rmse_plot <- rmse_plot +
    geom_ribbon(
      aes(ymin = rmse_ci_lower, ymax = rmse_ci_upper, fill = method_display),
      alpha = 0.18,
      linewidth = 0,
      color = NA,
      show.legend = FALSE
    )
}

rmse_plot <- rmse_plot +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.9) +
  facet_grid(parameter_display ~ scenario_display, scales = "free_y") +
  scale_color_manual(
    values = c("PANORAMIC" = "#0072B2", "Naive" = "#D55E00"),
    breaks = c("PANORAMIC", "Naive")
  ) +
  scale_fill_manual(
    values = c("PANORAMIC" = "#0072B2", "Naive" = "#D55E00"),
    breaks = c("PANORAMIC", "Naive")
  ) +
  scale_x_continuous(breaks = sort(unique(perf_use$n_patients))) +
  labs(
    x = "Number of patients",
    y = "RMSE",
    color = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "#F3F3F3", color = "#C7C7C7"),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 8)
  )

avg_cmp_long <- dplyr::bind_rows(
  breakpoint_estimates %>%
    dplyr::transmute(
      scenario, n_patients, replicate, feature_id,
      parameter = "mu", method = "panoramic",
      estimate = mu_hat, truth = true_mu
    ),
  breakpoint_estimates %>%
    dplyr::transmute(
      scenario, n_patients, replicate, feature_id,
      parameter = "mu", method = "naive",
      estimate = mu_naive, truth = true_mu
    ),
  breakpoint_estimates %>%
    dplyr::transmute(
      scenario, n_patients, replicate, feature_id,
      parameter = "tau2_patient", method = "panoramic",
      estimate = tau2_patient, truth = true_tau2_patient
    ),
  breakpoint_estimates %>%
    dplyr::transmute(
      scenario, n_patients, replicate, feature_id,
      parameter = "tau2_patient", method = "naive",
      estimate = tau2_patient_naive, truth = true_tau2_patient
    ),
  breakpoint_estimates %>%
    dplyr::transmute(
      scenario, n_patients, replicate, feature_id,
      parameter = "tau2_sample", method = "panoramic",
      estimate = tau2_sample, truth = true_tau2_sample
    ),
  breakpoint_estimates %>%
    dplyr::transmute(
      scenario, n_patients, replicate, feature_id,
      parameter = "tau2_sample", method = "naive",
      estimate = tau2_sample_naive, truth = true_tau2_sample
    )
) %>%
  dplyr::mutate(
    error = estimate - truth,
    abs_error = abs(error)
  )

avg_perf_by_replicate <- avg_cmp_long %>%
  dplyr::group_by(scenario, n_patients, parameter, method, replicate) %>%
  dplyr::summarise(
    bias = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mae = mean(abs_error, na.rm = TRUE),
    n_features = dplyr::n_distinct(feature_id),
    .groups = "drop"
  )

avg_perf_use <- avg_perf_by_replicate %>%
  dplyr::group_by(scenario, n_patients, parameter, method) %>%
  dplyr::summarise(
    bias = mean(bias, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
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
      pmax(0, rmse - rmse_ci_mult * rmse_se),
      rmse
    ),
    rmse_ci_upper = dplyr::if_else(
      n_reps > 1L,
      rmse + rmse_ci_mult * rmse_se,
      rmse
    )
  ) %>%
  dplyr::filter(
    parameter %in% c("mu", "tau2_patient", "tau2_sample"),
    method %in% c("panoramic", "naive"),
    scenario %in% scenario_levels
  ) %>%
  dplyr::mutate(
    scenario_display = factor(
      ifelse(
        scenario %in% names(scenario_pretty),
        scenario_pretty[scenario],
        scenario
      ),
      levels = unname(scenario_pretty[scenario_levels])
    ),
    parameter_display = factor(
      dplyr::case_when(
        parameter == "mu" ~ "Mean",
        parameter == "tau2_patient" ~ "Patient tau2",
        parameter == "tau2_sample" ~ "Sample tau2",
        TRUE ~ parameter
      ),
      levels = c("Mean", "Patient tau2", "Sample tau2")
    ),
    method_display = factor(
      ifelse(method == "panoramic", "PANORAMIC", "Naive"),
      levels = c("PANORAMIC", "Naive")
    )
  ) %>%
  dplyr::arrange(parameter_display, scenario_display, n_patients, method_display)

avg_rmse_plot <- ggplot(
  avg_perf_use,
  aes(x = n_patients, y = rmse, color = method_display, group = method_display)
)

if (all(c("rmse_ci_lower", "rmse_ci_upper") %in% colnames(avg_perf_use)) &&
    any(is.finite(avg_perf_use$rmse_ci_lower)) &&
    any(is.finite(avg_perf_use$rmse_ci_upper))) {
  avg_rmse_plot <- avg_rmse_plot +
    geom_ribbon(
      aes(ymin = rmse_ci_lower, ymax = rmse_ci_upper, fill = method_display),
      alpha = 0.18,
      linewidth = 0,
      color = NA,
      show.legend = FALSE
    )
}

avg_rmse_plot <- avg_rmse_plot +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.9) +
  facet_grid(parameter_display ~ scenario_display, scales = "free_y") +
  scale_color_manual(
    values = c("PANORAMIC" = "#0072B2", "Naive" = "#D55E00"),
    breaks = c("PANORAMIC", "Naive")
  ) +
  scale_fill_manual(
    values = c("PANORAMIC" = "#0072B2", "Naive" = "#D55E00"),
    breaks = c("PANORAMIC", "Naive")
  ) +
  scale_x_continuous(breaks = sort(unique(avg_perf_use$n_patients))) +
  labs(
    x = "Number of patients",
    y = "Average RMSE across cell pairs",
    color = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "#F3F3F3", color = "#C7C7C7"),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 8)
  )

print(multi_patient_plot)
print(progressive_plot)
print(rmse_plot)
print(avg_rmse_plot)

readr::write_csv(
  multi_patient_meta,
  file.path(out_tbl_dir, "intersample_publication_multi_patient_sample_metadata.csv")
)
readr::write_csv(
  progressive_meta,
  file.path(out_tbl_dir, "intersample_publication_progressive_degradation_metadata.csv")
)
readr::write_csv(
  perf_use %>%
    dplyr::select(
      dplyr::any_of(c(
        "scenario", "n_patients", "parameter", "method",
        "rmse", "rmse_ci_lower", "rmse_ci_upper", "n_reps",
        "bias", "mae", "n_features"
      ))
    ),
  file.path(out_tbl_dir, "intersample_publication_rmse_mu_tau2_data.csv")
)
readr::write_csv(
  avg_perf_use %>%
    dplyr::select(
      dplyr::any_of(c(
        "scenario", "n_patients", "parameter", "method",
        "rmse", "rmse_ci_lower", "rmse_ci_upper", "n_reps",
        "bias", "mae", "n_features"
      ))
    ),
  file.path(out_tbl_dir, "intersample_publication_rmse_mean_tau2_patient_tau2_sample_avg_data.csv")
)

save_plot_formats(
  plot_obj = multi_patient_plot,
  stem = file.path(out_fig_dir, "intersample_publication_multi_patient_samples"),
  width = 4.75,
  height = 3
)
save_plot_formats(
  plot_obj = rmse_plot,
  stem = file.path(out_fig_dir, "intersample_publication_rmse_mu_tau2_degradation_curves"),
  width = 6.5,
  height = 3
)
save_plot_formats(
  plot_obj = avg_rmse_plot,
  stem = file.path(out_fig_dir, "intersample_publication_rmse_mean_tau2_patient_tau2_sample_avg_ci"),
  width = 6.5,
  height = 4.5
)
save_plot_formats(
  plot_obj = progressive_plot,
  stem = file.path(out_fig_dir, "intersample_publication_progressive_degradation_single_patient"),
  width = 1.75,
  height = 3
)

message("Saved intersample publication figures and tables.")
message("  Figure dir: ", out_fig_dir)
message("  Table dir: ", out_tbl_dir)
