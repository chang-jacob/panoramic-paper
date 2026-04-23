# TMA Analysis Helper Script: Shared Runtime/Table Utilities
# Author: Jacob Chang
# Date: 2026-04-23
# Summary:
# - Shared helper functions for CRC/HNSCC TMA run scripts.
# - Provides runtime tracking, metadata attachment, result extraction, and LaTeX table rendering.

collapse_kv <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  nms <- names(x)
  if (is.null(nms)) nms <- rep("", length(x))
  parts <- vapply(
    seq_along(x),
    function(i) {
      key <- nms[[i]]
      val <- as.character(x[[i]])
      if (!nzchar(key)) val else paste0(key, "=", val)
    },
    FUN.VALUE = character(1)
  )
  paste(parts, collapse = "; ")
}

create_runtime_tracker <- function(run_prefix) {
  runtime_records <- list()

  summarize_se_dims <- function(obj) {
    if (inherits(obj, "SummarizedExperiment")) {
      return(list(
        n_features = as.integer(nrow(obj)),
        n_samples = as.integer(ncol(obj)),
        n_assays = as.integer(length(SummarizedExperiment::assays(obj)))
      ))
    }
    list(n_features = NA_integer_, n_samples = NA_integer_, n_assays = NA_integer_)
  }

  append_runtime_metric <- function(stage, status, started_at, ended_at,
                                    cache_path = NA_character_,
                                    params = NULL,
                                    obj = NULL) {
    dims <- summarize_se_dims(obj)
    runtime_records[[length(runtime_records) + 1L]] <<- data.frame(
      run_prefix = run_prefix,
      stage = as.character(stage),
      status = as.character(status),
      started_at = format(started_at, tz = "UTC", usetz = TRUE),
      ended_at = format(ended_at, tz = "UTC", usetz = TRUE),
      elapsed_sec = round(as.numeric(difftime(ended_at, started_at, units = "secs")), 3),
      cache_path = as.character(cache_path),
      params = collapse_kv(params),
      n_features = dims$n_features,
      n_samples = dims$n_samples,
      n_assays = dims$n_assays,
      stringsAsFactors = FALSE
    )
  }

  write_runtime_metrics <- function(path) {
    if (length(runtime_records) == 0L) return(invisible(NULL))
    runtime_df <- dplyr::bind_rows(runtime_records)
    readr::write_csv(runtime_df, path)
    message("Runtime metrics CSV: ", path)
    runtime_df
  }

  list(
    append_runtime_metric = append_runtime_metric,
    write_runtime_metrics = write_runtime_metrics
  )
}

get_unique <- function(spe, col) {
  vals <- unique(as.character(SummarizedExperiment::colData(spe)[[col]]))
  vals <- vals[!is.na(vals)]
  if (length(vals) != 1L) {
    stop(
      "Expected exactly one unique value for colData[[", col, "]] per sample. Got: ",
      paste(vals, collapse = ", ")
    )
  }
  vals
}

attach_sample_metadata <- function(
    se,
    sample_meta_df,
    patient_col = "patient",
    group_col = "group",
    region_col = "region",
    sample_col = "sample",
    sample_meta_sample_col = "sample",
    sample_meta_patient_col = NULL,
    sample_meta_group_col = NULL
) {
  if (is.null(sample_meta_patient_col)) {
    sample_meta_patient_col <- if (patient_col %in% names(sample_meta_df)) patient_col else "patient"
  }
  if (is.null(sample_meta_group_col)) {
    sample_meta_group_col <- if (group_col %in% names(sample_meta_df)) group_col else "group"
  }
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  idx <- match(cd[[sample_col]], sample_meta_df[[sample_meta_sample_col]])
  if (!is.null(region_col) && nzchar(region_col)) {
    cd[[region_col]] <- cd[[sample_col]]
  }
  cd[[patient_col]] <- sample_meta_df[[sample_meta_patient_col]][idx]
  cd[[group_col]] <- sample_meta_df[[sample_meta_group_col]][idx]
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd)
  se
}

add_coloc_direction <- function(df) {
  if (!is.data.frame(df)) {
    stop("Expected a data.frame in add_coloc_direction().")
  }
  if (!all(c("ct1", "ct2") %in% colnames(df))) {
    stop("add_coloc_direction() requires ct1/ct2 columns.")
  }
  df$coloc_source <- as.character(df$ct2)
  df$coloc_target <- as.character(df$ct1)
  df$coloc_direction <- paste0(df$coloc_source, " -> ", df$coloc_target)
  df
}

extract_spatialstats_table <- function(se) {
  yi <- SummarizedExperiment::assay(se, "yi")
  vi <- SummarizedExperiment::assay(se, "vi")
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  rid <- rep(seq_len(nrow(yi)), times = ncol(yi))
  cid <- rep(seq_len(ncol(yi)), each = nrow(yi))
  out <- cbind(
    rd[rid, , drop = FALSE],
    cd[cid, , drop = FALSE],
    yi = as.numeric(yi),
    vi = as.numeric(vi)
  )
  out <- add_coloc_direction(out)
  rownames(out) <- NULL
  out
}

extract_meta_table <- function(se) {
  out <- as.data.frame(SummarizedExperiment::rowData(se))
  out <- add_coloc_direction(out)
  rownames(out) <- NULL
  out
}

extract_contrast_table <- function(se) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  contrast_cols <- c("beta_diff", "se_diff", "z_diff", "p_diff", "fdr_diff")
  missing_cols <- setdiff(contrast_cols, colnames(rd))
  if (length(missing_cols) > 0L) {
    stop("Missing contrast columns in rowData(se): ", paste(missing_cols, collapse = ", "))
  }
  keep <- intersect(c("ct1", "ct2", "radius_um", "stat"), colnames(rd))
  out <- rd[, c(keep, contrast_cols), drop = FALSE]
  out <- add_coloc_direction(out)
  rownames(out) <- NULL
  out
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#$])", "\\\\\\1", x, perl = TRUE)
  x
}

format_group_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x <- gsub("met adjacent", "met-adjacent", x, ignore.case = TRUE)
  x
}

format_cell_pair_latex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s*->\\s*", " XXARROWTOKENXX ", x)
  x <- latex_escape(x)
  x <- gsub("XXARROWTOKENXX", "$\\\\to$", x, fixed = TRUE)
  x
}

build_panoramic_results_table <- function(meta_tbl, sig_only = TRUE, alpha = 0.05) {
  if (!is.data.frame(meta_tbl)) {
    stop("meta_tbl must be a data.frame.")
  }
  if (nrow(meta_tbl) == 0L) {
    stop("meta_tbl has 0 rows.")
  }

  mu_cols <- grep("_mu_hat$", names(meta_tbl), value = TRUE)
  if (length(mu_cols) != 2L) {
    stop(
      "Expected exactly two '*_mu_hat' columns in meta table. Found: ",
      paste(mu_cols, collapse = ", ")
    )
  }

  required_cols <- c("coloc_direction", "beta_diff", "p_diff", "fdr_diff")
  missing_cols <- setdiff(required_cols, names(meta_tbl))
  if (length(missing_cols) > 0L) {
    stop("Meta table missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  group_labels <- gsub("_mu_hat$", "", mu_cols)
  group_labels <- format_group_label(group_labels)
  out <- data.frame(
    "Cell type pair" = as.character(meta_tbl$coloc_direction),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out[[group_labels[[1]]]] <- as.numeric(meta_tbl[[mu_cols[[1]]]])
  out[[group_labels[[2]]]] <- as.numeric(meta_tbl[[mu_cols[[2]]]])
  out[["Effect size"]] <- as.numeric(meta_tbl$beta_diff)
  out[["p-value"]] <- as.numeric(meta_tbl$p_diff)
  out[["adjusted p-value"]] <- as.numeric(meta_tbl$fdr_diff)

  ord <- order(out[["adjusted p-value"]], out[["p-value"]], na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  if (isTRUE(sig_only)) {
    keep <- is.finite(out[["adjusted p-value"]]) & out[["adjusted p-value"]] <= alpha
    out <- out[keep, , drop = FALSE]
    rownames(out) <- NULL
  }

  out
}

render_panoramic_results_latex <- function(
    tbl,
    caption = "PANORAMIC colocalization results.",
    label = "tab:s_panoramic_results",
    first_col_width = "3.2cm",
    continued_text = "Continued on next page"
) {
  col_names <- names(tbl)
  if (length(col_names) != 6L) {
    stop("Expected 6 columns in results table; got ", length(col_names), ".")
  }

  col_align <- paste0(">{\\raggedright\\arraybackslash}p{", first_col_width, "}rrrrr")
  header <- paste(latex_escape(col_names), collapse = " & ")

  fmt_est <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 4), "NA")
  fmt_p <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 3), "NA")

  if (nrow(tbl) == 0L) {
    row_lines <- "No significant rows & NA & NA & NA & NA & NA\\\\"
  } else {
    row_lines <- vapply(seq_len(nrow(tbl)), function(i) {
      vals <- c(
        format_cell_pair_latex(tbl[[1]][i]),
        fmt_est(tbl[[2]][i]),
        fmt_est(tbl[[3]][i]),
        fmt_est(tbl[[4]][i]),
        fmt_p(tbl[[5]][i]),
        fmt_p(tbl[[6]][i])
      )
      paste0(paste(vals, collapse = " & "), "\\\\")
    }, FUN.VALUE = character(1))
  }

  caption_line <- paste0("\\caption{", latex_escape(caption), "}")
  label_line <- if (!is.null(label) && nzchar(as.character(label))) {
    paste0("\\label{", as.character(label), "}\\\\")
  } else {
    "\\\\"
  }

  lines <- c(
    paste0("\\begin{longtable}[t]{", col_align, "}"),
    caption_line,
    label_line,
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
    paste0("\\multicolumn{6}{r}{", latex_escape(continued_text), "}\\\\"),
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

write_panoramic_latex_results <- function(
    meta_tbl,
    out_path,
    sig_only = TRUE,
    alpha = 0.05,
    caption = "PANORAMIC colocalization results.",
    label = "tab:s_panoramic_results",
    first_col_width = "3.2cm",
    print_console = FALSE
) {
  tbl <- build_panoramic_results_table(meta_tbl = meta_tbl, sig_only = sig_only, alpha = alpha)
  latex_txt <- render_panoramic_results_latex(
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
