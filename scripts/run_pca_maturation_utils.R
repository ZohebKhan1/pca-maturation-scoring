# PCA Maturation Scoring: Internal Utils
# sourced by run_pca_maturation.R -- don't call these directly
#
# ================================================================
# Table of Contents
# ================================================================
#
# 1. Logging -- progress msgs and file output tracking
#    1.1 log_step
#    1.2 log_output
#
# 2. Validation -- input checks before analysis starts
#    2.1 check_required_packages
#    2.2 check_required_columns
#
# 3. Metadata -- row reordering to align w/ dds col order
#    3.1 reorder_by_sample_id
#
# 4. PCA and Scoring -- VST prep, projection, centroids, vector scoring
#    4.1 compute_vst_matrix
#    4.2 project_into_pca_space
#    4.3 compute_timepoint_centroids
#    4.4 compute_best_fit_direction
#    4.5 project_onto_scoring_vector
#
# 5. Colors and Layout -- auto palette, boxplot x-offsets
#    5.1 assign_group_colors
#    5.2 compute_group_x_offsets
#
# 6. Welch t-tests -- per-tp group comparisons on maturation scores
#    6.1 validate_t_test_pairs
#    6.2 compute_group_stats
#    6.3 run_welch_t_tests
#
# 7. Boxplot -- score range bars, p-val brackets, SVG export
#    7.0 format_axis_labels
#    7.1 compute_score_ranges
#    7.2 format_pvalue_label
#    7.3 compute_bracket_positions
#    7.4 build_and_save_boxplot
#
# 8. XLSX export -- plain publication tables for maturation scores
#    8.1 xlsx helpers
#    8.2 summarize_maturation_scores
#    8.3 save_maturation_score_xlsx


# ================================================================
# 1. Logging
# ================================================================

# 1.1 numbered progress line (e.g. "[2/5] running LRT")
log_step <- function(step_number, step_total, text) {
  message(sprintf('[%d/%d] %s', step_number, step_total, text))
}

# 1.2 report that a file was saved
log_output <- function(label, path) {
  message(paste0('  saved ', label, ' -> ', normalizePath(path, mustWork = FALSE)))
}


# ================================================================
# 2. Validation
# ================================================================

# 2.1 bail early if required pkgs aren't installed
check_required_packages <- function() {
  required <- c('DESeq2', 'SummarizedExperiment', 'ggplot2', 'svglite')
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) stop('missing required packages: ', paste(missing, collapse = ', '))
}

# 2.2 Stop if a df is missing expected cols
check_required_columns <- function(df, columns, df_name) {
  missing_cols <- setdiff(columns, colnames(df))
  if (length(missing_cols) > 0L) stop(df_name, ' is missing required columns: ', paste(missing_cols, collapse = ', '))
}


# ================================================================
# 3. Metadata
# ================================================================

# 3.1 reorder metadata rows to match a given vector of sample IDs.
# keeps metadata aligned w/ dds col order for DESeqDataSetFromMatrix.
reorder_by_sample_id <- function(meta, id_col, ids) {
  idx <- match(ids, meta[[id_col]])
  if (anyNA(idx)) stop('metadata is missing one or more requested sample IDs.')
  meta[idx, , drop = FALSE]
}


# ================================================================
# 4. PCA and Scoring
# ================================================================

# 4.1 Compute VST and optionally remove nuisance batch structure
# while protecting the supplied biological condition design.
compute_vst_matrix <- function(dds, metadata, sample_id_col_name,
                               condition_col_name,
                               batch_col_name = NULL,
                               batch2_col_name = NULL) {
  vst_mat <- SummarizedExperiment::assay(
    DESeq2::vst(dds, blind = FALSE)
  )

  if (is.null(batch_col_name) && is.null(batch2_col_name)) {
    return(list(
      matrix = vst_mat,
      applied_batch_cols = character(0),
      dropped_batch_cols = character(0)
    ))
  }

  if (!requireNamespace('limma', quietly = TRUE)) {
    stop('limma is required for optional VST batch correction.')
  }

  meta_aligned <- reorder_by_sample_id(
    metadata, sample_id_col_name, colnames(vst_mat)
  )
  protected_condition <- factor(
    as.character(meta_aligned[[condition_col_name]])
  )
  protected_design <- stats::model.matrix(
    ~protected_condition
  )

  batch <- NULL
  batch2 <- NULL
  applied_batch_cols <- character(0)
  dropped_batch_cols <- character(0)

  if (!is.null(batch_col_name)) {
    batch <- factor(as.character(meta_aligned[[batch_col_name]]))
    if (length(unique(stats::na.omit(batch))) < 2L) {
      batch <- NULL
      dropped_batch_cols <- c(dropped_batch_cols, batch_col_name)
    } else {
      applied_batch_cols <- c(applied_batch_cols, batch_col_name)
    }
  }

  if (!is.null(batch2_col_name)) {
    batch2 <- factor(as.character(meta_aligned[[batch2_col_name]]))
    if (length(unique(stats::na.omit(batch2))) < 2L) {
      batch2 <- NULL
      dropped_batch_cols <- c(dropped_batch_cols, batch2_col_name)
    } else {
      applied_batch_cols <- c(applied_batch_cols, batch2_col_name)
    }
  }

  if (is.null(batch) && is.null(batch2)) {
    return(list(
      matrix = vst_mat,
      applied_batch_cols = character(0),
      dropped_batch_cols = dropped_batch_cols
    ))
  }

  corrected <- limma::removeBatchEffect(
    vst_mat,
    batch = batch,
    batch2 = batch2,
    design = protected_design
  )

  list(
    matrix = corrected,
    applied_batch_cols = applied_batch_cols,
    dropped_batch_cols = dropped_batch_cols
  )
}

# 4.2 Project samples into an existing PCA space using the same
# centering/scaling parameters learned by prcomp().
project_into_pca_space <- function(vst_subset, pca_fit, n_pcs = 3L) {
  n_pcs <- as.integer(n_pcs[[1]])
  pc_cols <- paste0('PC', seq_len(n_pcs))
  rotation <- pca_fit[['rotation']]
  if (!identical(rownames(vst_subset), rownames(rotation))) {
    stop('vst_subset rows must match pca_fit rotation rows exactly.')
  }

  pca_scale <- pca_fit[['scale']]
  if (is.logical(pca_scale) && identical(pca_scale, FALSE)) {
    pca_scale <- FALSE
  } else if (is.numeric(pca_scale)) {
    pca_scale <- pca_scale[rownames(rotation)]
  }

  pca_center <- pca_fit[['center']]
  if (is.logical(pca_center) && identical(pca_center, FALSE)) {
    pca_center <- FALSE
  } else if (is.numeric(pca_center)) {
    pca_center <- pca_center[rownames(rotation)]
  }

  projected <- scale(
    t(vst_subset),
    center = pca_center,
    scale = pca_scale
  )
  coord_df <- as.data.frame(projected %*% rotation)
  for (pc in pc_cols) {
    if (!pc %in% colnames(coord_df)) coord_df[[pc]] <- 0
  }
  coord_df[, pc_cols, drop = FALSE]
}

# 4.3 mean PCA position for each ref tp across the requested scoring PCs.
# these define the scoring vector endpoints.
compute_timepoint_centroids <- function(score_df, timepoints, pc_cols = c('PC1', 'PC2', 'PC3')) {
  rows <- lapply(timepoints, function(tp) {
    sub <- score_df[score_df[['sample_set']] == 'reference' & score_df[['timepoint_numeric']] == tp, , drop = FALSE]
    if (nrow(sub) == 0L) stop('centroid could not be calculated for timepoint ', tp, '. no reference samples found.')
    row <- data.frame(
      timepoint_numeric = tp,
      stringsAsFactors = FALSE
    )
    for (pc in pc_cols) {
      row[[pc]] <- mean(sub[[pc]])
    }
    row
  })
  do.call(rbind, rows)
}

# 4.4 direction of best fit through ref samples in PCA space.
# fits a line via PCA on the ref PCA coords; PC1 gives the direction
# of greatest variance (maturation axis). oriented start -> end.
compute_best_fit_direction <- function(ref_coords, start_pt, end_pt, pc_cols = c('PC1', 'PC2', 'PC3')) {
  X <- as.matrix(ref_coords[, pc_cols])
  mu <- colMeans(X)
  Xc <- sweep(X, 2, mu)
  pca_line <- stats::prcomp(Xc, center = FALSE, scale. = FALSE)
  v <- pca_line[['rotation']][, 1]

  # orient so direction aligns w/ start -> end
  centroid_vec <- end_pt - start_pt
  if (sum(v * centroid_vec) < 0) v <- -v

  list(direction = v, origin = mu)
}

# 4.5 project samples onto scoring axis, compute normalized scores.
# two modes:
#   centroid_vector -- axis from start centroid to end centroid
#   best_fit_line   -- PC1 of ref samples in PCA space,
#                      normalized so start = 0, end = 1
project_onto_scoring_vector <- function(score_df, centroid_df, start_tp,
                                        end_tp,
                                        scoring_mode = 'centroid_vector',
                                        pc_cols = c('PC1', 'PC2', 'PC3')) {
  scoring_mode <- match.arg(scoring_mode, c('centroid_vector', 'best_fit_line'))

  start_row <- centroid_df[
    centroid_df[['timepoint_numeric']] == start_tp, pc_cols,
    drop = FALSE
  ]
  end_row <- centroid_df[
    centroid_df[['timepoint_numeric']] == end_tp, pc_cols,
    drop = FALSE
  ]
  start_pt <- as.numeric(start_row[1, ])
  end_pt <- as.numeric(end_row[1, ])
  if (anyNA(start_pt) || anyNA(end_pt)) {
    stop('could not find both scoring vector endpoints in the centroids.')
  }

  if (scoring_mode == 'centroid_vector') {
    vec_raw <- end_pt - start_pt
    vec_len_eff <- sqrt(sum(vec_raw^2))
    if (vec_len_eff == 0) {
      stop('scoring vector has zero length (identical centroids).')
    }
    vec_unit <- vec_raw / vec_len_eff
    origin <- start_pt
  } else {
    # best_fit_line: PC1 of ref samples, oriented start -> end
    is_ref <- score_df[['sample_set']] == 'reference'
    ref_rows <- score_df[is_ref, , drop = FALSE]
    bf <- compute_best_fit_direction(ref_rows, start_pt, end_pt, pc_cols = pc_cols)
    vec_unit <- bf[['direction']]
    origin <- bf[['origin']]
    vec_len_eff <- NULL
  }

  coords <- as.matrix(score_df[, pc_cols])
  offsets <- sweep(coords, 2, origin)
  raw_proj <- as.numeric(offsets %*% vec_unit)
  proj_coords <- matrix(
    origin,
    nrow = nrow(coords), ncol = length(pc_cols), byrow = TRUE
  ) + raw_proj %o% vec_unit

  if (scoring_mode == 'centroid_vector') {
    distance <- raw_proj
  } else {
    # shift so start centroid = 0; normalize by centroid separation
    start_proj <- sum((start_pt - origin) * vec_unit)
    end_proj <- sum((end_pt - origin) * vec_unit)
    vec_len_eff <- end_proj - start_proj
    if (abs(vec_len_eff) < 1e-12) {
      stop('centroids project to the same point on the best-fit line.')
    }
    distance <- raw_proj - start_proj
  }

  list(
    start_point = start_pt, end_point = end_pt,
    vector_unit = vec_unit, vector_length = abs(vec_len_eff),
    distance_along_vector = distance,
    maturation_score = distance / vec_len_eff,
    projected_coords = proj_coords
  )
}


# ================================================================
# 5. Colors and Layout
# ================================================================

# 5.1 resolve user-supplied colors or auto-generate a
# colorblind-friendly fill palette for groups.
assign_group_colors <- function(group_levels,
                                custom_colors = NULL) {
  if (!is.null(custom_colors)) {
    missing <- setdiff(group_levels, names(custom_colors))
    if (length(missing) > 0L) {
      stop(
        'group_colors is missing entries for: ',
        paste(missing, collapse = ', ')
      )
    }
    return(custom_colors[group_levels])
  }

  n <- length(group_levels)
  stats::setNames(
    grDevices::hcl.colors(max(n, 3), palette = 'Dark 3')[seq_len(n)],
    group_levels
  )
}

# 5.2 small x-axis offsets so multiple groups at the same tp
# don't overlap. evenly spaced, scales w/ group count.
compute_group_x_offsets <- function(group_levels,
                                    custom_offsets = NULL) {
  if (!is.null(custom_offsets)) {
    missing <- setdiff(group_levels, names(custom_offsets))
    if (length(missing) > 0L) {
      stop(
        'group_x_offsets is missing entries for: ',
        paste(missing, collapse = ', ')
      )
    }
    return(custom_offsets[group_levels])
  }
  n <- length(group_levels)
  if (n == 1L) return(stats::setNames(0, group_levels))
  max_off <- min(0.20, 0.45 / n)
  stats::setNames(
    seq(-max_off, max_off, length.out = n),
    group_levels
  )
}


# ================================================================
# 6. Welch t-tests
# ================================================================

# 6.1 Validate each t-test pair: must be a length-2 char vec of group
# names present in the scored data. deduplicates pairs.
validate_t_test_pairs <- function(pairs, available_groups) {
  if (!is.list(pairs) || length(pairs) == 0L) {
    stop('welch_t_test_comparison must be a non-empty list of length-2 character vectors.')
  }
  pair_list <- lapply(seq_along(pairs), function(i) {
    p <- as.character(pairs[[i]])
    p <- p[!is.na(p) & nzchar(p)]
    if (length(p) != 2L) stop('each t-test pair must contain exactly two group values.')
    if (!all(p %in% available_groups)) {
      stop('t-test group not in data: ', paste(setdiff(p, available_groups), collapse = ', '))
    }
    p
  })
  keys <- vapply(pair_list, paste, collapse = '__', character(1))
  pair_list[!duplicated(keys)]
}

# 6.2 mean, min, max for a numeric vec. returns all NA if n = 0.
compute_group_stats <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0L) return(list(mean = NA_real_, min = NA_real_, max = NA_real_))
  list(mean = mean(x), min = min(x), max = max(x))
}

# 6.3 Welch t-test at each shared tp for each group pair.
# returns df w/ one row per tp-comparison: group IDs, sample sizes,
# means, min/max, test stat, p-val, status.
run_welch_t_tests <- function(score_df, pairs) {
  make_row <- function(g1, g2, tp_num, n1, n2, s1, s2,
                       est_diff, diff_lo, diff_hi, stat, pval,
                       status, warn = NA_character_) {
    data.frame(
      group_1 = g1, group_2 = g2,
      comparison = paste(g1, 'vs', g2),
      timepoint_numeric = tp_num,
      n_group_1 = n1, n_group_2 = n2,
      mean_group_1 = s1[['mean']],
      min_group_1 = s1[['min']],
      max_group_1 = s1[['max']],
      mean_group_2 = s2[['mean']],
      min_group_2 = s2[['min']],
      max_group_2 = s2[['max']],
      estimate_difference = est_diff,
      difference_ci_low = diff_lo,
      difference_ci_high = diff_hi,
      statistic = stat, p_value = pval,
      status = status, warning_message = warn,
      stringsAsFactors = FALSE
    )
  }

  result_list <- vector('list', length(pairs))
  for (i in seq_along(pairs)) {
    g1 <- pairs[[i]][[1]]
    g2 <- pairs[[i]][[2]]
    g1_df <- score_df[score_df[['group_value']] == g1, , drop = FALSE]
    g2_df <- score_df[score_df[['group_value']] == g2, , drop = FALSE]
    shared_tp <- sort(intersect(
      unique(g1_df[['timepoint_numeric']]),
      unique(g2_df[['timepoint_numeric']])
    ))

    if (length(shared_tp) == 0L) {
      s1 <- compute_group_stats(g1_df[['maturation_score']])
      s2 <- compute_group_stats(g2_df[['maturation_score']])
      result_list[[i]] <- make_row(
        g1, g2, NA_real_, nrow(g1_df), nrow(g2_df), s1, s2,
        s1[['mean']] - s2[['mean']],
        NA_real_, NA_real_, NA_real_, NA_real_,
        'no_shared_timepoints'
      )
      next
    }

    tp_results <- lapply(shared_tp, function(tp) {
      x <- g1_df[g1_df[['timepoint_numeric']] == tp, 'maturation_score', drop = TRUE]
      y <- g2_df[g2_df[['timepoint_numeric']] == tp, 'maturation_score', drop = TRUE]
      s1 <- compute_group_stats(x)
      s2 <- compute_group_stats(y)

      if (length(x) < 2L || length(y) < 2L) {
        return(make_row(
          g1, g2, tp, length(x), length(y), s1, s2,
          s1[['mean']] - s2[['mean']],
          NA_real_, NA_real_, NA_real_, NA_real_,
          'insufficient_replicates'
        ))
      }

      warn <- NA_character_
      t_res <- tryCatch(
        withCallingHandlers(
          stats::t.test(x, y, var.equal = FALSE),
          warning = function(w) {
            warn <<- conditionMessage(w)
            invokeRestart('muffleWarning')
          }
        ),
        error = function(e) e
      )

      if (inherits(t_res, 'error')) {
        return(make_row(
          g1, g2, tp, length(x), length(y), s1, s2,
          s1[['mean']] - s2[['mean']],
          NA_real_, NA_real_, NA_real_, NA_real_,
          't_test_failed', conditionMessage(t_res)
        ))
      }

      make_row(
        g1, g2, tp, length(x), length(y), s1, s2,
        unname(t_res[['estimate']][[1]] - t_res[['estimate']][[2]]),
        unname(t_res[['conf.int']][1]),
        unname(t_res[['conf.int']][2]),
        unname(t_res[['statistic']]),
        unname(t_res[['p.value']]),
        'ok', warn
      )
    })
    result_list[[i]] <- do.call(rbind, tp_results)
  }
  out <- do.call(rbind, result_list)
  rownames(out) <- NULL
  out
}


# ================================================================
# 8. XLSX export
# ================================================================

# 8.1 Minimal plain-XLSX writer. Avoids optional workbook packages so
# project scripts can write simple publication tables from the renv baseline.
maturation_xlsx_col_ref <- function(col_idx) {
  letters <- character(0)
  while (col_idx > 0L) {
    rem <- (col_idx - 1L) %% 26L
    letters <- c(LETTERS[[rem + 1L]], letters)
    col_idx <- (col_idx - rem - 1L) %/% 26L
  }
  paste(letters, collapse = '')
}

maturation_xlsx_escape <- function(x) {
  x <- as.character(x)
  x <- gsub('&', '&amp;', x, fixed = TRUE)
  x <- gsub('<', '&lt;', x, fixed = TRUE)
  x <- gsub('>', '&gt;', x, fixed = TRUE)
  x <- gsub('"', '&quot;', x, fixed = TRUE)
  x
}

maturation_xlsx_safe_sheet_name <- function(x, existing = character()) {
  x <- gsub('[][\\\\/*?:]', '_', x)
  x <- substr(x, 1L, 31L)
  if (!nzchar(x)) {
    x <- 'Sheet'
  }
  candidate <- x
  suffix <- 1L
  while (candidate %in% existing) {
    suffix_text <- paste0('_', suffix)
    candidate <- paste0(substr(x, 1L, 31L - nchar(suffix_text)), suffix_text)
    suffix <- suffix + 1L
  }
  candidate
}

maturation_xlsx_cell_xml <- function(value, row_idx, col_idx, header = FALSE) {
  if (length(value) == 0L || is.na(value)) {
    return('')
  }
  ref <- paste0(maturation_xlsx_col_ref(col_idx), row_idx)
  style_id <- if (header) 1L else 0L
  if (inherits(value, 'Date')) {
    value <- as.character(value)
  }
  if (is.numeric(value) || is.integer(value)) {
    return(sprintf(
      '<c r="%s" s="%d"><v>%s</v></c>',
      ref,
      style_id,
      format(value, scientific = FALSE, trim = TRUE, digits = 15)
    ))
  }
  sprintf(
    '<c r="%s" s="%d" t="inlineStr"><is><t>%s</t></is></c>',
    ref,
    style_id,
    maturation_xlsx_escape(value)
  )
}

maturation_write_plain_xlsx <- function(sheets, output_path) {
  if (!is.list(sheets) || length(sheets) == 0L) {
    stop('sheets must be a non-empty named list.')
  }
  if (is.null(names(sheets)) || any(!nzchar(names(sheets)))) {
    stop('sheets must be named.')
  }

  output_path <- as.character(output_path)
  if (!grepl('^/', output_path)) {
    output_path <- file.path(getwd(), output_path)
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  output_path <- file.path(
    normalizePath(dirname(output_path), mustWork = TRUE),
    basename(output_path)
  )
  temp_root <- tempfile('maturation_xlsx_')
  dir.create(file.path(temp_root, '_rels'), recursive = TRUE)
  dir.create(file.path(temp_root, 'xl', '_rels'), recursive = TRUE)
  dir.create(file.path(temp_root, 'xl', 'worksheets'), recursive = TRUE)

  safe_names <- character(0)
  sheet_names <- vapply(names(sheets), function(sheet_name) {
    safe <- maturation_xlsx_safe_sheet_name(sheet_name, existing = safe_names)
    safe_names <<- c(safe_names, safe)
    safe
  }, character(1))

  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
      '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
      '<Default Extension="xml" ContentType="application/xml"/>',
      '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
      '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
      paste0(
        '<Override PartName="/xl/worksheets/sheet',
        seq_along(sheets),
        '.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
      ),
      '</Types>'
    ),
    file.path(temp_root, '[Content_Types].xml')
  )
  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
      '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>',
      '</Relationships>'
    ),
    file.path(temp_root, '_rels', '.rels')
  )
  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
      '<fonts count="2"><font><sz val="10"/><name val="Arial"/></font><font><b/><sz val="10"/><name val="Arial"/></font></fonts>',
      '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>',
      '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>',
      '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>',
      '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0"/></cellXfs>',
      '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>',
      '<dxfs count="0"/>',
      '<tableStyles count="0" defaultTableStyle="TableStyleMedium2" defaultPivotStyle="PivotStyleLight16"/>',
      '</styleSheet>'
    ),
    file.path(temp_root, 'xl', 'styles.xml')
  )

  workbook_sheets <- paste0(
    '<sheet name="',
    maturation_xlsx_escape(sheet_names),
    '" sheetId="',
    seq_along(sheet_names),
    '" r:id="rId',
    seq_along(sheet_names),
    '"/>'
  )
  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
      '<sheets>',
      workbook_sheets,
      '</sheets>',
      '</workbook>'
    ),
    file.path(temp_root, 'xl', 'workbook.xml')
  )
  workbook_rels <- paste0(
    '<Relationship Id="rId',
    seq_along(sheets),
    '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet',
    seq_along(sheets),
    '.xml"/>'
  )
  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
      workbook_rels,
      paste0(
        '<Relationship Id="rId',
        length(sheets) + 1L,
        '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
      ),
      '</Relationships>'
    ),
    file.path(temp_root, 'xl', '_rels', 'workbook.xml.rels')
  )

  for (sheet_idx in seq_along(sheets)) {
    sheet_df <- as.data.frame(sheets[[sheet_idx]], stringsAsFactors = FALSE)
    col_names <- colnames(sheet_df)
    header_cells <- vapply(seq_along(col_names), function(col_idx) {
      maturation_xlsx_cell_xml(col_names[[col_idx]], 1L, col_idx, header = TRUE)
    }, character(1))
    rows_xml <- c(
      sprintf('<row r="1">%s</row>', paste(header_cells, collapse = ''))
    )
    if (nrow(sheet_df) > 0L) {
      data_rows <- vapply(seq_len(nrow(sheet_df)), function(row_idx) {
        xlsx_row <- row_idx + 1L
        cells <- vapply(seq_along(col_names), function(col_idx) {
          maturation_xlsx_cell_xml(
            sheet_df[[col_idx]][[row_idx]],
            row_idx = xlsx_row,
            col_idx = col_idx
          )
        }, character(1))
        sprintf('<row r="%d">%s</row>', xlsx_row, paste(cells, collapse = ''))
      }, character(1))
      rows_xml <- c(rows_xml, data_rows)
    }
    dim_ref <- paste0(
      'A1:',
      maturation_xlsx_col_ref(max(1L, length(col_names))),
      max(1L, nrow(sheet_df) + 1L)
    )
    writeLines(
      c(
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
        paste0('<dimension ref="', dim_ref, '"/>'),
        '<sheetViews><sheetView workbookViewId="0" showGridLines="1"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>',
        '<sheetData>',
        rows_xml,
        '</sheetData>',
        '</worksheet>'
      ),
      file.path(temp_root, 'xl', 'worksheets', paste0('sheet', sheet_idx, '.xml'))
    )
  }

  if (file.exists(output_path)) {
    unlink(output_path)
  }
  zip_code <- paste0(
    'import os,sys,zipfile; ',
    'root=sys.argv[1]; out=sys.argv[2]; ',
    'zf=zipfile.ZipFile(out,"w",zipfile.ZIP_DEFLATED); ',
    '[zf.write(os.path.join(dp,f),os.path.relpath(os.path.join(dp,f),root)) ',
    'for dp,_,fs in os.walk(root) for f in fs]; ',
    'zf.close()'
  )
  zip_status <- system2(
    'python3',
    c(
      '-c',
      shQuote(zip_code),
      shQuote(temp_root),
      shQuote(output_path)
    )
  )
  unlink(temp_root, recursive = TRUE)
  if (!identical(zip_status, 0L)) {
    stop('failed to create XLSX file: ', output_path)
  }
  invisible(output_path)
}

# 8.2 Per-condition score summary.
summarize_maturation_scores <- function(score_df) {
  required_cols <- c(
    'condition',
    'genotype',
    'treatment',
    'timepoint',
    'timepoint_numeric',
    'group_value',
    'sample_set',
    'maturation_score'
  )
  check_required_columns(score_df, required_cols, 'score_df')
  split_key <- interaction(
    score_df[['condition']],
    score_df[['genotype']],
    score_df[['treatment']],
    score_df[['timepoint']],
    score_df[['timepoint_numeric']],
    score_df[['group_value']],
    score_df[['sample_set']],
    drop = TRUE
  )
  summary_df <- do.call(rbind, lapply(split(score_df, split_key), function(sub) {
    scores <- sub[['maturation_score']]
    data.frame(
      condition = as.character(sub[['condition']][[1]]),
      genotype = as.character(sub[['genotype']][[1]]),
      treatment = as.character(sub[['treatment']][[1]]),
      timepoint = as.character(sub[['timepoint']][[1]]),
      timepoint_numeric = sub[['timepoint_numeric']][[1]],
      group_value = as.character(sub[['group_value']][[1]]),
      sample_set = as.character(sub[['sample_set']][[1]]),
      n = sum(!is.na(scores)),
      mean_maturation_score = mean(scores, na.rm = TRUE),
      sd_maturation_score = stats::sd(scores, na.rm = TRUE),
      min_maturation_score = min(scores, na.rm = TRUE),
      max_maturation_score = max(scores, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_df) <- NULL
  summary_df[order(
    summary_df[['timepoint_numeric']],
    summary_df[['genotype']],
    summary_df[['treatment']],
    summary_df[['condition']]
  ), , drop = FALSE]
}

# 8.3 Save sample scores, condition means, and scoring metadata.
save_maturation_score_xlsx <- function(score_df,
                                       output_path,
                                       scoring_genes = NULL,
                                       lrt_results = NULL,
                                       pca_variance_explained = NULL,
                                       summary_table = NULL) {
  sample_cols <- intersect(
    c(
      'sampleid',
      'condition',
      'genotype',
      'treatment',
      'timepoint',
      'timepoint_numeric',
      'biological_replicate',
      'sequencing_batch',
      'group_value',
      'sample_set',
      'PC1',
      'PC2',
      'PC3',
      'maturation_score',
      'vector_PC1',
      'vector_PC2',
      'vector_PC3',
      'scoring_mode',
      'scoring_n_pcs',
      'n_lrt_genes'
    ),
    colnames(score_df)
  )
  sample_scores <- score_df[, sample_cols, drop = FALSE]
  sample_scores <- sample_scores[order(
    sample_scores[['timepoint_numeric']],
    sample_scores[['genotype']],
    sample_scores[['treatment']],
    sample_scores[['biological_replicate']],
    sample_scores[['sampleid']]
  ), , drop = FALSE]
  rownames(sample_scores) <- NULL

  group_means <- summarize_maturation_scores(score_df)

  sheets <- list(
    Sample_scores = sample_scores,
    Group_means = group_means
  )

  if (!is.null(pca_variance_explained)) {
    sheets[['PCA_variance']] <- data.frame(
      principal_component = names(pca_variance_explained),
      proportion_variance_explained = as.numeric(pca_variance_explained),
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(summary_table)) {
    sheets[['D9_summary']] <- as.data.frame(
      summary_table,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(scoring_genes)) {
    sheets[['Scoring_genes']] <- data.frame(
      gene_id = as.character(scoring_genes),
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(lrt_results)) {
    lrt_cols <- intersect(
      c('gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'),
      colnames(lrt_results)
    )
    sheets[['LRT_results']] <- lrt_results[, lrt_cols, drop = FALSE]
  }

  maturation_write_plain_xlsx(
    sheets = sheets,
    output_path = output_path
  )
}


# ================================================================
# 7. Boxplot
# ================================================================

# 7.0 format y-axis labels: whole numbers w/o decimals, fractional as-is.
format_axis_labels <- function(x, digits = NULL) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (!is.null(digits)) {
      label_digits <- if (abs((v * 2) - round(v * 2)) < 1e-9) {
        digits
      } else {
        max(digits, 2L)
      }
      return(formatC(v, format = 'f', digits = label_digits))
    }
    if (abs(v - round(v)) < 1e-9) {
      return(format(round(v), trim = TRUE))
    }
    format(v, trim = TRUE)
  }, character(1))
}

# 7.1 Min/max maturation score per group-tp combo.
# defines the vertical range bars on the boxplot.
compute_score_ranges <- function(score_df) {
  split_key <- interaction(
    score_df[['group_value']],
    score_df[['timepoint_numeric']],
    drop = TRUE
  )
  split_df <- split(score_df, split_key)
  rows <- lapply(split_df, function(sub) {
    scores <- sub[['maturation_score']]
    data.frame(
      group_value = as.character(sub[['group_value']][1]),
      timepoint_numeric = sub[['timepoint_numeric']][1],
      min_score = min(scores, na.rm = TRUE),
      max_score = max(scores, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# 7.2 format p-val as plotmath expression w/ italic p
format_pvalue_label <- function(pval) {
  if (is.na(pval)) return(NA_character_)
  if (pval < 0.001) {
    ex <- floor(log10(pval))
    cf <- signif(pval / 10^ex, 2)
    return(paste0(
      'italic(p)*"="*', cf,
      ' %*% 10^{', ex, '}'
    ))
  }
  p_str <- format(
    signif(pval, 2),
    scientific = FALSE, trim = TRUE
  )
  paste0('italic(p)*"="*', p_str)
}

# 7.3 Compute x/y positions for significance brackets.
# for each tp w/ a successful t-test, places a horiz bracket between
# the two group cols plus a p-val label above. stacks for multiple pairs.
compute_bracket_positions <- function(score_df, t_test_results, tp_breaks,
                                      y_limits = NULL) {
  empty <- data.frame(
    x_position = numeric(0), x_start = numeric(0), x_end = numeric(0),
    line_y = numeric(0), y_position = numeric(0), label = character(0),
    stringsAsFactors = FALSE
  )
  if (is.null(t_test_results) || nrow(t_test_results) == 0L) return(empty)

  if (is.null(y_limits)) {
    y_range <- range(score_df[['maturation_score']], na.rm = TRUE)
  } else {
    y_range <- range(y_limits)
  }
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span <= 0) y_span <- 0.1
  label_gap <- 0.012 * y_span

  rows <- lapply(seq_along(tp_breaks), function(i) {
    tp_val <- tp_breaks[[i]]
    sub_scores <- score_df[score_df[['timepoint_numeric']] == tp_val, , drop = FALSE]
    sub_tests <- t_test_results[
      t_test_results[['timepoint_numeric']] == tp_val & t_test_results[['status']] == 'ok', ,
      drop = FALSE
    ]
    if (nrow(sub_tests) == 0L) return(NULL)
    line_y_vals <- NULL
    if (is.null(y_limits)) {
      base_y <- max(sub_scores[['maturation_score']], na.rm = TRUE) + 0.04 * y_span
      line_y_vals <- base_y + (seq_len(nrow(sub_tests)) - 1L) * 0.10 * y_span
    } else {
      fixed_upper <- max(y_limits)
      available_top <- fixed_upper - 0.10 * y_span - label_gap
      lower_bound <- max(sub_scores[['maturation_score']], na.rm = TRUE) + 0.04 * y_span
      if (nrow(sub_tests) == 1L) {
        line_y_vals <- min(available_top, max(lower_bound, available_top))
      } else {
        step_target <- 0.10 * y_span
        start_y <- max(
          lower_bound,
          available_top - (nrow(sub_tests) - 1L) * step_target
        )
        if (start_y > available_top) {
          start_y <- available_top
        }
        step_eff <- min(
          step_target,
          (available_top - start_y) / (nrow(sub_tests) - 1L)
        )
        line_y_vals <- start_y + (seq_len(nrow(sub_tests)) - 1L) * step_eff
      }
    }
    do.call(rbind, lapply(seq_len(nrow(sub_tests)), function(j) {
      row <- sub_tests[j, , drop = FALSE]
      lx <- sub_scores[sub_scores[['group_value']] == row[['group_1']], 'x_position', drop = TRUE][1]
      rx <- sub_scores[sub_scores[['group_value']] == row[['group_2']], 'x_position', drop = TRUE][1]
      ly <- line_y_vals[[j]]
      y_pos <- ly + label_gap
      if (!is.null(y_limits)) {
        y_pos <- min(y_pos, max(y_limits) - 0.01 * y_span)
      }
      data.frame(
        x_position = mean(c(lx, rx)),
        x_start = min(lx, rx), x_end = max(lx, rx),
        line_y = ly, y_position = y_pos,
        label = format_pvalue_label(row[['p_value']]),
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(empty)
  rownames(out) <- NULL
  out
}

# 7.4 Assemble and save the maturation score dotplot as SVG.
# Colors auto-generated (colorblind-friendly), all points shape 21.
# layers: dashed lines at 0/1, min-max range bars, mean lines,
# sample points, significance brackets w/ p-val labels and hooks.
build_and_save_boxplot <- function(score_df, t_test_results,
                                   plot_title, output_path,
                                   group_colors = NULL,
                                   group_order = NULL,
                                   group_x_offsets = NULL,
                                   point_size = 1.1,
                                   point_stroke = 0.2,
                                   point_alpha = 1,
                                   y_limits = NULL,
                                   y_breaks = NULL,
                                   plot_style = 'dotplot',
                                   show_legend = TRUE,
                                   connect_group_medians = FALSE,
                                   show_dotplot_mean_bar = TRUE,
                                   point_jitter_width = NULL,
                                   point_jitter_height = 0.008,
                                   timepoint_gap = 0.85,
                                   box_width = 0.24,
                                   plot_width = 6.0,
                                   plot_height = 2.1,
                                   force_panel_height_in = NULL,
                                   axis_title_size = 9,
                                   axis_text_size = 7,
                                   reference_line_width = 0.3,
                                   axis_line_width = 0.25,
                                   range_line_width = 0.25,
                                   summary_line_width = 0.38,
                                   show_y_grid = TRUE,
                                   plot_margin = ggplot2::margin(4, 6, 4, 4),
                                   show_significance = TRUE,
                                   group_shapes = NULL,
                                   group_linetypes = NULL,
                                   proportional_timepoint_spacing = FALSE,
                                   y_axis_label = 'Maturation score',
                                   x_axis_tick_prefix = '',
                                   y_axis_label_digits = NULL,
                                   axis_tick_length = grid::unit(2, 'pt'),
                                   barplot_outline_colour = 'black',
                                   median_line_outline_width = 0,
                                   group_line_colors = NULL,
                                   reference_line_alpha = 1,
                                   x_expand_add = NULL) {
  plot_style <- match.arg(plot_style, c('dotplot', 'barplot', 'boxplot', 'lineplot'))
  group_levels <- unique(as.character(score_df[['group_value']]))
  if (!is.null(group_order)) {
    missing <- setdiff(group_levels, group_order)
    if (length(missing) > 0L) {
      stop(
        'boxplot_group_order is missing groups: ',
        paste(missing, collapse = ', ')
      )
    }
    group_levels <- group_order[group_order %in% group_levels]
  }
  group_colors <- assign_group_colors(
    group_levels, group_colors
  )
  legend_labels <- stats::setNames(
    gsub('_', ' ', group_levels), group_levels
  )
  if (!is.null(group_shapes)) {
    shape_vals <- group_shapes[group_levels]
  } else {
    shape_vals <- stats::setNames(
      rep(21L, length(group_levels)), group_levels
    )
  }
  if (!is.null(group_linetypes)) {
    linetype_vals <- group_linetypes[group_levels]
  } else {
    linetype_vals <- stats::setNames(
      rep('solid', length(group_levels)), group_levels
    )
  }

  bp_df <- score_df
  bp_df[['group_value']] <- factor(
    bp_df[['group_value']],
    levels = group_levels
  )
  tp_breaks <- sort(unique(bp_df[['timepoint_numeric']]))

  if (isTRUE(proportional_timepoint_spacing)) {
    tp_positions <- stats::setNames(
      as.numeric(tp_breaks), as.character(tp_breaks)
    )
  } else {
    tp_positions <- stats::setNames(
      1 + (seq_along(tp_breaks) - 1) * timepoint_gap,
      as.character(tp_breaks)
    )
  }
  offsets <- compute_group_x_offsets(
    group_levels,
    custom_offsets = group_x_offsets
  )
  if (identical(plot_style, 'lineplot')) {
    offsets[] <- 0
  }
  range_df <- compute_score_ranges(bp_df)

  bp_df[['x_position']] <- unname(
    tp_positions[as.character(bp_df[['timepoint_numeric']])]
  ) + unname(offsets[as.character(bp_df[['group_value']])])
  range_df[['x_position']] <- unname(
    tp_positions[as.character(range_df[['timepoint_numeric']])]
  ) + unname(offsets[as.character(range_df[['group_value']])])
  point_df <- bp_df[
    sample.int(nrow(bp_df)), ,
    drop = FALSE
  ]

  if (isTRUE(show_significance)) {
    ann_df <- compute_bracket_positions(
      bp_df, t_test_results, tp_breaks,
      y_limits = y_limits
    )
  } else {
    ann_df <- data.frame(
      x_position = numeric(0), x_start = numeric(0), x_end = numeric(0),
      line_y = numeric(0), y_position = numeric(0), label = character(0),
      stringsAsFactors = FALSE
    )
  }
  tp_labels <- paste0(x_axis_tick_prefix, tp_breaks)

  if (is.null(y_limits)) {
    score_range <- range(
      bp_df[['maturation_score']],
      na.rm = TRUE
    )
    y_span <- diff(score_range)
    if (y_span < 1e-6) y_span <- 0.1
    if (is.null(y_breaks)) {
      y_breaks <- pretty(score_range, n = 5)
    }
    y_lo <- min(y_breaks[1], score_range[1]) -
      y_span * 0.06
    y_hi <- max(
      y_breaks[length(y_breaks)],
      score_range[2]
    ) + y_span * 0.06
  } else {
    y_lo <- min(y_limits)
    y_hi <- max(y_limits)
    y_span <- y_hi - y_lo
    if (is.null(y_breaks)) {
      y_breaks <- pretty(c(y_lo, y_hi), n = 5)
    }
  }

  # hook segments (short downward caps at bracket ends)
  if (nrow(ann_df) > 0L) {
    hook_len <- y_span * 0.012
    hook_nudge <- y_span * 0.002
    hook_df <- data.frame(
      x = c(ann_df[['x_start']], ann_df[['x_end']]),
      y_top = rep(ann_df[['line_y']] + hook_nudge, 2),
      y_bot = rep(ann_df[['line_y']] - hook_len, 2),
      stringsAsFactors = FALSE
    )
  } else {
    hook_df <- data.frame(
      x = numeric(0), y_top = numeric(0),
      y_bot = numeric(0), stringsAsFactors = FALSE
    )
  }

  # per-group mean at each tp for mean bars
  split_key <- interaction(
    bp_df[['group_value']],
    bp_df[['timepoint_numeric']],
    drop = TRUE
  )
  mean_df <- do.call(rbind, lapply(
    split(bp_df, split_key),
    function(sub) {
      data.frame(
        group_value = as.character(sub[['group_value']][1]),
        timepoint_numeric = sub[['timepoint_numeric']][1],
        x_position = sub[['x_position']][1],
        mean_score = mean(
          sub[['maturation_score']],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(mean_df) <- NULL
  mean_hw <- 0.055

  if (identical(plot_style, 'lineplot')) {
    line_summary_df <- do.call(rbind, lapply(
      split(bp_df, split_key),
      function(sub) {
        data.frame(
          group_value = as.character(sub[['group_value']][1]),
          timepoint_numeric = sub[['timepoint_numeric']][1],
          x_position = sub[['x_position']][1],
          summary_score = mean(
            sub[['maturation_score']],
            na.rm = TRUE
          ),
          stringsAsFactors = FALSE
        )
      }
    ))
    rownames(line_summary_df) <- NULL
    d5_anchor <- line_summary_df[
      line_summary_df[['timepoint_numeric']] == min(tp_breaks),
      c('group_value', 'timepoint_numeric', 'x_position', 'summary_score'),
      drop = FALSE
    ]
    d5_anchor <- d5_anchor[
      d5_anchor[['group_value']] %in% c('Di21_ctrl', 'Tri21_ctrl'), ,
      drop = FALSE
    ]

    build_line_df <- function(line_group, anchor_group = NULL) {
      rows <- line_summary_df[
        line_summary_df[['group_value']] == line_group, ,
        drop = FALSE
      ]
      if (!is.null(anchor_group)) {
        anchor_row <- d5_anchor[
          d5_anchor[['group_value']] == anchor_group, ,
          drop = FALSE
        ]
        if (nrow(anchor_row) != 1L) {
          stop('lineplot D5 anchor could not be determined for ', line_group, '.')
        }
        anchor_row[['group_value']] <- line_group
        rows <- rbind(anchor_row, rows)
      }
      rows[['line_group']] <- line_group
      rows
    }

    # Direct groups: those with mean values at every shown timepoint.
    # Anchored groups (e.g. SAG) are drawn from the D5 ctrl anchor when
    # they start after the shared baseline. Skip any group not in the data.
    line_df_parts <- list()
    direct_groups <- intersect(
      c('Di21_ctrl', 'Tri21_ctrl'),
      group_levels
    )
    for (g in direct_groups) {
      line_df_parts[[g]] <- build_line_df(g)
    }
    anchored_group_map <- list(
      Di21_SAG = 'Di21_ctrl',
      Tri21_SAG = 'Tri21_ctrl'
    )
    for (g in names(anchored_group_map)) {
      if (!g %in% group_levels) next
      anchor_group <- anchored_group_map[[g]]
      if (!anchor_group %in% group_levels) next
      if (!anchor_group %in% d5_anchor[['group_value']]) next
      line_df_parts[[g]] <- build_line_df(
        g,
        anchor_group = anchor_group
      )
    }
    line_df <- do.call(rbind, line_df_parts)
    line_df[['group_value']] <- factor(line_df[['group_value']], levels = group_levels)
    line_df[['line_group']] <- factor(line_df[['line_group']], levels = group_levels)
    line_df <- line_df[order(
      line_df[['line_group']],
      line_df[['timepoint_numeric']]
    ), , drop = FALSE]
    line_type_values <- stats::setNames(
      ifelse(grepl('SAG$', group_levels), 'dashed', 'solid'),
      group_levels
    )
  } else {
    line_df <- NULL
    line_type_values <- NULL
  }

  median_df <- do.call(rbind, lapply(
    split(bp_df, split_key),
    function(sub) {
      data.frame(
        group_value = as.character(sub[['group_value']][1]),
        timepoint_numeric = sub[['timepoint_numeric']][1],
        x_position = sub[['x_position']][1],
        median_score = stats::median(
          sub[['maturation_score']],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(median_df) <- NULL
  median_df[['group_value']] <- factor(
    median_df[['group_value']],
    levels = group_levels
  )
  median_df <- median_df[order(
    median_df[['group_value']],
    median_df[['timepoint_numeric']]
  ), , drop = FALSE]

  if (is.null(point_jitter_width)) {
    point_jitter_width <- if (identical(plot_style, 'dotplot') || identical(plot_style, 'lineplot')) {
      0
    } else {
      0.018
    }
  }
  point_position <- if (point_jitter_width == 0 && point_jitter_height == 0) {
    ggplot2::position_identity()
  } else {
    ggplot2::position_jitter(
      width = point_jitter_width,
      height = point_jitter_height,
      seed = 42
    )
  }

  if (identical(plot_style, 'boxplot')) {
    boxplot_summary_df <- do.call(rbind, lapply(
      split(bp_df, split_key),
      function(sub) {
        stats <- grDevices::boxplot.stats(
          sub[['maturation_score']]
        )[['stats']]
        data.frame(
          group_value = as.character(sub[['group_value']][1]),
          timepoint_numeric = sub[['timepoint_numeric']][1],
          x_position = sub[['x_position']][1],
          whisker_low = stats[[1]],
          q1 = stats[[2]],
          median_score = stats[[3]],
          q3 = stats[[4]],
          whisker_high = stats[[5]],
          xmin = sub[['x_position']][1] - box_width / 2,
          xmax = sub[['x_position']][1] + box_width / 2,
          stringsAsFactors = FALSE
        )
      }
    ))
    rownames(boxplot_summary_df) <- NULL
    boxplot_summary_df[['group_value']] <- factor(
      boxplot_summary_df[['group_value']],
      levels = group_levels
    )
  } else {
    boxplot_summary_df <- NULL
  }

  p <- ggplot2::ggplot(bp_df, ggplot2::aes(
    x = .data[['x_position']],
    y = .data[['maturation_score']],
    fill = .data[['group_value']],
    shape = .data[['group_value']],
    group = interaction(
      .data[['x_position']],
      .data[['group_value']]
    )
  )) +
    ggplot2::geom_hline(
      yintercept = c(0, 1), linewidth = reference_line_width,
      colour = 'black', alpha = reference_line_alpha,
      linetype = 'dashed'
    )

  if (identical(plot_style, 'dotplot')) {
    p <- p +
      ggplot2::geom_segment(
        data = range_df,
        mapping = ggplot2::aes(
          x = .data[['x_position']],
          xend = .data[['x_position']],
          y = .data[['min_score']],
          yend = .data[['max_score']]
        ),
        inherit.aes = FALSE, linewidth = range_line_width,
        colour = 'black', lineend = 'round'
      ) +
      ggplot2::geom_point(
        data = point_df,
        size = point_size, stroke = point_stroke,
        colour = 'black', alpha = point_alpha,
        position = point_position
      )
    if (isTRUE(show_dotplot_mean_bar)) {
      p <- p +
        ggplot2::geom_segment(
          data = mean_df,
          mapping = ggplot2::aes(
            x = .data[['x_position']] - mean_hw,
            xend = .data[['x_position']] + mean_hw,
            y = .data[['mean_score']],
            yend = .data[['mean_score']]
          ),
          inherit.aes = FALSE, linewidth = range_line_width,
          colour = 'black'
        )
    }
  } else if (identical(plot_style, 'lineplot')) {
    p <- p +
      ggplot2::geom_line(
        data = line_df,
        mapping = ggplot2::aes(
          x = .data[['x_position']],
          y = .data[['summary_score']],
          colour = .data[['group_value']],
          linetype = .data[['line_group']],
          group = .data[['line_group']]
        ),
        inherit.aes = FALSE,
        linewidth = summary_line_width,
        lineend = 'round'
      ) +
      ggplot2::geom_point(
        data = point_df,
        size = point_size, stroke = point_stroke,
        colour = 'black', alpha = point_alpha,
        position = point_position
      )
  } else if (identical(plot_style, 'barplot')) {
    p <- p +
      ggplot2::geom_col(
        data = mean_df,
        mapping = ggplot2::aes(
          x = .data[['x_position']],
          y = .data[['mean_score']],
          fill = .data[['group_value']]
        ),
        inherit.aes = FALSE,
        position = ggplot2::position_identity(),
        width = box_width,
        alpha = 0.92,
        colour = barplot_outline_colour,
        linewidth = if (is.na(barplot_outline_colour)) 0 else axis_line_width
      ) +
      ggplot2::geom_point(
        data = point_df,
        size = point_size, stroke = point_stroke,
        colour = 'black', alpha = point_alpha,
        position = point_position
      )
  } else {
    p <- p +
      ggplot2::geom_rect(
        data = boxplot_summary_df,
        mapping = ggplot2::aes(
          xmin = .data[['xmin']],
          xmax = .data[['xmax']],
          ymin = .data[['q1']],
          ymax = .data[['q3']],
          fill = .data[['group_value']]
        ),
        inherit.aes = FALSE,
        colour = 'black',
        linewidth = range_line_width
      ) +
      ggplot2::geom_segment(
        data = boxplot_summary_df,
        mapping = ggplot2::aes(
          x = .data[['x_position']],
          xend = .data[['x_position']],
          y = .data[['whisker_low']],
          yend = .data[['q1']]
        ),
        inherit.aes = FALSE,
        colour = 'black',
        linewidth = range_line_width,
        lineend = 'butt'
      ) +
      ggplot2::geom_segment(
        data = boxplot_summary_df,
        mapping = ggplot2::aes(
          x = .data[['x_position']],
          xend = .data[['x_position']],
          y = .data[['q3']],
          yend = .data[['whisker_high']]
        ),
        inherit.aes = FALSE,
        colour = 'black',
        linewidth = range_line_width,
        lineend = 'butt'
      ) +
      ggplot2::geom_segment(
        data = boxplot_summary_df,
        mapping = ggplot2::aes(
          x = .data[['xmin']],
          xend = .data[['xmax']],
          y = .data[['median_score']],
          yend = .data[['median_score']]
        ),
        inherit.aes = FALSE,
        colour = 'black',
        linewidth = range_line_width,
        lineend = 'butt'
      ) +
      ggplot2::geom_point(
        data = point_df,
        size = point_size, stroke = point_stroke,
        colour = 'black', alpha = point_alpha,
        position = point_position
      )
  }

  if (isTRUE(connect_group_medians)) {
    dashed_groups <- names(linetype_vals)[linetype_vals != 'solid']
    solid_groups <- names(linetype_vals)[linetype_vals == 'solid']
    has_outlined_dashed <- median_line_outline_width > 0 &&
      length(dashed_groups) > 0L
    if (has_outlined_dashed) {
      # dashed groups: ggtrace outline (fill = line color,
      # color = outline color, stroke = outline width)
      outlined_df <- median_df[
        median_df[['group_value']] %in% dashed_groups, ,
        drop = FALSE
      ]
      if (nrow(outlined_df) > 0L) {
        p <- p +
          ggtrace::geom_line_trace(
            data = outlined_df,
            mapping = ggplot2::aes(
              x = .data[['x_position']],
              y = .data[['median_score']],
              fill = .data[['group_value']],
              linetype = .data[['group_value']],
              group = .data[['group_value']]
            ),
            inherit.aes = FALSE,
            color = 'black',
            stroke = median_line_outline_width,
            linewidth = summary_line_width,
            lineend = 'round'
          )
      }
      # solid groups: regular geom_line
      solid_df <- median_df[
        median_df[['group_value']] %in% solid_groups, ,
        drop = FALSE
      ]
      if (nrow(solid_df) > 0L) {
        p <- p +
          ggplot2::geom_line(
            data = solid_df,
            mapping = ggplot2::aes(
              x = .data[['x_position']],
              y = .data[['median_score']],
              colour = .data[['group_value']],
              linetype = .data[['group_value']],
              group = .data[['group_value']]
            ),
            inherit.aes = FALSE,
            linewidth = summary_line_width,
            lineend = 'round'
          )
      }
    } else {
      # no outline needed: all groups use regular geom_line
      p <- p +
        ggplot2::geom_line(
          data = median_df,
          mapping = ggplot2::aes(
            x = .data[['x_position']],
            y = .data[['median_score']],
            colour = .data[['group_value']],
            linetype = .data[['group_value']],
            group = .data[['group_value']]
          ),
          inherit.aes = FALSE,
          linewidth = summary_line_width,
          lineend = 'round'
        )
    }
  }

  if (isTRUE(show_significance)) {
    p <- p +
      ggplot2::geom_segment(
        data = ann_df,
        mapping = ggplot2::aes(
          x = .data[['x_start']],
          xend = .data[['x_end']],
          y = .data[['line_y']],
          yend = .data[['line_y']]
        ),
        inherit.aes = FALSE, linewidth = reference_line_width,
        colour = 'black', lineend = 'butt'
      ) +
      ggplot2::geom_segment(
        data = hook_df,
        mapping = ggplot2::aes(
          x = .data[['x']],
          xend = .data[['x']],
          y = .data[['y_top']],
          yend = .data[['y_bot']]
        ),
        inherit.aes = FALSE, linewidth = reference_line_width,
        colour = 'black', lineend = 'butt'
      ) +
      ggplot2::geom_text(
        data = ann_df,
        mapping = ggplot2::aes(
          x = .data[['x_position']],
          y = .data[['y_position']],
          label = .data[['label']]
        ),
        parse = TRUE,
        inherit.aes = FALSE, size = 2.2, vjust = 0,
        family = 'Arial', colour = 'black'
      )
  }

  if (is.null(x_expand_add)) {
    x_expand_add <- if (identical(plot_style, 'boxplot')) {
      c(0.14, 0.14)
    } else if (identical(plot_style, 'lineplot')) {
      c(0.30, 0.30)
    } else {
      c(0.20, 0.35)
    }
  }

  p <- p +
    ggplot2::scale_fill_manual(
      name = '', values = group_colors,
      labels = legend_labels, drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = '', values = shape_vals,
      labels = legend_labels, drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(tp_positions),
      labels = tp_labels, minor_breaks = NULL,
      expand = ggplot2::expansion(add = x_expand_add)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = function(x) {
        format_axis_labels(x, digits = y_axis_label_digits)
      }
    ) +
    ggplot2::coord_cartesian(
      ylim = c(y_lo, y_hi),
      clip = if (is.null(y_limits)) 'off' else 'on'
    ) +
    ggplot2::labs(
      x = 'Differentiation day',
      y = y_axis_label,
      title = plot_title
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(
          linetype = 0, stroke = 0.2,
          alpha = 1, size = 1.5
        ),
        keywidth = grid::unit(0.35, 'cm'),
        keyheight = grid::unit(0.30, 'cm')
      ),
      shape = ggplot2::guide_legend(
        keywidth = grid::unit(0.35, 'cm'),
        keyheight = grid::unit(0.30, 'cm')
      )
    ) +
    ggplot2::theme_classic(base_family = 'Arial') +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) 'top' else 'none',
      legend.direction = 'horizontal',
      legend.background = ggplot2::element_rect(
        fill = 'transparent', colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = 'transparent', colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = 'transparent', colour = NA
      ),
      legend.title = ggplot2::element_blank(),
      panel.grid.major.y = if (isTRUE(show_y_grid)) {
        ggplot2::element_line(
          color = 'grey95', linewidth = 0.18
        )
      } else {
        ggplot2::element_blank()
      },
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(
        colour = 'black', linewidth = axis_line_width
      ),
      axis.ticks = ggplot2::element_line(
        colour = 'black', linewidth = axis_line_width
      ),
      axis.ticks.length = axis_tick_length,
      axis.text.x = ggplot2::element_text(
        color = 'black', size = axis_text_size,
        lineheight = 0.9
      ),
      axis.text.y = ggplot2::element_text(
        color = 'black', size = axis_text_size
      ),
      axis.title = ggplot2::element_text(
        color = 'black', size = axis_title_size
      ),
      legend.text = ggplot2::element_text(
        color = 'black', size = 6.5
      ),
      legend.key.width = grid::unit(0.35, 'cm'),
      legend.key.height = grid::unit(0.30, 'cm'),
      legend.spacing.x = grid::unit(0.05, 'cm'),
      legend.spacing.y = grid::unit(0.01, 'cm'),
      legend.margin = ggplot2::margin(0, 0, 2, 0),
      plot.margin = plot_margin
    )

  line_colors <- if (!is.null(group_line_colors)) {
    group_line_colors
  } else {
    group_colors
  }
  if (isTRUE(connect_group_medians) || identical(plot_style, 'lineplot')) {
    p <- p +
      ggplot2::scale_color_manual(
        name = '', values = line_colors,
        labels = legend_labels, drop = FALSE
      )
    if (identical(plot_style, 'lineplot')) {
      p <- p +
        ggplot2::scale_linetype_manual(
          name = '', values = line_type_values,
          labels = legend_labels, drop = FALSE
        ) +
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            override.aes = list(
              shape = 21,
              linewidth = summary_line_width
            ),
            keywidth = grid::unit(0.35, 'cm'),
            keyheight = grid::unit(0.30, 'cm')
          ),
          linetype = ggplot2::guide_legend(
            override.aes = list(
              colour = group_colors[group_levels],
              linewidth = summary_line_width
            ),
            keywidth = grid::unit(0.35, 'cm'),
            keyheight = grid::unit(0.30, 'cm')
          )
        )
    } else {
      p <- p +
        ggplot2::scale_linetype_manual(
          name = '', values = linetype_vals,
          labels = legend_labels, drop = FALSE
        ) +
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            override.aes = list(
              shape = NA,
              linewidth = summary_line_width
            ),
            keywidth = grid::unit(0.35, 'cm'),
            keyheight = grid::unit(0.30, 'cm')
          ),
          linetype = 'none'
        )
    }
  }

  if (!is.null(force_panel_height_in)) {
    p <- p +
      ggh4x::force_panelsizes(
        rows = grid::unit(force_panel_height_in, 'in')
      )
  }

  ggplot2::ggsave(
    filename = output_path, plot = p,
    width = plot_width, height = plot_height, units = 'in',
    dpi = 600, device = svglite::svglite,
    bg = 'transparent',
    fix_text_size = FALSE
  )
  invisible(p)
}
