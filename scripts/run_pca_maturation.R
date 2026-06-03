# PCA Maturation Scoring
# Zoheb Khan, Moskowitz Lab
# v15 03/23/2026
#
# What this does:
#   1. computes VST from full dds (all samples, blind = FALSE)
#   2. identifies ref/ctrl samples via condition_col + reference_condition_names
#   3. runs DESeq2 LRT on ref/ctrl timepoints (or lrt_condition_names for cross-cohort)
#   4. trains PCA on ref samples using LRT-sig genes
#   5. identifies comparison samples via group_col + comparison_groups
#   6. Projects all samples (ref + comparison) into ref-defined PC space
#   7. scores each sample along a maturation axis (centroid_vector or best_fit_line)
#   8. optionally runs per-tp Welch t-tests (if welch_t_test_comparison provided)
#   9. saves 3 files: maturation scores CSV, LRT gene list CSV, boxplot SVG
#
# outputs use output_file_path for tables and optionally
# boxplot_output_path for the SVG:
#   {output_file_path}_maturation_scores.csv
#   {output_file_path}_lrt_genes.csv
#   {boxplot_output_path or output_file_path}_score_boxplot.svg
#
# Example:
#
#  # e.g. Di21_ctrl, Tri21_ctrl at days 0, 3, 7
#  # group_col_name = col w/ treatment-level groups like Di21_ctrl
#  # condition_col_name = col w/ per-tp condition labels like ctrl_Di21_D3
#
#   run_pca_maturation(
#     dds = my_dds,
#     metadata = my_metadata,
#     sample_id_col_name = 'sample_id',
#     timepoint_numeric_col_name = 'timepoint_numeric',
#     condition_col_name = 'condition',
#     group_col_name = 'treatment_group',
#     reference_condition_names = c('ctrl_Di21_D0', 'ctrl_Di21_D3', 'ctrl_Di21_D7'),
#     comparison_groups = c('Tri21_ctrl'),
#     vector_start_timepoint = 0,
#     vector_end_timepoint = 7,
#     output_file_path = 'results/maturation_score/ctrl_d0d7'
#   )

# source utilities from the same directory as this file when sourced,
# with a repo-root fallback for interactive sessions.
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (!is.na(this_file)) {
  source(file.path(dirname(this_file), 'run_pca_maturation_utils.R'))
} else {
  source('scripts/run_pca_maturation_utils.R')
}

# ================================================================
# Params
# ================================================================
#
# Required:
#
#' @param dds             [DESeqDataSet] Raw counts w/ all samples. VST computed from this once.
#' @param metadata        [data.frame]   One row per sample. Must contain the four col-name params below.
#' @param sample_id_col_name     [character(1)] Metadata col with sample IDs. Must match colnames(dds).
#' @param timepoint_numeric_col_name [character(1)] Metadata col with numeric tp values (e.g. 5, 7, 9).
#' @param condition_col_name     [character(1)] Metadata col with per-tp condition labels (e.g. 'ctrl_D7').
#' @param group_col_name         [character(1)] Metadata col for group-level identity (e.g. 'treatment_group').
#'   All ref samples must share one value in this col.
#' @param reference_condition_names [character vec] Values from condition_col defining the ref trajectory.
#'   These train the PCA and (by default) run the LRT. Must span >= 2 numeric tps.
#' @param comparison_groups      [character vec] Values from group_col for groups to project + score.
#'   All samples in these groups (all tps) included automatically. Can be character(0)
#'   when analysis_condition_names is used to select projection samples explicitly.
#' @param vector_start_timepoint [numeric(1)] Numeric tp for start of scoring vector (score = 0).
#' @param vector_end_timepoint   [numeric(1)] Numeric tp for end of scoring vector (score = 1).
#' @param output_file_path       [character(1)] Path prefix for all saved files. Dir created if needed.
#'   e.g. 'results/score/my_run' produces my_run_maturation_scores.csv and
#'   my_run_lrt_genes.csv.
#' @param boxplot_output_path    [character(1) or NULL] Optional path prefix for the SVG boxplot.
#'   Default NULL uses output_file_path, keeping all outputs in one directory.
#
# Optional:
#
#' @param lrt_condition_names    [character vec or NULL] Condition values for cross-cohort LRT. When set,
#'   LRT gene set comes from these samples instead of ref. Default NULL (LRT = ref).
#' @param lrt_dds                [DESeqDataSet or NULL] Optional external cohort for the LRT gene set.
#'   Default NULL uses `dds`.
#' @param lrt_metadata           [data.frame or NULL] Optional metadata matching `lrt_dds`.
#'   Default NULL uses `metadata`.
#' @param vst_batch_col_name     [character(1) or NULL] Optional metadata column to remove from
#'   VST space with `limma::removeBatchEffect()`, while protecting `condition_col_name`.
#' @param vst_batch2_col_name    [character(1) or NULL] Optional second nuisance column to remove
#'   from VST space alongside `vst_batch_col_name`.
#' @param lrt_block_col_names    [character vec or NULL] Optional nuisance columns to include in
#'   both the full and reduced LRT models. Use this for within-cohort blocking factors such as
#'   biological replicate.
#' @param lrt_tpm_matrix         [matrix/data.frame or NULL] Optional TPM matrix for post-LRT
#'   expression filtering on the LRT cohort. Row names must be gene IDs and column names
#'   must be sample IDs. Default NULL disables TPM filtering.
#' @param analysis_condition_names [character vec or NULL] Specific condition values to include in
#'   scoring/projection in addition to reference samples and any comparison_groups. Use this when
#'   projection samples should be restricted to a subset of timepoints instead of all samples from a group.
#' @param lrt_padj_cutoff        [numeric(1)] padj threshold for LRT gene selection. Default 0.05.
#' @param lrt_tpm_cutoff         [numeric(1) or NULL] Minimum TPM required for the optional LRT
#'   expression filter. Must be non-negative when set.
#' @param lrt_tpm_filter_mode    [character(1)] TPM filter summary mode. 'min_samples' requires
#'   `lrt_tpm_min_samples` samples with TPM >= cutoff. 'timepoint_mean_any' requires at least
#'   one LRT timepoint mean TPM >= cutoff. Default 'min_samples'.
#' @param lrt_tpm_min_samples    [integer(1) or NULL] Minimum number of LRT samples meeting
#'   `lrt_tpm_cutoff` when `lrt_tpm_filter_mode = 'min_samples'`. Must be >= 1 when used.
#' @param lrt_keep_direction     [character(1)] Which LRT genes to retain for scoring:
#'   'all' (default) keeps all LRT-significant genes; 'up' keeps only genes whose
#'   mean VST increases across the ordered LRT timepoints.
#' @param scoring_gene_ids       [character vec or NULL] Optional preselected gene IDs to use
#'   directly for PCA training/scoring. When set, the LRT is skipped and these genes define
#'   the PCA/scoring feature set after intersection with the VST matrix.
#' @param scoring_mode           [character(1)] Scoring axis method. 'centroid_vector' (default) draws
#'   axis between start/end tp centroids. 'best_fit_line' fits PC1 through ref samples in PCA space,
#'   normalized so start centroid = 0, end centroid = 1.
#' @param scoring_n_pcs          [integer(1)] Number of leading PCs to use for projection and scoring.
#'   Default 3 (PC1-PC3).
#' @param welch_t_test_comparison [list of length-2 char vecs or NULL] Group pairs for per-tp Welch
#'   t-tests on maturation scores. Each element c('group_A', 'group_B'). Default NULL.
#' @param group_colors           [named character vec or NULL] Optional colors keyed by group values.
#'   Must include every group shown in the boxplot. Default NULL (auto palette).
#' @param boxplot_group_order    [character vec or NULL] Optional explicit left-to-right group order
#'   within each timepoint. Must include every group shown in the boxplot.
#' @param boxplot_group_x_offsets [named numeric vec or NULL] Optional explicit x offsets keyed by
#'   group values for within-timepoint spacing. Must include every displayed group.
#' @param boxplot_point_size     [numeric(1)] Point size for the boxplot. Default 1.1.
#' @param boxplot_y_limits       [numeric(2) or NULL] Fixed y-axis limits for the boxplot.
#'   Default NULL uses the data range.
#' @param boxplot_y_breaks       [numeric vec or NULL] Explicit y-axis breaks for the boxplot.
#'   Default NULL uses pretty() breaks.
#' @param boxplot_title          [character or NULL] Title for boxplot. Default NULL (no title).
#' @param boxplot_plot_style     [character(1)] 'dotplot' (default), 'barplot',
#'   'boxplot', or 'lineplot'.
#' @param boxplot_show_legend    [logical(1)] Whether to show the group legend.
#'   Default TRUE.
#' @param boxplot_connect_group_medians [logical(1)] Whether to draw a same-color
#'   line through the per-group medians across displayed timepoints. Default FALSE.
#' @param boxplot_show_dotplot_mean_bar [logical(1)] Whether dotplots should draw
#'   the horizontal summary bar at the mean. Default TRUE.
#' @param boxplot_point_jitter_width [numeric(1) or NULL] Optional x jitter width
#'   for sample points. Default NULL uses the style-specific default.
#' @param boxplot_point_jitter_height [numeric(1)] Y jitter height for sample
#'   points. Default 0.008.
#' @param boxplot_timepoint_gap [numeric(1)] Spacing between adjacent displayed
#'   timepoints on the x-axis. Default 0.85.
#' @param boxplot_box_width     [numeric(1)] Width of each rendered box/bar.
#'   Default 0.24.
#' @param boxplot_width          [numeric(1)] Export width in inches. Default 6.0.
#' @param boxplot_height         [numeric(1)] Export height in inches. Default 2.1.
#' @param boxplot_axis_title_size [numeric(1)] Axis-title font size. Default 9.
#' @param boxplot_axis_text_size [numeric(1)] Axis tick-label font size. Default 7.
#' @param boxplot_x_axis_tick_prefix [character(1)] Prefix for x-axis tick labels.
#' @param boxplot_y_axis_label_digits [integer(1) or NULL] Decimal places for y-axis tick labels.
#' @param boxplot_axis_tick_length [grid unit] Length of axis tick marks.
#' @param boxplot_barplot_outline_colour [character(1)] Bar outline colour for barplot style.
#' @param show_significance      [logical(1)] Whether to draw significance brackets and p-value
#'   labels on the plot. Welch tests still run and remain available in outputs. Default TRUE.
#' @param display_timepoints     [numeric vec or NULL] Subset of tps shown in boxplot only.
#'   Scoring still uses all tps. Default NULL (show all).
#
# ================================================================

run_pca_maturation <- function(dds,
                               metadata,
                               sample_id_col_name,
                               timepoint_numeric_col_name,
                               condition_col_name,
                               group_col_name,
                               reference_condition_names,
                               comparison_groups,
                               vector_start_timepoint,
                               vector_end_timepoint,
                               output_file_path,
                               boxplot_output_path = NULL,
                               lrt_condition_names = NULL,
                               lrt_dds = NULL,
                               lrt_metadata = NULL,
                               vst_batch_col_name = NULL,
                               vst_batch2_col_name = NULL,
                               lrt_block_col_names = NULL,
                               lrt_tpm_matrix = NULL,
                               analysis_condition_names = NULL,
                               lrt_padj_cutoff = 0.05,
                               lrt_tpm_cutoff = NULL,
                               lrt_tpm_filter_mode = 'min_samples',
                               lrt_tpm_min_samples = NULL,
                               lrt_keep_direction = 'all',
                               scoring_gene_ids = NULL,
                               scoring_mode = 'centroid_vector',
                               scoring_n_pcs = 3L,
                               scoring_pca_scale = FALSE,
                               welch_t_test_comparison = NULL,
                               welch_t_test_bh_adjust = TRUE,
                               group_colors = NULL,
                               boxplot_group_order = NULL,
                               boxplot_group_x_offsets = NULL,
                               boxplot_point_size = 1.1,
                               boxplot_point_stroke = 0.2,
                               boxplot_point_alpha = 1,
                               boxplot_y_limits = NULL,
                               boxplot_y_breaks = NULL,
                               boxplot_title = NULL,
                               boxplot_plot_style = 'dotplot',
                               boxplot_show_legend = TRUE,
                               boxplot_connect_group_medians = FALSE,
                               boxplot_show_dotplot_mean_bar = TRUE,
                               boxplot_point_jitter_width = NULL,
                               boxplot_point_jitter_height = 0.008,
                               boxplot_timepoint_gap = 0.85,
                               boxplot_box_width = 0.24,
                               boxplot_width = 6.0,
                               boxplot_height = 2.1,
                               boxplot_force_panel_height_in = NULL,
                               boxplot_axis_title_size = 9,
                               boxplot_axis_text_size = 7,
                               boxplot_reference_line_width = 0.3,
                               boxplot_axis_line_width = 0.25,
                               boxplot_range_line_width = 0.25,
                               boxplot_summary_line_width = 0.38,
                               boxplot_show_y_grid = TRUE,
                               boxplot_plot_margin = ggplot2::margin(4, 6, 4, 4),
                               show_significance = TRUE,
                               display_timepoints = NULL,
                               boxplot_group_shapes = NULL,
                               boxplot_group_linetypes = NULL,
                               boxplot_proportional_timepoint_spacing = FALSE,
                               boxplot_y_axis_label = 'Maturation score',
                               boxplot_x_axis_tick_prefix = '',
                               boxplot_y_axis_label_digits = NULL,
                               boxplot_axis_tick_length = grid::unit(2, 'pt'),
                               boxplot_barplot_outline_colour = 'black',
                               boxplot_median_line_outline_width = 0,
                               boxplot_group_line_colors = NULL) {
  n_steps <- 5L

  # ================================================================
  # step 1: validate inputs
  # ================================================================
  log_step(1L, n_steps, 'validating inputs')
  check_required_packages()

  if (!inherits(dds, 'DESeqDataSet')) stop('dds must be a DESeqDataSet.')
  if (!is.data.frame(metadata)) stop('metadata must be a data frame.')
  if (is.null(lrt_dds)) {
    lrt_dds <- dds
  }
  if (is.null(lrt_metadata)) {
    lrt_metadata <- metadata
  }
  if (!inherits(lrt_dds, 'DESeqDataSet')) stop('lrt_dds must be a DESeqDataSet.')
  if (!is.data.frame(lrt_metadata)) stop('lrt_metadata must be a data frame.')
  if (!is.numeric(lrt_padj_cutoff) || length(lrt_padj_cutoff) != 1L || is.na(lrt_padj_cutoff)) stop('lrt_padj_cutoff must be a single numeric between 0 and 1.')
  if (lrt_padj_cutoff <= 0 || lrt_padj_cutoff >= 1) stop('lrt_padj_cutoff must be between 0 and 1.')
  lrt_tpm_filter_mode <- match.arg(
    lrt_tpm_filter_mode,
    c('min_samples', 'timepoint_mean_any')
  )
  use_lrt_tpm_filter <- !is.null(lrt_tpm_cutoff)
  if (use_lrt_tpm_filter) {
    if (is.null(lrt_tpm_matrix)) {
      stop('lrt_tpm_matrix must be provided when lrt_tpm_cutoff is set.')
    }
    if (!is.matrix(lrt_tpm_matrix) && !is.data.frame(lrt_tpm_matrix)) {
      stop('lrt_tpm_matrix must be a matrix or data frame.')
    }
    lrt_tpm_matrix <- as.matrix(lrt_tpm_matrix)
    if (!is.numeric(lrt_tpm_matrix)) stop('lrt_tpm_matrix must be numeric.')
    if (is.null(rownames(lrt_tpm_matrix)) || is.null(colnames(lrt_tpm_matrix))) {
      stop('lrt_tpm_matrix must have gene IDs as row names and sample IDs as column names.')
    }
    if (!is.numeric(lrt_tpm_cutoff) || length(lrt_tpm_cutoff) != 1L || is.na(lrt_tpm_cutoff) || lrt_tpm_cutoff < 0) {
      stop('lrt_tpm_cutoff must be a single non-negative numeric value.')
    }
    if (identical(lrt_tpm_filter_mode, 'min_samples')) {
      if (!is.numeric(lrt_tpm_min_samples) || length(lrt_tpm_min_samples) != 1L || is.na(lrt_tpm_min_samples)) {
        stop('lrt_tpm_min_samples must be a single positive integer value when lrt_tpm_filter_mode = min_samples.')
      }
      lrt_tpm_min_samples <- as.integer(lrt_tpm_min_samples)
      if (lrt_tpm_min_samples < 1L) stop('lrt_tpm_min_samples must be >= 1.')
    } else if (!is.null(lrt_tpm_min_samples)) {
      stop('lrt_tpm_min_samples must be NULL when lrt_tpm_filter_mode = timepoint_mean_any.')
    }
  } else if (!is.null(lrt_tpm_matrix)) {
    lrt_tpm_matrix <- as.matrix(lrt_tpm_matrix)
  } else if (!is.null(lrt_tpm_min_samples)) {
    stop('lrt_tpm_min_samples requires lrt_tpm_cutoff.')
  }
  if (!is.null(vst_batch_col_name)) vst_batch_col_name <- as.character(vst_batch_col_name[[1]])
  if (!is.null(vst_batch2_col_name)) vst_batch2_col_name <- as.character(vst_batch2_col_name[[1]])
  if (!is.numeric(scoring_n_pcs) || length(scoring_n_pcs) != 1L || is.na(scoring_n_pcs)) {
    stop('scoring_n_pcs must be a single positive integer.')
  }
  scoring_n_pcs <- as.integer(scoring_n_pcs[[1]])
  if (scoring_n_pcs < 1L) {
    stop('scoring_n_pcs must be >= 1.')
  }
  if (is.null(lrt_block_col_names)) {
    lrt_block_col_names <- character(0)
  } else {
    lrt_block_col_names <- as.character(lrt_block_col_names)
  }
  lrt_keep_direction <- match.arg(lrt_keep_direction, c('all', 'up'))
  if (is.null(scoring_gene_ids)) {
    scoring_gene_ids <- character(0)
  } else {
    scoring_gene_ids <- unique(as.character(scoring_gene_ids))
    scoring_gene_ids <- scoring_gene_ids[!is.na(scoring_gene_ids) & nzchar(scoring_gene_ids)]
  }
  scoring_mode <- match.arg(scoring_mode, c('centroid_vector', 'best_fit_line'))
  boxplot_plot_style <- match.arg(boxplot_plot_style, c('dotplot', 'barplot', 'boxplot', 'lineplot'))
  if (!is.logical(boxplot_show_legend) || length(boxplot_show_legend) != 1L || is.na(boxplot_show_legend)) {
    stop('boxplot_show_legend must be a single TRUE/FALSE value.')
  }
  if (!is.logical(boxplot_connect_group_medians) || length(boxplot_connect_group_medians) != 1L || is.na(boxplot_connect_group_medians)) {
    stop('boxplot_connect_group_medians must be a single TRUE/FALSE value.')
  }
  if (!is.logical(boxplot_show_dotplot_mean_bar) || length(boxplot_show_dotplot_mean_bar) != 1L || is.na(boxplot_show_dotplot_mean_bar)) {
    stop('boxplot_show_dotplot_mean_bar must be a single TRUE/FALSE value.')
  }
  if (!is.null(boxplot_point_jitter_width)) {
    if (!is.numeric(boxplot_point_jitter_width) ||
      length(boxplot_point_jitter_width) != 1L ||
      is.na(boxplot_point_jitter_width) ||
      boxplot_point_jitter_width < 0) {
      stop('boxplot_point_jitter_width must be NULL or a single non-negative numeric value.')
    }
  }
  if (!is.numeric(boxplot_point_jitter_height) ||
    length(boxplot_point_jitter_height) != 1L ||
    is.na(boxplot_point_jitter_height) ||
    boxplot_point_jitter_height < 0) {
    stop('boxplot_point_jitter_height must be a single non-negative numeric value.')
  }
  if (!is.numeric(boxplot_timepoint_gap) ||
    length(boxplot_timepoint_gap) != 1L ||
    is.na(boxplot_timepoint_gap) ||
    boxplot_timepoint_gap <= 0) {
    stop('boxplot_timepoint_gap must be a single positive numeric value.')
  }
  if (!is.numeric(boxplot_box_width) ||
    length(boxplot_box_width) != 1L ||
    is.na(boxplot_box_width) ||
    boxplot_box_width <= 0) {
    stop('boxplot_box_width must be a single positive numeric value.')
  }
  if (!is.numeric(boxplot_width) ||
    length(boxplot_width) != 1L ||
    is.na(boxplot_width) ||
    boxplot_width <= 0) {
    stop('boxplot_width must be a single positive numeric value.')
  }
  if (!is.numeric(boxplot_height) ||
    length(boxplot_height) != 1L ||
    is.na(boxplot_height) ||
    boxplot_height <= 0) {
    stop('boxplot_height must be a single positive numeric value.')
  }
  if (!is.numeric(boxplot_point_stroke) ||
    length(boxplot_point_stroke) != 1L ||
    is.na(boxplot_point_stroke) ||
    boxplot_point_stroke < 0) {
    stop('boxplot_point_stroke must be a single non-negative numeric value.')
  }
  if (!is.numeric(boxplot_axis_title_size) ||
    length(boxplot_axis_title_size) != 1L ||
    is.na(boxplot_axis_title_size) ||
    boxplot_axis_title_size <= 0) {
    stop('boxplot_axis_title_size must be a single positive numeric value.')
  }
  if (!is.numeric(boxplot_axis_text_size) ||
    length(boxplot_axis_text_size) != 1L ||
    is.na(boxplot_axis_text_size) ||
    boxplot_axis_text_size <= 0) {
    stop('boxplot_axis_text_size must be a single positive numeric value.')
  }
  if (!is.numeric(boxplot_reference_line_width) ||
    length(boxplot_reference_line_width) != 1L ||
    is.na(boxplot_reference_line_width) ||
    boxplot_reference_line_width < 0) {
    stop('boxplot_reference_line_width must be a single non-negative numeric value.')
  }
  if (!is.numeric(boxplot_axis_line_width) ||
    length(boxplot_axis_line_width) != 1L ||
    is.na(boxplot_axis_line_width) ||
    boxplot_axis_line_width < 0) {
    stop('boxplot_axis_line_width must be a single non-negative numeric value.')
  }
  if (!is.numeric(boxplot_range_line_width) ||
    length(boxplot_range_line_width) != 1L ||
    is.na(boxplot_range_line_width) ||
    boxplot_range_line_width < 0) {
    stop('boxplot_range_line_width must be a single non-negative numeric value.')
  }
  if (!is.numeric(boxplot_summary_line_width) ||
    length(boxplot_summary_line_width) != 1L ||
    is.na(boxplot_summary_line_width) ||
    boxplot_summary_line_width < 0) {
    stop('boxplot_summary_line_width must be a single non-negative numeric value.')
  }
  if (!is.logical(boxplot_show_y_grid) ||
    length(boxplot_show_y_grid) != 1L ||
    is.na(boxplot_show_y_grid)) {
    stop('boxplot_show_y_grid must be a single TRUE/FALSE value.')
  }
  if (!is.logical(show_significance) || length(show_significance) != 1L || is.na(show_significance)) {
    stop('show_significance must be a single TRUE/FALSE value.')
  }

  check_required_columns(metadata, c(sample_id_col_name, timepoint_numeric_col_name, condition_col_name, group_col_name), 'metadata')
  check_required_columns(lrt_metadata, c(sample_id_col_name, timepoint_numeric_col_name, condition_col_name, group_col_name), 'lrt_metadata')
  extra_metadata_cols <- unique(stats::na.omit(c(
    vst_batch_col_name, vst_batch2_col_name, lrt_block_col_names
  )))
  if (length(extra_metadata_cols) > 0L) {
    check_required_columns(metadata, extra_metadata_cols, 'metadata')
    check_required_columns(lrt_metadata, extra_metadata_cols, 'lrt_metadata')
  }
  if (!is.numeric(metadata[[timepoint_numeric_col_name]])) stop(timepoint_numeric_col_name, ' must be a numeric column.')
  if (!is.numeric(lrt_metadata[[timepoint_numeric_col_name]])) stop('lrt_metadata ', timepoint_numeric_col_name, ' must be a numeric column.')
  if (anyDuplicated(metadata[[sample_id_col_name]]) > 0L) stop('sample IDs in metadata must be unique.')
  if (anyDuplicated(lrt_metadata[[sample_id_col_name]]) > 0L) stop('sample IDs in lrt_metadata must be unique.')

  reference_condition_names <- as.character(reference_condition_names)
  comparison_groups <- as.character(comparison_groups)
  if (is.null(analysis_condition_names)) {
    analysis_condition_names <- character(0)
  } else {
    analysis_condition_names <- as.character(analysis_condition_names)
  }
  vector_start_timepoint <- as.numeric(vector_start_timepoint)
  vector_end_timepoint <- as.numeric(vector_end_timepoint)

  missing_ref <- setdiff(reference_condition_names, as.character(metadata[[condition_col_name]]))
  if (length(missing_ref) > 0L) stop('reference_condition_names not found in ', condition_col_name, ': ', paste(head(missing_ref, 5), collapse = ', '))

  if (length(comparison_groups) > 0L) {
    missing_comp <- setdiff(comparison_groups, as.character(metadata[[group_col_name]]))
    if (length(missing_comp) > 0L) stop('comparison_groups not found in ', group_col_name, ': ', paste(head(missing_comp, 5), collapse = ', '))
  }

  if (length(analysis_condition_names) > 0L) {
    missing_analysis <- setdiff(analysis_condition_names, as.character(metadata[[condition_col_name]]))
    if (length(missing_analysis) > 0L) stop('analysis_condition_names not found in ', condition_col_name, ': ', paste(head(missing_analysis, 5), collapse = ', '))
  }

  if (length(vector_start_timepoint) != 1L || is.na(vector_start_timepoint)) stop('vector_start_timepoint must be a single numeric value.')
  if (length(vector_end_timepoint) != 1L || is.na(vector_end_timepoint)) stop('vector_end_timepoint must be a single numeric value.')
  if (identical(vector_start_timepoint, vector_end_timepoint)) stop('vector_start_timepoint and vector_end_timepoint must differ.')

  # output dirs from path prefixes
  output_dir <- dirname(output_file_path)
  if (!nzchar(basename(output_file_path))) stop('output_file_path must end with a filename prefix.')
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(boxplot_output_path)) {
    boxplot_output_path <- output_file_path
  } else {
    if (!nzchar(basename(boxplot_output_path))) stop('boxplot_output_path must end with a filename prefix.')
    dir.create(dirname(boxplot_output_path), recursive = TRUE, showWarnings = FALSE)
  }

  if (!is.null(display_timepoints)) {
    display_timepoints <- as.numeric(display_timepoints)
    display_timepoints <- display_timepoints[!is.na(display_timepoints)]
  }
  if (!is.null(boxplot_group_order)) {
    boxplot_group_order <- as.character(boxplot_group_order)
  }
  if (!is.null(boxplot_group_x_offsets)) {
    offset_names <- names(boxplot_group_x_offsets)
    boxplot_group_x_offsets <- as.numeric(boxplot_group_x_offsets)
    names(boxplot_group_x_offsets) <- offset_names
    if (is.null(offset_names)) {
      stop('boxplot_group_x_offsets must be a named numeric vector.')
    }
  }
  if (!is.numeric(boxplot_point_size) ||
    length(boxplot_point_size) != 1L ||
    is.na(boxplot_point_size) ||
    boxplot_point_size <= 0) {
    stop('boxplot_point_size must be a single positive numeric value.')
  }
  if (!is.null(boxplot_y_limits)) {
    boxplot_y_limits <- as.numeric(boxplot_y_limits)
    if (length(boxplot_y_limits) != 2L || anyNA(boxplot_y_limits)) {
      stop('boxplot_y_limits must contain exactly 2 numeric values.')
    }
    if (boxplot_y_limits[[1]] >= boxplot_y_limits[[2]]) {
      stop('boxplot_y_limits must be increasing.')
    }
  }
  if (!is.null(boxplot_y_breaks)) {
    boxplot_y_breaks <- as.numeric(boxplot_y_breaks)
    boxplot_y_breaks <- boxplot_y_breaks[!is.na(boxplot_y_breaks)]
  }

  # ================================================================
  # Step 2: compute VST on all samples (once)
  # ================================================================
  log_step(2L, n_steps, 'computing VST on all samples (blind = FALSE)')
  vst_info <- compute_vst_matrix(
    dds = dds,
    metadata = metadata,
    sample_id_col_name = sample_id_col_name,
    condition_col_name = condition_col_name,
    batch_col_name = vst_batch_col_name,
    batch2_col_name = vst_batch2_col_name
  )
  vst_mat <- vst_info[['matrix']]
  if (length(vst_info[['applied_batch_cols']]) > 0L) {
    message(
      '  VST batch correction: removed ',
      paste(vst_info[['applied_batch_cols']], collapse = ', '),
      ' while protecting ', condition_col_name
    )
  }
  if (length(vst_info[['dropped_batch_cols']]) > 0L) {
    message(
      '  VST batch correction: skipped ',
      paste(vst_info[['dropped_batch_cols']], collapse = ', '),
      ' because they had < 2 levels.'
    )
  }

  # ================================================================
  # step 3: select sample subsets and run LRT
  # ================================================================
  log_step(3L, n_steps, 'selecting sample subsets and running LRT')

  # ref samples: rows where condition matches reference_condition_names
  ref_meta <- metadata[as.character(metadata[[condition_col_name]]) %in% reference_condition_names, , drop = FALSE]
  if (nrow(ref_meta) == 0L) stop('reference_condition_names matched no samples.')

  ref_timepoints <- sort(unique(ref_meta[[timepoint_numeric_col_name]]))
  if (length(ref_timepoints) < 2L) stop('reference must span at least 2 timepoints.')
  if (!vector_start_timepoint %in% ref_timepoints) stop('vector_start_timepoint (', vector_start_timepoint, ') not in reference timepoints.')
  if (!vector_end_timepoint %in% ref_timepoints) stop('vector_end_timepoint (', vector_end_timepoint, ') not in reference timepoints.')

  # LRT samples: default to ref, or separate cohort if lrt_condition_names specified
  if (is.null(lrt_condition_names)) {
    lrt_meta <- ref_meta
  } else {
    lrt_condition_names <- as.character(lrt_condition_names)
    missing_lrt <- setdiff(lrt_condition_names, as.character(lrt_metadata[[condition_col_name]]))
    if (length(missing_lrt) > 0L) stop('lrt_condition_names not found in lrt ', condition_col_name, ': ', paste(head(missing_lrt, 5), collapse = ', '))
    lrt_meta <- lrt_metadata[as.character(lrt_metadata[[condition_col_name]]) %in% lrt_condition_names, , drop = FALSE]
  }
  if (nrow(lrt_meta) == 0L) stop('LRT conditions matched no samples.')
  lrt_timepoints <- sort(unique(lrt_meta[[timepoint_numeric_col_name]]))
  if (length(lrt_timepoints) < 2L) stop('LRT needs at least 2 timepoints.')

  # ref must belong to a single group
  n_ref_groups <- length(unique(as.character(ref_meta[[group_col_name]])))
  if (n_ref_groups != 1L) stop('reference samples must share one value in ', group_col_name, '.')

  # comparison samples: all rows matching comparison_groups
  if (length(comparison_groups) > 0L) {
    comp_meta <- metadata[as.character(metadata[[group_col_name]]) %in% comparison_groups, , drop = FALSE]
  } else {
    comp_meta <- metadata[0, , drop = FALSE]
  }

  # optional extra analysis samples: specific conditions to score/project
  if (length(analysis_condition_names) > 0L) {
    extra_meta <- metadata[as.character(metadata[[condition_col_name]]) %in% analysis_condition_names, , drop = FALSE]
  } else {
    extra_meta <- metadata[0, , drop = FALSE]
  }

  if (nrow(comp_meta) == 0L && nrow(extra_meta) == 0L) {
    stop('comparison_groups and analysis_condition_names matched no samples.')
  }

  ref_ids <- as.character(ref_meta[[sample_id_col_name]])
  comp_ids <- as.character(comp_meta[[sample_id_col_name]])
  extra_ids <- as.character(extra_meta[[sample_id_col_name]])
  analysis_ids <- unique(c(ref_ids, comp_ids, extra_ids))
  lrt_ids <- as.character(lrt_meta[[sample_id_col_name]])

  missing_analysis_in_dds <- setdiff(analysis_ids, colnames(dds))
  if (length(missing_analysis_in_dds) > 0L) stop('analysis sample IDs not in dds: ', paste(head(missing_analysis_in_dds, 5), collapse = ', '))
  missing_lrt_in_dds <- setdiff(lrt_ids, colnames(lrt_dds))
  if (length(missing_lrt_in_dds) > 0L) stop('LRT sample IDs not in lrt_dds: ', paste(head(missing_lrt_in_dds, 5), collapse = ', '))
  if (use_lrt_tpm_filter) {
    missing_lrt_in_tpm <- setdiff(lrt_ids, colnames(lrt_tpm_matrix))
    if (length(missing_lrt_in_tpm) > 0L) {
      stop('LRT sample IDs not in lrt_tpm_matrix: ', paste(head(missing_lrt_in_tpm, 5), collapse = ', '))
    }
  }

  lrt_meta <- reorder_by_sample_id(lrt_meta, sample_id_col_name, lrt_ids)
  analysis_meta <- reorder_by_sample_id(
    metadata[as.character(metadata[[sample_id_col_name]]) %in% analysis_ids, , drop = FALSE],
    sample_id_col_name, analysis_ids
  )

  lrt_tpm_subset <- NULL
  lrt_expression_keep <- NULL
  if (use_lrt_tpm_filter) {
    lrt_tpm_subset <- lrt_tpm_matrix[, lrt_ids, drop = FALSE]
    if (identical(lrt_tpm_filter_mode, 'min_samples')) {
      lrt_expression_keep <- rowSums(
        lrt_tpm_subset >= lrt_tpm_cutoff,
        na.rm = TRUE
      ) >= lrt_tpm_min_samples
    } else {
      lrt_tpm_tp_means <- vapply(
        lrt_timepoints,
        function(tp) {
          cols <- lrt_ids[lrt_meta[[timepoint_numeric_col_name]] == tp]
          rowMeans(lrt_tpm_subset[, cols, drop = FALSE], na.rm = TRUE)
        },
        numeric(nrow(lrt_tpm_subset))
      )
      rownames(lrt_tpm_tp_means) <- rownames(lrt_tpm_subset)
      colnames(lrt_tpm_tp_means) <- as.character(lrt_timepoints)
      lrt_expression_keep <- apply(lrt_tpm_tp_means, 1, function(x) {
        any(x >= lrt_tpm_cutoff, na.rm = TRUE)
      })
    }
    if (!any(lrt_expression_keep)) {
      stop('no genes passed the pre-LRT TPM expression filter.')
    }
  }

  tp_col <- timepoint_numeric_col_name
  lrt_results <- NULL
  if (length(scoring_gene_ids) == 0L) {
    # run LRT: full ~ blocking factors + timepoint, reduced ~ blocking factors
    lrt_model_meta <- lrt_meta
    lrt_model_meta[[tp_col]] <- factor(
      lrt_model_meta[[tp_col]],
      levels = lrt_timepoints
    )
    kept_lrt_block_cols <- character(0)
    dropped_lrt_block_cols <- character(0)
    if (length(lrt_block_col_names) > 0L) {
      for (block_col in lrt_block_col_names) {
        lrt_model_meta[[block_col]] <- droplevels(
          factor(as.character(lrt_model_meta[[block_col]]))
        )
        if (length(unique(as.character(lrt_model_meta[[block_col]]))) < 2L) {
          dropped_lrt_block_cols <- c(dropped_lrt_block_cols, block_col)
        } else {
          kept_lrt_block_cols <- c(kept_lrt_block_cols, block_col)
        }
      }
    }
    full_terms <- c(kept_lrt_block_cols, tp_col)
    full_fml <- stats::as.formula(
      paste('~', paste(full_terms, collapse = ' + '))
    )
    if (length(kept_lrt_block_cols) > 0L) {
      reduced_fml <- stats::as.formula(
        paste('~', paste(kept_lrt_block_cols, collapse = ' + '))
      )
    } else {
      reduced_fml <- ~1
    }
    if (length(kept_lrt_block_cols) > 0L) {
      message(
        '  LRT blocking factors: ',
        paste(kept_lrt_block_cols, collapse = ', ')
      )
    }
    if (length(dropped_lrt_block_cols) > 0L) {
      message(
        '  LRT blocking factors skipped (< 2 levels): ',
        paste(dropped_lrt_block_cols, collapse = ', ')
      )
    }

    dds_lrt <- DESeq2::DESeqDataSetFromMatrix(
      countData = DESeq2::counts(
        lrt_dds,
        normalized = FALSE
      )[if (is.null(lrt_expression_keep)) {
        rownames(DESeq2::counts(lrt_dds, normalized = FALSE))
      } else {
        rownames(lrt_tpm_subset)[lrt_expression_keep]
      }, lrt_ids, drop = FALSE],
      colData = lrt_model_meta,
      design = full_fml
    )
    dds_lrt <- DESeq2::DESeq(
      dds_lrt,
      test = 'LRT',
      reduced = reduced_fml, quiet = TRUE
    )

    lrt_results <- as.data.frame(DESeq2::results(dds_lrt, alpha = lrt_padj_cutoff))
    lrt_results[['gene_id']] <- rownames(lrt_results)
    lrt_results[['pass_lrt_padj']] <- !is.na(lrt_results[['padj']]) & lrt_results[['padj']] < lrt_padj_cutoff

    if (identical(lrt_keep_direction, 'up')) {
      if (identical(lrt_dds, dds)) {
        lrt_vst <- vst_mat[, lrt_ids, drop = FALSE]
      } else {
        lrt_vst <- compute_vst_matrix(
          dds = lrt_dds,
          metadata = lrt_metadata,
          sample_id_col_name = sample_id_col_name,
          condition_col_name = condition_col_name,
          batch_col_name = vst_batch_col_name,
          batch2_col_name = vst_batch2_col_name
        )[['matrix']][, lrt_ids, drop = FALSE]
      }
      if (!is.null(lrt_expression_keep)) {
        lrt_vst <- lrt_vst[rownames(lrt_tpm_subset)[lrt_expression_keep], , drop = FALSE]
      }
      lrt_tp_means <- vapply(
        lrt_timepoints,
        function(tp) {
          cols <- lrt_ids[lrt_meta[[tp_col]] == tp]
          rowMeans(lrt_vst[, cols, drop = FALSE], na.rm = TRUE)
        },
        numeric(nrow(lrt_vst))
      )
      rownames(lrt_tp_means) <- rownames(lrt_vst)
      colnames(lrt_tp_means) <- as.character(lrt_timepoints)

      pass_up_direction <- apply(lrt_tp_means, 1, function(x) {
        if (length(x) < 2L) {
          return(FALSE)
        }
        all(diff(x) > 0)
      })

      lrt_results[['pass_lrt_direction']] <- pass_up_direction[lrt_results[['gene_id']]]
    } else {
      lrt_results[['pass_lrt_direction']] <- TRUE
    }

    if (use_lrt_tpm_filter) {
      lrt_results[['pass_lrt_expression']] <- TRUE
    } else {
      lrt_results[['pass_lrt_expression']] <- TRUE
    }

    lrt_results[['pass_lrt_gene_filter']] <- lrt_results[['pass_lrt_padj']] &
      lrt_results[['pass_lrt_direction']] &
      lrt_results[['pass_lrt_expression']]
    lrt_results <- lrt_results[
      order(lrt_results[['padj']], lrt_results[['pvalue']], na.last = TRUE),
      c(
        'gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat',
        'pvalue', 'padj', 'pass_lrt_padj', 'pass_lrt_direction',
        'pass_lrt_expression',
        'pass_lrt_gene_filter'
      ),
      drop = FALSE
    ]
    rownames(lrt_results) <- NULL
    scoring_genes <- lrt_results[['gene_id']][lrt_results[['pass_lrt_gene_filter']]]
    if (length(scoring_genes) == 0L) stop(
      'no genes passed the LRT gene filter (padj < ', lrt_padj_cutoff,
      ', direction = ', lrt_keep_direction, ').'
    )
    if (!all(scoring_genes %in% rownames(vst_mat))) stop('VST matrix is missing genes selected by the LRT.')

    message(
      '  LRT: ', length(scoring_genes),
      ' genes passed padj < ', lrt_padj_cutoff,
      ' with direction = ', lrt_keep_direction,
      if (use_lrt_tpm_filter) {
        if (identical(lrt_tpm_filter_mode, 'min_samples')) {
          paste0(
            ', TPM >= ', format(lrt_tpm_cutoff, trim = TRUE),
            ' in >= ', lrt_tpm_min_samples, ' LRT samples'
          )
        } else {
          paste0(
            ', mean TPM >= ', format(lrt_tpm_cutoff, trim = TRUE),
            ' at >= 1 LRT timepoint'
          )
        }
      } else {
        ''
      },
      ' (out of ', nrow(lrt_results), ' tested)'
    )
  } else {
    scoring_genes <- intersect(scoring_gene_ids, rownames(vst_mat))
    if (length(scoring_genes) == 0L) {
      stop('none of scoring_gene_ids were found in the VST matrix.')
    }
    missing_scoring_gene_ids <- setdiff(scoring_gene_ids, scoring_genes)
    if (length(missing_scoring_gene_ids) > 0L) {
      message(
        '  Preselected gene set: skipped ',
        length(missing_scoring_gene_ids),
        ' genes not present in the VST matrix.'
      )
    }
    lrt_results <- data.frame(
      gene_id = scoring_genes,
      baseMean = NA_real_,
      log2FoldChange = NA_real_,
      lfcSE = NA_real_,
      stat = NA_real_,
      pvalue = NA_real_,
      padj = NA_real_,
      pass_lrt_padj = NA,
      pass_lrt_direction = NA,
      pass_lrt_expression = TRUE,
      pass_lrt_gene_filter = TRUE,
      stringsAsFactors = FALSE
    )
    message(
      '  Preselected gene set: using ',
      length(scoring_genes),
      ' genes directly for PCA training/scoring.'
    )
  }

  # ================================================================
  # step 4: train PCA on ref samples, compute maturation scores
  # ================================================================
  log_step(4L, n_steps, 'training PCA and computing maturation scores')

  analysis_vst <- vst_mat[scoring_genes, analysis_ids, drop = FALSE]
  ref_vst <- analysis_vst[, ref_ids, drop = FALSE]
  if (isTRUE(scoring_pca_scale)) {
    ref_sd <- apply(ref_vst, 1L, stats::sd)
    keep_sd <- ref_sd > 0
    if (any(!keep_sd)) {
      message(
        '  Dropping ',
        sum(!keep_sd),
        ' zero-variance training genes before scaled PCA.'
      )
      ref_vst <- ref_vst[keep_sd, , drop = FALSE]
      analysis_vst <- analysis_vst[keep_sd, , drop = FALSE]
      scoring_genes <- scoring_genes[keep_sd]
    }
  }
  pca_fit <- stats::prcomp(t(ref_vst), center = TRUE, scale. = scoring_pca_scale)
  pca_importance <- summary(pca_fit)[['importance']]
  variance_pc_names <- intersect(
    paste0('PC', seq_len(max(5L, scoring_n_pcs))),
    colnames(pca_importance)
  )
  pca_variance_explained <- pca_importance[
    'Proportion of Variance',
    variance_pc_names,
    drop = TRUE
  ]
  pc_cols <- paste0('PC', seq_len(scoring_n_pcs))
  analysis_coords <- project_into_pca_space(
    analysis_vst,
    pca_fit,
    n_pcs = scoring_n_pcs
  )
  analysis_coords <- analysis_coords[
    analysis_meta[[sample_id_col_name]],
    pc_cols,
    drop = FALSE
  ]

  # build scored sample table w/ standardized internal cols
  score_df <- cbind(analysis_meta, analysis_coords)
  rownames(score_df) <- score_df[[sample_id_col_name]]
  score_df[['timepoint_numeric']] <- score_df[[timepoint_numeric_col_name]]
  score_df[['group_value']] <- as.character(score_df[[group_col_name]])
  score_df[['sample_set']] <- ifelse(
    as.character(score_df[[sample_id_col_name]]) %in% ref_ids,
    'reference', 'projection'
  )

  centroid_df <- compute_timepoint_centroids(
    score_df,
    ref_timepoints,
    pc_cols = pc_cols
  )
  vec_proj <- project_onto_scoring_vector(
    score_df, centroid_df,
    vector_start_timepoint, vector_end_timepoint,
    scoring_mode,
    pc_cols = pc_cols
  )

  score_df[['maturation_score']] <- vec_proj[['maturation_score']]
  for (pc_idx in seq_along(pc_cols)) {
    score_df[[paste0('vector_', pc_cols[[pc_idx]])]] <- vec_proj[['projected_coords']][, pc_idx]
  }
  score_df[['scoring_start']] <- vector_start_timepoint
  score_df[['scoring_end']] <- vector_end_timepoint
  score_df[['scoring_n_pcs']] <- scoring_n_pcs
  score_df[['lrt_padj_cutoff']] <- lrt_padj_cutoff
  score_df[['lrt_tpm_cutoff']] <- if (use_lrt_tpm_filter) lrt_tpm_cutoff else NA_real_
  score_df[['lrt_tpm_filter_mode']] <- if (use_lrt_tpm_filter) lrt_tpm_filter_mode else NA_character_
  score_df[['lrt_tpm_min_samples']] <- if (use_lrt_tpm_filter) lrt_tpm_min_samples else NA_integer_
  score_df[['lrt_keep_direction']] <- lrt_keep_direction
  score_df[['n_lrt_genes']] <- length(scoring_genes)
  score_df[['scoring_mode']] <- scoring_mode

  # per-tp Welch t-tests (optional)
  group_levels <- unique(score_df[['group_value']])
  t_test_results <- NULL

  if (!is.null(welch_t_test_comparison)) {
    validated_pairs <- validate_t_test_pairs(
      welch_t_test_comparison, group_levels
    )
    t_test_results <- run_welch_t_tests(
      score_df, validated_pairs
    )

    # optional BH correction within each comparison pair
    t_test_results[['p_adjusted']] <- NA_real_
    ok_idx <- t_test_results[['status']] == 'ok'
    if (isTRUE(welch_t_test_bh_adjust)) {
      comps <- unique(
        t_test_results[['comparison']][ok_idx]
      )
      for (comp in comps) {
        rows <- ok_idx &
          t_test_results[['comparison']] == comp
        t_test_results[['p_adjusted']][rows] <-
          stats::p.adjust(
            t_test_results[['p_value']][rows],
            method = 'BH'
          )
      }
    }

    ok_rows <- t_test_results[ok_idx, , drop = FALSE]
    if (nrow(ok_rows) > 0L) {
      for (comp in unique(ok_rows[['comparison']])) {
        comp_rows <- ok_rows[ok_rows[['comparison']] == comp, , drop = FALSE]
        tp_p <- vapply(seq_len(nrow(comp_rows)), function(r) {
          sprintf('D%g p=%s', comp_rows[['timepoint_numeric']][r], signif(comp_rows[['p_value']][r], 2))
        }, character(1))
        message('  Welch: ', comp, ' -- ', paste(tp_p, collapse = ', '))
      }
    }
  }

  # filter to display_timepoints for boxplot if set
  plotted_score_df <- score_df
  plotted_t_test <- t_test_results
  if (length(display_timepoints) > 0L) {
    tp_col <- 'timepoint_numeric'
    in_display <- function(df) {
      df[[tp_col]] %in% display_timepoints
    }
    missing_tp <- setdiff(
      display_timepoints,
      unique(score_df[[tp_col]])
    )
    if (length(missing_tp) > 0L) {
      stop(
        'display_timepoints not in scored data: ',
        paste(missing_tp, collapse = ', ')
      )
    }
    plotted_score_df <- score_df[
      in_display(score_df), ,
      drop = FALSE
    ]
    if (!is.null(t_test_results)) {
      plotted_t_test <- t_test_results[
        in_display(t_test_results), ,
        drop = FALSE
      ]
    }
  }

  # ================================================================
  # step 5: save outputs
  # ================================================================
  log_step(
    5L, n_steps,
    paste0(
      'saving outputs to ',
      normalizePath(output_dir, mustWork = FALSE)
    )
  )

  paths <- list(
    score_csv = paste0(output_file_path, '_maturation_scores.csv'),
    lrt_genes_csv = paste0(output_file_path, '_lrt_genes.csv'),
    boxplot_svg = paste0(boxplot_output_path, '_score_boxplot.svg')
  )

  utils::write.csv(
    score_df, paths[['score_csv']],
    row.names = FALSE
  )
  log_output('maturation scores', paths[['score_csv']])

  utils::write.csv(
    data.frame(
      gene_id = scoring_genes,
      stringsAsFactors = FALSE
    ),
    paths[['lrt_genes_csv']],
    row.names = FALSE
  )
  log_output('lrt gene list', paths[['lrt_genes_csv']])

  boxplot_gg <- build_and_save_boxplot(
    score_df = plotted_score_df,
    t_test_results = plotted_t_test,
    plot_title = boxplot_title,
    output_path = paths[['boxplot_svg']],
    group_colors = group_colors,
    group_order = boxplot_group_order,
    group_x_offsets = boxplot_group_x_offsets,
    point_size = boxplot_point_size,
    point_stroke = boxplot_point_stroke,
    point_alpha = boxplot_point_alpha,
    y_limits = boxplot_y_limits,
    y_breaks = boxplot_y_breaks,
    plot_style = boxplot_plot_style,
    show_legend = boxplot_show_legend,
    connect_group_medians = boxplot_connect_group_medians,
    show_dotplot_mean_bar = boxplot_show_dotplot_mean_bar,
    point_jitter_width = boxplot_point_jitter_width,
    point_jitter_height = boxplot_point_jitter_height,
    timepoint_gap = boxplot_timepoint_gap,
    box_width = boxplot_box_width,
    plot_width = boxplot_width,
    plot_height = boxplot_height,
    force_panel_height_in = boxplot_force_panel_height_in,
    axis_title_size = boxplot_axis_title_size,
    axis_text_size = boxplot_axis_text_size,
    reference_line_width = boxplot_reference_line_width,
    axis_line_width = boxplot_axis_line_width,
    range_line_width = boxplot_range_line_width,
    summary_line_width = boxplot_summary_line_width,
    show_y_grid = boxplot_show_y_grid,
    plot_margin = boxplot_plot_margin,
    show_significance = show_significance,
    group_shapes = boxplot_group_shapes,
    group_linetypes = boxplot_group_linetypes,
    proportional_timepoint_spacing = boxplot_proportional_timepoint_spacing,
    y_axis_label = boxplot_y_axis_label,
    x_axis_tick_prefix = boxplot_x_axis_tick_prefix,
    y_axis_label_digits = boxplot_y_axis_label_digits,
    axis_tick_length = boxplot_axis_tick_length,
    barplot_outline_colour = boxplot_barplot_outline_colour,
    median_line_outline_width = boxplot_median_line_outline_width,
    group_line_colors = boxplot_group_line_colors
  )
  log_output(
    if (identical(boxplot_plot_style, 'barplot')) {
      'score barplot'
    } else if (identical(boxplot_plot_style, 'lineplot')) {
      'score lineplot'
    } else if (identical(boxplot_plot_style, 'boxplot')) {
      'score boxplot'
    } else {
      'score dotplot'
    },
    paths[['boxplot_svg']]
  )

  log_step(5L, n_steps, 'done')

  invisible(list(
    lrt_results = lrt_results,
    score_table = score_df,
    t_test_results = t_test_results,
    centroid_table = centroid_df,
    pca_fit = pca_fit,
    pca_variance_explained = pca_variance_explained,
    scoring_pc_columns = pc_cols,
    vector_unit = vec_proj[['vector_unit']],
    vector_length = vec_proj[['vector_length']],
    output_paths = paths,
    vst_mat = vst_mat,
    scoring_genes = scoring_genes,
    boxplot_plot = boxplot_gg
  ))
}
