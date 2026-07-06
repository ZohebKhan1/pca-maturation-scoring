# Created:
# 2026-07-06
#
# Inputs:
# - expression_matrix: numeric matrix with genes as rows and samples as columns
# - metadata: data.frame with sample IDs, reference labels, and numeric timepoints
# - temporal_genes: character vector of genes used to train the timing axis
#
# Outputs:
# - score_differentiation_timing(): list with sample scores, fitted PCA object,
#   reference day-centroid polyline, and retained temporal genes
#
# Purpose:
# Score samples along a reference differentiation trajectory by projecting them
# into a PCA space trained on temporal genes and then onto the ordered
# day-centroid polyline. This function only performs polyline scoring.

score_differentiation_timing <- function(expression_matrix,
                                         metadata,
                                         temporal_genes,
                                         sample_id_col = 'sample_id',
                                         time_col = 'day_numeric',
                                         reference_col = NULL,
                                         reference_values = NULL,
                                         n_pcs = 3L) {
  if (!base::is.matrix(expression_matrix) && !base::is.data.frame(expression_matrix)) {
    base::stop('expression_matrix must be a numeric matrix or data.frame.', call. = FALSE)
  }
  expression_matrix <- base::as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- 'numeric'

  required_cols <- base::c(sample_id_col, time_col)
  missing_cols <- base::setdiff(required_cols, base::colnames(metadata))
  if (base::length(missing_cols) > 0L) {
    base::stop('metadata is missing required columns: ', base::paste(missing_cols, collapse = ', '), call. = FALSE)
  }
  if (base::anyDuplicated(metadata[[sample_id_col]]) > 0L) {
    base::stop('metadata sample IDs must be unique.', call. = FALSE)
  }
  if (base::anyDuplicated(base::colnames(expression_matrix)) > 0L) {
    base::stop('expression_matrix column names must be unique sample IDs.', call. = FALSE)
  }

  metadata <- metadata[base::match(base::colnames(expression_matrix), metadata[[sample_id_col]]), , drop = FALSE]
  if (base::anyNA(metadata[[sample_id_col]])) {
    base::stop('metadata must contain every expression_matrix sample column.', call. = FALSE)
  }
  if (!base::is.numeric(metadata[[time_col]])) {
    base::stop('time_col must be numeric.', call. = FALSE)
  }

  temporal_genes <- base::intersect(base::unique(temporal_genes), base::rownames(expression_matrix))
  if (base::length(temporal_genes) < 2L) {
    base::stop('At least two temporal genes must be present in expression_matrix.', call. = FALSE)
  }

  if (base::is.null(reference_col)) {
    reference_idx <- base::rep(TRUE, base::nrow(metadata))
  } else {
    if (!reference_col %in% base::colnames(metadata)) {
      base::stop('reference_col is not present in metadata.', call. = FALSE)
    }
    if (base::is.null(reference_values)) {
      base::stop('reference_values must be supplied when reference_col is used.', call. = FALSE)
    }
    reference_idx <- metadata[[reference_col]] %in% reference_values
  }
  if (base::sum(reference_idx) < 3L) {
    base::stop('At least three reference samples are required.', call. = FALSE)
  }

  reference_days <- base::sort(base::unique(metadata[[time_col]][reference_idx]))
  if (base::length(reference_days) < 2L) {
    base::stop('Reference samples must span at least two timepoints.', call. = FALSE)
  }
  missing_reference_days <- reference_days[
    !base::vapply(reference_days, function(day_value) {
      base::sum(reference_idx & metadata[[time_col]] == day_value) > 0L
    }, base::logical(1))
  ]
  if (base::length(missing_reference_days) > 0L) {
    base::stop('Reference centroid cannot be computed for every reference day.', call. = FALSE)
  }

  n_pcs <- base::as.integer(n_pcs[[1]])
  max_pcs <- base::min(base::sum(reference_idx) - 1L, base::length(temporal_genes))
  if (!base::is.finite(n_pcs) || n_pcs < 1L || n_pcs > max_pcs) {
    base::stop('n_pcs must be between 1 and ', max_pcs, ' for these inputs.', call. = FALSE)
  }

  reference_matrix <- base::t(expression_matrix[temporal_genes, reference_idx, drop = FALSE])
  pca_fit <- stats::prcomp(reference_matrix, center = TRUE, scale. = FALSE)
  rotation <- pca_fit$rotation[, base::seq_len(n_pcs), drop = FALSE]

  centered_matrix <- base::scale(
    base::t(expression_matrix[temporal_genes, , drop = FALSE]),
    center = pca_fit$center,
    scale = FALSE
  )
  pca_coordinates <- base::as.data.frame(centered_matrix %*% rotation)
  base::colnames(pca_coordinates) <- base::paste0('PC', base::seq_len(n_pcs))
  pca_coordinates[[sample_id_col]] <- base::rownames(pca_coordinates)
  pca_coordinates[[time_col]] <- metadata[[time_col]]
  pca_coordinates[['is_reference']] <- reference_idx

  pc_cols <- base::paste0('PC', base::seq_len(n_pcs))
  centroid_polyline <- stats::aggregate(
    pca_coordinates[reference_idx, pc_cols, drop = FALSE],
    by = base::list(timepoint = metadata[[time_col]][reference_idx]),
    FUN = base::mean
  )
  centroid_polyline <- centroid_polyline[base::order(centroid_polyline$timepoint), , drop = FALSE]

  projection <- .project_to_centroid_polyline(
    point_matrix = base::as.matrix(pca_coordinates[, pc_cols, drop = FALSE]),
    centroid_polyline = centroid_polyline,
    pc_cols = pc_cols
  )

  first_day <- centroid_polyline$timepoint[[1]]
  last_day <- centroid_polyline$timepoint[[base::nrow(centroid_polyline)]]
  day_span <- last_day - first_day
  if (!base::is.finite(day_span) || day_span <= 0) {
    base::stop('Reference timepoints must increase.', call. = FALSE)
  }

  scores <- base::data.frame(
    sample_id = metadata[[sample_id_col]],
    observed_time = metadata[[time_col]],
    is_reference = reference_idx,
    predicted_time = projection$predicted_time,
    differentiation_score = (projection$predicted_time - first_day) / day_span,
    nearest_segment_start = projection$segment_start,
    nearest_segment_end = projection$segment_end,
    segment_fraction = projection$segment_fraction,
    squared_distance = projection$squared_distance,
    stringsAsFactors = FALSE
  )

  base::list(
    scores = scores,
    pca_coordinates = pca_coordinates,
    pca_fit = pca_fit,
    centroid_polyline = centroid_polyline,
    temporal_genes = temporal_genes,
    reference_time_range = base::c(start = first_day, end = last_day)
  )
}

.project_to_centroid_polyline <- function(point_matrix, centroid_polyline, pc_cols) {
  if (base::nrow(centroid_polyline) < 2L) {
    base::stop('centroid_polyline must contain at least two timepoints.', call. = FALSE)
  }

  start_points <- base::as.matrix(centroid_polyline[-base::nrow(centroid_polyline), pc_cols, drop = FALSE])
  end_points <- base::as.matrix(centroid_polyline[-1L, pc_cols, drop = FALSE])
  segment_vectors <- end_points - start_points
  segment_lengths_sq <- base::rowSums(segment_vectors^2)
  keep_segments <- base::is.finite(segment_lengths_sq) & segment_lengths_sq > 0
  if (!base::any(keep_segments)) {
    base::stop('The centroid polyline has no nonzero-length segments.', call. = FALSE)
  }

  start_points <- start_points[keep_segments, , drop = FALSE]
  segment_vectors <- segment_vectors[keep_segments, , drop = FALSE]
  segment_lengths_sq <- segment_lengths_sq[keep_segments]
  segment_start <- centroid_polyline$timepoint[-base::nrow(centroid_polyline)][keep_segments]
  segment_end <- centroid_polyline$timepoint[-1L][keep_segments]

  projected <- base::lapply(base::seq_len(base::nrow(point_matrix)), function(i) {
    point <- point_matrix[i, ]
    point_repeated <- base::matrix(
      point,
      nrow = base::nrow(start_points),
      ncol = base::ncol(start_points),
      byrow = TRUE
    )
    segment_offset <- point_repeated - start_points
    segment_fraction <- base::rowSums(segment_offset * segment_vectors) / segment_lengths_sq
    segment_fraction <- base::pmin(base::pmax(segment_fraction, 0), 1)
    projected_points <- start_points + segment_fraction * segment_vectors
    squared_distance <- base::rowSums((point_repeated - projected_points)^2)
    best_segment <- base::which.min(squared_distance)
    predicted_time <- segment_start[[best_segment]] +
      segment_fraction[[best_segment]] *
        (segment_end[[best_segment]] - segment_start[[best_segment]])

    base::data.frame(
      predicted_time = predicted_time,
      segment_start = segment_start[[best_segment]],
      segment_end = segment_end[[best_segment]],
      segment_fraction = segment_fraction[[best_segment]],
      squared_distance = squared_distance[[best_segment]]
    )
  })

  base::do.call(rbind, projected)
}
