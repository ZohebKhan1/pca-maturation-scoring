# Author: Zoheb Khan
# Date: 2026-05-25
#
# Input file paths:
# - data/GSE122380_metadata.rds - Sample metadata
# - data/GSE122380_counts.rds - Raw count matrix
# - data/GSE122380_vst.rds - VST expression matrix
#
# Output file paths:
# - tmp/GSE122380_leave_one_line_out_temporal_trajectory.rds
# - tmp/GSE122380_leave_one_line_out_temporal_genes.rds
#
# Goal:
# Leave each cell line out of the LRT and temporal trajectory definition, then
# project the held-out line onto the ordered day-centroid polyline.
#
# Notes:
# This maintained optional validation script checkpoints each held-out line's
# temporal gene list as soon as it is generated, so interrupted all-line runs
# can resume without repeating LRTs.

# 0.0 set seed -----------------

base::set.seed(2026)

# 0.1 load required packages -----------------

required_packages <- base::c(
  'DESeq2',
  'edgeR'
)

missing_packages <- required_packages[
  !base::vapply(required_packages, requireNamespace, base::logical(1), quietly = TRUE)
]
if (base::length(missing_packages) > 0L) {
  base::stop(
    'Missing required packages: ',
    base::paste(missing_packages, collapse = ', '),
    call. = FALSE
  )
}

# 0.2 load data -----------------

metadata_path <- base::file.path('data', 'GSE122380_metadata.rds')
counts_path <- base::file.path('data', 'GSE122380_counts.rds')
vst_path <- base::file.path('data', 'GSE122380_vst.rds')
loo_cache_path <- base::file.path(
  'tmp',
  'GSE122380_leave_one_line_out_temporal_trajectory.rds'
)
temporal_gene_cache_path <- base::file.path(
  'tmp',
  'GSE122380_leave_one_line_out_temporal_genes.rds'
)

metadata <- base::readRDS(metadata_path)
counts <- base::readRDS(counts_path)
vst <- base::readRDS(vst_path)

base::stopifnot(!base::anyDuplicated(metadata$sample_id))
base::stopifnot(base::all(metadata$sample_id %in% base::colnames(counts)))
base::stopifnot(base::all(metadata$sample_id %in% base::colnames(vst)))

counts <- counts[, metadata$sample_id, drop = FALSE]
vst <- vst[, metadata$sample_id, drop = FALSE]
base::stopifnot(base::identical(base::as.character(metadata$sample_id), base::colnames(counts)))
base::stopifnot(base::identical(base::as.character(metadata$sample_id), base::colnames(vst)))

metadata$day_factor <- base::factor(
  metadata$day_numeric,
  levels = base::sort(base::unique(metadata$day_numeric))
)
metadata$cell_line <- base::droplevels(base::factor(metadata$cell_line))

# 0.3 define params -----------------

expression_cpm_cutoff <- 10
vst_dynamic_range_cutoff <- 0.6
lrt_padj_cutoff <- 1e-7
ordered_days <- base::sort(base::unique(metadata$day_numeric))
cell_lines <- base::sort(base::unique(base::as.character(metadata$cell_line)))
selected_cell_lines <- cell_lines
required_score_methods <- 'Polyline'
cache_key <- base::list(
  expression_cpm_cutoff = expression_cpm_cutoff,
  vst_dynamic_range_cutoff = vst_dynamic_range_cutoff,
  lrt_padj_cutoff = lrt_padj_cutoff,
  heldout_lines = selected_cell_lines,
  lrt_gene_set = 'genes passing training-set day-mean TMM CPM, LRT padj, and VST dynamic-range filters'
)

# 1.0 helper functions -----------------

read_temporal_gene_cache <- function() {
  if (base::file.exists(temporal_gene_cache_path)) {
    cache <- base::readRDS(temporal_gene_cache_path)
    if (!base::isTRUE(base::identical(cache$parameters, cache_key))) {
      cache <- base::list(
        temporal_genes = base::list(),
        lrt_results = base::list()
      )
    }
  } else {
    cache <- base::list(
      temporal_genes = base::list(),
      lrt_results = base::list()
    )
  }

  cache$parameters <- cache_key
  cache
}

write_temporal_gene_cache <- function(cache) {
  base::dir.create(base::dirname(temporal_gene_cache_path), recursive = TRUE, showWarnings = FALSE)
  base::saveRDS(cache, temporal_gene_cache_path)
}

read_trajectory_cache <- function() {
  if (!base::file.exists(loo_cache_path)) {
    return(NULL)
  }
  cache <- base::readRDS(loo_cache_path)
  if (!base::isTRUE(base::identical(cache$parameters, cache_key))) {
    return(NULL)
  }
  cache
}

write_trajectory_cache <- function(results) {
  score_table <- base::do.call(
    rbind,
    base::lapply(results, function(result) result$scores)
  )
  base::rownames(score_table) <- NULL

  summary_table <- base::do.call(
    rbind,
    base::lapply(results, function(result) {
      pca_summary <- result$pca_summary
      pca_summary$heldout_line <- result$heldout_line
      pca_summary$temporal_genes <- base::length(result$temporal_genes)
      pca_summary[
        ,
        base::c(
          'heldout_line',
          'gene_set',
          'gene_count',
          'temporal_genes',
          'pc1_percent',
          'pc2_percent',
          'pc3_percent'
        ),
        drop = FALSE
      ]
    })
  )
  base::rownames(summary_table) <- NULL

  loo_cache <- base::list(
    scores = score_table,
    summary = summary_table,
    results = results,
    selected_cell_lines = base::names(results),
    temporal_gene_cache_path = temporal_gene_cache_path,
    parameters = cache_key,
    notes = base::list(
      gene_set_split = paste(
        'Maturation genes have positive training-set Spearman correlation',
        'with day; progenitor genes have negative training-set Spearman correlation with day.'
      ),
      method_note = paste(
        'Polyline scores project held-out samples onto the nearest segment of',
        'the ordered day-centroid polyline.'
      )
    )
  )

  base::dir.create(base::dirname(loo_cache_path), recursive = TRUE, showWarnings = FALSE)
  base::saveRDS(loo_cache, loo_cache_path)
  loo_cache
}

get_training_expression_filter <- function(training_metadata) {
  training_counts <- counts[, training_metadata$sample_id, drop = FALSE]
  edge_dge <- edgeR::DGEList(counts = training_counts)
  edge_dge <- edgeR::calcNormFactors(edge_dge, method = 'TMM')
  tmm_cpm <- edgeR::cpm(edge_dge, normalized.lib.sizes = TRUE)

  mean_tmm_cpm_by_day <- base::sapply(ordered_days, function(day_value) {
    day_samples <- training_metadata$sample_id[training_metadata$day_numeric == day_value]
    base::rowMeans(tmm_cpm[, day_samples, drop = FALSE], na.rm = TRUE)
  })
  max_stage_mean_cpm <- base::apply(mean_tmm_cpm_by_day, 1, base::max)

  max_stage_mean_cpm >= expression_cpm_cutoff
}

filter_by_training_vst_dynamic_range <- function(gene_ids, training_metadata) {
  gene_day_means <- base::sapply(ordered_days, function(day_value) {
    day_samples <- training_metadata$sample_id[training_metadata$day_numeric == day_value]
    base::rowMeans(vst[gene_ids, day_samples, drop = FALSE], na.rm = TRUE)
  })
  dynamic_range <- base::apply(gene_day_means, 1, function(gene_values) {
    base::max(gene_values, na.rm = TRUE) - base::min(gene_values, na.rm = TRUE)
  })
  gene_ids[dynamic_range[gene_ids] >= vst_dynamic_range_cutoff]
}

run_training_lrt <- function(training_metadata, expression_pass) {
  filtered_genes <- base::names(expression_pass)[expression_pass]
  base::stopifnot(base::length(filtered_genes) > 1L)

  filtered_counts <- counts[filtered_genes, training_metadata$sample_id, drop = FALSE]
  lrt_metadata <- training_metadata
  base::rownames(lrt_metadata) <- lrt_metadata$sample_id

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = lrt_metadata,
    design = ~ cell_line + day_factor
  )
  dds <- DESeq2::DESeq(
    dds,
    test = 'LRT',
    reduced = ~cell_line,
    quiet = TRUE
  )

  lrt_results <- base::as.data.frame(DESeq2::results(dds, alpha = lrt_padj_cutoff))
  lrt_results$gene_id <- base::rownames(lrt_results)
  lrt_results[
    base::order(lrt_results$padj, lrt_results$pvalue, na.last = TRUE), ,
    drop = FALSE
  ]
}

project_with_pca <- function(pca_fit, projected_matrix) {
  centered_matrix <- base::scale(
    base::t(projected_matrix),
    center = pca_fit$center,
    scale = FALSE
  )
  projected_scores <- centered_matrix %*% pca_fit$rotation[, base::seq_len(3), drop = FALSE]
  base::as.data.frame(projected_scores)
}

calculate_centroid_polyline <- function(training_pca_df) {
  centroids <- stats::aggregate(
    base::cbind(PC1, PC2, PC3) ~ day_numeric,
    data = training_pca_df,
    FUN = base::mean
  )
  centroids[base::order(centroids$day_numeric), , drop = FALSE]
}

predict_day_polyline <- function(training_pca_df, heldout_pca_df) {
  centroids <- calculate_centroid_polyline(training_pca_df)
  pc_cols <- base::c('PC1', 'PC2', 'PC3')
  start_points <- base::as.matrix(centroids[-base::nrow(centroids), pc_cols, drop = FALSE])
  end_points <- base::as.matrix(centroids[-1L, pc_cols, drop = FALSE])
  segment_vectors <- end_points - start_points
  segment_lengths_sq <- base::rowSums(segment_vectors^2)
  keep_segments <- base::is.finite(segment_lengths_sq) & segment_lengths_sq > 0
  start_points <- start_points[keep_segments, , drop = FALSE]
  end_points <- end_points[keep_segments, , drop = FALSE]
  segment_vectors <- segment_vectors[keep_segments, , drop = FALSE]
  segment_lengths_sq <- segment_lengths_sq[keep_segments]
  segment_start_days <- centroids$day_numeric[-base::nrow(centroids)][keep_segments]
  segment_end_days <- centroids$day_numeric[-1L][keep_segments]

  heldout_points <- base::as.matrix(heldout_pca_df[, pc_cols])
  base::vapply(base::seq_len(base::nrow(heldout_points)), function(i) {
    point <- heldout_points[i, ]
    point_repeated <- base::matrix(
      point,
      nrow = base::nrow(start_points),
      ncol = base::length(pc_cols),
      byrow = TRUE
    )
    segment_offset <- point_repeated - start_points
    segment_fraction <- base::rowSums(segment_offset * segment_vectors) / segment_lengths_sq
    segment_fraction <- base::pmin(base::pmax(segment_fraction, 0), 1)
    projected_points <- start_points + segment_fraction * segment_vectors
    squared_distance <- base::rowSums((point_repeated - projected_points)^2)
    best_segment <- base::which.min(squared_distance)
    segment_start_days[[best_segment]] +
      segment_fraction[[best_segment]] *
        (segment_end_days[[best_segment]] - segment_start_days[[best_segment]])
  }, base::numeric(1))
}

get_cached_lrt_result <- function(heldout_line) {
  temporal_gene_cache <- read_temporal_gene_cache()
  if (!base::is.null(temporal_gene_cache$temporal_genes[[heldout_line]]) &&
    !base::is.null(temporal_gene_cache$lrt_results[[heldout_line]])) {
    return(base::list(
      temporal_genes = temporal_gene_cache$temporal_genes[[heldout_line]],
      lrt_results = temporal_gene_cache$lrt_results[[heldout_line]]
    ))
  }

  existing_cache <- read_trajectory_cache()
  if (base::is.null(existing_cache)) {
    return(NULL)
  }
  existing_result <- existing_cache$results[[heldout_line]]
  if (base::is.null(existing_result$lrt_results) || base::is.null(existing_result$temporal_genes)) {
    return(NULL)
  }
  existing_result
}

split_temporal_genes <- function(temporal_genes, training_metadata) {
  gene_day_means <- base::sapply(ordered_days, function(day_value) {
    day_samples <- training_metadata$sample_id[training_metadata$day_numeric == day_value]
    base::rowMeans(vst[temporal_genes, day_samples, drop = FALSE], na.rm = TRUE)
  })
  gene_day_cor <- base::apply(gene_day_means, 1, function(gene_values) {
    stats::cor(ordered_days, gene_values, method = 'spearman', use = 'pairwise.complete.obs')
  })
  gene_day_cor[!base::is.finite(gene_day_cor)] <- 0

  gene_sets <- base::list(
    `All temporal` = temporal_genes,
    Maturation = base::names(gene_day_cor)[gene_day_cor > 0],
    Progenitor = base::names(gene_day_cor)[gene_day_cor < 0]
  )
  base::stopifnot(base::all(base::vapply(gene_sets, base::length, base::integer(1)) > 1L))

  gene_sets
}

score_gene_set <- function(gene_ids, gene_set_name, heldout_line, training_metadata, heldout_metadata) {
  pca_fit <- stats::prcomp(
    base::t(vst[gene_ids, training_metadata$sample_id, drop = FALSE]),
    center = TRUE,
    scale. = FALSE
  )
  pca_var <- base::round(base::summary(pca_fit)$importance[2, 1:3] * 100, 2)

  training_scores <- base::as.data.frame(pca_fit$x[, base::seq_len(3), drop = FALSE])
  training_scores$sample_id <- base::rownames(training_scores)
  training_scores$day_numeric <- training_metadata$day_numeric[
    base::match(training_scores$sample_id, training_metadata$sample_id)
  ]
  training_scores$cell_line <- training_metadata$cell_line[
    base::match(training_scores$sample_id, training_metadata$sample_id)
  ]

  heldout_scores <- project_with_pca(
    pca_fit = pca_fit,
    projected_matrix = vst[gene_ids, heldout_metadata$sample_id, drop = FALSE]
  )
  heldout_scores$sample_id <- base::rownames(heldout_scores)
  heldout_scores$day_numeric <- heldout_metadata$day_numeric[
    base::match(heldout_scores$sample_id, heldout_metadata$sample_id)
  ]
  heldout_scores$cell_line <- heldout_line

  polyline_predicted_day <- predict_day_polyline(training_scores, heldout_scores)

  score_df <- base::data.frame(
    heldout_line = heldout_line,
    sample_id = heldout_scores$sample_id,
    cell_line = heldout_line,
    actual_day = heldout_scores$day_numeric,
    gene_set = gene_set_name,
    method = 'Polyline',
    predicted_day = polyline_predicted_day,
    residual = polyline_predicted_day - heldout_scores$day_numeric,
    stringsAsFactors = FALSE
  )

  base::list(
    gene_set = gene_set_name,
    gene_count = base::length(gene_ids),
    pca_variance_percent = pca_var,
    scores = score_df
  )
}

analyze_heldout_line <- function(heldout_line) {
  base::message('Analyzing held-out line: ', heldout_line)
  training_metadata <- metadata[metadata$cell_line != heldout_line, , drop = FALSE]
  heldout_metadata <- metadata[metadata$cell_line == heldout_line, , drop = FALSE]
  training_metadata$cell_line <- base::droplevels(training_metadata$cell_line)

  cached_result <- get_cached_lrt_result(heldout_line)
  if (base::is.null(cached_result)) {
    expression_pass <- get_training_expression_filter(training_metadata)
    lrt_results <- run_training_lrt(training_metadata, expression_pass)
    lrt_genes <- lrt_results$gene_id[
      !base::is.na(lrt_results$padj) & lrt_results$padj < lrt_padj_cutoff
    ]
    temporal_genes <- filter_by_training_vst_dynamic_range(
      gene_ids = lrt_genes,
      training_metadata = training_metadata
    )
    temporal_gene_cache <- read_temporal_gene_cache()
    temporal_gene_cache$temporal_genes[[heldout_line]] <- temporal_genes
    temporal_gene_cache$lrt_results[[heldout_line]] <- lrt_results
    temporal_gene_cache$updated_at <- base::format(base::Sys.time(), usetz = TRUE)
    write_temporal_gene_cache(temporal_gene_cache)
  } else {
    lrt_results <- cached_result$lrt_results
    temporal_genes <- cached_result$temporal_genes
  }
  base::stopifnot(base::length(temporal_genes) > 1L)

  gene_sets <- split_temporal_genes(temporal_genes, training_metadata)
  gene_set_results <- base::Map(
    score_gene_set,
    gene_sets,
    base::names(gene_sets),
    MoreArgs = base::list(
      heldout_line = heldout_line,
      training_metadata = training_metadata,
      heldout_metadata = heldout_metadata
    )
  )
  score_df <- base::do.call(rbind, base::lapply(gene_set_results, function(result) result$scores))
  pca_summary <- base::do.call(rbind, base::lapply(gene_set_results, function(result) {
    base::data.frame(
      gene_set = result$gene_set,
      gene_count = result$gene_count,
      pc1_percent = result$pca_variance_percent[[1]],
      pc2_percent = result$pca_variance_percent[[2]],
      pc3_percent = result$pca_variance_percent[[3]],
      stringsAsFactors = FALSE
    )
  }))

  base::list(
    heldout_line = heldout_line,
    temporal_genes = temporal_genes,
    lrt_results = lrt_results,
    gene_sets = gene_sets,
    pca_summary = pca_summary,
    scores = score_df
  )
}

# 2.0 run leave-one-line-out analysis -----------------

existing_trajectory_cache <- read_trajectory_cache()
if (base::is.null(existing_trajectory_cache)) {
  loo_results <- base::list()
} else {
  loo_results <- existing_trajectory_cache$results
  temporal_gene_cache <- read_temporal_gene_cache()
  for (heldout_line in base::intersect(base::names(loo_results), selected_cell_lines)) {
    if (base::is.null(temporal_gene_cache$temporal_genes[[heldout_line]]) &&
      !base::is.null(loo_results[[heldout_line]]$temporal_genes)) {
      temporal_gene_cache$temporal_genes[[heldout_line]] <- loo_results[[heldout_line]]$temporal_genes
      temporal_gene_cache$lrt_results[[heldout_line]] <- loo_results[[heldout_line]]$lrt_results
    }
  }
  temporal_gene_cache$updated_at <- base::format(base::Sys.time(), usetz = TRUE)
  write_temporal_gene_cache(temporal_gene_cache)
}

for (heldout_line in selected_cell_lines) {
  cached_methods <- if (base::is.null(loo_results[[heldout_line]])) {
    base::character(0)
  } else {
    base::unique(loo_results[[heldout_line]]$scores$method)
  }
  if (!base::is.null(loo_results[[heldout_line]]) &&
    base::all(required_score_methods %in% cached_methods)) {
    base::message('Using cached trajectory result for held-out line: ', heldout_line)
  } else {
    loo_results[[heldout_line]] <- analyze_heldout_line(heldout_line)
    loo_results <- loo_results[selected_cell_lines[selected_cell_lines %in% base::names(loo_results)]]
    loo_cache <- write_trajectory_cache(loo_results)
  }
}

loo_results <- loo_results[selected_cell_lines]
loo_cache <- write_trajectory_cache(loo_results)
summary_table <- loo_cache$summary

# 3.0 end-of-script summary -----------------

summary_table
base::c(
  loo_cache_path,
  temporal_gene_cache_path
)
