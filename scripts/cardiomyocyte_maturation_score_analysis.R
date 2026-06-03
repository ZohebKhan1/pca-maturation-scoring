# Created:
# 2026-05-25
#
# Inputs:
# - data/GSE122380_metadata.rds: QC-filtered GSE122380 sample metadata
# - data/GSE122380_counts.rds: filtered raw count matrix
# - data/GSE122380_vst.rds: VST expression matrix aligned to metadata
#
# Outputs:
# - R objects used by scripts/cardiomyocyte_maturation_score.Rmd to render the
#   bookdown site in docs/
#
# Purpose:
# Build the cardiomyocyte maturation scoring tutorial figures and summary
# tables from the active QC-processed GSE122380 inputs.
#
# Notes:
# This script is sourced by scripts/cardiomyocyte_maturation_score.Rmd and
# can also be run directly for validation. It does not write data, result, or
# figure files.

# 0.0 define local helpers -----------------

required_packages <- base::c(
  'DESeq2',
  'ComplexHeatmap',
  'circlize',
  'ggplot2',
  'gt',
  'patchwork',
  'plotly',
  'viridis'
)
missing_packages <- required_packages[
  !base::vapply(
    required_packages,
    base::requireNamespace,
    base::logical(1),
    quietly = TRUE
  )
]
if (base::length(missing_packages) > 0L) {
  base::stop('missing required packages: ', base::paste(missing_packages, collapse = ', '))
}

find_repo_root <- function() {
  candidates <- base::unique(base::normalizePath(
    base::c(base::getwd(), base::file.path(base::getwd(), '..'), base::file.path(base::getwd(), '../..')),
    mustWork = FALSE
  ))
  for (candidate in candidates) {
    if (base::file.exists(base::file.path(candidate, 'data', 'GSE122380_metadata.rds'))) {
      return(candidate)
    }
  }
  base::stop('Could not locate repository root containing data/GSE122380_metadata.rds.')
}

repo_root <- find_repo_root()
data_dir <- base::file.path(repo_root, 'data')
metadata_path <- base::file.path(data_dir, 'GSE122380_metadata.rds')
counts_path <- base::file.path(data_dir, 'GSE122380_counts.rds')
vst_path <- base::file.path(data_dir, 'GSE122380_vst.rds')

plot_font_family <- 'Inter'
lrt_padj_cutoff <- 0.05
n_heatmap_genes <- 1500L
n_temporal_clusters <- 4L
n_pc1_loading_genes_per_direction <- 10L

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = plot_font_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = 'black', family = plot_font_family),
      axis.text = ggplot2::element_text(color = 'black', family = plot_font_family),
      axis.title = ggplot2::element_text(color = 'black', family = plot_font_family),
      axis.line = ggplot2::element_line(color = 'black', linewidth = 0.3),
      axis.ticks = ggplot2::element_line(color = 'black', linewidth = 0.3),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(color = 'black', family = plot_font_family),
      legend.text = ggplot2::element_text(color = 'black', family = plot_font_family),
      strip.text = ggplot2::element_text(color = 'black', family = plot_font_family)
    )
}

label_ensembl_genes <- function(gene_ids) {
  clean_ids <- base::sub('\\..*$', '', gene_ids)

  if (!base::requireNamespace('AnnotationDbi', quietly = TRUE) ||
    !base::requireNamespace('org.Hs.eg.db', quietly = TRUE)) {
    return(gene_ids)
  }

  orgdb <- base::get('org.Hs.eg.db', envir = base::asNamespace('org.Hs.eg.db'))
  symbols <- AnnotationDbi::mapIds(
    orgdb,
    keys = clean_ids,
    keytype = 'ENSEMBL',
    column = 'SYMBOL',
    multiVals = 'first'
  )

  labels <- base::unname(symbols[clean_ids])
  missing_labels <- base::is.na(labels) | labels == ''
  labels[missing_labels] <- gene_ids[missing_labels]
  labels
}

orient_pca_fit_by_day <- function(pca_fit, sample_days) {
  pc1_day_cor <- stats::cor(
    pca_fit$x[, 'PC1'],
    sample_days[base::match(base::rownames(pca_fit$x), base::names(sample_days))],
    use = 'pairwise.complete.obs'
  )
  if (base::is.finite(pc1_day_cor) && pc1_day_cor < 0) {
    pca_fit$x[, 'PC1'] <- -pca_fit$x[, 'PC1']
    pca_fit$rotation[, 'PC1'] <- -pca_fit$rotation[, 'PC1']
  }
  pca_fit
}

get_pca_color_option <- function(set_label) {
  if (base::identical(set_label, 'C1+C2 (Early)')) {
    return('inferno')
  }
  if (base::identical(set_label, 'C3+C4 (Late)')) {
    return('mako')
  }
  'viridis'
}

format_lm_label <- function(x, y) {
  r_value <- stats::cor(x, y, use = 'pairwise.complete.obs')
  base::sprintf("r = %.2f\nR\u00b2 = %.2f", r_value, r_value^2)
}

# 1.0 load and validate inputs -----------------

base::message('CODEX_STEP load_data: loading active tutorial inputs')

meta <- base::readRDS(metadata_path)
counts <- base::readRDS(counts_path)
vst <- base::readRDS(vst_path)

base::stopifnot(!base::anyDuplicated(meta$sample_id))
base::stopifnot(base::all(meta$sample_id %in% base::colnames(counts)))
base::stopifnot(base::all(meta$sample_id %in% base::colnames(vst)))

counts <- counts[, meta$sample_id, drop = FALSE]
vst <- vst[, meta$sample_id, drop = FALSE]
base::stopifnot(base::identical(base::as.character(meta$sample_id), base::colnames(counts)))
base::stopifnot(base::identical(base::as.character(meta$sample_id), base::colnames(vst)))

meta$day_factor <- base::factor(
  meta$day_numeric,
  levels = base::sort(base::unique(meta$day_numeric))
)
meta$cell_line <- base::droplevels(base::factor(meta$cell_line))
sample_days <- stats::setNames(meta$day_numeric, meta$sample_id)

cohort_summary <- base::data.frame(
  samples = base::nrow(meta),
  cell_lines = base::length(base::unique(meta$cell_line)),
  genes = base::nrow(vst),
  first_day = base::min(meta$day_numeric),
  last_day = base::max(meta$day_numeric),
  stringsAsFactors = FALSE
)

samples_by_day <- base::as.data.frame(base::table(meta$day_numeric), stringsAsFactors = FALSE)
base::names(samples_by_day) <- base::c('day', 'samples')
samples_by_day$day <- base::as.integer(base::as.character(samples_by_day$day))

input_summary <- base::data.frame(
  timepoints = base::length(base::unique(meta$day_numeric)),
  samples_per_timepoint = base::paste0(
    base::min(samples_by_day$samples),
    '-',
    base::max(samples_by_day$samples)
  ),
  total_samples = base::nrow(meta),
  total_cell_lines = base::length(base::unique(meta$cell_line)),
  stringsAsFactors = FALSE
)

input_summary_table <- gt::gt(input_summary)
input_summary_table <- gt::cols_label(
  input_summary_table,
  timepoints = 'Timepoints',
  samples_per_timepoint = 'Samples per timepoint',
  total_samples = 'Total samples',
  total_cell_lines = 'Cell lines'
)
input_summary_table <- gt::fmt_integer(
  input_summary_table,
  columns = base::c('timepoints', 'total_samples', 'total_cell_lines')
)
input_summary_table <- gt::tab_options(
  input_summary_table,
  table.font.names = plot_font_family
)

cohort_summary_table <- gt::gt(cohort_summary)
cohort_summary_table <- gt::cols_label(
  cohort_summary_table,
  samples = 'Samples',
  cell_lines = 'Cell lines',
  genes = 'Genes',
  first_day = 'First day',
  last_day = 'Last day'
)
cohort_summary_table <- gt::fmt_integer(
  cohort_summary_table,
  columns = base::names(cohort_summary)
)
cohort_summary_table <- gt::tab_options(
  cohort_summary_table,
  table.font.names = plot_font_family
)

# 2.0 identify temporal genes -----------------

days <- base::sort(base::unique(meta$day_numeric))

lrt_meta <- meta
base::rownames(lrt_meta) <- lrt_meta$sample_id
lrt_meta$day_factor <- base::droplevels(lrt_meta$day_factor)
lrt_meta$cell_line <- base::droplevels(lrt_meta$cell_line)

base::message('CODEX_STEP lrt: selecting temporal genes with DESeq2 LRT')

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData = lrt_meta,
  design = ~ cell_line + day_factor
)

dds <- DESeq2::DESeq(
  dds,
  test = 'LRT',
  reduced = ~cell_line,
  quiet = TRUE
)

lrt_res <- base::as.data.frame(DESeq2::results(dds, alpha = lrt_padj_cutoff))
lrt_res$gene_id <- base::rownames(lrt_res)
lrt_res <- lrt_res[
  base::order(lrt_res$padj, lrt_res$pvalue, na.last = TRUE), ,
  drop = FALSE
]

sig_lrt <- lrt_res$gene_id[
  !base::is.na(lrt_res$padj) & lrt_res$padj < lrt_padj_cutoff
]

temporal_genes <- sig_lrt

lrt_summary <- base::data.frame(
  total_number_of_genes = base::nrow(lrt_res),
  lrt_significant_genes = base::length(sig_lrt),
  stringsAsFactors = FALSE
)

lrt_summary_table <- gt::gt(lrt_summary)
lrt_summary_table <- gt::cols_label(
  lrt_summary_table,
  total_number_of_genes = 'Total number of genes',
  lrt_significant_genes = 'Genes passing LRT filter of padj < 0.05'
)
lrt_summary_table <- gt::fmt_integer(
  lrt_summary_table,
  columns = base::names(lrt_summary)
)
lrt_summary_table <- gt::tab_options(
  lrt_summary_table,
  table.font.names = plot_font_family
)

# 4.0 create temporal heatmap and clusters -----------------

base::message('CODEX_STEP heatmap: clustering temporal gene trajectories')

base::set.seed(42)

heatmap_genes <- utils::head(temporal_genes, base::min(n_heatmap_genes, base::length(temporal_genes)))
hm_mat <- vst[heatmap_genes, meta$sample_id, drop = FALSE]
hm_z <- base::t(base::scale(base::t(hm_mat)))
hm_z[!base::is.finite(hm_z)] <- 0
hm_z[hm_z > 2] <- 2
hm_z[hm_z < -2] <- -2

order_one_day <- function(day_value) {
  day_samples <- meta$sample_id[meta$day_numeric == day_value]
  if (base::length(day_samples) <= 2L) {
    return(day_samples)
  }
  day_dist <- stats::dist(base::t(hm_z[, day_samples, drop = FALSE]))
  day_samples[stats::hclust(day_dist, method = 'average')$order]
}

column_order <- base::unlist(base::lapply(days, order_one_day), use.names = FALSE)
hm_z <- hm_z[, column_order, drop = FALSE]
hm_meta <- meta[base::match(column_order, meta$sample_id), , drop = FALSE]

row_km <- stats::kmeans(hm_z, centers = n_temporal_clusters, nstart = 30)
raw_cluster <- row_km$cluster

cluster_day_means <- base::do.call(rbind, base::lapply(base::sort(base::unique(raw_cluster)), function(cl) {
  genes <- base::names(raw_cluster)[raw_cluster == cl]
  day_values <- base::vapply(days, function(day) {
    samples <- hm_meta$sample_id[hm_meta$day_numeric == day]
    base::mean(hm_z[genes, samples, drop = FALSE], na.rm = TRUE)
  }, base::numeric(1))
  base::data.frame(
    raw_cluster = cl,
    peak_day = days[base::which.max(day_values)],
    mean_signal = base::max(day_values),
    stringsAsFactors = FALSE
  )
}))

cluster_order <- cluster_day_means$raw_cluster[
  base::order(cluster_day_means$peak_day, -cluster_day_means$mean_signal)
]
cluster_map <- stats::setNames(base::paste0('C', base::seq_along(cluster_order)), cluster_order)
gene_cluster <- base::factor(
  cluster_map[base::as.character(raw_cluster)],
  levels = base::paste0('C', base::seq_along(cluster_order))
)
base::names(gene_cluster) <- base::names(raw_cluster)

cluster_colors <- stats::setNames(
  base::c('#0072B2', '#009E73', '#E69F00', '#D55E00'),
  base::levels(gene_cluster)
)

col_fun <- circlize::colorRamp2(
  base::seq(-2, 2, length.out = 100),
  viridis::viridis(100)
)

cluster_annotation <- ComplexHeatmap::rowAnnotation(
  Cluster = gene_cluster[base::rownames(hm_z)],
  col = base::list(Cluster = cluster_colors),
  annotation_legend_param = base::list(
    Cluster = base::list(
      title_gp = grid::gpar(fontface = 'bold', fontfamily = plot_font_family),
      labels_gp = grid::gpar(fontfamily = plot_font_family)
    )
  ),
  show_annotation_name = FALSE,
  simple_anno_size = grid::unit(4, 'mm')
)

p_lrt_heatmap <- ComplexHeatmap::Heatmap(
  hm_z,
  name = 'Z-score',
  col = col_fun,
  left_annotation = cluster_annotation,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = gene_cluster[base::rownames(hm_z)],
  row_title = base::levels(gene_cluster),
  row_title_rot = 0,
  row_title_gp = grid::gpar(
    fontsize = 13,
    fontface = 'bold',
    fontfamily = plot_font_family,
    col = cluster_colors
  ),
  row_gap = grid::unit(0.75, 'mm'),
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = hm_meta$day_factor,
  column_gap = grid::unit(0, 'mm'),
  column_title = base::paste0('D', base::levels(hm_meta$day_factor)),
  column_title_gp = grid::gpar(fontsize = 11, fontfamily = plot_font_family),
  heatmap_legend_param = base::list(
    title_gp = grid::gpar(fontface = 'bold', fontfamily = plot_font_family),
    labels_gp = grid::gpar(fontfamily = plot_font_family)
  ),
  border = FALSE,
  use_raster = TRUE
)

cluster_genes <- base::names(gene_cluster)
cluster_vst <- vst[cluster_genes, meta$sample_id, drop = FALSE]

cluster_day_gene <- base::do.call(rbind, base::lapply(cluster_genes, function(gene) {
  means <- base::vapply(days, function(day) {
    samples <- meta$sample_id[meta$day_numeric == day]
    base::mean(cluster_vst[gene, samples], na.rm = TRUE)
  }, base::numeric(1))
  z <- base::as.numeric(base::scale(means))
  z[!base::is.finite(z)] <- 0
  base::data.frame(
    gene_id = gene,
    cluster = gene_cluster[[gene]],
    day_numeric = days,
    zscore = z,
    stringsAsFactors = FALSE
  )
}))

cluster_mean <- stats::aggregate(
  zscore ~ cluster + day_numeric,
  data = cluster_day_gene,
  FUN = base::mean
)

make_cluster_trajectory_plot <- function(cluster_name) {
  cluster_df <- cluster_day_gene[cluster_day_gene$cluster == cluster_name, ]
  mean_df <- cluster_mean[cluster_mean$cluster == cluster_name, ]
  cluster_n <- base::length(base::unique(cluster_df$gene_id))
  cluster_title <- base::paste0(base::sub('^C', 'Cluster ', cluster_name), ' (n=', cluster_n, ')')

  ggplot2::ggplot(cluster_df, ggplot2::aes(day_numeric, zscore, group = gene_id)) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = 'black',
      linetype = 'dashed',
      linewidth = 0.2
    ) +
    ggplot2::geom_line(color = 'gray72', linewidth = 0.18) +
    ggplot2::geom_line(
      data = mean_df,
      ggplot2::aes(day_numeric, zscore),
      inherit.aes = FALSE,
      color = cluster_colors[[cluster_name]],
      linewidth = 0.6
    ) +
    ggplot2::scale_x_continuous(breaks = days) +
    ggplot2::labs(title = cluster_title, x = 'Differentiation day', y = 'Mean VST z-score') +
    theme_pub(base_size = 8.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = cluster_colors[[cluster_name]],
        face = 'bold',
        family = plot_font_family
      )
    )
}

p_cluster_trajectories <- patchwork::wrap_plots(
  base::lapply(base::levels(gene_cluster), make_cluster_trajectory_plot),
  ncol = 2
) +
  patchwork::plot_annotation(tag_levels = 'a') &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = 'bold', family = plot_font_family, size = 11))

# 5.0 create sample correlation annotation dendrogram -----------------

temporal_vst <- vst[temporal_genes, meta$sample_id, drop = FALSE]
sample_correlation <- stats::cor(temporal_vst, method = 'pearson')
sample_distance <- stats::as.dist(1 - sample_correlation)
sample_hclust <- stats::hclust(sample_distance, method = 'average')
sample_dendrogram <- stats::as.dendrogram(sample_hclust)
annotation_carrier <- base::matrix(
  0,
  nrow = 1,
  ncol = base::ncol(sample_correlation),
  dimnames = base::list('hidden', base::colnames(sample_correlation))
)

day_color_function <- circlize::colorRamp2(
  base::seq(
    base::min(meta$day_numeric),
    base::max(meta$day_numeric),
    length.out = 256
  ),
  viridis::viridis(256, option = 'D')
)

column_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Day = meta$day_numeric,
  col = base::list(Day = day_color_function),
  annotation_name_gp = grid::gpar(
    fontfamily = plot_font_family,
    fontsize = 9,
    fontface = 'bold'
  ),
  annotation_legend_param = base::list(
    Day = base::list(
      at = base::c(1, 5, 10, 15),
      labels = base::c('1', '5', '10', '15'),
      legend_direction = 'horizontal',
      title_position = 'topcenter',
      title_gp = grid::gpar(fontfamily = plot_font_family, fontface = 'bold'),
      labels_gp = grid::gpar(fontfamily = plot_font_family)
    )
  ),
  simple_anno_size = grid::unit(11, 'mm')
)

p_correlation_annotation <- ComplexHeatmap::Heatmap(
  annotation_carrier,
  name = 'Hidden',
  top_annotation = column_annotation,
  cluster_columns = sample_dendrogram,
  cluster_rows = FALSE,
  column_split = 3,
  column_gap = grid::unit(1.5, 'mm'),
  show_column_dend = TRUE,
  column_dend_side = 'top',
  show_heatmap_legend = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  rect_gp = grid::gpar(col = NA, fill = NA),
  border = FALSE,
  width = grid::unit(230, 'mm'),
  height = grid::unit(2, 'mm')
)

# 6.0 fit PCA and validate loadings -----------------

base::message('CODEX_STEP pca: fitting reference PCA and maturation vector')

pca_input <- base::t(vst[temporal_genes, , drop = FALSE])
pca_fit <- stats::prcomp(pca_input, center = TRUE, scale. = FALSE)
pca_fit <- orient_pca_fit_by_day(pca_fit, sample_days)

var_pct <- base::round(base::summary(pca_fit)$importance[2, 1:3] * 100, 1)

pca_df <- base::data.frame(
  sample_id = base::rownames(pca_fit$x),
  PC1 = pca_fit$x[, 1],
  PC2 = pca_fit$x[, 2],
  PC3 = pca_fit$x[, 3],
  meta[base::match(base::rownames(pca_fit$x), meta$sample_id), ],
  row.names = NULL
)

p_pca_day <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = day_numeric)) +
  ggplot2::geom_point(size = 1.8, alpha = 1) +
  ggplot2::scale_color_viridis_c(name = "Day") +
  ggplot2::labs(
    x = base::paste0("PC1 (", var_pct[1], "%)"),
    y = base::paste0("PC2 (", var_pct[2], "%)")
  ) +
  theme_pub()

early_cluster_genes <- base::names(gene_cluster)[gene_cluster %in% base::c('C1', 'C2')]
late_cluster_genes <- base::names(gene_cluster)[gene_cluster %in% base::c('C3', 'C4')]

gene_sets <- base::list(
  All = temporal_genes,
  `C1+C2 (Early)` = early_cluster_genes,
  `C3+C4 (Late)` = late_cluster_genes
)

base::stopifnot(base::all(base::vapply(gene_sets, base::length, base::integer(1)) > 1L))
base::stopifnot(base::all(base::unlist(gene_sets, use.names = FALSE) %in% base::rownames(vst)))

calculate_pca <- function(gene_ids, set_label) {
  set_pca_input <- base::t(vst[gene_ids, meta$sample_id, drop = FALSE])
  set_pca_fit <- stats::prcomp(set_pca_input, center = TRUE, scale. = FALSE)
  set_pca_fit <- orient_pca_fit_by_day(set_pca_fit, sample_days)
  set_pca_var <- base::round(base::summary(set_pca_fit)$importance[2, 1:3] * 100, 2)

  set_pca_df <- base::data.frame(
    sample_id = base::rownames(set_pca_fit$x),
    PC1 = set_pca_fit$x[, 1],
    PC2 = set_pca_fit$x[, 2],
    PC3 = set_pca_fit$x[, 3],
    day_numeric = meta$day_numeric[
      base::match(base::rownames(set_pca_fit$x), meta$sample_id)
    ],
    cell_line = meta$cell_line[
      base::match(base::rownames(set_pca_fit$x), meta$sample_id)
    ],
    gene_set = set_label,
    stringsAsFactors = FALSE
  )

  base::list(
    data = set_pca_df,
    variance_percent = set_pca_var,
    gene_count = base::length(gene_ids)
  )
}

pca_results <- base::Map(calculate_pca, gene_sets, base::names(gene_sets))

make_pca_plot <- function(pca_result, plot_title) {
  plot_pca_df <- pca_result$data
  plot_pca_var <- pca_result$variance_percent
  color_option <- get_pca_color_option(plot_title)

  ggplot2::ggplot(plot_pca_df, ggplot2::aes(PC1, PC2, color = day_numeric)) +
    ggplot2::geom_point(size = 1.35, alpha = 0.9) +
    ggplot2::scale_color_viridis_c(name = 'Day', option = color_option) +
    ggplot2::labs(
      title = plot_title,
      x = base::paste0('PC1 (', plot_pca_var[[1]], '%)'),
      y = base::paste0('PC2 (', plot_pca_var[[2]], '%)')
    ) +
    theme_pub(base_size = 9) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5),
      legend.position = 'none'
    )
}

make_pc1_time_plot <- function(pca_result) {
  plot_pca_df <- pca_result$data
  color_option <- get_pca_color_option(plot_pca_df$gene_set[[1]])
  label_text <- format_lm_label(plot_pca_df$day_numeric, plot_pca_df$PC1)
  label_x <- base::min(plot_pca_df$day_numeric)
  label_y <- base::max(plot_pca_df$PC1, na.rm = TRUE)

  ggplot2::ggplot(plot_pca_df, ggplot2::aes(day_numeric, PC1)) +
    ggplot2::geom_point(ggplot2::aes(color = day_numeric), size = 1.2, alpha = 0.9) +
    ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = FALSE,
      color = 'black',
      linewidth = 0.45
    ) +
    ggplot2::annotate(
      'text',
      x = label_x,
      y = label_y,
      label = label_text,
      hjust = 0,
      vjust = 1,
      family = plot_font_family,
      size = 4
    ) +
    ggplot2::scale_x_continuous(breaks = days) +
    ggplot2::scale_color_viridis_c(name = 'Day', option = color_option) +
    ggplot2::labs(
      x = 'Differentiation day',
      y = 'PC1 position'
    ) +
    theme_pub(base_size = 10.5) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 14),
      legend.position = 'none'
    )
}

p_all_pca <- make_pca_plot(pca_results[['All']], 'All')
p_early_pca <- make_pca_plot(pca_results[['C1+C2 (Early)']], 'C1+C2 (Early)')
p_late_pca <- make_pca_plot(pca_results[['C3+C4 (Late)']], 'C3+C4 (Late)')

p_all_pc1 <- make_pc1_time_plot(pca_results[['All']])
p_early_pc1 <- make_pc1_time_plot(pca_results[['C1+C2 (Early)']])
p_late_pc1 <- make_pc1_time_plot(pca_results[['C3+C4 (Late)']])

p_pca_day_with_fit <- (p_all_pca | p_all_pc1) +
  patchwork::plot_layout(widths = base::c(1.05, 0.95))
p_pca_day <- p_pca_day_with_fit

p_pca_grid <- (
  p_early_pca | p_late_pca
) / (
  p_early_pc1 | p_late_pc1
)

make_line_trace <- function(line_df, trace_name) {
  base::list(
    x = line_df$PC1,
    y = line_df$PC2,
    z = line_df$PC3,
    name = trace_name
  )
}

calculate_best_fit_line <- function(fit_pca_df) {
  day_range <- base::range(fit_pca_df$day_numeric)
  line_days <- base::seq(day_range[[1]], day_range[[2]], length.out = 100)
  prediction_df <- base::data.frame(day_numeric = line_days)

  base::data.frame(
    day_numeric = line_days,
    PC1 = stats::predict(stats::lm(PC1 ~ day_numeric, data = fit_pca_df), prediction_df),
    PC2 = stats::predict(stats::lm(PC2 ~ day_numeric, data = fit_pca_df), prediction_df),
    PC3 = stats::predict(stats::lm(PC3 ~ day_numeric, data = fit_pca_df), prediction_df),
    stringsAsFactors = FALSE
  )
}

add_3d_arrow_cone <- function(plot_obj, line_df) {
  start_point <- line_df[1, base::c('PC1', 'PC2', 'PC3')]
  end_point <- line_df[base::nrow(line_df), base::c('PC1', 'PC2', 'PC3')]
  direction <- base::as.numeric(end_point - start_point)
  direction <- direction / base::sqrt(base::sum(direction^2))

  plotly::add_trace(
    plot_obj,
    type = 'cone',
    x = end_point$PC1,
    y = end_point$PC2,
    z = end_point$PC3,
    u = direction[[1]],
    v = direction[[2]],
    w = direction[[3]],
    anchor = 'tip',
    sizemode = 'absolute',
    sizeref = 8,
    colorscale = base::list(base::list(0, 'black'), base::list(1, 'black')),
    showscale = FALSE,
    name = 'D1 to D15 arrow',
    inherit = FALSE
  )
}

make_3d_pca_plot <- function(pca_result, line_df, plot_title, line_name) {
  plot_pca_df <- pca_result$data
  plot_pca_var <- pca_result$variance_percent
  line_trace <- make_line_trace(line_df, line_name)

  hover_text <- base::paste0(
    'Sample: ', plot_pca_df$sample_id,
    '<br>Line: ', plot_pca_df$cell_line,
    '<br>Day: ', plot_pca_df$day_numeric
  )

  plot_obj <- plotly::plot_ly(
    plot_pca_df,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~day_numeric,
    colors = viridis::viridis(100),
    type = 'scatter3d',
    mode = 'markers',
    marker = base::list(size = 4, opacity = 0.88),
    text = hover_text,
    hoverinfo = 'text',
    name = 'Samples',
    colorbar = base::list(title = 'Day')
  )

  plot_obj <- plotly::add_trace(
    plot_obj,
    x = line_trace$x,
    y = line_trace$y,
    z = line_trace$z,
    type = 'scatter3d',
    mode = 'lines',
    line = base::list(color = 'black', width = 7),
    name = line_trace$name,
    inherit = FALSE
  )
  plot_obj <- plotly::add_trace(
    plot_obj,
    x = line_trace$x[base::c(1, base::length(line_trace$x))],
    y = line_trace$y[base::c(1, base::length(line_trace$y))],
    z = line_trace$z[base::c(1, base::length(line_trace$z))],
    type = 'scatter3d',
    mode = 'markers+text',
    marker = base::list(color = 'black', size = base::c(4, 7)),
    text = base::c('  D1  ', '  D15  '),
    textposition = base::c('top left', 'top right'),
    name = 'Trajectory endpoints',
    inherit = FALSE
  )
  plot_obj <- add_3d_arrow_cone(plot_obj, line_df)

  plotly::layout(
    plot_obj,
    title = plot_title,
    legend = base::list(x = 0.02, y = 0.98),
    font = base::list(family = plot_font_family),
    scene = base::list(
      xaxis = base::list(title = base::paste0('PC1 (', plot_pca_var[[1]], '%)')),
      yaxis = base::list(title = base::paste0('PC2 (', plot_pca_var[[2]], '%)')),
      zaxis = base::list(title = base::paste0('PC3 (', plot_pca_var[[3]], '%)'))
    )
  )
}

all_pca_df <- pca_results[['All']]$data
best_fit_line <- calculate_best_fit_line(all_pca_df)

p_3d_best_fit <- make_3d_pca_plot(
  pca_result = pca_results[['All']],
  line_df = best_fit_line,
  plot_title = 'All Temporal Genes: 3D PCA Best-Fit Trajectory',
  line_name = 'Best-fit D1 to D15 trajectory'
)
p_pca_3d <- p_3d_best_fit

pc1_loadings <- pca_fit$rotation[, 'PC1']

top_pc1_pos <- base::names(
  base::sort(pc1_loadings, decreasing = TRUE)
)[base::seq_len(n_pc1_loading_genes_per_direction)]
top_pc1_neg <- base::names(
  base::sort(pc1_loadings, decreasing = FALSE)
)[base::seq_len(n_pc1_loading_genes_per_direction)]
top_pc1_genes <- base::c(top_pc1_neg, top_pc1_pos)

pc1_loading_df <- base::data.frame(
  gene_id = top_pc1_genes,
  gene_label = label_ensembl_genes(top_pc1_genes),
  loading = pc1_loadings[top_pc1_genes],
  direction = base::rep(
    base::c('PC1-', 'PC1+'),
    each = n_pc1_loading_genes_per_direction
  ),
  stringsAsFactors = FALSE
)

pc1_loading_df$gene_label <- base::factor(
  pc1_loading_df$gene_label,
  levels = pc1_loading_df$gene_label[base::order(pc1_loading_df$loading)]
)

p_pc1_loadings <- ggplot2::ggplot(pc1_loading_df, ggplot2::aes(gene_label, loading, fill = direction)) +
  ggplot2::geom_col(width = 0.72, color = 'black', linewidth = 0.18) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = base::c('PC1-' = '#2563EB', 'PC1+' = '#DC2626')) +
  ggplot2::labs(x = NULL, y = 'PC1 loading', fill = NULL) +
  theme_pub(base_size = 8) +
  ggplot2::theme(
    legend.position = 'top',
    axis.text.y = ggplot2::element_text(face = 'italic', family = plot_font_family)
  )

top_pc1_z <- base::t(base::scale(base::t(vst[top_pc1_genes, meta$sample_id, drop = FALSE])))
top_pc1_z[!base::is.finite(top_pc1_z)] <- 0

top_pc1_long <- base::data.frame(
  gene_id = base::rep(base::rownames(top_pc1_z), times = base::ncol(top_pc1_z)),
  sample_id = base::rep(base::colnames(top_pc1_z), each = base::nrow(top_pc1_z)),
  zscore = base::as.vector(top_pc1_z),
  stringsAsFactors = FALSE
)
top_pc1_long$day_numeric <- meta$day_numeric[
  base::match(top_pc1_long$sample_id, meta$sample_id)
]
top_pc1_long$direction <- pc1_loading_df$direction[
  base::match(top_pc1_long$gene_id, pc1_loading_df$gene_id)
]

top_pc1_mean <- stats::aggregate(
  zscore ~ direction + gene_id + day_numeric,
  data = top_pc1_long,
  FUN = base::mean
)

p_pc1_trajectories <- ggplot2::ggplot(top_pc1_mean, ggplot2::aes(day_numeric, zscore, group = gene_id)) +
  ggplot2::geom_line(ggplot2::aes(color = direction), linewidth = 0.4) +
  ggplot2::stat_summary(
    ggplot2::aes(group = direction, color = direction),
    fun = base::mean,
    geom = 'line',
    linewidth = 0.65
  ) +
  ggplot2::facet_wrap(~direction, nrow = 1) +
  ggplot2::scale_x_continuous(breaks = base::seq(base::min(days), base::max(days), 2)) +
  ggplot2::scale_color_manual(values = base::c('PC1-' = '#2563EB', 'PC1+' = '#DC2626')) +
  ggplot2::guides(color = 'none') +
  ggplot2::labs(x = 'Differentiation day', y = 'Mean z-scored VST') +
  theme_pub(base_size = 8)

p_pc1_validation <- p_pc1_loadings / p_pc1_trajectories +
  patchwork::plot_layout(heights = base::c(1.35, 1)) +
  patchwork::plot_annotation(
    tag_levels = 'a',
    theme = ggplot2::theme(
      text = ggplot2::element_text(family = plot_font_family, color = 'black'),
      plot.tag = ggplot2::element_text(face = 'bold', family = plot_font_family, size = 11)
    )
  )

# 7.0 score samples on the best-fit maturation vector -----------------

pc_cols <- base::c('PC1', 'PC2', 'PC3')

start_day <- base::min(days)
end_day <- base::max(days)
start_pt <- base::as.numeric(best_fit_line[1, pc_cols])
end_pt <- base::as.numeric(best_fit_line[base::nrow(best_fit_line), pc_cols])
base::names(start_pt) <- pc_cols
base::names(end_pt) <- pc_cols

vector_raw <- end_pt - start_pt
vector_length_sq <- base::sum(vector_raw^2)
base::stopifnot(vector_length_sq > 0)

score_raw <- base::as.matrix(pca_df[, pc_cols]) -
  base::matrix(start_pt, nrow = base::nrow(pca_df), ncol = base::length(pc_cols), byrow = TRUE)

pca_df$maturation_score <- base::as.numeric(score_raw %*% vector_raw) / vector_length_sq

projection_coords <- base::matrix(
  start_pt,
  nrow = base::nrow(pca_df),
  ncol = base::length(pc_cols),
  byrow = TRUE
) + pca_df$maturation_score %o% vector_raw

pca_df$vector_PC1 <- projection_coords[, 1]
pca_df$vector_PC2 <- projection_coords[, 2]
pca_df$vector_PC3 <- projection_coords[, 3]

p_vector <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2)) +
  ggplot2::geom_point(ggplot2::aes(color = day_numeric), size = 1.2, alpha = 1) +
  ggplot2::geom_path(
    data = best_fit_line,
    ggplot2::aes(PC1, PC2),
    inherit.aes = FALSE,
    color = 'black',
    linewidth = 0.45
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = start_pt[["PC1"]], y = start_pt[["PC2"]],
      xend = end_pt[["PC1"]], yend = end_pt[["PC2"]]
    ),
    inherit.aes = FALSE,
    linewidth = 0.9,
    color = "#059669",
    arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed")
  ) +
  ggplot2::annotate("text",
    x = start_pt[['PC1']], y = start_pt[['PC2']],
    label = base::paste0('Day ', start_day), vjust = -1, size = 3,
    family = plot_font_family
  ) +
  ggplot2::annotate("text",
    x = end_pt[['PC1']], y = end_pt[['PC2']],
    label = base::paste0('Day ', end_day), vjust = -1, size = 3,
    family = plot_font_family
  ) +
  ggplot2::scale_color_viridis_c(name = "Day") +
  ggplot2::labs(
    x = base::paste0("PC1 (", var_pct[1], "%)"),
    y = base::paste0("PC2 (", var_pct[2], "%)")
  ) +
  theme_pub()

score_summary <- stats::aggregate(
  maturation_score ~ day_numeric,
  data = pca_df,
  FUN = function(x) base::c(mean = base::mean(x), sd = stats::sd(x), n = base::length(x))
)

score_summary <- base::data.frame(
  day_numeric = score_summary$day_numeric,
  mean_score = score_summary$maturation_score[, "mean"],
  sd_score = score_summary$maturation_score[, "sd"],
  n = score_summary$maturation_score[, "n"],
  row.names = NULL
)

p_score_by_day <- ggplot2::ggplot(pca_df, ggplot2::aes(day_numeric, maturation_score)) +
  ggplot2::geom_hline(
    yintercept = base::c(0, 1), linetype = "dashed",
    linewidth = 0.3, color = "black"
  ) +
  ggplot2::geom_point(ggplot2::aes(color = base::factor(cell_line)), size = 1.3, alpha = 1) +
  ggplot2::stat_summary(
    fun = base::mean, geom = "line", linewidth = 0.8,
    color = "black"
  ) +
  ggplot2::scale_x_continuous(breaks = base::seq(base::min(days), base::max(days), 2)) +
  ggplot2::guides(color = "none") +
  ggplot2::labs(x = "Differentiation day", y = "Maturation score") +
  theme_pub()

make_panel_df <- function(x_pc, y_pc, panel_label) {
  base::data.frame(
    x = pca_df[[x_pc]],
    y = pca_df[[y_pc]],
    day_numeric = pca_df$day_numeric,
    panel = panel_label
  )
}

panel_df <- base::rbind(
  make_panel_df("PC1", "PC2", "PC1 vs PC2"),
  make_panel_df("PC1", "PC3", "PC1 vs PC3"),
  make_panel_df("PC2", "PC3", "PC2 vs PC3")
)

seg_df <- base::data.frame(
  x = base::c(start_pt[["PC1"]], start_pt[["PC1"]], start_pt[["PC2"]]),
  y = base::c(start_pt[["PC2"]], start_pt[["PC3"]], start_pt[["PC3"]]),
  xend = base::c(end_pt[["PC1"]], end_pt[["PC1"]], end_pt[["PC2"]]),
  yend = base::c(end_pt[["PC2"]], end_pt[["PC3"]], end_pt[["PC3"]]),
  panel = base::c("PC1 vs PC2", "PC1 vs PC3", "PC2 vs PC3")
)

p_projection_pairwise <- ggplot2::ggplot(panel_df, ggplot2::aes(x, y)) +
  ggplot2::geom_segment(
    data = seg_df,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    color = "#059669",
    linewidth = 0.8,
    arrow = grid::arrow(length = grid::unit(0.13, "cm"), type = "closed")
  ) +
  ggplot2::geom_point(ggplot2::aes(color = day_numeric), size = 1.1, alpha = 1) +
  ggplot2::facet_wrap(~panel, nrow = 1, scales = "free") +
  ggplot2::scale_color_viridis_c(name = "Day") +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_pub() +
  ggplot2::theme(legend.position = "top")

base::message('CODEX_DONE built tutorial figures and summary tables')
