# Differentiation timing score

These R functions estimate differentiation timing from bulk RNA-seq expression by comparing samples with a reference time course. The workflow selects time-dependent genes, fits a trajectory through the reference timepoints in principal component analysis (PCA) space, and projects each sample onto that trajectory. It returns a predicted reference time, a differentiation score from 0 to 1, and a squared distance that helps assess how closely the sample fits the reference.

The [online tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/) explains the method using a cardiomyocyte differentiation time course. It includes gene selection, trajectory fitting, figures, leave-one-cell-line-out validation, and guidance for applying the method to another study.

## Quick start

Download or clone this repository and open R in its root directory. Gene selection requires DESeq2 and edgeR from Bioconductor:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR"))
```

The example below uses the included counts, variance-stabilized (VST) expression, and metadata from [GSE122380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380). The loader aligns both matrices to the metadata sample order. Gene selection accounts for cell-line differences in this reference.

```r
source("data/load_GSE122380_data.R")
source("functions/get_temporal_genes.R")
source("functions/score_differentiation_timing.R")

reference_data <- load_GSE122380_data()
temporal_selection <- get_temporal_genes(
  raw_counts = reference_data$counts,
  vst_expression = reference_data$vst,
  metadata = reference_data$metadata,
  adjustment_covariates = "cell_line"
)

timing_fit <- score_differentiation_timing(
  expression_matrix = reference_data$vst,
  metadata = reference_data$metadata,
  temporal_genes = temporal_selection$temporal_genes
)

timing_fit$scores[, c(
  "sample_id",
  "observed_time",
  "predicted_time",
  "differentiation_score",
  "squared_distance"
)]
```

This example scores the reference samples themselves. To evaluate experimental samples, select genes from an independent reference cohort and pass its IDs through `reference_samples` when scoring. The tutorial shows how to select a reference group and assess predictions in held-out samples.

## How scoring works

[`get_temporal_genes()`](functions/get_temporal_genes.R) selects genes with sufficient expression, a significant time effect in a DESeq2 likelihood-ratio test, and sufficient variation across mean VST expression at reference timepoints. All three filters use the reference samples:

<p align="center"><img src="docs/assets/diagrams/temporal_gene_selection_v2.png" alt="Reference counts and metadata pass through an expression filter, DESeq2 likelihood-ratio test, and VST range filter to select temporal genes." width="560"></p>

[`score_differentiation_timing()`](functions/score_differentiation_timing.R) fits centered, unscaled PCA to the mean temporal-gene expression profile at each reference timepoint. It connects the resulting centroids in time order and projects each sample onto the nearest trajectory segment:

<p align="center"><img src="docs/assets/diagrams/reference_trajectory_v2.png" alt="Temporal genes and reference samples define a PCA trajectory. Sample projections give predicted time, differentiation score, and squared distance." width="700"></p>

[`run_loo_maturation()`](functions/run_loo_maturation.R) validates prediction by repeating gene selection and trajectory fitting while holding out each cell line in turn.

## Input requirements

- **Gene selection:** integer-like, unnormalized counts and VST expression with matching gene and sample identifiers. Both matrices have genes in rows and samples in columns. The function calculates TMM-normalized counts per million (CPM) internally.
- **Metadata:** one row per sample, matching sample IDs, and finite numeric timepoints for the reference samples. The default column names are `sample_id` and `day_numeric`.
- **Scoring:** normalized expression, matching metadata, and selected temporal genes. The reference requires at least two distinct, finite numeric timepoints; non-reference samples can have missing timepoints.
- **Validation:** aligned count and VST matrices, with `sample_id`, `cell_line`, and `day_numeric` in the metadata.

## Interpret the result

`timing_fit$scores` contains one row per sample:

- `predicted_time` is the estimated position on the original reference-time scale.
- `differentiation_score` maps the earliest and latest reference timepoints to 0 and 1.
- `squared_distance` is the squared distance to the nearest fitted trajectory segment in PCA space.

Interpret the score relative to the reference timepoints and compare samples scored with the same reference model. Read `squared_distance` alongside the score to identify samples whose expression fits the reference path poorly. See the [tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/) for interpretation, validation results, and references; its source is [`docs/tutorial/tutorial.Rmd`](docs/tutorial/tutorial.Rmd).
