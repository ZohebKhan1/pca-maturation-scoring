# Created:
# 2026-06-23
#
# Inputs:
# - report/_bookdown.yml: bookdown configuration
# - report/tutorial.Rmd: tutorial source document
# - scripts/01_build_tutorial_objects.R: sourced analysis
# - data/GSE122380_metadata.rds: QC-filtered GSE122380 sample metadata
# - data/GSE122380_counts.rds: filtered raw count matrix
# - data/GSE122380_vst.rds: VST expression matrix aligned to metadata
#
# Outputs:
# - report/index.html: rendered tutorial site entry point
#
# Purpose:
# Render the active GSE122380 bookdown tutorial in report/.

# 1.0 validate render dependencies -----------------

base::stopifnot(base::requireNamespace('bookdown', quietly = TRUE))

# 2.0 render bookdown site -----------------

bookdown::render_book(
  input = 'report',
  output_format = 'bookdown::gitbook',
  clean = TRUE)

base::stopifnot(base::file.exists('report/index.html'))
