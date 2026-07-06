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
# - docs/index.html: GitHub Pages copy of the rendered tutorial site
#
# Purpose:
# Render the active GSE122380 bookdown tutorial in report/ and copy the static
# site to docs/ for GitHub Pages.

# 1.0 validate render dependencies -----------------

base::stopifnot(base::requireNamespace('bookdown', quietly = TRUE))

# 2.0 render bookdown site -----------------

bookdown::render_book(
  input = 'report',
  output_format = 'bookdown::gitbook',
  clean = TRUE)

base::stopifnot(base::file.exists('report/index.html'))

# 3.0 refresh GitHub Pages copy -----------------

if (base::dir.exists('docs')) {
  base::unlink('docs', recursive = TRUE)
}
base::dir.create('docs', recursive = TRUE, showWarnings = FALSE)
base::file.copy(
  from = base::list.files('report', all.files = TRUE, no.. = TRUE, full.names = TRUE),
  to = 'docs',
  recursive = TRUE,
  copy.date = TRUE
)
base::file.create(base::file.path('docs', '.nojekyll'))
base::stopifnot(base::file.exists('docs/index.html'))
