# Maturation Scoring Tutorial

PCA-based maturation scoring for bulk RNA-seq time-course
experiments.

This repository contains a reusable, adaptable workflow using the
QC-processed GSE122380 iPSC-to-cardiomyocyte differentiation dataset.

## Active Data

The tutorial uses only the small set of QC-processed cardiomyocyte
inputs in `data/`.

- `GSE122380_metadata.rds`: QC-filtered sample metadata
- `GSE122380_counts.rds`: filtered raw counts for retained genes
  and samples
- `GSE122380_vst.rds`: VST matrix after cell-line effect removal
  while protecting nonlinear day effects

Raw GSE122380 files are intentionally not kept in this repository.

## Repository Structure

```text
data/                      Active QC-processed cardiomyocyte inputs
scripts/                   Rmd source, tutorial analysis, and reusable R functions
analysis/legacy_d21_t21/   Archived D21/T21 tutorial and helper scripts
docs/                      Rendered tutorial site for GitHub Pages
```

## Rendering

```bash
tools/r_codex_utils render scripts \
  --expect-output docs/index.html \
  --show-stdout \
  --stdout-tail 120
```

Rendered site output is written to `docs/`.

## Contact

**Author:** Zoheb Khan

**Affiliation:** Bioinformatician @ Moskowitz Lab at the University of Chicago Department of Pathology, Pediatrics, and Human Genetics

**Email:** zohebkhan600@gmail.com

**Website:** https://zohebkhan1.github.io/
