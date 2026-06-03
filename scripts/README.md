# Cardiomyocyte Maturation Score Tutorial

The active tutorial is maintained as one source file:

```text
scripts/cardiomyocyte_maturation_score.Rmd
```

The computational workflow is kept in:

```text
scripts/cardiomyocyte_maturation_score_analysis.R
```

The Rmd sources that script, then presents the method, input data,
figures, interpretation, reproducibility commands, and references in a
single maintainable document.

## Render

```bash
tools/r_codex_utils render scripts \
  --expect-output docs/index.html \
  --show-stdout \
  --stdout-tail 120
```

The rendered GitHub Pages site is written to `docs/index.html`.
