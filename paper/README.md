# GLBFP paper kit

This directory contains reproducible material for software-paper submission
tracks (The R Journal or JSS).

## Contents

- `rjournal_outline.qmd`: short-form paper structure for The R Journal.
- `jss_outline.qmd`: long-form paper structure for Journal of Statistical Software.
- `reproduce.R`: runs the benchmark pipeline and writes compact result tables.
- `results/`: generated CSV outputs used in manuscript tables/figures.

## Reproduce

From repository root:

```bash
Rscript paper/reproduce.R
quarto render paper/rjournal_outline.qmd
quarto render paper/jss_outline.qmd
```

The script only uses local package code and simulated data; no network access is required.
