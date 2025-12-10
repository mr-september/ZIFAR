# ZIFAR

**ZIFAR** is a native R implementation of Zero-Inflated Factor Analysis (ZIFA) for dimensionality reduction of single-cell RNA-seq data.

[![R-CMD-check](https://github.com/mr-september/ZIFAR/workflows/R-CMD-check/badge.svg)](https://github.com/mr-september/ZIFAR/actions)

## Overview

Single-cell RNA sequencing data often contains excess zeros due to "dropout" events where expressed genes fail to be detected. Standard dimensionality reduction methods (PCA, factor analysis) don't account for this, leading to suboptimal results.

ZIFA explicitly models the dropout process as:
- **Probability of dropout** ∝ exp(-λ × expression²)
- Lower expression → higher dropout probability

This R package is a faithful port of the [original Python implementation](https://github.com/epierson9/ZIFA) by Pierson & Yau (2015).

## Installation

```r
# Install from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
  
devtools::install_github("mr-september/ZIFAR")
```

## Quick Start

```r
library(ZIFAR)

# Simulate single-cell data (genes x cells)
set.seed(123)
counts <- matrix(rpois(1000 * 100, lambda = 5), nrow = 1000, ncol = 100)

# Preprocess: log-transform
Y <- preprocess_data(counts)

# Add realistic dropout
for (i in 1:nrow(Y)) {
  for (j in 1:ncol(Y)) {
    if (runif(1) < exp(-0.5 * Y[i, j]^2)) Y[i, j] <- 0
  }
}

# Fit ZIFA (genes x cells input)
result <- fit_zifa(Y, k = 2)

# View results
head(result$Z)  # Low-dimensional cell embeddings (cells x k)
```

## For Large Datasets

When working with >2000 genes, use `fit_block_zifa()` for efficiency:

```r
result <- fit_block_zifa(Y, k = 2, n_blocks = 4)
```

## Bioconductor Integration

ZIFAR works directly with `SummarizedExperiment` objects:

```r
library(SummarizedExperiment)
se <- SummarizedExperiment(assays = list(logcounts = Y))
result <- fit_zifa(se, k = 2)
```

## Visualization

```r
# Basic plot
plot_zifa(result$Z)

# With cell type labels
plot_zifa(result$Z, labels = cell_types)
```

## Output

`fit_zifa()` returns a list containing:
- `Z`: Cell embeddings matrix (cells × k dimensions)
- `model_params`: List with A (loadings), mus (means), sigmas (noise), decay_coef (λ)
- `n_iter`: Number of EM iterations

## Citation

If you use ZIFAR, please cite the original ZIFA paper:

> Pierson, E., & Yau, C. (2015). ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. *Genome Biology*, 16(1), 241. https://doi.org/10.1186/s13059-015-0805-z

## License

MIT License
