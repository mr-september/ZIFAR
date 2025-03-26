# ZIFAR

**ZIFAR** is a native R package implementing Zero-Inflated Factor Analysis (ZIFA) for dimensionality reduction, particularly tailored for zero-inflated data such as single-cell gene expression datasets.

## Features

- **Dimensionality Reduction:** Performs factor analysis adapted for datasets with a high proportion of zeros.
- **Block Processing:** Efficiently processes large datasets by dividing genes into blocks.
- **Custom Initialisation:** Utilises standard factor analysis for robust parameter initialisation.
- **SummarizedExperiment Compatibility:** Direct support for `SummarizedExperiment` objects from Bioconductor.

## Installation

You can install the development version from GitHub:

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
  
devtools::install_github("mr-september/ZIFAR")
install.packages("devtools")
devtools::install_github("mr-september/ZIFAR")
```

## Overview

Single-cell RNA sequencing data often contains many zero counts due to both biological zeros and technical dropouts. ZIFA models the dropout probability as a function of the expected expression level, providing more accurate dimensionality reduction than methods that don't account for this zero-inflation.

This package provides:
- A complete R implementation of the ZIFA algorithm
- Support for sparse matrices to efficiently handle large datasets
- Native compatibility with Bioconductor's SummarizedExperiment objects
- Visualization functions for exploring the reduced dimensions

## Usage

Basic usage with a count matrix:

```r
library(ZIFAR)

With a regular count matrix
results <- zifa(counts_matrix, k = 10)

Visualize the results
plot_zifa(results)

With a SummarizedExperiment object
library(SummarizedExperiment)
se_results <- zifa(summarized_experiment, k = 10)

```

## Parameters

- `data`: Count matrix (genes × cells) or SummarizedExperiment object
- `k`: Number of latent dimensions
- `n_iterations`: Maximum number of iterations for optimization
- `convergence_threshold`: Convergence criterion
- `dropout_model`: Model to use for dropout events (default: "exponential")

## References

This package implements the method described in:

Pierson, E., & Yau, C. (2015). ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. Genome biology, 16(1), 1-10.

## License

MIT
