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
```

## Usage
Fitting the ZIFA Model
``` r
# Simulate some data
set.seed(123)
true_Z <- matrix(rnorm(100), ncol = 2)
lambda <- matrix(runif(1000, -1, 1), nrow = 500)
mu <- true_Z %*% t(lambda)
Y <- mu + matrix(rnorm(500*50), nrow = 500)
# Introduce zero-inflation
mask <- matrix(rbinom(500*50, 1, prob = 0.3), nrow = 500)
Y[mask == 1] <- 0

# Fit the ZIFA model
result <- fit_zifa(Y, k = 2)
```

Block-wise ZIFA
```r
# Fit the block-wise ZIFA model
result_block <- fit_block_zifa(Y, k = 2)
```

Using SummarizedExperiment Objects
```r
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  library(SummarizedExperiment)
  # Create a SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = Y))
  # Fit the ZIFA model directly on the SummarizedExperiment
  result_se <- fit_zifa(se, k = 2)
}
```

## Dependencies
- R (>= 3.5)
- ggplot2 (for plotting, used in plot_zifa)
- SummarizedExperiment (optional, for enhanced compatibility)

## License
This project is licensed under the MIT License.

## Acknowledgements
This package is inspired by the original ZIFA algorithm and aims to provide a user-friendly implementation for the R community, as originally implemented in:
Pierson, E., & Yau, C. (2015). ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. Genome biology, 16(1), 1-10.
