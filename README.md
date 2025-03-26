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
