#' Preprocess Data for ZIFA
#'
#' Preprocesses raw count data for use with ZIFA by log-transforming it.
#' This is the recommended preprocessing step before running ZIFA.
#'
#' @param counts Matrix of raw count data (genes x cells)
#' @param base Base for logarithm transformation (default: 2)
#' @param pseudocount Pseudocount to add before log transformation (default: 1)
#'
#' @return Log-transformed data matrix of same dimensions as input
#'
#' @details
#' ZIFA expects log-transformed count data, not raw counts. This function

#' applies the transformation: log_base(counts + pseudocount).
#'
#' Common choices:
#' - base=2, pseudocount=1: log2(counts + 1), commonly used in scRNA-seq
#' - base=exp(1), pseudocount=1: natural log, ln(counts + 1)
#' - base=10, pseudocount=1: log10(counts + 1)
#'
#' @examples
#' # Generate raw count data
#' counts <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
#' 
#' # Preprocess data
#' Y <- preprocess_data(counts)
#' 
#' # Now Y is ready for ZIFA
#' # result <- fit_zifa(Y, k = 2)
#'
#' @export
preprocess_data <- function(counts, base = 2, pseudocount = 1) {
  if (any(counts < 0)) {
    stop("Count data cannot contain negative values")
  }
  log(counts + pseudocount, base = base)
}

#' Plot ZIFA Results
#'
#' Creates a scatter plot of ZIFA results in two dimensions.
#'
#' @param Z Matrix of ZIFA projections (cells x dimensions)
#' @param labels Optional vector of labels for coloring points (e.g., cell types)
#' @param dims Which dimensions to plot (default: first two dimensions)
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' # After running ZIFA
#' result <- fit_zifa(Y, k = 2)
#' 
#' # Basic plot
#' plot_zifa(result$Z)
#' 
#' # With cell type labels
#' cell_types <- c(rep("TypeA", 25), rep("TypeB", 25))
#' plot_zifa(result$Z, labels = cell_types)
#' 
#' # Plot dimensions 2 and 3
#' result_3d <- fit_zifa(Y, k = 3)
#' plot_zifa(result_3d$Z, dims = c(2, 3))
#' }
#'
#' @export
plot_zifa <- function(Z, labels = NULL, dims = c(1, 2)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install it with: install.packages('ggplot2')")
  }
  
  if (max(dims) > ncol(Z)) {
    stop(sprintf("Requested dimensions %s but Z only has %d dimensions", 
                 paste(dims, collapse = ", "), ncol(Z)))
  }
  
  # Extract dimensions to plot
  plot_data <- as.data.frame(Z[, dims])
  colnames(plot_data) <- c("Dim1", "Dim2")
  
  # Add labels if provided
  if (!is.null(labels)) {
    if (length(labels) != nrow(Z)) {
      stop("Length of labels must match number of rows in Z")
    }
    plot_data$Label <- as.factor(labels)
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Dim1, y = .data$Dim2, color = .data$Label)) +
      ggplot2::geom_point(alpha = 0.7, size = 2)
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Dim1, y = .data$Dim2)) +
      ggplot2::geom_point(alpha = 0.7, size = 2, color = "#3366CC")
  }
  
  p + ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = paste("ZIFA Dimension", dims[1]),
      y = paste("ZIFA Dimension", dims[2]),
      title = "ZIFA Projection",
      color = "Group"
    )
}
