#' Preprocess Data for ZIFA
#'
#' Preprocesses raw count data for use with ZIFA by log-transforming it.
#'
#' @param counts Matrix of raw count data
#' @param base Base for logarithm transformation (default: 2)
#' @param pseudocount Pseudocount to add before log transformation (default: 1)
#'
#' @return Log-transformed data matrix
#'
#' @examples
#' # Generate raw count data
#' counts <- matrix(rpois(1000, lambda = 5), nrow = 100)
#' # Preprocess data
#' Y <- preprocess_data(counts)
#'
#' @export
preprocess_data <- function(counts, base = 2, pseudocount = 1) {
  log(counts + pseudocount, base = base)
}

#' Plot ZIFA Results
#'
#' Creates a scatter plot of ZIFA results in two dimensions.
#'
#' @param Z Matrix of ZIFA projections
#' @param labels Optional vector of labels for coloring points
#' @param dims Dimensions to plot (default: first two dimensions)
#'
#' @return A ggplot2 object
#'
#' @examples
#' # Generate simulated data
#' set.seed(123)
#' Y <- matrix(rpois(1000, lambda = 5), nrow = 100)
#' Y[sample(length(Y), 500)] <- 0
#' result <- fit_zifa(Y, k = 2)
#' # Plot results
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_zifa(result$Z)
#' }
#'
#' @export
plot_zifa <- function(Z, labels = NULL, dims = c(1, 2)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install it.")
  }
  
  # Extract dimensions to plot
  plot_data <- as.data.frame(Z[, dims])
  colnames(plot_data) <- c("Dim1", "Dim2")
  
  # Add labels if provided
  if (!is.null(labels)) {
    plot_data$Label <- as.factor(labels)
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dim1, y = Dim2, color = Label)) +
      ggplot2::geom_point()
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dim1, y = Dim2)) +
      ggplot2::geom_point()
  }
  
  p + ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste("Dimension", dims[1]),
      y = paste("Dimension", dims[2]),
      title = "ZIFA Projection"
    )
}
