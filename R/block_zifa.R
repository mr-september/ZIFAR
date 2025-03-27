#' Fit Block ZIFA Model
#'
#' Implements Block Zero-Inflated Factor Analysis for more efficient processing
#' of large datasets by dividing genes into blocks.
#'
#' @param Y Matrix of expression data (genes as rows, cells as columns)
#' @param k Integer specifying the number of latent dimensions
#' @param n_blocks Number of blocks to divide genes into (default: genes/500)
#' @param ... Additional parameters passed to fit_zifa
#'
#' @return A list containing:
#'   \item{Z}{Matrix of low-dimensional projections (cells x k dimensions)}
#'   \item{model_params}{List of model parameters}
#'
#' @examples
#' # Generate large simulated dataset
#' set.seed(123)
#' true_Z <- matrix(rnorm(100), ncol = 2)
#' lambda <- matrix(runif(1000*2, -1, 1), nrow = 1000)
#' mu <- true_Z %*% t(lambda)
#' Y <- mu + matrix(rnorm(1000*50), nrow = 1000)
#' # Add zero-inflation
#' mask <- matrix(rbinom(1000*50, 1, prob = 0.3), nrow = 1000)
#' Y[mask == 1] <- 0
#'
#' # Run Block ZIFA
#' result <- fit_block_zifa(Y, k = 2, n_blocks = 2)
#'
#' @export
fit_block_zifa <- function(Y, k, n_blocks = NULL, ...) {

  if (inherits(Y, "SummarizedExperiment")) {
    Y <- SummarizedExperiment::assay(Y)
  }

  n_genes <- nrow(Y)
  n_samples <- ncol(Y)
  
  # Determine number of blocks if not specified
  if (is.null(n_blocks)) {
    n_blocks <- max(1, floor(n_genes / 500))
  }
  
  # Divide genes into blocks
  block_size <- ceiling(n_genes / n_blocks)
  blocks <- list()
  for (b in 1:n_blocks) {
    start_idx <- (b - 1) * block_size + 1
    end_idx <- min(b * block_size, n_genes)
    block_indices <- start_idx:end_idx  # Define block_indices here
    
    if (length(block_indices) == 0) {
      warning("Empty block encountered, skipping")
      next
    }
    
    blocks[[b]] <- block_indices
  }
  
  # Initialize results
  Z_sum <- matrix(0, nrow = n_samples, ncol = k)
  weights_sum <- matrix(0, nrow = n_samples, ncol = k)
  
  # Process each block
  for (b in 1:n_blocks) {
    cat("Processing block", b, "of", n_blocks, "\n")
    
    # Extract genes for this block
    Y_block <- Y[blocks[[b]], , drop = FALSE]
    
    # Run ZIFA on this block
    block_result <- fit_zifa(Y_block, k, ...)
    
    # Accumulate results (weighted by precision)
    Z_sum <- Z_sum + block_result$Z
    weights_sum <- weights_sum + 1
  }
  
  # Compute final Z as weighted average
  Z_final <- Z_sum / weights_sum
  
  # Final full run using Z_final as initialization
  final_result <- fit_zifa(Y, k, Z_init = Z_final, ...)
  
  return(final_result)
}

#' Fit Block ZIFA Model for SummarizedExperiment Objects
#'
#' This method extracts the primary assay data from a SummarizedExperiment
#' object and fits the block-wise ZIFA model.
#'
#' @param se A SummarizedExperiment object.
#' @param k Integer specifying the number of latent dimensions.
#' @param ... Additional parameters passed to fit_block_zifa.
#'
#' @return The result from fit_block_zifa.
#' @export
fit_block_zifa.SummarizedExperiment <- function(se, k, ...) {
  # Convert SummarizedExperiment input to a matrix if needed
  if (inherits(se, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package required but not installed")
    }
    Y <- SummarizedExperiment::assay(se) 
  } else if (!is.matrix(se)) {
    stop("Input se must be a matrix or SummarizedExperiment object")
  }
  
  # Fit the ZIFA model on the extracted data
  result <- fit_block_zifa(Y, k, ...)
  return(result)
}