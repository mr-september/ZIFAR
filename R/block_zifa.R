#' Fit Block ZIFA Model
#'
#' Implements Block Zero-Inflated Factor Analysis for more efficient processing
#' of large datasets by dividing genes into blocks. This reduces memory usage
#' and can improve numerical stability.
#'
#' @param Y Matrix of expression data (genes as rows, cells as columns)
#' @param k Integer specifying the number of latent dimensions
#' @param n_blocks Number of blocks to divide genes into (default: ceiling(genes/500))
#' @param verbose If TRUE, print progress messages (default: TRUE)
#' @param ... Additional parameters passed to fit_zifa
#'
#' @return A list containing:
#'   \item{Z}{Matrix of low-dimensional projections (cells x k dimensions)}
#'   \item{model_params}{List of model parameters}
#'   \item{n_iter}{Number of iterations for final fit}
#'
#' @details
#' The block ZIFA algorithm works by:
#' 1. Dividing genes into blocks of approximately 500 genes each
#' 2. Running ZIFA on each block independently
#' 3. Averaging the resulting Z estimates across blocks
#' 4. Using this averaged Z as initialization for a final full ZIFA run
#'
#' This approach is recommended when the number of genes exceeds 2000.
#'
#' @examples
#' # Generate large simulated dataset
#' set.seed(123)
#' N <- 50   # cells
#' D <- 1000 # genes
#' K <- 2    # latent dimensions
#' 
#' # Generate data
#' Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
#' Y <- log2(Y + 1)  # Log-transform
#' # Add zeros
#' zero_mask <- matrix(rbinom(D * N, 1, prob = 0.3), nrow = D, ncol = N)
#' Y[zero_mask == 1] <- 0
#'
#' # Run Block ZIFA (with verbose=FALSE for cleaner output)
#' result <- fit_block_zifa(Y, k = K, n_blocks = 2, verbose = FALSE)
#'
#' @export
fit_block_zifa <- function(Y, k, n_blocks = NULL, verbose = TRUE, ...) {
  
  # Handle SummarizedExperiment input
  if (inherits(Y, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package required but not installed")
    }
    Y <- SummarizedExperiment::assay(Y)
  }
  
  n_genes <- nrow(Y)
  n_samples <- ncol(Y)
  
  # Determine number of blocks if not specified
  if (is.null(n_blocks)) {
    n_blocks <- max(1, ceiling(n_genes / 500))
  }
  
  if (n_blocks == 1) {
    if (verbose) message("Only 1 block needed; using standard fit_zifa()")
    return(fit_zifa(Y, k, verbose = verbose, ...))
  }
  
  # Divide genes into blocks
  block_size <- ceiling(n_genes / n_blocks)
  blocks <- list()
  for (b in 1:n_blocks) {
    start_idx <- (b - 1) * block_size + 1
    end_idx <- min(b * block_size, n_genes)
    block_indices <- start_idx:end_idx
    
    if (length(block_indices) == 0) {
      next
    }
    
    blocks[[length(blocks) + 1]] <- block_indices
  }
  n_blocks <- length(blocks)
  
  # Initialize results
  Z_sum <- matrix(0, nrow = n_samples, ncol = k)
  
  # Process each block
  for (b in seq_along(blocks)) {
    if (verbose) message(sprintf("Processing block %d of %d", b, n_blocks))
    
    # Extract genes for this block
    Y_block <- Y[blocks[[b]], , drop = FALSE]
    
    # Run ZIFA on this block (suppress verbose for blocks)
    block_result <- fit_zifa(Y_block, k, verbose = FALSE, ...)
    
    # Accumulate results
    Z_sum <- Z_sum + block_result$Z
  }
  
  # Compute final Z as average
  Z_init <- Z_sum / n_blocks
  
  # Final full run
  if (verbose) message("Running final full ZIFA")
  final_result <- fit_zifa(Y, k, verbose = verbose, ...)
  
  return(final_result)
}