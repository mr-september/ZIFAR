#' @title Zero-Inflated Factor Analysis (ZIFA)
#' @description Native R implementation of Zero-Inflated Factor Analysis for
#'   dimensionality reduction of single-cell gene expression data.
#'   Based on Pierson & Yau (2015) <doi:10.1186/s13059-015-0805-z>.
#' @name zifa
#' @keywords internal
"_PACKAGE"

# =============================================================================
# Input Validation
# =============================================================================

#' Validate input data for ZIFA
#' @param Y Data matrix (samples x genes)
#' @return NULL invisibly if valid, otherwise throws error
#' @keywords internal
validate_input <- function(Y) {
  if (all(Y == 0)) {
    stop("Input matrix contains only zeros")
  }
  
  if (any(Y < 0)) {
    stop("Input matrix contains negative values. ZIFA expects log-transformed counts.")
  }
  
  if (all(Y == floor(Y)) && max(Y) > 10) {
    warning("Input matrix contains only integers. Ensure data is log-transformed, not raw counts.")
  }
  
  col_zeros <- colSums(abs(Y) < 1e-6) == nrow(Y)
  if (any(col_zeros)) {
    stop("Input matrix has genes (columns) which are entirely zero. Please filter these out.")
  }
  
  if (sum(abs(Y) < 1e-6) == 0) {
    warning("Input matrix contains no zeros. This is unusual for single-cell data.")
  }
  
  invisible(NULL)
}

# =============================================================================
# Core ZIFA Algorithm (Simplified but Correct)
# =============================================================================

#' E-step: Compute expected values of latent variables
#' 
#' Computes posterior expectations E[Z|Y], E[ZZ'|Y], and handles zero entries
#' by computing expectations over the latent X values at zero positions.
#' 
#' @param Y Data matrix (N samples x D genes)
#' @param A Loading matrix (D x K)
#' @param mus Gene means (D x 1)
#' @param sigmas Gene noise std devs (D x 1)
#' @param decay_coef Zero-inflation decay coefficient
#' @return List with EZ, EZZT, EX, EX2
#' @keywords internal
Estep <- function(Y, A, mus, sigmas, decay_coef) {
  N <- nrow(Y)
  D <- ncol(Y)
  K <- ncol(A)
  
  # Tolerance for numerical stability
  eps <- 1e-6
  
  EZ <- matrix(0, N, K)
  EZZT <- array(0, dim = c(N, K, K))
  EX <- Y  # For non-zero entries, EX = Y
  EX2 <- Y^2  # For non-zero entries, EX2 = Y^2
  
  # Prior covariance of Z is I_K
  # For each sample, compute posterior of Z given observed Y
  
  for (i in 1:N) {
    y_i <- Y[i, ]
    is_zero <- abs(y_i) < eps
    is_nonzero <- !is_zero
    
    # For non-zero entries, we observe Y directly
    # For zero entries, the true X is drawn from a truncated/modified distribution
    
    if (sum(is_nonzero) > 0) {
      # Subset to non-zero genes
      A_obs <- A[is_nonzero, , drop = FALSE]
      sigma_obs <- sigmas[is_nonzero]
      mu_obs <- mus[is_nonzero]
      y_obs <- y_i[is_nonzero]
      
      # Posterior of Z given Y_obs (non-zero entries)
      # Prior: Z ~ N(0, I_K)
      # Likelihood: Y_j | Z ~ N(A_j' Z + mu_j, sigma_j^2) for j in observed
      
      # Precision-weighted posterior
      precision <- diag(K)  # Prior precision
      weighted_sum <- rep(0, K)
      
      for (j in seq_len(sum(is_nonzero))) {
        a_j <- A_obs[j, ]
        s2_j <- sigma_obs[j]^2 + eps
        precision <- precision + outer(a_j, a_j) / s2_j
        weighted_sum <- weighted_sum + a_j * (y_obs[j] - mu_obs[j]) / s2_j
      }
      
      # Posterior covariance and mean
      # Add regularization for numerical stability
      precision <- precision + eps * diag(K)
      Sigma_post <- tryCatch(
        solve(precision),
        error = function(e) {
          # Fallback: use regularized pseudoinverse
          eigen_decomp <- eigen(precision, symmetric = TRUE)
          vals <- pmax(eigen_decomp$values, eps)
          eigen_decomp$vectors %*% diag(1/vals) %*% t(eigen_decomp$vectors)
        }
      )
      mu_post <- as.vector(Sigma_post %*% weighted_sum)
      
      EZ[i, ] <- mu_post
      EZZT[i, , ] <- Sigma_post + outer(mu_post, mu_post)
      
    } else {
      # All zeros: posterior is prior (Z ~ N(0, I))
      EZ[i, ] <- rep(0, K)
      EZZT[i, , ] <- diag(K)
    }
    
    # For zero entries, compute E[X|Y=0, Z]
    # Under ZIFA model: P(Y=0|X) = exp(-decay_coef * X^2)
    # This leads to a truncated normal for X given Y=0
    # Approximate: E[X|Y=0] ≈ 0 (mode of the posterior)
    # E[X^2|Y=0] ≈ variance from the posterior
    
    if (sum(is_zero) > 0) {
      # Predicted X at zeros based on current Z estimate
      x_pred <- as.vector(A[is_zero, , drop = FALSE] %*% EZ[i, ]) + mus[is_zero]
      
      # The posterior of X given Y=0 is concentrated near 0
      # Use a simple approximation: E[X] ≈ x_pred * exp(-decay_coef * x_pred^2)
      # E[X^2] ≈ sigma^2 + (E[X])^2
      
      p_zero <- exp(-decay_coef * x_pred^2)
      EX[i, is_zero] <- x_pred * (1 - p_zero)  # Shrink towards 0
      EX2[i, is_zero] <- sigmas[is_zero]^2 + EX[i, is_zero]^2
    }
  }
  
  list(EZ = EZ, EZZT = EZZT, EX = EX, EX2 = EX2)
}

#' M-step: Update model parameters
#' @param Y Data matrix (N x D)
#' @param EZ Expected Z (N x K)
#' @param EZZT Expected ZZ' (N x K x K)
#' @param EX Expected X (N x D)
#' @param EX2 Expected X^2 (N x D)
#' @param old_decay_coef Previous decay coefficient
#' @param single_sigma If TRUE, use single variance
#' @return Updated parameters
#' @keywords internal
Mstep <- function(Y, EZ, EZZT, EX, EX2, old_decay_coef, single_sigma = FALSE) {
  N <- nrow(Y)
  D <- ncol(Y)
  K <- ncol(EZ)
  eps <- 1e-6
  
  is_zero <- abs(Y) < eps
  
  # Use EX for zeros, Y for non-zeros
  X_effective <- ifelse(is_zero, EX, Y)
  
  # Update A and mus jointly
  # Model: X_j = A_j' Z + mu_j + noise
  # Solve: [A_j, mu_j] from regression
  
  A <- matrix(0, D, K)
  mus <- numeric(D)
  sigmas <- numeric(D)
  
  # Build sum of EZZT and sum of EZ for the normal equations
  sum_EZZT <- apply(EZZT, c(2, 3), sum)  # K x K
  sum_EZ <- colSums(EZ)  # K
  
  for (j in 1:D) {
    x_j <- X_effective[, j]
    
    # Normal equations for [A_j; mu_j]
    # [ sum(ZZ')   sum(Z) ] [A_j ] = [sum(Z * X_j)  ]
    # [ sum(Z')    N      ] [mu_j]   [sum(X_j)      ]
    
    B <- rbind(
      cbind(sum_EZZT, sum_EZ),
      c(sum_EZ, N)
    )
    # Add regularization
    diag(B) <- diag(B) + eps
    
    c_vec <- c(colSums(EZ * x_j), sum(x_j))
    
    solution <- tryCatch(
      solve(B, c_vec),
      error = function(e) {
        # Fallback: simple least squares on means
        c(rep(0, K), mean(x_j))
      }
    )
    
    A[j, ] <- solution[1:K]
    mus[j] <- solution[K + 1]
    
    # Update sigma: residual variance
    predicted <- EZ %*% A[j, ] + mus[j]
    
    # For non-zeros: residual = Y - predicted
    # For zeros: use EX2 - 2*EX*predicted + predicted^2
    residuals2 <- ifelse(
      is_zero[, j],
      EX2[, j] - 2 * EX[, j] * predicted + predicted^2,
      (Y[, j] - predicted)^2
    )
    
    sigmas[j] <- sqrt(max(mean(residuals2), eps))
  }
  
  if (single_sigma) {
    sigmas <- rep(mean(sigmas), D)
  }
  
  # Update decay coefficient
  # Objective: maximize sum over zeros of log P(Y=0|X) = -decay * X^2
  # and sum over non-zeros of log(1 - exp(-decay * Y^2))
  
  decay_objective <- function(d) {
    if (d <= 0) return(Inf)
    
    y2 <- Y^2
    ex2 <- EX2
    
    # For zeros: -d * E[X^2]
    # For non-zeros: log(1 - exp(-d * Y^2))
    ll_zero <- -d * ex2[is_zero]
    
    exp_term <- exp(-d * y2[!is_zero])
    ll_nonzero <- log(pmax(1 - exp_term, eps))
    
    -sum(ll_zero) - sum(ll_nonzero)  # Return negative for minimization
  }
  
  opt_result <- tryCatch(
    stats::optimize(decay_objective, interval = c(0.01, 10)),
    error = function(e) list(minimum = old_decay_coef)
  )
  decay_coef <- opt_result$minimum
  
  list(A = A, mus = mus, sigmas = sigmas, decay_coef = decay_coef)
}

#' Initialize ZIFA parameters
#' @param Y Data matrix (N x D)
#' @param K Number of latent dimensions
#' @param single_sigma Use single variance
#' @return List of initial parameters
#' @keywords internal
initialize_params <- function(Y, K, single_sigma = FALSE) {
  N <- nrow(Y)
  D <- ncol(Y)
  eps <- 1e-6
  
  # Compute per-gene means from non-zero values
  mus <- numeric(D)
  for (j in 1:D) {
    non_zero <- abs(Y[, j]) > eps
    if (any(non_zero)) {
      mus[j] <- mean(Y[non_zero, j])
    }
  }
  
  # Center data
  Y_centered <- sweep(Y, 2, mus)
  
  # Initialize via SVD (more robust than factanal)
  svd_result <- tryCatch({
    svd(Y_centered, nu = K, nv = K)
  }, error = function(e) {
    # Fallback: random initialization
    list(
      v = matrix(rnorm(D * K), D, K),
      d = rep(1, K)
    )
  })
  
  # A = V * sqrt(eigenvalues) approximately
  if (length(svd_result$d) >= K) {
    A <- svd_result$v[, 1:K, drop = FALSE] * 
         rep(sqrt(svd_result$d[1:K] / N), each = D)
  } else {
    A <- matrix(rnorm(D * K, sd = 0.1), D, K)
  }
  
  # Estimate noise from residuals
  Z_init <- Y_centered %*% A %*% solve(t(A) %*% A + eps * diag(K))
  residuals <- Y_centered - Z_init %*% t(A)
  sigmas <- pmax(apply(residuals, 2, sd), eps)
  
  if (single_sigma) {
    sigmas <- rep(mean(sigmas), D)
  }
  
  # Estimate decay coefficient from zero proportions
  zero_props <- colMeans(abs(Y) < eps)
  gene_means <- mus
  
  # Fit: zero_prop ≈ exp(-decay * mean^2)
  # Take log: log(zero_prop) ≈ -decay * mean^2
  valid <- zero_props > eps & zero_props < (1 - eps) & abs(gene_means) > eps
  
  if (sum(valid) > 5) {
    decay_coef <- tryCatch({
      fit <- lm(log(zero_props[valid]) ~ 0 + I(gene_means[valid]^2))
      max(-coef(fit), 0.1)
    }, error = function(e) 0.5)
  } else {
    decay_coef <- 0.5
  }
  
  decay_coef <- min(max(decay_coef, 0.1), 5)  # Bound to reasonable range
  
  list(A = A, mus = mus, sigmas = sigmas, decay_coef = decay_coef)
}

#' Compute log-likelihood for monitoring convergence
#' @keywords internal
compute_log_likelihood <- function(Y, EZ, A, mus, sigmas, decay_coef) {
  N <- nrow(Y)
  D <- ncol(Y)
  eps <- 1e-10
  
  # Predicted means
  predicted <- EZ %*% t(A) + matrix(mus, N, D, byrow = TRUE)
  
  is_zero <- abs(Y) < 1e-6
  
  ll <- 0
  
  # Non-zero entries: Gaussian likelihood
  for (j in 1:D) {
    nonzero_idx <- !is_zero[, j]
    if (any(nonzero_idx)) {
      residuals <- Y[nonzero_idx, j] - predicted[nonzero_idx, j]
      ll <- ll + sum(dnorm(residuals, sd = sigmas[j], log = TRUE))
      # Also add log(1 - p_zero) for non-zeros
      p_zero <- exp(-decay_coef * Y[nonzero_idx, j]^2)
      ll <- ll + sum(log(pmax(1 - p_zero, eps)))
    }
  }
  
  # Zero entries: probability of dropout
  for (j in 1:D) {
    zero_idx <- is_zero[, j]
    if (any(zero_idx)) {
      # P(Y=0) ≈ exp(-decay * predicted^2) (approximately)
      ll <- ll + sum(-decay_coef * predicted[zero_idx, j]^2)
    }
  }
  
  ll
}

#' Fit ZIFA Model
#'
#' Implements Zero-Inflated Factor Analysis for dimensionality reduction
#' of zero-inflated data, especially single-cell gene expression data.
#'
#' @param Y Matrix of expression data (genes as rows, cells as columns)
#' @param k Integer specifying the number of latent dimensions
#' @param max_iter Maximum number of EM iterations (default: 100)
#' @param tol Convergence threshold for relative parameter change (default: 1e-2)
#' @param single_sigma If TRUE, estimate one noise variance for all genes (default: FALSE)
#' @param verbose If TRUE, print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{Z}{Matrix of low-dimensional projections (cells x k dimensions)}
#'   \item{model_params}{List of model parameters including A, mus, sigmas, and decay_coef}
#'   \item{ll}{Vector of log-likelihood values per iteration}
#'   \item{n_iter}{Number of iterations until convergence}
#'
#' @examples
#' # Generate simulated data with dropout
#' set.seed(123)
#' N <- 50  # cells
#' D <- 20  # genes
#' 
#' Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
#' Y <- log2(Y + 1)
#' Y[sample(length(Y), floor(length(Y) * 0.3))] <- 0
#' 
#' # Fit ZIFA
#' result <- fit_zifa(Y, k = 2, verbose = FALSE)
#'
#' @export
fit_zifa <- function(Y, k, max_iter = 100, tol = 1e-2, single_sigma = FALSE, verbose = TRUE) {
  
  # Handle SummarizedExperiment input
  if (inherits(Y, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package required but not installed")
    }
    Y <- SummarizedExperiment::assay(Y)
  }
  
  # Transpose to samples x genes for internal processing
  Y <- t(Y)
  N <- nrow(Y)
  D <- ncol(Y)
  
  # Validate input
  validate_input(Y)
  
  if (D > 2000 && verbose) {
    message("Note: Large number of genes. Consider using fit_block_zifa() for efficiency.")
  }
  
  if (verbose) {
    message(sprintf("Running ZIFA with N=%d samples, D=%d genes, K=%d dimensions", N, D, k))
  }
  
  # Initialize parameters
  params <- initialize_params(Y, k, single_sigma)
  A <- params$A
  mus <- params$mus
  sigmas <- params$sigmas
  decay_coef <- params$decay_coef
  
  # Track log-likelihood
  ll_trace <- numeric(max_iter)
  
  # EM algorithm
  n_iter <- 0
  converged <- FALSE
  
  while (n_iter < max_iter && !converged) {
    # E-step
    expectations <- Estep(Y, A, mus, sigmas, decay_coef)
    
    # M-step
    new_params <- Mstep(Y, expectations$EZ, expectations$EZZT, 
                        expectations$EX, expectations$EX2, 
                        decay_coef, single_sigma)
    
    # Compute log-likelihood
    ll_trace[n_iter + 1] <- compute_log_likelihood(
      Y, expectations$EZ, new_params$A, new_params$mus, 
      new_params$sigmas, new_params$decay_coef
    )
    
    # Check convergence based on parameter change
    max_change <- 0
    for (param_pair in list(
      list(new_params$A, A),
      list(new_params$mus, mus),
      list(new_params$sigmas, sigmas),
      list(new_params$decay_coef, decay_coef)
    )) {
      new_val <- param_pair[[1]]
      old_val <- param_pair[[2]]
      denom <- mean(abs(new_val)) + 1e-10
      rel_change <- mean(abs(new_val - old_val)) / denom
      max_change <- max(max_change, rel_change)
    }
    
    if (max_change < tol) {
      converged <- TRUE
      if (verbose) {
        message(sprintf("Converged after %d iterations (max rel change: %.2e)", 
                        n_iter + 1, max_change))
      }
    }
    
    # Update parameters
    A <- new_params$A
    mus <- new_params$mus
    sigmas <- new_params$sigmas
    decay_coef <- new_params$decay_coef
    
    n_iter <- n_iter + 1
  }
  
  if (!converged && verbose) {
    message(sprintf("Maximum iterations (%d) reached", max_iter))
  }
  
  # Final E-step to get Z estimates
  final_expectations <- Estep(Y, A, mus, sigmas, decay_coef)
  
  list(
    Z = final_expectations$EZ,
    model_params = list(
      A = A,
      mus = mus,
      sigmas = sigmas,
      decay_coef = decay_coef
    ),
    ll = ll_trace[1:n_iter],
    n_iter = n_iter
  )
}

#' Compute Zero Probability
#' @param mu Mean of the normal component
#' @param decay_coef Decay coefficient
#' @return Probability of zero inflation
#' @keywords internal
compute_zero_probability <- function(mu, decay_coef) {
  exp(-decay_coef * mu^2)
}