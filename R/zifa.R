#' Fit ZIFA Model
#'
#' Implements Zero-Inflated Factor Analysis for dimensionality reduction
#' of zero-inflated data, especially single-cell gene expression data.
#'
#' @param Y Matrix of expression data (genes as rows, cells as columns)
#' @param k Integer specifying the number of latent dimensions
#' @param iterations Maximum number of iterations for the EM algorithm
#' @param convergence_criterion Convergence threshold for log-likelihood
#' @param lambda_init Initial value for lambda parameter
#' @param sigma_init Initial value for sigma parameter
#' @param decay_coef Initial value for zero decay coefficient
#'
#' @return A list containing:
#'   \item{Z}{Matrix of low-dimensional projections (cells x k dimensions)}
#'   \item{model_params}{List of model parameters including lambda, sigma, and decay_coef}
#'   \item{ll}{Vector of log-likelihood values per iteration}
#'
#' @examples
#' # Generate simulated data with dropout zeros
#' set.seed(123)
#' true_Z <- matrix(rnorm(100), ncol = 2)
#' lambda <- matrix(runif(20*2, -1, 1), nrow = 20)
#' mu <- true_Z %*% t(lambda)
#' Y <- mu + matrix(rnorm(20*50), nrow = 20)
#' # Add zero-inflation
#' mask <- matrix(rbinom(20*50, 1, prob = 0.3), nrow = 20)
#' Y[mask == 1] <- 0
#'
#' # Run ZIFA
#' result <- fit_zifa(Y, k = 2)
#'
#' @export
fit_zifa <- function(Y, k, iterations = 100, 
                    convergence_criterion = 1e-5,
                    lambda_init = NULL,
                    sigma_init = NULL,
                    decay_coef = 1.0) {
  
  n_genes <- nrow(Y)
  n_samples <- ncol(Y)
  
  # Initialize parameters using standard factor analysis if not provided
  if (is.null(lambda_init) || is.null(sigma_init)) {
    fa_result <- stats::factanal(t(Y), factors = k, scores = "regression")
    lambda_init <- t(fa_result$loadings)
    sigma_init <- rep(1, n_genes)
  }
  
  # Initialize model parameters
  lambda <- lambda_init
  sigma <- sigma_init
  Z <- matrix(0, nrow = n_samples, ncol = k)
  
  # Track log-likelihood
  ll_trace <- numeric(iterations)
  
  # EM algorithm
  for (iter in 1:iterations) {
    # E-step: estimate posterior of Z
    Z_params <- compute_posterior_Z(Y, lambda, sigma, decay_coef)
    Z <- Z_params$mean
    
    # M-step: update parameters
    params <- update_parameters(Y, Z, Z_params$covariance, decay_coef)
    lambda <- params$lambda
    sigma <- params$sigma
    decay_coef <- params$decay_coef
    
    # Calculate log-likelihood
    ll <- compute_log_likelihood(Y, Z, lambda, sigma, decay_coef)
    ll_trace[iter] <- ll
    
    # Check convergence
    if (iter > 1 && abs(ll_trace[iter] - ll_trace[iter-1]) < convergence_criterion) {
      break
    }
  }
  
  # Return results
  list(
    Z = Z,
    model_params = list(
      lambda = lambda,
      sigma = sigma,
      decay_coef = decay_coef
    ),
    ll = ll_trace[1:iter]
  )
}

#' Compute posterior distribution of Z given Y
#'
#' @param Y Data matrix
#' @param lambda Loading matrix
#' @param sigma Vector of variances
#' @param decay_coef Zero-inflation decay coefficient
#'
#' @return List with posterior mean and covariance
#' @keywords internal
compute_posterior_Z <- function(Y, lambda, sigma, decay_coef) {
  # Implementation of the computation of the posterior
  # of Z given the non-zero entries of Y
  n_samples <- ncol(Y)
  k <- nrow(lambda)
  
  Z_mean <- matrix(0, nrow = n_samples, ncol = k)
  Z_cov <- array(0, dim = c(k, k, n_samples))
  
  # Process each sample
  for (i in 1:n_samples) {
    # Get non-zero indices
    nonzero_idx <- which(Y[, i] != 0)
    
    if (length(nonzero_idx) > 0) {
      # Extract relevant parts of lambda and Y
      lambda_nz <- lambda[, nonzero_idx, drop = FALSE]
      sigma_nz <- sigma[nonzero_idx]
      Y_nz <- Y[nonzero_idx, i]
      
      # Compute posterior covariance
      precision_matrix <- diag(k) + 
        lambda_nz %*% (t(lambda_nz) / sigma_nz)
      
      Z_cov[,,i] <- solve(precision_matrix)
      
      # Compute posterior mean
      Z_mean[i,] <- Z_cov[,,i] %*% 
        (lambda_nz %*% (Y_nz / sigma_nz))
    } else {
      # All zeros case
      Z_cov[,,i] <- diag(k)
      Z_mean[i,] <- rep(0, k)
    }
  }
  
  list(mean = Z_mean, covariance = Z_cov)
}

#' Update model parameters
#'
#' @param Y Data matrix
#' @param Z_mean Posterior means of Z
#' @param Z_cov Posterior covariances of Z
#' @param decay_coef Current decay coefficient
#'
#' @return List of updated parameters
#' @keywords internal
update_parameters <- function(Y, Z_mean, Z_cov, decay_coef) {
  n_genes <- nrow(Y)
  n_samples <- ncol(Y)
  k <- ncol(Z_mean)
  
  # Update lambda
  lambda <- matrix(0, nrow = k, ncol = n_genes)
  
  for (j in 1:n_genes) {
    # Indices where Y is non-zero for gene j
    nonzero_idx <- which(Y[j,] != 0)
    
    if (length(nonzero_idx) > 0) {
      # Extract relevant parts
      Y_nz <- Y[j, nonzero_idx]
      Z_nz <- Z_mean[nonzero_idx, , drop = FALSE]
      
      # Create matrices for the update
      A <- matrix(0, nrow = k, ncol = k)
      b <- rep(0, k)
      
      for (i in nonzero_idx) {
        A <- A + Z_mean[i,] %*% t(Z_mean[i,]) + Z_cov[,,i]
        b <- b + Z_mean[i,] * Y[j,i]
      }
      
      # Update lambda for gene j
      lambda[,j] <- solve(A) %*% b
    }
  }
  
  # Update sigma
  sigma <- rep(0, n_genes)
  for (j in 1:n_genes) {
    nonzero_idx <- which(Y[j,] != 0)
    if (length(nonzero_idx) > 0) {
      Y_nz <- Y[j, nonzero_idx]
      Z_nz <- Z_mean[nonzero_idx, , drop = FALSE]
      lambda_j <- lambda[,j]
      
      residuals <- Y_nz - Z_nz %*% lambda_j
      sigma[j] <- mean(residuals^2)
    } else {
      sigma[j] <- 1.0  # Default when all observations are zero
    }
  }
  
  # Update decay coefficient using numerical optimization
  decay_opt <- optimize(function(d) {
    -decay_log_likelihood(Y, Z_mean, lambda, sigma, d)
  }, c(0.1, 10))
  
  updated_decay <- decay_opt$minimum
  
  list(
    lambda = lambda,
    sigma = sigma,
    decay_coef = updated_decay
  )
}

#' Compute log-likelihood for current parameter estimates
#'
#' @param Y Data matrix
#' @param Z Mean of latent variables
#' @param lambda Loading matrix
#' @param sigma Vector of variances
#' @param decay_coef Zero-inflation decay coefficient
#'
#' @return Log-likelihood value
#' @keywords internal
compute_log_likelihood <- function(Y, Z, lambda, sigma, decay_coef) {
  n_genes <- nrow(Y)
  n_samples <- ncol(Y)
  ll <- 0
  
  for (i in 1:n_samples) {
    for (j in 1:n_genes) {
      mu <- sum(Z[i,] * lambda[,j])
      if (Y[j,i] == 0) {
        # Zero component
        p_zero <- compute_zero_probability(mu, decay_coef)
        p_nz_zero <- dnorm(0, mean = mu, sd = sqrt(sigma[j]))
        ll <- ll + log(p_zero + (1 - p_zero) * p_nz_zero)
      } else {
        # Non-zero component
        p_zero <- compute_zero_probability(mu, decay_coef)
        ll <- ll + log(1 - p_zero) + 
          dnorm(Y[j,i], mean = mu, sd = sqrt(sigma[j]), log = TRUE)
      }
    }
  }
  
  return(ll)
}

#' Compute zero probability
#'
#' @param mu Mean of the normal component
#' @param decay_coef Decay coefficient
#'
#' @return Probability of zero inflation
#' @keywords internal
compute_zero_probability <- function(mu, decay_coef) {
  return(exp(-decay_coef * mu^2))
}

#' Log-likelihood for the decay coefficient
#'
#' @param Y Data matrix
#' @param Z Mean of latent variables
#' @param lambda Loading matrix
#' @param sigma Vector of variances
#' @param decay_coef Zero-inflation decay coefficient
#'
#' @return Log-likelihood value for decay coefficient
#' @keywords internal
decay_log_likelihood <- function(Y, Z, lambda, sigma, decay_coef) {
  n_genes <- nrow(Y)
  n_samples <- ncol(Y)
  ll <- 0
  
  for (i in 1:n_samples) {
    for (j in 1:n_genes) {
      mu <- sum(Z[i,] * lambda[,j])
      p_zero <- compute_zero_probability(mu, decay_coef)
      
      if (Y[j,i] == 0) {
        ll <- ll + log(p_zero)
      } else {
        ll <- ll + log(1 - p_zero)
      }
    }
  }
  
  return(ll)
}
