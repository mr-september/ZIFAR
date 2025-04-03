test_that("fit_zifa returns a list with expected elements", {
  set.seed(123)
  # Simulate data
  true_Z <- matrix(rnorm(100), ncol = 2)
  lambda <- matrix(runif(500 * 2, -1, 1), nrow = 500)
  mu <- lambda %*% t(true_Z)
  Y <- mu + matrix(rnorm(500 * 50, sd = 0.5), nrow = 500) + matrix(rnorm(500 * 50, mean = 0, sd = 1e-6), nrow = 500)
  # Introduce zero-inflation
  mask <- matrix(rbinom(500 * 50, 1, prob = 0.3), nrow = 500)
  Y[mask == 1] <- 0
  
  result <- fit_zifa(Y, k = 2)
  
  expect_true(is.list(result))
  expect_true("Z" %in% names(result))
  expect_true("model_params" %in% names(result))
  expect_true("ll" %in% names(result))
  expect_equal(ncol(result$Z), 2)
  
  # Enhanced verification
  expect_true(all(result$model_params$sigma > 0))
  expect_true(result$model_params$decay_coef > 0)
  expect_true(all(diff(result$ll) >= -1e-6)) # Allow small numerical decreases
})

test_that("properly handles zero-inflation", {
  set.seed(123)
  # Create data with known zero-inflation pattern
  true_Z <- matrix(rnorm(100), ncol = 2)
  lambda <- matrix(runif(500 * 2, -1, 1), nrow = 500)
  mu <- lambda %*% t(true_Z)
  Y <- mu + matrix(rnorm(500 * 50, sd = 0.5), nrow = 500)
  
  # Add structured zeros (higher probability where mu is small)
  zero_probs <- 1/(1 + exp(2*mu))
  mask <- matrix(rbinom(500 * 50, 1, prob = zero_probs), nrow = 500)
  Y[mask == 1] <- 0
  
  result <- fit_zifa(Y, k = 2)
  
  # Verify zeros affected the result
  cor_zero <- cor(as.vector(mask), abs(as.vector(result$Z %*% result$model_params$lambda)))
  expect_true(cor_zero > 0.3) # Some positive correlation expected
})

test_that("fit_zifa handles a matrix with all zeros", {
  Y <- matrix(0, nrow = 100, ncol = 20)
  expect_error(fit_zifa(Y, k = 2), "Input matrix contains only zeros")
})

test_that("fit_block_zifa returns results with expected dimensions", {
  set.seed(123)
  true_Z <- matrix(rnorm(100), ncol = 2)
  lambda <- matrix(runif(500 * 2, -1, 1), nrow = 500)
  mu <- lambda %*% t(true_Z)
  Y <- mu + matrix(rnorm(500 * 50, sd = 0.5), nrow = 500) + matrix(rnorm(500 * 50, mean = 0, sd = 1e-6), nrow = 500)
  mask <- matrix(rbinom(500 * 50, 1, prob = 0.3), nrow = 500)
  Y[mask == 1] <- 0
  
  result_block <- fit_block_zifa(Y, k = 2)
  expect_true(is.list(result_block))
  expect_equal(ncol(result_block$Z), 2)
})

test_that("S3 method for SummarizedExperiment works for fit_zifa", {
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    library(SummarizedExperiment)
    set.seed(123)
    true_Z <- matrix(rnorm(100), ncol = 2)
    lambda <- matrix(runif(500 * 2, -1, 1), nrow = 500)
    mu <- lambda %*% t(true_Z)
    Y <- mu + matrix(rnorm(500 * 50, sd = 0.5), nrow = 500) + matrix(rnorm(500 * 50, mean = 0, sd = 1e-6), nrow = 500)
    mask <- matrix(rbinom(500 * 50, 1, prob = 0.3), nrow = 500)
    Y[mask == 1] <- 0
    se <- SummarizedExperiment(assays = list(counts = Y))
    
    result_se <- fit_zifa(se, k = 2)
    expect_true(is.list(result_se))
    expect_equal(ncol(result_se$Z), 2)
  } else {
    skip("SummarizedExperiment package not installed")
  }
})

test_that("S3 method for SummarizedExperiment works for fit_block_zifa", {
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    library(SummarizedExperiment)
    set.seed(123)
    true_Z <- matrix(rnorm(100), ncol = 2)
    lambda <- matrix(runif(500 * 2, -1, 1), nrow = 500)
    mu <- lambda %*% t(true_Z)
    Y <- mu + matrix(rnorm(500 * 50, sd = 0.5), nrow = 500) + matrix(rnorm(500 * 50, mean = 0, sd = 1e-6), nrow = 500)
    mask <- matrix(rbinom(500 * 50, 1, prob = 0.3), nrow = 500)
    Y[mask == 1] <- 0
    se <- SummarizedExperiment(assays = list(counts = Y))
    
    result_se <- fit_block_zifa(se, k = 2)
    expect_true(is.list(result_se))
    expect_equal(ncol(result_se$Z), 2)
  } else {
    skip("SummarizedExperiment package not installed")
  }
})

test_that("compute_log_likelihood returns a numeric scalar", {
  set.seed(123)
  true_Z <- matrix(rnorm(100), ncol = 2)
  lambda <- matrix(runif(500 * 2, -1, 1), nrow = 500)
  mu <- lambda %*% t(true_Z)
  Y <- mu + matrix(rnorm(500 * 50, sd = 0.5), nrow = 500) + matrix(rnorm(500 * 50, mean = 0, sd = 1e-6), nrow = 500)
  mask <- matrix(rbinom(500 * 50, 1, prob = 0.3), nrow = 500)
  Y[mask == 1] <- 0
  
  result <- fit_zifa(Y, k = 2)
  ll_val <- compute_log_likelihood(Y, result$Z, result$model_params$lambda, result$model_params$sigma, result$model_params$decay_coef)
  
  expect_true(is.numeric(ll_val))
  expect_equal(length(ll_val), 1)
})

test_that("invalid input raises an error", {
  expect_error(fit_zifa("not a matrix", k = 2))
  expect_error(fit_block_zifa(12345, k = 2))
})
