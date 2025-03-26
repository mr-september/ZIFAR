test_that("ZIFA produces correct dimensionality reduction", {
  set.seed(123)
  
  # Generate simulated data
  n_genes <- 20
  n_samples <- 30
  k_true <- 2
  
  # Create true latent space
  Z_true <- matrix(rnorm(n_samples * k_true), ncol = k_true)
  
  # Create loading matrix
  lambda_true <- matrix(runif(n_genes * k_true, -1, 1), nrow = k_true)
  
  # Generate data
  mu <- Z_true %*% lambda_true
  Y <- t(mu) + matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  
  # Add zero-inflation
  mask <- matrix(rbinom(n_genes * n_samples, 1, prob = 0.2), nrow = n_genes)
  Y[mask == 1] <- 0
  
  # Fit ZIFA
  result <- fit_zifa(Y, k = k_true, iterations = 5)
  
  # Test output format
  expect_type(result, "list")
  expect_equal(names(result), c("Z", "model_params", "ll"))
  
  # Test dimensions
  expect_equal(dim(result$Z), c(n_samples, k_true))
  
  # Test that log-likelihood increases
  expect_true(result$ll[length(result$ll)] > result$ll[1])
})
