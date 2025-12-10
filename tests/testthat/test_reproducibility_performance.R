context("Reproducibility Tests")

test_that("results are seed reproducible", {
  N <- 30
  D <- 20
  K <- 2
  
  # Run 1
  set.seed(42)
  Y1 <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y1 <- log2(Y1 + 1)
  Y1[sample(length(Y1), floor(length(Y1) * 0.3))] <- 0
  
  set.seed(123)
  result1 <- fit_zifa(Y1, k = K, verbose = FALSE)
  
  # Run 2 with same seeds
  set.seed(42)
  Y2 <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y2 <- log2(Y2 + 1)
  Y2[sample(length(Y2), floor(length(Y2) * 0.3))] <- 0
  
  set.seed(123)
  result2 <- fit_zifa(Y2, k = K, verbose = FALSE)
  
  expect_equal(result1$Z, result2$Z, tolerance = 1e-6)
  expect_equal(result1$model_params$A, result2$model_params$A, tolerance = 1e-6)
})

test_that("block and non-block give similar subspaces", {
  skip_on_cran()  # Skip on CRAN due to time
  
  set.seed(123)
  N <- 30
  D <- 100
  K <- 2
  
  Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y <- log2(Y + 1)
  Y[sample(length(Y), floor(length(Y) * 0.3))] <- 0
  
  result_block <- fit_block_zifa(Y, k = K, n_blocks = 2, verbose = FALSE)
  result_normal <- fit_zifa(Y, k = K, verbose = FALSE)
  
  # Check subspace alignment via correlation
  # Note: signs may flip, so compare absolute correlation
  for (k in 1:K) {
    cor_val <- abs(cor(result_block$Z[, k], result_normal$Z[, k]))
    expect_true(cor_val > 0.5)  # Should have reasonable alignment
  }
})