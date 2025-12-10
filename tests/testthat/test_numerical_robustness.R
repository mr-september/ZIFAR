context("Numerical Robustness Tests")

test_that("handles sparse data", {
  set.seed(123)
  N <- 30
  D <- 20
  
  Y <- matrix(rpois(D * N, lambda = 2), nrow = D, ncol = N)
  Y <- log2(Y + 1)
  # Make very sparse (80% zeros)
  Y[sample(length(Y), floor(length(Y) * 0.8))] <- 0
  
  # Should still run without error
  expect_no_error(fit_zifa(Y, k = 2, max_iter = 20, verbose = FALSE))
})

test_that("converges with different tolerances", {
  set.seed(123)
  N <- 30
  D <- 20
  
  Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y <- log2(Y + 1)
  Y[sample(length(Y), floor(length(Y) * 0.3))] <- 0
  
  # Test loose convergence
  result_loose <- fit_zifa(Y, k = 2, tol = 1e-1, verbose = FALSE)
  # Test tight convergence 
  result_tight <- fit_zifa(Y, k = 2, tol = 1e-4, verbose = FALSE)
  
  # Loose should converge in fewer iterations

  expect_true(result_loose$n_iter <= result_tight$n_iter)
})

test_that("single_sigma option works", {
  set.seed(123)
  N <- 30
  D <- 20
  
  Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y <- log2(Y + 1)
  Y[sample(length(Y), floor(length(Y) * 0.3))] <- 0
  
  result <- fit_zifa(Y, k = 2, single_sigma = TRUE, verbose = FALSE)
  
  # All sigmas should be identical
  expect_equal(length(unique(result$model_params$sigmas)), 1)
})

test_that("handles small datasets", {
  set.seed(123)
  N <- 10
  D <- 5
  
  Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y <- log2(Y + 1)
  Y[sample(length(Y), floor(length(Y) * 0.3))] <- 0
  
  # Should work with k = 1

  result <- fit_zifa(Y, k = 1, verbose = FALSE)
  expect_equal(ncol(result$Z), 1)
})