test_that("fit_zifa returns a list with expected elements", {
  set.seed(123)
  # Simulate data (samples x genes, then transpose for fit_zifa)
  N <- 50
  D <- 20
  K <- 2
  
  # Generate data
  Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y <- log2(Y + 1)  # Log-transform
  Y[sample(length(Y), floor(length(Y) * 0.3))] <- 0  # Add zeros
  
  result <- fit_zifa(Y, k = K, verbose = FALSE)
  
  expect_true(is.list(result))
  expect_true("Z" %in% names(result))
  expect_true("model_params" %in% names(result))
  expect_true("n_iter" %in% names(result))
  expect_equal(ncol(result$Z), K)
  expect_equal(nrow(result$Z), N)
  
  # Check model parameters
  expect_true(all(result$model_params$sigmas > 0))
  expect_true(result$model_params$decay_coef > 0)
  expect_equal(length(result$model_params$mus), D)
})

test_that("fit_zifa handles a matrix with all zeros", {
  Y <- matrix(0, nrow = 100, ncol = 20)
  expect_error(fit_zifa(Y, k = 2, verbose = FALSE), "Input matrix contains only zeros")
})

test_that("fit_zifa rejects negative values", {
  Y <- matrix(rnorm(200), nrow = 20, ncol = 10)  # Contains negatives
  expect_error(fit_zifa(Y, k = 2, verbose = FALSE), "negative values")
})

test_that("fit_block_zifa returns results with expected dimensions", {
  set.seed(123)
  N <- 30
  D <- 100
  K <- 2
  
  Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y <- log2(Y + 1)
  Y[sample(length(Y), floor(length(Y) * 0.3))] <- 0
  
  result_block <- fit_block_zifa(Y, k = K, n_blocks = 2, verbose = FALSE)
  expect_true(is.list(result_block))
  expect_equal(ncol(result_block$Z), K)
  expect_equal(nrow(result_block$Z), N)
})

test_that("S3 method for SummarizedExperiment works for fit_zifa", {
  skip_if_not_installed("SummarizedExperiment")
  
  library(SummarizedExperiment)
  set.seed(123)
  N <- 30
  D <- 20
  K <- 2
  
  Y <- matrix(rpois(D * N, lambda = 5), nrow = D, ncol = N)
  Y <- log2(Y + 1)
  Y[sample(length(Y), floor(length(Y) * 0.3))] <- 0
  
  se <- SummarizedExperiment(assays = list(counts = Y))
  
  result_se <- fit_zifa(se, k = K, verbose = FALSE)
  expect_true(is.list(result_se))
  expect_equal(ncol(result_se$Z), K)
})

test_that("preprocess_data works correctly", {
  counts <- matrix(c(0, 1, 10, 100), nrow = 2)
  # Matrix layout:
  #      [,1] [,2]
  # [1,]    0   10
  # [2,]    1  100
  
  # Default: log2(counts + 1)
  result <- preprocess_data(counts)
  expect_equal(result[1, 1], 0)  # log2(0 + 1) = 0
  expect_equal(result[2, 1], 1)  # log2(1 + 1) = 1
  
  # Natural log
  result_ln <- preprocess_data(counts, base = exp(1))
  expect_equal(result_ln[1, 1], 0)  # ln(0 + 1) = 0
})

test_that("preprocess_data rejects negative counts", {
  counts <- matrix(c(-1, 1, 2, 3), nrow = 2)
  expect_error(preprocess_data(counts), "negative values")
})

test_that("input validation catches all-zero columns", {
  Y <- matrix(rpois(200, 5), nrow = 20, ncol = 10)
  Y <- log2(Y + 1)
  Y[5, ] <- 0  # Make ROW 5 (gene 5) all zeros - input is genes x samples
  
  expect_error(fit_zifa(Y, k = 2, verbose = FALSE), "entirely zero")
})

test_that("plot_zifa returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  
  Z <- matrix(rnorm(100), ncol = 2)
  p <- plot_zifa(Z)
  expect_s3_class(p, "ggplot")
  
  # With labels
  labels <- rep(c("A", "B"), each = 25)
  p_labeled <- plot_zifa(Z, labels = labels)
  expect_s3_class(p_labeled, "ggplot")
})

test_that("plot_zifa validates dimensions", {
  skip_if_not_installed("ggplot2")
  
  Z <- matrix(rnorm(100), ncol = 2)
  expect_error(plot_zifa(Z, dims = c(2, 3)), "only has 2 dimensions")
})
