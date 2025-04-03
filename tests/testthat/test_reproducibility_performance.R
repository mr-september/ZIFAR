context("Reproducibility and Performance Tests")

test_that("results are seed reproducible", {
  set.seed(123)
  Y1 <- matrix(rpois(500, 5), nrow=50)
  Y1[sample(length(Y1), 150)] <- 0
  result1 <- fit_zifa(Y1, k=2)
  
  set.seed(123)
  Y2 <- matrix(rpois(500, 5), nrow=50) 
  Y2[sample(length(Y2), 150)] <- 0
  result2 <- fit_zifa(Y2, k=2)
  
  expect_equal(result1$Z, result2$Z, tolerance=1e-6)
})

test_that("block and non-block give same results", {
  set.seed(123)
  Y <- matrix(rpois(1000, 5), nrow=100)
  Y[sample(length(Y), 300)] <- 0
  
  result_block <- fit_block_zifa(Y, k=2, n_blocks=2)
  result_normal <- fit_zifa(Y, k=2)
  
  expect_equal(result_block$Z, result_normal$Z, tolerance=1e-4)
})

test_that("handles large matrices", {
  skip_on_cran() # Skip on CRAN due to time
  set.seed(123)
  Y <- matrix(rpois(1e5, 5), nrow=500)
  Y[sample(length(Y), 3e4)] <- 0
  
  expect_silent(fit_zifa(Y, k=5, iterations=50))
})