context("Numerical Robustness Tests")

test_that("handles rank-deficient matrices", {
  set.seed(123)
  Y &lt;- matrix(rnorm(100), ncol=10)
  Y[,5:10] &lt;- Y[,1:4] %*% matrix(rnorm(24), ncol=6) # Create rank deficiency
  Y[sample(length(Y), 30)] &lt;- 0 # Add zeros
  expect_silent(fit_zifa(Y, k=4))
})

test_that("converges with different thresholds", {
  set.seed(123)
  Y &lt;- matrix(rpois(1000, 5), nrow=100)
  Y[sample(length(Y), 300)] &lt;- 0
  
  # Test loose convergence
  result_loose &lt;- fit_zifa(Y, k=2, convergence_criterion=1e-3)
  # Test tight convergence 
  result_tight &lt;- fit_zifa(Y, k=2, convergence_criterion=1e-6)
  
  expect_true(length(result_loose$ll) &lt; length(result_tight$ll))
})

test_that("monotonically increasing log-likelihood", {
  set.seed(123)
  Y &lt;- matrix(rpois(500, 5), nrow=50)
  Y[sample(length(Y), 150)] &lt;- 0
  result &lt;- fit_zifa(Y, k=2)
  
  # Check log-likelihood increases
  ll_diff &lt;- diff(result$ll)
  expect_true(all(ll_diff &gt;= 0 | abs(ll_diff) &lt; 1e-6)) # Allow small numerical decreases
})