test_that("ivgrf works", {
  n <- 2000
  p <- 5
  X <- matrix(rbinom(n * p, 1, 0.5), n, p)
  Z <- rbinom(n, 1, 0.5)
  Q <- rbinom(n, 1, 0.5)
  W <- Q * Z
  tau <-  X[, 1] / 2
  Y <- rowSums(X[, 1:3]) + tau * W + Q + rnorm(n)

  fit <- causal_iv_forest(X, Y, W, Z)
  pp <- predict(fit)$predictions
  scores <- grf::get_scores(fit)

  ate <- grf::average_treatment_effect(fit)
  rate <- grf::rank_average_treatment_effect(fit, pp)


  expect_true(TRUE)
})
