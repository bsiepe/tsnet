asym_beta <- function() {
  matrix(c(
    0.3, 0.5, 0.0,
    0.0, 0.2, 0.0,
    0.0, -0.3, 0.4
  ), nrow = 3, byrow = TRUE)
}

test_that("simulate_gvar returns a data frame with the truth attached", {
  Beta <- asym_beta()
  Sigma <- diag(3)
  set.seed(1)
  sim <- simulate_gvar(Beta, Sigma, n_t = 50)

  expect_s3_class(sim, "data.frame")
  expect_equal(dim(sim), c(50, 3))
  expect_named(sim, c("V1", "V2", "V3"))
  expect_equal(
    attr(sim, "truth"),
    list(Beta = Beta, Sigma = Sigma, mu = rep(0, 3))
  )

  colnames(Beta) <- c("a", "b", "c")
  expect_named(simulate_gvar(Beta, Sigma, n_t = 5), c("a", "b", "c"))
  # names of Sigma win, e.g. for a transposed beta_mu with ".l1" suffixes
  colnames(Beta) <- c("a.l1", "b.l1", "c.l1")
  dimnames(Sigma) <- list(c("a", "b", "c"), c("a", "b", "c"))
  expect_named(simulate_gvar(Beta, Sigma, n_t = 5), c("a", "b", "c"))
})

test_that("simulate_gvar uses the session RNG", {
  set.seed(2)
  sim_a <- simulate_gvar(asym_beta(), diag(3), n_t = 20)
  set.seed(2)
  sim_b <- simulate_gvar(asym_beta(), diag(3), n_t = 20)
  sim_c <- simulate_gvar(asym_beta(), diag(3), n_t = 20)

  expect_identical(sim_a, sim_b)
  expect_false(identical(sim_a, sim_c))
})

test_that("simulate_gvar validates its input", {
  Beta <- asym_beta()
  Sigma <- diag(3)
  bad_sigma <- Sigma
  bad_sigma[1, 2] <- 0.5
  not_pd <- matrix(0.99, 3, 3)
  diag(not_pd) <- c(1, 1, -1)

  expect_error(simulate_gvar(Beta[1:2, ], Sigma, n_t = 10), "square")
  expect_error(simulate_gvar(Beta, diag(2), n_t = 10), "dimensions")
  expect_error(simulate_gvar(Beta, bad_sigma, n_t = 10), "symmetric")
  expect_error(simulate_gvar(Beta, not_pd, n_t = 10), "positive definite")
  expect_error(simulate_gvar(diag(1.01, 3), Sigma, n_t = 10), "stationary")
  expect_error(simulate_gvar(Beta, Sigma, n_t = 0), "n_t")
  expect_error(simulate_gvar(Beta, Sigma, n_t = 10, burn_in = -1), "burn_in")
  expect_error(simulate_gvar(Beta, Sigma, n_t = 10, mu = 1:2), "mu")
})

test_that("simulate_gvar uses rows of Beta as outcomes", {
  Beta <- asym_beta()
  set.seed(3)
  sim <- as.matrix(simulate_gvar(Beta, diag(3), n_t = 20000))
  n <- nrow(sim)
  # row i of the estimate holds the regression of variable i on all lags
  est <- t(stats::coef(stats::lm(sim[-1, ] ~ 0 + sim[-n, ])))

  # the standard error of each coefficient is about 0.007
  expect_lt(max(abs(unname(est) - Beta)), 0.04)
  expect_gt(max(abs(unname(est) - t(Beta))), 0.3)
})

test_that("simulate_gvar starts from the stationary distribution", {
  # close to a unit root: 100 burn-in steps from zero would only reach 87%
  # of the stationary variance
  Beta <- matrix(0.99)
  Sigma <- matrix(1)
  set.seed(6)
  first <- replicate(4000, simulate_gvar(Beta, Sigma, n_t = 1)[1, 1])

  expect_equal(var(first), 1 / (1 - 0.99^2), tolerance = 0.08)

  # cross-covariances of the first observation
  Beta <- asym_beta()
  set.seed(7)
  first <- t(replicate(
    4000, unlist(simulate_gvar(Beta, diag(3), n_t = 1, burn_in = 2))
  ))
  expect_equal(
    unname(cov(first)),
    gvar_stationary(Beta, diag(3))$cov,
    tolerance = 0.08
  )
})

test_that("simulate_gvar respects process means", {
  mu <- c(10, -5, 2)
  set.seed(4)
  sim <- simulate_gvar(asym_beta(), diag(3), n_t = 20000, mu = mu)

  expect_equal(unname(colMeans(sim)), mu, tolerance = 0.01)
  expect_equal(attr(sim, "truth")$mu, mu)
})

test_that("gvar_stationary solves the Lyapunov equation", {
  Beta <- asym_beta()
  Sigma <- diag(c(1, 2, 0.5))
  Sigma[1, 2] <- Sigma[2, 1] <- 0.4
  res <- gvar_stationary(Beta, Sigma)

  expect_equal(res$cov, Beta %*% res$cov %*% t(Beta) + Sigma)

  sds <- sqrt(diag(res$cov))
  expect_equal(res$Beta_std, diag(1 / sds) %*% Beta %*% diag(sds))
  expect_equal(res$Sigma_std, diag(1 / sds) %*% Sigma %*% diag(1 / sds))
  # the standardized process has unit variances
  expect_equal(
    diag(gvar_stationary(res$Beta_std, res$Sigma_std)$cov),
    rep(1, 3)
  )

  Theta <- solve(Sigma)
  expect_equal(
    res$pcor[1, 2],
    -Theta[1, 2] / sqrt(Theta[1, 1] * Theta[2, 2])
  )
  expect_equal(diag(res$pcor), rep(0, 3))
  expect_error(gvar_stationary(diag(1, 3), Sigma), "stationary")
})

test_that("gvar_stationary matches the covariance of simulated data", {
  Beta <- asym_beta()
  Sigma <- diag(c(1, 2, 0.5))
  Sigma[1, 2] <- Sigma[2, 1] <- 0.4
  set.seed(5)
  sim <- simulate_gvar(Beta, Sigma, n_t = 100000)

  expect_equal(
    unname(stats::cov(sim)),
    gvar_stationary(Beta, Sigma)$cov,
    tolerance = 0.03
  )
})
