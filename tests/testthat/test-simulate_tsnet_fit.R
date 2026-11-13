fake_fit <- function(cnames = c("a", "b", "c")) {
  structure(
    list(stan_fit = NULL, arguments = list(n_t = 40, cnames = cnames)),
    class = c("tsnet_fit", "list")
  )
}

# 10 draws: draws 1 to 6 are stationary, draws 7 to 10 are explosive
fake_draws <- function() {
  beta <- array(0, c(10, 3, 3))
  sigma <- array(0, c(10, 3, 3))
  for (s in 1:10) {
    beta[s, , ] <- if (s <= 6) diag(0.1 * s / 2, 3) else diag(1.2, 3)
    beta[s, 1, 2] <- 0.3
    sigma[s, , ] <- diag(3)
  }
  list(Beta = beta, Sigma = sigma)
}

mock_draws <- function(draws = fake_draws(), env = parent.frame()) {
  local_mocked_bindings(
    .extract_draws = function(fit, par) draws[[par]],
    .env = env
  )
}

test_that("simulate.tsnet_fit simulates one data set per draw", {
  mock_draws()
  sims <- simulate(fake_fit(), nsim = 3, seed = 1)

  expect_length(sims, 3)
  expect_equal(dim(sims[[1]]), c(40, 3))
  expect_named(sims[[1]], c("a", "b", "c"))
  expect_equal(dim(simulate(fake_fit(), n_t = 15)[[1]]), c(15, 3))

  draw_info <- attr(sims, "draws")
  expect_equal(anyDuplicated(draw_info$draw), 0)
  expect_true(all(draw_info$stationary))
  # the truth is the draw in stan orientation
  truth <- attr(sims[[1]], "truth")$Beta
  expect_equal(unname(truth), fake_draws()$Beta[draw_info$draw[1], , ])
  expect_equal(truth[["a", "b"]], 0.3)
  expect_equal(truth[["b", "a"]], 0)

  expect_equal(
    attr(sims, "n_draws"),
    c(
      total = 10, stationary = 6, nonstationary = 4, invalid = 0,
      requested = 3, returned = 3
    )
  )
})

test_that("simulate.tsnet_fit follows the seed conventions of simulate", {
  mock_draws()
  set.seed(10)
  state <- .Random.seed
  sims_a <- simulate(fake_fit(), nsim = 2, seed = 123)
  expect_identical(.Random.seed, state)
  sims_b <- simulate(fake_fit(), nsim = 2, seed = 123)

  expect_identical(sims_a, sims_b)
  expect_equal(attr(sims_a, "seed")[1], 123)

  sims_c <- simulate(fake_fit(), nsim = 2)
  expect_false(identical(.Random.seed, state))
  expect_false(identical(sims_a[[1]], sims_c[[1]]))
})

test_that("simulate.tsnet_fit handles non-stationary draws", {
  mock_draws()

  expect_warning(
    sims <- simulate(fake_fit(), nsim = 8, seed = 1),
    "Only 6 usable draws"
  )
  expect_length(sims, 6)
  expect_true(all(attr(sims, "draws")$draw <= 6))
  expect_equal(attr(sims, "n_draws")[["returned"]], 6)

  kept <- simulate(fake_fit(),
    nsim = 10, seed = 1,
    nonstationary = "keep_flagged"
  )
  expect_length(kept, 10)
  expect_equal(sum(!attr(kept, "draws")$stationary), 4)
  expect_equal(
    attr(kept, "draws")$stationary,
    attr(kept, "draws")$spectral_radius < 1
  )

  expect_error(
    simulate(fake_fit(), nsim = 10, seed = 1, nonstationary = "error"),
    "non-stationary"
  )
  expect_warning(
    simulate(fake_fit(), nsim = 11, nonstationary = "keep_flagged"),
    "Only 10 usable draws"
  )
})

test_that("simulate.tsnet_fit never uses invalid draws", {
  draws <- fake_draws()
  draws$Sigma[1, , ] <- matrix(1, 3, 3)
  draws$Beta[2, , ] <- NaN
  draws$Beta[7, 1, 1] <- NaN
  mock_draws(draws)

  expect_warning(sims <- simulate(fake_fit(), nsim = 6, seed = 1))
  expect_false(any(1:2 %in% attr(sims, "draws")$draw))
  expect_equal(
    attr(sims, "n_draws")[c("stationary", "nonstationary", "invalid")],
    c(stationary = 4, nonstationary = 3, invalid = 3)
  )
})

test_that("simulate.tsnet_fit validates its input", {
  mock_draws()

  expect_error(simulate(fake_fit(), nsim = 0), "nsim")
  expect_error(simulate(fake_fit(), n_t = -1), "n_t")
  expect_error(simulate(fake_fit(), nonstationary = "other"))
})
