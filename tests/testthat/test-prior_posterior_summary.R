summary_fit <- function(prior_only = FALSE, priors = gvar_priors(3)) {
  structure(
    list(
      stan_fit = NULL,
      arguments = list(
        p = 3, n_t = 40, cnames = c("a", "b", "c"),
        priors = priors, prior_only = prior_only
      )
    ),
    class = c("tsnet_fit", "list")
  )
}

# posterior draws are narrow, prior draws are wide
summary_draws <- function(fit, par) {
  set.seed(1)
  n <- 4000
  sd_draws <- if (isTRUE(fit$arguments$prior_only)) 0.5 else 0.1
  draws <- array(rnorm(n * 9, sd = sd_draws), c(n, 3, 3))
  if (par == "Beta") {
    # effect of lagged b on a
    draws[, 1, 2] <- draws[, 1, 2] + 0.4
  }
  if (par == "Rho") {
    draws[, 1, 3] <- draws[, 3, 1] <- draws[, 1, 3] - 0.3
  }
  draws
}

test_that("prior_posterior_summary summarizes all edges", {
  local_mocked_bindings(
    .extract_draws = summary_draws,
    .cov_model = function(fit) "VAR_wishart"
  )
  res <- prior_posterior_summary(summary_fit(), summary_fit(prior_only = TRUE))

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 9 + 3)
  expect_equal(sum(res$type == "contemporaneous"), 3)

  # orientation: from = lagged predictor, to = outcome
  edge <- res[res$type == "temporal" & res$from == "b" & res$to == "a", ]
  expect_equal(edge$post_mean, 0.4, tolerance = 0.05)
  reverse <- res[res$type == "temporal" & res$from == "a" & res$to == "b", ]
  expect_equal(reverse$post_mean, 0, tolerance = 0.02)

  temporal <- res[res$type == "temporal", ]
  expect_equal(temporal$prior_mean, rep(0, 9))
  expect_equal(temporal$prior_sd, rep(0.5, 9))
  expect_equal(temporal$contraction, rep(1 - 0.1^2 / 0.5^2, 9), tolerance = 0.05)

  pcor <- res[res$type == "contemporaneous", ]
  expect_equal(pcor$from, c("a", "a", "b"))
  expect_equal(pcor$to, c("b", "c", "c"))
  expect_equal(pcor$post_mean, c(0, -0.3, 0), tolerance = 0.02)
  expect_equal(pcor$prior_sd, rep(0.5, 3), tolerance = 0.05)
  expect_equal(pcor$contraction, rep(0.96, 3), tolerance = 0.05)
})

test_that("prior_posterior_summary works without a prior-only fit", {
  local_mocked_bindings(.extract_draws = summary_draws)
  res <- prior_posterior_summary(summary_fit())
  pcor <- res[res$type == "contemporaneous", ]

  expect_true(all(is.na(pcor$prior_sd)))
  expect_true(all(is.na(pcor$contraction)))
  expect_false(anyNA(res[res$type == "temporal", ]))
})

test_that("prior_posterior_summary works for a single variable", {
  fit <- summary_fit(priors = gvar_priors(1))
  fit$arguments$p <- 1
  fit$arguments$cnames <- "a"
  local_mocked_bindings(
    .extract_draws = function(fit, par) array(rnorm(100, sd = 0.1), c(100, 1, 1))
  )
  res <- prior_posterior_summary(fit)

  expect_equal(nrow(res), 1)
  expect_equal(res$type, "temporal")
})

test_that("prior_posterior_summary validates its input", {
  local_mocked_bindings(
    .extract_draws = summary_draws,
    .cov_model = function(fit) "VAR_wishart"
  )
  no_priors <- summary_fit()
  no_priors$arguments$priors <- NULL
  other_priors <- summary_fit(
    prior_only = TRUE,
    priors = gvar_priors(3, prior_Rho_marginal = 0.1)
  )

  expect_error(prior_posterior_summary(list()), "tsnet_fit")
  expect_error(
    prior_posterior_summary(summary_fit(prior_only = TRUE)),
    "prior-only"
  )
  expect_error(prior_posterior_summary(no_priors), "return_all")
  expect_error(
    prior_posterior_summary(summary_fit(), summary_fit()),
    "prior_only = TRUE"
  )
  expect_warning(
    prior_posterior_summary(summary_fit(), other_priors),
    "same priors"
  )
})
