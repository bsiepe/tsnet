# Evidence for dev/prior_spec_check.md: what the priors on the contemporaneous
# network of stan_gvar() imply, based on prior-only fits.
# Run from the package root: Rscript dev/prior_spec_check.R
# Writes dev/output/prior_spec_check.txt and two figures.

devtools::load_all(quiet = TRUE)
source("dev/prior_helpers.R")
dir.create("dev/output", showWarnings = FALSE)
set.seed(2027)

rho_pooled <- function(res) {
  rho <- rstan::extract(res$fit$stan_fit, "Rho")$Rho
  K <- dim(rho)[2]
  x <- as.vector(apply(rho, 1, \(m) m[upper.tri(diag(K))]))
  x[is.finite(x)]
}

describe <- function(x) {
  c(
    mean = mean(x), sd = sd(x),
    quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)),
    "P(|rho|>0.5)" = mean(abs(x) > 0.5),
    "P(|rho|>0.9)" = mean(abs(x) > 0.9)
  )
}

## fact 6: LKJ model, beta_proportion(mu, kappa) --------------------------------
lkj_settings <- expand.grid(K = c(3, 6), kappa = c(sqrt(0.5), 2))
lkj_fits <- lapply(seq_len(nrow(lkj_settings)), function(r) {
  K <- lkj_settings$K[r]
  fit_prior(
    list(prior_Rho_scale = matrix(lkj_settings$kappa[r], K, K)),
    "LKJ",
    K = K, seed = r,
    # stuck chains run at the maximum tree depth and are very slow
    iter_sampling = 2000, iter_warmup = 1000
  )
})
lkj_rho <- lapply(lkj_fits, rho_pooled)
lkj_tab <- cbind(lkj_settings, do.call(rbind, lapply(lkj_rho, describe)))
lkj_diag <- do.call(rbind, lapply(seq_along(lkj_fits), function(r) {
  diagnostics(lkj_fits[[r]], sprintf(
    "LKJ K=%d kappa=%.3f", lkj_settings$K[r], lkj_settings$kappa[r]
  ))
}))

png("dev/output/prior_spec_lkj.png", width = 1600, height = 1400, res = 200)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (r in seq_len(nrow(lkj_settings))) {
  kappa <- lkj_settings$kappa[r]
  hist(lkj_rho[[r]],
    breaks = seq(-1, 1, by = 0.05), freq = FALSE,
    col = "grey80", border = "white", xlab = "partial correlation",
    main = sprintf(
      "K = %d, kappa = %s%s", lkj_settings$K[r], format(kappa, digits = 3),
      if (kappa < 1) " (default up to 0.2.0)" else " (new default)"
    )
  )
  # beta factor on its own, and the uniform density
  curve(dbeta((x + 1) / 2, 0.5 * kappa, 0.5 * kappa) / 2,
    from = -0.995, to = 0.995, add = TRUE, lty = 2
  )
  abline(h = 0.5, lty = 3)
}
dev.off()

## fact 7: Wishart model ---------------------------------------------------------
marginals <- c(0.1, 0.25, 0.5)
iw_fits <- lapply(marginals, function(m) {
  fit_prior(list(prior_Rho_marginal = m), "IW", seed = 10)
})
iw_rho <- lapply(iw_fits, rho_pooled)
iw_tab <- data.frame(
  prior_Rho_marginal = marginals,
  delta = 1 / marginals - 1,
  var_rho_draws = vapply(iw_rho, var, numeric(1)),
  sd_rho_draws = vapply(iw_rho, sd, numeric(1)),
  sd_rho_analytic = sqrt(marginals)
)
iw_diag <- do.call(rbind, lapply(seq_along(iw_fits), function(r) {
  diagnostics(iw_fits[[r]], sprintf("IW marginal=%.2f", marginals[r]))
}))

png("dev/output/prior_spec_wishart.png", width = 1600, height = 800, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(iw_rho[[2]],
  breaks = seq(-1, 1, by = 0.05), freq = FALSE,
  col = "grey80", border = "white", xlab = "partial correlation",
  main = "prior_Rho_marginal = 0.25 (default)", ylim = c(0, 1.7)
)
# implied: delta = 3, Beta(1.5, 1.5); documented: SD = 0.25, Beta(7.5, 7.5)
curve(dbeta((x + 1) / 2, 1.5, 1.5) / 2, add = TRUE, lwd = 2)
curve(dbeta((x + 1) / 2, 7.5, 7.5) / 2, add = TRUE, lty = 2)
legend("topleft",
  lty = c(1, 2), lwd = c(2, 1), bty = "n", cex = 0.7,
  legend = c("Beta(1.5, 1.5): SD 0.5", "Beta(7.5, 7.5): SD 0.25")
)
plot(marginals, iw_tab$sd_rho_draws,
  pch = 19, xlim = c(0, 0.55), ylim = c(0, 0.8),
  xlab = "prior_Rho_marginal", ylab = "SD of partial correlations",
  main = "S = I"
)
curve(sqrt(x), add = TRUE, lwd = 2)
abline(0, 1, lty = 2)
legend("topleft",
  lty = c(1, 2), lwd = c(2, 1), bty = "n", cex = 0.7,
  legend = c("sqrt(value) (code)", "value (documentation)")
)
dev.off()

## report -------------------------------------------------------------------------
options(width = 150, digits = 3)
sink("dev/output/prior_spec_check.txt")
cat("tsnet prior specification check,", format(Sys.time()), "\n")
cat("\n== Fact 6: LKJ model, pooled partial correlations ==\n")
cat("beta factor: Beta(0.5 * kappa, 0.5 * kappa) on (rho + 1) / 2\n")
print(lkj_tab, row.names = FALSE)
print(lkj_diag, row.names = FALSE)
cat("\n== Fact 7: Wishart model, S = I, K = 3 ==\n")
print(iw_tab, row.names = FALSE)
print(iw_diag, row.names = FALSE)
cat("\n== Fact 8: default of cov_prior ==\n")
print(formals(stan_gvar)$cov_prior)
sink()
