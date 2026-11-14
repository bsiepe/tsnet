# Checks that prior-only fits of stan_gvar() sample from the documented priors.
# Run from the package root: Rscript dev/prior_checks.R
# Results are written to dev/output/prior_checks.txt and prior_checks.rds.
#
# Quantiles of the Stan draws are compared with reference quantiles through
# z = (q_stan - q_ref) / sqrt(mcse_stan^2 + se_ref^2), where mcse_stan comes
# from posterior::mcse_quantile() and se_ref is zero for closed-form
# references and based on the Kish effective sample size otherwise.

devtools::load_all(quiet = TRUE)
source("dev/prior_helpers.R")
dir.create("dev/output", showWarnings = FALSE)
set.seed(2026)

K <- 3
prs <- pairs_upper(K)
n_ref <- 4e5
results <- list()
diag_tab <- list()

add <- function(case, tab) {
  tab$case <- case
  results[[length(results) + 1]] <<- tab
}

## 0. share of stationary draws under the default Beta prior ------------------
stationary_share <- vapply(c(3, 4, 6, 8), function(k) {
  mean(replicate(20000, {
    max(Mod(eigen(matrix(rnorm(k * k, 0, 0.5), k), only.values = TRUE)$values))
  }) < 1)
}, numeric(1))
names(stationary_share) <- paste0("K=", c(3, 4, 6, 8))

## 1. validate the R sampler for LKJ(eta) --------------------------------------
lkj_validation <- do.call(rbind, lapply(
  list(c(3, 1), c(3, 2), c(6, 1), c(6, 3)),
  function(par) {
    k <- par[1]
    eta <- par[2]
    S <- rlkj(2e5, k, eta)
    shape <- eta - 1 + k / 2
    idx <- pairs_upper(k)
    do.call(rbind, lapply(seq_len(nrow(idx)), function(r) {
      x <- S[, idx[r, 1], idx[r, 2]]
      q_ref <- 2 * qbeta(probs, shape, shape) - 1
      se <- sqrt(probs * (1 - probs) / length(x)) /
        (dbeta((q_ref + 1) / 2, shape, shape) / 2)
      data.frame(
        label = sprintf("K=%d eta=%d r[%d,%d]", k, eta, idx[r, 1], idx[r, 2]),
        prob = probs,
        sampler = unname(quantile(x, probs)),
        exact = q_ref,
        z = (unname(quantile(x, probs)) - q_ref) / se
      )
    }))
  }
))

## 2. Beta: Normal(loc, scale) elementwise -------------------------------------
beta_loc <- matrix(c(
  0.4, -0.3, 0.1,
  0.0, 0.2, 0.5,
  -0.6, 0.0, -0.1
), nrow = K, byrow = TRUE)
beta_scale <- matrix(c(
  0.1, 0.2, 0.3,
  0.5, 0.7, 1.0,
  0.05, 1.5, 0.4
), nrow = K, byrow = TRUE)

for (cov_prior in c("IW", "LKJ")) {
  # flat beta factor on the partial correlations (ignored by the IW model):
  # prior-only chains of the LKJ model get stuck with the default, see (c)
  res <- fit_prior(
    list(
      prior_Beta_loc = beta_loc, prior_Beta_scale = beta_scale,
      prior_Rho_scale = matrix(2, K, K)
    ),
    cov_prior
  )
  arr <- stan_array(res, "Beta")
  for (i in 1:K) {
    for (j in 1:K) {
      add(paste("Beta", cov_prior), compare_quantiles(
        arr[, , sprintf("Beta[%d,%d]", i, j)],
        list(q = qnorm(probs, beta_loc[i, j], beta_scale[i, j])),
        sprintf("Beta[%d,%d]", i, j)
      ))
    }
  }
  diag_tab[[length(diag_tab) + 1]] <- diagnostics(res, paste("Beta", cov_prior))
}

## 3. LKJ models ----------------------------------------------------------------
compare_rho <- function(res, case, ref_fun) {
  arr <- stan_array(res, "Rho")
  for (r in seq_len(nrow(prs))) {
    i <- prs[r, 1]
    j <- prs[r, 2]
    add(case, compare_quantiles(
      arr[, , sprintf("Rho[%d,%d]", i, j)],
      ref_fun(i, j),
      sprintf("Rho[%d,%d]", i, j)
    ))
  }
  diag_tab[[length(diag_tab) + 1]] <<- diagnostics(res, case)
}

corr_from_chol <- function(L) {
  apply(L, 1, tcrossprod, simplify = FALSE) |>
    simplify2array() |>
    aperm(c(3, 1, 2))
}

# importance weights of LKJ(eta) draws for the joint prior with beta factors
lkj_reference <- function(eta, loc, kappa) {
  rho <- pcor_draws(rlkj(n_ref, K, eta))
  log_w <- rep(0, n_ref)
  for (r in seq_len(nrow(prs))) {
    i <- prs[r, 1]
    j <- prs[r, 2]
    log_w <- log_w + dbeta(
      rho[, i, j] / 2 + 0.5,
      loc[i, j] * kappa[i, j], (1 - loc[i, j]) * kappa[i, j],
      log = TRUE
    )
  }
  w <- exp(log_w - max(log_w))
  list(
    rho = rho, w = w, ess = kish_ess(w),
    pareto_k = suppressWarnings(loo::psis(log_w, r_eff = NA))$diagnostics$pareto_k
  )
}
is_info <- list()

# (a) flat beta factor: (rho + 1) / 2 ~ Beta(eta, eta)
for (eta in c(1, 2)) {
  res <- fit_prior(
    list(
      prior_Rho_loc = matrix(0.5, K, K),
      prior_Rho_scale = matrix(2, K, K),
      prior_Eta = eta
    ),
    "LKJ"
  )
  compare_rho(
    res, sprintf("LKJ (a) kappa=2 eta=%d", eta),
    \(i, j) list(q = 2 * qbeta(probs, eta, eta) - 1)
  )
  if (eta == 1) {
    # Rho recomputed in R from the L_Theta draws
    draws_flat <- rstan::extract(res$fit$stan_fit, c("L_Theta", "Rho"))
    rho_recompute_flat <- max(abs(
      pcor_draws(corr_from_chol(draws_flat$L_Theta)) - draws_flat$Rho
    ))
  }
}

# (b) asymmetric locations, kappa = 10, importance-sampling reference
rho_loc <- matrix(0.5, K, K)
rho_loc[1, 2] <- rho_loc[2, 1] <- 0.7
rho_loc[2, 3] <- rho_loc[3, 2] <- 0.35
kappa10 <- matrix(10, K, K)
for (eta in c(1, 2)) {
  case <- sprintf("LKJ (b) kappa=10 eta=%d", eta)
  res <- fit_prior(
    list(prior_Rho_loc = rho_loc, prior_Rho_scale = kappa10, prior_Eta = eta),
    "LKJ"
  )
  ref <- lkj_reference(eta, rho_loc, kappa10)
  is_info[[case]] <- c(kish_ess = ref$ess, pareto_k = ref$pareto_k)
  compare_rho(res, case, \(i, j) list(x = ref$rho[, i, j], w = ref$w))

  # the same draws against a reference with mirrored locations: the check has
  # to fail clearly if the sign of the partial correlations were wrong
  ref_flip <- lkj_reference(eta, 1 - rho_loc, kappa10)
  arr <- stan_array(res, "Rho")
  for (r in seq_len(nrow(prs))) {
    i <- prs[r, 1]
    j <- prs[r, 2]
    add(paste(case, "SIGN-FLIPPED reference"), compare_quantiles(
      arr[, , sprintf("Rho[%d,%d]", i, j)],
      list(x = ref_flip$rho[, i, j], w = ref_flip$w),
      sprintf("Rho[%d,%d]", i, j)
    ))
  }
}

# (c) package defaults: smoke test, no pass/fail
res_default <- fit_prior(NULL, "LKJ")
ref_default <- lkj_reference(
  1, matrix(0.5, K, K), matrix(sqrt(0.5), K, K)
)
is_info[["LKJ (c) defaults"]] <- c(
  kish_ess = ref_default$ess, pareto_k = ref_default$pareto_k
)
compare_rho(
  res_default, "LKJ (c) defaults (smoke test)",
  \(i, j) list(x = ref_default$rho[, i, j], w = ref_default$w)
)
draws_default <- rstan::extract(res_default$fit$stan_fit, c("L_Theta", "Rho"))
rho_in_r <- pcor_draws(corr_from_chol(draws_default$L_Theta))
rho_recompute_maxdiff <- suppressWarnings(
  max(abs(rho_in_r - draws_default$Rho), na.rm = TRUE)
)
default_summary <- rstan::summary(
  res_default$fit$stan_fit,
  pars = c("Rho[1,2]", "Rho[1,3]", "Rho[2,3]", "sigma_theta")
)$summary

## 4. Wishart models ------------------------------------------------------------
wishart_reference <- function(delta, S) {
  nu <- delta + K - 1
  Theta <- array(
    apply(rWishart(n_ref, nu, solve(S)), 3, solve),
    c(K, K, n_ref)
  ) |> aperm(c(3, 1, 2))
  list(rho = pcor_draws(Theta, is_precision = TRUE), Theta = Theta)
}

S_nonid <- matrix(c(
  1.0, 0.5, 0.2,
  0.5, 2.0, -0.3,
  0.2, -0.3, 0.5
), nrow = K)

wishart_cases <- list(
  list(marginal = 0.25, S = diag(K), name = "S=I marginal=0.25 (default)"),
  list(marginal = 0.1, S = diag(K), name = "S=I marginal=0.1"),
  list(marginal = 0.25, S = S_nonid, name = "S non-identity marginal=0.25")
)
wishart_moments <- list()
for (wc in wishart_cases) {
  delta <- 1 / wc$marginal - 1
  res <- fit_prior(
    list(prior_Rho_marginal = wc$marginal, prior_S = wc$S), "IW"
  )
  ref <- wishart_reference(delta, wc$S)
  compare_rho(
    res, paste("IW", wc$name, "vs inverse-Wishart precision draws"),
    \(i, j) list(x = ref$rho[, i, j])
  )
  if (identical(wc$S, diag(K))) {
    arr <- stan_array(res, "Rho")
    for (r in seq_len(nrow(prs))) {
      nm <- sprintf("Rho[%d,%d]", prs[r, 1], prs[r, 2])
      add(
        paste("IW", wc$name, "vs Beta(delta/2, delta/2)"),
        compare_quantiles(
          arr[, , nm],
          list(q = 2 * qbeta(probs, delta / 2, delta / 2) - 1), nm
        )
      )
    }
  }
  # precision matrix itself, on the log scale for the diagonal
  theta_stan <- rstan::extract(res$fit$stan_fit, "Theta", permuted = FALSE)
  add(
    paste("IW", wc$name, "log Theta[1,1]"),
    compare_quantiles(
      log(theta_stan[, , "Theta[1,1]"]),
      list(x = log(ref$Theta[, 1, 1])), "log Theta[1,1]"
    )
  )
  rho_stan <- rstan::extract(res$fit$stan_fit, "Rho")$Rho
  wishart_moments[[wc$name]] <- c(
    prior_Rho_marginal = wc$marginal,
    delta = delta,
    sd_rho_stan = mean(
      apply(rho_stan, c(2, 3), sd, na.rm = TRUE)[upper.tri(diag(K))]
    ),
    sd_if_marginal_were_sd = wc$marginal,
    sd_beta_formula = sqrt(1 / (delta + 1))
  )
}

## 5. report ---------------------------------------------------------------------
all_res <- do.call(rbind, results)
all_diag <- do.call(rbind, diag_tab)
real <- !grepl("SIGN-FLIPPED|smoke test", all_res$case)
by_case <- do.call(rbind, lapply(split(all_res, all_res$case), function(d) {
  data.frame(
    case = d$case[1],
    comparisons = nrow(d),
    n_abs_z_gt_3 = sum(abs(d$z) > 3, na.rm = TRUE),
    max_abs_z = max(abs(d$z), na.rm = TRUE),
    rms_z = sqrt(mean(d$z^2, na.rm = TRUE))
  )
}))
rownames(by_case) <- NULL

saveRDS(
  list(
    results = all_res, diagnostics = all_diag, by_case = by_case,
    is_info = is_info, lkj_validation = lkj_validation,
    stationary_share = stationary_share, wishart_moments = wishart_moments,
    rho_recompute_maxdiff = rho_recompute_maxdiff,
    rho_recompute_flat = rho_recompute_flat,
    default_summary = default_summary
  ),
  "dev/output/prior_checks.rds"
)

options(width = 150, digits = 4)
sink("dev/output/prior_checks.txt")
cat("tsnet prior checks,", format(Sys.time()), "\n")
cat("K = 3, 4 chains x 5000 draws per prior-only fit,", n_ref, "reference draws\n")
cat("\n== Share of stationary draws, Beta_ij ~ N(0, 0.5) ==\n")
print(stationary_share)
cat("\n== LKJ sampler in R vs Beta(eta - 1 + K/2) marginals ==\n")
cat(
  "comparisons:", nrow(lkj_validation),
  " |z| > 3:", sum(abs(lkj_validation$z) > 3),
  " max |z|:", max(abs(lkj_validation$z)), "\n"
)
cat("\n== Summary per case ==\n")
print(by_case, row.names = FALSE)
cat(
  "\nAll pass/fail comparisons:", sum(real), " |z| > 3:",
  sum(abs(all_res$z[real]) > 3, na.rm = TRUE), "\n"
)
cat("\n== Importance-sampling references (Kish ESS, Pareto k) ==\n")
print(do.call(rbind, is_info))
cat("\n== Sampler diagnostics of the prior-only fits ==\n")
print(all_diag, row.names = FALSE)
cat("\n== Wishart models: SD of partial correlations ==\n")
print(do.call(rbind, wishart_moments))
cat("\n== Rho recomputed in R from L_Theta ==\n")
cat(
  "kappa = 2, eta = 1: max |Rho from L_Theta in R - Rho from Stan|:",
  rho_recompute_flat, "\n"
)
cat("\n== LKJ defaults ==\n")
cat("max |Rho from L_Theta in R - Rho from Stan|:", rho_recompute_maxdiff, "\n")
cat(
  "draws with a singular correlation matrix L_Theta L_Theta':",
  sum(is.na(rho_in_r[, 1, 2])), "of", dim(rho_in_r)[1], "\n"
)
print(default_summary)
cat("\nwarnings of the default LKJ fit:\n")
print(unique(res_default$warnings))
cat("\n== All quantile comparisons ==\n")
print(all_res[, c("case", "label", "prob", "stan", "ref", "mcse_stan", "se_ref", "z")],
  row.names = FALSE
)
cat("\n== LKJ sampler validation, all rows ==\n")
print(lkj_validation, row.names = FALSE)
sink()
