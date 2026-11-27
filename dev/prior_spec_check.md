# Prior specification check

Three places where the documentation of `stan_gvar()` in tsnet 0.2.0 and the
code disagree about the priors on the contemporaneous network. Line numbers of
the documentation refer to `R/stan_gvar.R` at `origin/main` (78be303), those
of the code to this branch.

Evidence was produced by `dev/prior_spec_check.R` (figures and
`dev/output/prior_spec_check.txt`) and `dev/prior_checks.R`
(`dev/output/prior_checks.txt`), both of which use prior-only fits
(`prior_only = TRUE`). Only the documentation was changed on this branch.
Whether the code or the old documentation reflects the intended prior is a
decision for the maintainer; the options are listed per issue.

**Update after review of this report.** The maintainer decided to change the
default of `prior_Rho_scale` to 2 (issue 1, second option). "Default" in
section 1, in the figures and in `dev/output/*.txt` refers to the default up
to version 0.2.0, `sqrt(0.5)`; these outputs were produced before the change.
Issues 2 and 3 remain documentation-only.

## 1. LKJ model: `prior_Rho_scale` is a precision, the default is U-shaped

**Documentation (0.2.0, l. 128-134).** "a proportional (mean-variance
parameterized) beta prior", and "The variance parameter of sqrt(0.5) implies a
uniform distribution of partial correlations."

**Code.** `inst/stan/VAR_LKJ.stan:66` (and `VAR_LKJ_beep.stan:75`):

```stan
target += beta_proportion_lpdf(Rho[i,j] / 2 + 0.5 | prior_Rho_loc[i,j], prior_Rho_scale[i,j]);
```

for `i < j`, in addition to `lkj_corr_cholesky_lpdf(L_Theta | prior_Eta)`.
The default is `prior_Rho_scale = sqrt(0.5)` (`R/gvar_priors.R:57`).

**Implied prior.** Stan's `beta_proportion(mu, kappa)` is
Beta(mu * kappa, (1 - mu) * kappa); kappa is a precision. The joint prior
density of the correlation matrix R is

    p(R) ∝ LKJ(R | eta) * prod_{i<j} Beta((rho_ij + 1) / 2 | mu_ij kappa_ij, (1 - mu_ij) kappa_ij)

* With mu = 0.5 and kappa = 2 the beta factor is Beta(1, 1), i.e., flat. The
  prior is then LKJ(eta), under which each partial correlation given all
  other variables has the marginal (rho + 1) / 2 ~ Beta(eta, eta) for any K.
  This is uniform for eta = 1. So the value that gives the documented
  "uniform distribution" is 2, not sqrt(0.5).
* With the default kappa = sqrt(0.5) the factor is Beta(0.354, 0.354):
  U-shaped, with poles at partial correlations of -1 and 1.
* The default is not only U-shaped; for K >= 3 the joint prior appears to be
  improper. Since 1 - rho_ij^2 = det(R) det(R_{-ij}) / (det(R_{-i}) det(R_{-j})),
  every beta factor behaves like det(R)^(a - 1) with a = 0.354 near the
  boundary of the set of correlation matrices, so that
  p(R) ~ det(R)^(eta - 1 - 0.646 K (K - 1) / 2). For K = 3 and eta = 1 the
  exponent is -1.94, and det(R)^c is only integrable at 0 for c > -1. (K = 2
  is proper, exponent -0.646.) This derivation is mine and worth a second
  look. With data the posterior is proper, because the likelihood vanishes at
  the boundary, but the prior pulls towards strong partial correlations and
  nearly singular innovation covariance matrices, more so for larger K and
  shorter time series.

**Evidence from prior-only draws** (4 chains, pooled over all pairs;
`dev/output/prior_spec_lkj.png`):

| K | kappa | SD(rho) | P(abs(rho) > 0.5) | P(abs(rho) > 0.9) | divergent | at max. treedepth | max. Rhat |
|---|-------|---------|-------------------|-------------------|-----------|-------------------|-----------|
| 3 | sqrt(0.5) (default) | 0.95 | 1.00 | 1.00 | 2491 / 8000 | 5509 / 8000 | not defined (156 NaN draws) |
| 6 | sqrt(0.5) (default) | 0.99 | 1.00 | 1.00 | 1308 / 8000 | 6692 / 8000 | 7.6e14 |
| 3 | 2 | 0.576 | 0.497 | 0.100 | 0 | 0 | 1.00 |
| 6 | 2 | 0.574 | 0.493 | 0.097 | 0 | 0 | 1.00 |

With the default, every chain runs into the boundary and stays there: all
draws of the partial correlations are numerically -1 or 1, and the
correlation matrix `L_Theta L_Theta'` is singular in 20000 of 20000 draws of
the K = 3 fit in `dev/prior_checks.R`. These draws are not a sample from
anything and only show that the mass sits at the boundary. An
importance-sampling reference for the default prior fails for the same
reason (Kish ESS 1.4 of 4e5 draws, Pareto k 2.0).

With kappa = 2 the draws match the closed form: for a uniform distribution
SD = 0.577, P(abs(rho) > 0.5) = 0.5 and P(abs(rho) > 0.9) = 0.1, and all 42
quantile comparisons for eta = 1 and eta = 2 in `dev/prior_checks.R` are
within abs(z) < 1.9. Informative factors (kappa = 10, locations 0.7, 0.5,
0.35) also match their importance-sampling reference (42 comparisons,
max abs(z) = 2.9, Pareto k < 0), so the mechanism itself, including the
sign of the partial correlations, works as intended.

**Options.**

* *Documentation only (done on this branch, commit e460fdc).* Results of
  existing analyses do not change. Users of `cov_prior = "LKJ"` with default
  priors keep a prior that favors strong partial correlations; prior-only
  fits of this model are unusable (now with a warning, commit cda5389).
* *Change the default to 2* (and/or rename the argument, since it is a
  precision). This gives the documented uniform prior. Results for
  `cov_prior = "LKJ"` with default priors change: partial correlations are
  shrunk less towards -1 and 1, i.e., estimates move towards zero compared
  with 0.2.0, mostly for short time series. Users who set `prior_Rho_scale`
  themselves, believing it to be a variance-like quantity, are affected
  under both options and should be pointed to the parameterization.

## 2. Wishart model: the prior is on the precision matrix, and `prior_Rho_marginal` is a variance

**Documentation (0.2.0, l. 140-153).** "The Inverse-Wishart prior is a
distribution on the innovation covariance matrix Sigma: Sigma ~ IW(nu, S)",
with delta = 1 / s_rho - 1, where s_rho "is the standard deviation of the
implied marginal beta distribution of the partial correlations", and "The
user can set a prior on the expected standard deviation of the partial
correlations by specifying a `prior_Rho_marginal` parameter. The default
value is 0.25". The code comment (0.2.0, l. 243) says
"SD = 1/(delta+1), delta = (1 / SD) - 1".

**Code.** `inst/stan/VAR_wishart.stan:56` (and `VAR_wishart_beep.stan:58`):

```stan
cov_matrix[K] Theta;                       // parameters
matrix[K,K] Sigma = inverse_spd(Theta);    // transformed parameters
target += inv_wishart_lpdf(Theta | prior_delta + K - 1, prior_S);
```

and `R/stan_gvar.R:280`: `prior_delta <- (1 / priors$prior_Rho_marginal) - 1`.

**Implied prior.** The inverse Wishart is on the precision matrix
Theta = Sigma^-1, which is equivalent to Sigma ~ Wishart(nu, S^-1) with
nu = delta + K - 1. For Theta ~ IW(delta + K - 1, I), the correlations of
Theta marginally follow (r + 1) / 2 ~ Beta(delta / 2, delta / 2), and the
partial correlations are rho_ij = -r_ij, so they have the same symmetric
distribution, with Var(rho) = 1 / (delta + 1). Inserting
delta = 1 / m - 1 gives Var(rho) = m. `prior_Rho_marginal` is therefore the
prior variance of the partial correlations; the SD is its square root. The
default 0.25 means delta = 3, a Beta(1.5, 1.5) shape and SD 0.5, not 0.25.
An SD of 0.25 would need `prior_Rho_marginal = 0.0625` (delta = 15). The code
comment has the same error. The parameterization as such, a prior on the
partial correlations with marginal Beta(delta / 2, delta / 2), is the usual
matrix-F / BGGM-style setup; only the description of the argument is off.

**Evidence from prior-only draws** (K = 3, S = I, 20000 draws;
`dev/output/prior_spec_wishart.png`):

| `prior_Rho_marginal` | delta | Var(rho), draws | SD(rho), draws | sqrt(value) |
|----------------------|-------|-----------------|----------------|-------------|
| 0.10 | 9 | 0.099 | 0.315 | 0.316 |
| 0.25 (default) | 3 | 0.253 | 0.503 | 0.500 |
| 0.50 | 1 | 0.511 | 0.715 | 0.707 |

In `dev/prior_checks.R`, the partial correlations match both
Beta(delta / 2, delta / 2) and partial correlations computed from
`solve(rWishart(1, nu, solve(prior_S))[, , 1])` used as precision matrix,
including one non-identity `prior_S` (105 quantile comparisons, none with
abs(z) > 3, max 2.2). `log(Theta[1, 1])` matches the same reference (21
comparisons, max abs(z) = 1.75), which confirms that the inverse Wishart is
on the precision matrix.

Two side observations: `prior_delta` is declared `real<lower=1>`, so
`prior_Rho_marginal` has to be at most 0.5 (now documented and checked by
`gvar_priors()`); and at that limit (delta = 1, a U-shaped Beta(0.5, 0.5))
the prior-only fit had 2418 of 20000 divergent transitions.

**Options.**

* *Documentation only (done on this branch, commit 6a0d75b).* No results
  change. The default remains a fairly wide prior (SD 0.5) that the old
  documentation described as SD 0.25. If the simulation study behind "has
  worked well" used the code, its conclusions hold for SD 0.5.
* *Change the code so that the argument is an SD*
  (`prior_delta <- 1 / prior_Rho_marginal^2 - 1`). With the default 0.25 this
  gives delta = 15 instead of 3, a much more informative prior that shrinks
  partial correlations towards zero more strongly. All results with
  `cov_prior = "IW"`, which is the default of `stan_gvar()`, would change,
  noticeably so for short time series, and would no longer correspond to the
  published simulation results. A middle way that keeps results unchanged
  would be to interpret the argument as an SD and change the default to 0.5.

## 3. Default of `cov_prior`

**Documentation (0.2.0, l. 38).** "Default: `"LKJ"`".

**Code.** `R/stan_gvar.R:218`: `cov_prior = "IW"`; confirmed by
`formals(stan_gvar)$cov_prior`. The README and the examples that set
`cov_prior` use `"IW"`.
The existing tests that do not set `cov_prior` therefore test the Wishart
model.

**Options.** Documentation only (done, commit 757ceb9). Changing the code
default to `"LKJ"` instead would switch every call without an explicit
`cov_prior` to the model with the problematic default prior from section 1
and is not recommended.
