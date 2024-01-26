// Save this file as inst/stan/test.stan
// for experimental pkg development purposes
data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real intercept;
  real beta;
  real<lower=0> sigma;
}
model {
  // ... priors, etc.

  y ~ normal(intercept + beta * x, sigma);
}
