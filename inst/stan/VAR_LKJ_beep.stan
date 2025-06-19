////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> T; // number of time points
  array[T] int beep; // beep number
  array[T] vector[K] Y; // responses
  // Priors
  matrix[K,K] prior_Beta_loc; // locations for priors on Beta matrix
  matrix[K,K] prior_Beta_scale;  // scales for priors on Beta matrix
  matrix[K,K] prior_Rho_loc;  // locations for priors on partial correlations
  matrix[K,K] prior_Rho_scale;   // scales for priors on partial correlations
  int<lower=1> prior_Eta; // prior for LKJ
  // Forecast
  int<lower=0> ahead; // forecasted time points
  int<lower = 0, upper = 1> compute_log_lik; // compute log likelihood?
  array[ahead] vector[K] Y_future; // future responses for loglik computation

}
////////////////////////////////////////////////////////////////////////////////
transformed data{
  int first_beep = min(beep);
}
////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  matrix[K,K] Beta_raw; //
  //real mu_Beta;
  //real<lower=0> sigma_Beta;

  // Contemporaneous
  cholesky_factor_corr[K] L_Theta;
  vector[K] sigma_theta;
}
////////////////////////////////////////////////////////////////////////////////
transformed parameters{
  // Non-centered parameterization for Beta matrix
  matrix[K,K] Beta = Beta_raw .* prior_Beta_scale + prior_Beta_loc;
  //matrix[K,K] Beta = Beta_raw * sigma_Beta + mu_Beta;

  // Covariance matrix from cholesky corr matrix and SDs
  matrix[K,K] Sigma = diag_pre_multiply(exp(sigma_theta), L_Theta) *
                      diag_pre_multiply(exp(sigma_theta), L_Theta)';
  // Partial correlation matrix
  matrix[K,K] Rho;
  {
    // Precision matrix
    matrix[K,K] Theta = inverse_spd(Sigma);
    for(i in 1:K){
      for(j in 1:K){
        if(i != j){
          Rho[i,j] = -Theta[i,j] / sqrt(Theta[i,i] * Theta[j,j]);
        }else{
          Rho[i,j] = 0;
        } // end else
      } // end j
    } // end i
  }
}
////////////////////////////////////////////////////////////////////////////////
model {
  // Priors
  target+=   std_normal_lpdf(to_vector(Beta_raw));    // prior on Beta

  target+= lkj_corr_cholesky_lpdf(L_Theta | prior_Eta);
  target+=   student_t_lpdf(sigma_theta | 3,0,2);   // prior on sigma_theta
  // Priors on partial correlations
  for(i in 1:K){
    for(j in 1:K){
      if(i < j){
        // Scaled beta prior on partial correlations
        // (Rho[i,j] / 2 + 0.5 is a scaling to the unit interval)
        target+= beta_proportion_lpdf(
          Rho[i,j] / 2 + 0.5 | prior_Rho_loc[i,j], prior_Rho_scale[i,j]);
        }
      }
    }
  {
    // Cholesky decomposition of the covariance matrix
    matrix[K, K] Sigma_chol = diag_pre_multiply(exp(sigma_theta), L_Theta);
    for(t in 2:T){
      if(beep[t] > first_beep){
        vector[K] mu = Beta * Y[t-1,];
        target += multi_normal_cholesky_lpdf(Y[t,] | mu, Sigma_chol);
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
generated quantities {
  int min_beep = first_beep;
  vector[T + ahead] log_lik;

  // Initialize log_lik with zero for all time points
  for (t in 1:(T + ahead)) {
    log_lik[t] = 0;
  }

  // Cholesky decomposition of the covariance matrix
  matrix[K, K] Sigma_chol = diag_pre_multiply(exp(sigma_theta), L_Theta);
  for (t in 2:T) {
    if (beep[t] > first_beep) {
      vector[K] mu = Beta * Y[t-1,];
      // use t here meaning that ll of timepoint two is stored in log_lik[2]
      log_lik[t] = multi_normal_cholesky_lpdf(Y[t, ] | mu, Sigma_chol);
    }
  }
  if(ahead > 0){
    // Forecasting ahead steps
    array[ahead] vector[K] Y_forecast; // forecasted responses
    vector[K] current_Y = Y[T]; // initialize current_Y to the last observed value
    for (s in 1:ahead) {
      vector[K] mu = Beta * current_Y; // mu for current step
      Y_forecast[s] = multi_normal_cholesky_rng(mu, Sigma_chol); // generate forecast
      print(Y_forecast[s]);
      current_Y = Y_forecast[s]; // update current_Y to predicted value
      // forecasted log-likelihood
      if (compute_log_lik == 1) {
        log_lik[T + s] = multi_normal_cholesky_lpdf(Y_future[s] | mu, Sigma_chol);
    }
  }
  }

}
