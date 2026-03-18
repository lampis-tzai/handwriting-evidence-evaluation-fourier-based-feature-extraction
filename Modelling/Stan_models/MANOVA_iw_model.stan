data {
  int<lower=1> N; // number of observations
  int<lower=1> P; // number of dimensions
  int<lower=1> L; // number of levels
  vector[P] mu; // prior mean for theta
  matrix[P, P] B; // prior covariance for theta
  matrix[L, P] beta_mu; // prior means for beta
  cov_matrix[P] beta_cov[L]; // prior covariances for beta
  int<lower=1, upper=L> letters[N]; // categorical predictor
  matrix[N, P] y; // observed data
  real nu;  // degrees of freedom for inverse Wishart
  matrix[P, P] U; // scale matrix for inverse Wishart
}

parameters {
  vector[P] theta; // mean vector
  matrix[P, L-1] beta_raw; // raw coefficients for levels 2 to L
  cov_matrix[P] Sigma; // covariance matrix
}

transformed parameters {
  matrix[P, L] beta;
  beta = append_col(rep_vector(0, P), beta_raw); // set beta[,1] to 0 and append the rest
}

model {
  // Priors
  target += inv_wishart_lpdf(Sigma | nu, U);
  target += multi_normal_lpdf(theta | mu, B);
  for (l in 2:L) {
    target += multi_normal_lpdf(beta_raw[, l-1] | beta_mu[l]', beta_cov[l]);
  }

  // Likelihood
  for (n in 1:N) {
    target += multi_normal_lpdf(y[n] | theta + beta[, letters[n]], Sigma);
  }
}
