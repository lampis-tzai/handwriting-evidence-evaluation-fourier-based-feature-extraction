data {
  int<lower=1> N;           // Number of observations
  int<lower=1> P;           // Number of variables
  matrix[N, P] y;           // Data matrix
  vector[P] mu;             // Prior mean vector
  matrix[P, P] B;           // Between writers variability
  matrix[P, P] U;           // Prior scale matrix for Inverse Wishart
  real nu;                  // Prior degrees of freedom for Inverse Wishart
}
parameters {
  vector[P] theta;          // Mean vector
  cov_matrix[P] Sigma;      // Covariance matrix
}
model {
  // Normal-Inverse-Wishart prior
  target += inv_wishart_lpdf(Sigma | nu, U);
  target += multi_normal_lpdf(theta | mu, B);

  // Likelihood
  for (n in 1:N) {
    target += multi_normal_lpdf(y[n] | theta, Sigma);
  }
}
