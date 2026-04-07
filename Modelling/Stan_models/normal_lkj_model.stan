
data {
int<lower=1> N;  // Number of observations
int<lower=1> P;  // Number of variables
matrix[N, P] y;  // Data matrix
vector[P] mu;    // Prior mean vector
matrix[P, P] B;  // Prior covariance matrix for the mean vector
vector[P] loc;   // real loc
vector[P] sc;    // real sc
real eta;  // eta lkj
}
parameters {
vector[P] theta;  // Mean vector
vector<lower=0>[P] sigma;  // Standard deviations
cholesky_factor_corr[P] L;  // Cholesky factor of the correlation matrix
}
transformed parameters {
//matrix[P, P] Sigma;  // Covariance matrix
//matrix[P, P] prec;   // Precision matrix
matrix[P, P] diag_L;   // Cholesky factor with scaling
diag_L = diag_pre_multiply(sigma, L);
//Sigma = diag_L * diag_L';
//prec = inverse(Sigma);
}
model {
// Priors
target += multi_normal_lpdf(theta | mu, B);
  //target += cauchy_lpdf(sigma | loc, sc);
  target += lognormal_lpdf(sigma | loc, sc);
  target += lkj_corr_cholesky_lpdf(L | eta);

  // Likelihood
  for (n in 1:N) {
    target += multi_normal_cholesky_lpdf(y[n] | theta, diag_L);
  }
}

