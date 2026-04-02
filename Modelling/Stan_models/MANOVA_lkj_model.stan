
data {
int<lower=1> N; // number of observations
int<lower=1> P; // number of dimensions
int<lower=1> L; // number of levels
vector[P] mu; // prior mean for theta
matrix[P, P] B; // prior covariance for theta
matrix[L,P] beta_mu; // prior means for beta
cov_matrix[P] beta_cov[L]; // prior covariances for beta
int<lower=1, upper=L> letters[N]; // categorical predictor
matrix[N, P] y; // observed data
real loc;   // 
real sc;    // 
real eta;  // eta lkj
}

parameters {
vector[P] theta; // mean vector
matrix[P, L-1] beta_raw; // raw coefficients for levels 2 to L
cholesky_factor_corr[P] Omega; // Cholesky factor of correlation matrix
vector<lower=0>[P] sigma; // scale parameters
}

transformed parameters {
//matrix[P, P] Sigma; // covariance matrix
//matrix[P, P] diag_Omega;   // Cholesky factor with scaling
//diag_Omega = diag_pre_multiply(sigma, Omega);
//Sigma = diag_Omega * diag_Omega';
matrix[P, L] beta;
matrix[P, P] L_Sigma;                     // Cholesky factor of covariance

beta = append_col(rep_vector(0, P), beta_raw);  // set baseline (level 1) to 0
L_Sigma = diag_pre_multiply(sigma, Omega);      // Cholesky factor of full covariance
}

model {
  // Priors
  //target += cauchy_lpdf(sigma | loc, sc);
  target += lognormal_lpdf(sigma | loc, sc);
  target += lkj_corr_cholesky_lpdf(Omega | eta);
  target += multi_normal_lpdf(theta | mu, B);
  for (l in 2:L) {
    target += multi_normal_lpdf(beta_raw[, l-1] | beta_mu[l]', beta_cov[l]);
  }

  // Likelihood
  for (n in 1:N) {
    target += multi_normal_cholesky_lpdf(y[n] | theta + beta[, letters[n]], L_Sigma);
  }
}


