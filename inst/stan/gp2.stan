data {
  int<lower=1> d;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> L;
  matrix[(2*d), K] b;
  matrix[L>0 ? (2*d) : 0, L] frac;
  real t[N];
  matrix[d, N] X;
  matrix[d, N] X_sd;
  real<lower=0> sigma;
  real<lower=0> rho;
}
transformed data {
  matrix[d, N] X_vars = (square(X_sd) + square(sigma)) / 2;
  matrix[d, N] X_means;
  matrix[N, N] G = gp_exp_quad_cov(t, 1.0, rho);
  matrix[N, N] LL[d];
  vector[N] means2[d];
  matrix[N, N] covs1[d];
  matrix[N, N] covs2[d];

  for(i in 1:d) X_means[i] = rep_row_vector(mean(X[i]), N);
  for(i in 1:N) G[i,i] += 1e-9;
  for(i in 1:d) {
    LL[i] = cholesky_decompose(add_diag(variance(X[i]) * G, X_vars[i]));
    means2[i] = mdivide_left_tri_low(LL[i], X_means[i]');
    covs1[i] = mdivide_left_tri_low(LL[i], diag_matrix(X_vars[i]'));
    covs2[i] = mdivide_left_tri_low(LL[i], variance(X[i]) * G);
  }
}
parameters {
  simplex[K] f[N];
  matrix<lower=0, upper=1>[L, N] r;
  matrix<lower=-0.5, upper=0.5>[d, K] baux;
  matrix[L>0 ? d:0, L] A;
}
transformed parameters {
  matrix[d, K] S;
  matrix[d, N] mu;

  S = b[1:d,] + baux .* b[(d+1):(2*d),];
  for (i in 1:N) mu[,i] = S * f[i];
  if (L > 0) mu += A * log(r);
}
model {
  vector[N] mean1;

  for(l in 1:L) target += normal_lpdf(A[,l] | frac[1:d,l], frac[(d+1):(2*d),l]);
  for(i in 1:d) {
    mean1 = mdivide_left_tri_low(LL[i], mu[i]');
    target += multi_normal_lpdf(X[i] | covs2[i]' * mean1 + covs1[i]' * means2[i], add_diag(covs1[i]' * covs2[i], X_vars[i]));
    target += multi_normal_cholesky_lpdf(mu[i] | X_means[i], LL[i]);
  }
}
