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
  vector<lower=0>[d] eta;
  real<lower=0> rho;
}
transformed data {
  matrix[d, N] X_mean;
  matrix[d, N] X_scale;
  matrix[d, N] sdev;

  matrix[N, N] LL;
  matrix[N, N] G = gp_exp_quad_cov(t, 1.0, rho);
  for(i in 1:N) G[i,i] += 1e-9;
  LL = cholesky_decompose(G);

  for(i in 1:d) X_mean[i] = rep_row_vector(mean(X[i]), N);
  for(i in 1:d) X_scale[i] = rep_row_vector(sd(X[i]), N);
  for(i in 1:d) sdev[i] = sqrt(square(X_sd[i]) + square(eta[i])) / sqrt(2);
}
parameters {
   simplex[K] f[N];
  matrix<lower=0, upper=1>[L, N] r;
  matrix[d, N] E;
  matrix<lower=-0.5, upper=0.5>[d, K] baux;
  matrix[L>0 ? d:0, L] A;
}
transformed parameters {
  matrix[d, N] W;
  matrix[d, K] S;
  matrix[d, N] mu;

  W = X_mean + X_scale .* E * LL';
  S = b[1:d,] + baux .* b[(d+1):(2*d),];
  for (i in 1:N) mu[,i] = S * f[i];
  if (L > 0) mu += A * log(r);
}
model {
  for(i in 1:d) target += std_normal_lpdf(E[i]);
  for(l in 1:L) target += normal_lpdf(A[,l] | frac[1:d,l], frac[(d+1):(2*d),l]);
  for(i in 1:d) target += normal_lpdf(W[i] | mu[i], sdev[i]);
  for(i in 1:d) target += normal_lpdf(X[i] | W[i], sdev[i]);
}
