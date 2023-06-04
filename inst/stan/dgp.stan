functions {
  matrix simplex_base(int K) {
    matrix[K, (K-1)] B;
    for(k in 1:(K-1)) {
      for(i in 1:k) B[i,k] = sqrt(1.0 / (k*(k+1)));
      B[k+1,k] = -sqrt(1.0 * k/(k+1));
      for(i in (k+2):K) B[i,k] = 0.0;
    }
    return(B);
  }
}
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
  real<lower=0> rho_r;
}
transformed data {
  matrix[d, N] sdev = sqrt(square(X_sd)+ square(sigma));
  matrix[K, (K-1)] SB = simplex_base(K);
  matrix[N, N] LL_f;
  matrix[N, N] LL_r;
  matrix[N, N] G_f = gp_exp_quad_cov(t, 1.0, rho);
  matrix[N, N] G_r = gp_exp_quad_cov(t, 1.0, rho_r);

  for(i in 1:N) G_f[i,i] += 1e-9;
  LL_f = cholesky_decompose(G_f);

  for(i in 1:N) G_r[i,i] += 1e-9;
  LL_r = cholesky_decompose(G_r);
}
parameters {
  matrix[(K-1), N] Z_f;
  matrix[L, N] Z_r;
  matrix<lower=-0.5, upper=0.5>[d, K] baux;
  matrix[d, L] A;
}
transformed parameters {
  matrix[d, K] S;
  matrix[d, N] mu;
  matrix[(K-1), N] W;
  matrix[K, N] f;
  matrix[L, N] r;

  S = b[1:d,] + baux .* b[(d+1):(2*d),];
  W = Z_f * LL_f';
  r = inv_logit(Z_r * LL_r');
  for(i in 1:N) f[,i] = softmax(SB * W[,i]);
  for (i in 1:N) mu[, i] = S * f[, i] + A * log(r[, i]);
}
model {
  for(j in 1:(K-1)) target += std_normal_lpdf(Z_f[j]);
  for(l in 1:L) target += std_normal_lpdf(Z_r[l]);
  for(l in 1:L) target += normal_lpdf(A[,l] | frac[1:d,l], frac[(d+1):(2*d),l]);
  for(i in 1:d) target += normal_lpdf(X[i] | mu[i], sdev[i]);
}
