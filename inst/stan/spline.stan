functions {
  matrix simplex_base(int K) {
    matrix[K, (K-1)] B;
    for(k in 1:(K-1)) {
      for(i in 1:k) {
        B[i,k] = sqrt(1.0 / (k*(k+1)));
      }
      B[k+1,k] = -sqrt(1.0 * k/(k+1));
      for(i in (k+2):K) {
        B[i,k] = 0.0;
      }
    }
    return(B);
  }
}
data {
  int<lower=1> d;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> L;
  int<lower=1> M;
  int<lower=0> M_r;
  matrix[(2*d), K] b;
  matrix[L>0 ? (2*d) : 0, L] frac;
  real t[N];
  matrix[d, N] X;
  matrix[d, N] X_sd;
  real<lower=0> sigma;
  matrix[M, N] basis;
  matrix[M_r, N] basis_r;
}
transformed data {
  matrix[K, (K-1)] SB = simplex_base(K);
  matrix[d, N] sdev;
  for(i in 1:d) sdev[i] = sqrt(square(X_sd[i]) + square(sigma));
}
parameters {
  matrix[(K-1), M] faux;
  matrix[L, M_r] raux;
  matrix<lower=-0.5, upper=0.5>[d, K] baux;
  matrix[d, L] A;
}
transformed parameters {
  matrix[K, N] f;
  matrix[L, N] r;
  matrix[d, K] S;
  matrix[d, N] mu;

  matrix[K, N] feval = SB * faux * basis;
  matrix[L, N] reval = raux * basis_r;
  for(i in 1:N) f[,i] = softmax(feval[,i]);
  for(i in 1:N) r[,i] = inv_logit(reval[,i]);
  S = b[1:d,] + baux .* b[(d+1):(2*d),];
  for(i in 1:N) mu[,i] = S * f[,i];
  if (L > 0) mu += A * log(r);
}
model {
  for(k in 1:(K-1)) target += normal_lpdf(faux[k] | 0, 1);
  for(l in 1:L) target += normal_lpdf(raux[l] | 0, 1);
  for(l in 1:L) target += normal_lpdf(A[,l] | frac[1:d,l], frac[(d+1):(2*d),l]);
  for(i in 1:d) target += normal_lpdf(X[i] | mu[i], sdev[i]);
}
