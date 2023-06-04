data {
  int<lower=1> d;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> L;
  matrix[(2*d), K] b;
  matrix[L>0 ? (2*d) : 0, L] frac;
  matrix[d, N] X;
  matrix[d, N] X_sd;
  real<lower=0> sigma;
  real<lower=0> alpha;
  real<lower=0> alpha_r;
}
transformed data {
  matrix[d, N] sdev = sqrt(square(X_sd) + square(sigma));
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
  if(alpha != 1) for(i in 1:N) f[i] ~ dirichlet(rep_vector(alpha, K));
  if(alpha_r != 1) for(l in 1:L) r[l] ~ gamma(alpha_r, 1);
  for(l in 1:L) A[,l] ~ normal(frac[1:d,l], frac[(d+1):(2*d),l]);
  for(i in 1:d) X[i] ~ normal(mu[i], sdev[i]);
}
