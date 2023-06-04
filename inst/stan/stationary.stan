data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> L;
  matrix[2,N] x;
  matrix[2,K] b;
  matrix[2,K] sdev;
  vector[2] sigma;
  matrix[L > 0 ? 2 : 0, L > 0 ? L : 0] frac;
  matrix[L > 0 ? 2 : 0, L > 0 ? L : 0] frac_sdev;
}
parameters {
  simplex[K] f;
  vector<lower=0, upper=1>[L] r;
  matrix<lower=-0.5, upper=0.5>[2,K] s;
  matrix[L > 0 ? 2 : 0, L > 0 ? L : 0] a;
}
transformed parameters {
  vector[2] mu;
  mu = (b + s .* sdev) * f;
  if(L > 0) mu += a * log(r);
}
model {
  if (L > 0) for (d in 1:2) a[d] ~ normal(frac[d], frac_sdev[d]);
  for (d in 1:2) x[d] ~ normal(mu[d], sigma[d]);
}
