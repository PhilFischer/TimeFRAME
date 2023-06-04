data {
  int<lower=1> d;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> L;
  matrix[d,N] x;
  matrix[d,K] b;
  matrix[d,K] sdev;
  vector[d] sigma;
  matrix[L > 0 ? d : 0, L > 0 ? L : 0] frac;
  matrix[L > 0 ? d : 0, L > 0 ? L : 0] frac_sdev;
}
parameters {
  simplex[K] f;
  vector<lower=0, upper=1>[L] r;
  matrix<lower=-0.5, upper=0.5>[d,K] s;
  matrix[L > 0 ? d : 0, L > 0 ? L : 0] a;
}
transformed parameters {
  vector[d] mu;
  mu = (b + s .* sdev) * f;
  if(L > 0) mu += a * log(r);
}
model {
  if (L > 0) for (i in 1:d) a[i] ~ normal(frac[i], frac_sdev[i]);
  for (i in 1:d) x[i] ~ normal(mu[i], sigma[i]);
}
