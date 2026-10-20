// Presence-only-with-false-positives occupancy model for detection-only data.
//   z_i ~ Bernoulli(psi_i)            latent occupancy
//   q                                 probability a site is visited/recorded
//   theta11 = q * P(y=1 | z=1)        theta10 = q * P(y=1 | z=0)
data {
  int<lower=1> N;                         // number of sites
  int<lower=1> K;                         // number of covariates
  array[N] int<lower=0, upper=1> y;       // detection-only data
  matrix[N, K] x;                         // covariates
}
parameters {
  real beta0;                             // intercept
  vector[K] beta;                         // regression coefficients
  real<lower=0.001, upper=1> q;           // visit probability (uniform prior)
  real<lower=0, upper=q> theta11;         // Pr(y=1, visited | z=1)
  real<lower=0, upper=q> theta10;         // Pr(y=1, visited | z=0)
}
transformed parameters {
  real p11 = theta11 / q;                 // Pr(y=1 | z=1, visited)
  real p10 = theta10 / q;                 // Pr(y=1 | z=0, visited)
  vector[N] psi = inv_logit(beta0 + x * beta);
  vector[N] p_obs = theta11 * psi + theta10 * (1 - psi);
}
model {
  // Priors
  beta0 ~ normal(0, 2);
  beta ~ normal(0, 1);
  // Induced priors on p11 and p10 (with Jacobian for theta = p * q)
  target += beta_lpdf(p11 | 5, 2);
  target += beta_lpdf(p10 | 2, 20);
  target += -2 * log(q);
  // Soft constraint: p11 > p10
  {
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  // Likelihood
  y ~ bernoulli(p_obs);
}
generated quantities {
  array[N] int<lower=0, upper=1> ypred;   // posterior predictive draws
  vector[N] Pr_z_eq_1;                    // Pr(z_i = 1 | y_i)
  for (n in 1:N) {
    ypred[n] = bernoulli_rng(p_obs[n]);
    Pr_z_eq_1[n] = y[n] == 1
                   ? psi[n] * theta11 / p_obs[n]
                   : psi[n] * (1 - theta11) / (1 - p_obs[n]);
  }
}
