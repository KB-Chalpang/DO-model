// Occupancy model with misclassification (DN-1-M, one visit per site) for the New Zealand case study.
//   z_i ~ Bernoulli(psi_i)            latent occupancy
//   p11 = Pr(y=1 | z=1)               sensitivity
//   p10 = Pr(y=1 | z=0)               false-positive rate
// The likelihood sums over the latent state z_i, so the visits at a site share the same z_i.
data {
  int<lower=1> N;                                  // number of sites
  int<lower=1> K;                                  // number of covariates
  int<lower=1> T_visits;                           // visits per site
  array[N, T_visits] int<lower=0, upper=1> y;      // detection (1) / non-detection (0) at each visit
  matrix[N, K] x;                                  // covariates
}
parameters {
  real beta0;                                      // intercept
  vector[K] beta;                                  // regression coefficients
  real<lower=0, upper=1> p11;
  real<lower=0, upper=1> p10;
}
transformed parameters {
  vector[N] psi = inv_logit(beta0 + x * beta);     // occupancy probability
}
model {
  // Priors
  beta0 ~ normal(0, 2);
  beta ~ normal(0, 1);
  p11 ~ beta(5, 2);                                // high sensitivity
  p10 ~ beta(2, 20);                               // low false positives
  // Soft constraint: p11 > p10
  {
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  // Likelihood: mixture over z_i = 1 (occupied) and z_i = 0 (unoccupied)
  for (i in 1:N) {
    real log_z1 = log(psi[i]);
    real log_z0 = log1m(psi[i]);
    for (t in 1:T_visits) {
      if (y[i, t] == 1) {
        log_z1 += log(p11);
        log_z0 += log(p10);
      } else {
        log_z1 += log1m(p11);
        log_z0 += log1m(p10);
      }
    }
    target += log_sum_exp(log_z1, log_z0);
  }
}
generated quantities {
  vector[N] Pr_z_eq_1;                             // posterior Pr(z_i = 1 | y_i)
  for (i in 1:N) {
    real log_z1 = log(psi[i]);
    real log_z0 = log1m(psi[i]);
    for (t in 1:T_visits) {
      if (y[i, t] == 1) {
        log_z1 += log(p11);
        log_z0 += log(p10);
      } else {
        log_z1 += log1m(p11);
        log_z0 += log1m(p10);
      }
    }
    Pr_z_eq_1[i] = inv_logit(log_z1 - log_z0);
  }
}
