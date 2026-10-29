# Compile Stan models for the simulation study.
#
# Usage (compile models before fitting or running simulation replicates):
#   Rscript compile_stan_models.R
#   Rscript compile_stan_models.R /path/to/model/output
# The output directory defaults to the current working directory. Fitting
# scripts must read the saved model objects from the chosen output directory.
# Keep compiled executables in place: saved CmdStanModel objects refer to them.
# Compilation requires cmdstanr, a configured CmdStan installation and a C++ toolchain.
#
# Model families:
#   true: Bernoulli logistic regression; supply true z as y for an oracle fit.
#         Supplying observed records instead fits a GLM to those records.
#   RL:   occupancy model with false positives and repeated visits. T_visits is
#         supplied in the data; these models are not restricted to one visit.
#   do:   detection-only model with unvisited sites represented by zero records.
#
# Misspecification scenarios:
#   1: two predictors (important predictor omitted).
#   2: three predictors (wrong predictor supplied by the data preparation code).
#      The Stan formula alone does not implement replacement of a predictor.
#   3: four predictors (additional irrelevant predictor).
#
# Model assumptions:
#   RL scenario 3 assigns no explicit prior to beta4 (an implicit flat prior).
#   Detection-only models estimate q rather than treating it as known data.
#   Their likelihood depends on theta11 and theta10 and cannot separately
#   identify q and conditional detection probabilities without additional information.

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Install and configure cmdstanr before running this script.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript compile_stan_models.R [output_dir]")
output_dir <- if (length(args)) path.expand(args[[1L]]) else getwd()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

# Save each model immediately after compiling, so earlier successes are retained
# if compilation of a later model fails. Existing source/RDS files are replaced.
compile_and_save <- function(stan_code, model_name) {
  stan_file <- file.path(output_dir, paste0(model_name, ".stan"))
  writeLines(stan_code, stan_file)
  message("Compiling ", model_name)
  model <- cmdstanr::cmdstan_model(stan_file)
  saveRDS(model, file.path(output_dir, paste0(model_name, ".rds")))
  model
}

# 01. Bernoulli logistic regression: misspec1
# ------------------------------------------------------------------------
stan_model_true_misspec1 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> y;
  vector[N] x1;
  vector[N] x2;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
}

// Prior contributions followed by the observation likelihood.
model {
  vector[N] p_obs = inv_logit(beta0 + beta1 * x1 + beta2 * x2);
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  // One binary outcome per site with the model-specific marginal probability.
  y ~ bernoulli(p_obs);
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] p_obs_out;
  p_obs_out = inv_logit(beta0 + beta1 * x1 + beta2 * x2);
}
"

model_true_misspec1 <- compile_and_save(
  stan_model_true_misspec1, "model_true_misspec1"
)

# 02. Bernoulli logistic regression: misspec2
# ------------------------------------------------------------------------
stan_model_true_misspec2 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> y;
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
}

// Prior contributions followed by the observation likelihood.
model {
  vector[N] p_obs = inv_logit(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3);
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  // One binary outcome per site with the model-specific marginal probability.
  y ~ bernoulli(p_obs);
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] p_obs_out;
  p_obs_out = inv_logit(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3);
}
"

model_true_misspec2 <- compile_and_save(
  stan_model_true_misspec2, "model_true_misspec2"
)

# 03. Bernoulli logistic regression: misspec3
# ------------------------------------------------------------------------
stan_model_true_misspec3 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> y;
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;
  vector[N] x4;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;
}

// Prior contributions followed by the observation likelihood.
model {
  vector[N] p_obs = inv_logit(beta0 +
  beta1 * x1 +
  beta2 * x2 +
  beta3 * x3 +
  beta4 * x4);
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  beta4 ~ normal(0, 1);
  // One binary outcome per site with the model-specific marginal probability.
  y ~ bernoulli(p_obs);
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] p_obs_out;
  p_obs_out = inv_logit(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + beta4 * x4);
}
"

model_true_misspec3 <- compile_and_save(
  stan_model_true_misspec3, "model_true_misspec3"
)

# 04. Repeated-visit occupancy model: misspec1
# ------------------------------------------------------------------------
stan_modelRL_misspec1 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=1> N;
  int<lower=1> T_visits;
  // Include sites with observed visits only; these data cannot contain NA.
  array[N, T_visits] int<lower=0, upper=1> y;
  // Rows are sites; columns follow the predictor order in the regression.
  matrix[N, 2] x;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
  // Conditional on a visit: p11 is sensitivity; p10 is the false-positive rate.
  real<lower=0, upper=1> p11;
  real<lower=0, upper=1> p10;
}

// Deterministic probabilities derived from the sampled parameters.
transformed parameters {
  // Probability that each site is occupied, before conditioning on records.
  vector[N] psi;
  // Marginal probability of one positive record at each site.
  vector[N] p_obs;
  psi = inv_logit(beta0
      + beta1 * x[, 1]
      + beta2 * x[, 2]);
  p_obs = psi * p11 + (1 - psi) * p10;
}

// Prior contributions followed by the observation likelihood.
model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  p11 ~ beta(8, 2);
  p10 ~ beta(1, 19);
  {
    // Truncation excludes p11 < p10; this is a hard ordering restriction.
    // The normal density additionally weights the nonnegative difference.
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  // All visits at a site share the same latent occupancy state.
  // Multiply visit probabilities within each state, then sum over states.
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
    // Sum over the two possible occupancy states on the log scale.
    target += log_sum_exp(log_z1, log_z0);
  }
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] log_lik;
  vector[N] Pr_z_eq_1;
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
    // Joint log probability of this site's complete visit history.
    log_lik[i] = log_sum_exp(log_z1, log_z0);
    // Occupancy probability conditional on this site's observed history.
    Pr_z_eq_1[i] = inv_logit(log_z1 - log_z0);
  }
}
"

mystanmodelRL_misspec1 <- compile_and_save(
  stan_modelRL_misspec1, "mystanmodelRL_misspec1"
)

# 05. Repeated-visit occupancy model: misspec2
# ------------------------------------------------------------------------
stan_modelRL_misspec2 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=1> N;
  int<lower=1> T_visits;
  // Include sites with observed visits only; these data cannot contain NA.
  array[N, T_visits] int<lower=0, upper=1> y;
  // Rows are sites; columns follow the predictor order in the regression.
  matrix[N, 3] x;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  // Conditional on a visit: p11 is sensitivity; p10 is the false-positive rate.
  real<lower=0, upper=1> p11;
  real<lower=0, upper=1> p10;
}

// Deterministic probabilities derived from the sampled parameters.
transformed parameters {
  // Probability that each site is occupied, before conditioning on records.
  vector[N] psi;
  // Marginal probability of one positive record at each site.
  vector[N] p_obs;
  psi = inv_logit(beta0
      + beta1 * x[, 1]
      + beta2 * x[, 2]
      + beta3 * x[, 3]);
  p_obs = psi * p11 + (1 - psi) * p10;
}

// Prior contributions followed by the observation likelihood.
model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  p11 ~ beta(8, 2);
  p10 ~ beta(1,19);
  {
    // Truncation excludes p11 < p10; this is a hard ordering restriction.
    // The normal density additionally weights the nonnegative difference.
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  // All visits at a site share the same latent occupancy state.
  // Multiply visit probabilities within each state, then sum over states.
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
    // Sum over the two possible occupancy states on the log scale.
    target += log_sum_exp(log_z1, log_z0);
  }
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] log_lik;
  vector[N] Pr_z_eq_1;
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
    // Joint log probability of this site's complete visit history.
    log_lik[i] = log_sum_exp(log_z1, log_z0);
    // Occupancy probability conditional on this site's observed history.
    Pr_z_eq_1[i] = inv_logit(log_z1 - log_z0);
  }
}
"

mystanmodelRL_misspec2 <- compile_and_save(
  stan_modelRL_misspec2, "mystanmodelRL_misspec2"
)

# 06. Repeated-visit occupancy model: misspec3
# ------------------------------------------------------------------------
stan_modelRL_misspec3 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=1> N;
  int<lower=1> T_visits;
  // Include sites with observed visits only; these data cannot contain NA.
  array[N, T_visits] int<lower=0, upper=1> y;
  // Rows are sites; columns follow the predictor order in the regression.
  matrix[N, 4] x;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  // Conditional on a visit: p11 is sensitivity; p10 is the false-positive rate.
  real<lower=0, upper=1> p11;
  real<lower=0, upper=1> p10;
}

// Deterministic probabilities derived from the sampled parameters.
transformed parameters {
  // Probability that each site is occupied, before conditioning on records.
  vector[N] psi;
  // Marginal probability of one positive record at each site.
  vector[N] p_obs;
  psi = inv_logit(beta0
      + beta1 * x[, 1]
      + beta2 * x[, 2]
      + beta3 * x[, 3]
      + beta4 * x[, 4]);
  p_obs = psi * p11 + (1 - psi) * p10;
}

// Prior contributions followed by the observation likelihood.
model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  // beta4 has an implicit flat prior; no prior density is added for this coefficient.
  p11 ~ beta(8, 2);
  p10 ~ beta(1, 19);
  {
    // Truncation excludes p11 < p10; this is a hard ordering restriction.
    // The normal density additionally weights the nonnegative difference.
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  // All visits at a site share the same latent occupancy state.
  // Multiply visit probabilities within each state, then sum over states.
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
    // Sum over the two possible occupancy states on the log scale.
    target += log_sum_exp(log_z1, log_z0);
  }
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] log_lik;
  vector[N] Pr_z_eq_1;
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
    // Joint log probability of this site's complete visit history.
    log_lik[i] = log_sum_exp(log_z1, log_z0);
    // Occupancy probability conditional on this site's observed history.
    Pr_z_eq_1[i] = inv_logit(log_z1 - log_z0);
  }
}
"

mystanmodelRL_misspec3 <- compile_and_save(
  stan_modelRL_misspec3, "mystanmodelRL_misspec3"
)

# 07. Detection-only model: misspec1
# ------------------------------------------------------------------------
stan_model_do_misspec1 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=0> N;
  // One record per site; zero includes unvisited sites in the DO-M analysis.
  array[N] int<lower=0, upper=1> y;
  vector[N] x1;
  vector[N] x2;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
  // Site-visit probability; the lower bound avoids division by zero.
  real<lower=0.001, upper=1> q;
  // Record probabilities including visitation: theta11 = q*p11, theta10 = q*p10.
  real<lower=0, upper=q> theta11;
  real<lower=0, upper=q> theta10;
}

// Deterministic probabilities derived from the sampled parameters.
transformed parameters {
  // Conditional detection probabilities, given a visit and occupancy state.
  real p11 = theta11 / q;
  real p10 = theta10 / q;
  vector[N] psi = inv_logit(
  beta0 + beta1 * x1 + beta2 * x2
  );
  vector[N] p_obs = theta11 . * psi + theta10 . * (1 - psi);
}

// Prior contributions followed by the observation likelihood.
model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  q ~ beta(1, 1);
  target += beta_lpdf(theta11 / q | 8, 2);
  target += beta_lpdf(theta10 / q | 1, 19);
  // Jacobian for (theta11, theta10) -> (p11, p10), with q unchanged.
  // Each division by q contributes -log(q) to the log density.
  target += -2 * log(q);
  {
    // Truncation excludes p11 < p10; this is a hard ordering restriction.
    // The normal density additionally weights the nonnegative difference.
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  // One binary outcome per site with the model-specific marginal probability.
  y ~ bernoulli(p_obs);
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] log_lik;
  vector[N] Pr_z_eq_1;
  array[N] int<lower=0, upper=1> ypred;
  for (n in 1:N) {
    log_lik[n] = bernoulli_lpmf(y[n] | p_obs[n]);
    // Replicated record for posterior predictive checks.
    ypred[n] = bernoulli_rng(p_obs[n]);
    {
      real log_prob_z1 = log(psi[n]);
      real log_prob_z0 = log1m(psi[n]);
      if (y[n] == 1) {
        log_prob_z1 += log(theta11);
        log_prob_z0 += log(theta10);
      } else {
        log_prob_z1 += log1m(theta11);
        log_prob_z0 += log1m(theta10);
      }
      // Bayes' rule using the two occupancy-specific record probabilities.
      Pr_z_eq_1[n] = inv_logit(log_prob_z1 - log_prob_z0);
    }
  }
}
"

mystanmodel_do_misspec1 <- compile_and_save(
  stan_model_do_misspec1, "mystanmodel_do_misspec1"
)

# 08. Detection-only model: misspec2
# ------------------------------------------------------------------------
stan_model_do_misspec2 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=0> N;
  // One record per site; zero includes unvisited sites in the DO-M analysis.
  array[N] int<lower=0, upper=1> y;
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  // Site-visit probability; the lower bound avoids division by zero.
  real<lower=0.001, upper=1> q;
  // Record probabilities including visitation: theta11 = q*p11, theta10 = q*p10.
  real<lower=0, upper=q> theta11;
  real<lower=0, upper=q> theta10;
}

// Deterministic probabilities derived from the sampled parameters.
transformed parameters {
  // Conditional detection probabilities, given a visit and occupancy state.
  real p11 = theta11 / q;
  real p10 = theta10 / q;
  vector[N] psi = inv_logit(
  beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3
  );
  vector[N] p_obs = theta11 . * psi + theta10 . * (1 - psi);
}

// Prior contributions followed by the observation likelihood.
model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  q ~ beta(1, 1);
  target += beta_lpdf(theta11 / q | 8, 2);
  target += beta_lpdf(theta10 / q | 1, 19);
  // Jacobian for (theta11, theta10) -> (p11, p10), with q unchanged.
  // Each division by q contributes -log(q) to the log density.
  target += -2 * log(q);
  {
    // Truncation excludes p11 < p10; this is a hard ordering restriction.
    // The normal density additionally weights the nonnegative difference.
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  // One binary outcome per site with the model-specific marginal probability.
  y ~ bernoulli(p_obs);
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] log_lik;
  vector[N] Pr_z_eq_1;
  array[N] int<lower=0, upper=1> ypred;
  real mean_z;
  real mean_p_obs;
  for (n in 1:N) {
    log_lik[n] = bernoulli_lpmf(y[n] | p_obs[n]);
    // Replicated record for posterior predictive checks.
    ypred[n] = bernoulli_rng(p_obs[n]);
    {
      real log_prob_z1 = log(psi[n]);
      real log_prob_z0 = log1m(psi[n]);
      if (y[n] == 1) {
        log_prob_z1 += log(theta11);
        log_prob_z0 += log(theta10);
      } else {
        log_prob_z1 += log1m(theta11);
        log_prob_z0 += log1m(theta10);
      }
      // Bayes' rule using the two occupancy-specific record probabilities.
      Pr_z_eq_1[n] = inv_logit(log_prob_z1 - log_prob_z0);
    }
  }
  mean_z     = mean(Pr_z_eq_1);
  mean_p_obs = mean(p_obs);
}
"

mystanmodel_do_misspec2 <- compile_and_save(
  stan_model_do_misspec2, "mystanmodel_do_misspec2"
)

# 09. Detection-only model: misspec3
# ------------------------------------------------------------------------
stan_model_do_misspec3 <- "
// Observed outcomes and predictor values supplied as Stan data.
data {
  int<lower=0> N;
  // One record per site; zero includes unvisited sites in the DO-M analysis.
  array[N] int<lower=0, upper=1> y;
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;
  vector[N] x4;
}

// Unknown quantities estimated from the data.
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  // Site-visit probability; the lower bound avoids division by zero.
  real<lower=0.001, upper=1> q;
  // Record probabilities including visitation: theta11 = q*p11, theta10 = q*p10.
  real<lower=0, upper=q> theta11;
  real<lower=0, upper=q> theta10;
}

// Deterministic probabilities derived from the sampled parameters.
transformed parameters {
  // Conditional detection probabilities, given a visit and occupancy state.
  real p11 = theta11 / q;
  real p10 = theta10 / q;
  vector[N] psi = inv_logit(
  beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + beta4 * x4
  );
  vector[N] p_obs = theta11 . * psi + theta10 . * (1 - psi);
}

// Prior contributions followed by the observation likelihood.
model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  beta4 ~ normal(0, 1);
  q ~ beta(1, 1);
  target += beta_lpdf(theta11 / q | 8, 2);
  target += beta_lpdf(theta10 / q | 1, 19);
  // Jacobian for (theta11, theta10) -> (p11, p10), with q unchanged.
  // Each division by q contributes -log(q) to the log density.
  target += -2 * log(q);
  {
    // Truncation excludes p11 < p10; this is a hard ordering restriction.
    // The normal density additionally weights the nonnegative difference.
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  // One binary outcome per site with the model-specific marginal probability.
  y ~ bernoulli(p_obs);
}

// Derived outputs evaluated for each posterior draw; they do not affect fitting.
generated quantities {
  vector[N] log_lik;
  vector[N] Pr_z_eq_1;
  array[N] int<lower=0, upper=1> ypred;
  real mean_z;
  real mean_p_obs;
  for (n in 1:N) {
    log_lik[n] = bernoulli_lpmf(y[n] | p_obs[n]);
    // Replicated record for posterior predictive checks.
    ypred[n] = bernoulli_rng(p_obs[n]);
    {
      real log_prob_z1 = log(psi[n]);
      real log_prob_z0 = log1m(psi[n]);
      if (y[n] == 1) {
        log_prob_z1 += log(theta11);
        log_prob_z0 += log(theta10);
      } else {
        log_prob_z1 += log1m(theta11);
        log_prob_z0 += log1m(theta10);
      }
      // Bayes' rule using the two occupancy-specific record probabilities.
      Pr_z_eq_1[n] = inv_logit(log_prob_z1 - log_prob_z0);
    }
  }
  mean_z     = mean(Pr_z_eq_1);
  mean_p_obs = mean(p_obs);
}
"

mystanmodel_do_misspec3 <- compile_and_save(
  stan_model_do_misspec3, "mystanmodel_do_misspec3"
)

message("All 9 models compiled and saved in: ", output_dir)
