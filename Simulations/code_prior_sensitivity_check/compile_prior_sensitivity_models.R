# Compile six detection-only models with alternative priors on p11 and p10.
# Run once before fitting: Rscript compile_prior_sensitivity_models.R
# All paths are relative to the working directory used for submission.

library(cmdstanr)

# 1. Output directory and prior specifications 

output_dir <- "prior_sensitivity_compiled_models"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create directory: ", output_dir)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

priors <- data.frame(
  scenario = letters[1:6],
  p11_alpha = c(8, 4, 6, 18, 2, 2),
  p11_beta  = c(2, 1, 4, 2, 2, 3),
  p10_alpha = c(1, 0.5, 2, 1, 1, 2),
  p10_beta  = c(19, 9.5, 18, 39, 9, 5)
)

# 2. Shared Stan model

# Four placeholders are replaced with the selected Beta shape parameters.
# All other model terms are identical across the six scenarios.
stan_template <- r"---(
data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> y;
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;
}

parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  // The lower bound avoids division by zero in p11 and p10.
  real<lower=0.001, upper=1> q;
  // Record probabilities include visitation: theta11=q*p11, theta10=q*p10.
  real<lower=0, upper=q> theta11;
  real<lower=0, upper=q> theta10;
}

transformed parameters {
  // Detection probabilities conditional on visitation and occupancy state.
  real p11 = theta11 / q;
  real p10 = theta10 / q;
  vector[N] psi = inv_logit(
    beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3
  );
  vector[N] p_obs = theta11 .* psi + theta10 .* (1 - psi);
}

model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  q ~ beta(1, 1);

  // Scenario-specific Beta contributions for conditional detection probabilities.
  target += beta_lpdf(theta11 / q | @P11_ALPHA@, @P11_BETA@);
  target += beta_lpdf(theta10 / q | @P10_ALPHA@, @P10_BETA@);
  // Jacobian for (theta11, theta10) -> (p11, p10), holding q fixed.
  target += -2 * log(q);
  {
    // Require p11 >= p10 and additionally weight their nonnegative difference.
    real delta = p11 - p10;
    delta ~ normal(0.2, 0.5) T[0, ];
  }
  y ~ bernoulli(p_obs);
}

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
      // Occupancy probability conditional on this site's observed record.
      Pr_z_eq_1[n] = inv_logit(log_prob_z1 - log_prob_z0);
    }
  }
}
)---"

# 3. Write, compile and save each scenario ------------------------------------

for (i in seq_len(nrow(priors))) {
  prior <- priors[i, ]
  replacements <- c(
    "@P11_ALPHA@" = prior$p11_alpha, "@P11_BETA@" = prior$p11_beta,
    "@P10_ALPHA@" = prior$p10_alpha, "@P10_BETA@" = prior$p10_beta
  )
  stan_code <- stan_template
  for (token in names(replacements)) {
    stan_code <- gsub(token, as.character(replacements[[token]]),
                      stan_code, fixed = TRUE)
  }
  model_name <- paste0("model_do_prior_", prior$scenario)
  stan_file <- file.path(output_dir, paste0(model_name, ".stan"))
  rds_file <- file.path(output_dir, paste0(model_name, ".rds"))
  writeLines(stan_code, stan_file)

  # Unique source and executable names prevent one scenario from replacing
  # another. Compile sequentially to keep memory use predictable.
  cat("Compiling scenario", prior$scenario, "\n")
  model <- cmdstanr::cmdstan_model(stan_file, force_recompile = TRUE)
  saveRDS(model, rds_file)
  cat("Saved", rds_file, "\n")
}

# A readable record of the Beta shape parameters used to generate the sources.
write.csv(priors, file.path(output_dir, "prior_specifications.csv"), row.names = FALSE)
cat("Finished compiling all six models in", output_dir, "\n")
