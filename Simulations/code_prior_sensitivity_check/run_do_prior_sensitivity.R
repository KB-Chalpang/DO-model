# Detection-only model: sensitivity to priors on p11 and p10.
# Usage: Rscript run_do_prior_sensitivity.R <scenario> <replicate>
# Example: Rscript run_do_prior_sensitivity.R a 1
# Scenarios a-f use separate, previously compiled Stan models (see README).

library(cmdstanr)
library(posterior)

# 1. Scenario and replicate

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript run_do_prior_sensitivity.R <a-f> <replicate>")
}
scenario <- args[1]
task_id <- suppressWarnings(as.numeric(args[2]))
if (!scenario %in% letters[1:6]) stop("Scenario must be a, b, c, d, e or f.")
if (!is.finite(task_id) || task_id < 1 || task_id != floor(task_id) ||
    task_id > .Machine$integer.max - 9000) {
  stop("Replicate must be a positive integer within the supported seed range.")
}
task_id <- as.integer(task_id)

# The table records the intended priors. The priors used for fitting are
# defined inside each compiled Stan model, not by this table.
prior_settings <- data.frame(
  scenario = letters[1:6],
  p11_alpha = c(8, 4, 6, 18, 2, 2),
  p11_beta  = c(2, 1, 4, 2, 2, 3),
  p10_alpha = c(1, 0.5, 2, 1, 1, 2),
  p10_beta  = c(19, 9.5, 18, 39, 9, 5)
)
prior <- prior_settings[prior_settings$scenario == scenario, , drop = FALSE]

# 2. Input and output paths 

# Paths beginning with ~ refer to your home directory. Other paths are
# relative to the directory from which this script is run.
covariate_file <- path.expand("~/covariates_fixed.rds")
model_dir <- "prior_sensitivity_compiled_models"
model_file <- file.path(model_dir, paste0("model_do_prior_", scenario, ".rds"))
scenario_dir <- file.path("results_do_prior_sensitivity", paste0("scenario_", scenario))
results_dir <- file.path(scenario_dir, "results")
settings_dir <- file.path(scenario_dir, "settings")

if (!file.exists(covariate_file)) stop("Covariate file not found: ", covariate_file)
if (!file.exists(model_file)) stop("Compiled model not found: ", model_file)
x <- readRDS(covariate_file)
required_columns <- c("x1", "x2", "x3")
if (!is.data.frame(x) || !all(required_columns %in% names(x)) || nrow(x) == 0L) {
  stop("Covariates must be a nonempty data frame containing x1, x2 and x3.")
}
if (!all(vapply(x[required_columns], function(v) is.numeric(v) && all(is.finite(v)),
                logical(1)))) {
  stop("Covariates x1, x2 and x3 must contain finite numeric values.")
}
model <- readRDS(model_file)
for (path in c(results_dir, settings_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Could not create directory: ", path)
}

# 3. Data-generating settings 

beta0 <- 0
beta1 <- 1.2
beta2 <- 1
beta3 <- -1
p11 <- 0.80                   # Positive-record probability at an occupied site.
p10 <- 0.05                   # False-positive probability at an unoccupied site.
q <- 0.8                      # Probability that a site is visited.
T_visits <- 2L
N <- nrow(x)

# Matching replicate IDs use the same simulated data across prior scenarios.
# These seeds also make simulation and random initial values reproducible.
simulation_seed <- 4567L + task_id
initialization_seed <- 9000 + task_id
stan_seed <- 123L + task_id

# 4. Simulate occupancy and first-visit records 

simulate_records <- function(x, beta0, beta1, beta2, beta3, p11, p10, q, T_visits) {
  N <- nrow(x)
  psi <- plogis(beta0 + beta1 * x$x1 + beta2 * x$x2 + beta3 * x$x3)
  z <- rbinom(N, size = 1, prob = psi)
  visited <- rbinom(N, size = 1, prob = q)

  # A zero is stored for every unvisited site. For visited sites, the model
  # uses only the first visit's misclassified record.
  y <- integer(N)
  for (i in which(visited == 1L)) {
    for (t in seq_len(T_visits)) {
      # Retain the error-free draw and later visits to preserve the simulation
      # random-number sequence, although their outputs are not used for fitting.
      invisible(rbinom(1, size = 1, prob = z[i] * p11))
      if (z[i] == 1L) {
        record <- rbinom(1, size = 1, prob = p11)
      } else {
        record <- rbinom(1, size = 1, prob = p10)
      }
      if (t == 1L) y[i] <- record
    }
  }
  y
}

set.seed(simulation_seed)
y <- simulate_records(x, beta0, beta1, beta2, beta3, p11, p10, q, T_visits)
stan_data <- list(N = N, y = y, x1 = x$x1, x2 = x$x2, x3 = x$x3)

# 5. Fit the selected model 

# CmdStanR calls this function for each chain. q exceeds theta11 and theta10
# over these ranges, and theta11 exceeds theta10.
init_fun_do <- function() {
  list(
    beta0 = rnorm(1, 0, 1),
    beta1 = rnorm(1, 0, 1),
    beta2 = rnorm(1, 0, 1),
    beta3 = rnorm(1, 0, 1),
    q = runif(1, 0.5, 0.8),
    theta11 = runif(1, 0.15, 0.45),
    theta10 = runif(1, 0.001, 0.1)
  )
}

cat("Scenario:", scenario, " Replicate:", task_id, " Model:", model_file, "\n")
set.seed(initialization_seed)
start_time <- Sys.time()
fit <- model$sample(
  data = stan_data,
  seed = stan_seed,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  refresh = 0,
  max_treedepth = 15,
  adapt_delta = 0.99,
  init = init_fun_do
)
elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# 6. Summarise posterior draws 

params <- c("beta0", "beta1", "beta2", "beta3", "p11", "p10", "theta11", "theta10")
# Generating values stay fixed across scenarios; prior means do not replace them.
true_vals <- c(beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3,
               p11 = p11, p10 = p10, theta11 = q * p11, theta10 = q * p10)
draws <- fit$draws(variables = params)
summary_stats <- posterior::summarise_draws(draws, mean, sd, rhat, ess_bulk, ess_tail)
quantiles <- posterior::summarise_draws(
  draws, ~quantile(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
)
# Match by parameter name so the saved table has a consistent row order.
s <- match(params, summary_stats$variable)
u <- match(params, quantiles$variable)
results <- data.frame(
  scenario = scenario,
  replicate = task_id,
  parameter = params,
  true_value = unname(true_vals[params]),
  post_mean = summary_stats$mean[s],
  post_sd = summary_stats$sd[s],
  q2.5 = quantiles$`2.5%`[u],
  q25 = quantiles$`25%`[u],
  q50 = quantiles$`50%`[u],
  q75 = quantiles$`75%`[u],
  q97.5 = quantiles$`97.5%`[u],
  rhat = summary_stats$rhat[s],
  ess_bulk = summary_stats$ess_bulk[s],
  ess_tail = summary_stats$ess_tail[s]
)

# 7. Save one summary and one settings file per replicate 

saveRDS(results, file.path(results_dir, sprintf("post_means_%03d.rds", task_id)))
saveRDS(
  list(
    scenario = scenario, replicate = task_id, intended_priors = prior,
    generating_values = c(true_vals, q = q), N = N, T_visits = T_visits,
    simulation_seed = simulation_seed, initialization_seed = initialization_seed,
    stan_seed = stan_seed, elapsed_seconds = elapsed_seconds,
    covariate_file = normalizePath(covariate_file),
    model_file = normalizePath(model_file), session_info = sessionInfo()
  ),
  file.path(settings_dir, sprintf("settings_%03d.rds", task_id))
)
cat("Finished scenario", scenario, "replicate", task_id,
    "in", round(elapsed_seconds, 1), "seconds.\n")
