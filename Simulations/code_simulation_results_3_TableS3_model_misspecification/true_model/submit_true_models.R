# Simulate and fit one replicate of a true-state logistic regression model.
#
# Usage:
#   Rscript simulate_true_models.R <task_id> <scenario> [covariates] [models] [results]
# Examples:
#   Rscript simulate_true_models.R 1 misspec1
#   Rscript simulate_true_models.R 1 misspec2 covariates_fixed.rds compiled_models results_true
#
# Each invocation fits one scenario and saves one summary table. For SLURM
# arrays, pass SLURM_ARRAY_TASK_ID as task_id and allocate four CPUs per task.
# Relative paths are resolved from the working directory.
#
# Scenario   Fitted predictors   Interpretation
# misspec1   x1, x2              Omit the relevant predictor x3.
# misspec2   x1, x2, x4          Replace x3 with x4 (Stan's x3 input receives x4).
# misspec3   x1, x2, x3, x4      Include the additional predictor x4.
#
# The response is the simulated true occupancy state z at every site. No visit
# or detection process is needed because this model family observes z directly.
# The same task_id uses the same simulated z across all three scenarios when
# the covariates and generating coefficients are identical.

# 1. Dependencies and command-line arguments ---------------------------------
for (package in c("cmdstanr", "posterior")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required package is not installed: ", package)
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 5L) {
  stop(paste(
    "Usage: Rscript simulate_true_models.R <task_id> <scenario>",
    "[covariates] [models] [results]"
  ))
}

task_id <- suppressWarnings(as.numeric(args[[1L]]))
if (length(task_id) != 1L || !is.finite(task_id) || task_id < 1 ||
    task_id != floor(task_id) || task_id > .Machine$integer.max - 4569) {
  stop("task_id must be a positive integer within the supported seed range.")
}
task_id <- as.integer(task_id)
scenario <- match.arg(args[[2L]], c("misspec1", "misspec2", "misspec3"))
covariate_file <- path.expand(if (length(args) >= 3L) args[[3L]] else "~/covariates_fixed.rds")
model_dir <- path.expand(if (length(args) >= 4L) args[[4L]] else "compiled_models")
results_root <- path.expand(if (length(args) >= 5L) args[[5L]] else "results_true")

# 2. Simulation and sampling settings ---------------------------------------
# Coefficients define the data-generating process, independently of which
# predictors are retained in the fitted model. x4 has no generating effect.
generating_beta <- c(beta0 = 0, beta1 = 1.2, beta2 = 1, beta3 = -1, beta4 = 0)
simulation_seed <- 4569L + task_id
sampling_seed <- 123L + task_id
chains <- 4L
iter_warmup <- 2000L
iter_sampling <- 2000L

# Respect a SLURM CPU allocation; otherwise run the four chains concurrently.
allocated_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
if (is.na(allocated_cpus) || allocated_cpus < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
}
parallel_chains <- min(chains, allocated_cpus)

# 3. Load fixed covariates and the compiled model -----------------------------
# Reusing covariates holds the environmental conditions fixed across replicates.
# Only occupancy states and MCMC sampling vary with task_id.
x <- readRDS(covariate_file)
required_columns <- if (scenario == "misspec1") c("x1", "x2", "x3") else paste0("x", 1:4)
if (!is.data.frame(x) || nrow(x) < 1L || !all(required_columns %in% names(x))) {
  stop("Covariates must be a nonempty data frame containing: ",
       paste(required_columns, collapse = ", "))
}
valid_columns <- vapply(x[required_columns], function(column) {
  is.numeric(column) && all(is.finite(column))
}, logical(1))
if (!all(valid_columns)) stop("Required covariates must be numeric with no missing or infinite values.")

model_file <- file.path(model_dir, paste0("model_true_", scenario, ".rds"))
if (!file.exists(model_file)) stop("Compiled model not found: ", model_file)
model <- readRDS(model_file)

# Create the scenario folder before sampling so output-path errors fail early.
# Reusing the same scenario and task_id replaces that replicate's summary file.
results_dir <- file.path(results_root, scenario, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(results_dir) || file.access(results_dir, 2L) != 0L) {
  stop("Results directory is not writable: ", results_dir)
}

# 4. Simulate true occupancy -------------------------------------------------
# inv_logit(eta), implemented by plogis(), maps the linear predictor to [0, 1].
# Each site then receives one Bernoulli occupancy state with probability psi.
set.seed(simulation_seed)
eta <- generating_beta[["beta0"]] +
  generating_beta[["beta1"]] * x$x1 +
  generating_beta[["beta2"]] * x$x2 +
  generating_beta[["beta3"]] * x$x3
psi <- plogis(eta)
z <- rbinom(nrow(x), size = 1L, prob = psi)

# 5. Build the scenario-specific Stan data -----------------------------------
# Names on the left are Stan inputs; values on the right are source columns.
# Keeping this mapping explicit distinguishes Stan's beta3 in misspec2
# (the coefficient of x4) from beta3 in the data-generating process.
predictor_map <- switch(
  scenario,
  misspec1 = c(x1 = "x1", x2 = "x2"),
  misspec2 = c(x1 = "x1", x2 = "x2", x3 = "x4"),
  misspec3 = c(x1 = "x1", x2 = "x2", x3 = "x3", x4 = "x4")
)
stan_data <- c(
  list(N = nrow(x), y = z),
  setNames(lapply(unname(predictor_map), function(column) x[[column]]), names(predictor_map))
)
parameters <- paste0("beta", 0:length(predictor_map))

# 6. Fit the compiled model --------------------------------------------------
# Warmup adapts the sampler; only the subsequent iterations enter summaries.
# The compiled executable referenced by model must remain accessible.
start_time <- Sys.time()
fit <- model$sample(
  data = stan_data,
  seed = sampling_seed,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  refresh = 0
)
elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# 7. Summarize regression coefficients ---------------------------------------
# Extract only the coefficients fitted by the selected Stan model. Keep
# iteration and chain dimensions separate for R-hat and ESS diagnostics.
summary_parameters <- parameters
draws <- fit$draws(variables = parameters, format = "draws_array")

# Named quantiles provide stable column names in the saved table.
summary_table <- posterior::summarise_draws(
  draws,
  mean = base::mean,
  sd = stats::sd,
  rhat = posterior::rhat,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail,
  ~setNames(
    stats::quantile(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE),
    c("q2.5", "q25", "q50", "q75", "q97.5")
  )
)
summary_table <- summary_table[match(summary_parameters, summary_table$variable), ]

# These reference values belong to the generating model. Under omitted or
# substituted predictors, fitted coefficients need not converge to these
# values: they are not the pseudo-true coefficients of a misspecified model.
# In misspec2, Stan's beta3 refers to x4, whose generating coefficient is zero.
generating_effects <- c(
  intercept = generating_beta[["beta0"]],
  x1 = generating_beta[["beta1"]],
  x2 = generating_beta[["beta2"]],
  x3 = generating_beta[["beta3"]],
  x4 = generating_beta[["beta4"]]
)
predictors <- c("intercept", unname(predictor_map))
reference_values <- unname(generating_effects[predictors])
results_rep <- data.frame(
  replicate = task_id,
  scenario = scenario,
  parameter = summary_parameters,
  predictor = predictors,
  true_value = reference_values,
  post_mean = summary_table$mean,
  post_sd = summary_table$sd,
  q2.5 = summary_table$q2.5,
  q25 = summary_table$q25,
  q50 = summary_table$q50,
  q75 = summary_table$q75,
  q97.5 = summary_table$q97.5,
  rhat = summary_table$rhat,
  ess_bulk = summary_table$ess_bulk,
  ess_tail = summary_table$ess_tail,
  elapsed_seconds = elapsed_seconds
)

# 8. Save one summary table per replicate ------------------------------------
# The attribute records settings without changing the tabular result structure.
# Raw Stan CSV output uses CmdStanR's temporary directory; this file retains
# summaries only. Re-fitting requires the fixed covariates and compiled model.
attr(results_rep, "simulation_settings") <- list(
  generating_beta = generating_beta,
  simulation_seed = simulation_seed,
  sampling_seed = sampling_seed,
  covariate_file = normalizePath(covariate_file),
  model_file = normalizePath(model_file),
  chains = chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  session_info = utils::sessionInfo()
)
result_file <- file.path(results_dir, sprintf("results_%03d.rds", task_id))
saveRDS(results_rep, result_file)
message("Finished ", scenario, ", replicate ", task_id, ": ", result_file)

