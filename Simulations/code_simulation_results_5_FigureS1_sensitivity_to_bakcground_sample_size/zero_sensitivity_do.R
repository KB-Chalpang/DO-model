# Sensitivity to the number of zero records: detection-only (DO) model.
# Each array task fits one zero-count setting on the same simulated dataset.
# Submit tasks 1-10 from. Input paths and output names can be changed in Section 2. 
#See README_ZERO_SENSITIVITY.md for other intercept settings and submission examples.

library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)

# 1. Identify the zero-count setting 
# SLURM is checked first; SGE is used when no SLURM array ID is available.
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
if (is.na(task_id)) stop("No array task ID found.")

# 2. Input files, parameters and output paths 
x <- readRDS("~/covariates_fixed4.rds")
model <- readRDS("~/compiled_models/mystanmodel_do_misspec2.rds")

beta0 <- 0
beta1 <- 1.2
beta2 <- 1
beta3 <- -1
p11 <- 0.80
p10 <- 0.05
# q is the probability that a site is visited.
# With q = 0.1, about 10% of sites are visited and 90% are unvisited on average.
# Unvisited sites receive zero records, creating a large background pool for
# varying the number of zeros while retaining every detection record. Visited
# nondetections also enter this pool. A zero is not a confirmed absence: the
# species may occupy an unvisited site or be missed at a visited site.
q <- 0.1
T_visits <- 2L

results_dir <- "results_sens_do"

# 3. Simulate occupancy, visitation and detection records 
# Only the first misclassified detection is used in fitting. 
# Therefore T_visits cane be set equal to one 
simulate_records <- function(x, T_visits, q, beta0, beta1, beta2, beta3, p11, p10) {
  N <- nrow(x)
  eta <- beta0 + beta1 * x$x1 + beta2 * x$x2 + beta3 * x$x3
  psi <- plogis(eta)
  z <- rbinom(N, size = 1L, prob = psi)
  v <- rbinom(N, size = 1L, prob = q)

  # NA indicates an unvisited site.
  first_record <- rep(NA_real_, N)
  for (i in seq_len(N)) {
    if (v[i] == 1L) {
      for (t in seq_len(T_visits)) {
        invisible(rbinom(1L, size = 1L, prob = z[i] * p11))
        if (z[i] == 1L) {
          record <- rbinom(1L, size = 1L, prob = p11)
        } else {
          record <- rbinom(1L, size = 1L, prob = p10)
        }
        if (t == 1L) first_record[i] <- record
      }
    }
  }
  first_record
}

# Fixed simulation seed: tasks vary the zero count, not the underlying data.
set.seed(2223)
first_record <- simulate_records(x, T_visits, q, beta0, beta1, beta2, beta3, p11, p10)

# 4. Define the zero-count settings 
# Visited nondetections and unvisited sites both become background data (pseudo-absences) in the response.
y_all <- first_record
y_all[is.na(y_all)] <- 0
one_idx <- which(y_all == 1)
zero_idx <- which(y_all == 0)

zero_sizes <- sort(c(
  round(seq(1000, length(zero_idx), length.out = 6)),
  5000, 10000, 20000, 40000
))
# Report invalid requests 
if (task_id < 1L || task_id > length(zero_sizes)) {
  stop("Array task ID must be between 1 and ", length(zero_sizes), ".")
}
k <- zero_sizes[task_id]
if (k < 1 || k > length(zero_idx)) stop("Requested zero count exceeds the available pool or is less than 1.")

# 5. Select all positive records and the requested zeros 
# Shuffle all zeros once with the same seed in every task; take the first k.
# This produces nested zero samples within the intercept scenario.
set.seed(2223)
zero_perm <- sample(zero_idx)
keep_zeros <- zero_perm[1:k]
keep_idx <- c(one_idx, keep_zeros)
stan_data <- list(
  N = length(keep_idx), y = y_all[keep_idx],
  x1 = x$x1[keep_idx], x2 = x$x2[keep_idx], x3 = x$x3[keep_idx]
)

# 6. Fit the compiled model
# Fixed starts satisfy theta10 < theta11 < q. Stan derives p11 and p10
# from theta11/q and theta10/q; they are not initialized separately.
init_fun <- function() {
  list(
    q = 0.6,
    theta11 = 0.3,
    theta10 = 0.01,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    beta3 = 0
  )
}
fit <- model$sample(
  data = stan_data,
  seed = 1000 + task_id,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  init = init_fun,
  refresh = 0
)

# 7. Save posterior draws 
# One row per parameter per posterior draw. Only post-warmup draws are saved.
# Column names are parameter, value, zeros and task for downstream plotting.
parameters <- c("beta0", "beta1", "beta2", "beta3", "q", "p11", "p10", "theta11", "theta10")
beta_df <- as.data.frame(as_draws_df(fit$draws(variables = parameters))) |>
  select(all_of(parameters)) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(zeros = k, task = task_id)

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(beta_df, file = file.path(results_dir, sprintf("betas_%03d.rds", task_id)))
cat("Finished task", task_id, "with", k, "zeros\n")
