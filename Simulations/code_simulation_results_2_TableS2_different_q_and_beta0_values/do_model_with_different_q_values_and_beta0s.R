# This is the code used for the detection-only model (DO-M)
# with different levels of q. The code is run on an HPC cluster as a SLURM job
# array: each array task is one simulation replicate (see submit_ajob_do.sh).
#
# Usage: Rscript submit_ajob_do.R <task_id> <q> <beta0> <results_dir>
#   task_id      simulation replicate (1, 2, ..., 500)
#   q            probability that a site is visited (default 0.5)
#   beta0        intercept of the occupancy model (default 0)
#   results_dir  folder for the results (default: named after q and beta0)

# Install and load the following libraries
library(cmdstanr)
library(dplyr)       # provides tibble()
library(posterior)   # provides summarise_draws()

# create the data frame covariates_fixed, so that the environmental conditions
# do not change for the different replications. Therefore, the code below should be
# run separately, once, before the simulations.
#
# set.seed(4567)
#N <- 5000
#x <- data.frame(
#  x1 = rnorm(N),
#  x2 = rbinom(N, 1, 0.5),
#  x3 = rnorm(N),
#  x4 = rnorm(N)
#)
#saveRDS(x, "covariates_fixed.rds")
#
# load the covariates_fixed data frame
x <- readRDS("~/covariates_fixed.rds")
x1 <- x$x1
x2 <- x$x2
x3 <- x$x3

# The task id, q and beta0 are passed on by the submission script
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1]) # each task_id is a simulation replicate or run
q       <- if (length(args) >= 2) as.numeric(args[2]) else 0.5  # probability that a site is visited
beta0   <- if (length(args) >= 3) as.numeric(args[3]) else 0    # intercept

# A different seed for each replicate, so that replicates are independent
# but every run can be reproduced exactly
set.seed(4567 + task_id)

# ------------ parameters --------------
N <- nrow(x)
# three intercept scenarios are considered in the simulation, beta0 = 0,
# beta0 = 1.5 and beta0 = -1.5; the value of beta0 is the third argument of the script
beta1 <- 1.5
beta2 <- 1
beta3 <- -1
p11 <- .85  # Pr(y = 1 | z = 1): probability of recording an occupied site
p10 <- .05  # Pr(y = 1 | z = 0): probability of a false positive at an unoccupied site
# run simulation for each value of q for {0.5, 0.6, 0.7, 0.8, 0.9, 1};
# the value of q is the second argument of the script
T_visits <- 2 # number of visits

# save and compile the stan model separately, so that each replicate or run
# will use the object mystanmodel_do.rds so the stan model will not be compiled
# in each replicate
#
#mystanmodel_do <- readRDS("~/mystanmodel_do.rds")
mystanmodel_do <- readRDS("~/mystanmodel_do.rds")

# data simulation function
# Simulates, for N sites and T_visits visits per site:
#   z      latent occupancy state, z ~ Bernoulli(psi), psi = plogis(beta0 + beta1*x1 + beta2*x2 + beta3*x3)
#   v      whether the site is visited, v ~ Bernoulli(q); unvisited sites have no records (NA)
#   y_true detections without misclassification (only occupied sites can be detected)
#   y_miss detections with misclassification (occupied sites: detected with probability p11;
#          unoccupied sites: false positive with probability p10)
# The random numbers are drawn in a fixed order (z, then v, then site by site), so the
# same seed always gives the same data.
sim <- function(N, T_visits, q,
                beta0, beta1, beta2, beta3,
                p11, p10, x1, x2, x3) {
  
  # True occupancy probability
  eta <- beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3
  psi <- plogis(eta)
  # Latent occupancy state (1 = occupied, 0 = unoccupied)
  z <- rbinom(N, 1, psi)
  
  # Visit process (whether sites are visited or not)
  v <- rbinom(N, 1, q)
  
  # Simulate detection data: one row per site, one column per visit.
  # Entries stay NA for sites that are not visited.
  y_true <- matrix(NA_integer_, nrow = N, ncol = T_visits,
                   dimnames = list(NULL, paste0("y_true_", 1:T_visits)))
  y_miss <- matrix(NA_integer_, nrow = N, ncol = T_visits,
                   dimnames = list(NULL, paste0("y_miss_", 1:T_visits)))
  for (i in which(v == 1)) {          # only visited sites
    for (t in 1:T_visits) {
      # True detection
      y_true[i, t] <- rbinom(1, 1, z[i] * p11)
      # Misclassification: true positive if occupied (p11), false positive if not (p10)
      y_miss[i, t] <- rbinom(1, 1, if (z[i] == 1) p11 else p10)
    }
  }
  
  # Combine into main dataset
  # (y_binom_* aggregate the detections over visits; NAs of unvisited sites count as 0)
  data <- data.frame(
    site = 1:N,
    z = z,
    v = v,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    y_binom_true = rowSums(y_true, na.rm = TRUE),
    y_binom_miss = rowSums(y_miss, na.rm = TRUE),
    T_visits_actual = v * T_visits,   # number of visits actually made (T_visits or 0)
    y_true, y_miss
  )
  
  # Detection-only data: extend to all T visits
  # The detection-only records are the misclassified detections, renamed y_do_1, y_do_2, ...
  y_do <- y_miss
  colnames(y_do) <- paste0("y_do_", 1:T_visits)
  # version with NA's
  data_do_NAs <- cbind(data[c("site", "z", "v", "x1", "x2", "x3")], y_do)
  
  # version where NA --> 0 (pseudo-absence)
  data_do_0s <- data_do_NAs
  data_do_0s[is.na(data_do_0s)] <- 0
  
  list(
    data = data,
    # Subset for occupancy model (visited sites only)
    data_occ_true = subset(data, v == 1),
    data_do_NAs = data_do_NAs,
    data_do_0s = data_do_0s,
    betas = c(beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3),
    observation_probs = c(p11, p10),
    true_psi = psi,
    true_eta = eta
  )
}

# simulate data
simdata <- sim(N, T_visits, q,
               beta0, beta1, beta2, beta3,
               p11, p10,
               x1, x2, x3)

# Stan data for detection-only model with correction for misclassification
# (the DO-M is fitted to the record of the first visit, y_do_1, with 0 for unvisited sites)
stan_data_do <- list(
  N = nrow(simdata$data_do_0s),
  y = simdata$data_do_0s$y_do_1,
  x1 = simdata$data_do_0s$x1,
  x2 = simdata$data_do_0s$x2,
  x3 = simdata$data_do_0s$x3
)

# initialisation (random starting values for each chain)
init_fun_do <- function() {
  list(
    q = runif(1, 0.5, 0.8),
    theta11 = runif(1, 0.15, 0.4), theta10 = runif(1, 0.01, 0.1),
    beta0 = rnorm(1, 0, 1), beta1 = rnorm(1, 0, 1), beta2 = rnorm(1, 0, 1),
    beta3 = rnorm(1, 0, 1)
  )
}

# fitting the stan model
start.time <- Sys.time()
fit <- mystanmodel_do$sample(
  data = stan_data_do,
  seed = 123 + task_id,
  chains = 4,
  parallel_chains = 4,   # one core per chain (matches --cpus-per-task in the job script)
  iter_warmup = 2000,
  iter_sampling = 2000,
  refresh = 0,
  init = init_fun_do
)
end.time <- Sys.time()
computational.time <- end.time - start.time
computational.time

# Parameters of interest
params <- c("beta0", "beta1", "beta2", "beta3", "p11", "p10", "theta11", "theta10")
# vector of true values; theta = q * p is the probability of a record
# (a site is visited and the species is recorded)
true_vals <- c(simdata$betas,
               p11 = p11, p10 = p10,
               theta11 = q * p11, theta10 = q * p10)

# Summarise posterior draws with desired quantiles
quantiles <- summarise_draws(
  fit,
  ~quantile(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
)

# Compute posterior mean, sd, and the diagnostics rhat, ess_bulk, ess_tail
sum_df <- summarise_draws(fit, mean, sd, rhat, ess_bulk, ess_tail)

# Construct a tidy table for this replicate
results_rep <- tibble(
  replicate  = task_id,
  parameter  = params,
  true_value = unname(true_vals[match(params, names(true_vals))]),
  post_mean  = sum_df$mean[match(params, sum_df$variable)],
  post_sd    = sum_df$sd[match(params, sum_df$variable)],
  q2.5       = quantiles$`2.5%`[match(params, quantiles$variable)],
  q25        = quantiles$`25%`[match(params, quantiles$variable)],
  q50        = quantiles$`50%`[match(params, quantiles$variable)],
  q75        = quantiles$`75%`[match(params, quantiles$variable)],
  q97.5      = quantiles$`97.5%`[match(params, quantiles$variable)],
  rhat       = sum_df$rhat[match(params, sum_df$variable)],
  ess_bulk   = sum_df$ess_bulk[match(params, sum_df$variable)],
  ess_tail   = sum_df$ess_tail[match(params, sum_df$variable)]
)
results_rep
# save results
# one folder per scenario, named after q and beta0, for example results_do05_beta0_0
# (q = 0.5, beta0 = 0), results_do06_beta0_m1.5 (q = 0.6, beta0 = -1.5; "m" stands for minus)
results_dir <- if (length(args) >= 4) args[4] else
  sprintf("results_do%02d_beta0_%s", round(q * 10), gsub("-", "m", as.character(beta0)))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Save each replicate separately
saveRDS(results_rep,
        file = sprintf("%s/post_means_%03d.rds", results_dir, task_id))

cat("Finished replicate", task_id, "\n")