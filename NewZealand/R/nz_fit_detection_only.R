# New Zealand case study: detection-only model (DO-M)
#
# Fits the detection-only model to the presence-absence records of one species,
# treated as detection-only data, and saves the objects used by
# nz_case_study.Rmd:
# y_obs_do_<species>.rds (observed 0/1 records)
# yrep_do_<species>.rds (posterior predictive draws (draws x sites))
# results_rep_po_do1_<species>.rds  (posterior summary of the parameters)
#
# all these files are provided in the Results folder, since 
# the sample size is large and it will take long time to run.
#
# Usage: run from the project folder (the one containing R/ and stan/). Set `species_id`
# below ("nz35" or "nz42"), run the script, then repeat for the other species.
#
# Requires: CmdStan (via cmdstanr), stan/do_model_nz.stan, R/prepare_nz_data.R

library(cmdstanr)
library(posterior)
library(dplyr)

path <- "insert path"
source(file.path(path, "prepare_nz_data.R"))

# Settings 
species_id <- "nz42"       # species to fit: "nz35" or "nz42"
out_dir    <- path    # folder where the results are saved

stan_file  <- file.path(path, "do_model_nz.stan")
seed       <- 123

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)       # seeds the random initial values

# Data 
nz <- prepare_nz_data(species_id)

stan_data <- list(N = nrow(nz$x), K = ncol(nz$x), y = nz$y, x = nz$x)


# Fit the model 
model <- cmdstan_model(stan_file)

init_fun <- function() {
  list(q = runif(1, 0.5, 0.8),
       theta11 = runif(1, 0.15, 0.45), theta10 = runif(1, 0.001, 0.1),
       beta0 = rnorm(1), beta = rnorm(stan_data$K))
}

start_time <- Sys.time()
fit <- model$sample(
  data = stan_data,
  seed = seed,
  chains = 4, parallel_chains = 4,
  iter_warmup = 2000, iter_sampling = 2000,
  max_treedepth = 13, adapt_delta = 0.99,
  refresh = 500,
  init = init_fun
)
computational_time <- Sys.time() - start_time
print(computational_time)


# Posterior summary 
# Parameter names follow the Rmd: beta0, beta1, ..., beta11, theta11, theta10
results_rep_po_do1 <- fit$summary(
  variables = c("beta0", "beta", "theta11", "theta10"),
  mean, sd,
  ~quantile2(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
  rhat, ess_bulk, ess_tail
) |>
  mutate(parameter = sub("\\[(\\d+)\\]", "\\1", variable)) |> 
  select(parameter, post_mean = mean, post_sd = sd,
         q2.5, q25, q50, q75, q97.5, rhat, ess_bulk, ess_tail)

print(results_rep_po_do1, n = Inf)


# Objects for the posterior predictive checks 
y_obs_do <- nz$y

# Stored as an integer matrix (draws x sites) to halve the file size
yrep_draws <- fit$draws("ypred", format = "matrix")
yrep_do    <- matrix(as.integer(yrep_draws), nrow = nrow(yrep_draws))
rm(yrep_draws)

# Save 
saveRDS(y_obs_do,           file.path(out_dir, sprintf("y_obs_do_%s.rds",            species_id)))
saveRDS(yrep_do,            file.path(out_dir, sprintf("yrep_do_%s.rds",             species_id)))
saveRDS(results_rep_po_do1, file.path(out_dir, sprintf("results_rep_po_do1_%s.rds",  species_id)))
