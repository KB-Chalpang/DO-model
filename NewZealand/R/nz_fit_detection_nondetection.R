# New Zealand case study: detection/non-detection model with one visit (DN-1-M)
#
# Fits the occupancy model with misclassification (sensitivity p11 and
# false-positive rate p10) to the presence-absence records of one species, and
# saves the posterior summary used by nz_case_study.Rmd:
# results_rep_occ_<species>.rds 
#
# all these files are provided in the Results folder, since 
# the sample size is large and it will take long time to run.
#
# Usage: run from the project folder (the one containing R/ and stan/). Set `species_id`
# below ("nz35" or "nz42"), run the script, then repeat for the other species.
#
# Requires: CmdStan (via cmdstanr), stan/occupancy_model_nz.stan, R/prepare_nz_data.R

library(cmdstanr)
library(posterior)
library(dplyr)

path <- "C:/Users/ahw881/Downloads/nz_files/R"
source(file.path(path, "prepare_nz_data.R"))


# Settings 
species_id <- "nz35"       # species to fit: "nz35" or "nz42"
out_dir    <- path    # folder where the results are saved

stan_file  <- file.path(path, "occupancy_model_nz.stan")
T_visits   <- 1L     # one visit per site
seed       <- 123

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Data 
nz <- prepare_nz_data(species_id)

stan_data <- list(N = nrow(nz$x), K = ncol(nz$x), T_visits = T_visits,
                  y = as.matrix(nz$y),   # number of detections per site (0/1 for one visit)
                  x = nz$x)

# Fit the model 
model <- cmdstan_model(stan_file)

init_fun <- function() {
  list(beta0 = 0, beta = rep(0, stan_data$K), p11 = 0.01, p10 = 0.01)
}

start_time <- Sys.time()
fit <- model$sample(
  data = stan_data,
  seed = seed,
  chains = 4, parallel_chains = 4,
  iter_warmup = 2000, iter_sampling = 2000,
  refresh = 500,
  init = init_fun
)
computational_time <- Sys.time() - start_time
print(computational_time)


#  Posterior summary 
# Parameter names follow the Rmd: beta0, beta1, ..., beta11, p11, p10
results_rep_occ <- fit$summary(
  variables = c("beta0", "beta", "p11", "p10"),
  mean, sd,
  ~quantile2(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
  rhat, ess_bulk, ess_tail
) |>
  mutate(parameter = sub("\\[(\\d+)\\]", "\\1", variable)) |>   # "beta[3]" -> "beta3"
  select(parameter, post_mean = mean, post_sd = sd,
         q2.5, q25, q50, q75, q97.5, rhat, ess_bulk, ess_tail)

print(results_rep_occ, n = Inf)

# Save 
saveRDS(results_rep_occ, file.path(out_dir, sprintf("results_rep_occ_%s.rds", species_id)))
