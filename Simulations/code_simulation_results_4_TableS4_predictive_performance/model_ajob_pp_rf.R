# Five-fold cross-validation for random forest.
# Usage: Rscript model_ajob_pp_rf.R [q]
# q may be 0.5 or 1; it defaults to the value in Section 2.
# One SLURM/SGE array task is one simulation replicate. Outside an array, the
# replicate defaults to 1. All model scripts use matching data and outer folds
# for identical replicate IDs, covariates, simulation settings, etc.

# 1. Replicate and dependencies 
task_text <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (!nzchar(task_text)) task_text <- Sys.getenv("SGE_TASK_ID", "")
if (!nzchar(task_text)) task_text <- "1"
task_number <- suppressWarnings(as.numeric(task_text))
if (!is.finite(task_number) || task_number < 1 || task_number != floor(task_number) ||
    task_number > .Machine$integer.max - 1000) stop("Invalid replicate ID.")
task_id <- as.integer(task_number)
set.seed(1000 + task_id)
library(pROC)
library(randomForest)

# 2. Inputs, simulation settings and output directory 
q <- 1
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript model_ajob_pp_rf.R [0.5|1]")
if (length(args)) q <- suppressWarnings(as.numeric(args[[1L]]))
if (!is.finite(q) || !q %in% c(0.5, 1)) stop("q must be 0.5 or 1.")

covariate_file <- "~/covariates_fixed2.rds"
x <- readRDS(covariate_file)
predictors <- c("x1", "x2", "x3")
if (!is.data.frame(x) || !all(predictors %in% names(x))) {
  stop("Covariates must be a data frame with x1, x2 and x3.")
}
if (!all(vapply(x[predictors], function(a) is.numeric(a) && all(is.finite(a)), logical(1)))) {
  stop("Predictors must be numeric without missing or infinite values.")
}

beta0 <- 0
beta1 <- 1.2
beta2 <- 1
beta3 <- -1
p11 <- 0.80
p10 <- 0.05
T_visits <- 3L
K <- 5L
if (nrow(x) < K) stop("At least K sites are needed for K-fold CV.")
q_tag <- if (q == 0.5) "05" else "1"
results_dir <- paste0("summaries_pp_rf_q", q_tag)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
cat("Running random forest replicate", task_id, "with q =", q, "\n")

# 3. Simulate occupancy and observed records 
# q controls site visitation. Unvisited sites become zero records; occupied
# sites can be missed, and unoccupied sites can generate false positives.
# Retain the order of all detection draws across two visits, even though
# only the first misclassified record is fitted.

simulate_data <- function(x, T_visits, q, beta0, beta1, beta2, beta3, p11, p10) {
  N <- nrow(x)
  psi <- plogis(beta0 + beta1 * x$x1 + beta2 * x$x2 + beta3 * x$x3)
  z <- rbinom(N, 1, psi)
  v <- rbinom(N, 1, q)
  first_record <- rep(NA_real_, N)
  for (i in seq_len(N)) {
    if (v[i] == 1L) {
      for (t in seq_len(T_visits)) {
        invisible(rbinom(1L, 1, z[i] * p11))
        if (z[i] == 1L) {
          record <- rbinom(1L, 1, p11)
        } else {
          record <- rbinom(1L, 1, p10)
        }
        if (t == 1L) first_record[i] <- record
      }
    }
  }
  first_record[is.na(first_record)] <- 0
  data.frame(z = z, y = first_record, x1 = x$x1, x2 = x$x2, x3 = x$x3)
}
data_cv <- simulate_data(x, T_visits, q, beta0, beta1, beta2, beta3, p11, p10)

# 4. Assign outer cross-validation folds 
# Balanced fold sizes are randomized with the same fixed seed in every script.
# Training uses observed records y; evaluation uses held-out true occupancy z.
set.seed(123)
data_cv$fold <- sample(rep(1:K, length.out = nrow(data_cv)))
cv_results <- vector("list", K)

# 5. Fit and evaluate each outer fold 

start_time <- Sys.time()
for (k in seq_len(K)) {
  cat("Fold", k, "\n")
  train_data <- data_cv[data_cv$fold != k, , drop = FALSE]
  test_data <- data_cv[data_cv$fold == k, , drop = FALSE]
  # A factor response selects classification rather than regression.
  train_data$y <- as.factor(train_data$y)
  set.seed(346)
  model_fit <- randomForest::randomForest(
    y ~ x1 + x2 + x3, data = train_data, ntree = 500
  )
  p_hat <- predict(model_fit, test_data, type = "prob")[, "1"]
  # Scores use true occupancy z, while models are trained on observed y.
  # Clipping avoids log(0). AUC retains pROC's automatic direction selection.
  p_hat <- pmin(pmax(p_hat, 1e-6), 1 - 1e-6)
  y_test <- test_data$z
  cv_results[[k]] <- list(
    elpd = sum(y_test * log(p_hat) + (1 - y_test) * log(1 - p_hat)),
    brier = mean((y_test - p_hat)^2),
    auc = as.numeric(pROC::auc(response = y_test, predictor = p_hat, quiet = TRUE))
  )
}
elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# 6. Summarize fold metrics 
# Each fold contributes one score per metric. Bounds are the 2.5% and
# 97.5% quantiles of the five fold scores; they are not a CV confidence interval.
metric_names <- c("elpd", "brier", "auc")
summary_results <- data.frame(
  Metric = c("ELPD", "Brier", "AUC"),
  Mean = vapply(metric_names, function(metric) {
    mean(vapply(cv_results, function(result) result[[metric]], numeric(1)))
  }, numeric(1)),
  CI_lower = vapply(metric_names, function(metric) {
    quantile(vapply(cv_results, function(result) result[[metric]], numeric(1)), 0.025)
  }, numeric(1)),
  CI_upper = vapply(metric_names, function(metric) {
    quantile(vapply(cv_results, function(result) result[[metric]], numeric(1)), 0.975)
  }, numeric(1)),
  row.names = NULL
)

# 7. Save the summary and the fold-level values used to construct it 
# Summary column names remain Metric, Mean, CI_lower and CI_upper. See the
# README for their exact calculation; the labels do not define a new interval.
saveRDS(summary_results, file.path(results_dir, paste0("summary_rf_rep_", task_id, ".rds")))
saveRDS(cv_results, file.path(results_dir, paste0("fold_metrics_rf_rep_", task_id, ".rds")))
saveRDS(list(
  model = "rf", replicate = task_id, q = q, beta0 = beta0,
  beta1 = beta1, beta2 = beta2, beta3 = beta3, p11 = p11, p10 = p10,
  T_visits = T_visits, K = K, simulation_seed = 1000 + task_id,
  fold_seed = 123, fold_assignment = data_cv$fold,
  covariate_file = path.expand(covariate_file),
  elapsed_seconds = elapsed_seconds, session_info = utils::sessionInfo()
), file.path(results_dir, paste0("settings_rf_rep_", task_id, ".rds")))
cat("Finished random forest replicate", task_id, "in", elapsed_seconds, "seconds\n")
