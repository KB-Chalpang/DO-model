# Sensitivity to the number of zeros or pseudo-absences (backround data)

This study retains every positive detection record and varies the number of
pseudo-absences used in fitting. Each intercept scenario uses one simulated dataset;
array tasks are zero-count settings, not independent simulation replicates.

## Files and model paths

- `zero_sensitivity_naive.R`: naive logistic regression, default `beta0 = 0`.
- `zero_sensitivity_do.R`: detection-only model, default `beta0 = 0`.
- `submit_zero_sensitivity_naive.sh`: submission script for the naive model.
- `submit_zero_sensitivity_do.sh`: submission script for the DO model.


### Generate fixed covariates

Run this code once before the simulations. Reuse the saved file to keep covariates 
consistent across replicates and model scenarios.

```r
# Set the seed for reproducibility.
set.seed(4567)

# Number of sites.
N <- 150000

# Generate continuous and binary predictors.
x <- data.frame(
  x1 = rnorm(N),
  x2 = rbinom(N, size = 1, prob = 0.5),
  x3 = rnorm(N),
  x4 = rnorm(N)
)

# Save the covariates in your working directory.
saveRDS(x, file = "covariates_fixed4.rds")
```

The predictors `x1`, `x3`, and `x4` follow standard normal distributions.
The predictor `x2` is binary, with probability 0.5 of taking the value 1.

Set the covariate path in each simulation script to the location of this file.


Both R scripts (for naive and do) load `~/covariates_fixed4.rds` with numeric `x1`, `x2`, `x3`,

Set the model paths in Section 2 to your actual files:

| Model | Default model path |
|---|---|
| Naive | `~/model_true_misspec2.rds` |
| DO | `~/mystanmodel_do_misspec2.rds` |

If models are in `compiled_models/`, use that directory instead of `~/`. If
an existing DO object still uses `mystanmodel_po_misspec2.rds`, point to that file
only after confirming it contains the intended model. The referenced compiled
executables must remain accessible. Both fits use the actual x1, x2 and x3.

## Why use q = 0.1?

Here, q is the probability that a site is visited.
Setting q to 0.1 means that approximately 10% of sites are visited
and 90% are unvisited on average.

Unvisited sites have no detection record and are coded as zero. This creates
a large pool of background sites that can be sampled as pseudo-absences. 
A background zero is not a confirmed species absence, it is a pseudo-absence. 
Visited sites with no positive record also enter that pool. The analysis can
therefore vary the amount of background data over a broad range while retaining
all positive records from the same underlying dataset.

Both models receive the same type of pool of background data. The naive logistic model
ignores visitation and misclassification. The DO model includes them through
q, p11 and p10 (and their products theta11 and theta10).

## Simulation and zero selection

Defaults are slopes 1.2, 1 and -1, p11 = 0.80, p10 = 0.05 and two visits.
Only the first misclassified visit is fitted.
Each task resets the simulation seed to 2223.

The sizes of the count of zeros or pseudo-absences are set as follows:

```r
zero_sizes <- sort(c(
  round(seq(1000, length(zero_idx), length.out = 6)),
  5000, 10000, 20000, 40000
))
```

Both sample without replacement and retain all positive records. DO samples are
nested.

## Submit the beta0 = 0 runs

From the directory containing the scripts:

```bash
mkdir -p logs_sens_naive logs_sens_do
sbatch submit_zero_sensitivity_naive.sh
sbatch submit_zero_sensitivity_do.sh
```

Create log directories before submission: SLURM opens log files before running
the shell script. The R scripts create `results_sens_naive/` and `results_sens_do/`
themselves.


The scripts load R/4.4.1. Install cmdstanr, posterior, dplyr and tidyr in an
accessible R library. If needed, add `export R_LIBS_USER="$HOME/.local/R_libs"`
after `module load` in the submission scripts.

## Run the other intercepts

Change beta0 and results_dir in the R script, then update the SLURM job name
and both log directives. Suggested names are:

| Model | beta0 | results_dir | Log directory |
|---|---:|---|---|
| Naive | 0 | `results_sens_naive` | `logs_sens_naive` |
| Naive | 1.5 | `results_sens_naive_beta0_1p5` | `logs_sens_naive_beta0_1p5` |
| Naive | -1.5 | `results_sens_naive_beta0_m1p5` | `logs_sens_naive_beta0_m1p5` |
| DO | 0 | `results_sens_do` | `logs_sens_do` |
| DO | 1.5 | `results_sens_do_beta0_1p5` | `logs_sens_do_beta0_1p5` |
| DO | -1.5 | `results_sens_do_beta0_m1p5` | `logs_sens_do_beta0_m1p5` |

For example, for DO with beta0 = 1.5, set in the R script:

```r
beta0 <- 1.5
results_dir <- "results_sens_do_beta0_1p5"
```

Update the submission script:

```bash
#SBATCH --job-name=sens_N_do_beta0_1p5
#SBATCH --output=logs_sens_do_beta0_1p5/%A_%a.out
#SBATCH --error=logs_sens_do_beta0_1p5/%A_%a.err
```

Then run:

```bash
mkdir -p logs_sens_do_beta0_1p5
sbatch submit_zero_sensitivity_do.sh
```

Changing a job name or log path does not change beta0. Leave the compiled model,
seeds, zero-count grid and initialization unchanged. Keep the array at 1-10.
Wait for existing jobs before editing their shared R script. For simultaneous
intercept runs, make separate script copies and change the Rscript filename in
each submission script accordingly.

## Saved results and validation

Each results directory contains `betas_001.rds` through `betas_010.rds`, with
columns `parameter`, `value`, `zeros`, `task`. Naive files contain beta0-beta3;
DO files also contain q, p11, p10, theta11 and theta10. Re-running a task overwrites its file.
