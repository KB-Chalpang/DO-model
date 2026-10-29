# Cross-validation of GLM, Maxent, BRT, RF and DO models

This bundle contains five independent R scripts, five SLURM submission scripts
and one launcher, `submit_all_pp.sh`, that submits the model arrays together.
Each R script simulates one replicate, assigns five outer folds, fits its model
five times and evaluates the held-out predictions. The two visitation scenarios
are q = 0.5 and q = 1; they use separate logs and output directories.

## Files

| Method | R script | SLURM script |
|---|---|---|
| Weighted GLM | `model_ajob_pp_glm.R` | `submit_pp_glm.sh` |
| Maxent | `model_ajob_pp_maxent.R` | `submit_pp_maxent.sh` |
| Boosted regression trees | `model_ajob_pp_brt.R` | `submit_pp_brt.sh` |
| Random forest | `model_ajob_pp_rf.R` | `submit_pp_rf.sh` |
| Detection-only | `model_ajob_pp_do.R` | `submit_pp_do.sh` |

Each R script is self-contained. No shared R helper, task table or preparation
stage is required.

## Prerequisites and settings

Create the fixed covariates once before submitting jobs:

```r
# Set the seed for reproducibility.
set.seed(4567)

# Number of sites.
N <- 15000

# Generate continuous and binary predictors.
x <- data.frame(
  x1 = rnorm(N),
  x2 = rbinom(N, size = 1, prob = 0.5),
  x3 = rnorm(N)
)

# Save in the location read by the simulation scripts.
saveRDS(x, file = "~/covariates_fixed2.rds")
```

All scripts read `~/covariates_fixed2.rds`, containing numeric x1, x2 and x3.
The DO script also reads `~/mystanmodel_do_misspec2.rds`. Edit these paths in
Section 2 to match your files.

The default generating values are beta0 = 0, beta1 = 1.2, beta2 = 1, beta3 = -1,
p11 = 0.80 and p10 = 0.05. All model formulas use the actual x1, x2 and x3.

Package requirements:

| Method | R packages |
|---|---|
| GLM | pROC |
| Maxent | pROC, maxnet |
| BRT | pROC, dismo, gbm |
| RF | pROC, randomForest |
| DO | cmdstanr, pROC; configured CmdStan environment |

Submission scripts load R/4.4.1. If required by your installation, add
`export R_LIBS_USER="$HOME/.local/R_libs"` after the module command.

## Submit with SLURM

Place `submit_all_pp.sh` alongside the five R scripts and five model-specific
SLURM scripts. Run the launcher on your cluster's login/submission node using
**bash**, because it calls `sbatch` for each model itself.

Submit all five models for **both q values** with one command:

```bash
bash submit_all_pp.sh
```

Alternatively, submit all five models for just one q value:

```bash
bash submit_all_pp.sh 0.5
```

```bash
bash submit_all_pp.sh 1
```

The default (or the argument `all`) submits ten job arrays: five models for
each of two q values. Selecting one q submits five job arrays. Each array
contains 500 replicate tasks, so these commands submit 5,000 or 2,500 tasks,
respectively. Each task fits five outer folds. SLURM schedules the jobs as
resources become available; submitting them together does not require them
to start simultaneously.

For a small initial run across all models and both q values:

```bash
ARRAY=1-2 bash submit_all_pp.sh
```

For a small initial run at q = 0.5 only:

```bash
ARRAY=1-2 bash submit_all_pp.sh 0.5
```

The launcher prints the model, q value and `sbatch` submission response,
including the job ID. It checks that all required scripts exist before
submitting anything. If a submission fails, it stops; earlier submitted
jobs remain queued or running. Review the printed job IDs before retrying
so that you do not submit the same arrays twice.

The launcher submits from its own directory, where logs and results will be
created. No log directories need to be created beforehand: SLURM initially
writes to that directory, and each model script creates and redirects to its
model/q-specific log folder before running R.

Alos, individual submission remains available from the scripts directory:

```bash
sbatch submit_pp_glm.sh 0.5
sbatch --array=1-2 submit_pp_do.sh 1
```

Omitting q from an individual model submission selects q = 1. Omitting q from
the all-model launcher selects **both** q values. The launcher uses each
model's existing resource requests, R script and output directories.


## Simulation and reproducibility

The simulation seed is 1000 + replicate ID. Occupancy is Bernoulli with
logistic probability based on x1, x2 and x3. Visitation is Bernoulli(q).
For each visited site, positive records have probability p11 if occupied and
p10 if unoccupied. Unvisited sites are assigned pseudo-absence records.


The fold seed is 123. The fold assignment samples a repeated sequence 1:5,
giving balanced fold sizes, without class stratification. Matching input data,
replicate IDs, settings and RNG (Random Number Generator) behavior give matching simulation and fold
assignments across methods. Training uses observed records y; evaluation uses
true occupancy z for every method.


## Models and scores as implemented

- GLM uses a binomial logit model with positive weight 1 and pseudo-absences weight
  n_positive/n_background, calculated within the training fold. This equalizes
  total positive and background weight.
- Maxent uses maxnet's default feature classes, regularization multiplier 1,
  and cloglog predictions. Other maxnet options retain their defaults.
- BRT uses the same positive/background weights, tree complexity 1 when fewer
  than 50 training positives are present and 5 otherwise, learning rate 0.001,
  and bag fraction 0.75. gbm.step performs internal cross-validation within
  the outer training fold to select its tree count.
- RF uses a factor response, 500 trees and class-1 probabilities. No weighting
  is supplied.
- DO uses four chains, 2,000 warmup and 2,000 sampling iterations.

All predictions are clipped to [0.000001, 0.999999] for score calculations.
AUC retains pROC's automatic direction selection. No direction or class-level
override has been added.

The GLM, Maxent, BRT and RF scripts compute, for each test fold:

- `ELPD`: the sum of Bernoulli log scores evaluated against z.
- `Brier`: the mean squared difference between the prediction and z.
- `AUC`: the ROC area against z using the package's automatic direction.

For DO, each posterior draw produces psi from its regression coefficients,
then `p_obs = theta11 * psi + theta10 * (1 - psi)`. These record probabilities
are evaluated against z. 


## Summary columns

The summary schema remains `Metric`, `Mean`, `CI_lower`, `CI_upper` with rows
ELPD, Brier and AUC.

For GLM, Maxent, BRT and RF, Mean averages the five fold scores. CI_lower and
CI_upper are the 2.5% and 97.5% quantiles of those five scores. For DO, Mean
averages the five within-fold posterior score means; each bound averages the
corresponding within-fold posterior quantile.

ELPD uses fold sums and then averages across folds. 

## Outputs and changing scenarios

The q argument determines both output paths and the value passed to R:

| q | Log directory | Results directory |
|---|---|---|
| 0.5 | `logs_pp_<model>_q05/` | `summaries_pp_<model>_q05/` |
| 1 | `logs_pp_<model>_q1/` | `summaries_pp_<model>_q1/` |

Replace `<model>` with glm, maxent, brt, rf or do.
RF output follows the actual q setting, including q1 when q = 1.

Each replicate writes:

- The summary RDS: `summary_<model>_rep_<id>.rds`.
- `fold_metrics_<model>_rep_<id>.rds`: the fold-score list. DO contains vectors
  of draw-wise scores within each fold; other methods contain scalar scores.
- `settings_<model>_rep_<id>.rds`: generating settings, fold assignment, elapsed
  time, input covariate path and R session information.

When combining summaries, select only the `summary_` filenames rather than
reading every RDS file in the directory. Re-running the same model/q/replicate
overwrites its RDS outputs. Log filenames include the job ID to distinguish
submissions. Changing beta0 or any other generating setting requires a new
results path or copied study directory to avoid overwriting earlier results.

Startup logs are `pp_<model>_setup_<job>_<task>.out/.err` in the submission
directory. After redirection, logs are written under the matching model/q folder.
If a run fails, inspect both locations. Degenerate folds or model-fitting
failures are not silently skipped, resampled or replaced with another method.
