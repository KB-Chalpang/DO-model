# Prior sensitivity for the detection-only model

This study compares six specifications for the priors on p11 and p10:
one baseline and five alternatives. The data-generating parameters remain
fixed; only the fitted model's prior specification changes.

## Files

- `compile_prior_sensitivity_models.R`: writes and compiles the six Stan models.
- `submit_compile_prior_sensitivity.sh`: runs compilation as one SLURM job.
- `run_do_prior_sensitivity.R`: simulates one replicate, fits one scenario and saves posterior summaries.
- `submit_do_prior_sensitivity.sh`: submits all six scenarios as one SLURM array.

## Prior scenarios

The mean of a Beta(alpha, beta) distribution is alpha / (alpha + beta).
Here, "center" refers to this mean, not necessarily the mode.

| Scenario | Prior for p11 | Mean | Prior for p10 | Mean |
|---|---|---:|---|---:|
| a: baseline | Beta(8, 2) | 0.80 | Beta(1, 19) | 0.05 |
| b | Beta(4, 1) | 0.80 | Beta(0.5, 9.5) | 0.05 |
| c | Beta(6, 4) | 0.60 | Beta(2, 18) | 0.10 |
| d | Beta(18, 2) | 0.90 | Beta(1, 39) | 0.03 |
| e | Beta(2, 2) | 0.50 | Beta(1, 9) | 0.10 |
| f | Beta(2, 3) | 0.40 | Beta(2, 5) | 0.29 |

Means in the table are rounded to two decimal places. The exact means for
p10 in scenarios d and f are 0.025 and 2/7, respectively.


Scenario b retains the baseline means with less concentration. The other
alternatives also change the prior means. This is to assess whether effect estimates
are sensitive to to the different prior specifications.

## Prepare the compiled models

Run compilation once from the directory containing all five scripts:

```bash
sbatch submit_compile_prior_sensitivity.sh
```

The compiler creates **`prior_sensitivity_compiled_models/`** and writes six
uniquely named `.stan` sources, their compiled executables and six `.rds`
CmdStanR model objects:

| Scenario | Compiled model object |
|---|---|
| a | `prior_sensitivity_compiled_models/model_do_prior_a.rds` |
| b | `prior_sensitivity_compiled_models/model_do_prior_b.rds` |
| c | `prior_sensitivity_compiled_models/model_do_prior_c.rds` |
| d | `prior_sensitivity_compiled_models/model_do_prior_d.rds` |
| e | `prior_sensitivity_compiled_models/model_do_prior_e.rds` |
| f | `prior_sensitivity_compiled_models/model_do_prior_f.rds` |

It also saves `prior_specifications.csv` with the Beta shape parameters.
The fitting script reads these model paths automatically.

The shared Stan template is inside `compile_prior_sensitivity_models.R`.
Only the four Beta shape parameters vary between scenarios.

Keep the directory and executables accessible at the same paths on compute
nodes: the RDS objects store executable paths, not self-contained executables.
Compile in the final project location. Recompile if that location changes.
Do not recompile while sampling jobs are using these executables. A compilation
rerun overwrites the scenario sources, executables and RDS objects.

Compilation is sequential, requesting one CPU, 8 GB total memory and one hour.
Adjust these requests for your cluster. Compilation logs are
`compile_do_priors_<job>.out/.err` in the submission directory. There is no
pre-submission log directory to create. On a node where compilation is allowed,
you may instead run `Rscript compile_prior_sensitivity_models.R` directly.

## Covariates and software

The R script reads `~/covariates_fixed.rds`: a nonempty data frame with finite,
numeric columns x1, x2 and x3. Extra columns are ignored. 
Every scenario uses the same covariates.

Install cmdstanr and posterior and configure CmdStan on the cluster.
The submission script loads `R/4.4.1`; adjust that module for your cluster.
If needed, add `export R_LIBS_USER="$HOME/.local/R_libs"` after loading R.

## Submit all scenarios

After compilation completes successfully, run from the same directory:

```bash
sbatch submit_do_prior_sensitivity.sh
```

You can also submit compilation and fitting together with an SLURM dependency:

```bash
compile_job=$(sbatch --parsable submit_compile_prior_sensitivity.sh) &&
  sbatch --dependency=afterok:"${compile_job%%;*}" submit_do_prior_sensitivity.sh
```

The fitting array starts only after compilation exits successfully. If
compilation fails, inspect its logs and submit a new compilation job; a fresh
fitting submission must depend on the new job ID. Use either this workflow or
the separate submissions above to avoid duplicate runs.

The array has 3,000 tasks: six scenarios with 500 replicates each.

| Array IDs | Scenario | Replicates |
|---|---|---|
| 1-500 | a | 1-500 |
| 501-1000 | b | 1-500 |
| 1001-1500 | c | 1-500 |
| 1501-2000 | d | 1-500 |
| 2001-2500 | e | 1-500 |
| 2501-3000 | f | 1-500 |

To test replicate 1 of every scenario:

```bash
sbatch --array=1,501,1001,1501,2001,2501 submit_do_prior_sensitivity.sh
```

To submit only scenario b:

```bash
sbatch --array=501-1000 submit_do_prior_sensitivity.sh
```

Resource requests retain the supplied settings: 30 minutes, four CPUs and
3 GB per CPU (12 GB per task). Adjust these for your cluster and observed
runtimes. If the cluster limits array size below 3,000, submit the six ranges
above separately. To change the number of replicates per scenario, change both
`replicates=500` and the `#SBATCH --array` upper bound (six times that number).

No log directory must be created beforehand. Initial SLURM output is written
to `do_prior_setup_<array>_<task>.out/.err` in the submission directory. Once
the script creates its scenario log directory, output is redirected there.
Inspect the startup logs if a job stops before its scenario logs appear.

For a single fit on a node where R execution is permitted:

```bash
Rscript run_do_prior_sensitivity.R a 1
```

## Simulation and reproducibility

The generating values are beta0 = 0, beta1 = 1.2, beta2 = 1, beta3 = -1,
p11 = 0.80, p10 = 0.05 and q = 0.8. Occupancy follows a Bernoulli distribution
with logistic probability based on the three predictors. Visitation follows
Bernoulli(q). At visited sites, positive records occur with probability p11
when occupied and p10 when unoccupied. Unvisited sites are coded zero.

Two visits are drawn, since visits are assumed to be independent, 
only the first visit's misclassified record is used for model fitting.  

For replicate r, the seeds are:

- Simulation: `4567 + r`.
- Random initial values: `9000 + r`.
- Stan sampling: `123 + r`.

Using the same replicate ID across scenarios gives the same data and initial
values under the same R/package/RNG configuration. Different priors still give
different posterior distributions. RNG means random number generator.
These seeds make new runs reproducible; the original unseeded simulation
cannot be reconstructed without its saved data or random-number state.

Initialisation uses independent normal draws for beta0-beta3, q from
Uniform(0.5, 0.8), theta11 from Uniform(0.15, 0.45), and theta10 from
Uniform(0.001, 0.1). Check compatibility with your
Stan model's constraints. Each fit uses four parallel chains, 2,000 warmup and
2,000 sampling iterations per chain, adapt_delta = 0.99 and max_treedepth = 15.

## Output files

All outputs are under `results_do_prior_sensitivity`. Each scenario has:

- `scenario_a/logs/`: task stdout and stderr, including the array job ID.
- `scenario_a/results/post_means_001.rds`: posterior summary for replicate 1.
- `scenario_a/settings/settings_001.rds`: intended priors, generating values,
  seeds, input paths, elapsed time and R session information.

Replace `a` with b-f and `001` with the replicate number as appropriate.
The settings directory is separate so that summary files can be combined
without accidentally reading metadata files.

Summary columns include scenario, replicate, parameter, generating true value,
posterior mean and SD, five posterior quantiles, Rhat, bulk ESS and tail ESS.
The parameters are beta0-beta3, p11, p10, theta11 and theta10. 

Re-running the same scenario/replicate overwrites its RDS files. Use a separate
study directory when changing generating settings or revising the priors.
