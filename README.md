# DO-Model
# Simulation studies and case studies

This project contains simulation studies evaluating occupancy models under
alternative observation processes, model specifications and prior assumptions.
This README outlines the workflow. Each study has its own README with the
required files, settings, submission commands and output descriptions.

## 1. Prepare the covariates

Generate and save the fixed covariates before running the simulations. Use the
covariate file specified in each study's README: the model-comparison and prior
sensitivity studies use `covariates_fixed.rds`, cross-validation uses
`covariates_fixed2.rds`, and zero-count sensitivity uses `covariates_fixed4.rds`.
Check the input paths and required R packages before submitting jobs.

## 2. Compile the Stan models

Run `compile_stan_models.R` using its submission script to prepare the True,
RL (occupancy models) and detection-only (DO-M) models. Keep the compiled executables accessible
alongside the saved model objects. Prior sensitivity has a separate compiler,
described in Step 6. Complete compilation before submitting model-fitting jobs.

## 3. Compare model specifications

Run the following model families:

- **True:** fits the simulated occupancy states directly.
- **RL:** fits detection/non-detection records, with one and two visits evaluated separately.
- **DO:** fits first-visit records across all sites, using detection-only/background data.

The misspecification scenarios omit a relevant predictor, replace a relevant
predictor with an irrelevant predictor, or add an irrelevant predictor.
The correctly specified comparison uses the actual predictors `x1`, `x2`
and `x3`. Follow the corresponding study READMEs to submit the arrays and
locate the results for each model and scenario.

## 4. Assess sensitivity to the number of zeros

Run the naive (GLM) and DO-M size of background data sensitivity scripts separately. 
Each analysis retains all positive records while varying the number of sampled background records.
The visitation probability is `q = 0.1`, providing a large pool of unvisited
sites that contribute background data; these zeros are not confirmed absences, hence, are pseudo-absences.

The baseline scripts use an occupancy intercept of zero. Follow the study
README to repeat the analysis with intercepts of 1.5 and -1.5, using separate
output directories.

## 5. Evaluate predictive performance

Run five-fold cross-validation for GLM, Maxent, boosted regression trees,
random forest and DO at `q = 0.5` and `q = 1`. Models are fitted to observed
records and evaluated against simulated occupancy using the implemented
log-score, Brier-score and AUC calculations.

See `README_CROSS_VALIDATION.md` for more details and model settings.

## 6. Assess sensitivity to prior specifications

Compare six prior scenarios for the true-positive probability `p11` and the
false-positive probability `p10`: the baseline specification and five
alternatives. The data-generating parameters remain fixed across scenarios.

First run `compile_prior_sensitivity_models.R` through its compilation
submission script. It creates `prior_sensitivity_compiled_models/`. After
successful compilation, submit `submit_do_prior_sensitivity.sh` to fit all
six scenarios. See `README_DO_PRIOR_SENSITIVITY.md` for the prior table,
compilation dependency command and output structure.

## Case studies

The project also includes case studies in the United Kingdom (UK) and
New Zealand. Their analysis code and data are available in the folders UK and New Zealand.

