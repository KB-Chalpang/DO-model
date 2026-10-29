#!/bin/bash
#SBATCH --job-name=sens_N_naive_beta0_0
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --array=1-10
#SBATCH --output=logs_sens_naive/%A_%a.out
#SBATCH --error=logs_sens_naive/%A_%a.err

# One task fits one zero-count setting; tasks 1-10 use the R script's grid.
# Four CPUs run four MCMC chains. Memory is 32 GB per CPU (128 GB per task).
# Submit from the folder containing zero_sensitivity_naive.R:
#   mkdir -p logs_sens_naive
#   sbatch submit_sensitivity_naive.sh
# SLURM opens its output files before executing this script, so creating the
# log folder inside the job would be too late.
#
# Other intercepts: set beta0 to 1.5 or -1.5 in zero_sensitivity_naive.R.
# Use descriptive output names to keep runs separate, for example:
#   results_sens_naive_beta0_1p5 and logs_sens_naive_beta0_1p5
#   results_sens_naive_beta0_m1p5 and logs_sens_naive_beta0_m1p5
# Update results_dir in R and both output/error directives above, then create
# the corresponding log folder before submission. Change job-name for clarity.
# The array remains 1-10. See README_ZERO_SENSITIVITY.md.
# Finish existing jobs before editing shared scripts. For simultaneous runs,
# make intercept-specific copies and update the Rscript filename below.

set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"
module load R/4.4.1

# The R script reads SLURM_ARRAY_TASK_ID directly and creates results_dir.
Rscript zero_sensitivity_naive.R
