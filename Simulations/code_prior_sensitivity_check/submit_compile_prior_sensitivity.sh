#!/bin/bash
#SBATCH --job-name=compile_do_priors
#SBATCH --output=compile_do_priors_%j.out
#SBATCH --error=compile_do_priors_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Submit once from the directory containing the R scripts:
#   sbatch submit_compile_prior_sensitivity.sh
# This job compiles six models sequentially; it does not run simulations.
# Adjust time and memory if compilation requires more resources on your cluster.
set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:?Submit this script using sbatch}"
module load R/4.4.1
# Add R_LIBS_USER here if your packages use a custom library directory.
Rscript compile_prior_sensitivity_models.R
