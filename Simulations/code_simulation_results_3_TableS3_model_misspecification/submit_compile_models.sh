#!/bin/bash
#SBATCH --job-name=compile_stan_models
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Compile and save all nine Stan models in a single job.
# Submit from the directory containing compile_stan_models.R:
#   mkdir -p logs
#   sbatch submit_compile_models.sh
# The logs directory must exist before submission because SLURM opens the
# output files before this script starts.

set -euo pipefail

cd -- "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"

# Load the R environment used by the simulation jobs. CmdStan and a compatible
# C++ toolchain must already be installed and available in this environment.
module load R/4.4.1
export R_LIBS_USER="${HOME}/.local/R_libs"

# Compile sequentially within the single allocated CPU.
export MAKEFLAGS="-j1"

if [[ ! -f compile_stan_models.R ]]; then
  echo "Error: compile_stan_models.R is missing from ${PWD}." >&2
  exit 1
fi

# Keep the RDS objects, Stan sources and executables together. Saved model
# objects refer to these executables, so retain this folder after compilation.
output_dir="${SLURM_SUBMIT_DIR}/compiled_models"
mkdir -p -- "${output_dir}"

echo "Compiling Stan models into ${output_dir}"
Rscript compile_stan_models.R "${output_dir}"
echo "Compilation completed successfully."
