#!/bin/bash
#SBATCH --job-name=correct_models
#SBATCH --array=1-800
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/dev/null
#SBATCH --error=slurm_correct_setup_%A_%a.err

# Run all four configurations with one submission:
#   sbatch submit_correct_models.sh
# Submit from the directory containing simulate_correct_models.R and compiled_models/.
# Each configuration has logs/ and results/ subfolders under results_correct/.
# Folders are created by the job before redirecting its output. SLURM setup
# errors remain in slurm_correct_setup_<job>_<task>.err in the submission directory.
#
# Tasks 1-200: true; 201-400: rl_T1; 401-600: rl_T2; 601-800: do.
# Each configuration uses replicate IDs 1-200, giving matching simulation seeds
# across configurations. To change the replicate count, update both
# replicates_per_configuration and the array upper bound (4 times that count).

set -euo pipefail

cd -- "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"

replicates_per_configuration=200
configurations=(true rl_T1 rl_T2 do)
array_task_id="${SLURM_ARRAY_TASK_ID:?This script requires a SLURM array job}"

if (( array_task_id < 1 || array_task_id > 4 * replicates_per_configuration )); then
  echo "Error: array task ID is outside the configured configuration range." >&2
  exit 1
fi

# Convert the global array ID into a configuration and a within-configuration replicate.
configuration_index=$(( (array_task_id - 1) / replicates_per_configuration ))
replicate_id=$(( (array_task_id - 1) % replicates_per_configuration + 1 ))
configuration="${configurations[configuration_index]}"

# Scenario folders are created at runtime, so submission needs no setup commands.
configuration_dir="${SLURM_SUBMIT_DIR}/results_correct/${configuration}"
mkdir -p -- "${configuration_dir}/logs" "${configuration_dir}/results"
log_stem="${configuration_dir}/logs/replicate_${replicate_id}_${SLURM_ARRAY_JOB_ID}_${array_task_id}"
exec >"${log_stem}.out" 2>"${log_stem}.err"

module load R/4.4.1
export R_LIBS_USER="${HOME}/.local/R_libs"

# Four allocated CPUs allow the four MCMC chains to run concurrently.
# Results are saved as results_correct/<configuration>/results/results_<replicate>.rds.
echo "Running ${configuration}, replicate ${replicate_id} (array task ${array_task_id})"
Rscript simulate_correct_models.R \
  "${replicate_id}" \
  "${configuration}" \
  "${HOME}/covariates_fixed.rds" \
  "${SLURM_SUBMIT_DIR}/compiled_models" \
  "${SLURM_SUBMIT_DIR}/results_correct"
