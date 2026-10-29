#!/bin/bash
#SBATCH --job-name=rl_models
#SBATCH --array=1-600
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/dev/null
#SBATCH --error=slurm_rl_setup_%A_%a.err

# Run all three scenarios with one submission:
#   sbatch submit_rl_models.sh
# Submit from the directory containing simulate_rl_models.R and compiled_models/.
# Each scenario has logs/ and results/ subfolders under results_rl/.
# Folders are created by the job before redirecting its output. SLURM setup
# errors remain in slurm_rl_setup_<job>_<task>.err in the submission directory.
#
# Array tasks 1-200 run misspec1, 201-400 run misspec2, and 401-600 run misspec3.
# Each scenario uses replicate IDs 1-200, giving matching simulation seeds
# across scenarios. To change the replicate count, update both
# replicates_per_scenario and the array upper bound (3 times that count).

# To run the two-visit scenarios:
# 1. Set T_visits <- 2L in simulate_rl_models.R.
# 2. Replace results_rl with results_rl_T2 throughout this submission script.
# 3. Submit with sbatch submit_rl_models.sh; all three scenarios will use two visits.
# Finish existing jobs before editing their shared R script, or use separate
# script copies for simultaneous one-visit and two-visit runs.

set -euo pipefail

cd -- "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"

replicates_per_scenario=200
scenarios=(misspec1 misspec2 misspec3)
array_task_id="${SLURM_ARRAY_TASK_ID:?This script requires a SLURM array job}"

if (( array_task_id < 1 || array_task_id > 3 * replicates_per_scenario )); then
  echo "Error: array task ID is outside the configured scenario range." >&2
  exit 1
fi

# Convert the global array ID into a scenario and a within-scenario replicate.
scenario_index=$(( (array_task_id - 1) / replicates_per_scenario ))
replicate_id=$(( (array_task_id - 1) % replicates_per_scenario + 1 ))
scenario="${scenarios[scenario_index]}"

# Scenario folders are created at runtime, so submission needs no setup commands.
scenario_dir="${SLURM_SUBMIT_DIR}/results_rl/${scenario}"
mkdir -p -- "${scenario_dir}/logs" "${scenario_dir}/results"
log_stem="${scenario_dir}/logs/replicate_${replicate_id}_${SLURM_ARRAY_JOB_ID}_${array_task_id}"
exec >"${log_stem}.out" 2>"${log_stem}.err"

module load R/4.4.1
export R_LIBS_USER="${HOME}/.local/R_libs"

# Four allocated CPUs allow the four MCMC chains to run concurrently.
# Results are saved as results_rl/<scenario>/results/results_<replicate>.rds.
echo "Running ${scenario}, replicate ${replicate_id} (array task ${array_task_id})"
Rscript simulate_rl_models.R \
  "${replicate_id}" \
  "${scenario}" \
  "${HOME}/covariates_fixed.rds" \
  "${SLURM_SUBMIT_DIR}/compiled_models" \
  "${SLURM_SUBMIT_DIR}/results_rl"
