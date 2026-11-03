#!/bin/bash
#SBATCH --job-name=do_prior_sensitivity
#SBATCH --array=1-3000
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=4
#SBATCH --output=do_prior_setup_%A_%a.out
#SBATCH --error=do_prior_setup_%A_%a.err

# Run from the directory containing run_do_prior_sensitivity.R:
#   sbatch submit_do_prior_sensitivity.sh
# Six scenarios each have 500 replicates. The array task selects both the
# scenario and the replicate; the R script uses the replicate for its seeds.
# Startup logs use the submission directory, which exists before the job starts.
set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:?Submit this script using sbatch}"

replicates=500
scenarios=(a b c d e f)
array_id="${SLURM_ARRAY_TASK_ID:?Missing array task ID}"
if (( array_id < 1 || array_id > replicates * ${#scenarios[@]} )); then
  echo "Array task ID is outside the configured scenario/replicate range." >&2
  exit 1
fi
scenario_index=$(( (array_id - 1) / replicates ))
replicate=$(( (array_id - 1) % replicates + 1 ))
scenario="${scenarios[$scenario_index]}"

# Redirect each task's logs after creating the matching scenario directory.
log_dir="results_do_prior_sensitivity/scenario_${scenario}/logs"
mkdir -p -- "$log_dir"
log_stem="${log_dir}/${SLURM_ARRAY_JOB_ID}_${replicate}"
exec >"${log_stem}.out" 2>"${log_stem}.err"

module load R/4.4.1
# Add R_LIBS_USER here if your cluster requires a custom R package directory.
Rscript run_do_prior_sensitivity.R "$scenario" "$replicate"
