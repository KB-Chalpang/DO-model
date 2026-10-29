#!/bin/bash
#SBATCH --job-name=model_pp_rf
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=1-500
#SBATCH --output=pp_rf_setup_%A_%a.out
#SBATCH --error=pp_rf_setup_%A_%a.err

# Submit from the directory containing model_ajob_pp_rf.R:
#   sbatch submit_pp_rf.sh 0.5
#   sbatch submit_pp_rf.sh 1
# Default q is 1. Each array task is one replicate with five CV folds.
# Logs and summaries use separate q05/q1 folders. Startup failures before
# redirection are recorded in pp_rf_setup_<job>_<task>.err.
# --mem is total memory per task; resource limits may need cluster-specific tuning.

set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"
q="${1:-1}"
case "${q}" in
  0.5) q_tag=05 ;;
  1) q_tag=1 ;;
  *) echo "q must be 0.5 or 1" >&2; exit 1 ;;
esac
if (( $# > 1 )); then echo "Usage: sbatch submit_pp_rf.sh [0.5|1]" >&2; exit 1; fi

# Create directories before redirecting output; SLURM startup logs use the
# submission directory, so no pre-submission mkdir command is required.
log_dir="logs_pp_rf_q${q_tag}"
mkdir -p -- "${log_dir}"
log_stem="${log_dir}/model_pp_rf_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
exec >"${log_stem}.log" 2>"${log_stem}.err"
module load R/4.4.1
# If needed for your installation, set R_LIBS_USER after loading the module.
Rscript model_ajob_pp_rf.R "${q}"
