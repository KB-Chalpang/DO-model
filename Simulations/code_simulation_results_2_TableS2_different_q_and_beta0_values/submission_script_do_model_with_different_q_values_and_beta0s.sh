#!/bin/bash
# ------------------------------------------------------------------------------
# One script for the whole simulation study. A scenario is one value of beta0 and
# one value of q. Each scenario is run as its own job array of 500 replicates, with
# its own results folder and its own log folder:
#
#   results_do<q>_beta0_<beta0>/   e.g. results_do07_beta0_1.5/   (post_means_001.rds, ...)
#   logs/do<q>_beta0_<beta0>/      e.g. logs/do07_beta0_1.5/       (<jobid>_<task>.out / .err)
#
# Usage, from the folder containing submit_ajob_do.R (either command submits ALL scenarios):
#   sbatch submit_ajob_do.sh
#   bash submit_ajob_do.sh
# To submit only some scenarios, give the values to use (space-separated):
#   BETA0_VALUES="1.5" Q_VALUES="0.7 0.8" bash submit_ajob_do.sh
# ("m" in a folder name stands for minus, e.g. beta0 = -1.5 gives beta0_m1.5)
# ------------------------------------------------------------------------------

# The lines below are for the launcher only (the part that submits the arrays); the
# resources of the replicate jobs are set in the sbatch command further down.
#SBATCH --job-name=do-launcher
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

# Scenarios and number of replicates
BETA0_VALUES=${BETA0_VALUES:-"-1.5 0 1.5"}
Q_VALUES=${Q_VALUES:-"0.5 0.7 0.9 1"}
N_REPLICATES=10

# Scenario name from q and beta0, e.g. q = 0.7, beta0 = -1.5 -> do07_beta0_m1.5
scenario_name() {
  local q_tag
  q_tag=$(printf "%02d" "$(awk -v q="$1" 'BEGIN{printf "%d", q * 10 + 0.5}')")
  echo "do${q_tag}_beta0_${2//-/m}"
}

# ---- Submit mode: not inside an array task -> submit one job array per scenario ----
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  for BETA0 in $BETA0_VALUES; do
    for Q in $Q_VALUES; do
      SCENARIO=$(scenario_name "$Q" "$BETA0")
      mkdir -p "results_${SCENARIO}" "logs/${SCENARIO}"   # SLURM needs the log folder to exist
      sbatch --array=1-$N_REPLICATES \
             --time=01:00:00 --cpus-per-task=4 --mem-per-cpu=8G \
             --job-name="${SCENARIO}" \
             --output="logs/${SCENARIO}/%A_%a.out" \
             --error="logs/${SCENARIO}/%A_%a.err" \
             --export=ALL,Q="$Q",BETA0="$BETA0" \
             "$0"
    done
  done
  exit 0
fi

# ---- Job mode: one array task = one replicate of one scenario ------------------------
Q=${Q:-0.5}
BETA0=${BETA0:-0}
SCENARIO=$(scenario_name "$Q" "$BETA0")

module load R/4.4.1
export R_LIBS_USER=$HOME/.local/R_libs         # adjust if your R packages are installed elsewhere

Rscript submit_ajob_dom_q_beta_testing.R $SLURM_ARRAY_TASK_ID $Q $BETA0 "results_${SCENARIO}"

