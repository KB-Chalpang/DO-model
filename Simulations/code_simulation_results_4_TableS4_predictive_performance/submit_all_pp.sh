#!/bin/bash
# Submit all five model arrays for both q values, or for one selected q.
# Run on the login/submission node:
#   bash submit_all_pp.sh
#   bash submit_all_pp.sh 0.5
#   bash submit_all_pp.sh 1
# This launcher calls sbatch; it is not itself a SLURM batch job.
set -euo pipefail

if (( $# > 1 )); then
  echo 'Usage: bash submit_all_pp.sh [all|0.5|1]' >&2
  exit 1
fi
case "${1:-all}" in
  all) q_values=(0.5 1) ;;
  0.5|1) q_values=("$1") ;;
  *) echo 'Usage: bash submit_all_pp.sh [all|0.5|1]' >&2; exit 1 ;;
esac

# Submit from the scripts directory so relative R and output paths resolve there.
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd -- "$script_dir"
command -v sbatch >/dev/null || { echo 'sbatch is not available.' >&2; exit 1; }
models=(glm maxent brt rf do)

# Check every required script before submitting the first job.
for model in "${models[@]}"; do
  for file in "submit_pp_${model}.sh" "model_ajob_pp_${model}.R"; do
    [[ -r "$file" ]] || { echo "Missing or unreadable file: $file" >&2; exit 1; }
  done
done

# Each existing submission script supplies its own resources and 500-task array.
# Optional ARRAY lets every model use a smaller array for an initial check.
options=()
if [[ -n "${ARRAY:-}" ]]; then options+=("--array=$ARRAY"); fi
for q in "${q_values[@]}"; do
  for model in "${models[@]}"; do
    echo "Submitting $model with q=$q"
    if ! sbatch "${options[@]}" "submit_pp_${model}.sh" "$q"; then
      echo "Submission failed for $model, q=$q. Stopping; earlier jobs remain submitted." >&2
      echo 'Check the job IDs printed above before retrying to avoid duplicate submissions.' >&2
      exit 1
    fi
  done
done
