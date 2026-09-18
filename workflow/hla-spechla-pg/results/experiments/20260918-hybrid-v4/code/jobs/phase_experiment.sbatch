#!/bin/bash
#SBATCH --job-name=spg_phase_v2
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/leechuck/hla/spechla-pg/logs/phase-v2-%A_%a.log
set -euo pipefail
root=/home/leechuck/hla/spechla-pg
: "${RELEASE:?Set RELEASE to an immutable experiment directory}"
read -r donor fold cram < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$root/source/donors_folds.tsv")
for arm in ${ARMS:-A-replay A-fixed AE-ipd365 AE-panel}; do
  "$root/../mm/envs/asian50-spechla/bin/python3" "$RELEASE/code/experiment.py" \
    --root "$root" --release "$RELEASE" phase --donor "$donor" --fold "$fold" --arm "$arm" \
    --threads "$SLURM_CPUS_PER_TASK"
done
