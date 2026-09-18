#!/bin/bash
#SBATCH --job-name=spg_hybrid
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/spechla-pg/logs/hybrid-v2-%A_%a.log
set -euo pipefail
root=/home/leechuck/hla/spechla-pg
: "${RELEASE:?Set RELEASE to the new hybrid code/output directory}"
: "${PHASE_RELEASE:?Set PHASE_RELEASE to completed paired phase experiments}"
read -r donor fold cram < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$root/source/donors_folds.tsv")
python="$root/../mm/envs/asian50-spechla/bin/python3"
"$python" "$RELEASE/code/graph_diagnostic.py" --root "$root" --donor "$donor" --fold "$fold" \
  --gene DRB1 --pairs 1000000 --threads "$SLURM_CPUS_PER_TASK" --timeout 300 \
  --out "$RELEASE/diagnostics/$donor/DRB1"
for arm in A-pairing AE-ipd365-pairing; do
  case "$arm" in A-pairing) output=AG-homozygous; control=AF-full-DRB1;; AE-ipd365-pairing) output=AEG-homozygous; control=AEF-full-DRB1;; esac
  "$python" "$RELEASE/code/structural_overlay.py" --donor "$donor" --naming-only \
    --source "$PHASE_RELEASE/runs/$donor/$arm" --graph "$RELEASE/diagnostics/$donor" \
    --out "$RELEASE/runs/$donor/$control" --threads "$SLURM_CPUS_PER_TASK"
  "$python" "$RELEASE/code/structural_overlay.py" --donor "$donor" \
    --source "$PHASE_RELEASE/runs/$donor/$arm" --graph "$RELEASE/diagnostics/$donor" \
    --out "$RELEASE/runs/$donor/$output" --threads "$SLURM_CPUS_PER_TASK"
done
