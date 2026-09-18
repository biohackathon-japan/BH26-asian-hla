#!/bin/bash
#SBATCH --job-name=spg_desig
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=/home/leechuck/hla/spechla-pg/logs/redesignate-%j.log
# Re-run step F (scripts/designate.py) on every existing run directory, without re-running A-E.
# Used after a change to the designation logic or database. MODES / DONORS env vars restrict the sweep.
set -uo pipefail
shopt -s extglob
cd /home/leechuck/hla/spechla-pg
export PATH=/home/leechuck/hla/mm/envs/hlatyper/bin:$PATH
while read -r donor fold cram; do
  [ "$donor" == donor ] && continue
  if [ -n "${DONORS:-}" ] && [[ " $DONORS " != *" $donor "* ]]; then continue; fi
  for mode in ${MODES:-native A AD ADEF}; do
    out=runs/$donor/$mode
    [ -e "$out/hla.allele.1.HLA_A.fasta" ] || continue
    python3 scripts/designate.py --sample "$donor" --indir "$out" --fold "$fold" --root "$(pwd)" \
      > "$out/designate.log" 2>&1 && echo "DESIGNATED $donor $mode" || echo "FAILED $donor $mode"
  done
done < source/donors_folds.tsv
echo ALLDESIGNATED $(date)
