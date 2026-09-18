#!/bin/bash
#SBATCH --job-name=spg_type
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=/home/leechuck/hla/spechla-pg/logs/type-%A_%a.log
# Type one dev donor with every configuration: native SpecHLA (baseline) and SpecHLA-PG modes A, AD, ADEF.
# Array index i -> data line i+1 of source/donors_folds.tsv. MODES / GRAPH env vars override the defaults.
set -euo pipefail
cd /home/leechuck/hla/spechla-pg
read -r donor fold cram < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" source/donors_folds.tsv)
test -n "$donor"; test -e reads/$donor/COMPLETE
r1=$(readlink -f reads/$donor/r1.fq.gz); r2=$(readlink -f reads/$donor/r2.fq.gz)
graph=${GRAPH:-pggb}
for mode in ${MODES:-native A AD ADEF}; do
  out=runs/$donor/$mode; [ $mode != native ] && [ $graph != pggb ] && out=runs/$donor/$mode.$graph
  if [ -e $out/COMPLETE ]; then echo "ALREADY $donor $mode"; continue; fi
  rm -rf $out; mkdir -p $out; t0=$(date +%s)
  if [ $mode == native ]; then
    export PATH=/home/leechuck/hla/mm/envs/asian50-spechla/bin:$PATH CONDA_PREFIX=/home/leechuck/hla/mm/envs/asian50-spechla
    /usr/bin/time -v spechla -n $donor -1 $r1 -2 $r2 -o runs/$donor/native_tmp -j 4 -u 0 -p Unknown > $out/spechla.log 2> $out/time.log || echo "native FAILED" >> $out/spechla.log
    mv runs/$donor/native_tmp/$donor/* $out/ 2>/dev/null || true; rm -rf runs/$donor/native_tmp
    printf 'stage\tseconds\ntotal\t%d\n' $(( $(date +%s) - t0 )) > $out/timing.tsv
    grep -E "Maximum resident|Elapsed" $out/time.log || true
    python3 scripts/designate.py --sample $donor --indir $out --fold $fold --root $(pwd) > $out/designate.log 2>&1 || echo "designate FAILED" >> $out/spechla.log
    touch $out/COMPLETE
  else
    /usr/bin/time -v bash scripts/spechla_pg.sh -n $donor -1 $r1 -2 $r2 -o $out -k $fold -m $mode -g $graph -j 4 > $out/driver.log 2> $out/time.log || echo "$mode FAILED" >> $out/driver.log
    grep -E "Maximum resident|Elapsed" $out/time.log || true
  fi
  echo "DONE $donor $mode $(( $(date +%s) - t0 ))s"
done
echo ALLDONE $donor
