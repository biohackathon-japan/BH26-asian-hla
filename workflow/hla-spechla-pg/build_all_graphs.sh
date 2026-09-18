#!/bin/bash
# Build every per-fold, per-gene graph sequentially with one builder: build_all_graphs.sh pggb|mc [THREADS] [FOLDS]
# (used on the login node in the background to keep the shared compute node free for typing jobs; skips finished graphs)
set -uo pipefail
builder=$1; threads=${2:-4}; folds=${3:-"0 1 2 3 4"}
cd /home/leechuck/hla/spechla-pg
for fold in $folds; do for gene in ${GENES:-A B C DPA1 DPB1 DQA1 DQB1 DRB1}; do
  bash scripts/build_graph_$builder.sh $fold $gene $threads > logs/$builder-$fold-$gene.log 2>&1 || echo "FAILED $builder $fold $gene"
  tail -n 1 logs/$builder-$fold-$gene.log
done; done
echo ALLGRAPHS $builder $(date)
