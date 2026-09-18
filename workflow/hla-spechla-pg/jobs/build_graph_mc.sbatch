#!/bin/bash
#SBATCH --job-name=spg_mc
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/leechuck/hla/spechla-pg/logs/mc-%A_%a.log
# Per-fold, per-gene graph via scripts/build_graph_mc.sh. Array index = fold*8 + gene index (A B C DPA1 DPB1 DQA1 DQB1 DRB1).
set -euo pipefail
genes=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
i=${SLURM_ARRAY_TASK_ID}; fold=$((i/8)); gene=${genes[$((i%8))]}
bash /home/leechuck/hla/spechla-pg/scripts/build_graph_mc.sh $fold $gene 4
