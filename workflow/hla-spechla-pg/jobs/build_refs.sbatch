#!/bin/bash
#SBATCH --job-name=spg_refs
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/home/leechuck/hla/spechla-pg/logs/build_refs-%j.log
# Build the five fold-specific, leakage-free reference sets (binning db + bowtie2 index, E/F allele db + blast db, graph inputs).
set -euo pipefail
cd /home/leechuck/hla/spechla-pg
export PATH=/home/leechuck/hla/mm/envs/asian50-spechla/bin:$PATH
python3 scripts/build_refs.py --threads 4 --folds "${FOLDS:-0,1,2,3,4}"
echo REFSDONE $(date)
