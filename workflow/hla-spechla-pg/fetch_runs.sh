#!/bin/bash
# Fetch the small per-run result files from the cluster into results/runs/ (tar streamed over the two-hop ssh).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/runs
ssh -o BatchMode=yes ddbj "ssh -o BatchMode=yes -o IdentitiesOnly=yes -i ~/.ssh/id_ed25519 a001 'cd /home/leechuck/hla/spechla-pg/runs && find . -maxdepth 3 \( -name \"hla.result*.txt\" -o -name \"hla.result.pg.details.tsv\" -o -name \"hla.allele.*.fasta\" -o -name \"timing.tsv\" -o -name \"time.log\" -o -name \"bin_counts.tsv\" -o -name \"merge.flagstat.txt\" -o -name \"*_freq.txt\" -o -name \"COMPLETE\" -o -name \"driver.log\" -o -name \"*.giraffe.log\" -o -name \"bowtie2.log\" -o -name \"assign.log\" \) -print0 | tar czf - --null -T -'" | tar xzf - -C results/runs
find results/runs -name COMPLETE | sort | sed 's#results/runs/##' | tr '\n' ' '; echo
