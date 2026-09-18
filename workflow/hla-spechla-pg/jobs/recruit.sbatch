#!/bin/bash
#SBATCH --job-name=spg_recruit
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/home/leechuck/hla/spechla-pg/logs/recruit-%A_%a.log
# HLA read recruitment from a 1000G 30x CRAM (S3), identical logic to hla-typer/stageA/recruit.sbatch:
# chr6 MHC +/-100 kb (= graph ctg1), every chr6 alt, all 525 HLA-* contigs, complete pairs (-P), plus unmapped pairs.
# Array index i -> data line i+1 of source/donors_folds.tsv (donor, fold, cram URL). Reuses the main session's
# recruitment read-only when its COMPLETE sentinel exists.
set -euo pipefail
cd /home/leechuck/hla/spechla-pg
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
read -r donor fold cram < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" source/donors_folds.tsv)
test -n "$donor"
out=reads/$donor; mkdir -p "$out"
if [ -e "$out/COMPLETE" ]; then echo ALREADY "$donor"; exit 0; fi
other=/home/leechuck/hla/hla-typer/reads/$donor
for i in $(seq 1 60); do [ -e "$other/COMPLETE" ] && break; sleep 30; done   # shared recruitment is normally finished; wait up to 30 min
if [ -e "$other/COMPLETE" ] && [ -s "$other/r1.fq.gz" ]; then
  ln -sf "$other/r1.fq.gz" "$out/r1.fq.gz"; ln -sf "$other/r2.fq.gz" "$out/r2.fq.gz"
  cp "$other/flagstat.txt" "$out/flagstat.txt" 2>/dev/null || true; cp "$other/contig_classes.tsv" "$out/" 2>/dev/null || true
  echo "reused $other" > "$out/COMPLETE"; echo REUSED "$donor"; exit 0
fi
bed=reads/regions.bed
if [ ! -s "$bed" ]; then
 { printf 'chr6\t28410700\t33580488\n'; awk -v OFS='\t' '$1~/^chr6_.*_alt$/ || $1~/^HLA-/{print $1,0,$2}' "$ref.fai"; } > "$bed.$$.tmp"
 mv -f "$bed.$$.tmp" "$bed"
fi
t0=$(date +%s)
samtools view -@4 -T "$ref" -P -M -L "$bed" -C -o "$out/regions.cram" "$cram"
samtools view -@4 -T "$ref" -C -o "$out/unmapped.cram" "$cram" '*'
samtools cat -o "$out/recruited.cram" "$out/regions.cram" "$out/unmapped.cram"
rm "$out/regions.cram" "$out/unmapped.cram"
samtools collate -@4 -u -O --reference "$ref" "$out/recruited.cram" "$out/collate" |
 samtools fastq -@4 -F 0x900 -n -1 "$out/r1.fq.gz" -2 "$out/r2.fq.gz" -0 /dev/null -s "$out/singletons.fq.gz" -
samtools flagstat -@4 "$out/recruited.cram" > "$out/flagstat.txt"
samtools view -@4 -F 0x900 "$out/recruited.cram" | awk -F'\t' '{c=($3=="*")?"unmapped":($3~/^HLA-/)?"HLA_contig":($3~/_alt$/)?"chr6_alt":($3=="chr6")?"chr6":"other_mate";n[c]++} END{for(k in n)print k"\t"n[k]}' | sort > "$out/contig_classes.tsv"
rm -f "$out/recruited.cram"
echo "seconds $(( $(date +%s) - t0 ))" > "$out/COMPLETE"
echo COMPLETE "$donor"
