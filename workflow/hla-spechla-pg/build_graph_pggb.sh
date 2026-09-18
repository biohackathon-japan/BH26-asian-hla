#!/bin/bash
# PGGB-style graph (wfmash + seqwish + smoothxg/gfaffix/odgi normalisation)
# plus giraffe indexes for one fold and gene: build_graph_pggb.sh FOLD GENE THREADS
# Input graphs/fold$FOLD/HLA_$GENE.in.fa (SpecHLA#0#HLA_<gene> reference + fold training panel sequences).
set -euo pipefail
fold=$1; gene=$2; threads=${3:-4}
cd "${SPECHLA_PG_ROOT:-/home/leechuck/hla/spechla-pg}"
export PATH=/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin:/home/leechuck/hla/mm/envs/pggb053/bin:/home/leechuck/hla/mm/envs/asian50-spechla/bin:$PATH
g=HLA_$gene; in=graphs/fold$fold/$g.in.fa; out=graphs/fold$fold/pggb_$g; P=graphs/fold$fold/$g.pggb
if [ -e $P.giraffe.gbz ] && [ -e $P.min ] && [ "${FORCE:-0}" != 1 ]; then echo ALREADY $fold $g; exit 0; fi
n=$(grep -c '>' $in); echo "fold $fold $g haplotypes=$n $(date)"
mkdir -p $out
t0=$(date +%s)
gfa=$(ls $out/*.seqwish.gfa 2>/dev/null | head -1 || true)
if [ -z "$gfa" ]; then
  rm -rf $out; mkdir -p $out
  bgzip -c -@$threads $in > $out/$g.fa.gz; samtools faidx $out/$g.fa.gz
  wfmash -s 5000 -l 25000 -p 90 -n $((n-1)) -k 19 -H 0.001 -X -t $threads --tmp-base $out $out/$g.fa.gz --approx-map > $out/$g.mappings.paf 2> $out/wfmash.log
  wfmash -s 5000 -l 25000 -p 90 -n $((n-1)) -k 19 -H 0.001 -X -t $threads --tmp-base $out $out/$g.fa.gz -i $out/$g.mappings.paf --invert-filtering > $out/$g.alignments.paf 2>> $out/wfmash.log
  seqwish -s $out/$g.fa.gz -p $out/$g.alignments.paf -k 19 -f 0 -g $out/$g.seqwish.gfa -B 10000000 -t $threads --temp-dir $out -P > $out/seqwish.log 2>&1
  gfa=$out/$g.seqwish.gfa
fi
t1=$(date +%s); echo "wfmash+seqwish seconds=$((t1-t0)) gfa=$gfa"
# PGGB's normalisation stage. Without it the raw seqwish graph keeps every redundant parallel walk, and
# giraffe's gapless-extension step explodes on the DRB1 tangle (measured: 461 s and 650 MB for a single
# read pair). smoothxg + gfaffix + odgi make the graph linear enough for short-read mapping.
norm=$out/$g.smooth.final.gfa
if [ ! -s $norm ]; then
  smoothxg -t $threads -T $threads -g $gfa -r $n --base $out --chop-to 100 -I 0.9 -R 0 -j 0 -e 0 \
    -l 700,900,1100 -p 1,19,39,3,81,1 -O 0.001 -Y $((100*n)) -d 0 -D 0 -V -o $out/$g.smooth.gfa > $out/smoothxg.log 2>&1
  gfaffix $out/$g.smooth.gfa -o $out/$g.smooth.fix.gfa > $out/$g.affixes.tsv 2> $out/gfaffix.log
  odgi build -t $threads -g $out/$g.smooth.fix.gfa -o - -O | odgi unchop -t $threads -i - -o - |
    odgi sort -p Ygs --temp-dir $out -t $threads -i - -o $out/$g.final.og
  odgi view -i $out/$g.final.og -g > $norm
fi
gfa=$norm
t1b=$(date +%s); echo "smoothxg+gfaffix+odgi seconds=$((t1b-t1)) gfa=$gfa"
vg gbwt -G $gfa --gbz-format -g $out/$g.raw.gbz --num-threads $threads
vg gbwt -Z $out/$g.raw.gbz --set-reference SpecHLA --gbz-format -g $P.giraffe.gbz
vg index -j $P.dist -t $threads $P.giraffe.gbz
vg minimizer -t $threads -d $P.dist -z $P.zipcodes -o $P.min $P.giraffe.gbz
t2=$(date +%s); echo "index seconds=$((t2-t1))"
vg stats -lz $P.giraffe.gbz > $P.stats.txt; cat $P.stats.txt
printf 'builder\tfold\tgene\thaplotypes\tbuild_seconds\tindex_seconds\n' > $P.timing.tsv
printf 'pggb\t%s\t%s\t%s\t%s\t%s\n' $fold $g $n $((t1-t0)) $((t2-t1)) >> $P.timing.tsv
rm -f $out/$g.raw.gbz $out/*.paf $out/$g.fa.gz* $out/$g.smooth.gfa $out/$g.smooth.fix.gfa $out/$g.affixes.tsv
echo GRAPHDONE $fold $g $(date)
