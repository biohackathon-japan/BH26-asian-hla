#!/bin/bash
# Minigraph-Cactus graph + giraffe indexes for one fold and gene: build_graph_mc.sh FOLD GENE THREADS
set -euo pipefail
fold=$1; gene=$2; threads=${3:-4}
cd /home/leechuck/hla/spechla-pg
C=/home/leechuck/hla/cactus/cactus-bin-v3.3.0
. $C/venv-cactus-v3.3.0/bin/activate
export PATH=$C/bin:/home/leechuck/hla/mm/envs/asian50-spechla/bin:$PATH PYTHONPATH=$C/lib:${PYTHONPATH:-}
g=HLA_$gene; in=graphs/fold$fold/$g.in.fa; W=$(pwd)/graphs/fold$fold/mc_$g; P=graphs/fold$fold/$g.mc
if [ -e $P.giraffe.gbz ] && [ -e $P.min ]; then echo ALREADY $fold $g; exit 0; fi
rm -rf $W; mkdir -p $W/fa $W/tmp; export TMPDIR=$W/tmp
# one FASTA per haplotype named sample.hap; the SpecHLA reference is sample 'SpecHLA' (no haplotype suffix)
python3 - "$in" "$W" <<'PY'
import sys
inp,W=sys.argv[1],sys.argv[2]
seqfile=open(W+'/seqfile.txt','w');name=None;buf=[]
def flush():
    if name is None:return
    s,h,c=name.split('#')
    nm='SpecHLA' if s=='SpecHLA' else s+'.'+h
    open(f'{W}/fa/{nm}.fa','w').write('>'+c+'\n'+''.join(buf)+'\n');seqfile.write(f'{nm}\t{W}/fa/{nm}.fa\n')
for line in open(inp):
    if line.startswith('>'):flush();name=line[1:].split()[0];buf=[]
    else:buf.append(line.strip())
flush();seqfile.close()
PY
n=$(wc -l < $W/seqfile.txt); echo "fold $fold $g haplotypes=$n $(date)"
t0=$(date +%s)
cactus-pangenome $W/js $W/seqfile.txt --outDir $W/out --outName $g --reference SpecHLA \
  --batchSystem single_machine --maxCores $threads --maxMemory 40G --workDir $TMPDIR \
  --mapCores $threads --consCores $threads --indexCores $threads --gfa full --gbz full --giraffe full --logFile $W/cactus.log > $W/run.log 2>&1
t1=$(date +%s); echo "cactus seconds=$((t1-t0))"
ls $W/out
cp $W/out/$g.full.gbz $P.giraffe.gbz; cp $W/out/$g.full.dist $P.dist
cp $W/out/$g.full.shortread.withzip.min $P.min 2>/dev/null || cp $W/out/$g.full.min $P.min
cp $W/out/$g.full.shortread.zipcodes $P.zipcodes 2>/dev/null || true
vg stats -lz $P.giraffe.gbz > $P.stats.txt; cat $P.stats.txt
printf 'builder\tfold\tgene\thaplotypes\tbuild_seconds\tindex_seconds\n' > $P.timing.tsv
printf 'mc\t%s\t%s\t%s\t%s\t0\n' $fold $g $n $((t1-t0)) >> $P.timing.tsv
rm -rf $W/js $W/tmp
echo GRAPHDONE $fold $g $(date)
