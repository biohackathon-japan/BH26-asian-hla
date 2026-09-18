#!/bin/bash
# SpecHLA-PG driver: SpecHLA's workflow (Wang et al. 2023) with pangenome references substituted at
# steps A (read extraction database), D (graph alignment for variant calling), E (block-linking allele
# database) and F (designation database + panel frequency prior); B (binning rule) and C (BWA/fermikit
# local assembly and realignment) are SpecHLA's own scripts, run unchanged.
#
# usage: spechla_pg.sh -n donor -1 r1.fq.gz -2 r2.fq.gz -o outdir -k fold -m MODE [-g pggb|mc] [-j threads]
#   MODE   A     : pangenome binning db; linear BWA alignment; native E db; native F (+ pangenome F for comparison)
#          AD    : + giraffe alignment to the fold's per-gene graph, surjected onto the SpecHLA reference path
#          ADEF  : + fold pangenome allele db for block linking (E) and pangenome designation (F)
# The binning step (A+B) is shared: if $outdir/../bin/COMPLETE exists its gene FASTQs are reused.
set -euo pipefail
threads=4; graph=pggb
while getopts "n:1:2:o:k:m:g:j:" opt; do case $opt in
 n) sample=$OPTARG;; 1) fq1=$OPTARG;; 2) fq2=$OPTARG;; o) outdir=$OPTARG;; k) fold=$OPTARG;; m) mode=$OPTARG;; g) graph=$OPTARG;; j) threads=$OPTARG;;
esac; done
ROOT=${SPECHLA_PG_ROOT:-/home/leechuck/hla/spechla-pg}
ENV=/home/leechuck/hla/mm/envs/asian50-spechla
export PATH=$ENV/bin:/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin:$PATH
export CONDA_PREFIX=$ENV LD_LIBRARY_PATH=$ENV/lib
NATIVE=$ENV/share/spechla/db; SCRIPT=$ENV/share/spechla/script
export SPECHLA_DB=$NATIVE SPECHLA_SCRIPT=$SCRIPT
FOLDDB=$ROOT/db/fold$fold
V2=hla_gen.format.filter.extend.DRB.no26789.v2.fasta
hlaref=$NATIVE/ref/hla.ref.extend.fa
group='@RG\tID:'$sample'\tSM:'$sample
hlas=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
case $mode in A) align=linear; edb=$NATIVE;; AD) align=graph; edb=$NATIVE;; ADEF) align=graph; edb=$FOLDDB;; *) echo bad mode; exit 1;; esac
mkdir -p $outdir $(dirname $outdir)/bin
outdir=$(readlink -f $outdir); bindir=$(readlink -f $(dirname $outdir)/bin)
fq1=$(readlink -f $fq1); fq2=$(readlink -f $fq2)
log(){ echo "[$(date +%H:%M:%S)] $*"; }
t_start=$(date +%s)

# ---------------- steps A+B: extract HLA reads against the pangenome-augmented allele db, bin with SpecHLA's rule
if [ ! -e $bindir/COMPLETE ]; then
  log "A: uniq read names"
  python3 $SCRIPT/uniq_read_name.py $fq1 $bindir/$sample.uniq.name.R1.gz
  python3 $SCRIPT/uniq_read_name.py $fq2 $bindir/$sample.uniq.name.R2.gz
  log "A: bowtie2 --very-sensitive -k 30 against fold$fold binning db"
  bowtie2 --very-sensitive -p $threads -k 30 -x $FOLDDB/ref/$V2 -1 $bindir/$sample.uniq.name.R1.gz -2 $bindir/$sample.uniq.name.R2.gz 2> $bindir/bowtie2.log |
    samtools view -bS - | samtools sort -@2 - > $bindir/$sample.map_database.bam
  samtools index $bindir/$sample.map_database.bam
  log "B: assign_reads_to_genes.py (SpecHLA, unchanged: NM<=2, no soft clip, both ends, score gap 0.1)"
  python3 $SCRIPT/assign_reads_to_genes.py -1 $bindir/$sample.uniq.name.R1.gz -2 $bindir/$sample.uniq.name.R2.gz -n $ENV/bin -o $bindir -d 0.1 -b $bindir/$sample.map_database.bam -nm 2 > $bindir/assign.log
  for hla in ${hlas[@]}; do echo -e "$hla\t$(( $(zcat $bindir/$hla.R1.fq.gz | wc -l) / 4 ))"; done > $bindir/bin_counts.tsv
  rm -f $bindir/$sample.uniq.name.R[12].gz $bindir/$sample.map_database.bam*
  echo "seconds $(( $(date +%s) - t_start ))" > $bindir/COMPLETE
fi
for hla in ${hlas[@]}; do ln -sf $bindir/$hla.R1.fq.gz $outdir/$hla.R1.fq.gz; ln -sf $bindir/$hla.R2.fq.gz $outdir/$hla.R2.fq.gz; done
t_bin=$(date +%s)

# ---------------- step C (linear) or D-alignment (graph): per-gene alignment in SpecHLA reference coordinates
{ printf '@HD\tVN:1.6\tSO:coordinate\n'; awk -v OFS='\t' '{print "@SQ","SN:"$1,"LN:"$2}' $hlaref.fai; printf "$group\n"; } > $outdir/header.sam
for hla in ${hlas[@]}; do
  if [ $align == linear ]; then
    bwa mem -t $threads -U 10000 -L 10000,10000 -R "$group" $NATIVE/HLA/HLA_$hla/HLA_$hla.fa $outdir/$hla.R1.fq.gz $outdir/$hla.R2.fq.gz 2>/dev/null |
      samtools view -bS -F 0x800 - | samtools sort - > $outdir/$hla.bam
  else
    G=$ROOT/graphs/fold$fold/HLA_$hla.$graph
    zip=""; [ -s $G.zipcodes ] && zip="-z $G.zipcodes"
    vg giraffe -Z $G.giraffe.gbz -d $G.dist -m $G.min $zip -f $outdir/$hla.R1.fq.gz -f $outdir/$hla.R2.fq.gz -t $threads \
      -o BAM --ref-name SpecHLA -N $sample -R $sample 2> $outdir/$hla.giraffe.log |
      samtools view -bS -F 0x800 - | samtools sort - > $outdir/$hla.raw.bam
    # surjected contig is named after the reference path; normalise to SpecHLA's HLA_<gene>
    samtools view -h $outdir/$hla.raw.bam | sed -E 's/SpecHLA#0#(HLA_[A-Z0-9]+)/\1/g' | samtools view -b -o $outdir/$hla.bam - && rm $outdir/$hla.raw.bam
  fi
  samtools index $outdir/$hla.bam
done
samtools merge -f -h $outdir/header.sam $outdir/$sample.merge.bam $outdir/A.bam $outdir/B.bam $outdir/C.bam $outdir/DPA1.bam $outdir/DPB1.bam $outdir/DQA1.bam $outdir/DQB1.bam $outdir/DRB1.bam
samtools index $outdir/$sample.merge.bam
samtools flagstat $outdir/$sample.merge.bam > $outdir/merge.flagstat.txt
log "C: fermikit local assembly + realignment (SpecHLA run.assembly.realign.sh, unchanged)"
bash $SCRIPT/run.assembly.realign.sh $sample $outdir/$sample.merge.bam $outdir 70 $SCRIPT/whole/select.region.txt $threads
t_aln=$(date +%s)

# ---------------- step D calling: freebayes -p 3 (SpecHLA, unchanged) + low-depth masking
log "D: freebayes"
freebayes -a -f $hlaref -p 3 $outdir/$sample.realign.sort.bam > $outdir/$sample.realign.vcf
bgzip -f $outdir/$sample.realign.vcf; tabix -f $outdir/$sample.realign.vcf.gz
zcat $outdir/$sample.realign.vcf.gz | grep "#" > $outdir/$sample.realign.filter.vcf
bcftools filter -t HLA_A:1000-4503,HLA_B:1000-5081,HLA_C:1000-5304,HLA_DPA1:1000-10775,HLA_DPB1:1000-12468,HLA_DQA1:1000-7492,HLA_DQB1:1000-8480,HLA_DRB1:1000-12229 \
  $outdir/$sample.realign.vcf.gz | grep -v "#" >> $outdir/$sample.realign.filter.vcf || true
bam=$outdir/$sample.realign.sort.bam; vcf=$outdir/$sample.realign.filter.vcf
samtools depth -aa $bam > $bam.depth
python3 $SCRIPT/mask_low_depth_region.py -c $bam.depth -o $outdir -w 20 -d 5 -f False > /dev/null
t_call=$(date +%s)

# ---------------- step E: SpecHap phasing; unlinked blocks linked through the allele database ($edb)
log "E: phase_variants.py with db=$edb"
for hla in ${hlas[@]}; do
  python3 $SCRIPT/phase_variants.py -o $outdir -b $bam -s nothing -v $vcf --fq1 $outdir/$hla.R1.fq.gz --fq2 $outdir/$hla.R2.fq.gz \
    --gene HLA_$hla --freq_bias 0.05 --snp_qual 0.01 --snp_dp 5 --ref $edb/ref/HLA_$hla.fa --tgs NA --nanopore NA --hic_fwd NA --hic_rev NA --tenx NA \
    --sa $sample --weight_imb 0 --exon 0 --thread_num $threads --use_database 1 --trio None --db $edb > $outdir/phase.HLA_$hla.log 2>&1
done
t_phase=$(date +%s)

# ---------------- step F: designation. Native (annoHLA.pl + G group) and pangenome (designate.py) both run on the same sequences.
log "F: native annoHLA.pl"
perl $SCRIPT/whole/annoHLA.pl -s $sample -i $outdir -p Unknown -d $NATIVE/HLA -r whole > $outdir/annoHLA.log 2>&1
python3 $SCRIPT/whole/g_group_annotation.py -s $sample -i $outdir -p Unknown -j $threads --db $NATIVE > $outdir/ggroup.log 2>&1
log "F: pangenome designate.py"
python3 $ROOT/scripts/designate.py --sample $sample --indir $outdir --fold $fold --root $ROOT --threads $threads > $outdir/designate.log 2>&1
t_end=$(date +%s)
printf 'stage\tseconds\nbin\t%d\nalign_assembly\t%d\ncall\t%d\nphase\t%d\ndesignate\t%d\ntotal\t%d\n' $((t_bin-t_start)) $((t_aln-t_bin)) $((t_call-t_aln)) $((t_phase-t_call)) $((t_end-t_phase)) $((t_end-t_start)) > $outdir/timing.tsv
rm -f $outdir/*.bam.depth $outdir/extract.fa $outdir/rematch.bam $outdir/assembly.fa* $outdir/$sample.realign.bam $outdir/[A-Z]*[0-9].bam $outdir/[ABC].bam $outdir/*.bam.bai.tmp
touch $outdir/COMPLETE
log "done mode=$mode align=$align graph=$graph"
