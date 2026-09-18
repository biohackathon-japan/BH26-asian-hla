#!/usr/bin/env python3
"""Build per-fold, leakage-free pangenome references for SpecHLA-PG (runs on the cluster).

For fold k every haplotype listed for fold k in excluded_paths.tsv (the 8 dev donors of the fold
plus pedigree relatives) is removed from every pangenome-derived reference:
  db/fold{k}/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta  step A binning index
        = SpecHLA lite IMGT db (6172 alleles, 3.38) + panel gene+-2kb sequences of 11 loci,
          panel records renamed '<GENE>*pg|<sample>#<hap>' so SpecHLA's binning (gene = name.split('*')[0])
          works unchanged. bowtie2 index built alongside.
  db/fold{k}/HLA/whole/HLA_<G>.fasta  step E/F allele database (block linking + designation)
        = IPD-IMGT/HLA 3.65.0 genomic alleles (<G>_gen.fasta) + panel records (PanSN names). makeblastdb.
  graphs/fold{k}/HLA_<G>.in.fa  step D graph input: SpecHLA representative reference 'SpecHLA#0#HLA_<G>'
        (gene + 1 kb flanks, the coordinate system of every downstream SpecHLA script) + panel records.
Everything else in db/fold{k} is a symlink to the native SpecHLA db (step C references, exon db, frequencies).
"""
import argparse,csv,os,subprocess,sys
from pathlib import Path

GENES=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1']
BIN_GENES=GENES+['DRB3','DRB4','DRB5']
SPECHLA_DB=Path('/home/leechuck/hla/mm/envs/asian50-spechla/share/spechla/db')
PANEL=Path('/home/leechuck/hla/hla-typer/source/genes')
V2='hla_gen.format.filter.extend.DRB.no26789.v2.fasta'

def read_fasta(path):
    name=None;seqs=[];out=[]
    for line in open(path):
        if line.startswith('>'):
            if name is not None:out.append((name,''.join(seqs)))
            name=line[1:].strip();seqs=[]
        else:seqs.append(line.strip())
    if name is not None:out.append((name,''.join(seqs)))
    return out

def write_fasta(path,records,width=80):
    with open(path,'w') as f:
        for name,seq in records:
            f.write('>'+name+'\n')
            for i in range(0,len(seq),width):f.write(seq[i:i+width]+'\n')

def excluded_haps(source,fold):
    """hap_ids excluded for a fold: excluded_paths.tsv rows of the fold, plus every haplotype whose
    validation group (relatives) contains one of those donors (belt and braces)."""
    rows=[r for r in csv.DictReader(open(source/'excluded_paths.tsv'),delimiter='\t') if int(r['fold'])==fold]
    haps={r['hap_id'] for r in rows};donors={r['donor'] for r in rows}|{r['sample'] for r in rows}
    groups=list(csv.DictReader(open(source/'validation_groups.tsv'),delimiter='\t'))
    hit_groups={g['group'] for g in groups if g['donor'] in donors or g['hap_id'].split('#')[0] in donors}
    extra={g['hap_id'] for g in groups if g['group'] in hit_groups}
    return haps|extra,sorted(extra-haps)

def panel_records(gene,excluded):
    recs=read_fasta(PANEL/f'HLA-{gene}.fa');kept=[];dropped=[]
    for name,seq in recs:
        hap='#'.join(name.split('#')[:2])
        (dropped if hap in excluded else kept).append((name.split()[0],seq.upper()))
    return kept,dropped

def link_tree(src,dst,skip):
    """Symlink every entry of src into dst except names in skip."""
    dst.mkdir(parents=True,exist_ok=True)
    for p in src.iterdir():
        if p.name in skip:continue
        t=dst/p.name
        if not t.exists() and not t.is_symlink():t.symlink_to(p)

def run(cmd):
    print('+',cmd,flush=True);subprocess.run(cmd,shell=True,check=True)

def build_fold(root,source,fold,threads,index):
    excluded,extra=excluded_haps(source,fold)
    log={'fold':fold,'excluded_haps':len(excluded),'extra_from_groups':extra}
    db=root/'db'/f'fold{fold}';gdir=root/'graphs'/f'fold{fold}';gdir.mkdir(parents=True,exist_ok=True)
    # --- step A: binning database
    link_tree(SPECHLA_DB/'ref',db/'ref',skip={p.name for p in (SPECHLA_DB/'ref').iterdir() if p.name.startswith(V2)})
    recs=[(n.split()[0],s) for n,s in read_fasta(SPECHLA_DB/'ref'/V2)]
    n_imgt=len(recs);counts={}
    for gene in BIN_GENES:
        kept,dropped=panel_records(gene,excluded)
        counts[gene]={'kept':len(kept),'dropped':len(dropped)}
        recs+=[(f'{gene}*pg|'+'#'.join(n.split('#')[:2]),s) for n,s in kept]
    write_fasta(db/'ref'/V2,recs)
    log['binning_db']={'imgt_lite':n_imgt,'panel':counts,'total':len(recs)}
    # --- step E/F database and step D graph inputs
    link_tree(SPECHLA_DB/'HLA',db/'HLA',skip={'whole'})
    link_tree(SPECHLA_DB/'HLA'/'whole',db/'HLA'/'whole',skip={p.name for p in (SPECHLA_DB/'HLA'/'whole').iterdir() if not p.name.startswith('HLA_DRB1.exon')})
    ref=dict(read_fasta(SPECHLA_DB/'ref'/'hla.ref.extend.fa'))
    for gene in GENES:
        kept,dropped=panel_records(gene,excluded)
        imgt=[(n.split()[1],s.upper()) for n,s in read_fasta(source/'imgt_gen'/f'{gene}_gen.fasta')]
        write_fasta(db/'HLA'/'whole'/f'HLA_{gene}.fasta',imgt+kept)
        write_fasta(gdir/f'HLA_{gene}.in.fa',[(f'SpecHLA#0#HLA_{gene}',ref[f'HLA_{gene}'].upper())]+kept)
        log.setdefault('ef_db',{})[gene]={'imgt_365_gen':len(imgt),'panel':len(kept)}
    if index:
        if not (db/'ref'/(V2+'.rev.2.bt2')).exists():run(f'cd {db}/ref && bowtie2-build --threads {threads} -q {V2} {V2}')
        for gene in GENES:
            run(f'cd {db}/HLA/whole && makeblastdb -in HLA_{gene}.fasta -dbtype nucl -out HLA_{gene} > /dev/null && samtools faidx HLA_{gene}.fasta')
    return log

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default='/home/leechuck/hla/spechla-pg');p.add_argument('--folds',default='0,1,2,3,4')
    p.add_argument('--threads',type=int,default=4);p.add_argument('--no-index',action='store_true');a=p.parse_args()
    root=Path(a.root);source=root/'source';import json
    logs=[build_fold(root,source,int(k),a.threads,not a.no_index) for k in a.folds.split(',')]
    (root/'db'/'build_log.json').write_text(json.dumps(logs,indent=1))
    print(json.dumps(logs,indent=1))
