#!/usr/bin/env python3
"""Exon (CDS) coordinates for every designation-database sequence (runs locally; output uploaded to source/exons.tsv).

Panel records (gene +-2 kb): CDS features of the haplotype's Immuannot GTF mapped into record coordinates;
validated against the catalogue cds_sha256. IPD-IMGT/HLA 3.65.0 genomic alleles: the allele's CDS (<gene>_nuc.fasta)
is split by one of the exon-length tuples observed in the panel for that gene and each exon is located in the
genomic sequence (exact search, in order); alleles whose CDS matches no observed tuple get no exon rows and are
scored on the whole gene only. Output columns: name, gene, exon, start0, end0 (record coordinates, 0-based half-open).
"""
import collections,csv,gzip,hashlib,re,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parent
GENES=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1']
SRC=REPO/'hla-analysis/source';IMGT=SRC/'imgt/fasta'

def read_fasta(path,key=0):
    out=[];name=None;buf=[]
    op=gzip.open if str(path).endswith('.gz') else open
    for line in op(path,'rt'):
        if line.startswith('>'):
            if name is not None:out.append((name,''.join(buf)))
            name=line[1:].split()[key];buf=[]
        else:buf.append(line.strip())
    if name is not None:out.append((name,''.join(buf)))
    return out

def panel_exons(gene,cat):
    """Yield (name, [(start0,end0),...]) for panel records; count validation failures."""
    rows=list(csv.DictReader(open(SRC/f'HLA-{gene}.regions.tsv'),delimiter='\t'))
    seqs=dict(read_fasta(SRC/f'HLA-{gene}.fa'))
    ok=bad=0;out=[]
    for r in rows:
        name=r['name'];sample,hap=r['sample'],r['haplotype']
        m=re.match(r'^(.+):(\d+)-(\d+)\((.)\)$',r['region']) or re.match(r'^(.+):(\d+)-(\d+)$',r['region'])
        contig,rs,re_=m[1],int(m[2]),int(m[3]);strand=r['strand']
        gtf=SRC/'gtf_all'/f'{sample}_{hap}.gtf.gz'
        if not gtf.exists():continue
        cds=[]
        for line in gzip.open(gtf,'rt'):
            if line.startswith('#'):continue
            f=line.rstrip('\n').split('\t')
            if len(f)<9 or f[2] not in ('CDS','stop_codon') or f[0]!=contig or f'gene_name "HLA-{gene}"' not in f[8]:continue
            s,e=int(f[3]),int(f[4])
            if s<rs or e>re_:continue
            cds.append((s,e))
        if not cds:continue
        L=re_-rs+1
        if strand=='+':iv=sorted((s-rs,e-rs+1) for s,e in cds)
        else:iv=sorted((re_-e,re_-s+1) for s,e in cds)
        merged=[]  # the stop codon is a separate GTF feature adjacent to (or inside) the last CDS
        for a,b in iv:
            if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
            else:merged.append((a,b))
        iv=merged
        seq=seqs[name].upper();concat=''.join(seq[a:b] for a,b in iv)
        if hashlib.sha256(concat.encode()).hexdigest()==cat[name]['cds_sha256']:ok+=1;out.append((name,iv))
        else:bad+=1
    return out,ok,bad

def imgt_exons(gene,tuples):
    gen=read_fasta(IMGT/f'{gene}_gen.fasta',key=1);nuc=dict(read_fasta(IMGT/f'{gene}_nuc.fasta',key=1))
    out=[];n_ok=n_nocds=n_nolen=n_nofind=0
    for name,g in gen:
        g=g.upper();cds=nuc.get(name,'').upper()
        if not cds:n_nocds+=1;continue
        found=None
        for t in tuples:
            if sum(t)!=len(cds):continue
            iv=[];pos=0;p=0;good=True
            for L in t:
                ex=cds[p:p+L];p+=L;i=g.find(ex,pos)
                if i<0:good=False;break
                iv.append((i,i+L));pos=i+L
            if good:found=iv;break
        if found is None:
            if all(sum(t)!=len(cds) for t in tuples):n_nolen+=1
            else:n_nofind+=1
            continue
        n_ok+=1;out.append((name,found))
    return out,n_ok,n_nocds,n_nolen,n_nofind

def main():
    cat={r['name']:r for r in csv.DictReader(open(REPO/'hla-analysis/results/sequence_catalogue.tsv'),delimiter='\t')}
    rows=[];report=[]
    for gene in GENES:
        pan,ok,bad=panel_exons(gene,cat)
        tuples=collections.Counter(tuple(b-a for a,b in iv) for _,iv in pan)
        order=[t for t,_ in tuples.most_common()]
        im,n_ok,n_nocds,n_nolen,n_nofind=imgt_exons(gene,order)
        report.append(f'{gene}\tpanel_ok={ok}\tpanel_bad={bad}\texon_tuples={len(order)}\ttop={order[0] if order else None}\timgt_ok={n_ok}\timgt_no_cds={n_nocds}\timgt_len_unmatched={n_nolen}\timgt_not_found={n_nofind}')
        for name,iv in pan+im:
            for k,(a,b) in enumerate(iv,1):rows.append((name,gene,k,a,b))
    with open(HERE/'source/exons.tsv','w') as f:
        f.write('name\tgene\texon\tstart0\tend0\n')
        for r in rows:f.write('\t'.join(map(str,r))+'\n')
    print('\n'.join(report));print('rows',len(rows))

if __name__=='__main__':main()
