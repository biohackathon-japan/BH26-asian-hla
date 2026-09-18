#!/usr/bin/env python3
"""Step F with the pangenome: designate SpecHLA's reconstructed haplotype sequences.

SpecHLA (annoHLA.pl) blasts each reconstructed haplotype against IPD-IMGT/HLA genomic alleles and
picks the highest identity, breaking ties with a population allele-frequency table. Here the
database is the fold's leakage-free set: IPD-IMGT/HLA 3.65.0 full-length genomic alleles plus the
training panel gene sequences (gene +-2 kb, labelled through the sequence catalogue), and the
score is exact edit distance (edlib, infix mode): first the sum over the database sequence's exons
(CDS segments from source/exons.tsv, which SpecHLA also privileges: its class II designation blasts
exon windows only), then, among exon ties, the whole gene body inside the reconstruction (the
pangenome's full-length haplotypes make this second stage meaningful). Remaining ties are broken by
(1) having a nomenclature label, (2) the number of training-panel haplotypes carrying the same
two-field label (the pangenome's own frequency prior), (3) name. Masked bases (N) in the
reconstruction count as mismatches, as in SpecHLA. Database sequences without exon annotation
(IPD alleles with an unusual CDS length, 1-3% per gene) rank after all annotated ones.

A panel haplotype with no exact IPD match at any level (31.6% of panel gene sequences; catalogue
assessment `unresolved_no_exact_match`) can win the sequence-level comparison but carries no
nomenclature label. Reporting it as `novel:<record>` would be a no-call, so the nomenclature call
falls back to the closest *labelled* record under the same ranking (`label_source=nearest:<record>`)
while the sequence-level winner stays in the `winner` column; the extra distance of that fallback is
reported as `nearest_exon_edits` / `nearest_gene_edits`.

Outputs in the run directory: hla.result.pg.txt (SpecHLA hla.result.txt layout) and
hla.result.pg.details.tsv (per haplotype: winner, label, edit distance, ties, N count).
"""
import argparse,collections,csv,os,re,sys
from pathlib import Path
GENES=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1']
FLANK=2000  # panel records are gene +- 2000 bp (verified against catalogue gene_sha256)

def read_fasta(path):
    out=[];name=None;buf=[]
    for line in open(path):
        if line.startswith('>'):
            if name is not None:out.append((name,''.join(buf)))
            name=line[1:].split()[0];buf=[]
        else:buf.append(line.strip())
    if name is not None:out.append((name,''.join(buf)))
    return out

def is_panel(name):return '#' in name

def gene_body(name,seq):
    """IMGT genomic alleles are already the gene body; panel records carry 2 kb flanks."""
    return seq[FLANK:-FLANK] if is_panel(name) else seq

def load_catalogue(path):
    cat={}
    for r in csv.DictReader(open(path),delimiter='\t'):cat[r['name']]=r
    return cat

def label_for(name,cat):
    """Nomenclature label of a database sequence: IMGT allele name, or the catalogue's exact genomic
    (else exact CDS) allele of a panel record; None for a panel record without any exact match."""
    if not is_panel(name):return name
    r=cat.get(name)
    if r is None:return None
    for col in ('exact_genomic_alleles','exact_cds_alleles'):
        v=r.get(col,'')
        if v:return v.split(';')[0]
    return None

def two_field(label):
    m=re.match(r'^[A-Z0-9]+\*(\d+):(\d+)',label or '')
    return m[0] if m else None

def filter_full_length(records):
    """Drop partial IMGT genomic entries (< 90% of the median full-length gene body)."""
    lens=sorted(len(gene_body(n,s)) for n,s in records)
    med=lens[len(lens)//2]
    return [(n,s) for n,s in records if len(gene_body(n,s))>=0.9*med]

def frequency_prior(records,cat):
    """Two-field label -> number of training-panel haplotypes carrying it."""
    c=collections.Counter()
    for n,_ in records:
        if is_panel(n):
            tf=two_field(label_for(n,cat))
            if tf:c[tf]+=1
    return c

def edit_distance(query,target):
    import edlib
    return edlib.align(query,target,mode='HW',task='distance')['editDistance']

def load_exons(path):
    """name -> [(start0,end0),...] in the record's own coordinates (panel: gene +-2 kb record; IPD: genomic allele)."""
    ex=collections.defaultdict(list)
    for r in csv.DictReader(open(path),delimiter='\t'):ex[r['name']].append((int(r['start0']),int(r['end0'])))
    return ex

INF=10**9
def choose(recon,records,cat,prior,exons=None,dist=None):
    """Return (winner_name, label, (exon_edits, gene_edits), ties, nearest) for one reconstructed haplotype.

    Rank by exon edit distance, then gene-body edit distance (ties = records equal on both). If the overall
    winner is an unlabelled (novel) panel haplotype, the nomenclature call falls back to the closest
    *labelled* record, ranked the same way, and `nearest` = (name, label, exon_edits, gene_edits); the
    returned label is then that record's. `nearest` is None when the winner itself carries a label."""
    dist=dist or edit_distance;exons=exons or {};recon=recon.upper();cache={}
    def d(q):
        if q not in cache:cache[q]=dist(q,recon)
        return cache[q]
    scored=[]
    for n,s in records:
        s=s.upper();iv=exons.get(n)
        e=sum(d(s[a:b]) for a,b in iv) if iv else INF
        scored.append((e,n,s))
    def key(n):
        lab=label_for(n,cat);tf=two_field(lab)
        return (0 if lab else 1,-prior.get(tf,0) if tf else 0,n)
    def rank(cands):
        best_e=min(e for e,_,_ in cands)
        stage=[(d(gene_body(n,s)),n) for e,n,s in cands if e==best_e]
        best_g=min(g for g,_ in stage)
        ties=sorted((n for g,n in stage if g==best_g),key=key)
        return ties[0],(best_e if best_e<INF else None,best_g),ties
    win,edits,ties=rank(scored)
    label=label_for(win,cat);nearest=None
    if label is None:
        labelled=[(e,n,s) for e,n,s in scored if label_for(n,cat)]
        if labelled:
            wn,(ne,ng),_=rank(labelled)
            label=label_for(wn,cat);nearest=(wn,label,ne,ng)
    return win,label,edits,ties,nearest

FIELDS=['gene','hap','winner','label','label_source','exon_edits','gene_edits','nearest_exon_edits','nearest_gene_edits','n_masked','recon_len','ties','tie_labels']
def designate_gene(indir,gene,records,cat,prior,exons=None):
    rows=[]
    for hap in (1,2):
        fa=Path(indir)/f'hla.allele.{hap}.HLA_{gene}.fasta'
        if not fa.exists() or os.path.getsize(fa)==0:
            rows.append({f:'' for f in FIELDS}|{'gene':gene,'hap':hap,'recon_len':0});continue
        recon=''.join(s for _,s in read_fasta(fa))
        win,label,(ee,ge),ties,near=choose(recon,records,cat,prior,exons)
        tie_labels=sorted({label_for(t,cat) or 'novel:'+t for t in ties})
        rows.append({'gene':gene,'hap':hap,'winner':win,'label':label or 'novel:'+win,
                     'label_source':'nearest:'+near[0] if near else ('exact' if label else 'none'),
                     'exon_edits':ee,'gene_edits':ge,
                     'nearest_exon_edits':near[2] if near else '','nearest_gene_edits':near[3] if near else '',
                     'n_masked':recon.upper().count('N'),
                     'recon_len':len(recon),'ties':len(ties),'tie_labels':';'.join(tie_labels)})
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--sample',required=True);p.add_argument('--indir',required=True);p.add_argument('--fold',type=int,required=True)
    p.add_argument('--root',default='/home/leechuck/hla/spechla-pg');p.add_argument('--threads',type=int,default=1);p.add_argument('--genes',default=','.join(GENES))
    a=p.parse_args();root=Path(a.root)
    cat=load_catalogue(root/'source'/'catalogue_slim.tsv');exons=load_exons(root/'source'/'exons.tsv')
    result={};details=[]
    for gene in a.genes.split(','):
        records=filter_full_length(read_fasta(root/'db'/f'fold{a.fold}'/'HLA'/'whole'/f'HLA_{gene}.fasta'))
        prior=frequency_prior(records,cat)
        rows=designate_gene(a.indir,gene,records,cat,prior,exons)
        details+=rows
        result[gene]=[r['label'] if r['winner'] else '-' for r in rows]
        print(gene,[(r['label'],r['exon_edits'],r['gene_edits'],r['ties']) for r in rows],flush=True)
    with open(Path(a.indir)/'hla.result.pg.txt','w') as f:
        f.write('# database: IPD-IMGT/HLA 3.65.0 genomic + panel fold%d; rank = exon edits, gene edits, label, panel frequency\n'%a.fold)
        f.write('Sample\t'+'\t'.join(f'HLA_{g}_{i}' for g in GENES for i in (1,2))+'\n')
        f.write(a.sample+'\t'+'\t'.join(result.get(g,['-','-'])[i] for g in GENES for i in (0,1))+'\n')
    with open(Path(a.indir)/'hla.result.pg.details.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');w.writeheader();w.writerows(details)

if __name__=='__main__':main()
