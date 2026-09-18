#!/usr/bin/env python3
"""Score SpecHLA-PG and native SpecHLA predictions on the development donors (runs locally on fetched outputs).

Inputs: results/runs/<donor>/<mode>/{hla.result.txt, hla.result.g.group.txt, hla.result.pg.txt, hla.allele.*.fasta}
Truth A: the donor's two catalogue haplotypes (hla-analysis/results/sequence_catalogue.tsv, IPD-IMGT/HLA 3.65.0)
Truth B: Gourraud 2014 experimental typings (hla-bench/truth.py)
Scoring: hla-bench/names.py (compare; unordered pair; no-call = error) at two_field and g_group; plus an exact
three-field (CDS) level and, for sequences, the edit distance of the truth gene body inside each reconstruction.
Outputs: results/predictions.tsv, results/scores.tsv, results/summary.tsv, results/sequence_scores.tsv,
results/errors.tsv (one row per wrong or missing call, with the designation evidence behind it)
"""
import argparse,collections,csv,itertools,re,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parent
sys.path.insert(0,str(REPO/'hla-bench'))
from names import Nomenclature,pred_slot,truth_slot,compare,split_name
GENES=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1']
FLANK=2000
# method = (mode, result file, designator label)
METHODS=[('native','hla.result.txt','F-native'),('native','hla.result.g.group.txt','F-native-G'),('native','hla.result.pg.txt','F-pg'),
         ('A','hla.result.txt','F-native'),('A','hla.result.g.group.txt','F-native-G'),('A','hla.result.pg.txt','F-pg'),
         ('AD','hla.result.txt','F-native'),('AD','hla.result.g.group.txt','F-native-G'),('AD','hla.result.pg.txt','F-pg'),
         ('ADEF','hla.result.txt','F-native'),('ADEF','hla.result.g.group.txt','F-native-G'),('ADEF','hla.result.pg.txt','F-pg')]

def read_result(path):
    """SpecHLA result layout -> {gene: [allele1, allele2]} ('-' / 'no_match' / '' -> None)."""
    lines=[l.rstrip('\n') for l in open(path) if l.strip() and not l.startswith('#')]
    if len(lines)<2:return {}
    head=lines[0].split('\t');vals=lines[1].split('\t');row=dict(zip(head,vals));out={}
    for g in GENES:
        pair=[]
        for i in (1,2):
            v=row.get(f'HLA_{g}_{i}','').strip()
            pair.append(None if v in ('','-','no_match') or v.startswith('novel:') else v)
        out[g]=pair
    return out

def three_field(name):
    s=split_name(name)
    return s[0]+'*'+':'.join(s[1][:3]) if s and len(s[1])>=3 else None

def read_fasta(path):
    out=[];name=None;buf=[]
    for line in open(path):
        if line.startswith('>'):
            if name is not None:out.append((name,''.join(buf)))
            name=line[1:].split()[0];buf=[]
        else:buf.append(line.strip())
    if name is not None:out.append((name,''.join(buf)))
    return out

def read_details(path):
    """hla.result.pg.details.tsv -> {(gene,hap): row}; empty when the pangenome designator did not run."""
    if not Path(path).exists():return {}
    return {(r['gene'],r['hap']):r for r in csv.DictReader(open(path),delimiter='\t')}

def truth_a(catalogue,donors):
    """{(donor,gene): [record1, record2]} with the catalogue rows of the donor's two haplotypes."""
    t=collections.defaultdict(dict)
    for r in csv.DictReader(open(catalogue),delimiter='\t'):
        sample,hap,g=r['name'].split('#');gene=g.replace('HLA-','')
        if sample in donors and gene in GENES:t[(sample,gene)][hap]=r
    return {k:[v.get('1'),v.get('2')] for k,v in t.items() if '1' in v and '2' in v}

def truth_slots(recs,level,nom,gene):
    """Truth slots from exact_cds_alleles ('/'-joined for names.truth_slot); None if a haplotype has no exact CDS allele."""
    slots=[]
    for r in recs:
        names=[x for x in r['exact_cds_alleles'].split(';') if x]
        if not names:return None
        if level=='three_field':
            s=frozenset(t for t in (three_field(n) for n in names) if t)
            slots.append(s or None)
        else:slots.append(truth_slot('/'.join(names),level,nom,gene))
    return None if None in slots else tuple(slots)

def pred_pair(pair,level,nom,gene):
    if level=='three_field':return [[frozenset([three_field(v)])] if v and three_field(v) else None for v in pair]
    return [pred_slot(v,level,nom,gene) if v else None for v in pair]

def sequence_score(run_dir,gene,recs,panel_seqs):
    """Best assignment of the two reconstructions to the two truth gene bodies: (edits1, edits2, exact, n_masked)."""
    import edlib
    recon=[]
    for hap in (1,2):
        fa=run_dir/f'hla.allele.{hap}.HLA_{gene}.fasta'
        recon.append(''.join(s for _,s in read_fasta(fa)).upper() if fa.exists() else '')
    truths=[panel_seqs[r['name']][FLANK:-FLANK].upper() for r in recs]
    def d(t,q):return edlib.align(t,q,mode='HW',task='distance')['editDistance'] if q else len(t)
    best=min(((d(truths[0],recon[i]),d(truths[1],recon[j])) for i,j in ((0,1),(1,0))),key=sum)
    return best[0],best[1],int(best[0]==0)+int(best[1]==0),sum(r.count('N') for r in recon)

def main():
    p=argparse.ArgumentParser();p.add_argument('--runs',default=str(HERE/'results/runs'));p.add_argument('--out',default=str(HERE/'results'))
    p.add_argument('--catalogue',default=str(REPO/'hla-analysis/results/sequence_catalogue.tsv'));p.add_argument('--panel',default=str(REPO/'hla-analysis/source'))
    a=p.parse_args();runs=Path(a.runs);out=Path(a.out);out.mkdir(exist_ok=True)
    nom=Nomenclature();donors=sorted(d.name for d in runs.iterdir() if d.is_dir())
    ta=truth_a(a.catalogue,set(donors))
    import truth as T;tb={lvl:T.gourraud(lvl,nom)[0] for lvl in ('two_field','g_group')}
    panel={g:dict(read_fasta(Path(a.panel)/f'HLA-{g}.fa')) for g in GENES}
    preds=[];scores=[];seqs=[];details={};pred_by={}
    for donor in donors:
        for mode,fname,desig in METHODS:
            f=runs/donor/mode/fname
            if not f.exists():continue
            res=read_result(f);method=f'{mode}/{desig}'
            if (donor,mode) not in details:details[(donor,mode)]=read_details(runs/donor/mode/'hla.result.pg.details.tsv')
            for g in GENES:
                pair=res.get(g,[None,None])
                preds.append({'donor':donor,'method':method,'gene':g,'allele1':pair[0] or 'NO_CALL','allele2':pair[1] or 'NO_CALL'})
                pred_by[(donor,method,g)]=pair
                for lvl in ('two_field','g_group','three_field'):
                    if (donor,g) in ta:
                        ts=truth_slots(ta[(donor,g)],lvl,nom,g)
                        if ts is not None:
                            called,correct,m=compare(pred_pair(pair,lvl,nom,g),ts)
                            scores.append({'donor':donor,'method':method,'gene':g,'truth':'A','level':lvl,'called':called,'correct':correct,'allele_matches':m,
                                           'truth1':'/'.join(sorted(ts[0]))[:60],'truth2':'/'.join(sorted(ts[1]))[:60]})
                    if lvl!='three_field' and (donor,g) in tb[lvl]:
                        called,correct,m=compare(pred_pair(pair,lvl,nom,g),tb[lvl][(donor,g)])
                        scores.append({'donor':donor,'method':method,'gene':g,'truth':'B','level':lvl,'called':called,'correct':correct,'allele_matches':m,
                                       'truth1':'/'.join(sorted(tb[lvl][(donor,g)][0]))[:60],'truth2':'/'.join(sorted(tb[lvl][(donor,g)][1]))[:60]})
        for mode in ('native','A','AD','ADEF'):
            if not (runs/donor/mode).exists():continue
            for g in GENES:
                if (donor,g) not in ta:continue
                e1,e2,exact,nmask=sequence_score(runs/donor/mode,g,ta[(donor,g)],panel[g])
                seqs.append({'donor':donor,'mode':mode,'gene':g,'edits_hap1':e1,'edits_hap2':e2,'exact_haplotypes':exact,'n_masked':nmask})
    # every wrong or missing call, with the designation evidence that produced it
    errors=[]
    for s in scores:
        if s['correct']:continue
        donor,method,g=s['donor'],s['method'],s['gene'];mode=method.split('/')[0]
        pair=pred_by.get((donor,method,g),[None,None]);det=details.get((donor,mode),{})
        seq=[x for x in seqs if x['donor']==donor and x['mode']==mode and x['gene']==g]
        errors.append({'donor':donor,'method':method,'gene':g,'truth_level':s['level'],'truth_source':s['truth'],
                       'truth1':s['truth1'],'truth2':s['truth2'],'pred1':pair[0] or 'NO_CALL','pred2':pair[1] or 'NO_CALL',
                       'matched_alleles':s['allele_matches'],
                       'pg_winner1':det.get((g,'1'),{}).get('winner',''),'pg_winner2':det.get((g,'2'),{}).get('winner',''),
                       'pg_label_source1':det.get((g,'1'),{}).get('label_source',''),'pg_label_source2':det.get((g,'2'),{}).get('label_source',''),
                       'exon_edits1':det.get((g,'1'),{}).get('exon_edits',''),'exon_edits2':det.get((g,'2'),{}).get('exon_edits',''),
                       'n_masked':seq[0]['n_masked'] if seq else '','seq_edits':f"{seq[0]['edits_hap1']}/{seq[0]['edits_hap2']}" if seq else ''})
    with open(out/'errors.tsv','w') as f:
        cols=['donor','method','gene','truth_level','truth_source','truth1','truth2','pred1','pred2','matched_alleles',
              'pg_winner1','pg_winner2','pg_label_source1','pg_label_source2','exon_edits1','exon_edits2','n_masked','seq_edits']
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(errors)
    with open(out/'predictions.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=['donor','method','gene','allele1','allele2'],delimiter='\t');w.writeheader();w.writerows(preds)
    with open(out/'scores.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=['donor','method','gene','truth','level','called','correct','allele_matches','truth1','truth2'],delimiter='\t');w.writeheader();w.writerows(scores)
    with open(out/'sequence_scores.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=['donor','mode','gene','edits_hap1','edits_hap2','exact_haplotypes','n_masked'],delimiter='\t');w.writeheader();w.writerows(seqs)
    # summary: per method x truth x level, per gene and pooled: correct/n
    agg=collections.defaultdict(lambda:[0,0,0]);rows=[]
    for s in scores:
        for gene in (s['gene'],'ALL'):
            k=(s['method'],s['truth'],s['level'],gene);agg[k][0]+=s['correct'];agg[k][1]+=1;agg[k][2]+=s['called']
    for (method,truth,level,gene),(c,n,called) in sorted(agg.items()):
        rows.append({'method':method,'truth':truth,'level':level,'gene':gene,'correct':c,'n':n,'called':called,'accuracy':f'{c/n:.4f}'})
    sagg=collections.defaultdict(lambda:[0,0,0,0])
    for s in seqs:
        for gene in (s['gene'],'ALL'):
            k=(s['mode'],gene);sagg[k][0]+=s['exact_haplotypes'];sagg[k][1]+=2;sagg[k][2]+=s['edits_hap1']+s['edits_hap2'];sagg[k][3]+=int(s['exact_haplotypes']==2)
    for (mode,gene),(ex,n,ed,both) in sorted(sagg.items()):
        rows.append({'method':f'{mode}/sequence','truth':'A','level':'exact_gene_haplotypes','gene':gene,'correct':ex,'n':n,'called':both,'accuracy':f'{ex/n:.4f}'})
        rows.append({'method':f'{mode}/sequence','truth':'A','level':'mean_edit_distance','gene':gene,'correct':ed,'n':n,'called':both,'accuracy':f'{ed/n:.1f}'})
    with open(out/'summary.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=['method','truth','level','gene','correct','n','called','accuracy'],delimiter='\t');w.writeheader();w.writerows(rows)
    for r in rows:
        if r['gene']=='ALL':print('\t'.join(str(r[k]) for k in ('method','truth','level','correct','n','accuracy')))

if __name__=='__main__':main()
