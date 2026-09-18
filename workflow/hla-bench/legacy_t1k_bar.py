#!/usr/bin/env python3
"""Comparator bar from existing results only: legacy T1K (1000G_MHC run, --preset hla, IPD of
2026-09-15, MHC CRAM reads incl. HLA-* contigs, no unmapped reads) vs Gourraud 2014.
Cohorts from source/cohorts.tsv. No new tool runs. Quality > 0 rule as in hla-asian50."""
import csv,tarfile,collections,json
from pathlib import Path
from names import Nomenclature,pred_slot,compare
from truth import gourraud,GOURRAUD_GENES
B=Path(__file__).resolve().parent

def t1k_calls(tar_path):
    calls={}
    with tarfile.open(tar_path) as tar:
        for m in tar.getmembers():
            donor=m.name.split('/')[0]
            for line in tar.extractfile(m).read().decode().splitlines():
                v=line.split('\t');g=v[0].removeprefix('HLA-')
                if g not in GOURRAUD_GENES:continue
                n=int(v[1]);a1=v[2] if n>=1 and float(v[4])>0 else None
                a2=(a1 if n==1 else v[5] if n==2 and float(v[7])>0 else None)
                calls[donor,g]=(a1,a2) if n else (None,None)
    return calls

def main():
    nom=Nomenclature();cohorts=list(csv.DictReader(open(B/'source/cohorts.tsv'),delimiter='\t'))
    calls=t1k_calls(B/'source/t1k_legacy/t1k_genotypes.tar.gz')
    out=[];summary=collections.defaultdict(lambda:collections.Counter())
    for level in ['two_field','g_group']:
        truth,reasons=gourraud(level,nom)
        for c in cohorts:
            if not c['gourraud'] or c['gourraud']=='0':continue
            for g in GOURRAUD_GENES:
                k=(c['donor'],g)
                if k not in truth:
                    for grp in [c['superpopulation'],'ALL']:summary[level,c['cohort'],grp]['ineligible']+=1
                    continue
                if k[0] not in {d for d,_ in calls}:
                    for grp in [c['superpopulation'],'ALL']:summary[level,c['cohort'],grp]['no_t1k_output']+=1
                    continue
                raw=calls.get(k,(None,None))
                pair=[pred_slot(x,level,nom,g) if x else None for x in raw]
                called,correct,matches=compare(pair,truth[k])
                out.append(dict(level=level,cohort=c['cohort'],donor=c['donor'],superpopulation=c['superpopulation'],gene=g,prediction=';'.join(x or 'NO_CALL' for x in raw),called=called,correct=correct,allele_matches=matches))
                for grp in [c['superpopulation'],'ALL']:
                    s=summary[level,c['cohort'],grp];s['genotypes']+=1;s['called']+=called;s['correct']+=correct;s['alleles_correct']+=matches
    with open(B/'results/legacy_t1k_scores.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
    rows=[]
    for (level,cohort,grp),s in sorted(summary.items()):
        n=s['genotypes'];rows.append(dict(level=level,cohort=cohort,superpopulation=grp,genotypes=n,ineligible=s['ineligible'],no_t1k_output=s['no_t1k_output'],
            call_rate_pct=round(100*s['called']/n,2) if n else '',genotype_accuracy_pct=round(100*s['correct']/n,2) if n else '',allele_accuracy_pct=round(50*s['alleles_correct']/n,2) if n else ''))
    with open(B/'results/legacy_t1k_summary.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    for r in rows:print(*r.values(),sep='\t')
if __name__=='__main__':main()
