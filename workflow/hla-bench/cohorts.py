#!/usr/bin/env python3
"""Freeze benchmark cohorts from repo metadata only (no predictions are read).

test_A_loo      panel donors with an assembly in the panel and a NYGC 30x CRAM
test_B_gourraud Gourraud 2014 donors with a 30x CRAM, not in the panel, no panel donor in their pedigree family
leak_gourraud   Gourraud donors with a 30x CRAM that share a pedigree family with a panel donor (not in panel)
dev             the 40 hla-asian50 pilot donors (inspected; excluded from primary test statistics)
"""
import csv,json,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent/'source'

def main():
    ped={}
    for r in csv.DictReader(open(P/'hla-pilot/source/1000G.ped'),delimiter=' '):ped[r['SampleID']]=r
    cram={}
    for line in open(P/'hla-targeted/source/C4Investigator/resources/1000Genomes_resources/tgp_full_30x.tsv'):
        s,path=line.rstrip('\n').split('\t')[:2];cram[s]=path
    panel={r['donor_id'] for r in csv.DictReader(open(P/'hla-audit/2026-09-16/panel_manifest.tsv'),delimiter='\t') if r['cohort']!='REF'}
    panel_fams={ped[d]['FamilyID'] for d in panel if d in ped}
    gourraud=[l.split()[0].strip('"') for l in open(P/'1000g_ground_truth/1000G_2014/20140702_hla_diversity.txt')][1:]
    dev={r['donor'] for r in csv.DictReader(open(P/'hla-asian50/run/source/donors.tsv'),delimiter='\t')}
    rows=[]
    def add(cohort,d):
        r=ped[d];path=cram[d]
        local=f'/usr/local/shared_data/public-human-genomes/GRCh38/1000Genomes/CRAM/{d}/{d}.cram' if path.startswith('1000G_2504_high_coverage/data/') else 'https://s3.amazonaws.com/1000genomes/'+path
        rows.append(dict(cohort=cohort,donor=d,family=r['FamilyID'],population=r['Population'],superpopulation=r['Superpopulation'],
                         dev=int(d in dev),gourraud=int(d in gourraud),cram=local))
    for d in sorted(panel):
        if d in cram and d in ped:add('test_A_loo',d)
    for d in sorted(set(gourraud)):
        if d not in cram or d not in ped or d in panel:continue
        add('leak_gourraud' if ped[d]['FamilyID'] in panel_fams else 'test_B_gourraud',d)
    OUT.mkdir(exist_ok=True)
    with open(OUT/'cohorts.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    counts={}
    for r in rows:
        c=counts.setdefault(r['cohort'],{});c[r['superpopulation']]=c.get(r['superpopulation'],0)+1;c['ALL']=c.get('ALL',0)+1
        if r['dev']:c['dev']=c.get('dev',0)+1
    summary=dict(counts=counts,sha256=hashlib.sha256((OUT/'cohorts.tsv').read_bytes()).hexdigest())
    (OUT/'cohorts.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n');print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=='__main__':main()
