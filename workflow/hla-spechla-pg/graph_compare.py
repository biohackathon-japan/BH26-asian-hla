#!/usr/bin/env python3
"""Summarise graph builds (results/graphs/*.timing.tsv, *.stats.txt fetched from the cluster) into results/graph_builds.tsv."""
import csv,re,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
rows=[]
for t in sorted((HERE/'results/graphs').glob('*.timing.tsv')):
    r=list(csv.DictReader(open(t),delimiter='\t'))[0]
    stats=t.with_name(t.name.replace('.timing.tsv','.stats.txt'))
    s={}
    if stats.exists():
        for line in open(stats):
            m=re.match(r'(nodes|edges|length):\s+(\d+)',line)
            if m:s[m[1]]=int(m[2])
    r.update({'nodes':s.get('nodes',''),'edges':s.get('edges',''),'length':s.get('length','')});rows.append(r)
with open(HERE/'results/graph_builds.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=['builder','fold','gene','haplotypes','build_seconds','index_seconds','nodes','edges','length'],delimiter='\t');w.writeheader();w.writerows(rows)
for r in rows:print('\t'.join(str(r[k]) for k in ['builder','fold','gene','haplotypes','build_seconds','index_seconds','nodes','edges','length']))
