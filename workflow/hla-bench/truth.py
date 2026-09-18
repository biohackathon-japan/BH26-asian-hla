#!/usr/bin/env python3
"""Truth loaders. Eligibility is decided from truth alone, never from predictions."""
from pathlib import Path
from names import truth_slot
P=Path(__file__).resolve().parent.parent
GOURRAUD=P/'1000g_ground_truth/1000G_2014/20140702_hla_diversity.txt'
GOURRAUD_GENES=['A','B','C','DRB1','DQB1']

def gourraud(level,nom,path=GOURRAUD):
    """{(donor,gene): (slot1,slot2)} for eligible genotypes, plus {(donor,gene): reason} for ineligible ones.
    A member needs at least two numeric fields; '0000' and one-field typings are ineligible."""
    lines=[l.split() for l in open(path)];head=[h.strip('"') for h in lines[0]]
    truth={};reasons={}
    for v in lines[1:]:
        r=dict(zip(head,[x.strip('"') for x in v]))
        for g in GOURRAUD_GENES:
            raw=[r[g],r[g+'.1']];k=(r['id'],g)
            if any(m.count(':')<1 for x in raw for m in x.split('/')):reasons[k]='missing_or_one_field';continue
            slots=[truth_slot(x,level,nom,g) for x in raw]
            if None in slots:reasons[k]='unmappable_name';continue
            truth[k]=tuple(slots)
    return truth,reasons
