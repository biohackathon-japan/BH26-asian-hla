#!/usr/bin/env python3
"""HLA allele-name equivalence for scoring, pinned to one IPD-IMGT/HLA release.

Levels:
  two_field  exact numeric first two fields (as hla-asian50/refined/score_hla.py `norm`)
  g_group    ARS-level equivalence via hla_nom_g.txt. A name of any resolution expands to
             the G groups of every current allele it prefixes (field-boundary prefix).
A predicted slot is a list of alternatives (a tool's ambiguity list). It matches a truth
slot only if EVERY alternative is compatible, so ambiguity is never rewarded.
A truth slot is a set of acceptable values (assay ambiguity): any member may match.
"""
import re,collections
from pathlib import Path
IMGT=Path(__file__).resolve().parent.parent/'hla-analysis/source/imgt'
HISTORY=Path(__file__).resolve().parent/'source/imgt/Allelelist_history.txt'  # IPD 3.65.0, commit 5b915f27
NAME=re.compile(r'^(?:HLA-)?([A-Z0-9]+)\*(\d+(?::\d+)*)([A-Z]?)$')

def split_name(x,gene=None):
    """'HLA-A*02:01:01:01N' -> ('A', ('02','01','01','01')). Bare '02:01' needs gene."""
    x=x.strip()
    if '*' not in x:
        if gene is None or not re.fullmatch(r'\d+(?::\d+)*[A-Z]?',x):return None
        x=gene+'*'+x
    m=NAME.match(x)
    return (m[1],tuple(m[2].split(':'))) if m else None

def norm(x):
    x=x.split('*')[-1];m=re.match(r'^(\d+):(\d+)',x)
    return ':'.join(m.groups()) if m and int(m[2]) else None

class Nomenclature:
    def __init__(self,imgt=IMGT):
        self.release=next((l.split('version:')[1].strip() for l in open(Path(imgt)/'wmda/hla_nom_g.txt') if l.startswith('# version:')),'')
        self.group={}
        for line in open(Path(imgt)/'wmda/hla_nom_g.txt'):
            if line.startswith('#') or not line.strip():continue
            gene,alleles,gname=line.rstrip('\n').split(';')
            gene=gene.rstrip('*')
            for a in alleles.split('/'):
                fields=tuple(re.sub('[A-Z]$','',a).split(':'))
                self.group[gene,fields]=gname or ':'.join(fields[:3])
        self.prefix=collections.defaultdict(set)
        for (gene,fields),gname in self.group.items():
            for k in range(1,len(fields)+1):self.prefix[gene,fields[:k]].add(gname)
    def g_groups(self,name,gene=None):
        s=split_name(name,gene)
        return frozenset(self.prefix.get(s,())) if s else frozenset()
    def current_names(self,name,gene):
        """Historical name (any release since 3.0.0) -> current full names via Allelelist_history.
        Only used for names that no longer exist; deleted alleles map to nothing."""
        if not hasattr(self,'history'):
            import csv
            self.history=collections.defaultdict(set)
            rows=csv.reader(l for l in open(HISTORY) if not l.startswith('#'))
            next(rows)
            for r in rows:
                if r[1]=='NA':continue
                for old in set(r[2:]):
                    if old and old!='NA' and old!=r[1]:
                        s=split_name(old)
                        if s:
                            for k in range(1,len(s[1])+1):self.history[s[0],s[1][:k]].add(r[1])
        s=split_name(name,gene)
        return sorted(self.history.get(s,())) if s else []

def equiv(value,level,nom,gene):
    """Map one allele name to its comparison key set at `level` (empty set = unparseable)."""
    if level=='two_field':
        n=norm(value);return frozenset([n]) if n else frozenset()
    if level=='g_group':return nom.g_groups(value,gene)
    raise ValueError(level)

def truth_slot(text,level,nom,gene):
    """'32:01:01/32:01:02' -> union of keys over members. A member whose name no longer
    exists is translated through Allelelist_history; deleted members are dropped.
    None if no member can be mapped."""
    keys=[]
    for v in text.split('/'):
        if nom.g_groups(v,gene):keys.append(equiv(v,level,nom,gene));continue
        for cur in nom.current_names(v,gene):keys.append(equiv(cur,level,nom,gene))
    keys=[k for k in keys if k]
    return frozenset().union(*keys) if keys else None

def pred_slot(text,level,nom,gene):
    """'A*01:01,A*01:02' -> list of alternative key sets; None for no-call/unmappable."""
    if not text or text in {'.','NO_CALL'}:return None
    alts=[equiv(v,level,nom,gene) for v in re.split('[,;]',text) if v]
    return alts if alts and all(alts) else None

def compare(pair,truth):
    """pair: two pred slots (list of key sets or None); truth: two key sets.
    Returns (called, correct, allele_matches) like hla-asian50/refined/score_hla.py."""
    def ok(p,t):return p is not None and all(a & t for a in p)
    called=all(p is not None for p in pair)
    orders=[pair,pair[::-1]]
    correct=bool(called and any(ok(p0,truth[0]) and ok(p1,truth[1]) for p0,p1 in orders))
    matches=max(ok(p0,truth[0])+ok(p1,truth[1]) for p0,p1 in orders)
    return int(called),int(correct),matches
