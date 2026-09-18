#!/usr/bin/env python3
"""Fetch compact experiment results and logs, never reads/BAM/GAM/databases."""
import argparse
import io
import json
from pathlib import Path
import tarfile
from remote import remote, ROOT

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('release',help='Experiment directory basename')
p.add_argument('--out',type=Path)
a=p.parse_args()
if Path(a.release).name != a.release or a.release in ('.','..'):
    p.error('release must be one directory name')
out=a.out or Path(__file__).parent/'results/experiments'/a.release
root=ROOT+'/experiments/'+a.release
script='''import sys,tarfile
from pathlib import Path
root=Path(ROOT_VALUE)
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as tar:
 for sub in [root/'PREPARED.json',root/'code',root/'runs',*root.glob('diagnostics*')]:
  paths=[sub] if sub.is_file() else sub.rglob('*')
  for f in paths:
   if f.is_symlink() or not f.is_file(): continue
   if sub.name=='runs' and len(f.relative_to(sub).parts)!=3: continue
   if sub.name=='code' and f.suffix in ('.py','.sbatch','.json'):
    tar.add(f,arcname=str(f.relative_to(root)),recursive=False)
    continue
   if (f.name in ('COMPLETE','timing.tsv','structural_candidates.tsv') or f.suffix in ('.json','.log')
       or (f.name.startswith('hla.allele.') and f.suffix=='.fasta')
       or (f.name.startswith('hla.result.') and f.suffix=='.txt')
       or f.name=='gam.stats.txt'):
    tar.add(f,arcname=str(f.relative_to(root)),recursive=False)
'''.replace('ROOT_VALUE',repr(root))
r=remote('python3 -',input=script.encode(),capture_output=True)
out.mkdir(parents=True,exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(r.stdout),mode='r:gz') as tar:
    tar.extractall(out,filter='data')
print(out)
