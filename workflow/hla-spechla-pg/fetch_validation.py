#!/usr/bin/env python3
"""Fetch a compact validation snapshot without reads, alignments or graph indexes."""
import argparse
import io
from pathlib import Path
import tarfile
from remote import remote, ROOT


def fetch(release, out):
    if Path(release).name != release or release in ('.','..'):
        raise ValueError('release must be a basename')
    root=ROOT+'/validations/'+release
    script='''import sys,tarfile
from pathlib import Path
root=Path(ROOT_VALUE)
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as tar:
 for sub in ('runs','final/runs','phase/runs','diagnostics'):
  for f in (root/sub).glob('*/*/*'):
   if f.is_symlink() or not f.is_file(): continue
   if (f.name in ('COMPLETE','timing.tsv','structural_candidates.tsv','gam.stats.txt') or
       f.suffix == '.json' or f.name.startswith('hla.allele.') and f.suffix=='.fasta' or
       f.name.startswith('hla.result.') and f.suffix=='.txt'):
    tar.add(f,arcname=str(f.relative_to(root)),recursive=False)
 for name in ('code/MANIFEST.json','phase/PREPARED.json','cohort.tsv','frozen-method.json','native-long-smoke/STATUS.json','native-long-smoke/configuration.json'):
  f=root/name
  if f.is_file(): tar.add(f,arcname=name,recursive=False)
'''.replace('ROOT_VALUE',repr(root))
    r=remote('python3 -',input=script.encode(),capture_output=True)
    out.mkdir(parents=True,exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(r.stdout),mode='r:gz') as tar:
        tar.extractall(out,filter='data')
    print(out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('release'); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); fetch(a.release,a.out)
