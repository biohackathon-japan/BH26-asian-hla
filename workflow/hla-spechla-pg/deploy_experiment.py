#!/usr/bin/env python3
"""Deploy a new immutable code release and submit reference preparation to Slurm."""
import argparse
import io
import hashlib
import json
import shlex
from pathlib import Path
import tarfile
from remote import ROOT, remote

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('release',help='New release basename; existing directories are refused')
p.add_argument('--folds',default='0')
p.add_argument('--hybrid-phase-release',help='Completed phase-release basename; submit the hybrid pilot instead of preparing phase databases')
a=p.parse_args()
if Path(a.release).name != a.release or a.release in ('.','..'):
    p.error('release must be one directory name')
folds=','.join(str(int(x)) for x in a.folds.split(','))
release=ROOT+'/experiments/'+a.release
files=['experiment.py','phase_linkage.py','build_refs.py','graph_diagnostic.py','structural_overlay.py',
       'jobs/phase_experiment.sbatch','jobs/hybrid_experiment.sbatch']
buf=io.BytesIO()
with tarfile.open(fileobj=buf,mode='w') as tar:
    for name in files:
        tar.add(Path(__file__).parent/name,arcname='code/'+name)
    hashes={n:hashlib.sha256((Path(__file__).parent/n).read_bytes()).hexdigest() for n in files}
    data=(json.dumps(hashes,indent=2)+'\n').encode()
    info=tarfile.TarInfo('code/MANIFEST.json'); info.size=len(data)
    tar.addfile(info,io.BytesIO(data))
if a.hybrid_phase_release and (Path(a.hybrid_phase_release).name!=a.hybrid_phase_release or a.hybrid_phase_release in ('.','..')):
    p.error('hybrid phase release must be one directory name')
remote('mkdir -p '+shlex.quote(ROOT+'/experiments')+' && mkdir '+shlex.quote(release)+' && tar -xf - -C '+shlex.quote(release),input=buf.getvalue())
if a.hybrid_phase_release:
    phase=ROOT+'/experiments/'+a.hybrid_phase_release
    r=remote(shlex.join(['sbatch','--parsable','--array=0-7','--export=ALL,RELEASE='+release+',PHASE_RELEASE='+phase,
                        release+'/code/jobs/hybrid_experiment.sbatch']),capture_output=True,text=True)
    print('release='+release+'\nhybrid_job='+r.stdout.strip())
    raise SystemExit(0)
cmd=shlex.join(['/home/leechuck/hla/mm/envs/asian50-spechla/bin/python3',release+'/code/experiment.py','--release',release,'prepare','--folds',folds])
r=remote(shlex.join(['sbatch','--parsable','--partition=asianhla-c32','--account=asianhla-group','--cpus-per-task=1','--mem=4G','--time=00:15:00','--job-name=spg_prepare','--output='+ROOT+'/logs/prepare-%j.log','--wrap='+cmd]),capture_output=True,text=True)
print('release='+release+'\nprepare_job='+r.stdout.strip())
