#!/usr/bin/env python3
"""Deploy isolated DogoHLA validation and submit graph, preparation and donor jobs."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import shlex
import subprocess
import tarfile
from remote import ROOT, remote


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('release')
    a = p.parse_args()
    if Path(a.release).name != a.release or a.release in ('.', '..'):
        p.error('release must be a new basename')
    local = Path(__file__).parent
    root = ROOT + '/validations/' + a.release
    files = ['dogohla.py', 'experiment.py', 'phase_linkage.py', 'build_refs.py', 'graph_diagnostic.py',
             'structural_overlay.py', 'spechla_pg.sh', 'build_graph_pggb.sh', 'designate.py', 'VERSION']
    record = dict(commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                  files={n: hashlib.sha256((local/n).read_bytes()).hexdigest() for n in files})
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode='w') as tar:
        for name in files:
            tar.add(local/name, arcname='code/'+name)
        for name in ('cohort.tsv', 'PROTOCOL.md', 'frozen-method.json'):
            tar.add(local/'validation/20260918'/name, arcname=name)
        data = (json.dumps(record, indent=2)+'\n').encode()
        info = tarfile.TarInfo('code/MANIFEST.json'); info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    remote('mkdir -p '+shlex.quote(ROOT+'/validations')+' && mkdir '+shlex.quote(root)+
           ' && tar -xf - -C '+shlex.quote(root), input=archive.getvalue())
    setup = f'''from pathlib import Path
root=Path({root!r}); old=Path({ROOT!r})
for name in ('source','db'): (root/name).symlink_to(old/name)
(root/'reads').symlink_to('/home/leechuck/hla/hla-typer/reads')
(root/'scripts').symlink_to(root/'code')
(root/'logs').mkdir()
(root/'phase').mkdir(); (root/'phase/code').symlink_to(root/'code')
for fold in range(1,5):
    d=root/'graphs'/f'fold{{fold}}'; d.mkdir(parents=True)
    source=old/'graphs'/f'fold{{fold}}'
    (d/'HLA_DRB1.in.fa').symlink_to(source/'HLA_DRB1.in.fa')
    w=d/'pggb_HLA_DRB1'; w.mkdir()
    raw=list((source/'pggb_HLA_DRB1').glob('*.seqwish.gfa'))
    assert len(raw)==1, raw
    (w/raw[0].name).symlink_to(raw[0])
'''
    remote('python3 -', input=setup.encode())
    base = ['sbatch', '--parsable', '--partition=asianhla-c32', '--account=asianhla-group']
    python = '/home/leechuck/hla/mm/envs/asian50-spechla/bin/python3'
    def submit(name, command, cpus, mem, duration, extra=()):
        args = base + ['--job-name='+name, '--cpus-per-task='+str(cpus), '--mem='+mem,
                       '--time='+duration, '--output='+root+'/logs/'+name+'-%A_%a.log'] + list(extra) + ['--wrap='+command]
        return remote(shlex.join(args), capture_output=True, text=True).stdout.strip().split(';')[0]
    prep = submit('dogo_prepare', shlex.join([python, root+'/code/experiment.py', '--root', root,
                  '--release', root+'/phase', 'prepare', '--folds', '1,2,3,4']), 1, '4G', '00:20:00')
    graphs = submit('dogo_graph', 'export SPECHLA_PG_ROOT='+shlex.quote(root)+'; bash '+shlex.quote(root+'/code/build_graph_pggb.sh')+
                    ' "$SLURM_ARRAY_TASK_ID" DRB1 4', 4, '24G', '00:45:00', ['--array=1-4'])
    donors = submit('dogohla32', shlex.join([python, root+'/code/dogohla.py', '--root', root, '--paired',
                     '--threads', '4', '--cohort-index'])+' "$SLURM_ARRAY_TASK_ID"', 4, '24G', '01:30:00',
                     ['--array=0-31', '--dependency=afterok:'+prep+':'+graphs])
    record.update(root=root, jobs=dict(prepare=prep, graphs=graphs, donors=donors))
    (local/'validation/20260918/launch.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__': main()
