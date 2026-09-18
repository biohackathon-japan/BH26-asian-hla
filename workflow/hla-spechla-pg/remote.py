#!/usr/bin/env python3
"""DDBJ SSH transport (copied from hla-asian50/run/remote.py with a new root)."""
import argparse,shlex,subprocess
from pathlib import Path
ROOT='/home/leechuck/hla/spechla-pg'
def remote(command,**kwargs):
    hop='ssh -o ConnectTimeout=20 -o BatchMode=yes -o IdentitiesOnly=yes -i ~/.ssh/id_ed25519 a001 '+shlex.quote(command)
    return subprocess.run(['ssh','-o','ConnectTimeout=20','-o','BatchMode=yes','ddbj',hop],check=True,**kwargs)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',nargs='?');p.add_argument('--upload');p.add_argument('--download');p.add_argument('--dest');a=p.parse_args()
    if a.upload:
        dest=a.dest or ROOT+'/'+a.upload
        remote('mkdir -p '+shlex.quote(str(Path(dest).parent)))
        remote('cat > '+shlex.quote(dest),input=Path(a.upload).read_bytes())
    elif a.download:
        r=remote('cat '+shlex.quote(a.download),capture_output=True);Path(a.dest).write_bytes(r.stdout)
    else:remote(a.command)
