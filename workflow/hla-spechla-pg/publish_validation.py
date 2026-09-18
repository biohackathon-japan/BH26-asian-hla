#!/usr/bin/env python3
"""Publish only the final compact validation report, refusing unrelated staged work."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parent.parent
REPORT='hla-spechla-pg/results/validation-20260918'

def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()


def main():
    for attempt in range(31):
        record=json.loads((ROOT/REPORT/'monitor.json').read_text())
        if record['complete'] or record['deadline_reached']:
            break
        if attempt < 30:
            time.sleep(10)
    if not record['complete'] and not record['deadline_reached']:
        raise RuntimeError('Monitor has not produced a terminal report')
    if git('branch','--show-current') != 'hla-experiments':
        raise RuntimeError('Repository branch changed; refusing automated commit')
    if git('diff','--cached','--name-only'):
        raise RuntimeError('Unrelated staged work; refusing automated commit')
    files=[REPORT+'/monitor.json']+[str(p.relative_to(ROOT)) for p in sorted((ROOT/REPORT/'analysis').iterdir()) if p.is_file()]
    git('add','--',*files)
    git('-c','core.whitespace=cr-at-eol','diff','--cached','--check')
    if git('diff','--cached','--name-only'):
        git('commit','-m','Report frozen DogoHLA 32-donor validation outcomes')
    git('push','origin','hla-experiments')
    print(datetime.now(timezone.utc).isoformat(),git('rev-parse','HEAD'))


if __name__=='__main__': main()
