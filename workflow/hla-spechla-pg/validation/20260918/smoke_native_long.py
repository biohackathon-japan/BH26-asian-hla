#!/usr/bin/env python3
"""Feasibility gate for SpecHLA's optional ScanIndel mode on a pilot donor."""
import os
from pathlib import Path
import shutil
import socket
import sys
sys.path.insert(0, str(Path(sys.argv[1]) / 'code'))
from experiment import ENV, command, environment, json_write, sha, validate_outputs
root = Path(sys.argv[1])
out = root / 'native-long-smoke'
out.mkdir()
script = out / 'vendor/script'
shutil.copytree(ENV / 'share/spechla/script', script)
pipeline = script / 'whole/SpecHLA.sh'
text = pipeline.read_text()
old = 'port=$(date +%N|cut -c5-9)'
assert text.count(old) == 1
# Execution-only fix: upstream generates ports outside the valid TCP range.
with socket.socket() as s:
    s.bind(('127.0.0.1', 0))
    port = s.getsockname()[1]
pipeline.write_text(text.replace(old, f'port={port}'))
env = environment(script)
reads = root / 'reads/HG02155'
args = ['timeout', '--kill-after=30s', '1200', 'spechla', '-n', 'HG02155',
        '-1', reads/'r1.fq.gz', '-2', reads/'r2.fq.gz', '-o', out/'runs', '-j', '4', '-u', '0', '-p', 'Unknown', '-v', 'True']
json_write(out/'configuration.json',dict(donor='HG02155', port=port, script_sha256=sha(pipeline),
           modification='Only replace invalid random TCP port with OS-selected free valid port', command=list(map(str,args))))
try:
    command(args,out/'spechla.log',env)
    result=out/'runs/HG02155'
    for log in [out/'spechla.log', *result.rglob('*.log')]:
        content=log.read_text(errors='replace')
        if any(x in content for x in ('Traceback (most recent call last)', 'command not found', 'Segmentation fault')):
            raise ValueError(f'Nested failure: {log}')
    outputs=validate_outputs(result,'HG02155')
    json_write(out/'STATUS.json',dict(status='complete',outputs=outputs))
except BaseException as exc:
    json_write(out/'STATUS.json',dict(status='failed',error=repr(exc)))
    raise
