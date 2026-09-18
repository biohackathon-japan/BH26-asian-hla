#!/usr/bin/env python3
"""DōgoHLA: run the frozen AsianHLA/SpecHLA reconstruction on a prepared root."""
import argparse
import csv
from pathlib import Path
import shutil
import time
from experiment import ENV, command, environment, failed, finish, phase_run, sha, start, validate_outputs
from graph_diagnostic import diagnose
from structural_overlay import run as overlay


def checked_legacy(root, donor, fold, threads, native=False, long_indels=False):
    arm = 'native-long' if long_indels else ('native' if native else 'A')
    out = root / 'runs' / donor / arm
    reads = root / 'reads' / donor
    if not (reads / 'COMPLETE').exists():
        raise ValueError(f'Incomplete read recruitment: {reads}')
    inputs = [reads / f'r{i}.fq.gz' for i in (1, 2)]
    script = ENV / 'share/spechla/script'
    config = dict(donor=donor, fold=fold, arm=arm, threads=threads,
                  inputs={str(p): sha(p) for p in inputs}, driver_sha256=sha(__file__),
                  pipeline_sha256=sha(script / 'whole/SpecHLA.sh') if native else sha(root / 'code/spechla_pg.sh'))
    if not start(out, config):
        return
    env = environment(script)
    env['SPECHLA_PG_ROOT'] = str(root)
    t0 = time.monotonic()
    try:
        if native:
            tmp = root / 'runs' / donor / (arm + '_work')
            args = ['timeout', '--kill-after=30s', '3600', 'spechla', '-n', donor,
                    '-1', inputs[0], '-2', inputs[1], '-o', tmp, '-j', threads, '-u', '0', '-p', 'Unknown']
            if long_indels:
                args += ['-v', 'True']
            command(args, out / 'spechla.log', env)
            for p in (tmp / donor).iterdir():
                shutil.move(str(p), out / p.name)
        else:
            command(['timeout', '--kill-after=30s', '3600', 'bash', root / 'code/spechla_pg.sh',
                     '-n', donor, '-1', inputs[0], '-2', inputs[1], '-o', out, '-k', fold,
                     '-m', 'A', '-j', threads], out / 'driver.log', env)
            # Replace the legacy unstructured marker only after validating outputs.
            (out / 'COMPLETE').unlink(missing_ok=True)
        for log in out.rglob('*.log'):
            text = log.read_text(errors='replace')
            if any(x in text for x in ('Traceback (most recent call last)', 'command not found', 'Segmentation fault')):
                raise ValueError(f'Nested command failure: {log}')
        outputs = validate_outputs(out, donor)
        (out / 'timing.tsv').write_text(f'stage\tseconds\ntotal\t{time.monotonic()-t0:.3f}\n')
        finish(out, config, outputs)
    except BaseException as exc:
        (out / 'COMPLETE').unlink(missing_ok=True)
        failed(out, exc)
        raise


def run(root, donor, fold, threads):
    checked_legacy(root, donor, fold, threads)
    phase_run(root, root / 'phase', donor, fold, 'AE-ipd365-pairing', threads)
    source = root / 'phase/runs' / donor / 'AE-ipd365-pairing'
    graph = root / 'diagnostics' / donor
    overlay(source, graph, donor, root / 'final/runs' / donor / 'DogoHLA-no-graph', [], threads)
    diagnose(root, graph / 'DRB1', fold, 'DRB1', donor, 'pggb', 1000000, threads, 300)
    overlay(source, graph, donor, root / 'final/runs' / donor / 'DogoHLA', ['DRB1'], threads)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version='DogoHLA 0.1.0')
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--donor')
    p.add_argument('--fold', type=int)
    p.add_argument('--cohort-index', type=int)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--native-only', action='store_true')
    p.add_argument('--long-indels', action='store_true')
    p.add_argument('--paired', action='store_true', help='Also run native SpecHLA, retaining either arm on failure')
    a = p.parse_args()
    root = a.root.resolve()
    if a.cohort_index is not None:
        rows = list(csv.DictReader(open(root / 'cohort.tsv'), delimiter='\t'))
        row = rows[a.cohort_index]
        a.donor, a.fold = row['donor'], int(row['fold'])
    if not a.donor or a.fold is None or a.threads <= 0:
        p.error('Provide donor/fold (or cohort-index) and positive threads')
    errors = []
    if a.native_only or a.paired:
        try:
            checked_legacy(root, a.donor, a.fold, a.threads, native=True, long_indels=a.long_indels)
        except Exception as exc:
            errors.append(('native', repr(exc)))
    if not a.native_only:
        try:
            run(root, a.donor, a.fold, a.threads)
        except Exception as exc:
            errors.append(('DogoHLA', repr(exc)))
    if errors:
        raise RuntimeError(errors)


if __name__ == '__main__':
    main()
