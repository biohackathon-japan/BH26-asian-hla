#!/usr/bin/env python3
"""Prepare and execute isolated, provenance-checked SpecHLA experiments on DDBJ.

The phase arms reuse EXACTLY the A-arm BAM/VCF/read bins, changing only E.
All output lives under a new experiments/<release>/ directory.
"""
import argparse
import ast
import csv
import hashlib
import inspect
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from build_refs import GENES, excluded_haps, read_fasta, write_fasta
from phase_linkage import Linkage, read_blast

ROOT = Path('/home/leechuck/hla/spechla-pg')
ENV = Path('/home/leechuck/hla/mm/envs/asian50-spechla')
VG = Path('/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin')
ARMS = ('A-replay', 'A-fixed', 'AE-ipd365', 'AE-panel')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for data in iter(lambda: f.read(1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def json_write(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def patch_block_linker(source):
    """Replace two inspected definitions, failing closed if the API changes."""
    tree = ast.parse(source)
    linkage = next(x for x in tree.body if isinstance(x, ast.ClassDef) and x.name == 'Linkage')
    analyze = next(x for x in tree.body if isinstance(x, ast.ClassDef) and x.name == 'Analyze_map')
    reader = next(x for x in analyze.body if isinstance(x, ast.FunctionDef) and x.name == 'read_blast')
    lines = source.splitlines(keepends=True)
    replacements = [
        (linkage.lineno, linkage.end_lineno, inspect.getsource(Linkage)),
        (reader.lineno, reader.end_lineno,
         ''.join('    ' + line if line.strip() else line for line in inspect.getsource(read_blast).splitlines(keepends=True))),
    ]
    for start, end, text in sorted(replacements, reverse=True):
        lines[start - 1:end] = [text + '\n']
    result = ''.join(lines)
    compile(result, '<patched block linker>', 'exec')
    return result


def environment(script):
    env = dict(os.environ)
    env.update(PATH=f'{ENV}/bin:{VG}:' + env.get('PATH', ''),
               CONDA_PREFIX=str(ENV), LD_LIBRARY_PATH=str(ENV / 'lib'),
               SPECHLA_DB=str(ENV / 'share/spechla/db'), SPECHLA_SCRIPT=str(script),
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               PYTHONHASHSEED='0')
    return env


def command(args, log, env, timeout=None, stdout=None):
    """Every direct subprocess must succeed; logs remain available on failure."""
    with open(log, 'ab') as err:
        err.write(('COMMAND ' + repr([str(x) for x in args]) + '\n').encode())
        err.flush()
        subprocess.run([str(x) for x in args], env=env, stdout=stdout or err,
                       stderr=err, check=True, timeout=timeout)


def validate_outputs(out, sample, genes=GENES):
    paths = []
    for gene in genes:
        for hap in (1, 2):
            path = out / f'hla.allele.{hap}.HLA_{gene}.fasta'
            records = read_fasta(path)
            if len(records) != 1 or not records[0][1] or set(records[0][1].upper()) - set('ACGTN'):
                raise ValueError(f'Invalid sequence output: {path}')
            paths.append(path)
    for filename in ('hla.result.txt', 'hla.result.g.group.txt'):
        path = out / filename
        lines = [line.rstrip().split('\t') for line in path.read_text().splitlines()
                 if line.strip() and not line.startswith('#')]
        if len(lines) != 2 or len(lines[0]) != len(lines[1]) or lines[1][0] != sample:
            raise ValueError(f'Invalid result table: {path}')
        expected = {f'HLA_{gene}_{hap}' for gene in GENES for hap in (1, 2)}
        if not expected.issubset(lines[0]):
            raise ValueError(f'Missing result columns: {path}')
        paths.append(path)
    # Several upstream scripts use os.system, which can hide subprocess failure.
    for log in out.glob('phase.*.log'):
        content = log.read_text(errors='replace')
        if any(marker in content for marker in ('Traceback (most recent call last)', 'command not found',
                                               'No such file or directory', 'Segmentation fault')):
            raise ValueError(f'Nested SpecHLA command failed: {log}')
    return {str(p.relative_to(out)): sha(p) for p in paths}


def completed(out, config):
    marker = out / 'COMPLETE'
    if not marker.exists():
        return False
    record = json.loads(marker.read_text())
    if record['configuration_sha256'] != identity(config):
        raise ValueError(f'Refusing stale completion/configuration: {out}')
    for name, digest in record['outputs'].items():
        if sha(out / name) != digest:
            raise ValueError(f'Completed output changed: {out / name}')
    return True


def start(out, config):
    if completed(out, config):
        print('VERIFIED COMPLETE', out, flush=True)
        return False
    # Never silently mix a failed or partial run with a fresh attempt.
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'Partial run exists; inspect and choose a new output directory: {out}')
    out.mkdir(parents=True, exist_ok=True)
    json_write(out / 'manifest.json', dict(configuration=config, status='running', started=time.time()))
    return True


def finish(out, config, outputs):
    json_write(out / 'COMPLETE', dict(configuration_sha256=identity(config), outputs=outputs,
                                     finished=time.time()))
    record = json.loads((out / 'manifest.json').read_text())
    record.update(status='complete', finished=time.time())
    json_write(out / 'manifest.json', record)


def failed(out, exc):
    path = out / 'manifest.json'
    if path.exists():
        record = json.loads(path.read_text())
        record.update(status='failed', error=repr(exc), finished=time.time())
        json_write(path, record)


def link(source, dest):
    dest.symlink_to(source.resolve())


def prepare(root, release, folds):
    if release.exists() and any(p.name != 'code' for p in release.iterdir()):
        raise FileExistsError(f'Use a fresh release directory: {release}')
    release.mkdir(parents=True, exist_ok=True)
    installed = ENV / 'share/spechla'
    hashes = {}
    for variant in ('original', 'fixed'):
        target = release / 'vendor' / variant
        shutil.copytree(installed / 'script', target / 'script', ignore=shutil.ignore_patterns('__pycache__'))
        link(installed / 'db', target / 'db')
        link(ENV / 'bin', target / 'bin')
        block = target / 'script/whole/map_block2_database.py'
        hashes[variant + '_before'] = sha(block)
        if variant == 'fixed':
            block.write_text(patch_block_linker(block.read_text()))
        hashes[variant + '_after'] = sha(block)
        # Both arms share deterministic eigenvector initialisation, independent of the score repair.
        phaser = target / 'script/phase_unlinked_block.py'
        text = phaser.read_text()
        old = "eigsh(L, k=2, which='SM')"
        if text.count(old) != 1:
            raise ValueError('Unexpected SpecHLA eigensolver source; inspect before patching')
        phaser.write_text(text.replace(old, "eigsh(L, k=2, which='SM', v0=np.random.RandomState(20260918).normal(size=L.shape[0]))"))
    db_manifest = {}
    for fold in folds:
        excluded, _ = excluded_haps(root / 'source', fold)
        for kind in ('ipd365', 'panel'):
            db = release / 'db' / f'fold{fold}' / kind
            (db / 'HLA/whole').mkdir(parents=True)
            link(installed / 'db/ref', db / 'ref')
            for p in (installed / 'db/HLA').iterdir():
                if p.name != 'whole':
                    link(p, db / 'HLA' / p.name)
            for p in (installed / 'db/HLA/whole').iterdir():
                link(p, db / 'HLA/whole' / p.name)
            for gene in GENES:
                records = [(name.split()[1], seq.upper()) for name, seq in
                           read_fasta(root / 'source/imgt_gen' / f'{gene}_gen.fasta')]
                n_ipd = len(records)
                if kind == 'panel':
                    # Same IPD records in both arms; only training-panel gene bodies differ.
                    panel = read_fasta(root / 'db' / f'fold{fold}' / 'HLA/whole' / f'HLA_{gene}.fasta')
                    leaked = [n for n, _ in panel if '#'.join(n.split('#')[:2]) in excluded]
                    if leaked:
                        raise ValueError(f'Held-out haplotypes in phase database: {leaked}')
                    records += [(name, seq[2000:-2000]) for name, seq in panel if '#' in name]
                    if any(not seq for _, seq in records):
                        raise ValueError(f'Empty gene body after flank removal: {gene}')
                target = db / 'HLA/whole' / f'HLA_{gene}.fasta'
                if target.is_symlink():
                    target.unlink()
                write_fasta(target, records)
                db_manifest[str(target.relative_to(release))] = dict(sha256=sha(target), ipd_records=n_ipd,
                                                                     panel_records=len(records) - n_ipd)
    source_hashes = {str(p.relative_to(release / 'vendor')): sha(p)
                     for p in sorted((release / 'vendor').glob('*/script/**/*')) if p.is_file()}
    json_write(release / 'PREPARED.json', dict(block_linker=hashes, databases=db_manifest,
              vendor_files=source_hashes, software={name: sha(Path(__file__).parent / name)
              for name in ('experiment.py', 'phase_linkage.py', 'build_refs.py')}, folds=folds,
              exclusion_manifest_sha256=sha(root / 'db/build_log.json'), seed=20260918))
    print('PREPARED', release, flush=True)


def verify_prepared(root, release, donor, fold):
    record = json.loads((release / 'PREPARED.json').read_text())
    if fold not in record['folds']:
        raise ValueError(f'Fold {fold} was not prepared')
    rows = list(csv.DictReader(open(root / 'source/donors_folds.tsv'), delimiter='\t'))
    if not any(r['donor'] == donor and int(r['fold']) == fold for r in rows):
        raise ValueError(f'Incorrect donor/fold assignment: {donor}/{fold}')
    checks = [(root / 'db/build_log.json', record['exclusion_manifest_sha256'])]
    checks += [(release / 'vendor' / p, h) for p, h in record['vendor_files'].items()]
    checks += [(release / p, info['sha256']) for p, info in record['databases'].items()]
    checks += [(Path(__file__).parent / p, h) for p, h in record['software'].items()]
    for path, expected in checks:
        if sha(path) != expected:
            raise ValueError(f'Prepared input changed: {path}')


def phase_run(root, release, donor, fold, arm, threads, output=None):
    if arm not in ARMS:
        raise ValueError(arm)
    verify_prepared(root, release, donor, fold)
    original = ENV / 'share/spechla/db'
    script = release / 'vendor' / ('original' if arm == 'A-replay' else 'fixed') / 'script'
    edb = original if arm in ('A-replay', 'A-fixed') else release / 'db' / f'fold{fold}' / ('ipd365' if arm == 'AE-ipd365' else 'panel')
    source = root / 'runs' / donor / 'A'
    bin_dir = root / 'runs' / donor / 'bin'
    if not (source / 'COMPLETE').exists():
        raise ValueError(f'Source A arm incomplete: {source}')
    # Validate the legacy source before adopting it as a frozen input.
    source_outputs = validate_outputs(source, donor)
    required = [source / f'{donor}.realign.sort.bam', source / f'{donor}.realign.sort.bam.bai',
                source / f'{donor}.realign.filter.vcf', source / 'low_depth.bed']
    required += [bin_dir / f'{g}.R{h}.fq.gz' for g in GENES for h in (1, 2)]
    config = dict(donor=donor, fold=fold, arm=arm, threads=threads, source_outputs=source_outputs,
                  inputs={str(p): sha(p) for p in required}, prepared_sha256=sha(release / 'PREPARED.json'),
                  driver_sha256=sha(__file__), phase_database=str(edb), naming_database=str(original))
    out = output or release / 'runs' / donor / arm
    if not start(out, config):
        return
    env = environment(script)
    t0 = time.monotonic()
    try:
        # Inputs are immutable; SpecHLA writes only new files into this run directory.
        for p in required:
            link(p, out / p.name)
        bam = out / f'{donor}.realign.sort.bam'
        vcf = out / f'{donor}.realign.filter.vcf'
        command(['samtools', 'quickcheck', '-v', bam], out / 'input-check.log', env)
        for gene in GENES:
            command(['python3', script / 'phase_variants.py', '-o', out, '-b', bam, '-s', 'nothing', '-v', vcf,
                     '--fq1', out / f'{gene}.R1.fq.gz', '--fq2', out / f'{gene}.R2.fq.gz', '--gene', f'HLA_{gene}',
                     '--freq_bias', '.05', '--snp_qual', '.01', '--snp_dp', '5', '--ref', edb / 'ref' / f'HLA_{gene}.fa',
                     '--tgs', 'NA', '--nanopore', 'NA', '--hic_fwd', 'NA', '--hic_rev', 'NA', '--tenx', 'NA',
                     '--sa', donor, '--weight_imb', '0', '--exon', '0', '--thread_num', threads,
                     '--use_database', '1', '--trio', 'None', '--db', edb], out / f'phase.HLA_{gene}.log', env,
                    timeout=1800)
        command(['perl', script / 'whole/annoHLA.pl', '-s', donor, '-i', out, '-p', 'Unknown',
                 '-d', original / 'HLA', '-r', 'whole'], out / 'annoHLA.log', env)
        command(['python3', script / 'whole/g_group_annotation.py', '-s', donor, '-i', out, '-p', 'Unknown',
                 '-j', threads, '--db', original], out / 'ggroup.log', env)
        outputs = validate_outputs(out, donor)
        (out / 'timing.tsv').write_text(f'stage\tseconds\nphase_and_designate\t{time.monotonic()-t0:.3f}\n')
        finish(out, config, outputs)
        print('COMPLETE', donor, arm, flush=True)
    except BaseException as exc:
        failed(out, exc)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--release', type=Path, required=True)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--folds', default='0')
    run = sub.add_parser('phase')
    run.add_argument('--donor', required=True)
    run.add_argument('--fold', type=int, required=True)
    run.add_argument('--arm', choices=ARMS, required=True)
    run.add_argument('--threads', type=int, default=4)
    run.add_argument('--output', type=Path)
    a = parser.parse_args()
    if a.action == 'prepare':
        prepare(a.root.resolve(), a.release.resolve(), [int(f) for f in a.folds.split(',')])
    else:
        phase_run(a.root.resolve(), a.release.resolve(), a.donor, a.fold, a.arm, a.threads, a.output)


if __name__ == '__main__':
    main()
