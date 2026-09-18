#!/usr/bin/env python3
"""Bounded graph preservation/mapping/calling gate; run inside a Slurm allocation.

A prefix of read pairs is a performance smoke test, never an accuracy benchmark.
Original GAM and graph-native variant calls are retained alongside surjected BAM.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import re
import statistics
import time

from build_refs import excluded_haps, read_fasta
from experiment import ROOT, ENV, VG, command, environment, failed, finish, json_write, sha, start


def canonical(name):
    # vg may render an unsplit full-length haplotype with a zero-start suffix.
    name = re.sub(r'\[0\]$', '', name.split()[0])
    # GBWT adds phase block #0 to PanSN haplotype paths (but not references).
    # https://github.com/vgteam/vg/wiki/Changing-References
    if len(name.split('#')) == 4 and name.endswith('#0'):
        name = name[:-2]
    return name


def verify_paths(expected, observed, excluded):
    source = dict(expected)
    if len(source) != len(expected):
        raise ValueError('Duplicate input path names')
    actual = {}
    for name, seq in observed:
        name = canonical(name)
        if name in actual:
            raise ValueError(f'Duplicate exported path: {name}')
        actual[name] = seq.upper()
    leaked = [n for n in source if '#'.join(n.split('#')[:2]) in excluded]
    missing = sorted(set(source) - set(actual))
    extra = sorted(set(actual) - set(source))
    changed = sorted(n for n, s in source.items() if n in actual and s.upper() != actual[n])
    return dict(passed=not (leaked or missing or extra or changed), input_paths=len(source),
                exported_paths=len(actual), leaked=leaked, missing=missing, extra=extra, changed=changed)


def subset_pairs(first, second, out1, out2, limit):
    n = 0
    with gzip.open(first, 'rt') as a, gzip.open(second, 'rt') as b, open(out1, 'w') as c, open(out2, 'w') as d:
        while n < limit:
            x, y = [a.readline() for _ in range(4)], [b.readline() for _ in range(4)]
            if not x[0] and not y[0]:
                break
            if not all(x) or not all(y):
                raise ValueError('Truncated or unequal paired FASTQ')
            names = [re.sub(r'/[12]$', '', z[0].split()[0]) for z in (x, y)]
            if names[0] != names[1] or any(not z[0].startswith('@') or not z[2].startswith('+') or len(z[1].strip()) != len(z[3].strip()) for z in (x, y)):
                raise ValueError('Invalid or mismatched paired FASTQ')
            c.writelines(x); d.writelines(y); n += 1
    if not n:
        raise ValueError('No read pairs')
    return n


def bam_metrics(path):
    import pysam
    result = dict(primary=0, mapped=0, mapq20=0, proper_pair=0, softclip_bases=0, query_bases=0)
    with pysam.AlignmentFile(path, 'rb') as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary:
                continue
            result['primary'] += 1
            result['mapped'] += not r.is_unmapped
            result['mapq20'] += not r.is_unmapped and r.mapping_quality >= 20
            result['proper_pair'] += r.is_proper_pair
            result['query_bases'] += r.query_length or 0
            result['softclip_bases'] += sum(n for op, n in r.cigartuples or [] if op == 4)
    return result


def fragment_distribution(path):
    """Estimate insert length from confident A-arm pairs, never from truth paths."""
    import pysam
    with pysam.AlignmentFile(path, 'rb') as bam:
        lengths = [abs(r.template_length) for r in bam if r.is_read1 and r.is_proper_pair
                   and not r.is_secondary and not r.is_supplementary and r.mapping_quality >= 20
                   and 0 < abs(r.template_length) < 2000]
    if len(lengths) < 30:
        raise ValueError('Insufficient confident pairs for fragment-length control')
    median = statistics.median(lengths)
    mad = statistics.median(abs(x-median) for x in lengths)
    trimmed = [x for x in lengths if abs(x-median) <= max(50, 5*mad)]
    return dict(mean=statistics.mean(trimmed), stdev=max(1,statistics.stdev(trimmed)), pairs=len(trimmed))


def graph_support(reference_json, reads_json):
    with open(reference_json) as f:
        reference = [json.loads(line) for line in f]
    nodes = {str(m['position']['node_id']) for a in reference for m in a['path']['mapping']}
    result = dict(reads=0, mapped=0, mapq20=0, aligned_bases=0, off_reference_node_bases=0)
    with open(reads_json) as stream:
        for line in stream:
            a=json.loads(line); mappings=a.get('path',{}).get('mapping',[])
            result['reads']+=1; result['mapped']+=bool(mappings)
            result['mapq20']+=bool(mappings) and a.get('mapping_quality',0)>=20
            for m in mappings:
                length=sum(e.get('to_length',0) for e in m.get('edit',[]) if e.get('from_length',0)>0)
                result['aligned_bases']+=length
                if str(m['position']['node_id']) not in nodes:
                    result['off_reference_node_bases']+=length
    return result


def diagnose(root, out, fold, gene, donor, kind, pairs, threads, timeout):
    with open(root / 'source/donors_folds.tsv') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if not any(r['donor'] == donor and int(r['fold']) == fold for r in rows):
        raise ValueError(f'Incorrect donor/fold assignment: {donor}/{fold}')
    prefix = root / 'graphs' / f'fold{fold}' / f'HLA_{gene}.{kind}'
    paths = {suffix: Path(str(prefix) + suffix) for suffix in ('.giraffe.gbz', '.dist', '.min', '.zipcodes')}
    inp = root / 'graphs' / f'fold{fold}' / f'HLA_{gene}.in.fa'
    reads = [root / 'runs' / donor / 'bin' / f'{gene}.R{i}.fq.gz' for i in (1, 2)]
    exclusions = root / 'source/excluded_paths.tsv'
    groups = root / 'source/validation_groups.tsv'
    baseline = root / 'runs' / donor / 'A' / f'{donor}.realign.sort.bam'
    linear_ref = ENV / 'share/spechla/db/HLA' / f'HLA_{gene}' / f'HLA_{gene}.fa'
    config = dict(fold=fold, gene=gene, donor=donor, kind=kind, pairs=pairs, threads=threads,
                  timeout_seconds=timeout, inputs={str(p): sha(p) for p in [*paths.values(), inp, *reads, exclusions, groups, baseline, linear_ref]},
                  driver_sha256=sha(__file__))
    if not start(out, config):
        return
    env = environment(ENV / 'share/spechla/script')
    env['OMP_NUM_THREADS'] = str(threads)
    log = out / 'commands.log'
    timings = {}
    def run(label, args, output=None):
        t = time.monotonic()
        try:
            if output:
                with open(out / output, 'wb') as f:
                    command(args, log, env, timeout=timeout, stdout=f)
            else:
                command(args, log, env, timeout=timeout)
        finally:
            timings[label] = time.monotonic() - t
            json_write(out / 'timing.json', timings)
    try:
        gbz = paths['.giraffe.gbz']
        run('export_paths', ['vg', 'paths', '-x', gbz, '-F', '-t', threads], 'paths.fa')
        expected, actual = read_fasta(inp), read_fasta(out / 'paths.fa')
        excluded, _ = excluded_haps(root / 'source', fold)
        preservation = verify_paths(expected, actual, excluded)
        json_write(out / 'preservation.json', preservation)
        if not preservation['passed']:
            raise ValueError('Graph path preservation/leakage gate failed; see preservation.json')
        ref = next(n for n, _ in actual if canonical(n) == f'SpecHLA#0#HLA_{gene}')
        fqs = [out / f'reads.R{i}.fq' for i in (1, 2)]
        n = subset_pairs(*reads, *fqs, pairs)
        fragment = fragment_distribution(baseline)
        json_write(out / 'fragment.json', fragment)
        run('giraffe', ['vg', 'giraffe', '-Z', gbz, '-d', paths['.dist'], '-m', paths['.min'],
                        '-z', paths['.zipcodes'], '-f', fqs[0], '-f', fqs[1], '-t', threads,
                        '--fragment-mean', fragment['mean'], '--fragment-stdev', fragment['stdev'],
                        '-N', donor, '-o', 'gam'], 'reads.gam')
        run('alignment_stats', ['vg', 'stats', '-a', out / 'reads.gam'], 'gam.stats.txt')
        run('reference_walk', ['vg','paths','-x',gbz,'-Q','SpecHLA#0#','-X','-t',threads], 'reference.gam')
        run('reference_json', ['vg','view','-aj',out/'reference.gam'], 'reference.jsonl')
        run('reads_json', ['vg','view','-aj',out/'reads.gam'], 'reads.jsonl')
        support = graph_support(out/'reference.jsonl', out/'reads.jsonl')
        run('surject', ['vg', 'surject', '-x', gbz, '-p', ref, '-i', '-b', '-t', threads,
                        '--read-length', 'short', out / 'reads.gam'], 'projected.bam')
        run('bam_check', ['samtools', 'quickcheck', '-v', out / 'projected.bam'])
        metrics = bam_metrics(out / 'projected.bam')
        if metrics['primary'] != n * 2:
            raise ValueError(f'Read loss during projection: {metrics["primary"]} vs {n*2}')
        run('linear_bwa', ['bwa', 'mem', '-t', threads, '-U', '10000', '-L', '10000,10000', linear_ref, *fqs], 'linear.sam')
        run('linear_bam', ['samtools', 'view', '-b', out / 'linear.sam'], 'linear.bam')
        linear_metrics = bam_metrics(out / 'linear.bam')
        run('pack', ['vg', 'pack', '-x', gbz, '-g', out / 'reads.gam', '-o', out / 'reads.pack', '-Q', '5', '-t', threads])
        run('call', ['vg', 'call', gbz, '-k', out / 'reads.pack', '-s', donor, '-z', '-a', '-p', ref, '-t', threads], 'graph.vcf')
        import pysam
        with pysam.VariantFile(str(out / 'graph.vcf')) as vcf:
            records = list(vcf)
            if donor not in vcf.header.samples:
                raise ValueError('Graph VCF has wrong sample')
            called_alt = [r for r in records if any(a is not None and a > 0 for a in (r.samples[donor]['GT'] or ()))]
            structural = []
            for r in called_alt:
                call = r.samples[donor]
                gt = call['GT']
                for allele in sorted({i for i in gt if i is not None and i>0}):
                    alt = r.alleles[allele]
                    if abs(len(alt)-len(r.ref)) < 50: continue
                    structural.append(dict(chrom=r.chrom, pos=r.pos, ref=r.ref, alt=alt,
                                           dosage=gt.count(allele), genotype='/'.join('.' if x is None else str(x) for x in gt),
                                           quality=r.qual, depth=call.get('DP'), genotype_quality=call.get('GQ'),
                                           allele_depth=call.get('AD')[allele] if call.get('AD') else None))
        with open(out/'structural_candidates.tsv','w') as f:
            w=csv.DictWriter(f,fieldnames=['chrom','pos','ref','alt','dosage','genotype','quality','depth','genotype_quality','allele_depth'],delimiter='\t')
            w.writeheader(); w.writerows(structural)
        json_write(out / 'summary.json', dict(read_pairs=n, subset='all binned pairs' if n < pairs else 'first pairs; performance gate only',
                   graph_support=support, projection=metrics, linear=linear_metrics, fragment_distribution=fragment,
                   graph_variant_records=len(records), nonreference_records=len(called_alt),
                   nonreference_indel50_alleles=len(structural), timing_seconds=timings))
        finish(out, config, {name: sha(out / name) for name in ('preservation.json', 'summary.json', 'reads.gam', 'graph.vcf', 'projected.bam','structural_candidates.tsv')})
    except BaseException as exc:
        failed(out, exc)
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--fold', type=int, default=0)
    p.add_argument('--gene', choices=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1'], required=True)
    p.add_argument('--donor', default='HG00706')
    p.add_argument('--kind', choices=['pggb','mc'], default='pggb')
    p.add_argument('--pairs', type=int, default=256)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--timeout', type=int, default=300)
    a=p.parse_args()
    if a.pairs <= 0 or a.threads <= 0 or a.timeout <= 0:
        p.error('pairs, threads and timeout must be positive')
    diagnose(a.root, a.out, a.fold, a.gene, a.donor, a.kind, a.pairs, a.threads, a.timeout)
