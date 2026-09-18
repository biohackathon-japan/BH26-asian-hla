#!/usr/bin/env python3
"""Score a fixed, fully completed paired cohort. Never silently drop failed arms.

Name accuracy is diploid-genotype accuracy. Whole-gene exactness uses global edit
alignment; masked N bases count as errors. Legacy infix matching is also reported
but cannot hide terminal sequence differences in the new primary endpoint.
"""
import argparse
import collections
import csv
import json
from pathlib import Path
import edlib
from experiment import ARMS, completed, validate_outputs
from score import (GENES, HERE, REPO, Nomenclature, compare, pred_pair, read_fasta,
                   read_result, truth_a, truth_slots)


def gene_score(recon, truths):
    matrix = [[edlib.align(t, q, mode='NW', task='distance')['editDistance'] if q else len(t)
               for q in recon] for t in truths]
    assignment = min(((0,1),(1,0)), key=lambda p: sum(matrix[i][p[i]] for i in (0,1)))
    distances = [matrix[i][assignment[i]] for i in (0,1)]
    legacy = min(sum(edlib.align(truths[i], recon[p[i]], mode='HW', task='distance')['editDistance']
                     if recon[p[i]] else len(truths[i]) for i in (0,1)) for p in ((0,1),(1,0)))
    return dict(global_edits=sum(distances), exact_genes=sum(d == 0 for d in distances),
                truth_bases=sum(map(len, truths)), reconstructed_bases=sum(map(len,recon)),
                masked_bases=sum(s.count('N') for s in recon), legacy_infix_edits=legacy)


def verify_run(path, donor):
    manifest = json.loads((path / 'manifest.json').read_text())
    if manifest['status'] != 'complete' or not completed(path, manifest['configuration']):
        raise ValueError(f'Incomplete run: {path}')
    config = manifest['configuration']
    if config['donor'] != donor or config['arm'] != path.name:
        raise ValueError(f'Wrong donor or arm: {path}')
    validate_outputs(path, donor)


def evaluate(runs, donors, arms, out, panel, catalogue):
    # Validate the entire predeclared cohort before producing any aggregate.
    for donor in donors:
        for arm in arms:
            verify_run(runs / donor / arm, donor)
    nom = Nomenclature()
    truth = truth_a(catalogue, set(donors))
    sequences = {g: dict(read_fasta(panel / f'HLA-{g}.fa')) for g in GENES}
    rows = []
    for donor in donors:
        for arm in arms:
            path = runs / donor / arm
            names = read_result(path / 'hla.result.txt')
            gnames = read_result(path / 'hla.result.g.group.txt')
            for gene in GENES:
                recs = truth.get((donor,gene))
                if not recs:
                    raise ValueError(f'Missing assembly truth: {donor}/{gene}')
                ts = [sequences[gene][r['name']][2000:-2000].upper() for r in recs]
                rs = [read_fasta(path / f'hla.allele.{h}.HLA_{gene}.fasta')[0][1].upper() for h in (1,2)]
                row = dict(donor=donor, arm=arm, gene=gene, **gene_score(rs,ts))
                for level, predictions in [('two_field',names),('g_group',gnames),('three_field',names)]:
                    slots = truth_slots(recs,level,nom,gene)
                    row[level+'_eligible'] = int(slots is not None)
                    row[level+'_correct'] = compare(pred_pair(predictions[gene],level,nom,gene),slots)[1] if slots is not None else 0
                rows.append(row)
    out.mkdir(parents=True,exist_ok=True)
    with open(out / 'paired_scores.tsv','w') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    aggregate = collections.defaultdict(collections.Counter)
    for row in rows:
        for g in (row['gene'],'ALL'):
            a=aggregate[(row['arm'],g)]; a['gene_pairs']+=1
            a.update({k:v for k,v in row.items() if isinstance(v,(int,float))})
    summary=[dict(arm=arm,gene=g,**counts) for (arm,g),counts in sorted(aggregate.items())]
    with open(out / 'summary.tsv','w') as stream:
        w=csv.DictWriter(stream,fieldnames=list(summary[0]),delimiter='\t'); w.writeheader(); w.writerows(summary)
    (out / 'cohort.json').write_text(json.dumps(dict(donors=donors,arms=arms,complete=True),indent=2)+'\n')
    for r in summary:
        if r['gene']=='ALL':
            print(f"{r['arm']}: exact genes {r['exact_genes']}/{2*r['gene_pairs']}; global edits {r['global_edits']}; two-field {r['two_field_correct']}/{r['two_field_eligible']}")


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--donors',required=True,help='Comma-separated fixed cohort; all must complete')
    p.add_argument('--arms',default=','.join(ARMS))
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--panel',type=Path,default=REPO/'hla-analysis/source')
    p.add_argument('--catalogue',type=Path,default=REPO/'hla-analysis/results/sequence_catalogue.tsv')
    a=p.parse_args(); donors=a.donors.split(','); arms=a.arms.split(',')
    if len(set(donors))!=len(donors) or len(set(arms))!=len(arms): p.error('Duplicate donors or arms')
    evaluate(a.runs,donors,arms,a.out,a.panel,a.catalogue)
