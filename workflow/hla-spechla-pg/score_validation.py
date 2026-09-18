#!/usr/bin/env python3
"""Intention-to-test DogoHLA validation report, including missing/failed donors."""
import argparse
import collections
import csv
import json
from pathlib import Path
import random
from score_experiment import gene_score, verify_run
from score import GENES, REPO, Nomenclature, compare, pred_pair, read_fasta, read_result, truth_a, truth_slots
from truth import gourraud

ARMS = ('native', 'native-long', 'DogoHLA-no-graph', 'DogoHLA')


def locate(root, donor, arm):
    return root / ('runs' if arm.startswith('native') else 'final/runs') / donor / arm


def status(path, donor):
    if not path.exists():
        return 'not_started', 'No output directory'
    if not (path / 'COMPLETE').exists():
        try:
            m = json.loads((path / 'manifest.json').read_text())
            return m.get('status', 'incomplete'), m.get('error', '')
        except (OSError, ValueError):
            return 'incomplete', 'No validated completion'
    try:
        verify_run(path, donor)
        return 'complete', ''
    except Exception as exc:
        return 'invalid', repr(exc)


def write_tsv(path, rows):
    if not rows:
        return
    with open(path, 'w') as f:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def bootstrap(values, seed=20260918):
    if not values:
        return None
    rng = random.Random(seed)
    n = len(values)
    means = sorted(sum(rng.choices(values, k=n))/n for _ in range(10000))
    return dict(mean=sum(values)/n, lower_95=means[249], upper_95=means[9749], donor_count=n)


def evaluate(root, cohort, out):
    donors = list(csv.DictReader(open(cohort), delimiter='\t'))
    if len({d['donor'] for d in donors}) != len(donors) or len({d['family'] for d in donors}) != len(donors):
        raise ValueError('This bootstrap requires unique donors and families')
    nom = Nomenclature()
    truth = truth_a(REPO/'hla-analysis/results/sequence_catalogue.tsv', {d['donor'] for d in donors})
    sequences = {g: dict(read_fasta(REPO/'hla-analysis/source'/f'HLA-{g}.fa')) for g in GENES}
    experimental, _ = gourraud('two_field', nom)
    rows, statuses = [], []
    for d in donors:
        donor = d['donor']
        for arm in ARMS:
            path = locate(root, donor, arm)
            state, reason = status(path, donor)
            statuses.append(dict(donor=donor, arm=arm, status=state, reason=reason))
            valid = state == 'complete'
            names = read_result(path/'hla.result.txt') if valid else {}
            groups = read_result(path/'hla.result.g.group.txt') if valid else {}
            for gene in GENES:
                records = truth.get((donor, gene))
                eligible = records is not None
                r = dict(donor=donor, family=d['family'], stratum=d['stratum'], population=d['population'],
                         arm=arm, gene=gene, status=state, assembly_eligible=int(eligible),
                         completed=int(valid), exact_genes=0, global_edits=None, masked_bases=None,
                         truth_bases=None, reconstructed_bases=None, legacy_infix_edits=None)
                if eligible and valid:
                    ts = [sequences[gene][v['name']][2000:-2000].upper() for v in records]
                    rs = [read_fasta(path/f'hla.allele.{h}.HLA_{gene}.fasta')[0][1].upper() for h in (1,2)]
                    r.update(gene_score(rs, ts))
                for level, prediction in [('two_field', names), ('g_group', groups), ('three_field', names)]:
                    slots = truth_slots(records, level, nom, gene) if eligible else None
                    r[level+'_eligible'] = int(slots is not None)
                    r[level+'_correct'] = int(compare(pred_pair(prediction.get(gene,[None,None]), level, nom, gene), slots)[1]) if slots is not None else 0
                slots = experimental.get((donor, gene))
                r['experimental_eligible'] = int(slots is not None)
                r['experimental_correct'] = int(compare(pred_pair(names.get(gene,[None,None]), 'two_field', nom, gene), slots)[1]) if slots is not None else 0
                r['two_field_called'] = int(valid and all(names.get(gene, [None,None])))
                rows.append(r)
    aggregate = collections.defaultdict(collections.Counter)
    for row in rows:
        for stratum, gene in [('ALL','ALL'), (row['stratum'],'ALL'), ('ALL',row['gene']), (row['stratum'],row['gene'])]:
            a = aggregate[(row['arm'], stratum, gene)]
            a['planned_gene_pairs'] += 1
            for key, value in row.items():
                if isinstance(value, (int,float)):
                    a[key] += value
            a['edit_evaluable_pairs'] += int(row['global_edits'] is not None)
    summaries = [dict(arm=a, stratum=s, gene=g, **v) for (a,s,g),v in sorted(aggregate.items())]
    paired = []
    for baseline in ('native', 'native-long', 'DogoHLA-no-graph'):
        for d in donors:
            donor = d['donor']
            arms = {a: [r for r in rows if r['donor']==donor and r['arm']==a] for a in (baseline,'DogoHLA')}
            valid = all(r['status']=='complete' and r['assembly_eligible'] for rs in arms.values() for r in rs)
            paired.append(dict(donor=donor, stratum=d['stratum'], baseline=baseline,
                both_complete=int(valid), exact_gain=sum(r['exact_genes'] for r in arms['DogoHLA'])-sum(r['exact_genes'] for r in arms[baseline]),
                edit_reduction=(sum(r['global_edits'] for r in arms[baseline])-sum(r['global_edits'] for r in arms['DogoHLA'])) if valid else None,
                two_field_gain=sum(r['two_field_correct'] for r in arms['DogoHLA'])-sum(r['two_field_correct'] for r in arms[baseline]),
                experimental_gain=sum(r['experimental_correct'] for r in arms['DogoHLA'])-sum(r['experimental_correct'] for r in arms[baseline])))
    intervals = {}
    for baseline in ('native','native-long','DogoHLA-no-graph'):
        for stratum in ('ALL','EAS','SAS'):
            rs=[r for r in paired if r['baseline']==baseline and (stratum=='ALL' or r['stratum']==stratum)]
            intervals[baseline+'/'+stratum]={k:bootstrap([r[k] for r in rs if r[k] is not None])
                for k in ('exact_gain','edit_reduction','two_field_gain','experimental_gain')}
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out/'status.tsv',statuses); write_tsv(out/'gene_scores.tsv',rows)
    write_tsv(out/'summary.tsv',summaries); write_tsv(out/'paired_donors.tsv',paired)
    (out/'bootstrap.json').write_text(json.dumps(intervals,indent=2)+'\n')
    complete = all(s['status']=='complete' for s in statuses)
    metadata=dict(planned_donors=len(donors), arms=ARMS, complete=complete,
                  statuses=dict(collections.Counter(s['status'] for s in statuses)), seed=20260918,
                  bootstrap_replicates=10000, note='Exact/name gains use all planned donors; edit reductions use paired complete donors. Intervals are exploratory, not multiplicity-adjusted.')
    (out/'analysis.json').write_text(json.dumps(metadata,indent=2)+'\n')
    text=['# DōgoHLA 0.1.0 validation', '', '**'+('Complete' if complete else 'INCOMPLETE — do not interpret interim differences as final validation')+'**', '',
          '| Method | Completed donors | Exact haplotypes / eligible | Two-field genotypes | Experimental genotypes |',
          '|---|---:|---:|---:|---:|']
    for arm in ARMS:
        r=next(v for v in summaries if v['arm']==arm and v['stratum']=='ALL' and v['gene']=='ALL')
        text.append(f"| {arm} | {sum(s['arm']==arm and s['status']=='complete' for s in statuses)}/{len(donors)} | {r['exact_genes']}/{2*r['assembly_eligible']} | {r['two_field_correct']}/{r['two_field_eligible']} | {r['experimental_correct']}/{r['experimental_eligible']} |")
    text += ['', 'All planned donors remain in the exact/name denominators. Missing outputs receive no correct credit. Global edits require paired completed runs; see `paired_donors.tsv`. Donor/family bootstrap intervals and ancestry strata are retained in the accompanying tables. These are additional development-cohort donors, not the locked validation cohort.', '']
    (out/'REPORT.md').write_text('\n'.join(text))
    print('\n'.join(text))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--cohort',type=Path,default=Path(__file__).parent/'validation/20260918/cohort.tsv')
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); evaluate(a.root,a.cohort,a.out)
