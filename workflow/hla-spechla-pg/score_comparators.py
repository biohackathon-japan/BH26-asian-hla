#!/usr/bin/env python3
"""Supplemental matched T1K comparison and strict four-numeric-field naming accuracy.

Four-field truth uses exact genomic matches, never CDS matches. Underspecified
predictions are no-calls at four fields; alternative predictions must all match.
The original frozen primary validation scoring is not changed.
"""
import argparse
import csv
import json
import re
from pathlib import Path
from experiment import sha
from score import GENES, REPO, Nomenclature, compare, pred_pair, read_result, split_name, truth_a, truth_slots
from score_validation import ARMS, locate, status, write_tsv

FIRST_SIX = ('HG00658','HG02132','HG02178','HG03831','HG03874','HG04199')
T1K = REPO/'hla-typer/results/comparators/t1k'


def four_field(name):
    parsed=split_name(name)
    return parsed[0]+'*'+':'.join(parsed[1]) if parsed and len(parsed[1])==4 else None


def genomic_truth(records):
    slots=[]
    for r in records:
        names=[n for n in r['exact_genomic_alleles'].split(';') if n]
        # Strict endpoint: both haplotypes must have genomic identity established
        # at four fields. Do not pad shorter names or use a CDS-compatible label.
        if not names or any(four_field(n) is None for n in names): return None
        slots.append(frozenset(four_field(n) for n in names))
    return tuple(slots)


def four_pair(pair):
    result=[]
    for name in pair:
        alternatives=[four_field(n.strip()) for n in re.split('[,;]',name or '')]
        result.append([frozenset([n]) for n in alternatives] if alternatives and all(alternatives) else None)
    return result


def t1k_calls(path):
    calls={}
    for line in path.read_text().splitlines():
        v=line.split('\t'); gene=v[0].removeprefix('HLA-')
        if gene not in GENES: continue
        n=int(v[1]); a=v[2].removeprefix('HLA-') if n>=1 and float(v[4])>0 else None
        b=a if n==1 else v[5].removeprefix('HLA-') if n==2 and float(v[7])>0 else None
        calls[gene]=[a,b]
    return calls


def evaluate(root,cohort,out):
    donors=list(csv.DictReader(open(cohort),delimiter='\t'))
    donor_ids={d['donor'] for d in donors}
    truth=truth_a(REPO/'hla-analysis/results/sequence_catalogue.tsv',donor_ids)
    nom=Nomenclature(); rows=[]; states={}; input_hashes={}
    methods=(*ARMS,'T1K')
    for d in donors:
        donor=d['donor']
        for method in methods:
            if method=='T1K':
                p=T1K/donor/f'{donor}_genotype.tsv'
                valid=p.exists(); names=t1k_calls(p) if valid else {}
                state='archived_output' if valid else 'missing'
                if valid: input_hashes[str(p.relative_to(REPO))]=sha(p)
            else:
                p=locate(root,donor,method);state,_=status(p,donor)
                valid=state=='complete'; names=read_result(p/'hla.result.txt') if valid else {}
            states[donor,method]=valid
            for gene in GENES:
                records=truth.get((donor,gene))
                for level in ('two_field','four_field'):
                    slots=(genomic_truth(records) if level=='four_field' else truth_slots(records,level,nom,gene)) if records else None
                    pair=names.get(gene,[None,None])
                    predictions=four_pair(pair) if level=='four_field' else pred_pair(pair,level,nom,gene)
                    called,correct,matches=compare(predictions,slots) if slots is not None else (0,0,0)
                    rows.append(dict(donor=donor,stratum=d['stratum'],method=method,gene=gene,level=level,status=state,
                        eligible=int(slots is not None),called=called,correct=correct,allele_matches=matches,
                        allele1=pair[0] or 'NO_CALL',allele2=pair[1] or 'NO_CALL'))
    complete=[d['donor'] for d in donors if all(states[d['donor'],m] for m in ('native','DogoHLA-no-graph','DogoHLA','T1K'))]
    cohorts={'first_six':list(FIRST_SIX),'completed_matched':complete,'planned_32':sorted(donor_ids)}
    summary=[]
    for scope,ds in cohorts.items():
        for method in methods:
            for level in ('two_field','four_field'):
                for gene in ('ALL',*GENES):
                    rs=[r for r in rows if r['donor'] in ds and r['method']==method and r['level']==level and (gene=='ALL' or gene==r['gene'])]
                    summary.append(dict(scope=scope,donors=len(ds),completed_donors=sum(states[d,method] for d in ds),method=method,level=level,gene=gene,
                        **{k:sum(r[k] for r in rs) for k in ('eligible','called','correct','allele_matches')}))
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/'comparator_gene_scores.tsv',rows);write_tsv(out/'comparator_summary.tsv',summary)
    meta=dict(cohorts=cohorts,t1k_output_sha256=input_hashes,truth_source='exact_genomic_alleles for four fields; exact_cds_alleles for two fields',
        endpoint='Strict four numeric fields, expression suffix ignored; both truth haplotypes require exact genomic matches with four-field names. Shorter predictions are unresolved, never expanded into correct calls.',
        provenance='Archived T1K 1.0.6, hla-wgs, IPD 3.65, shared hla-typer/reads recruitment; historical outputs lack per-run FASTQ hashes. Native/Dogo naming uses IPD 3.38. No claim of equal database versions.',
        frozen_primary_unchanged=True)
    (out/'comparator_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    report=['# Matched T1K comparison and four-field performance','',
        'Supplemental analysis requested after the first six outcomes were available. Primary frozen scoring is unchanged. Figures are diploid genotype accuracy, not individual alleles.','']
    for scope in ('first_six','completed_matched'):
        ds=cohorts[scope]
        report += [f'## {scope}: {len(ds)} donors','',
            '| Method | Completed | Two-field correct | Four-field correct | Four-field resolved |',
            '|---|---:|---:|---:|---:|']
        for method in methods:
            a=next(r for r in summary if r['scope']==scope and r['method']==method and r['level']=='two_field' and r['gene']=='ALL')
            b=next(r for r in summary if r['scope']==scope and r['method']==method and r['level']=='four_field' and r['gene']=='ALL')
            four=(f"{b['correct']}/{b['eligible']}" if b['called'] else 'Not resolved') if a['completed_donors'] else 'Pending'
            report.append(f"| {method} | {a['completed_donors']}/{len(ds)} | {a['correct']}/{a['eligible']} | {four} | {b['called']}/{b['eligible']} |")
        report.append('')
    report += ['Four-field truth requires an exact genomic match for both haplotypes. Unmatched/novel truth sequences are excluded for every method and are counted in the separate whole-gene endpoint. Shorter allele names are not silently expanded. T1K\'s archived output may stop at three fields; that is an output-resolution limit, not proof its reconstructed sequence is wrong.','',
        'T1K uses IPD 3.65 whereas native/DōgoHLA allele naming uses IPD 3.38. These are comparisons of the deployed configurations, not an isolated algorithm comparison. No matched OptiType, HLA*LA or other established short-read caller outputs were found in the local project.','']
    (out/'COMPARATORS.md').write_text('\n'.join(report))
    print('\n'.join(report))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--cohort',type=Path,default=Path(__file__).parent/'validation/20260918/cohort.tsv')
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();evaluate(a.root,a.cohort,a.out)
