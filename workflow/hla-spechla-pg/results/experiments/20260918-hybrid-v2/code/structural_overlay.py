#!/usr/bin/env python3
"""Conservative homozygous graph-allele reconstruction on SpecHLA gene bodies.

This prototype accepts PASS homozygous alternate indels >=50 bp with GQ>=20,
DP>=10 and alternate depth>=5. It preserves the source reconstruction outside
nonoverlapping replacement intervals. Heterozygous/uncertain graph alleles are
reported but never assigned to phase arbitrarily. No truth is used by inference.
"""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil

import edlib
from build_refs import GENES, read_fasta, write_fasta
from experiment import ENV, ROOT, command, environment, failed, finish, json_write, sha, start, validate_outputs, completed

# Upstream SpecHLA's 0-based half-open reconstructed gene-body coordinates.
BOUNDS = dict(A=(1000,4503),B=(1000,5081),C=(1000,5304),DPA1=(1000,10775),
              DPB1=(1000,12468),DQA1=(1000,7492),DQB1=(1000,8480),DRB1=(1000,12229))


def project_interval(sequence, reference, start, end, anchor=30):
    """Find boundaries in a reconstructed haplotype, with contiguous flank checks.

CIGAR query=sequence, target=reference: I consumes query, D consumes target.
Require uniquely located, near-identical flanks and no boundary-crossing indel.
"""
    if not (anchor <= start < end <= len(reference)-anchor):
        raise ValueError('Replacement lacks flanking sequence')
    cigar=edlib.align(sequence,reference,mode='NW',task='path')['cigar']
    mapping={}; q=t=0
    for n, op in re.findall(r'(\d+)([=XID])',cigar):
        n=int(n)
        if op in '=X':
            for i in range(n): mapping[t+i]=q+i
            q+=n; t+=n
        elif op=='I': q+=n
        else: t+=n
    if start not in mapping or end not in mapping:
        raise ValueError('Replacement boundary is deleted')
    left,right=mapping[start],mapping[end]
    if right <= left:
        raise ValueError('Non-monotonic replacement coordinates')
    for a,b in ((start-anchor,start),(end,end+anchor)):
        indices=[mapping.get(i) for i in range(a,b)]
        if any(i is None for i in indices) or indices!=list(range(indices[0],indices[0]+anchor)):
            raise ValueError('Indel at replacement flank')
        fragment=sequence[indices[0]:indices[-1]+1]
        ref=reference[a:b]
        if 'N' in fragment or sum(x!=y for x,y in zip(fragment,ref))>3:
            raise ValueError('Uncertain replacement flank')
        if sequence.count(fragment)!=1 or reference.count(ref)!=1:
            raise ValueError('Non-unique replacement flank')
    return left,right


def apply_alleles(sequence,reference,variants):
    intervals=[]
    for start,ref,alt in sorted(variants):
        end=start+len(ref)
        if reference[start:end]!=ref:
            raise ValueError('Graph VCF REF does not match SpecHLA reference')
        if intervals and start<intervals[-1][0]:
            raise ValueError('Overlapping structural calls')
        left,right=project_interval(sequence,reference,start,end)
        intervals.append((end,left,right,alt))
    for _,left,right,alt in reversed(intervals):
        sequence=sequence[:left]+alt+sequence[right:]
    return sequence


def candidates(path,donor,gene):
    import pysam
    accepted=[]; audit=[]
    offset,stop=BOUNDS[gene]
    with pysam.VariantFile(str(path)) as vcf:
        if donor not in vcf.header.samples: raise ValueError('Wrong graph VCF donor')
        for r in vcf:
            call=r.samples[donor]; gt=call.get('GT') or ()
            called={i for i in gt if i is not None and i>0}
            if not any(abs(len(r.alleles[i])-len(r.ref))>=50 for i in called): continue
            reason='accepted'
            if len(gt)!=2 or gt[0]!=gt[1] or gt[0] in (None,0): reason='not_confident_homozygous'
            elif set(r.filter)!={'PASS'} or (call.get('GQ') or 0)<20 or (call.get('DP') or 0)<10 or (r.qual or 0)<20: reason='low_confidence'
            elif not call.get('AD') or call['AD'][gt[0]]<5: reason='low_allele_support'
            elif r.chrom!=f'HLA_{gene}' or r.start<offset or r.stop>stop: reason='outside_gene_body'
            elif set(r.ref+''.join(r.alts or []))-set('ACGT'): reason='non_explicit_allele'
            row=dict(gene=gene,pos=r.pos,genotype='/'.join('.' if i is None else str(i) for i in gt),
                     gq=call.get('GQ'),depth=call.get('DP'),ref_length=len(r.ref),
                     alt_lengths=','.join(str(len(r.alleles[i])) for i in sorted(called)),decision=reason)
            audit.append(row)
            if reason=='accepted': accepted.append((r.start-offset,r.ref,r.alleles[gt[0]]))
    return accepted,audit


def run(source,graph,donor,out,genes,threads):
    manifest=json.loads((source/'manifest.json').read_text())
    if manifest['status']!='complete' or not completed(source,manifest['configuration']):
        raise ValueError('Source phasing run is incomplete or changed')
    if manifest['configuration']['donor']!=donor: raise ValueError('Wrong source donor')
    inputs=validate_outputs(source,donor)
    for gene in genes:
        d=graph/gene
        m=json.loads((d/'manifest.json').read_text())
        if m['configuration']['donor']!=donor or m['configuration']['gene']!=gene:
            raise ValueError('Wrong graph diagnosis donor/gene')
        if m['configuration']['fold']!=manifest['configuration']['fold']:
            raise ValueError('Graph and phasing folds differ')
        source_inputs={Path(p).name:h for p,h in manifest['configuration']['inputs'].items()}
        graph_inputs={Path(p).name:h for p,h in m['configuration']['inputs'].items()}
        for mate in (1,2):
            name=f'{gene}.R{mate}.fq.gz'
            if source_inputs[name]!=graph_inputs[name]:
                raise ValueError('Graph and phasing read inputs differ')
        if m['status']!='complete' or not completed(d,m['configuration']):
            raise ValueError('Graph run incomplete or changed')
        s=json.loads((d/'summary.json').read_text())
        if s['subset']!='all binned pairs': raise ValueError('Structural inference requires all binned pairs')
    refpath=ENV/'share/spechla/db/ref/hla.ref.extend.fa'
    config=dict(donor=donor,arm=out.name,source=str(source),source_outputs=inputs,
                source_manifest_sha256=sha(source/'manifest.json'),reference_sha256=sha(refpath),
                graph_inputs={g:sha(graph/g/'graph.vcf') for g in genes},driver_sha256=sha(__file__),
                genes=genes,threads=threads,min_gq=20,min_dp=10,min_ad=5,min_qual=20,anchor=30,max_anchor_mismatches=3)
    if not start(out,config): return
    env=environment(ENV/'share/spechla/script')
    try:
        refs=dict(read_fasta(refpath)); audit=[]
        for path in source.glob('hla.allele.*.fasta'): shutil.copy2(path,out/path.name)
        for gene in genes:
            variants,rows=candidates(graph/gene/'graph.vcf',donor,gene)
            reference=refs[f'HLA_{gene}'][slice(*BOUNDS[gene])].upper()
            if variants:
                rewritten=[]
                try:
                    for hap in (1,2):
                        path=out/f'hla.allele.{hap}.HLA_{gene}.fasta'
                        name,sequence=read_fasta(path)[0]
                        rewritten.append((path,name,apply_alleles(sequence.upper(),reference,variants)))
                except ValueError as exc:
                    for row in rows:
                        if row['decision']=='accepted': row['decision']='withheld: '+str(exc)
                else:
                    for path,name,seq in rewritten: write_fasta(path,[(name,seq)])
            audit+=rows
        json_write(out/'structural_audit.json',audit)
        db=ENV/'share/spechla/db'; script=ENV/'share/spechla/script'
        command(['perl',script/'whole/annoHLA.pl','-s',donor,'-i',out,'-p','Unknown','-d',db/'HLA','-r','whole'],out/'annoHLA.log',env,timeout=600)
        command(['python3',script/'whole/g_group_annotation.py','-s',donor,'-i',out,'-p','Unknown','-j',threads,'--db',db],out/'ggroup.log',env,timeout=600)
        outputs=validate_outputs(out,donor); outputs['structural_audit.json']=sha(out/'structural_audit.json')
        finish(out,config,outputs)
        print(json.dumps(audit),flush=True)
    except BaseException as exc:
        failed(out,exc); raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--graph',type=Path,required=True,help='Graph diagnostic directory for this donor')
    p.add_argument('--donor',required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--genes',default='DRB1')
    p.add_argument('--threads',type=int,default=4)
    a=p.parse_args(); genes=a.genes.split(',')
    if set(genes)-set(GENES): p.error('Unknown gene')
    run(a.source,a.graph,a.donor,a.out,genes,a.threads)
