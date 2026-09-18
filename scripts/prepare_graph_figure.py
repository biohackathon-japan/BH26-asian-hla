#!/usr/bin/env python3
"""Regenerate pgr-tk bundles locally and prepare sequence-free Figure 1 data.

Inputs: final-panel classII.fa, HLA-A.fa and HLA-DPB1.fa; construction MHC.sv.gfa.gz;
haplotype_donors.tsv. Requires pgr-pbundle-decomp, NetworkX, SciPy, NumPy and Graphviz sfdp.
"""
from pathlib import Path
import argparse, collections, csv, gzip, hashlib, json, os, shlex, subprocess
import numpy as np
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def layout(nodes, edges):
    graph=nx.Graph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    if len(nodes)>100:
        pos=nx.kamada_kawai_layout(graph,weight=None)
    else:
        dot='graph G {\ngraph [overlap=prism, start=42];\nnode [shape=point];\n'
        dot+=''.join(f'"{n}";\n' for n in nodes)
        dot+=''.join(f'"{u}" -- "{v}";\n' for u,v in edges)+'}\n'
        result=subprocess.run(['sfdp','-Tplain'],input=dot,text=True,capture_output=True,check=True)
        pos={}
        for line in result.stdout.splitlines():
            v=shlex.split(line)
            if v and v[0]=='node':pos[v[1]]=[float(v[2]),float(v[3])]
    arr=np.array(list(pos.values())); arr-=arr.mean(0)
    _,_,rotation=np.linalg.svd(arr,full_matrices=False);arr=arr@rotation.T
    arr/=np.ptp(arr[:,0])
    return {node:xy.tolist() for node,xy in zip(pos,arr)}

def bundle_graph(name, work, out):
    nodes={};edges=collections.defaultdict(set);paths=collections.defaultdict(list)
    for line in (work/f'{name}.bed').open():
        if line.startswith('#'):continue
        f=line.rstrip().split('\t');bid=f[3].split(':')[0];start,end=int(f[1]),int(f[2])
        node=nodes.setdefault(bid,{'haplotypes':set(),'lengths':[],'ref_intervals':[]})
        node['haplotypes'].add(f[0]);node['lengths'].append(end-start);paths[f[0]].append(bid)
        if f[0].startswith('GRCh38#'):node['ref_intervals'].append([start,end])
    for hap,path in paths.items():
        for a,b in zip(path,path[1:]):
            if a!=b:edges[tuple(sorted((a,b)))].add(hap)
    for node in nodes.values():
        node['carriers']=len(node.pop('haplotypes'))
        node['length']=float(np.median(node.pop('lengths')))
    pos=layout(nodes,edges)
    data={'name':name,'paths':len(paths),
          'nodes':[{'id':k,**v,'xy':pos[k]} for k,v in nodes.items()],
          'edges':[{'source':u,'target':v,'carriers':len(h)} for (u,v),h in edges.items()]}
    (out/f'{name}.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    # Retain decomposition intervals, omitting transient filesystem paths in headers.
    with gzip.open(out/f'{name}.bed.gz','wt') as handle:
        handle.write('# contig\tstart0\tend0\tbundle_descriptor\n')
        for line in (work/f'{name}.bed').open():
            if not line.startswith('#'):handle.write(line)
    print(name,len(paths),'paths;',len(nodes),'bundles;',len(edges),'adjacencies',flush=True)
    return {'paths':len(paths),'bundles':len(nodes),'adjacencies':len(edges)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,required=True)
    p.add_argument('--bundle-dir',type=Path,required=True)
    p.add_argument('--structural-gfa',type=Path,required=True)
    p.add_argument('--donors',type=Path,required=True)
    p.add_argument('--out',type=Path,default=ROOT/'paper/data/graph')
    p.add_argument('--reuse-decompositions',action='store_true',help='Package existing local decompositions with these input files and parameters.')
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True);args.bundle_dir.mkdir(parents=True,exist_ok=True)
    provenance={'pgr_tk':subprocess.check_output(['pgr-pbundle-decomp','--version'],text=True).strip(),
                'graphviz':subprocess.run(['sfdp','-V'],capture_output=True,text=True,check=True).stderr.strip(),
                'networkx':nx.__version__, 'layout':'Kamada-Kawai (unweighted) for >100 nodes; otherwise sfdp, overlap=prism, start=42; principal-axis rotation and isotropic scaling', 'inputs':{},'bundles':{}}
    previous = args.out/'provenance.json'
    if previous.exists():
        prior = json.loads(previous.read_text())
        for key in ['source_locations', 'map']:
            if key in prior: provenance[key] = prior[key]
    for name in ['classII','HLA-A','HLA-DPB1']:
        source=args.input_dir/f'{name}.fa'
        options=['-w','48','-k','56','-r','2','--min-span', '8' if name=='classII' else '12',
                 '--min-cov','0','--min-branch-size','8','--bundle-length-cutoff','500' if name=='classII' else '200',
                 '--bundle-merge-distance','2000' if name=='classII' else '1000']
        cmd=['pgr-pbundle-decomp',*options,str(source),str(args.bundle_dir/name)]
        if not args.reuse_decompositions:
            subprocess.run(cmd,check=True,env={**os.environ,'RAYON_NUM_THREADS':'4'})
        provenance['inputs'][source.name]={'sha256':digest(source),'bytes':source.stat().st_size}
        summary=bundle_graph(name,args.bundle_dir,args.out)
        summary['command']=' '.join(['pgr-pbundle-decomp',*options,name+'.fa',name])
        provenance['bundles'][name]=summary
    nodes={};edges=[]
    with gzip.open(args.structural_gfa,'rt') as f:
        for line in f:
            v=line.rstrip().split('\t')
            if v[0]=='S':
                tags={t.split(':',2)[0]:t.split(':',2)[2] for t in v[3:]}
                nodes[v[1]]={'length':int(tags['LN']),'rank':int(tags['SR']),'source':tags['SN'],'start':int(tags['SO'])}
            elif v[0]=='L':edges.append((v[1],v[3]))
    pos=layout(nodes,edges)
    (args.out/'MHC.structural.json').write_text(json.dumps({'nodes':[{'id':k,**v,'xy':pos[k]} for k,v in nodes.items()],'edges':edges},separators=(',',':'))+'\n')
    provenance['structural_backbone']={'segments':len(nodes),'links':len(edges),'sha256':digest(args.structural_gfa)}
    donors=list(csv.DictReader(args.donors.open(),delimiter='\t'))
    cohorts=[]
    for cohort in sorted({r['cohort'] for r in donors}):
        sub=[r for r in donors if r['cohort']==cohort]
        cohorts.append({'cohort':cohort,'haplotypes':len(sub),'samples':len({r['donor_id'] for r in sub})})
    data={'cohorts':cohorts,'haplotypes':len(donors),'distinct_donors':len({r['donor_id'] for r in donors if r['cohort']!='REF'}),
          'sample_entries':sum(c['samples'] for c in cohorts if c['cohort']!='REF')}
    assert (data['haplotypes'],data['distinct_donors'],data['sample_entries'])==(754,371,376)
    (args.out/'cohorts.json').write_text(json.dumps(data,indent=2)+'\n')
    provenance['inputs']['haplotype_donors.tsv']={'sha256':digest(args.donors)}
    provenance['outputs']={p.name:digest(p) for p in sorted(args.out.iterdir()) if p.suffix in ['.gz','.json','.geojson'] and p.name!='provenance.json'}
    (args.out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Prepared Figure 1 data with provenance.',flush=True)

if __name__=='__main__':main()
