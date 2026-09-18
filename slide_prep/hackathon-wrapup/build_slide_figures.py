#!/usr/bin/env python3
"""Export vector slide panels from the frozen manuscript graph and paired scores."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
p=argparse.ArgumentParser();p.add_argument('--graph-data',type=Path,required=True);a=p.parse_args()
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'figures/paper-panels'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
def load(n):return json.loads((a.graph_data/n).read_text())
def save(fig,name):
 for ext in ['pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',bbox_inches='tight',pad_inches=.06)
 plt.close(fig)
def network(ax,name,bundle=False):
 d=load(name);pos={n['id']:n['xy'] for n in d['nodes']}
 edges=[(e['source'],e['target']) for e in d['edges']] if bundle else d['edges']
 ax.add_collection(LineCollection([[pos[u],pos[v]] for u,v in edges],color='#87939b',linewidths=.6,alpha=.8))
 xy=np.array([n['xy'] for n in d['nodes']]);colors=[n['carriers'] for n in d['nodes']] if bundle else ['#32789d' if n['rank']==0 else '#c77838' for n in d['nodes']]
 opts={'cmap':'viridis','vmin':0,'vmax':754} if bundle else {}
 sc=ax.scatter(*xy.T,c=colors,s=9 if bundle else 8,linewidths=0,**opts)
 ax.autoscale();ax.margins(.02,.10);ax.axis('off');return d,sc
f,ax=plt.subplots(figsize=(7.2,2.65));d,_=network(ax,'MHC.structural.json')
anchors=load('reference_anchors.json')
for gene,label,offset in [('HLA-A','HLA-A',(0,18)),('HLA-DRB1','DRB1',(0,-24)),('HLA-DPB1','DPB1',(0,20))]:
 iv=anchors['genes'][gene];target=(iv['start1']+iv['end1'])/2-1
 node=min([n for n in d['nodes'] if n['rank']==0 and n['source'].startswith('GRCh38')],key=lambda n:abs(n['start']+n['length']/2-target))
 ax.annotate(label,node['xy'],xytext=offset,textcoords='offset points',ha='center',fontsize=11,arrowprops={'arrowstyle':'-','lw':.6},bbox={'fc':'white','ec':'none','pad':1})
f.subplots_adjust(left=.01,right=.99,top=.92,bottom=.08);save(f,'mhc_backbone_large')
f,ax=plt.subplots(figsize=(6.5,2.8));d,sc=network(ax,'classII.json',True)
iv=anchors['genes']['HLA-DRB1'];target=(iv['start1']+iv['end1'])/2-anchors['classII_start1']
node=min([(abs((s+e)/2-target),n) for n in d['nodes'] for s,e in n['ref_intervals']],key=lambda v:v[0])[1]
ax.annotate('DRB1',node['xy'],xytext=(-20,-17),textcoords='offset points',fontsize=11,arrowprops={'arrowstyle':'-','lw':.6},bbox={'fc':'white','ec':'none','pad':1})
f.subplots_adjust(left=.01,right=.99,top=.98,bottom=.05);save(f,'classII_drb_large')
f,ax=plt.subplots(figsize=(7,2.7))
for feature in load('natural_earth_countries.geojson')['features']:
 g=feature['geometry'];polys=g['coordinates'] if g['type']=='MultiPolygon' else [g['coordinates']]
 for polygon in polys:
  xy=np.array(polygon[0]);chosen=feature['properties']['ADMIN'] in ['Saudi Arabia','United Arab Emirates','China','South Korea','Japan']
  ax.add_patch(Polygon(xy,fc='#dce7ee' if chosen else '#f3f3f3',ec='#b0b6ba',lw=.4))
for label,point,xy in [('UAE',(54,24),(68,12)),('Saudi Arabia',(45,24),(35,6)),('China',(104,35),(98,54)),('Korea',(128,36),(130,58)),('Japan\nHPRC + JaSaPaGe',(139,36),(153,16))]:
 ax.plot(*point,'o',ms=4,color='#32789d');ax.annotate(label,point,xytext=xy,fontsize=11,ha='center',arrowprops={'arrowstyle':'-','lw':.6},bbox={'fc':'white','ec':'none','pad':1})
ax.set(xlim=(25,177),ylim=(0,66));ax.axis('off');f.subplots_adjust(left=0,right=1,top=1,bottom=0);save(f,'cohort_map_hprc')
rows=list(csv.DictReader((ROOT/'data/variant_paired.tsv').open(),delimiter='\t'))
f,axs=plt.subplots(1,2,figsize=(12,4.0));blue='#245f85'
for ax,heading,universe,kind,limits,ticks in [(axs[0],'Structural variants','all_truth','truth_SV_length',(-1,6),[0,2,4,6]),(axs[1],'SNVs','frozen_threeway_SNV','SNV',(-.025,.10),[0,.025,.05,.075,.10])]:
 ax.set_title(heading,loc='left',fontsize=16,pad=20)
 ax.axvline(0,color='#7a7a7a',lw=1)
 for y,stratum in [(2,'EAS'),(1,'SAS')]:
  r=next(r for r in rows if r['stratum']==stratum and r['universe']==universe and r['variant_class']==kind and r['endpoint']=='all' and r['baseline']=='hprc')
  v,lo,hi=[float(r[k]) for k in ['full_minus_baseline_pp','ci95_low_pp','ci95_high_pp']]
  ax.errorbar(v,y,xerr=[[v-lo],[hi-v]],fmt='o',color=blue,ms=8,lw=2,capsize=4)
  ax.annotate(f'+{v:.2f}' if kind=='SNV' else f'+{v:.1f}',(v,y),xytext=(0,14),textcoords='offset points',ha='center',fontsize=12)
 ax.text(.5,0,'Not yet evaluated',transform=ax.get_yaxis_transform(),ha='center',va='center',color='#737373',fontsize=12)
 ax.set(ylim=(-.5,2.7),xlim=limits,xticks=ticks,yticks=[2,1,0],yticklabels=['East Asian (20)','South Asian (20)','Non-Asian'])
 ax.tick_params(axis='y',length=0,pad=12);ax.tick_params(axis='x',labelsize=10)
 ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color('#aaaaaa')
 ax.set_xlabel('Improvement over HPRC (percentage points)',labelpad=12,fontsize=11)
f.subplots_adjust(left=.16,right=.98,bottom=.20,top=.84,wspace=.62);save(f,'variant_improvement')
