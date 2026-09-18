#!/usr/bin/env python3
"""Plot Figure 1 from prepared graph topology, local pgr-tk bundles and cohort counts."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'paper/data/graph'
OUT=ROOT/'paper/figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
                    'pdf.fonttype':42,'savefig.dpi':260})
TEAL='#007C83';ORANGE='#C77735'

def load(name):return json.loads((DATA/name).read_text())
def title(ax,text):ax.set_title(text,loc='left',fontweight='bold',pad=8)
def network(ax,data,bundle=False,size=4):
    nodes=data['nodes'];pos={n['id']:np.array(n['xy']) for n in nodes}
    pairs=[(e['source'],e['target']) for e in data['edges']] if bundle else data['edges']
    ax.add_collection(LineCollection([[pos[a],pos[b]] for a,b in pairs],colors='#929DA3',linewidths=.45,alpha=.65,zorder=1))
    xy=np.array([n['xy'] for n in nodes])
    if bundle:
        ax.scatter(*xy.T,c=[n['carriers'] for n in nodes],s=size,cmap='viridis',vmin=0,vmax=754,linewidths=0,zorder=2)
    else:
        ax.scatter(*xy.T,c=[TEAL if n['rank']==0 else ORANGE for n in nodes],s=size,linewidths=0,zorder=2)
    ax.autoscale();ax.margins(.07,.2);ax.axis('off')
    return pos

def ref_anchor(nodes,target,bundle=False):
    if bundle:
        choices=[(abs((a+b)/2-target),n) for n in nodes for a,b in n['ref_intervals']]
    else:
        choices=[(abs(n['start']+n['length']/2-target),n) for n in nodes if n['rank']==0 and n['source'].startswith('GRCh38')]
    return min(choices,key=lambda x:x[0])[1]

def annotate(ax,node,label,offset):
    ax.scatter(*node['xy'],s=25,facecolors='none',edgecolors='#111111',linewidths=.7,zorder=4)
    ax.annotate(label,xy=node['xy'],xytext=offset,textcoords='offset points',fontsize=9,
                ha='center',va='center',bbox={'fc':'white','ec':'none','pad':1,'alpha':.9},
                arrowprops={'arrowstyle':'-','color':'#333333','lw':.6},zorder=5)

fig=plt.figure(figsize=(8,8.6))
gs=fig.add_gridspec(3,1,height_ratios=[2.5,1.7,3.5],hspace=.27)
top=gs[0].subgridspec(1,2,width_ratios=[3.2,1],wspace=.12)
mapax=fig.add_subplot(top[0]);info=fig.add_subplot(top[1]);info.axis('off')
title(mapax,'A   Source cohorts')
cohorts=load('cohorts.json');counts={x['cohort']:x['samples'] for x in cohorts['cohorts']}
selected={'Saudi Arabia','United Arab Emirates','China','South Korea','Japan'}
for feature in load('natural_earth_countries.geojson')['features']:
    geom=feature['geometry'];polygons=geom['coordinates'] if geom['type']=='MultiPolygon' else [geom['coordinates']]
    for polygon in polygons:
        xy=np.array(polygon[0])
        if (xy[:,0].max()<25 or xy[:,0].min()>167 or xy[:,1].max()<4 or xy[:,1].min()>63):continue
        chosen=feature['properties']['ADMIN'] in selected
        mapax.add_patch(Polygon(xy,facecolor='#D5E5E5' if chosen else '#F0F1F1',edgecolor='#A8B0B3',linewidth=.35))
locations=[
 ('APR','UAE · APR',(54,24),(70,13)),
 ('JaSaPaGe-Saudi','Saudi Arabia\nJaSaPaGe',(45,24),(39,9)),
 ('CPC-Chinese','China · CPC',(104,35),(93,55)),
 ('KPanRef-Korean','Korea · K-PanRef',(128,36),(132,57)),
 ('JaSaPaGe-Japanese','Japan · JaSaPaGe',(139,36),(150,20))]
for cohort,label,point,xytext in locations:
    mapax.scatter(*point,s=32,c=TEAL,edgecolors='white',linewidths=.6,zorder=5)
    mapax.annotate(f'{label}\n{counts[cohort]} samples',xy=point,xytext=xytext,
        fontsize=8.8,ha='center',va='center',zorder=6,
        arrowprops={'arrowstyle':'-','lw':.6,'color':TEAL},
        bbox={'fc':'white','ec':'none','alpha':.87,'pad':1.8})
mapax.set(xlim=(25,169),ylim=(1,65));mapax.set_aspect(1.15);mapax.axis('off')
info.text(0,1,'754 haplotype entries',fontsize=12,weight='bold',va='top',color=TEAL)
info.text(0,.84,'371 distinct donors',fontsize=10,va='top')
info.text(0,.66,'HPRC release 2\n232 samples, worldwide',fontsize=9.2,va='top')
info.text(0,.44,'Mapped cohorts\n144 sample entries',fontsize=9.2,va='top')
info.text(0,.23,'+ GRCh38 and CHM13',fontsize=9.2,va='top')
info.text(0,-.02,'5 donors occur in both\nHPRC and JaSaPaGe',fontsize=8.5,va='bottom',color='#555555')

ax=fig.add_subplot(gs[1]);data=load('MHC.structural.json')
network(ax,data,size=3.2)
title(ax,'B   MHC structural backbone')
ax.text(1,1.06,f"{len(data['nodes']):,} segments · {len(data['edges']):,} links",transform=ax.transAxes,ha='right',fontsize=9)
anchors=load('reference_anchors.json')
for label,gene,offset in [('HLA-A','HLA-A',(0,18)),('DRB1','HLA-DRB1',(0,-22)),('HLA-DPB1','HLA-DPB1',(0,21))]:
    interval=anchors['genes'][gene]
    target=(interval['start1']+interval['end1'])/2-1
    annotate(ax,ref_anchor(data['nodes'],target),label,offset)
ax.legend(handles=[Line2D([],[],marker='o',ls='',ms=4,color=TEAL,label='GRCh38 backbone'),Line2D([],[],marker='o',ls='',ms=4,color=ORANGE,label='Added sequence')],loc='lower center',bbox_to_anchor=(.5,-.17),ncol=2,frameon=False,fontsize=8.5)

bottom=gs[2].subgridspec(2,2,width_ratios=[1.6,1],hspace=.45,wspace=.25)
c=fig.add_subplot(bottom[:,0]);d=fig.add_subplot(bottom[0,1]);e=fig.add_subplot(bottom[1,1])
for ax,name,heading,size in [(c,'classII','C   Class II / DRB region',8),(d,'HLA-A','D   HLA-A',17),(e,'HLA-DPB1','E   HLA-DPB1',17)]:
    data=load(name+'.json');network(ax,data,bundle=True,size=size)
    title(ax,heading)
    ax.text(0,1.0,f"{data['paths']} haplotypes · {len(data['nodes'])} bundles",transform=ax.transAxes,fontsize=9,va='top')
    if name=='classII':
        # Midpoints in the original GRCh38 MHC extraction, minus class II start.
        anchors=load('reference_anchors.json')
        for label,gene,offset in [('DRA','HLA-DRA',(-5,12)),('DRB1','HLA-DRB1',(-15,-23)),('DQ','HLA-DQA1',(13,19)),('DMA','HLA-DMA',(0,14))]:
            interval=anchors['genes'][gene]
            target=(interval['start1']+interval['end1'])/2-anchors['classII_start1']
            annotate(ax,ref_anchor(data['nodes'],target,bundle=True),label,offset)
cbax=fig.add_axes([.59,.035,.29,.014]);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,754),cmap='viridis'),cax=cbax,orientation='horizontal')
cb.solids.set_rasterized(False)
cb.set_ticks([0,250,500,754]);cb.ax.tick_params(labelsize=8)
cb.set_label('Haplotypes carrying each bundle (C–E)',fontsize=9,labelpad=3)
fig.subplots_adjust(left=.055,right=.97,top=.96,bottom=.1)
for suffix in ['png','pdf','svg']:
    fig.savefig(OUT/f'fig1_asian_graph.{suffix}',bbox_inches='tight',metadata={'Creator':'AsianPGR figure script'})
print('Generated Figure 1: cohort map, structural backbone and three pgr-tk regional views.')

# Presentation derivatives use the same topology and map data, with larger labels.
# Export panels directly from plotting objects; never crop text out of a bitmap.
slide_out=OUT/'slide_panels'
slide_out.mkdir(exist_ok=True)
labels=['UAE\nAPR: 53','Saudi Arabia\nJaSaPaGe: 9','China\nCPC: 58','Korea\nK-PanRef: 14','Japan\nJaSaPaGe: 10']
for text,label in zip(mapax.texts,labels):
    text.set_text(label);text.set_fontsize(18)
mapax.set_title('',loc='left')
for panel in list(fig.axes):
    if panel is not mapax: panel.remove()
fig.set_size_inches(8,2.85)
mapax.set_position([.02,.04,.96,.92])
mapax.set_aspect('auto')
fig.canvas.draw()
mapbox=mapax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
for ext in ['pdf','svg']:
    fig.savefig(slide_out/f'cohort_map.{ext}',bbox_inches=mapbox,pad_inches=.04,facecolor='white')

backbone=load('MHC.structural.json')
f,a=plt.subplots(figsize=(8,1.1));network(a,backbone,size=4)
f.subplots_adjust(left=.01,right=.99,bottom=.03,top=.97)
for ext in ['pdf','svg']:
    f.savefig(slide_out/f'mhc_backbone.{ext}',bbox_inches='tight',pad_inches=.02,transparent=True)
plt.close(f)

f,a=plt.subplots(figsize=(3.5,2.6));regional=load('classII.json');network(a,regional,bundle=True,size=10)
anchors=load('reference_anchors.json');interval=anchors['genes']['HLA-DRB1']
annotate(a,ref_anchor(regional['nodes'],(interval['start1']+interval['end1'])/2-anchors['classII_start1'],bundle=True),'DRB1',(-20,-15))
f.subplots_adjust(left=.02,right=.98,bottom=.02,top=.98)
for ext in ['pdf','svg']:
    f.savefig(slide_out/f'classII_drb.{ext}',bbox_inches='tight',pad_inches=.06,facecolor='white')
plt.close(f)
print('Exported map, backbone and class II/DRB panels for slides.')
