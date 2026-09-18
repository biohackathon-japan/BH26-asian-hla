#!/usr/bin/env python3
"""Rebuild manuscript figures from the frozen, repository-local score tables."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
DATA, OUT = ROOT / 'paper/data', ROOT / 'paper/figures'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.titlesize': 10, 'pdf.fonttype': 42,
    'savefig.dpi': 240, 'font.family': 'DejaVu Sans'})
FULL, HPRC = '#007C83', '#666666'

def read(name):
    with (DATA / (name + '.tsv')).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def select(rows, **kwargs):
    return [r for r in rows if all(r[k] == v for k, v in kwargs.items())]

def one(rows, **kwargs):
    result = select(rows, **kwargs)
    assert len(result) == 1, (kwargs, len(result))
    return result[0]

def save(fig, name):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / f'{name}.{ext}', bbox_inches='tight', metadata={'Creator': 'Manuscript figure script'})
    plt.close(fig)

scores = read('named_HLA_scores')
genes = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
fig, axes = plt.subplots(1, 2, figsize=(6.4, 4.2), layout='constrained', sharex=True, sharey=True)
for ax, ancestry, panel in zip(axes, ['EAS', 'SAS'], ['A', 'B']):
    for i, gene in enumerate(genes):
        values = []
        for method, color in [('Locityper_hprc', HPRC), ('Locityper_full', FULL)]:
            rows = select(scores, stratum=ancestry, gene=gene, method=method, endpoint='assembly_exact_CDS_twofield')
            n = len(rows)
            covered = sum(int(r['truth_labels_in_panel']) for r in rows)
            values.append(100 * covered / n)
            ax.scatter(values[-1], i, color=color, marker='o' if method.endswith('full') else 's', s=45, zorder=3)
        ax.plot(values, [i, i], color='#bbbbbb', zorder=1)
        ax.text(102, i, f'{n}', va='center', fontsize=9)
    ax.set_title(f'{panel}   {ancestry} (20 donors)')
    ax.set_xlim(65, 110)
    ax.set_xticks([70, 80, 90, 100])
    ax.set_xlabel('Genotypes with both labels (%)')
    ax.set_yticks(range(8), ['HLA-' + g for g in genes])
    ax.text(102, -0.8, 'n', fontsize=9, ha='left')
    ax.grid(axis='x', color='#eeeeee')
axes[0].invert_yaxis()
fig.legend(handles=[Line2D([], [], marker='s', ls='', color=HPRC, label='HPRC-only'), Line2D([], [], marker='o', ls='', color=FULL, label='Full panel')], loc='outside upper center', ncol=2, frameon=False)
save(fig, 'figS1_reference_coverage')

rows, paired, summary = read('variant_per_donor'), read('variant_paired'), read('variant_summary')
fig, axes = plt.subplots(2, 2, figsize=(6.4, 5.8), layout='constrained')
for ax, ancestry, panel in zip(axes[0], ['EAS', 'SAS'], ['A', 'B']):
    sub = select(rows, stratum=ancestry, universe='all_truth', variant_class='truth_SV_length')
    donors = sorted({r['donor'] for r in sub})
    assert len(donors) == 20
    for donor in donors:
        pair = [one(sub, donor=donor, arm=arm) for arm in ['hprc', 'full']]
        assert pair[0]['n'] == pair[1]['n']
        values = [100 * int(r['correct']) / int(r['n']) for r in pair]
        ax.plot([0, 1], values, color='#b7c7c7', lw=0.8, zorder=1)
        ax.scatter([0, 1], values, color=[HPRC, FULL], s=17, zorder=2)
    for x, arm, color in [(0, 'hprc', HPRC), (1, 'full', FULL)]:
        r = one(summary, stratum=ancestry, universe='all_truth', variant_class='truth_SV_length', arm=arm)
        ax.scatter(x, float(r['accuracy_pct']), marker='D', s=65, color=color, edgecolor='white', zorder=4)
        ax.text(x, 94, f"{r['correct']}/{r['n']}", ha='center', va='top', fontsize=9)
    ax.set(title=f'{panel}   SV recovery: {ancestry}', xticks=[0, 1], xticklabels=['HPRC-only', 'Full panel'], xlim=(-.3, 1.3), ylim=(55, 96), ylabel='Exact genotype recovery (%)')
    ax.grid(axis='y', color='#eeeeee')
for ax, universe, vclass, panel, title, limits in [
    (axes[1,0], 'all_truth', 'truth_SV_length', 'C', 'SV-bearing genotypes', (-1, 7)),
    (axes[1,1], 'frozen_threeway_SNV', 'SNV', 'D', 'Shared SNV genotypes', (-.025, .10))]:
    for y, ancestry in enumerate(['EAS', 'SAS']):
        r = one(paired, stratum=ancestry, universe=universe, variant_class=vclass, endpoint='all', baseline='hprc')
        x, lo, hi = [float(r[k]) for k in ['full_minus_baseline_pp', 'ci95_low_pp', 'ci95_high_pp']]
        ax.errorbar(x, y, xerr=[[x-lo], [hi-x]], fmt='o', color=FULL, capsize=4)
        ax.text(x, y-.20, f'{x:+.2f}' if panel == 'C' else f'{x:+.3f}', ha='center', fontsize=9)
    ax.axvline(0, color='#999999', lw=1, ls='--')
    ax.set(title=f'{panel}   {title}', yticks=[0,1], yticklabels=['EAS','SAS'], ylim=(1.5,-.5), xlim=limits, xlabel='Full − HPRC-only (pp)')
save(fig, 'fig2_variant_recovery')

summary = read('named_HLA_summary')
fig, axes = plt.subplots(3, 1, figsize=(6.4, 6), layout='constrained', sharex=True, sharey=True)
methods = ['Locityper_full', 'Locityper_hprc', 'T1K', 'SpecHLA']
labels = ['Full + Locityper', 'HPRC-only + Locityper', 'T1K', 'SpecHLA']
colors = [FULL, '#D99064', '#D9D9D9']
for ax, ancestry, endpoint, title in zip(axes, ['EAS', 'EAS', 'SAS'], ['experimental_twofield', 'assembly_exact_CDS_twofield', 'assembly_exact_CDS_twofield'], ['A   Experimental EAS: 8 donors, 5 loci', 'B   Assembly EAS: 20 donors, 8 loci', 'C   Assembly SAS: 20 donors, 8 loci']):
    for y, method in enumerate(methods):
        r = one(summary, stratum=ancestry, endpoint=endpoint, method=method, gene='ALL')
        n, called, correct = [int(r[k]) for k in ['genotypes', 'called', 'correct']]
        counts = [correct, called-correct, n-called]
        assert sum(counts) == n and min(counts) >= 0
        left = 0
        for count, color in zip(counts, colors):
            width = 100*count/n
            ax.barh(y, width, left=left, color=color, height=.6)
            left += width
        ax.text(2, y, f'{correct}/{n}', color='white', va='center', fontsize=10, weight='bold')
    ax.set(title=title, xlim=(0,100), xticks=[0,25,50,75,100], xlabel='', yticks=range(4), yticklabels=labels)
axes[0].invert_yaxis()
axes[-1].set_xlabel('Eligible genotypes (%)')
fig.legend(handles=[Patch(color=c, label=l) for c,l in zip(colors, ['Correct pair', 'Incorrect pair', 'Unresolved'])], loc='outside lower center', ncol=3, frameon=False)
save(fig, 'fig3_hla_typing')
print('Generated Figures 2–3 and supplementary Figure S1 as PNG, vector PDF and SVG from frozen tables.')
