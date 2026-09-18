---
title: 'Evaluating an Asian-enriched MHC pangenome for variant genotyping and HLA typing'
title_short: 'BH26JP: Asian-enriched MHC pangenome'
tags:
  - Pangenomics
  - HLA
  - Population genetics
  - Genome assembly
  - Genotyping
authors: []
affiliations: []
date: 18 September 2026
cito-bibliography: paper.bib
event: BH26JP
biohackathon_name: "DBCLS BioHackathon 2026"
biohackathon_url:   "https://2026.biohackathon.org/"
biohackathon_location: "Matsuyama, Japan, 2026"
group: asian-hla
# URL to project git repo --- should contain the actual paper.md:
git_url: https://github.com/biohackathon-japan/BH26-asian-hla
# This is the short authors description that is used at the
# bottom of the generated paper (typically the first two authors):
authors_short: ""
---

# Abstract

Does expanding an HPRC reference panel with Asian and Arab haplotypes improve short-read
MHC genotyping? We constructed AsianPGR from 754 MHC haplotype entries and compared full
and HPRC-only inference panels in 20 East Asian and 20 South Asian donors. Test families
were excluded from the inference panels; both panels used a graph constructed from all
input assemblies. With PanGenie, the full panel increased exact recovery of assembly-derived
SV-bearing genotypes by 3.89 percentage points in East Asian donors and 1.97 points in South
Asian donors. Shared-SNV concordance increased by 0.047 and 0.023 points, respectively.
With Locityper, the full panel matched 31 of 39 experimental HLA genotypes, compared with
28 for HPRC-only, 34 for T1K and 36 for SpecHLA. The paired interval for the HLA panel gain
included zero. These results support expanded haplotype panels for MHC variant recovery
within a shared graph and motivate a larger independent evaluation of HLA typing.

# Introduction

HLA typing from short reads requires assigning reads to closely related genes and resolving
the two inherited alleles. HLA allele frequencies vary among populations [@Gourraud2014],
so a reference panel's population composition may affect which genotypes it can recover.
The surrounding major histocompatibility complex (MHC) also contains structural variation,
providing a second test of the value of additional reference haplotypes.

Different inference methods use reference sequences in different ways. PanGenie combines
short-read k-mer counts with a panel of haplotypes to genotype represented variants
[@Ebler2022PanGenie]. Locityper uses read alignment and depth to select pairs of locus
haplotypes [@Prodanov2025Locityper]. Specialist HLA methods include T1K, which estimates
allele abundances from read alignments [@Song2023T1K], and SpecHLA, which reconstructs
phased HLA sequences [@DeepOmicsSpecHLA]. An expanded panel can supply additional alleles,
but its effect on variant recovery and named HLA typing requires separate measurements.

We ask whether adding Asian and Arab haplotypes to an HPRC reference panel improves MHC
variant genotyping and HLA typing in East and South Asian donors. We constructed AsianPGR,
an MHC graph with 754 input haplotype entries, and compared full and HPRC-only inference
panels using the same reads, donor folds and genotyping software. We first measured
representation of held-out HLA labels, then exact variant-genotype recovery, and finally
HLA typing against experimental and assembly-derived labels. This design estimates the
combined effect of increasing panel size and changing its composition within a shared graph.

# Methods

## Haplotype panel and sequence annotation

We extended a 610-haplotype panel comprising the Arab Pangenome Reference (APR; 106 haplotypes)
[@Nassir2025APR], HPRC release 2 (464), JaSaPaGe Saudi (18), JaSaPaGe Japanese (20), GRCh38
and CHM13. We added 28 haplotypes from 14 K-PanRef individuals [@Shin2026KPanRef] and 116 from
58 Chinese Pangenome Consortium (CPC) individuals [@Wang2026CPC]. These additions used paths
from the projects' Minigraph-Cactus graphs. We extracted K-PanRef haplotype paths with
`vg paths -F`. For CPC, we used `odgi` to identify nodes in the CHM13 chromosome 6 interval
28–34 Mb and extracted each haplotype walk spanning those nodes. We retained segments with
at least 50 kb of matching bases and mapping quality (MAPQ) ≥20 in a minimap2 asm20 alignment
to the GRCh38 MHC [@Li2018minimap2]. Three CPC samples already represented in HPRC were
excluded from the additions.

The resulting panel contains 754 haplotype entries (Table \ref{tableCohorts}). Excluding the two
reference entries and five duplicate donor assembly pairs gives 742 haplotypes from 371
name-reconciled donors for donor-level analyses. Population metadata divide HPRC into Japanese
(JPT), other East Asian (CHB, CHS, CDX, KHV and HG005), HG002 (Ashkenazi), and remaining samples.
HG002 contributes two haplotypes from one individual. K-PanRef and CPC sequences inherit the
clipping applied during construction of their source graphs; this can shorten private insertions.

Table: Input haplotype entries used to construct AsianPGR. Counts precede donor reconciliation and fold-specific exclusions. \label{tableCohorts}

| Source | Haplotypes | Input |
| --- | ---: | --- |
| HPRC release 2 | 464 | Assemblies |
| Arab Pangenome Reference | 106 | Assemblies |
| JaSaPaGe Saudi / Japanese | 18 / 20 | Assemblies |
| K-PanRef | 28 | Graph paths |
| Chinese Pangenome Consortium | 116 | Graph paths |
| GRCh38 / CHM13 | 2 | References |
| **Total** | **754** | |

We used the Common Workflow Language workflow `hla_pangenome.cwl` to extract the extended MHC
with pgr-tk `pgr-query` [@Chin2023pgrtk] and annotate gene structures and allele identities with
Immuannot v3 [@Zhou2024Immuannot]. The extraction interval was GRCh38
chr6:28,510,120–33,480,577 with 100 kb flanks. Original annotations used IPD-IMGT/HLA v3.55,
IPD-KIR v2.13 and RefSeq C4. We regenerated per-gene sequences, pgr-tk bundle decompositions,
and pggb/odgi graphs for the enlarged panel [@Garrison2024pggb; @Guarracino2022odgi].
The subsequent coding-sequence audit used frozen IPD-IMGT/HLA v3.65.0 labels. Complete coding
sequences (CDS) required valid strand and interval extraction and passed completeness checks;
ambiguous or unmatched labels remained unresolved.

## MHC graph construction and visualisation

We built AsianPGR with Cactus v3.3.0 and the Minigraph-Cactus algorithm
[@Hickey2024MinigraphCactus], using all 754 entries and GRCh38 and CHM13 as reference paths.
The build requested clipped and full GFA/GBZ graphs, Giraffe indexes, VCF and odgi output,
and used 32 cores and 240 GB on the NIG supercomputer. A mapping check used MHC-recruited
HG00096 short reads and vg Giraffe v1.76.1 with eight threads [@Siren2021Giraffe].

For Figure \ref{figGraph}, we visualised the structural backbone from the construction
GFA and regenerated regional principal bundles with pgr-tk v0.6.0. The class II input
spans DRA to DMA with 20-kb flanks; HLA-A and HLA-DPB1 inputs include 2-kb gene flanks.
All three decompositions used SHIMMER parameters w=48, k=56 and r=2. Class II used
minimum span 8, bundle-length cutoff 500 bp and merge distance 2,000 bp; the gene
views used minimum span 12, cutoff 200 bp and merge distance 1,000 bp. Minimum
coverage was 0 and minimum branch size 8. We connected bundle identifiers observed
consecutively in each sequence and counted distinct carrying sequences per bundle.
The input sets contain 753 class II sequences, 753 HLA-A sequences and 754 HLA-DPB1
sequences. Frozen layout coordinates, input hashes and commands accompany the figure.

## Full versus HPRC-only panels in 40 donors

We selected 20 East Asian (EAS) and 20 South Asian (SAS) 1000 Genomes donors and assigned them
to five population-balanced folds with known families kept together. For each fold, we excluded
test donors, aliases and known relatives from both a full reference panel and an HPRC-only panel.
The panels contained approximately 363 and 224–225 training donors, respectively. Both used the
same whole-MHC graph, whose construction included the test assemblies. This design evaluates
held-out haplotype inference conditional on a shared graph topology.

We decomposed the graph into top-level variant sites relative to GRCh38 with `vg deconstruct`
and genotyped identical MHC-recruited reads with PanGenie v4.2.1. The initial panel builder
removed any site with a missing training genotype. We revised this rule to retain known phased
alleles and encode unknown alleles as missing, using PanGenie's native support for phased
missingness. Alternate alleles entered a panel only when observed in its training haplotypes.
Both panel arms used the revised rule and the original donor folds.

Variant truth comprised the test donors' assembly-derived allele sequences at graph sites within
GRCh38 chr6:28,510,121–33,480,577 (1-based, inclusive). We excluded missing, conflicting and
non-ACGT truth and required each allele's reference span to fit the interval. A correct genotype
required an exact match of the unordered pair of allele sequences. Missing sites and no-calls
counted as failures in recovery across all eligible truth genotypes. An SV-bearing genotype
contained at least one allele whose length differed from the reference by ≥50 bp; success required
the complete bubble allele pair, including embedded small variants. We also compared pure
single-nucleotide variant (SNV) sites shared by both original graph outputs and the published
NYGC GRCh38 PASS callset, keeping that comparison universe fixed after the repair. The NYGC
arm used the October 2020 filtered, phased release from joint whole-genome calling.

To obtain named HLA genotypes, we used Locityper v1.7.4 [@Prodanov2025Locityper] on intact
source locus sequences from the same training panels, with MHC-recruited reads and regional
background calibration. Predictions received a numeric two-field name when their frozen exact-CDS
labels resolved to one type. The assembly endpoint covered eight genes (A, B, C, DPA1, DPB1,
DQA1, DQB1 and DRB1). Experimental labels supplied 39 eligible genotypes at five genes in eight
EAS donors. T1K v1.0.6 and SpecHLA retained their frozen predictions and native databases;
SpecHLA used IPD-IMGT/HLA v3.38.0, while T1K used the run's current-IPD reference. T1K required
positive call quality and ambiguity that collapsed to a single numeric two-field allele.
Unresolved predictions counted as failures. We estimated paired differences with 10,000 donor
bootstrap resamples within each ancestry stratum, conditional on the fixed sites and folds.
The random seed was derived from SHA256 of the analysis date (20260917), endpoint and baseline;
intervals were unadjusted for multiple comparisons.

Table: Evaluation endpoints. All variant comparisons use 20 East Asian (EAS) and 20 South Asian (SAS) donors; experimental HLA donors form a subset of EAS. The graph topology includes test assemblies. \label{tableEndpoints}

| Endpoint | Truth and denominator | Success criterion |
| --- | --- | --- |
| SV-bearing genotypes | Assembly sequences; 1,465 EAS / 1,577 SAS genotypes | Exact diploid graph-bubble allele pair |
| Shared SNVs | Assembly genotypes at fixed three-way shared sites; 1,007,636 comparisons per stratum | Exact diploid genotype |
| Experimental HLA | Experimental typing; 39 genotypes, five loci, eight EAS donors | Exact unordered two-field pair |
| Assembly HLA | Complete-CDS labels; 159 EAS / 160 SAS genotypes, eight loci | Exact unordered two-field pair |

# Results

## An Asian-enriched graph of the MHC

AsianPGR contained 417,896 nodes and 579,137 edges, representing 5.89 Mb of sequence.
Of 752 non-reference input haplotypes, 747 covered at least 99% of the GRCh38 MHC interval;
720 were extracted as a single segment. Construction took 3 h 43 min. In the HG00096
mapping check, 654,059 of 657,052 recruited read pairs aligned, including 625,590 at
MAPQ ≥20. The input cohorts span the Arabian Peninsula and East Asia alongside the worldwide HPRC
panel (Figure \ref{figGraph}A). Their 376 sample entries represent 371 distinct donors;
five donors occur in both HPRC and JaSaPaGe. The two reference genomes bring the construction
panel to 754 haplotype entries.

Figure \ref{figGraph}B shows the Minigraph structural backbone used in graph construction.
To examine regional sequence organisation, we regenerated pgr-tk principal bundles locally
from the final-panel class II, HLA-A and HLA-DPB1 sequences. Each bundle represents a shared
sequence path, and edges connect bundles that occur consecutively along an input sequence.
The class II view contains 890 bundles and shows the organisation around DRB; HLA-A and
HLA-DPB1 provide complementary class I and class II views (Figure \ref{figGraph}C–E).

```{=latex}
\begin{figure}[p]
\makebox[\linewidth][r]{\includegraphics[width=1.25\linewidth,height=0.73\textheight,keepaspectratio]{figures/fig1_asian_graph.pdf}}
\caption{AsianPGR input cohorts and sequence-graph structure. A: sample counts by source for the mapped Asian and Arab cohorts, with the worldwide HPRC contribution shown separately. Five donors overlap HPRC and JaSaPaGe; counts by source therefore sum to 376 sample entries from 371 donors. GRCh38 and CHM13 contribute the two reference entries. Map: Natural Earth. B: the Minigraph structural backbone, with 1,677 sequence segments and 2,381 links; colours distinguish GRCh38 segments from added sequence. This backbone is refined into the 417,896-node Minigraph-Cactus graph. C–E: locally regenerated pgr-tk principal-bundle graphs for the class II region (DRA to DMA, including DRB), HLA-A and HLA-DPB1. Nodes are bundles; edges represent observed consecutive bundles, drawn without direction; colour gives the number of haplotype sequences carrying each bundle. Counts include reference sequences. Gene labels mark GRCh38 anchors. Graph layouts show connectivity; distances have no genomic scale.}
\label{figGraph}
\end{figure}
```

## Added haplotypes increase reference coverage of HLA labels

We first asked whether the added haplotypes supplied HLA labels absent from the HPRC-only
training panels. Both assembly-derived truth labels were represented for 156/159 eligible
East Asian genotypes with the full panel, compared with 151/159 for HPRC-only. Corresponding
South Asian counts were 155/160 and 149/160. Gains varied by locus, with the largest South
Asian gain at DPB1 (Suppl. Fig. \ref{figCoverage}). These counts measure two-field label coverage;
exact sequence recovery was assessed separately.

## The full panel improves recovery of SV-bearing genotypes

We next asked whether the added reference haplotypes improved inference from the same reads.
Both panels retained sites with partially missing training genotypes. This rule restored
19,043–19,072 sites per fold in the full panel and 221–227 in HPRC-only compared with the
initial complete-case filter. The comparison below uses the revised rule in both arms.

The full panel increased exact recovery of SV-bearing diploid genotypes from 72.35% to
76.25% in East Asian donors and from 75.33% to 77.30% in South Asian donors
(Figure \ref{figVariants}A–C). Recovery improved in 18/20 and 14/20 donors, respectively.
The paired gains were 3.89 percentage points (95% donor bootstrap interval 2.40–5.29) and
1.97 points (0.32–3.47). These results support additional haplotypes for recovering complete
allele pairs at SV-bearing graph sites.

Shared-SNV differences were smaller (Figure \ref{figVariants}D). On 1,007,636 fixed
comparisons per stratum, full-panel concordance was 99.8217% in East Asian donors and
99.8272% in South Asian donors. Gains over HPRC-only were 0.0467 points (0.0148–0.0828)
and 0.0233 points (−0.0060–0.0573). The published linear callset reached 99.6907% and
99.6405% on these same comparisons. This arm provides context from an existing joint-calling
pipeline; the paired panel contrast holds inference software and read inputs fixed.

![Variant recovery with full and HPRC-only panels. A–B: each line connects one donor's exact recovery of eligible SV-bearing genotypes; diamonds show pooled correct/eligible proportions, with counts above. An eligible genotype contains an allele differing in length from the reference by at least 50 bp. Success requires the complete unordered allele-sequence pair, including embedded small variants. Missing sites and no-calls count as failures. C–D: pooled full-minus-HPRC differences with 95% intervals from 10,000 paired donor bootstrap resamples. The shared-SNV endpoint uses 1,007,636 comparisons per ancestry stratum. Horizontal scales differ between C and D; pp denotes percentage points. EAS: East Asian; SAS: South Asian. \label{figVariants}](./figures/fig2_variant_recovery.pdf)

## HLA typing gains are small and depend on the truth source

We tested whether greater reference coverage also improved named HLA genotypes.
An initial whole-locus PanGenie representation yielded predominantly unresolved predictions;
its marker selection retained few informative k-mers across complete alleles. We therefore
used Locityper to compare intact locus haplotypes from the same training panels.

Against experimental typing in eight East Asian donors, full-panel Locityper recovered
31/39 genotypes, compared with 28/39 for HPRC-only, 34/39 for T1K and 36/39 for SpecHLA
(Figure \ref{figHLA}A). The full-minus-HPRC gain was 7.69 percentage points, with a wide
paired donor interval (−2.78 to 20.00). SpecHLA had the highest concordance on this endpoint.

Against assembly-derived labels, full-panel Locityper recovered 137/159 East Asian and
143/160 South Asian genotypes (Figure \ref{figHLA}B–C). Gains over HPRC-only were 1.89
points (−1.87 to 5.66) and 2.50 points (0.00 to 5.62). The full panel left seven and two
eligible genotypes unresolved because selected sequences lacked unambiguous complete-CDS
labels. DRB1 remained difficult: full-panel Locityper matched 13/20 East Asian and 16/20
South Asian assembly genotypes, while T1K matched 20/20 in each stratum.

The experimental and assembly endpoints therefore give different method rankings. Historical
experimental labels and current assembly labels were compatible in 36/38 overlapping eligible
genotypes; two DRB1 pairs disagreed. The separate truth sources and explicit unresolved-call
counts identify where additional typing and sequence validation would be informative.

![HLA typing against experimental and assembly-derived truth. Each bar uses all eligible genotypes within its panel and separates correct allele pairs, incorrect pairs and unresolved predictions. White labels give correct/eligible counts. Success requires an exact unordered numeric two-field pair. A: experimental typing at A, B, C, DRB1 and DQB1 in eight East Asian donors. B–C: assembly labels at eight loci in 20 East Asian and 20 South Asian donors. T1K and SpecHLA retain their frozen native reference databases. The experimental donors are a subset of the East Asian cohort. \label{figHLA}](./figures/fig3_hla_typing.pdf)

# Discussion

Expanding the HPRC reference panel with Asian and Arab haplotypes improved recovery of
assembly-derived SV-bearing genotypes in both tested ancestry strata. Shared-SNV gains were
small, and the named-HLA comparison gave uncertain incremental gains over HPRC-only.
The distinction matters for reference design: additional allele representation can improve
variant recovery while HLA typing still depends on locus inference and allele designation.

The comparison holds reads, inference software and folds fixed. Its scope is inference within
a graph constructed using all input assemblies, including test donors. Panel size and
population composition also change together. A size-matched comparison with graphs built
exclusively from training donors would isolate the contribution of ancestry composition and
measure performance on unseen graph variation. The present variant endpoint measures exact
assembly-derived bubble genotypes; independently validated SV events would provide a
complementary assessment of structural-variant detection.

Experimental HLA truth covers eight East Asian donors and five loci. A larger cohort with
experimental labels in South Asian and Arab donors would extend the population coverage and
narrow uncertainty in the panel comparison. Harmonising reference databases and read
recruitment would also clarify the differences among Locityper, T1K and SpecHLA.
Graph-derived source sequences inherit clipping, and unresolved complete-CDS labels limit
some named predictions. These are concrete targets for improving the reference and its
annotation before broader evaluation.

## Data and code availability

AsianPGR graph files and the frozen analysis archive are linked from the project release at
<https://github.com/biohackathon-japan/BH26-asian-hla/releases/tag/v0.1.0>.
The manuscript directory contains the score tables used for the figures; `scripts/plot_paper.py`
regenerates Figures 2–3 and Suppl. Fig. S1; `scripts/plot_graph_figure.py` regenerates Figure 1. Appendix A identifies the analysis sources and archived diagnostics.

## Acknowledgements

We thank the organisers of DBCLS BioHackathon Japan 2026 and the NIG supercomputer team;
the APR, K-PanRef, CPC, HPRC and JaSaPaGe teams for sequence resources; and the developers
of the annotation, graph and typing tools.

# References

```{=latex}
\AtEndDocument{%
```

# Appendix A: Analysis provenance

The versioned graph archive includes clipped and full GBZ graphs, Giraffe short-read indexes,
graph variants, contig-name mappings and the original build command. The frozen analysis
archive preserves extraction and annotation under `hla/`, graph construction details in
`hla/docs/MHC_GRAPH_README.md`, and the completed panel comparison in
`hla-asian50/refined/REPORT.md`. Donor manifests, exclusions, software parameters,
per-donor scores and machine-readable checksums accompany these analyses.

The six benchmark tables in `paper/data/` are unchanged copies of the corresponding refined analysis
results. They provide per-donor variant counts, variant summaries and paired intervals,
per-genotype HLA outcomes, locus summaries and paired HLA intervals. Plotting changes presentation
only; reported intervals are taken from the frozen analysis. Figure 1 has a separate
`paper/data/graph/` directory containing regenerated bundle intervals, topology,
layout coordinates, cohort counts and input provenance. Its inputs and regeneration
procedure are documented there.

Extraction diagnostics, assembly annotation audits and initial graph-consensus typing runs
remain in the analysis archive. The graph-consensus experiment used Giraffe, `vg call`,
WhatsHap [@Patterson2015WhatsHap] and consensus annotation with Immuannot. Its aggregate
comparisons use method-specific resolved-call subsets. A common donor–locus eligibility
set would be needed to compare overall typing success with the benchmark in Figure \ref{figHLA}.
The original whole-locus PanGenie HLA outputs and callable-site filtering diagnostics also
remain archived. Results use completed analyses available on 18 September 2026.

# Supplementary figure

```{=latex}
\setcounter{figure}{0}
\renewcommand{\thefigure}{S\arabic{figure}}
```

![Representation of held-out HLA labels in the training panels. At each locus, the points show the proportion of eligible donor genotypes for which both assembly-derived two-field labels occur in the corresponding training panel. Lines connect the two panel conditions; overlapping markers indicate equal coverage. Right-hand counts give eligible genotypes. A: East Asian donors. B: South Asian donors. Test donors and known relatives were excluded from the panels. \label{figCoverage}](./figures/figS1_reference_coverage.pdf)

```{=latex}
}
```
