---
title: 'Evaluating Asian-specific pangenomes for HLA typing'
title_short: 'BH26JP: Asian pangenomes for HLA typing'
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

Human leukocyte antigen (HLA) typing from short reads requires distinguishing closely related
genes and resolving highly polymorphic alleles. We present DōgoHLA, a new population-specific
method for HLA sequence reconstruction built around an Asian- and Arab-enriched pangenome.
We constructed a 754-haplotype major histocompatibility complex (MHC) panel and used its
sequences for read collection and its locus graphs for structural inference. DōgoHLA combines
these pangenome inputs with read-backed phasing and guarded reconstruction of graph-supported
noncoding indels. In eight development donors, DōgoHLA reconstructed 57 of 128 whole-gene
haplotypes exactly, compared with 36 for native SpecHLA, and reduced total sequence edit
distance by 22.8%. To assess the reference resource separately, we compared full and HPRC-only
panels in 40 East and South Asian donors. After correcting a callable-site filter, PanGenie
with the full panel recovered more assembly-derived SV-bearing genotypes, with gains of 3.89
and 1.97 percentage points in the two strata. For named HLA genotypes, full-panel Locityper
matched 31 of 39 eligible experimental genotypes, compared with 28 for HPRC-only, 34 for T1K
and 36 for SpecHLA. These complementary experiments assess the population-specific reference
and the reconstruction method. Evaluation of the frozen DōgoHLA method on additional donors
was ongoing at the report cutoff.

# Introduction

Human leukocyte antigen (HLA) genes encode molecules that present peptides to the immune system.
HLA variation contributes to transplant compatibility and susceptibility to immune-mediated disease.
Allele frequencies vary among populations, making population representation relevant to the design
and evaluation of HLA reference resources [@Gourraud2014]. The major histocompatibility complex
(MHC) also contains duplicated genes and structural variation. An HLA typing method must assign
reads to the correct locus, distinguish the two inherited alleles and, for full-sequence typing,
reconstruct each allele across coding and noncoding regions.

Existing methods address different parts of this problem. HLA\*LA projects read alignments onto a
population reference graph to infer HLA types [@Dilthey2019HLALA]. T1K estimates HLA and KIR
allele abundances from reads aligned to allele references [@Song2023T1K]. SpecHLA assigns reads
to HLA loci, uses local assembly to improve alignment in divergent regions, and phases variants
to reconstruct diploid gene sequences [@DeepOmicsSpecHLA]. Pangenome methods provide a
complementary use of assembled haplotypes: PanGenie combines k-mer counts with a reference
haplotype panel to genotype variants [@Ebler2022PanGenie], while Locityper selects pairs of
locus haplotypes using read alignment and depth [@Prodanov2025Locityper]. These approaches
make the choice of reference sequences and their representation part of the inference problem.
A larger panel can add informative alleles while also changing which sites remain callable.

We developed DōgoHLA as a new, natively pangenome-based method for population-specific HLA
typing. Its reference represents Asian and Arab haplotypes, and population-specific sequences
enter both read collection and structural inference. DōgoHLA reconstructs the two gene
haplotypes by combining read-backed coding variation with supported noncoding changes from
locus graphs. This connects the population reference to sequence reconstruction at two stages
of the method.

First, we combined assemblies and graph-derived sequences from several pangenome projects,
annotated HLA genes and constructed AsianPGR, the population-specific MHC reference. We evaluated alternative uses
of this resource through graph-consensus typing, full versus HPRC-only panel comparisons in
40 East and South Asian donors, and RCCX/C4 structural typing. We then evaluated DōgoHLA's
whole-gene reconstruction in eight development donors, with controls for phasing, reference
database and graph inference. This report presents the completed panel comparisons and
DōgoHLA development results. Evaluation of the frozen method on 32 additional donors was
ongoing at the analysis cutoff on 18 September 2026. DōgoHLA and AsianPGR are released together at
<https://github.com/biohackathon-japan/BH26-asian-hla/releases/tag/v0.1.0>, with method source,
graph indexes, checksums and the frozen analysis artifacts described in the appendix.

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

Table: Haplotypes in the panel, grouped by source and population stratum. \label{tableCohorts}

| Cohort / stratum | Haplotypes | Population | Source |
| --- | ---: | --- | --- |
| APR | 106 | Arab (UAE) | Assembly |
| HPRC r2, Japanese | 32 | Japanese (1000G JPT) | Assembly |
| HPRC r2, HG002 | 2 | Ashkenazi, one donor | Assembly |
| HPRC r2, other East Asian | 70 | CHB/CHS/CDX/KHV, HG005 | Assembly |
| HPRC r2, remaining samples | 360 | Mixed | Assembly |
| JaSaPaGe, Saudi | 18 | Arab (Saudi Arabia) | Assembly |
| JaSaPaGe, Japanese | 20 | Japanese (1000G JPT) | Assembly |
| K-PanRef | 28 | Korean | Graph path |
| CPC | 116 | Chinese | Graph path |
| GRCh38 / CHM13 | 2 | Reference | Reference |
| **Total** | **754** | | |

We used the Common Workflow Language workflow `hla_pangenome.cwl` to extract the extended MHC
with pgr-tk `pgr-query` [@Chin2023pgrtk] and annotate gene structures and allele identities with
Immuannot v3 [@Zhou2024Immuannot]. The extraction interval was GRCh38
chr6:28,510,120–33,480,577 with 100 kb flanks. Original annotations used IPD-IMGT/HLA v3.55,
IPD-KIR v2.13 and RefSeq C4. We regenerated per-gene sequences, pgr-tk bundle decompositions,
and pggb/odgi graphs for the enlarged panel [@Garrison2024pggb; @Guarracino2022odgi].
The subsequent coding-sequence audit used frozen IPD-IMGT/HLA v3.65.0 labels. Complete coding
sequences (CDS) required valid strand and interval extraction and passed completeness checks;
ambiguous or unmatched labels remained unresolved.

## Whole-MHC graph and consensus-based typing

We built AsianPGR, a whole-MHC graph from all 754 entries with Cactus v3.3.0 and the Minigraph-Cactus
algorithm [@Hickey2024MinigraphCactus], using GRCh38 and CHM13 as reference paths. The command
requested clipped and full GFA/GBZ graphs, Giraffe indexes, VCF, full odgi output and
visualisations. The build used 32 cores and 240 GB on the NIG supercomputer. We assessed mapping
on MHC-recruited short reads from HG00096 using vg Giraffe v1.76.1 with eight threads
[@Siren2021Giraffe].

For an initial typing comparison, we aligned 1000 Genomes reads with Giraffe, called variants
with `vg call`, phased them with WhatsHap [@Patterson2015WhatsHap], and generated two haplotype
consensus sequences with `bcftools consensus`. We annotated these sequences with Immuannot
using IPD-IMGT/HLA v3.55.0 and v3.65.0. Two separate T1K runs provided direct-read comparisons.
The analysis used the completed outputs available in the archived snapshot of the 2,504-sample
cohort; 955 cohort members overlap the experimental HLA resource of Gourraud et al.
[@Gourraud2014].

We compared unordered allele pairs at one-field and numeric two-field resolution. Pairwise
comparisons covered A, B, C, DRA, DRB1, DQA1, DQB1, DPA1 and DPB1; comparisons with experimental
typing covered A, B, C, DRB1 and DQB1. We removed Immuannot's unresolved `:new` suffix before
truncating names and included only calls with sufficient numeric fields in conditional
concordance. Each method therefore has its own denominator in this initial analysis. For the
assembly annotation benchmark, we retained the experimental ambiguity lists and required two
eligible current-CDS labels per donor and locus.

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

## RCCX/C4 structural typing

We evaluated two short-read prototypes for the RCCX region, which contains C4A/C4B copy-number
and long/short variation. The first fitted pairs of reference paths to sampled canonical 31-mer
counts in 18 pilot and 88 additional donors. Markers required support in at least three training
families, MHC locus specificity and base quality ≥20; overlapping mates contributed once per
fragment. A correction for total-copy dosage was estimated on the 18 pilot donors and applied
to the 88 additional donors. Calibrated ordinary reference depth supplied a simpler comparator.
We scored coarse DRB gene-content signatures separately.

The second prototype used 299 targeted probes after development on the preceding 106 donors
and was evaluated on 24 additional donors. Test families were excluded from probe discovery,
training statistics and candidate paths. The probes measured C4 A/B and long/short marginal
dosages; a reference-pair model then inferred complete structural signatures. Comparators were
the earlier sketch, calibrated reference depth and C4Investigator [@Marin2024C4Investigator].
C4Investigator's total-copy calibration used the 18 pilot donors. Structural truth comprised
assembly annotations, including provisional CYP21/TNX prototype assignments.

## DōgoHLA: population-specific pangenome inference

DōgoHLA v0.1.0 takes short reads and a population-specific panel of annotated MHC haplotypes
as inputs and returns reconstructed diploid HLA gene sequences and allele names. We built its
reference from the Asian- and Arab-enriched panel described above, excluding test donors and
known relatives. Panel sequences guide read collection, and locus graphs supply the represented
structural alleles used during reconstruction. These are integral inputs to DōgoHLA's inference.

We implemented read alignment, small-variant calling and initial sequence reconstruction using
components from SpecHLA v1.0.12 [@DeepOmicsSpecHLA]. We corrected phase-alternative scoring
to group complementary haplotype scores under the same phase alternative and used
IPD-IMGT/HLA v3.65 phase references. For DRB1 structural inference, we preserved population
haplotype paths in a normalised pggb graph, aligned reads with Giraffe and genotyped represented
alleles. We decomposed confident homozygous graph alleles into primitive long indels with
affine-gap alignment and applied eligible noncoding changes while preserving read-phased coding
variants. Graph calls required PASS status, genotype and variant quality ≥20, depth ≥10 and
at least five reads supporting the alternate allele. Candidate indels were ≥50 bp and required
valid reference sequence, mappable boundaries and unique flanks. The naming step used the
reconstructed sequences with the fixed 11-kb DRB1 query crop removed. Frozen thresholds and
code snapshots are retained with the experiments.

## DōgoHLA evaluation

Development used eight fold-0 donors, with their families excluded from reference inputs.
We compared reconstructed sequences with assembly-derived truth at eight genes per donor,
assigning the two haplotypes to minimise total global edit distance. Exactness required a global
whole-gene match; masked N bases counted as mismatches. We separately scored eligible two-field
genotypes. Controls compared phasing changes, phase-reference databases, graph reconstruction
and naming-only changes. The prospective evaluation froze the method before inference on the
remaining 32 donors, with native SpecHLA and a DōgoHLA arm without graph reconstruction as
comparators. This evaluation was incomplete at the report cutoff.

# Results

## MHC recovery and panel composition

We first assessed how much of the reference MHC interval was recovered from each input.
Of 752 non-reference haplotypes, 747 reached ≥0.99 coverage and 720 were represented by a single
extracted segment (Figure \ref{figExtraction}). Four CPC haplotypes and one JaSaPaGe Saudi
haplotype had lower coverage, with minima of 0.92 and 0.96, respectively. These measurements
quantify coverage of the GRCh38 interval. The source clipping inherited by CPC and K-PanRef
also affects the representation of sequence outside that interval's alignment.

![MHC extraction across 752 non-reference haplotypes. Left: coverage of the GRCh38 MHC interval, with counts at coverage ≥0.99. Right: extracted segments per haplotype. CPC and JaSaPaGe Saudi contain the five lower-coverage haplotypes; CPC and K-PanRef are the graph-derived input cohorts. \label{figExtraction}](./figures/fig1_mhc_extraction_coverage.png)

Allele frequencies differed among the nine source/population strata (Figure \ref{figPopulation}).
For example, HLA-A\*11:01 occurred in 29% of CPC and 26% of other East Asian haplotypes,
while A\*24:02 occurred in 35–41% of the Japanese strata. These are frequencies among sampled
panel haplotypes. At DRB1, 57–79% of annotated copies across strata lacked an exact full-length
match to IPD-IMGT/HLA v3.55, motivating the subsequent current-database and CDS audits.

![HLA annotations across nine strata. Top: HLA-A, HLA-B and HLA-DRB1 allele frequencies. Bottom: secondary DRB genes, C4A/C4B long/short forms and the fraction of copies lacking a full-length database match. The HPRC-Jewish column contains HG002 alone. \label{figPopulation}](./figures/fig2_population_hla.png)

We used external 1000 Genomes SNP calls to distinguish near-homozygosity from duplicated
assembly content (Figure \ref{figHomozygosity}). NA18976 and NA19909 carried fewer than 100
heterozygous SNPs per 100 kb across the MHC, supporting near-homozygosity. In contrast,
NA18952 had heterozygous read-derived calls across the region but only ten substitutions between
its two JaSaPaGe haplotypes. This supports duplication of one haplotype during assembly; the
donor-level panel retained the HPRC assembly for this individual.

![MHC heterozygosity from 1000 Genomes data compared with assembly haplotypes. Low external heterozygosity supports near-homozygosity in NA18976 and NA19909; discordance between external heterozygosity and the two JaSaPaGe sequences identifies NA18952 for assembly QC. \label{figHomozygosity}](./figures/fig3_mhc_homozygosity.png)

## Coding-sequence audit and experimental HLA concordance

We checked assembly-derived HLA labels against current database sequences and experimental
typing. HG02717's HLA-DQB1 sequence, previously labelled as a candidate new allele, matched
registered DQB1\*02:180. The NA20346 HLA-DPA1 candidate had local support from 24 of 51 reads in the
initial analysis. The broader complete-CDS audit identified 47 distinct protein
candidates absent from the frozen IPD-IMGT/HLA v3.65.0 database. Three candidates, NA18620
HLA-C and HG02976 and NA19159 HLA-DRB1, had at least five supporting fragments at each tested
distinguishing SNP. These observations support the tested local bases; full-allele phasing and
orthogonal sequence validation remain steps in candidate assessment.

Current exact-CDS labels agreed with the Gourraud experimental genotypes in 62/64 HLA-A,
64/64 HLA-B, 58/63 HLA-C, 59/63 HLA-DRB1 and 58/61 HLA-DQB1 comparisons. Ten of the
14 discrepancies had possible antigen-binding-exon compatibility under frozen G-group
definitions. Among the remaining four, NA18943 A/DRB1 and NA19007 A had independent local
read or cross-assembly support for the assembly bases; NA18608 DRB1 remained unresolved because
several probes were non-unique. A complementary comparison between FuFiHLA
[@Hu2026FuFiHLA] and Immuannot assessed agreement of annotations on the same assemblies
(Figure \ref{figTypingConcordance}). Together, these analyses establish the eligible sequence
labels and identify discrepancies for further review.

![Agreement between FuFiHLA, Immuannot and T1K at one- to four-field resolution on overlapping assembled donors. Denominators depend on the method pair, locus and available resolution. \label{figTypingConcordance}](./figures/fig7_typing_concordance.png)

## Graph-consensus typing resolves a limited subset of eligible genotypes

The whole-MHC graph contained 417,896 nodes and 579,137 edges, representing 5.89 Mb of sequence,
and took 3 h 43 min to construct. In the HG00096 mapping test, 654,059 of 657,052 read pairs
aligned (99.5%), including 625,590 at MAPQ ≥20 (95.2%), in 26 s. This established a working
mapping resource for the typing experiments.

We then compared HLA calls from graph-derived consensus sequences with direct-read calls.
Pairwise agreement varied with method, locus and resolution (Figure \ref{figGraphOverall}).
The two database annotations of identical consensus sequences also differed, particularly at
DRB1, showing that allele designation depends on the database as well as the reconstructed
sequence. These comparisons motivate separate evaluation of sequence reconstruction and naming.

![Pairwise HLA genotype agreement between two T1K runs and Immuannot annotations of graph-derived consensus sequences using two database releases. The plotted comparison pools nine HLA genes. \label{figGraphOverall}](./figures/fig5_gc_overall_concordance.png)

Against experimental typing, T1K had 98.8% genotype concordance among 1,843 checkable calls,
the separate T1K run had 96.5% among 4,308, and each Immuannot database annotation
had 98.5% among 613 checkable calls
(Figure \ref{figGraphTruth}). The archived analysis reports two-field resolution for approximately
90% of T1K's eligible calls and 33% of the graph-consensus calls. These conditional concordances
refer to method-specific subsets. The practical finding is the larger number of resolved genotypes
from T1K, which motivated evaluation of alternative uses of the haplotype panel.

![Conditional two-field concordance with published experimental typing at five genes. Each bar uses that method's resolved calls: 1,843 for the pipeline T1K run, 4,308 for the separate T1K run and 613 for each Immuannot annotation. \label{figGraphTruth}](./figures/fig6_gc_truth_concordance.png)

## Retaining callable sites improves recovery with the full panel

The initial PanGenie analysis exposed a loss of callable sites when the panel builder required
complete training genotypes. Allowing phased missing alleles restored 19,043–19,072 sites per
fold in the full panel and 221–227 in HPRC-only. Every previously retained site remained, and
every repaired HPRC-only site was present in the repaired full panel.

With missing sites and no-calls counted as failures, the repaired full panel recovered more
SV-bearing assembly genotypes than repaired HPRC-only in both ancestry strata
(Table \ref{tableVariants}). The paired gain was 3.89 percentage points in EAS (95% donor
bootstrap interval 2.40–5.29) and 1.97 points in SAS (0.32–3.47). This supports the use of
additional reference haplotypes for exact bubble-genotype recovery within the shared graph.

Table: Exact recovery of SV-bearing diploid bubble genotypes after callable-site repair. \label{tableVariants}

| Stratum | Eligible genotypes | HPRC-only | Full panel |
| --- | ---: | ---: | ---: |
| EAS | 1,465 | 1,060 (72.35%) | 1,117 (76.25%) |
| SAS | 1,577 | 1,188 (75.33%) | 1,219 (77.30%) |

On the fixed 1,007,636 shared SNV comparisons per stratum, concordance was 99.8217% for the
full panel, 99.7749% for HPRC-only and 99.6907% for the published linear callset in EAS;
the corresponding SAS values were 99.8272%, 99.8039% and 99.6405%. The full-minus-HPRC
paired intervals were 0.0148–0.0828 percentage points in EAS and −0.0060–0.0573 in SAS.
These small shared-site differences complement the larger recovery gains from retaining
previously omitted sites.

## Named HLA performance depends on the reference and truth source

We evaluated whether the expanded panel also improved named HLA genotypes. An initial
whole-locus PanGenie representation produced predominantly unresolved calls because its marker
selection left few informative k-mers shared across complete alleles. Locityper provided a
usable locus-haplotype inference method for the same panels. With Locityper, the full panel
recovered more eligible HLA genotypes than HPRC-only on both assembly-derived and experimental
labels (Table \ref{tableHLA}).

Table: Exact unordered numeric two-field genotypes. No-calls count as failures. Experimental truth covers five loci in eight EAS donors; assembly truth covers eight loci in 20 donors per stratum. \label{tableHLA}

| Method | Experimental EAS | Assembly EAS | Assembly SAS |
| --- | ---: | ---: | ---: |
| Full panel + Locityper | 31/39 | 137/159 | 143/160 |
| HPRC-only + Locityper | 28/39 | 134/159 | 139/160 |
| T1K | 34/39 | 134/159 | 120/160 |
| SpecHLA | 36/39 | 126/159 | 110/160 |

SpecHLA had the highest concordance on the experimental endpoint. The full panel's gain over
HPRC-only was 7.69 percentage points, with a paired donor interval of −2.78 to 20.00. The
assembly-derived endpoint favoured full-panel Locityper, with smaller gains over HPRC-only:
1.89 points in EAS (−1.87 to 5.66) and 2.50 in SAS (0.00 to 5.62). The truth-source audit
found compatible historical and current assembly labels in 36 of 38 overlapping eligible
genotypes, with two DRB1 disagreements. The differing method rankings therefore need to be
interpreted with the endpoint and reference database specified.

## C4 dosage is more readily recovered than complete structural pairs

The first structural prototype recovered coarse DRB gene content in all 106 donors. Its
uncorrected RCCX path and binary-marker variants recovered 38/106 and 64/106 complete
structural pairs, respectively. After pilot-only dosage calibration, the corrected method
recovered total copy number in 88/88 additional donors and complete signature pairs in 73/88.
Calibrated ordinary reference depth also recovered total copy number in 88/88.

In the subsequent 24-donor experiment, the targeted assay, calibrated reference depth and
pilot-calibrated C4Investigator each recovered total C4 copy number in 24/24 donors.
The targeted assay recovered 21/24 complete structural pairs, compared with 22/24 for the
earlier sketch. These results distinguish accurate total dosage from the more difficult
assignment of C4 forms to complete reference structures. The observed errors involved incorrect
A/B dosage, ambiguous linkage between A/B and long/short forms, and a structure absent from
the reference panel.

## DōgoHLA improves whole-gene reconstruction in development donors

We evaluated whether DōgoHLA improved whole-gene reconstruction by comparing it with native
SpecHLA on the same eight development donors. The selected DōgoHLA candidate increased exact
whole-gene matches from
36/128 in the saved native SpecHLA runs to 57/128. Total global edit distance decreased from
97,965 to 75,634 (22.8%), and correct two-field genotypes increased from 61/63 to 62/63.
Edit distance improved in all eight donors. These values use global whole-gene alignment
throughout.

The controls localised the improvements. With panel-assisted read collection held fixed,
complete phasing repair and updated phase references increased exact matches from 43/128 to
57/128. Adding panel sequences to the phase-reference database left exactness at 57/128.
The final guarded graph-indel step reduced total edit distance from 76,924 to 75,634 while
retaining 57/128 exact sequences and 62/63 correct genotypes. The reduction came from supported
noncoding DRB1 insertions in two donors. An earlier whole-block replacement had overwritten
coding differences and reduced two-field correctness to 60/63; decomposing supported changes
into eligible noncoding indels preserved the read-phased coding alleles. The selected method
therefore combines phasing repair with a narrower structural reconstruction step. Its prospective
32-donor evaluation will assess these choices with the method fixed before inference.

# Discussion

DōgoHLA introduces a population-specific HLA reconstruction method in which the pangenome
supplies both read-collection sequences and structural alleles. We built an Asian- and
Arab-enriched MHC reference, combined its locus graphs with read-backed haplotype inference,
and evaluated the resulting sequences against assembly-derived truth. In the eight development
donors, DōgoHLA increased exact whole-gene recovery and reduced sequence edit distance relative
to native SpecHLA. This provides an initial evaluation of the method with its population-specific
reference.

The component controls identify how the method produced these gains. With panel-assisted read
collection held fixed, phasing repair and updated phase references increased exact whole-gene
matches. The guarded graph step then reduced residual sequence errors in two DRB1 donors while
preserving coding variation and allele names. DōgoHLA combines these operations in a pangenome-based
inference workflow, using SpecHLA components for read alignment, small-variant calling and
initial reconstruction. A complete paired evaluation of DōgoHLA and its graph ablation will
quantify the structural step's contribution on additional donors.

The separate PanGenie and Locityper comparisons assess the reference resource. The strongest
evidence for the full panel came from recovery of SV-bearing assembly genotypes after
callable-site repair. Added haplotypes contribute useful evidence when the panel representation
preserves the sites at which they differ. The smaller and less certain named-HLA gains also
show why variant recovery, allele designation and whole-gene reconstruction require distinct
endpoints. These comparisons connect the design of the population reference to the different
inference tasks that use it.

The 40-donor panel comparison is conditional on a graph built with the test assemblies. The
full panel is also larger than HPRC-only, so its observed effect combines panel size and
composition. A size-matched augmentation experiment using graphs constructed exclusively from
training donors would separate these effects. Experimental HLA labels are available for only
eight EAS donors in this benchmark; additional experimental typing in SAS and Arab donors
would extend the ancestry coverage. Harmonised allele databases and read recruitment would
also make comparisons among typing methods easier to interpret. The published linear arm
provides a shared-SNV comparison under its original joint-calling and filtering design.

Sequence quality and truth definition impose further constraints. Graph-derived inputs inherit
clipping, some assembly haplotypes require QC, and historical HLA labels can resolve a different
sequence region from a complete-CDS annotation. The candidate catalogue therefore retains
local read evidence separately from complete-allele validation. Similarly, C4 dosage and
reference-signature imputation describe different levels of structural information. Long reads
spanning diagnostic sites and module boundaries would provide an independent test of the
inferred RCCX arrangements. These follow-up experiments can establish when the expanded panel
improves inference beyond the fixed-panel comparisons reported here.

## Acknowledgements

We thank the organisers of DBCLS BioHackathon Japan 2026 and the NIG supercomputer team for
compute access; the APR, K-PanRef, CPC, HPRC and JaSaPaGe teams for the sequence resources;
and the developers of the annotation, graph and typing tools used in this work.

# References

```{=latex}
\AtEndDocument{%
```

# Appendices

## Analysis provenance and reproduction

DōgoHLA v0.1.0 and AsianPGR v0.1.0 are available from the BioHackathon project repository
(<https://github.com/biohackathon-japan/BH26-asian-hla>) and its versioned release
(<https://github.com/biohackathon-japan/BH26-asian-hla/releases/tag/v0.1.0>). The release provides
DōgoHLA source and development evidence. The AsianPGR graph and mapping indexes and the
frozen analysis archive are hosted at <https://bio2vec.net/data/asianpgr/v0.1.0/>, alongside
the JaSaPaGe data, and linked from the GitHub release. Machine-readable checksums identify the
individual files and release archives. The software release preserves the prepared-environment
implementation used in the experiments; its README specifies dependencies and input preparation.

The AsianPGR archive contains clipped and full GBZ graphs, Giraffe short-read indexes, graph
variants, contig-name mappings, input sizes and the original build command. The source manifest
links these files to the 15 September 2026 graph build. The graph archive preserves original
filenames for compatibility with the analysis scripts.

The analysis archive retains the original directory layout. Extraction and annotation are in
`hla/`; graph construction is documented in `hla/docs/MHC_GRAPH_README.md`; the initial
consensus comparison uses `hla/analysis/gc_*.py` and `hla/results/tables/gc_*.tsv`.
The assembly annotation audit is in `hla-analysis/FINAL_REPORT.md` and its reproduction guide.
The completed 40-donor comparison is in `hla-asian50/refined/REPORT.md`; the original comparison
is retained in `hla-asian50/accuracy/REPORT.md`. Structural protocols and measurements are in
`hla-structural/` and `hla-targeted/`. DōgoHLA source, development controls and the prospective
protocol are in `workflow/hla-spechla-pg/` in the repository and source archive.

The analysis directories retain donor manifests, exclusions, parameter settings, per-donor
scores and input/output checksums. External allele databases, raw reads and large source
assemblies are obtained from their original providers using the provenance and preparation
instructions. Results in this report use the completed analyses available at the 18 September
2026 cutoff. The release retains the frozen eight-donor DōgoHLA development results and the
protocol for its subsequent evaluation.

```{=latex}
}
```
