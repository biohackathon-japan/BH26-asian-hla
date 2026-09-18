# Talk track (~3 minutes)

## 1 — Whose HLA is in your reference?

“We combined five projects into an MHC graph with 754 haplotype entries. The table
shows each project's sample and haplotype contributions. These include 462
Asian or Arab haplotype entries: 266 East Asian, 72 South Asian and 124 Arab.
We also have a reproducible extraction and annotation workflow and MHC reads
from the 2,504-sample 1000 Genomes cohort.”

Counts are source entries. APR contributes 53 samples, HPRC 232, JaSaPaGe 19
(nine Saudi and ten Japanese), K-PanRef 14 and CPC 58. They total 376 source
sample entries from 371 distinct donors, because five donors occur in both
HPRC and JaSaPaGe. Each sample entry contributes two haplotypes; GRCh38 and
CHM13 add two references. The former 390 Asian/Arab count omitted 72 South
Asian HPRC haplotypes. The current breakdown includes them explicitly.

## 2 — More haplotypes. Better SV calls.

“We evaluated the same reads with PanGenie and two reference panels in 20 East
Asian and 20 South Asian donors. The eligible graph truth includes 79,844 SNV
sites and 298 sites containing a structural allele per donor. Among the
SV-bearing truth genotypes, exact recovery increased from 72.4% to 76.2% in
East Asian donors and from 75.3% to 77.3% in South Asian donors.

There are 1,465 and 1,577 eligible SV-bearing donor genotypes, respectively.
The separate fixed shared-SNV comparison has 1,007,636 donor-site comparisons
per ancestry group. These are pilot results against assembly-derived truth.
Panel size and ancestry composition change together.”

SV-containing sites have at least one graph allele differing in length from
reference by at least 50 bp. An SV-bearing truth genotype has a long allele in
that donor. The figure scores the latter, so its denominator differs from all
298 sites multiplied by donor count. Exact success requires the whole unordered
allele-sequence pair, including embedded small variants. Missing calls count
as failures. The confidence intervals resample paired donors. Test families
were excluded from inference panels; test assemblies remain in graph topology.

Non-Asian sample and variant counts will be added after the ongoing comparison
completes and its denominators and outputs are verified. No interim non-Asian
numbers enter this deck.

## 3 — DogoHLA

“DogoHLA is a pangenome-native, population-specific HLA typer. Population panel
sequences guide read collection. We call and phase read-supported variants,
then use supported noncoding indels from locus graphs to refine the two gene
sequences and assign HLA names. The implementation uses SpecHLA components for
alignment and small-variant inference.

In eight development donors, exact global whole-gene matches increased from
36 of 128 with native SpecHLA to 57 of 128 with the selected DogoHLA candidate.
Total sequence edit distance decreased by 22.8%, and correct two-field genotypes
increased from 61 of 63 to 62 of 63. These are selected development results.
Further independent analysis is needed; the larger evaluation is ongoing.”

The comparison uses assembly truth, optimal pairing of the two haplotypes and
global whole-gene exactness, with N counted as a mismatch. These numbers replace
the earlier slide's 41/51 infix-matching comparison, which assessed a different
metric and an earlier candidate. Controls attribute most exactness gains to
phasing repair and updated phase references. The guarded graph step reduces
residual noncoding DRB1 errors in two donors; exact-gene counts remain 57/128
relative to the repaired-phasing control. See the evidence links in README.md.
