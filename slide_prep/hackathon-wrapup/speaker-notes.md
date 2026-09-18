# Speaker notes

## 1. Asian HLA pangenome

The construction panel combines HPRC, APR, JaSaPaGe, K-PanRef and CPC.
HPRC contributes 232 donors, including 16 Japanese JPT donors (32 haplotypes).
JaSaPaGe contributes 10 Japanese and 9 Saudi donors. Five donors occur in both
HPRC and JaSaPaGe: the table totals 376 project entries, or 371 distinct donors.
The 752 donor haplotype entries plus GRCh38 and CHM13 give 754 entries.
The Asian/Arab total is 462: 266 East Asian, 72 South Asian and 124 Arab.
Donor location labels describe source populations, not current residence.

The map, construction backbone and class II/DRB bundle graph use the manuscript
Figure 1 data. The construction backbone has 1,677 segments and 2,381 links;
the final Minigraph-Cactus graph has 417,896 nodes and 579,137 edges.
The regional graph was regenerated locally with pgr-tk from 753 sequences.
Layouts show adjacency; graph distances are not genomic coordinates.
The separate 1000 Genomes read resource contains 2,504 samples.

## 2. Variant recovery compared with HPRC

The same PanGenie caller and reads were used with expanded and HPRC-only panels
in 20 East Asian and 20 South Asian donors. Points show the pooled accuracy
difference, with 95% intervals from paired donor bootstrap resampling.
SV-bearing genotypes improve by 3.89 and 1.97 percentage points. Shared SNVs
improve by 0.047 and 0.023 points. The SNV axis is expanded; its gains are much
smaller, and the South Asian interval includes zero.

The SV endpoint requires exact recovery of the unordered whole allele pair
for donors carrying an allele at least 50 bp different in length from reference.
Its eligible denominators are 1,465 EAS and 1,577 SAS genotypes. The fixed
shared-SNV endpoint has 1,007,636 donor-site comparisons per ancestry group.
The broader truth includes 79,844 SNV sites and 298 SV-containing sites per
donor; these counts are kept here to leave the slide focused on the comparison.

Non-Asian variant recovery has not yet been evaluated. The non-Asian row has
no estimate. The ongoing DogoHLA typing experiment measures another endpoint.
Panel size and ancestry composition change together. Test families are excluded
from panels, while test assemblies remain in the graph topology.

## 3. DogoHLA

Dawn presents the method. This slide contains only her original DogoHLA image,
retained at its native resolution from the org repository.
