# Figure and table argument

The question is whether expanding an HPRC reference panel with Asian and Arab
haplotypes improves short-read MHC genotyping in East and South Asian donors.
The study estimates the combined effect of panel size and composition within a
shared graph. Independent graph construction and size-matched panels remain
follow-up experiments.

## Main figures

| Figure | Question | Evidence and display | Inference |
| --- | --- | --- | --- |
| 1. Asian HLA graph | What populations and sequence structures does the resource represent? | Asia map with source-specific sample counts; actual Minigraph structural backbone; locally regenerated class II/DRB, HLA-A and HLA-DPB1 pgr-tk bundle graphs | Establishes the resource, its population sampling and regional graph structure |
| 2. Variant recovery | Does the full panel improve inference on the same reads? | Paired donor SV recovery; pooled counts; full-minus-HPRC effects and donor bootstrap intervals for SV-bearing genotypes and shared SNVs | Shows consistent SV-bearing genotype gains and their much smaller SNV counterpart |
| 3. HLA typing | Does reference expansion yield more correct HLA calls? | Correct, incorrect and unresolved genotype fractions; common denominator within each truth-source panel; experimental truth first | Shows the modest panel gain, specialist-method comparison and dependence on truth source |

Figure 1 describes the resource. Supplementary Figure S1 measures held-out HLA
label representation. Figure 2 tests exact allele-sequence pairs at graph sites. Figure 3 tests numeric two-field HLA names. Each evaluation endpoint has its own truth source and success rule. The panels are connected by the scientific
question, with those distinctions retained in the captions and Methods.

Benchmark source tables and bootstrap intervals are frozen. Figure 1 uses
locally regenerated pgr-tk bundles from final-panel sequences. Its input hashes,
commands, sequence counts and layouts are retained with the figure data. Figure PDFs retain selectable text so labels
can be checked alongside the manuscript.

## Main tables

- **Table 1: Reference composition.** Source, input type and haplotype count.
  Distinguish the 754 graph entries from 371 reconciled donors and the smaller
  fold-specific training panels. Remove population subcolumns whose only role
  was to support the old descriptive heatmaps.
- **Table 2: Evaluation definitions.** Truth source, denominator and success
  criterion for SV-bearing genotypes, shared SNVs, experimental HLA and assembly
  HLA. This is the reader's guide to what each numerical comparison means.

Exact per-locus and per-donor outcomes remain in the machine-readable source
  tables. Repeating the figure totals in additional main-text tables would add
  length without a new inference.

## Disposition of previous material

- Remove the population-frequency/C4/database-novelty composite. Its panels
  answer different questions, and one column represents a single donor.
- Remove the aggregate method-agreement charts. Agreement between methods and
  conditional accuracy on different resolved-call subsets cannot supply the
  common-denominator comparison needed here.
- Keep extraction, heterozygosity and annotation QC in the archived analyses;
  retain the construction checks needed to understand the resource in Methods.
- Remove the old candidate-new-allele figure, whose labels predate the CDS audit.
- Exclude DogoHLA development and RCCX/C4 method comparisons from this paper.
- Remove all superseded images from the manuscript figure directory. The Git
  history retains the previous versions.

## Figure hygiene

Use method names and defined scientific conditions in labels. Remove usernames,
local filesystem paths, machine names, run nicknames and conversational labels.
Check both visible figure text and searchable PDF text/metadata. Preserve public
sample accessions where scientifically necessary in the underlying data.

## Evidence gaps retained in the draft

- Test assemblies contribute to graph topology; inference panels exclude test
  families.
- Panel size and ancestry composition change together.
- Experimental HLA truth includes only eight East Asian donors; 39 genotypes
  are correlated measurements across five loci.
- Existing HLA methods use different frozen databases and preprocessing.
- The exact CPC graph accession and historical T1K database release still need
  to be recovered from source manifests before submission.

## Figure 1 panel specification

- **A: source cohorts.** Map APR (53), JaSaPaGe Saudi (9), JaSaPaGe Japanese
  (10), K-PanRef (14) and CPC (58) samples. Show HPRC's 232 worldwide samples
  separately. Across projects, 376 sample entries represent 371 distinct donors;
  five donors occur in both HPRC and JaSaPaGe. Add GRCh38 and CHM13 to obtain
  754 haplotype entries. Markers indicate cohort geography.
- **B: MHC structural backbone.** Draw the 1,677 segments and 2,381 links from
  the actual construction GFA. The caption distinguishes this Minigraph
  backbone from the final 417,896-node base-level graph.
- **C: class II / DRB region.** Regenerate the 753-sequence DRA-to-DMA pgr-tk
  decomposition locally. The resulting 890 bundles and 1,469 undirected
  adjacencies reproduce the scale of the final-panel slide graph.
- **D–E: HLA-A and HLA-DPB1.** Regenerate both locally from the final panel:
  753 / 754 sequences and 42 / 23 bundles. These give class I and class II
  regional examples under the same gene-level bundle parameters.

Regional node colour represents the number of carrying haplotype sequences.
Gene anchors come from the frozen GRCh38 annotation. Layout distances have no
physical genomic scale; comparisons of bundle counts must account for the
region lengths and decomposition parameters.
