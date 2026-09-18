# Manuscript revision, 18 September 2026

This revision follows the request to focus on AsianPGR, variant genotyping and
HLA typing. It supersedes the earlier method-centred manuscript review.

The paper now asks whether an expanded reference panel improves inference over
HPRC-only panels in East and South Asian donors. DogoHLA and RCCX/C4 development
results have been removed from the manuscript. The reference construction and
completed 40-donor comparison supply the argument.

Three main figures show graph construction and topology, paired variant
recovery and HLA typing outcomes. Supplementary Figure S1 shows held-out HLA
label coverage. Figure 1 combines an Asia map, the construction backbone and
locally regenerated pgr-tk views of class II/DRB, HLA-A and HLA-DPB1. Two tables describe reference composition and
define evaluation endpoints. The figure rationale is in `FIGURE_PLAN.md`.
Six unchanged score tables and their provenance accompany the plotting script.
All numerical intervals are taken from the completed analysis.

Superseded figure assets have been removed from the manuscript directory,
including charts containing internal run labels. Replacement figure labels use
scientific names. The previous assets remain in Git history. Repository software
and historical releases have not been altered by this manuscript edit.

## Remaining submission details

- Authors, affiliations and the author footer remain empty as previously requested.
- Recover the exact CPC graph accession and historical T1K database release.
- The current paired comparison uses test-family exclusions from inference
  panels within a graph that includes test assemblies. Panel size and ancestry
  composition change together. Both limits remain explicit.
- Experimental HLA labels cover eight East Asian donors. Existing specialist
  methods retain their original reference databases and read processing.
- The public release still describes the earlier combined software/resource
  scope. Its release page and archives are historical artifacts; update public
  manuscript presentation when this revision is ready to publish.

## Build

The local working PDF uses a temporary copy of the BioHackrXiv generator with
its author minimum disabled and an empty author block supported. The repository
CI generator still requires author metadata. Figure and manuscript verification
results are recorded below after the build.

## Verification

- Rebuilt all four figures with Matplotlib 3.10.1 and visually checked them.
- Verified all variant and HLA summary counts against the detailed source rows.
- Verified all six copied TSV files byte-for-byte against the original results.
- Checked citation keys, image paths and all figure/table references.
- Built the revised PDF and inspected the title, tables and figure pages.
- Checked the final manuscript PDF and all figure PDFs for internal usernames,
  local paths, conversational run labels and obsolete method text, including
  PDF metadata. No matches or unresolved reference markers remained.
- `git diff --check` passed. Vector PDF and SVG exports accompany all figures and reused slide panels.
