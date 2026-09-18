# Manuscript review, 18 September 2026

Edited `paper.md` against the completed analyses in
`/home/leechuck/Public/software/pangenome` (tracked source snapshot
`01483dee04b6da61c313bba22f3f7752f3be7686`). The validation monitor was active
and its local outputs were incomplete when reviewed. No interim validation
scores enter the revised manuscript. Analysis code and result files were not changed.

## Scientific changes

- Replaced the outdated claim that the 40-donor benchmark awaited scoring with
  the completed callable-site repair, PanGenie variant comparison and Locityper
  named-HLA comparison. Sources: `hla-asian50/refined/REPORT.md`,
  `results/named_HLA_summary.tsv`, and the original `accuracy/REPORT.md`.
- Defined the graph-derived SV-bearing bubble-genotype endpoint, no-call rules,
  shared-SNV comparison, donor bootstrap, and experimental versus assembly HLA
  endpoints. Made inclusion of test assemblies in graph construction and the
  panel-size/composition confounding explicit.
- Added DōgoHLA's selected eight-donor development result and controls from
  `hla-spechla-pg/RESULTS-2026-09-18.md`. The native comparison is 36/128 exact
  genes; 43/128 belongs to an already panel-augmented control. The graph step
  reduces edit distance while the exact-gene count remains 57/128.
- Removed the unsupported statistical-equivalence claim for graph-consensus
  versus T1K. Restored the separate T1K run (96.5%, 4,308 calls) from
  `hla/results/tables/gc_truth_concordance_overall.tsv`. Conditional concordances
  use different resolved-call subsets. Removed the causal attribution of
  DRB1 disagreement exclusively to phasing: annotations of the same sequence
  change with the database.
- Corrected the conflation of 106/106 coarse DRB content calls with RCCX calls.
  Added the calibrated 88-donor RCCX results and the calibrated C4Investigator
  comparator from `hla-structural/REPORT.md` and `hla-targeted/REPORT.md`.
- Replaced confirmed-novel-allele language with evidence-ranked candidates from
  `hla-analysis/FINAL_REPORT.md`. Removed the old candidate figure from the
  manuscript because its embedded labels assert novelty and include an APR
  sequence rejected by the later complete-CDS audit. The original image remains
  in the repository. Local candidate read support is retained in prose.
- Corrected the extraction narrative: CPC and K-PanRef are graph-derived;
  JaSaPaGe Saudi is assembly-derived. Removed the unsupported tripling claim.
- Replaced the SpecHLA software-only citation with the published paper, added
  Locityper, and checked DOI metadata for existing references. Corrected
  inaccurate author lists and titles where publisher-deposited metadata resolved.

## Items requiring author input or further evidence

1. Author names, affiliations, ORCIDs and contributions are intentionally empty
   at Robert’s request. Complete these fields before submission.
2. The exact CPC graph release/accession used for the 116 retained paths needs
   to be linked to its originating release. The supplied citation resolves to
   the 2026 *1000 Chinese Pangenome* paper; the original draft called the input
   Phase 1. The revision describes the observed input without asserting a phase.
3. The graph-consensus snapshot should receive a fixed analysis date and exact
   eligible-call denominators behind the approximate 90%/33% resolution rates.
   The archived slides provide these percentages, while the checked TSV gives
   the resolved-call counts. Existing pairwise plots and aggregate TSVs pool
   different gene sets; the revision removes the disputed pooled percentages.
4. Exact historical versions for the initial consensus pipeline and the T1K
   database in the 40-donor comparison should be recovered from run manifests
   before submission. The current text preserves the documented versions and
   explicitly identifies native database differences.
5. The public repository exists, but unauthenticated checks returned HTTP 404 for
   `hla-asian50/refined/REPORT.md`, `hla-spechla-pg/RESULTS-2026-09-18.md` and the
   DōgoHLA validation protocol. Deposit/version these artifacts before submission;
   the manuscript states this availability gap. No push was performed.
6. Replace template publication metadata and settle the preprint licence before
   submission. The repository README describes the current template as CC0.

## Editorial approach

Added an abstract and organised the paper around reference construction,
variant recovery, HLA naming and sequence reconstruction. Used direct explanatory
prose, sentence-case headings and positive statements of what each endpoint
measures. Removed promotional wording, imagined-reader disclaimers, em dashes,
and claims of mechanism or general performance beyond the controls.

## Verification

- All 19 cited bibliography keys resolve; DOI metadata was retrieved for all
  DOI-bearing entries, including the added Locityper and SpecHLA publications.
- Metadata, image paths and figure/table labels pass structural checks.
- The BioHackrXiv generator completed its Pandoc, LuaLaTeX and Biber build.
  The regenerated PDF has 14 pages and no unresolved reference markers.
- Visually checked the title/abstract page and the new named-HLA table page.
- `git diff --check` passes. Manuscript and bibliography are edited locally;
  `paper/paper.pdf` has been regenerated. Nothing was committed or pushed.

## Author-requested follow-up

Authors, affiliations and the author footer are now empty. The title is
“Evaluating Asian-specific pangenomes for HLA typing”. The anonymous draft PDF
uses a temporary local copy of the generator with its two-author minimum
disabled and empty-author layout supported. The repository's upstream CI
generator still requires authors; fill in authorship before using that build.
