# DōgoHLA 0.1.0: prospective 32-donor extension

Frozen before new inference on 18 September 2026. Deadline: **09:02:07 UTC**
(12:02:07 Riyadh), five hours after the request. Cohort and method digests are
in `frozen-method.json`; `cohort.tsv` lists every planned donor.

## Design

Evaluate all 32 remaining fold 1–4 donors: 16 East Asian and 16 South Asian,
32 distinct families, none overlapping the eight-donor method-development
pilot. These are an extension of the existing development cohort, **not** the
untouched 228/946 locked benchmark. Use existing complete 1000 Genomes MHC
read recruitments, identical input FASTQs for each paired comparison.
Exclude the entire held-out fold and related family haplotypes from panel
references and graphs. Verify graph path sequences and exclusions before use.

Run native SpecHLA 1.0.12 with its default small-variant workflow, DōgoHLA
0.1.0, and DōgoHLA without graph reconstruction (same repaired phasing,
IPD 3.65 phase database and DRB1 naming correction). The last comparison
isolates the structural graph contribution. Test native SpecHLA's optional
long-indel mode on a pilot donor first; if executable within the time budget,
run it as an additional comparator. Report availability/failures explicitly.

## Outcomes and analysis

Primary: exact whole-gene haplotypes against assembly-derived truth across
8 genes × 2 haplotypes × 32 donors (512 planned sequences), using global
edit distance and minimum-distance unordered diploid assignment. Report
truth availability explicitly. Secondary: total global edits, masked Ns,
call/completion rates, assembly-derived two-field/G-group calls, and the
available independent experimental Gourraud genotypes. Infix scores are
secondary and must never be mixed with global scores.

Report every donor, gene and ancestry stratum, not only improvements. Missing
outputs and no-calls receive no exact/correct credit in planned denominators;
edit-distance comparisons require paired valid outputs and are labelled
conditional on completion. Missing truth is separately accounted for and
excluded identically for every method. Report paired donor differences and
95% percentile bootstrap intervals (10,000 resamples, seed 20260918), sampling
whole donors/families rather than treating correlated alleles as independent.
Show the graph ablation even if it is neutral or worse.

## Execution and stopping rules

Freeze inference code and thresholds before new results. Isolate outputs and
record code, reference and input hashes, Slurm IDs, timing, exit status and
validated completion markers. Do not silently reuse partial runs. Execution
repairs must be logged; changing inference after viewing results creates a
new exploratory version and cannot retain the frozen primary label.

Use Slurm on DDBJ, four CPUs per donor, and reuse existing read/reference
assets. Build the four missing normalised DRB1 graphs from cached raw graphs.
At the deadline, report all completed, failed and pending donors. Do not claim
32-donor validation if fewer completed. No outcome-dependent cohort extension
or early stopping. Native long-indel feasibility is an execution gate, not an
accuracy-based selection. This short validation can support further study;
it cannot establish general superiority or clinical validity.

## Launch record

Frozen method/documentation: `3f791b8`; isolated execution wrapper: `9a6c7d7`.
DDBJ release: `/home/leechuck/hla/spechla-pg/validations/dogohla-v0.1.0-20260918`.
Slurm jobs: preparation **20663742**, graph array **20663743_[1-4]**, paired
validation **20663744_[0-31]**, native long-indel feasibility **20663769**.
Preparation succeeded. The long-indel feasibility job uses pilot HG02155,
not a new validation donor. Its only upstream execution patch replaces the
invalid random TCP-port generator with an OS-selected valid free port;
its exact command and patched-script hash are recorded on DDBJ.

Truth-only audit before inference: all **256 diploid gene pairs / 512
haplotypes** have assembly truth; all 256 have eligible two-field truth.
Independent experimental truth is available for 29 gene pairs.

A local monitor fetches compact snapshots every three minutes and automatically
scores the fixed cohort on completion or at the deadline. Reproduce with:

```sh
python3 hla-spechla-pg/monitor_validation.py dogohla-v0.1.0-20260918 \
  --out hla-spechla-pg/results/validation-20260918
```

Read `analysis/REPORT.md`, `status.tsv`, `paired_donors.tsv` and `bootstrap.json`
inside that output directory. Bootstrap units are donors (each from a distinct
family); positive gains favour DōgoHLA. Exact/name gains retain failures in the
planned denominator. Edit reductions are explicitly conditional on both arms
completing. These exploratory intervals do not adjust for multiple endpoints.

The monitor runs as the local user service `dogohla-validation-20260918.service`
so it survives the interactive session. A separate deadline publisher runs
`publish_validation.py` after the terminal snapshot: it pushes only compact
analysis tables and the timestamped monitor record, and refuses to commit if
the branch changed or unrelated work is staged. No raw reads or alignments
are published. The generated report explicitly labels incomplete runs.

## Long-indel feasibility gate passed

Pilot HG02155 completed native `-v True` in approximately 11 minutes, with
valid sequence/naming outputs. ScanIndel reports three fully processed genes
and five explicit no-soft-clipped/unmapped-read exits; these latter exits use
code 1 upstream but are its expected no-candidate path, not hidden tool errors.
Consequently the prespecified supplemental **native-long** arm is included for
all 32 donors. It runs the complete native workflow on identical FASTQs.

Its isolated wrapper fixes two installation issues: invalid TCP ports and an
absent `db/ref/HLA_G.fa` phase-reference path (redirected to the installed
`db/HLA/HLA_G/HLA_G.fa`). It also makes ScanIndel's hardcoded eight BWA threads
respect the four allocated CPUs. Allele databases, thresholds, phasing and
variant inference remain native. The wrapper rejects nested tool errors and
requires a recognised terminal outcome for each of the eight ScanIndel stages.
This is labelled a compatibility-corrected native long-indel control, separately
from untouched native default SpecHLA. DōgoHLA's frozen inference is unchanged.

## User-requested supplemental T1K and four-field analysis

Added after viewing the first six outcomes, at the user's request; this is
supplemental and does not replace or alter the frozen primary endpoints.
`score_comparators.py` compares matched donors with T1K 1.0.6 (`hla-wgs`, IPD
3.65) archived in `hla-typer/results/comparators/t1k`. All 32 local genotype
files match completed DDBJ outputs by SHA256; see `t1k-provenance.json`.
The T1K job script uses the same shared recruited-read directory. Its legacy
completion markers lack per-run FASTQ hashes, so exact historical input
identity is supported by the workflow, not independently hash-certified.

Four-field scoring uses **exact_genomic_alleles**, not CDS-compatible labels.
Both truth haplotypes must have exact genomic matches with four numeric fields.
Missing/novel genomic labels and shorter truth labels are excluded equally for
all methods; eligible denominators are reported. Predictions with fewer than
four fields are unresolved, never padded or expanded into correct calls.
Comma/semicolon ambiguity must be compatible for every listed alternative.
The endpoint is unordered diploid genotype correctness; expression suffixes
are not additional numeric fields. Whole-gene accuracy remains separate.

Reports retain the original first-six cohort, the growing matched completed
cohort and the fixed planned 32-donor cohort, including per-gene scores.
T1K's archived calls do not provide four-field resolution on the eligible
completed cases; this is reported as unresolved, not evidence that DōgoHLA
has better sequence reconstruction. Native/DōgoHLA naming uses IPD 3.38,
whereas T1K uses 3.65, so this compares deployed configurations rather than
isolating algorithm effects. No superiority over other established short-read
callers is claimed without matched runs.

The persistent monitor now also writes `analysis/COMPARATORS.md`,
`comparator_summary.tsv`, `comparator_gene_scores.tsv` and
`comparator_metadata.json`; these are included in the deadline publication.
