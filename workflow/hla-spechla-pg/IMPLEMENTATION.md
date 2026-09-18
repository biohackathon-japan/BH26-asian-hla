# AsianHLA–SpecHLA implementation — 2026-09-18

The eight-donor development candidate improves exact whole-gene reconstruction
from 43/128 to 57/128 and two-field typing from 61/63 to 62/63. Guarded graph
indels remove a further 1,290 sequence edits without losing correct names.
This remains a development result, not independent validation.

DōgoHLA is a population-specific, pangenome-based method using the family-excluded
AsianHLA panel for read collection and graph inference. The implementation uses
SpecHLA components for alignment, small-variant calling and initial reconstruction. Results and limitations are in [RESULTS-2026-09-18.md](RESULTS-2026-09-18.md).

## Implemented

1. **Repair phase-block linkage.** `phase_linkage.py` retains the maximum shared
   allele score and the best BLAST HSP per subject. `experiment.py` additionally
   repairs the complementary-pair combination: upstream groups by `x == y`,
   mixing parallel and crossed phase alternatives; the correct grouping is
   `x == 0`. Original, score-repaired, and fully repaired scripts are separate
   copies. All comparisons share a fixed eigensolver initialisation.
2. **Isolate database effects.** Native IPD, IPD 3.65 alone, and IPD 3.65 plus
   training-panel gene bodies are compared on identical BAM, VCF, masked regions
   and read bins. Panel flanks are removed before phase-block matching. Fold
   donors and relatives are excluded, and the donor/fold assignment is checked.
   The native final naming database is held fixed.
3. **Retain graph evidence.** `graph_diagnostic.py` verifies every graph path
   against its input sequence, checks family exclusions, retains Giraffe GAM,
   makes projected BAM, and calls represented alleles with `vg pack`/`vg call`.
   It records off-reference-node support, timing, mapped reads, and a matched
   linear BWA control. Fragment length is estimated from confident read pairs
   in the existing A-arm BAM; a small locus-specific subset otherwise produced
   an invalid zero-length fragment estimate. No assembly truth is used by mapping.
4. **Reconstruct confident structural alleles.** `structural_overlay.py` accepts
   explicit PASS homozygous graph indels of at least 50 bp with GQ ≥20, depth ≥10,
   alternate depth ≥5 and QUAL ≥20. Affine alignment decomposes complex alleles;
   only long indels outside annotated coding exons are eligible, preserving
   read-phased SNPs and small indels. It verifies REF and gene coordinates, maps
   boundaries to both reconstructed haplotypes, requires contiguous, unique
   30-base flanks (≤3 mismatches, no Ns), rejects overlaps, and preserves all
   sequence outside the replacement. It retains the original sequence when
   confidence or coordinate checks fail. The naming tools run again with the same databases; an isolated DRB1-only
   patch removes its fixed 11-kb query crop. Naming-only controls measure this
   change on unchanged sequences.
   Heterozygous structural variants are deliberately left unresolved; this is
   a conservative prototype, not a general structural phaser.
5. **Score complete paired cohorts.** `score_experiment.py` verifies completion
   manifests and output hashes for every predeclared donor/arm before scoring.
   Whole-gene global edit distance and exactness count N as an error; legacy
   infix distance is kept as a secondary metric. Two-field, G-group and
   three-field naming use the existing catalogue truth and nomenclature scorer.

Existing environments and legacy runs are preserved. New output lives beneath
`/home/leechuck/hla/spechla-pg/experiments/` on DDBJ. Failed or stale runs are
never silently reused. The runner checks direct command failures, expected
outputs and known nested SpecHLA failure messages, because upstream scripts
sometimes ignore subprocess exit status.

## Reproduce the development comparisons

Prerequisites: the existing DDBJ SpecHLA environment, source data and completed
fold-0 A-arm runs. Local scoring requires `edlib` and this repository's downloaded
nomenclature and assembly inputs. These commands use the repository's two-hop
DDBJ transport and submit compute through Slurm.

```sh
python3 -m unittest discover -s hla-spechla-pg -p 'test_*.py' -v
python3 hla-spechla-pg/deploy_experiment.py NEW_RELEASE --folds 0
```

The deploy command prints the preparation job ID and refuses an existing release.
After preparation succeeds, submit `code/jobs/phase_experiment.sbatch` with
`--array=0-7`, `RELEASE` set to the absolute release path, and `ARMS` set to the
space-separated arms to compare. The fully repaired comparisons are
`A-pairing AE-ipd365-pairing AE-panel-pairing`. The four original/partial-repair
controls are `A-replay A-fixed AE-ipd365 AE-panel`.

The graph/structural workflow is `jobs/hybrid_experiment.sbatch`. Submit it with:

```sh
python3 hla-spechla-pg/deploy_experiment.py NEW_HYBRID_RELEASE \
  --hybrid-phase-release COMPLETED_PHASE_RELEASE
```

The phase release must contain completed `A-pairing` and `AE-ipd365-pairing`
runs. The deploy command submits `--array=0-7`. Each job validates/maps all
DRB1 read pairs, calls graph alleles, and creates matched `AG-homozygous` and
`AEG-homozygous` sequence reconstructions plus `AF-full-DRB1` and
`AEF-full-DRB1` naming-only controls. Structural decomposition uses Biopython. Different reads or folds between the
phasing and graph runs are rejected. More than one million pairs triggers the
subset gate rather than silently treating a partial dataset as complete.

```sh
python3 hla-spechla-pg/fetch_experiment.py RELEASE
python3 hla-spechla-pg/score_experiment.py \
  --runs hla-spechla-pg/results/experiments/RELEASE/runs \
  --donors HG00706,HG02155,HG03239,HG03742,HG03804,NA18620,NA19087,NA21093 \
  --arms A-pairing,AE-ipd365-pairing,AE-panel-pairing \
  --out hla-spechla-pg/results/experiments/RELEASE/paired-score
```

Use the appropriate arm list for each release; missing/failed arms stop scoring.
Fetched native sequences and logs are ignored by Git. Aggregate tables, cohort
lists, reference/code hashes and graph diagnostics are kept for review. Large
sequence databases, reads and alignment files remain on DDBJ.

## Still needed before calling this a validated typer

- Freeze a candidate and test the remaining 32 development donors, then the
  locked validation cohorts; the current eight donors are development data.
- Measure graph-aware heterozygous structural phasing and residual novel-sequence
  discovery. `vg call` genotypes represented alleles; it does not discover new SVs.
- Compare original local realignment, long-indel-enabled native SpecHLA, and an
  adaptive linear reference on matched inputs. The conservative overlay is not
  a substitute for these controls.
- Expand locus coverage beyond the current DRB1 structural intervention and
  measure locus-assignment/paralog errors and graph-absent alleles.
- Revisit naming and confidence after sequence inference is stable. Population
  enrichment must not become an ancestry-based allele whitelist.
