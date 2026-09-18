# hla-spechla-pg: SpecHLA with pangenome references

An independent short-read HLA typer that keeps SpecHLA's architecture (Wang et al. 2023, Figure 1 steps A-F)
and substitutes our Asian MHC pangenome at steps A (read extraction database), D (graph alignment for variant
calling), E (block-linking allele database) and F (designation database + frequency prior). Steps B and C
are SpecHLA's own scripts. See `DESIGN.md` for the step-by-step mapping and `REPORT.md` for results.

Cluster root: `/home/leechuck/hla/spechla-pg` (NIG/DDBJ, partition `asianhla-c32`). Local = scripting, scoring.

## Layout

| Path | Purpose |
|---|---|
| `remote.py` | two-hop ssh transport (`ssh ddbj` -> `a001`); `--upload FILE --dest PATH`, `--download PATH --dest FILE`, or a command |
| `push.sh` | upload every `*.py`, `*.sh` and `jobs/*` to `scripts/` and `jobs/` on the cluster |
| `fetch_runs.sh` | download the small per-run result files into `results/runs/<donor>/<mode>/` |
| `build_refs.py` | (cluster) build the five fold-specific, leakage-free references: binning db + bowtie2 index, E/F allele db + blast db, graph input FASTAs |
| `jobs/build_refs.sbatch` | Slurm wrapper for `build_refs.py` |
| `jobs/build_graph_pggb.sbatch` | array job: per fold x gene PGGB graph + giraffe indexes (`graphs/fold{k}/HLA_<gene>.pggb.*`) |
| `jobs/build_graph_mc.sbatch` | array job: the Minigraph-Cactus alternative (`*.mc.*`) used for the builder comparison |
| `jobs/recruit.sbatch` | array job over `source/donors_folds.tsv`: HLA read recruitment from the 1000G S3 CRAM (same logic as `hla-typer/stageA/recruit.sbatch`) |
| `spechla_pg.sh` | (cluster) the pipeline driver: `-m A|AD|ADEF`, `-g pggb|mc`, `-k fold` |
| `designate.py` | (cluster) step F with the pangenome database and frequency prior |
| `jobs/run_donor.sbatch` | array job: native SpecHLA baseline + modes A, AD, ADEF for one donor (`MODES`, `GRAPH` env overrides) |
| `jobs/redesignate.sbatch` | re-run step F only (`designate.py`) over every existing run directory, after a change to the designation logic or database |
| `score.py` | (local) predictions/scores/summary/errors against Truth A (assembly catalogue) and Truth B (Gourraud 2014) |
| `test_designate.py`, `test_score.py` | unit tests: `python3 -m unittest` |
| `source/` | fold definitions, donor list with CRAM URLs, IPD-IMGT/HLA 3.65.0 genomic FASTAs + nomenclature, slim catalogue |
| `results/` | `predictions.tsv`, `scores.tsv`, `summary.tsv`, `sequence_scores.tsv`, `errors.tsv` (one row per wrong or missing call), fetched run outputs |

## Running

```
./push.sh                                                    # upload scripts
python3 remote.py 'cd /home/leechuck/hla/spechla-pg && sbatch jobs/build_refs.sbatch'
python3 remote.py 'cd /home/leechuck/hla/spechla-pg && sbatch --array=0-39%2 jobs/build_graph_pggb.sbatch'
python3 remote.py 'cd /home/leechuck/hla/spechla-pg && sbatch --array=0-39%2 jobs/recruit.sbatch'     # index = data line of source/donors_folds.tsv
python3 remote.py 'cd /home/leechuck/hla/spechla-pg && sbatch --array=0-39%2 jobs/run_donor.sbatch'
python3 remote.py 'cd /home/leechuck/hla/spechla-pg && sbatch jobs/redesignate.sbatch'   # step F only, after a designate.py change
./fetch_runs.sh && python3 score.py                          # -> results/*.tsv
```

Cluster outputs: `reads/<donor>/`, `runs/<donor>/{bin,native,A,AD,ADEF}/` (`hla.result.txt` = native designation,
`hla.result.g.group.txt`, `hla.result.pg.txt` = pangenome designation, `hla.allele.{1,2}.HLA_<gene>.fasta`
reconstructed haplotypes, `timing.tsv`, `time.log`), logs in `logs/`.

Environments: `asian50-spechla` (SpecHLA 1.0.12, bowtie2, bwa, freebayes, fermi2, blast, SpecHap, edlib),
`pggb053` (pggb 0.5.3 with pinned wfmash 0.10, seqwish 0.7, smoothxg 0.7, odgi 0.8; the unpinned bioconda solve pairs pggb 0.5.3 with wfmash 0.24 whose CLI is incompatible), Cactus 3.3.0 bundle (vg 1.76.1, cactus-pangenome).
Package versions: `source/*.packages.json`.
