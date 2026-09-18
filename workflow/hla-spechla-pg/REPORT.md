# SpecHLA-PG results on the development donors

All numbers below were produced by `score.py` from `results/runs/`, which is fetched verbatim from
`/home/leechuck/hla/spechla-pg/runs/`. Every table states its own denominator. Scoring uses
`hla-bench/names.py` (`Nomenclature`, `truth_slot`, `pred_slot`, `compare`), imported unmodified.

## 1. What is being compared

| Arm | A (read extraction db) | D (alignment for calling) | E (block-linking db) | F (designation) |
|---|---|---|---|---|
| `native/F-native` | SpecHLA IMGT lite (6,172 alleles, IMGT 3.38) | bwa mem | SpecHLA IMGT | `annoHLA.pl` |
| `native/F-native-G` | as above | bwa mem | SpecHLA IMGT | `g_group_annotation.py` |
| `native/F-pg` | as above | bwa mem | SpecHLA IMGT | `designate.py` (pangenome) |
| `A/F-native` | **pangenome** (IMGT lite + fold training panel) | bwa mem | SpecHLA IMGT | `annoHLA.pl` |
| `A/F-pg` | **pangenome** | bwa mem | SpecHLA IMGT | **pangenome** `designate.py` |
| `AD/*`, `ADEF/*` | **pangenome** | **giraffe on the fold's per-gene graph** | ADEF: **pangenome** | see section 5 |

Steps B (binning rule) and C (fermikit local assembly + realignment) are SpecHLA's own scripts in every
arm. `native` is native SpecHLA 1.0.12 run exactly as asked
(`spechla -n <donor> -1 r1 -2 r2 -o <dir> -j 4 -u 0 -p Unknown`) on the same recruited reads, so it is a
clean comparator arm; the `hla-asian50` SpecHLA outputs are not reused.

Donors: the 8 fold-0 development donors (HG00706, HG02155, HG03239, HG03742, HG03804, NA18620, NA19087,
NA21093). Leave-one-out: for fold 0 every reference built from the panel excludes those donors' own
haplotypes and their `validation_groups.tsv` relatives (`source/excluded_paths.tsv`,
`db/build_log.json`). IPD-IMGT/HLA itself is public and is not held out.

Truth A = the donor's own two assembly haplotypes in `sequence_catalogue.tsv` (IPD 3.65.0 exact CDS
alleles; genes with no exact CDS allele on either haplotype are dropped, which is why the denominator is
63 and not 64 genotypes). Truth B = Gourraud 2014 experimental typings, 10 genotypes across A, B, C,
DQB1 and DRB1 for the donors it covers.

## 2. Genotype accuracy, Truth A, 8 donors x 8 genes

Unordered pair, no-call counted as an error.

| Arm | two-field | three-field (CDS) | G group |
|---|---|---|---|
| `native/F-native` | **61/63 = 96.8%** | 55/62 = 88.7% | 58/63 = 92.1% |
| `native/F-native-G` | 56/63 = 88.9% | 52/62 = 83.9% | **62/63 = 98.4%** |
| `native/F-pg` | 55/63 = 87.3% | 54/62 = 87.1% | 58/63 = 92.1% |
| `A/F-native` | **61/63 = 96.8%** | 55/62 = 88.7% | 58/63 = 92.1% |
| `A/F-native-G` | 56/63 = 88.9% | 52/62 = 83.9% | **62/63 = 98.4%** |
| `A/F-pg` | 58/63 = 92.1% | **56/62 = 90.3%** | 58/63 = 92.1% |

Truth B (Gourraud 2014, 10 genotypes): `native/F-native` 10/10 two-field and 10/10 G group;
`A/F-pg` 9/10 two-field (a DRB1*15:01 case, see below), 10/10 G group.

Per gene, two-field, Truth A (correct/n):

| Arm | A | B | C | DPA1 | DPB1 | DQA1 | DQB1 | DRB1 |
|---|---|---|---|---|---|---|---|---|
| `native/F-native` | 8/8 | 8/8 | 7/7 | 7/8 | 7/8 | 8/8 | 8/8 | 8/8 |
| `A/F-native` | 8/8 | 8/8 | 7/7 | 7/8 | 7/8 | 8/8 | 8/8 | 8/8 |
| `native/F-pg` | 8/8 | 7/8 | 7/7 | 8/8 | 6/8 | 7/8 | 8/8 | 4/8 |
| `A/F-pg` | 8/8 | 8/8 | 7/7 | 8/8 | 6/8 | 8/8 | 8/8 | 5/8 |

Reading: **replacing step A alone does not change the nomenclature call at all** (`A/F-native` is
identical to `native/F-native`, genotype for genotype). It does change the reconstructed sequences, and
there it helps (section 3). The pangenome designator (F) is 3 genotypes behind `annoHLA.pl` at
two-field, entirely on DRB1 and DPB1, but it is ahead at three-field, where `annoHLA.pl`'s
frequency-table tie-break sends it to the wrong sub-allele.

## 3. Sequence accuracy, Truth A, 8 donors x 8 genes x 2 haplotypes

Each reconstruction is compared with the donor's own assembled gene body (edlib infix, best of the two
pairings). "Exact" = edit distance 0. Masked bases (N) count as mismatches, as in the designation step.
n = 16 haplotypes per gene (8 donors x 2); pooled n = 128.

| Gene | native exact/16 | mode A exact/16 | native mean edits | mode A mean edits |
|---|---|---|---|---|
| A | 13 | **15** | 4.6 | **0.5** |
| B | 3 | **8** | 33.1 | **7.8** |
| C | 10 | **13** | 3.8 | **0.8** |
| DPA1 | 7 | 7 | 24.2 | **23.7** |
| DPB1 | 1 | 1 | 42.8 | **34.8** |
| DQA1 | 7 | 7 | 164.1 | **141.6** |
| DQB1 | 0 | 0 | 477.7 | **409.2** |
| DRB1 | 0 | 0 | 4272.9 | **3155.1** |
| **pooled** | **41/128 (32.0%)** | **51/128 (39.8%)** | **627.9** | **471.7** |

This is the clearest result so far: the pangenome read-extraction database at step A improves the
reconstructed sequence of every gene, and improves class I markedly (HLA-B goes from 3/16 to 8/16 exact,
mean edit distance 33.1 -> 7.8). The mechanism is the one the design predicted: 31.6% of panel gene
sequences have no exact IMGT 3.38 genomic match, and their reads fail SpecHLA's binning rule (both mates
aligned, no soft clip, NM <= 2) against the IMGT-only database. Adding the panel haplotypes as binning
targets recruits them.

The improvement does not reach the nomenclature call because `annoHLA.pl` already tolerates those
residual mismatches: a reconstruction 33 edits from truth still blasts to the right two-field allele.
The gain pays off for a caller that reports sequences, and at three-field/four-field resolution.

## 4. Failure analysis

Every wrong or missing call is in `results/errors.tsv` with the evidence that produced it (winner record,
`label_source`, per-exon edit distance, masked bases, sequence edit distances).

**`A/F-native` and `native/F-native`, 2 wrong genotypes of 63 (two-field):**

| Donor | Gene | Truth | Called | Cause |
|---|---|---|---|---|
| HG03742 | DPA1 | 01:03 / 02:02 | 02:02 / 02:02 | allele dropout: the 01:03 haplotype's reconstruction is far from truth and both haplotypes blast to 02:02 |
| NA19087 | DPB1 | 03:01 / 05:01 | 03:01 / 03:01 | same pattern; DPB1 exon 2 depth is the limit |

**`A/F-pg`, 5 wrong genotypes of 63 (two-field): 3 DRB1 + 2 DPB1.**

| Donor | Gene | Truth | Called | Evidence |
|---|---|---|---|---|
| HG00706 | DRB1 | 15:01 / 15:01 | 11:34 / 11:34 | sequence winner is the unlabelled panel haplotype `RY06#1#HLA-DRB1` at 29 exon edits; nearest labelled record is DRB1*11:34 |
| HG02155 | DRB1 | 15:01 / 15:01 | 11:34 / 11:34 | same |
| HG03804 | DRB1 | 15:01 / 15:04 | 13:24 / 11:34 | same |
| HG00706 | DPB1 | 02:02 / 14:01 | 1493:01 / 651:01 | exon edits 2 and 0, 178 masked bases; the winners are rare IPD alleles exon-identical to the N-corrupted reconstruction |
| NA19087 | DPB1 | 03:01 / 05:01 | 135:01 / 104:01 | exon edits 0 and 0, no masking; pure exon degeneracy |

Two distinct failure modes, neither a scoring artefact:

1. **DRB1*15:01 reconstructions are wrong, not the designator.** All three DRB1 errors are donors whose
   truth contains DRB1*15:01. SpecHLA's DRB1 reference contig is a DR1-group allele; DR2 haplotypes
   differ from it by kilobase-scale intron 1 and 3' indels, so `bwa mem -U 10000 -L 10000,10000`
   end-to-end alignment forces DR2 reads onto the wrong coordinates. The mean DRB1 reconstruction is
   3,155 edits from truth even in mode A. `annoHLA.pl` survives this because for DRB1 it blasts exon
   windows only and scores identity per HSP; `designate.py` sums exact edit distance and is therefore
   honest about a reconstruction that is genuinely broken. This is precisely what step D exists to fix
   and is the strongest argument for the AD arm.
2. **DPB1 exon degeneracy.** IPD contains hundreds of DPB1 alleles whose exons are identical or nearly
   identical and which differ only in introns. When the reconstruction carries a few errors or masked
   bases, an arbitrary rare allele (DPB1*1493:01, DPB1*135:01) reaches exon distance 0 before the common
   one does. The panel frequency prior does not rescue this because the rare winners are IPD records,
   not panel records, and the tie-break only fires on exact ties. Fix: break ties on the gene body
   before the frequency prior, or drop IPD records whose two-field label has zero panel support.

Two designation bugs were found and fixed during this work; both are generic and worth recording.

* **Infix length bias.** Ranking by `edlib` HW distance of a whole database gene body inside the
  reconstruction systematically favours the *shortest* database record, since HW charges nothing for the
  unmatched tail of the target. In the first run this made `DRB1*13:370N` (a short null allele) the
  winner for every DRB1 haplotype of every donor: DRB1 scored 0/7. The exon-first ranking stage removes
  the bias because exon lengths are near-constant across records.
* **Unlabelled panel winners.** 31.6% of panel gene sequences have no exact IPD match at any level. When
  one of them wins, reporting `novel:<record>` is a no-call, i.e. an error. `designate.py` now falls back
  to the closest *labelled* record for the nomenclature call (`label_source=nearest:<record>`) while
  keeping the novel record as the sequence-level winner. Before the fix DRB1 scored 0/7.

## 5. Step D: graph alignment

Status: **not yet measured.** The per-gene PGGB graphs for fold 0 were built for all 8 genes, and
`vg giraffe` maps and surjects 7 of them in about 2 minutes per donor in total. The DRB1 graph is the
exception and it broke the first AD run: on the raw seqwish graph (6,688 nodes, 9,167 edges, 38,887 bp of
sequence for a 13 kb locus) giraffe's gapless-extension step took up to 461 seconds and 650 MB of memory
*for a single read pair* (`runs/HG00706/AD/DRB1.giraffe.log`).

Cause: the first build ran wfmash + seqwish only and skipped PGGB's normalisation stage, so every
redundant parallel walk through the DRB1 intron repeat survived in the graph. `build_graph_pggb.sh` now
runs the full PGGB sequence (smoothxg -> gfaffix -> odgi build/unchop/sort) before indexing, which is
exactly the stage that exists to make such graphs mappable. The fold-0 DRB1 graph is being rebuilt with
it; AD and ADEF results for the 8 fold-0 donors belong in this section once it completes.

Graph complexity of the raw seqwish fold-0 graphs, for orientation:

| Gene | nodes | edges | graph bp | reference contig bp |
|---|---|---|---|---|
| A | 1,542 | 2,088 | 8,568 | 4,503 |
| B | 1,761 | 2,366 | 9,695 | 5,081 |
| C | 1,576 | 2,121 | 9,581 | 5,304 |
| DPA1 | 1,621 | 2,169 | 16,533 | 10,775 |
| DPB1 | 1,660 | 2,237 | 17,215 | 12,468 |
| DQA1 | 2,793 | 3,796 | 16,379 | 7,492 |
| DQB1 | 3,434 | 4,668 | 21,174 | 8,480 |
| DRB1 | 6,688 | 9,167 | 38,887 | 12,229 |

### Graph builder: PGGB vs Minigraph-Cactus

PGGB is the builder in use. The head-to-head test (fold 0, HLA-A and HLA-DRB1, identical inputs) was
started twice but `cactus-pangenome` on 738 haplotypes did not finish inside the compute available on the
shared node and was cancelled to keep the critical path moving, so **the comparison is not completed**
and the choice currently rests on weaker evidence:

* PGGB is reference-free; Minigraph-Cactus anchors on the reference path and clips sequence that does not
  align to the minigraph backbone. At this locus the divergent material (DRB1 intron 1, the DR haplotype
  groups) is exactly what must be kept, since it is where the reads that motivate the whole exercise
  come from.
* Minigraph-Cactus' advantage is an SV backbone built by minigraph, which at 5-15 kb gene scale has
  nothing to contribute.
* The one measurement obtained is a negative one about *parameters* rather than builders: an
  unnormalised graph is unusable for giraffe regardless of who built it. Minigraph-Cactus applies its own
  normalisation and would not have shown this failure; that is a point in its favour on robustness, not
  on representation.

This is the weakest part of the study and should be closed before the method is written up.

The global 754-haplotype whole-MHC Cactus graph (`/home/asianhla/data/upload/HLA/mhc_graph/`) is not used
anywhere in this pipeline: it contains the development donors' own assemblies and would leak.

## 6. Runtime and memory per donor (4 threads, `asianhla-c32`)

| Stage | native SpecHLA | SpecHLA-PG mode A |
|---|---|---|
| A+B read extraction and binning (computed once, shared by all PG modes) | included below | 555-768 s, median 608 s |
| pipeline after binning | 495-657 s, median 538 s | 91-1001 s, median 138 s |
| peak RSS | 0.8 GB | 0.3 GB (NA21093: 1.7 GB) |

End to end the two arms are comparable, roughly 700-1,700 s per donor. Binning dominates the PG cost
because `bowtie2 --very-sensitive -k 30` runs against a database about 10% larger. Of the recruited
reads, 7.8-8.6% align to the pangenome-augmented allele database, and 6,490-9,653 read pairs survive
SpecHLA's binning rule per donor.

## 7. Assessment: can this route beat T1K and SpecHLA?

**At two-field resolution on these 8 genes, no, and not by a small margin.** T1K on the identical reads
scores 98.4% over the 40 development donors. Native SpecHLA here scores 96.8% (61/63) on 8 donors, the
pangenome-at-A arm scores exactly the same, and the pangenome designator is behind at 92.1%. No
configuration measured so far overtakes T1K, and the honest reading is that two-field typing of the
classical genes is already saturated by a well-built allele-matching method.

**Where this route is genuinely better is sequence reconstruction**, which is a different product. Mode A
reconstructs 39.8% of gene haplotypes exactly against 32.0% for native SpecHLA and cuts the mean edit
distance by a quarter (628 -> 472), with the largest gains at HLA-B (33.1 -> 7.8 edits) and HLA-A
(4.6 -> 0.5). T1K does not produce sequences at all. If the goal is a novel-allele-capable caller for an
Asian cohort, that is the axis to develop, and the metric is three-field/four-field accuracy and exact
CDS recovery, where `A/F-pg` is already the best arm measured (56/62 = 90.3% at three-field, against
55/62 for `annoHLA.pl` and 52/62 for the G-group script).

**Three things must happen before any stronger claim.** (i) The AD arm has to run: every DRB1 error here
traces to a reconstruction that linear alignment to a DR1-group reference cannot produce, and step D is
the only part of the design that addresses it. (ii) DRB3/4/5 are not in this pipeline at all and are
exactly where T1K collapses (48/110); SpecHLA's architecture has no copy-number step, so adding them is
new work, not a configuration change. (iii) The denominator is 8 donors and 63 genotypes; a 3-genotype
difference is one donor. The 32 remaining development donors are staged and queued.

## 8. Reproducing

See `README.md`. `results/` is the exact output of `./fetch_runs.sh && python3 score.py` against the
cluster state described above.
