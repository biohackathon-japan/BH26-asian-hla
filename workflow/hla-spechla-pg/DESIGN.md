# SpecHLA-PG: SpecHLA's workflow with the Asian MHC pangenome substituted where it matters

Target: a short-read (Illumina 2x150, 30x WGS) HLA typer that keeps the architecture of SpecHLA
(Wang et al. 2023, Cell Reports Methods 3:100589, Figure 1 steps A-F) and replaces its IPD-IMGT/HLA-only
references by our 754-haplotype Asian-enriched MHC pangenome at steps A, D, E and F. Steps B and C are
SpecHLA's own scripts, run unchanged. Everything below refers to SpecHLA 1.0.12 as installed in
`/home/leechuck/hla/mm/envs/asian50-spechla` (`share/spechla/script/`, `share/spechla/db/`); the
IMGT/HLA release in its database is 3.38 (paper) / `hla_nom_g.txt` header says 3.51.0 for the G-group file.

## 1. SpecHLA as implemented (Figure 1 A-F)

| Step | What SpecHLA does | Reference used | Source |
|---|---|---|---|
| A. Extract HLA reads | `ExtractHLAread.sh` pulls MHC-region, HLA-contig and unmapped reads from an hg38 BAM. `SpecHLA.sh` then aligns all of them with `bowtie2 --very-sensitive -k 30` (novoalign if licensed) to the "lite" IMGT/HLA database. | `db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta`: 6,172 genomic alleles of 39 HLA genes/pseudogenes (IMGT 3.38), each extended with flanking sequence; records named `A*01:01:01:01` etc. | `script/ExtractHLAread.sh`, `script/whole/SpecHLA.sh` lines "assign the reads to original gene" |
| B. Reads binning | `assign_reads_to_genes.py`: per read pair and per allele, both mates must align to the same allele with no soft clip and NM <= 2 each; score = matched bases / mapped length; per gene take the best allele; assign to the top gene if its score >= 0.1 and exceeds the second gene by >= 0.1 (`-d 0.1`), else discard. Writes `<gene>.R{1,2}.fq.gz` for A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1; gene = allele name before `*`. | same bowtie2 BAM as A | `script/assign_reads_to_genes.py` |
| C. Alignment + local assembly | Per gene `bwa mem -U 10000 -L 10000,10000` (end-to-end) of the binned reads to the IMGT representative reference (first recorded allele per gene + 1 kb flanks: HLA_A 5,503 bp ... HLA_DRB1 13,229 bp), merged into `<sample>.merge.bam`. `run.assembly.realign.sh` then assembles the reads of 25 hypervariable windows (`whole/select.region.txt`) with `fermi2.pl unitig`, blasts the unitigs to the gene reference, realigns reads through the unitigs (`rematchblast.pl`, `realignblast.py`) and writes `<sample>.realign.sort.bam`. | `db/ref/hla.ref.extend.fa` (8 contigs `HLA_<gene>`), `db/HLA/HLA_<gene>/` (bwa + blast indexes) | `SpecHLA.sh`, `script/run.assembly.realign.sh` |
| D. Variant calling | `freebayes -a -p 3` on the realigned BAM; keep variants inside `HLA_<gene>:1000-(len-1000)`; `mask_low_depth_region.py` masks 20 bp windows with mean depth < 5 (`-k 5`). Long indels (> 150 bp) via ScanIndel only with `-v True` (off by default; off in all our runs). Inside `phase_variants.py::read_vcf` triploid calls are reduced to diploid, alleles with MAF < 0.05 (`-r`) are dropped, and variants whose two top alleles sum to < 0.7 of depth are discarded. | `hla.ref.extend.fa` | `SpecHLA.sh`, `script/phase_variants.py` |
| E. Phasing | `phase_variants.py`: ExtractHAIRs-style fragments from the BAM (`MNP_linkage`), `SpecHap --ncs --protocols ngs,matrix --window_size 15000` (spectral graph bipartition). Phase blocks that reads cannot link are linked through the database: `whole/map_block2_database.py` writes each block's two haplotypes with `bcftools consensus`, blasts them (`-subject db/HLA/whole/HLA_<gene>.fasta`), scores every 00/11 vs 01/10 block pairing by identity x length of the best shared allele, and `phase_unlinked_block.py` bipartitions the block graph by its Fiedler vector. `Share_reads.split_seg` then writes `hla.allele.{1,2}.HLA_<gene>.fasta` (whole reference contig incl. 1 kb flanks, consensus of the rephased VCF, N over masked windows; the DRB1 3,898-4,400 repeat is re-typed by blast against `ref/DRB1_dup_extract_ref.fasta`). | `db/HLA/whole/HLA_<gene>.fasta` (IMGT genomic alleles: A 2,474, B 3,106, C 2,924, DPA1 86, DPB1 441, DQA1 140, DQB1 239, DRB1 797) | `script/phase_variants.py`, `script/whole/map_block2_database.py`, `script/phase_unlinked_block.py` |
| F. Designation | `annoHLA.pl -r whole`: each haplotype's "diversity region" (e.g. `HLA_A_0:100-3300`, DQB1 two windows, DRB1 100-11000) is blasted against `db/HLA/whole/HLA_<gene>` (DRB1: `HLA_DRB1.exon`), score = 100 x (1 - (mismatches+gaps)/aligned length) summed over HSPs; alleles with zero frequency in `HLA_FREQ_HLA_I_II.txt` are skipped unless `-p nonuse`; ties -> highest population frequency; two hard-coded SNP fixes (DRB1*14:01/14:54, C*07:01/07:18). `g_group_annotation.py` blasts the haplotypes against `hla_exons.fasta` for a G-group call. Output `hla.result.txt`, `hla.result.g.group.txt`. | `db/HLA/whole/`, `db/HLA/HLA_FREQ_HLA_I_II.txt`, `db/HLA/hla_nom_g.txt`, `db/HLA/hla_exons.fasta` | `script/whole/annoHLA.pl`, `script/whole/g_group_annotation.py` |

## 2. What the pangenome replaces

The unit of the pangenome used here is the per-gene sequence set: for 11 loci (A, B, C, DRB1, DRB3, DRB4,
DRB5, DQA1, DQB1, DPA1, DPB1) the gene +-2 kb of every panel haplotype (`hla-analysis/source/HLA-<gene>.fa`,
records `sample#hap#HLA-<gene>`, labelled against IPD-IMGT/HLA 3.65.0 in
`hla-analysis/results/sequence_catalogue.tsv`). SpecHLA's coordinate system (`HLA_<gene>` = representative
allele + 1 kb flanks) is kept everywhere, so every SpecHLA script downstream of the alignment runs unchanged.

| Step | Substitution | Built by |
|---|---|---|
| A | Read-extraction/binning database = SpecHLA lite IMGT db (6,172 alleles, kept as decoys for the 31 non-target genes and pseudogenes) + the fold's training panel gene sequences of the 11 loci, renamed `<GENE>*pg\|<sample>#<hap>` so B's `name.split('*')[0]` still yields the gene. bowtie2 index. Reads from alleles absent from IMGT 3.38 (intron variants, novel Asian alleles: 31.6% of panel gene sequences have no exact IPD genomic match) can now satisfy B's NM <= 2 / no-soft-clip rule. | `build_refs.py` -> `db/fold{k}/ref/hla_gen...v2.fasta` |
| B | unchanged (`assign_reads_to_genes.py`, `-nm 2 -d 0.1`) | - |
| C | unchanged in mode A. In graph modes the linear `bwa mem` step is replaced by `vg giraffe` (see D) but the fermikit local assembly + realignment (`run.assembly.realign.sh`) still runs on the surjected BAM. | - |
| D | Alignment for variant calling goes through a per-fold, per-gene pangenome graph whose reference path is SpecHLA's `HLA_<gene>` contig (`SpecHLA#0#HLA_<gene>`) and whose other paths are the fold's training panel sequences. `vg giraffe -o BAM --ref-name SpecHLA` surjects onto that path, so the BAM is in SpecHLA coordinates and `freebayes -p 3`, masking and everything after run unchanged. Reads of divergent alleles are placed through their own haplotype path instead of being forced (`-L 10000,10000`) onto the representative allele. Graph builder: PGGB by default; Minigraph-Cactus available for the comparison (section 4). | `build_refs.py` (inputs), `jobs/build_graph_pggb.sbatch`, `jobs/build_graph_mc.sbatch` |
| E | Block-linking database `db/fold{k}/HLA/whole/HLA_<gene>.fasta` = IPD-IMGT/HLA 3.65.0 genomic alleles + training panel sequences. `map_block2_database.py` is unchanged; it now links unlinked blocks through Asian haplotypes that IMGT lacks. | `build_refs.py` |
| F | `designate.py`: each reconstructed haplotype is compared by exact edit distance (edlib infix mode) with every full-length gene body of the same database (IMGT 3.65 genomic alleles, partial entries < 90% of median length dropped, + training panel). Ranking is two-stage, as in SpecHLA's own class II designation which blasts exon windows only: first the summed edit distance over the database record's exons (`source/exons.tsv`, built by `build_exons.py` and validated against all 5,027 panel records with 0 `cds_sha256` mismatches; it annotates 97-99% of IPD genomic alleles per gene), then, among exon ties, the whole gene body. The exon stage also removes the length bias of infix alignment, which in the first run made the shortest surviving DRB1 record (`DRB1*13:370N`) win every DRB1 haplotype. Remaining ties -> labelled over novel, then the pangenome frequency prior (number of training haplotypes carrying the same two-field label; replaces SpecHLA's population frequency table), then name. A panel winner is reported with its catalogue label (exact genomic, else exact CDS allele). 31.6% of panel gene sequences have no exact IPD match at any level and so carry no label; when such a record wins, the nomenclature call falls back to the closest *labelled* record under the same ranking (`label_source=nearest:<record>`, extra distance in `nearest_*_edits`) while the novel record stays the sequence-level winner. Without the fallback every donor whose closest panel haplotype is an unlabelled Asian allele is a no-call: that alone cost DRB1 0/7 in the first run. The native `annoHLA.pl` + G-group scripts are also run on the same sequences, so F-native vs F-pangenome is a pure post-hoc comparison. | `designate.py` |

Population frequency database (Figure 1 top): SpecHLA's `HLA_FREQ_HLA_I_II.txt` (Caucasian/Black/Asian
two-field frequencies) is replaced in F by counts over the fold's training panel haplotypes.
Pedigree / long-read inputs (Figure 1 centre) are not used.

## 3. Configurations and ablations

`spechla_pg.sh -m MODE` (all share the A+B output of a donor in `runs/<donor>/bin/`):

| Mode | A | B | C | D alignment | E db | F |
|---|---|---|---|---|---|---|
| native (`spechla`, baseline) | IMGT lite | SpecHLA | SpecHLA | bwa | IMGT | annoHLA |
| A | pangenome | SpecHLA | SpecHLA | bwa | IMGT | annoHLA and designate.py |
| AD | pangenome | SpecHLA | SpecHLA (assembly on surjected BAM) | giraffe on fold graph | IMGT | annoHLA and designate.py |
| ADEF | pangenome | SpecHLA | SpecHLA | giraffe on fold graph | IMGT 3.65 + panel | annoHLA and designate.py (the ADEF call is designate.py) |

The A-only, A+D and A+D+E+F ablations asked for are mode A / F-native, mode AD / F-native and mode ADEF /
F-pangenome. Because F runs both designators on every mode, F's own effect (A vs A+F, AD vs AD+F) is also
measured for free.

## 4. Graph builder: PGGB vs Minigraph-Cactus

Per gene and fold the graph input is ~600-750 sequences of 5-15 kb. Both builders are run on fold 0 for HLA-A
(easy) and HLA-DRB1 (hardest, most divergent) with identical inputs; compared on build time, node/edge
counts, the fraction of binned reads that giraffe aligns with identity 1.0, and the downstream typing of the
fold-0 dev donors (mode AD, both graphs). Result and choice: see REPORT.md section "Graph builder".
Rationale before the test: PGGB is reference-free and does not clip haplotype-specific sequence, which
matters for the DRB1 intron 1 / 3' region where alleles differ by kilobase-scale indels; Minigraph-Cactus
anchors on the reference path and, at gene scale, its minigraph SV backbone offers no advantage.

The global 754-haplotype whole-MHC Cactus graph (`/home/asianhla/data/upload/HLA/mhc_graph/`) is NOT used
anywhere in this pipeline: it contains the test assemblies and would leak.

## 5. Leakage handling

For fold k (folds from `hla-asian50/run/source/donors.tsv`, 8 dev donors each) every haplotype listed for
fold k in `excluded_paths.tsv` (the donors' own two haplotypes plus pedigree relatives, 16-18 haplotypes per
fold) is removed from the binning database, the graph, the block-linking database, the designation database
and the frequency prior. `build_refs.py` additionally removes any haplotype whose `validation_groups.tsv`
group contains a fold donor and records what that added (`db/build_log.json`, field `extra_from_groups`).
A donor is always typed with `db/fold{k}` and `graphs/fold{k}` of its own fold. IPD-IMGT/HLA itself is a
public database and is not held out.

## 6. Reads

Recruitment is the main session's logic (chr6:28,410,700-33,580,488 + all chr6 alts + all 525 `HLA-*`
contigs with complete pairs + all unmapped reads), run by `jobs/recruit.sbatch` into `reads/<donor>/`
(re-using `/home/leechuck/hla/hla-typer/reads/<donor>/` when its COMPLETE sentinel exists). The defective
`hla-asian50/run/source/mhc_regions.bed` is not used.

## 7. Evaluation

Truth A: the donor's two catalogue haplotypes per gene (`exact_cds_alleles` -> two-field / three-field /
G-group; `gene_sha256` and the gene body for sequence identity). Truth B: Gourraud 2014 experimental typings
for the 8 dev donors it covers (A, B, C, DRB1, DQB1). Scoring with `hla-bench/names.py` (`compare`, levels
two_field and g_group; unordered pair; a no-call is an error). Sequence-level: edit distance of the truth
gene body inside the reconstructed haplotype (edlib infix), exact if 0. See `score.py`.
