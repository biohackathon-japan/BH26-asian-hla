# AsianPGR v0.1.0

AsianPGR is the whole-MHC pangenome graph built from 754 haplotype entries:
106 APR, 464 HPRC release 2, 18 JaSaPaGe Saudi, 20 JaSaPaGe Japanese,
28 K-PanRef, 116 CPC, GRCh38 and CHM13. The reference is enriched for Asian
and Arab haplotypes. It covers the extended MHC, using the GRCh38 interval
chr6:28,510,120–33,480,577 plus 100 kb flanks for extraction.

[Download the graph archive](https://bio2vec.net/data/asianpgr/v0.1.0/AsianPGR-v0.1.0.tar.gz).
[Release notes](https://github.com/biohackathon-japan/BH26-asian-hla/releases/tag/v0.1.0)
and source code are maintained in the BioHackathon repository.
Extract `AsianPGR-v0.1.0.tar.gz` and run `sha256sum -c SHA256SUMS` inside
its directory. The original `MHC.*` filenames preserve compatibility with the
build and mapping scripts.

## Files

| File | Purpose |
| --- | --- |
| `MHC.gbz` | Clipped graph and haplotype paths for Giraffe mapping |
| `MHC.full.gbz` | Full graph and haplotype paths |
| `MHC.dist`, `MHC.shortread.withzip.min`, `MHC.shortread.zipcodes` | Matching short-read mapping indexes |
| `MHC.snarls`, `MHC.full.snarls` | Variant-site decomposition indexes |
| `MHC.vcf.gz`, `MHC.vcf.gz.tbi` | Graph-derived variants relative to GRCh38 |
| `contig_names.tsv` | Source haplotype-to-renamed-contig mapping |
| `seqfile.txt`, `mhc_mc_graph.sbatch` | Original input manifest and build command |
| `MHC.input-contig-sizes.tsv.gz`, `MHC.stats.txt` | Build input sizes and graph statistics |
| `BUILD_README.md` | Original build and mapping notes |
| `SOURCE.json`, `SHA256SUMS` | Source identity and file checksums |

The clipped graph has 417,896 nodes, 579,137 edges and 5,885,602 sequence bases.
Cactus v3.3.0 built the graph on 15 September 2026 using the Minigraph-Cactus
algorithm, with GRCh38 and CHM13 reference paths. Full GFA, HAL and odgi derivatives
are not included in this compact release; the full GBZ retains the full graph.

## Map paired reads

With vg v1.76.1, from the extracted graph directory:

```sh
vg giraffe -t 8 -Z MHC.gbz -d MHC.dist \
  -m MHC.shortread.withzip.min -z MHC.shortread.zipcodes \
  -f reads_R1.fastq.gz -f reads_R2.fastq.gz -o gaf > reads.gaf
```

If vg rejects index timestamp ordering after extraction, preserve the files and
run `touch MHC.shortread.withzip.min MHC.shortread.zipcodes`. Reference paths use
`GRCh38#0#ctgN` and `CHM13#0#ctgN`; consult `contig_names.tsv` for source coordinates.
The graph is an MHC resource. Whole-genome read recruitment and competing paralogs
need to be handled by the chosen genotyping workflow.

## Provenance and interpretation

K-PanRef and CPC inputs were extracted from their source graphs and inherit
upstream clipping. The 754 inputs contain duplicate donor assembly sources;
donor-level analyses use 742 haplotypes from 371 reconciled donors after removing
the references and five duplicate pairs. The global graph includes assemblies
used as test donors in the paper. Family-excluded inference panels and the
DōgoHLA locus graphs are separate derived artifacts, described in the analysis
archive and method scripts.

Original `seqfile.txt` paths identify the build environment; reconstruct those
inputs from the source project data when rebuilding. Source projects and terms
are listed in [THIRD_PARTY.md](../THIRD_PARTY.md). Cite these projects, the
Minigraph-Cactus and Giraffe papers, and the BioHackathon manuscript when using
this resource.
