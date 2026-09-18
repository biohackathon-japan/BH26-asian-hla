# AsianPGR Figure 1 data

Figure 1 uses the final construction panel: 754 haplotype entries, comprising
376 sample entries across source projects and two reference sequences. Five
cross-project duplicate donors give 371 distinct donors. Cohort counts are
computed from `hla-typer/source/haplotype_donors.tsv` and stored in `cohorts.json`.
The map locates the five added Asian/Arab cohorts; the worldwide HPRC contribution
is shown separately. Markers indicate cohort geography, with counts by source.

## Rebuild the figure

From the repository root:

```sh
python3 scripts/plot_graph_figure.py
```

Requires Python, NumPy and Matplotlib. The script reads frozen topology and layout
coordinates, and writes `paper/figures/fig1_asian_graph.pdf` and its PNG preview.
These files include no internal run labels or local filesystem paths.

## Regenerate bundles locally

The bundles were regenerated on the local workstation from final-panel sequences
with pgr-tk v0.6.0 (build 102bdf5). The input FASTAs are classII.fa (DRA to DMA,
20-kb flanks), HLA-A.fa and HLA-DPB1.fa (2-kb flanks). Their byte counts and hashes,
exact decomposition options and source descriptions are in `provenance.json`.
The regional FASTAs derive from the same extraction/annotation workflow described
in the manuscript; large input FASTAs are kept outside this manuscript directory.

Place those three FASTAs in `inputs/`, obtain the construction structural backbone
`MHC.sv.gfa.gz` and construction donor manifest, then run locally:

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/prepare_graph_figure.py \
  --input-dir inputs \
  --bundle-dir bundles \
  --structural-gfa inputs/MHC.sv.gfa.gz \
  --donors inputs/haplotype_donors.tsv
python3 scripts/plot_graph_figure.py
```

This requires pgr-pbundle-decomp, NetworkX, SciPy, NumPy and Graphviz sfdp.
The preparer limits pgr-tk to four Rayon threads. `--reuse-decompositions` packages
already regenerated local BED outputs; omit it to run the decompositions.
The graph and regional FASTAs were retrieved/read as inputs; all bundle generation
and layout computation were performed locally.

## Graph representations

- `MHC.structural.json`: the actual Minigraph structural backbone from construction,
  1,677 sequence segments and 2,381 links. Segment metadata and connectivity are
  retained without sequence strings. GRCh38 segments have rank 0. The final
  Minigraph-Cactus graph has 417,896 nodes and 579,137 edges.
- `classII.json`, `HLA-A.json`, `HLA-DPB1.json`: bundle graphs reconstructed from
  locally regenerated BED files. An undirected edge joins distinct consecutive
  bundle IDs observed along at least one sequence. Node colour uses the number
  of distinct carrying sequences. Multiple occurrences within a sequence count
  once. Bundle orientations are omitted from these connectivity views.
- `*.bed.gz`: regenerated decomposition intervals, with transient command paths
  removed from the header. The exact options are retained in `provenance.json`.
- `reference_anchors.json`: GRCh38 gene intervals in the extracted MHC coordinate
  frame and the class II extraction interval. Labels anchor to the nearest
  reference segment/bundle midpoint; graph layouts have no genomic distance scale.

The class II / HLA-A / HLA-DPB1 inputs contain 753 / 753 / 754 sequences and yield
890 / 42 / 23 bundles. The first two input sets each lack one of the construction
haplotypes. Counts include reference sequences. Large graphs use unweighted
Kamada–Kawai layouts; small graphs use Graphviz sfdp with start=42. Coordinates
are rotated to their principal axes and scaled isotropically before plotting.
The independently chosen panel aspect ratios are for display.

## Geographic data

Country polygons are frozen in `natural_earth_countries.geojson`, from Natural
Earth's 1:110m countries dataset. The source URL, retrieval date and SHA256 are
in `provenance.json`. Natural Earth data are [public domain](https://www.naturalearthdata.com/about/terms-of-use/).
The figure uses the dataset for geographic context and depicts cohort locations
schematically. Gene/reference coordinates and cohort counts come from project
manifests, not from the map dataset.
