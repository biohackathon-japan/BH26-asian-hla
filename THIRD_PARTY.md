# Sources, attribution and third-party terms

The repository's existing CC0 dedication covers original project material to
which its authors can apply it. It does not replace third-party licences or
source-data terms. Release archives preserve the notices below.

## Software

DōgoHLA uses components from **SpecHLA v1.0.12**, by the DeepOmics Lab, under
the MIT licence. The release includes its notice at
`workflow/hla-spechla-pg/LICENSE-SpecHLA`; the inspected upstream copies and
phase-linkage test fixture retain this attribution. Cite Wang et al.,
[SpecHLA enables full-resolution HLA typing from sequencing data](https://doi.org/10.1016/j.crmeth.2023.100589).

External tools are installed separately under their own licences, including
[Cactus/Minigraph-Cactus](https://doi.org/10.1038/s41587-023-01793-w),
[vg/Giraffe](https://doi.org/10.1126/science.abg8871),
[pggb](https://doi.org/10.1038/s41592-024-02430-3),
[odgi](https://doi.org/10.1093/bioinformatics/btac308),
[pgr-tk](https://doi.org/10.1038/s41592-023-01914-y),
[Immuannot](https://doi.org/10.1101/gr.278985.124),
[PanGenie](https://doi.org/10.1038/s41588-022-01043-w) and
[Locityper](https://doi.org/10.1038/s41588-025-02362-4).
The manuscript bibliography supplies the remaining methodological references.

## HLA nomenclature and databases

IPD-IMGT/HLA v3.65.0 is pinned to commit
`5b915f27f7f620361cf83cb626eeac8a03c0247c` in the
[official repository](https://github.com/ANHIG/IMGTHLA).
`scripts/fetch_nomenclature.py` retrieves unmodified nomenclature files and the
provider's licence directly. Raw IPD sequence databases and nomenclature files
are excluded from the release archives; project-derived annotation labels and
benchmark tables identify their database version. The provider's licence and
citation instructions apply to downloaded files. The method's external
SpecHLA/IPD databases must be obtained and prepared separately.

## AsianPGR sequence sources

| Source | Included haplotypes | Attribution and provenance |
| --- | ---: | --- |
| Arab Pangenome Reference | 106 | Nassir et al., [Nature Communications 2025](https://doi.org/10.1038/s41467-025-61645-w); [APR dataset](https://zenodo.org/records/13752609), CC BY 4.0 |
| HPRC release 2 | 464 | [Human Pangenome Reference Consortium](https://humanpangenome.org/); source identifiers are retained in the contig mapping |
| JaSaPaGe | 38 | JaSaPaGe project, 18 Saudi and 20 Japanese haplotypes supplied for the BioHackathon; assembly identifiers are retained in the contig mapping |
| K-PanRef | 28 | Shin et al., [2026 preprint](https://doi.org/10.64898/2026.07.06.26357367); [source graph deposit](https://zenodo.org/records/20810335) |
| Chinese Pangenome Consortium | 116 | Wang et al., [The 1000 Chinese Pangenome empowers medical and population genetics](https://doi.org/10.1038/s41586-026-10315-y); graph-derived paths, with source names retained |
| GRCh38 and CHM13 | 2 | Genome Reference Consortium and Telomere-to-Telomere reference paths |

AsianPGR is a derived MHC graph. Source datasets retain their original terms;
this release assigns no additional licence to third-party sequences. Attribution,
source identifiers, clipping flags and build commands accompany the graph.
Raw individual sequencing reads and original whole-genome assemblies are not
included. The release preserves the graph constructed for the reported study.
