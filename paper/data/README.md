# Figure source data

These six TSV files are unchanged copies of the completed 17 September 2026
panel comparison under `hla-asian50/refined/results/` in the versioned analysis
archive linked from the project release. `provenance.json` gives original
relative paths and SHA256 hashes.

From the repository root, run `python3 scripts/plot_paper.py` with Python and
Matplotlib (tested with the installed version recorded in `paper/REVIEW.md`).
This produces Figures 2–3 and supplementary Figure S1 as PNG previews and vector
PDFs in `paper/figures/`. Figure 1 has separate inputs and instructions in
[`graph/README.md`](graph/README.md).
The manuscript uses the PDFs. No inference or bootstrap calculation is rerun.

Source method identifiers remain unchanged in these data; the plotting script
maps them to publication labels. Public donor accessions are retained for
reproducibility. Variant rows cover multiple comparison universes; the main
figure uses `all_truth/truth_SV_length` and `frozen_threeway_SNV/SNV`, with the
`all` endpoint and repaired `full` and `hprc` arms.
