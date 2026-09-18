# BioHackrXiv Publication Template — DBCLS BioHackathon 2026

Template for a [BioHackrXiv](https://biohackrxiv.org/) publication reporting work done at the
[DBCLS BioHackathon 2026](https://2026.biohackathon.org/) (BH26JP), held 13–19 September 2026
in Matsuyama, Japan. The event metadata is already filled in — see
[For BH26 participants](#for-bh26-participants) below for what you still need to change.

Project pages: [Projects](https://github.com/dbcls/bh26/wiki/Projects) ·
[Schedule](https://github.com/dbcls/bh26/wiki/Schedule) ·
[Participants](https://github.com/dbcls/bh26/wiki/Participants) ·
[Wiki](https://github.com/dbcls/bh26/wiki/)

## For BH26 participants

**Ask the organizers for a repository.** Repository creation under the
[`biohackathon-japan`](https://github.com/biohackathon-japan) organization is handled by the
BioHackathon organizers. Give them your project name and your GitHub account, and you will get a
repository named `BH26-<your-project-name>` created from this template, with you as an
administrator. Please do not use the green "Use this template" button to create a repository
elsewhere — reports collected under one organization are easier to find and to archive.

Once you have your repository, three things need your attention:

1. **Leave the event metadata alone.** BH26 is registered in the BioHackrXiv index as
   [`BH26JP`](https://index.biohackrxiv.org/tag/BH26JP). The `event`, `biohackathon_name`,
   `biohackathon_url` and `biohackathon_location` fields in `paper/paper.md` are already correct;
   changing them will detach your report from the event.
2. **Point `git_url` at your own repository.** It still points at this template.
3. **Replace the license.** This template is CC0. Change `LICENSE` to the license of your preprint
   so that you can submit it to BioHackrXiv as CC-BY.

## Step 1: Configuring the Markdown

The publication Markdown is found in the `paper/paper.md` file. At the top you can edit the
YAML code with metadata. It is important to get this part correct, because otherwise the PDF
generation will fail. The metadata looks like this:

```yaml
title: 'DBCLS BioHackathon 2026 report: Template for the very long title'
title_short: 'BioHackJP26: How we found breakfast'
tags:
  - Semantic web
  - Ontologies
  - Workflows
authors:
  - name: First Author
    affiliation: 1
    role: Writing – original draft
  - name: Last Author
    orcid: 0000-0000-0000-0000
    affiliation: 2
    role: Conceptualization, Writing – review & editing
affiliations:
  - name: First Affiliation
    index: 1
  - name: ELIXIR Europe
    ror: 044rwnt51
    index: 2
date: 18 September 2026
cito-bibliography: paper.bib
event: BH26JP
biohackathon_name: "DBCLS BioHackathon 2026"
biohackathon_url:   "https://2026.biohackathon.org/"
biohackathon_location: "Matsuyama, Japan, 2026"
group: YOUR-PROJECT-NAME-GOES-HERE
# URL to project git repo --- should contain the actual paper.md:
git_url: https://github.com/biohackathon-japan/bh26-bhxiv-template
# This is the short authors description that is used at the
# bottom of the generated paper (typically the first two authors):
authors_short: First Author \emph{et al.}
```

### Which metadata to update?

#### To change

The following fields should be changed:

* title
* title_short
* tags
* authors (name, affiliation, and optionally their ORCID identifier and CRediT role)
* affiliations
* date
* group — your project name
* authors_short

Particularly important to update is the following field, which should point to
your own repository, instead of the template:

* git_url: https://github.com/biohackathon-japan/bh26-bhxiv-template

See [paper/paper.md](paper/paper.md) itself for the details on ORCID identifiers, ROR
identifiers for affiliations, and CRediT contributor roles.

#### Not to change

These fields describe the event and are already registered with BioHackrXiv. Leave them as they are:

* event: BH26JP
* biohackathon_name: "DBCLS BioHackathon 2026"
* biohackathon_url:   "https://2026.biohackathon.org/"
* biohackathon_location: "Matsuyama, Japan, 2026"

## Step 2: Writing the article

A full Markdown example is given in [paper/paper.md](paper/paper.md). This includes instructions how to include
figures, tables, and annotate citations with the Citation Typing Ontology.

## Step 3: Previewing the paper as PDF

This repository builds the PDF for you. The `Generate PDF` GitHub Action:

* checks that the required metadata fields are present in `paper/paper.md`;
* on a push to `main`, builds `paper/paper.pdf` and commits it back to the repository;
* on a pull request, builds the PDF, uploads it as a downloadable artifact, and comments on the
  pull request with a link to it.

Only one build runs at a time per branch: if you push again while a build is still going, the
older build is superseded, so the committed PDF always matches the newest `paper.md`.

So the current PDF is always at `paper/paper.pdf`, and you can review changes as a PDF before
merging them. Alternatively, the BioHackrXiv [Preview Server](http://preview.biohackrxiv.org/)
will build a PDF from a repository URL.

## Troubleshooting

### The first page is badly formatted

Sometimes the list of authors plus affiliations runs over the page. We are working on a fix, but in the mean time you can try to shorten the affiliations. If that does not work move the affiliations into a repo and put the affiliations on a web page and use something like

```yaml
affiliations:
  - name: For remaining affiliations see \url{https://github.com/project/etc} \vspace{0.2in}
    index: \*
```
