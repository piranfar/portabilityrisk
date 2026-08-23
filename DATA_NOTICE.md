# Data notice — attribution for the derived observation table

Short companion to `DATA_LICENSE.md`. If you use the derived observation table,
this is what you owe and to whom.

## The one-line version

Cite this repository **and** the source study for every value you use. The DOI
is in the row.

## Provenance of the derived table

`audit/contracts/assay_aware_emergence/v1/ingestion/observations_v1.tsv`
contains 192 observations extracted from **9 peer-reviewed studies**. Each row
names its own origin; the table below lists the studies, the number of rows
taken from each, and the document the values were read from.

| study (PMID) | rows | source document referenced by the manifest |
|---|---:|---|
| 29431605 | 80 | `29431605_mgen000158.pdf` |
| 34826621 | 27 | `34826621_elsevier_fulltext.xml` |
| 36845973 | 26 | `36845973_PMC9948630.xml` |
| 32439569 | 23 | `32439569_1-s2.0-S2213716520301302-main.pdf` |
| 34653681 | 14 | `34653681_elsevier_fulltext.xml` |
| 38946899 | 9 | `38946899_PMC11211256.xml` |
| 30265709 | 8 | `30265709_pone.0204936.pdf` |
| 32021334 | 4 | `32021334_240404.docx` |
| 40239923 | 1 | `40239923_PMC12198106.xml` |

**None of those documents is stored in this repository.** The filenames appear
in `source_manifest.tsv` alongside a SHA-256 so that a reader who obtains the
document from its publisher can confirm they are looking at the same bytes the
extraction was made from. Resolve the `doi` field of any row to reach the
article.

The 8 rows from 30265709 are classified `held`, not `anchor` — their BioSample
link is absent from the source record. They are published as held, with the
reason recorded, rather than repaired by inference.

## What attribution requires

**When you use one or more values:**

- cite this repository (`CITATION.cff`); and
- cite the source study for each value, by the DOI carried in its row.

**Wherever rows go, the provenance goes with them.** `doi`, `source_sha256` and
`source_locator` must stay attached to the values. Those three columns are what
make a value checkable; a number separated from them is an assertion, not a
measurement. Note that this repository grants no permission to redistribute the
table — see `DATA_LICENSE.md` — so this is a condition on any onward use you have
an independent basis for, not an invitation to one.

**When you publish a figure or table built from this data,** say which studies
it draws on. "Derived from 9 published studies, listed in DATA_NOTICE.md" plus
the repository citation is sufficient for a methods section; a supplementary
table should carry the per-row DOIs.

## What attribution does not do

Citing this repository does not give you any right in the source articles, does
not make the maintainers a licensor of them, and does not by itself create a
permission to redistribute or relicense the derived table. Attribution is a
condition, not a grant. See `DATA_LICENSE.md`.

## Genome accessions

Every row carries `biosample_accession` and `assembly_accession`, resolvable at
NCBI. Those are public identifiers pointing at public records. No sequence data
is stored here, and the depositing authors are credited through the study
citation.

## Corrections and takedown

Extraction errors are defects — open an issue with the DOI and the cell.

If you are a rights holder and believe any part of this repository reproduces
more of your work than verification requires, open an issue or use the private
reporting route in `SECURITY.md`. The material will be shortened, replaced with
a citation, or removed.
