# Data licence and rights boundary

**The Apache License 2.0 in `LICENSE` covers code and schemas only. It does not
license data, and it grants no rights whatsoever over the source publications
this repository cites.**

This file exists because a single root `LICENSE` on a repository that ships both
software and extracted measurements is routinely read as licensing everything in
the tree. It does not, and reading it that way would be wrong in a way that
could mislead a downstream user into republishing material the maintainers have
no right to license.

## What is covered by Apache-2.0

- `audit/ingest/assay_aware_emergence/v1/` — the builder, the tests, the
  assay-attribution registry.
- `audit/contracts/assay_aware_emergence/v1/DATA_CONTRACT.md`, the JSON Schemas,
  `controlled_vocabularies.json`, `raw_data_manifest.template.tsv`.
- `tools/`, and the repository documentation.

## What is NOT covered by Apache-2.0

### 1. The source publications

The observations in this repository were read out of peer-reviewed articles.
Those articles, their tables, their figures, their supplements and their
full-text XML remain the property of their authors and publishers under whatever
terms those parties set.

**No right in any source publication is granted, sublicensed, transferred or
implied by this repository or by its Apache-2.0 LICENSE.** If you want a source
document, obtain it from the publisher under your own access rights.

### 2. Third-party source data

Any reference to an external database, accession, or resource is a citation.
Their terms are theirs. Nothing here relicenses them.

### 3. The derived observation table

`audit/contracts/assay_aware_emergence/v1/ingestion/observations_v1.tsv` and its
companion ledger, manifest and receipt are a **curated derived dataset**. Their
status is set out below.

## The derived observation table

### What it is

192 antimicrobial susceptibility observations extracted from 9 published
studies, normalised against the data contract, with the provenance of each value
retained.

**Every row carries all three of:**

| field | what it is |
|---|---|
| `doi` (with `study_group`) | the citation for the study the value came from |
| `source_sha256` | the SHA-256 of the exact source document the value was read from |
| `source_locator` | the table and cell — e.g. `Table 1, row Colistin, column "Isolate 724942", printed "1"` |

Those three fields exist so that any value can be checked against its origin by
a third party who obtains the source independently. They are not decoration:
wherever a row goes, they go with it.

### What is deliberately absent

**Source documents are not redistributed.** No publisher PDF, DOCX, supplement,
or full-text XML is stored in, served by, or obtainable from this repository.
`source_manifest.tsv` references each source by filename, digest and locator so
you can fetch it yourself and verify the digest matches. That is the entire
mechanism: **cite and hash, never mirror.**

### Its licence status: deliberately unassigned

**No open data licence is applied to the derived table. Not CC0, not CC-BY, not
ODbL, not any other.**

This is a decision, not an oversight, and it should not be resolved by anyone
other than the repository owner:

- Measurements are facts, and facts are generally not subject to copyright in
  the United States and are thin at best elsewhere. On that view a licence is
  unnecessary.
- But a *curated selection and arrangement* of facts can attract database or
  sui generis rights in some jurisdictions, and the extraction here is
  substantively curated — studies screened, drugs attributed per instrument,
  rows admitted or held against a contract.
- And the underlying articles are third-party works. Stamping CC0 or CC-BY on a
  table derived from them risks asserting a freedom over material the
  maintainers do not own, which is worse than saying nothing.

**Do not add CC0, CC-BY, or any other data licence to this table without the
explicit approval of the repository owner.** A pull request that does so will be
declined. If you need a specific licence for a specific use, open an issue and
ask.

### What is and is not permitted in the meantime

The table is published **for inspection, verification, and reproducibility
assessment. No licence or permission for redistribution, relicensing, or
commercial reuse is granted by this repository.**

That is a statement about what this repository grants, which is nothing. It is
not an assertion that the underlying measurements are owned — facts generally
are not — and it is not a prohibition dressed up as one. Whether any particular
use is lawful for you depends on your jurisdiction and your purpose, and this
repository is silent on that question rather than resolving it in your favour.
If you need a permission, ask for it: open an issue describing the use.

Whatever your legal basis for using the data, these obligations stand:

1. **Cite this repository** — see `CITATION.cff`.
2. **Cite the underlying studies.** Required, not optional. A value in this
   table is a measurement made by the authors of the source article; they did
   the work. Every row carries the DOI needed to credit them, and quoting a
   value without citing its source study misattributes someone's result.
3. **Keep the provenance columns together with the values.** `doi`,
   `source_sha256` and `source_locator` must travel with any rows that leave
   this repository, under whatever basis they leave it. A value stripped of its
   provenance cannot be verified and should not be trusted, including by you.
4. **Do not present it as a survey.** It is a strict-inclusion literature
   extract, not a representative sample. See `MODEL_CARD.md`.
5. **Do not use it clinically.** See `MODEL_CARD.md`.

## Quotation in the assay-attribution registry

`assay_attribution_registry.json` contains short methods sentences quoted from
source articles, each attributed by DOI. They are there so a reader can verify
which instrument measured which drug — the single most error-prone step in this
field — and they are quoted rather than paraphrased because a paraphrase cannot
be checked against the original.

These are brief factual method statements used for verification and criticism,
with attribution. If you are a rights holder and consider any quotation longer
than necessary, open an issue and it will be shortened or reduced to a citation.

## Corrections

If you believe any row misreads its source, misattributes an assay, or
reproduces more of a source than is warranted, that is a defect. Open an issue
with the DOI and the cell. Corrections are treated with the same seriousness as
code bugs.

## Summary

| material | licence |
|---|---|
| code, schemas, contract text, docs | Apache-2.0 (`LICENSE`) |
| derived observation table and its receipts | **no licence assigned.** Published for inspection, verification and reproducibility assessment; no permission for redistribution, relicensing or commercial reuse is granted here. Cite this repository and the source studies; provenance columns must travel with the data |
| source publications and supplements | **not licensed here, not redistributed**; obtain from the publisher |
| third-party databases and accessions | governed by their own terms |
