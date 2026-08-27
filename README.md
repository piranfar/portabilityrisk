# PortabilityRisk

**Replicon-resolved evidence that plasmid fraction underestimates the mobile resistome.**

Antimicrobial-resistance surveillance records whether a resistance gene is present. It does not
record whether that gene is positioned to move. This repository holds the manuscript, protocols,
results, data resource, validation record and code for a study that measures the second thing.

In 6,288 closed complete Gram-negative ESKAPE genomes, each of **74,349 acquired resistance-gene
occurrences** is assigned to a documented chromosome or plasmid using NCBI molecule designations
rather than prediction — **100 % resolution, 0 unmatched, 0 ambiguous, 0 missing coordinates** —
and receives one of five evidence-ranked portability classes.

Three results follow:

- **164 gene families occur in both compartments.** Portability is a property of the occurrence,
  not of the gene, and a gene-presence table has no field in which the difference could be written
  down.
- **46.39 % of chromosomal occurrences lie within 10 kb of a mobile-element marker** — in
  *Acinetobacter baumannii* at a median distance of **647 bp**. Restricting the endpoint to
  structurally complete insertion sequences reproduces the same host ordering across a
  21,955-block census and corroborates **73.80 %** of the homology calls.
- **Plasmid fraction and chromosomal mobile-element association are largely orthogonal.** Fitted
  on eight confirmation species, plasmid share explains **2.75 %** of the variance in chromosomal
  MGE association and the slope cannot be distinguished from zero.

**No transfer event is reported.** The claim is about genomic architecture, not about demonstrated
mobilization. The prohibited-claims list is in [`MANUSCRIPT.md`](manuscript/MANUSCRIPT.md) and is
binding on any revision.

---

## What is here

| directory | contents |
|---|---|
| [`manuscript/`](manuscript/) | the full draft, eight figures in SVG/PDF/300 dpi/600 dpi, and per-module interpretation documents |
| [`protocols/`](protocols/) | every frozen protocol and numbered amendment, the frozen class definitions, the pinned environment lock |
| [`results/`](results/) | result tables, summary tables the manuscript quotes, and a signed receipt per module |
| [`validation/`](validation/) | independent verification reports, result statuses, the validation ledger, the claim status matrix, the error and deviation register |
| [`data/`](data/) | the publication data resource: a **273-column data dictionary at 100 % curation**, table catalogue, denominator dictionary, evidence-level dictionary, provenance map, claim-to-evidence crosswalk, and a JSON Schema per table |
| [`code/`](code/) | analysis code: `validation/` for the NM modules and the figure generator, `context/` for the portability-context pipeline |

`MANIFEST.tsv` lists every file with its size and full 64-character SHA-256.

## Start here

- **The one number readers get wrong** — 36.586 % is `3,569 / 9,755` genome-collapsed events, not
  a mean of per-genome percentages (35.932 %), and not the occurrence-weighted share (52.736 %).
  [`data/PORTABILITYRISK_PUBLIC_DATA_README_V1.md`](data/PORTABILITYRISK_PUBLIC_DATA_README_V1.md).
- **What each evidence level does and does not establish** —
  [`data/PORTABILITYRISK_EVIDENCE_LEVEL_DICTIONARY_V1.tsv`](data/PORTABILITYRISK_EVIDENCE_LEVEL_DICTIONARY_V1.tsv).
- **Where any headline number comes from** —
  [`data/PORTABILITYRISK_CLAIM_TO_EVIDENCE_CROSSWALK_V1.tsv`](data/PORTABILITYRISK_CLAIM_TO_EVIDENCE_CROSSWALK_V1.tsv),
  which gives the table, the field and the digest for each.

## How the work is meant to be checked

Each module follows the same three-script shape, and the shape is the point:

- **freeze** — writes a protocol and hashes it *before* any outcome column is read;
- **score** — verifies those digests and refuses to run on a mismatch;
- **verify** — re-derives the published numbers independently and reports its own disagreement
  count.

The verifiers reported **0 disagreements**: 35 checks on the structural census, 43 on the
consolidated tables, 15 on the denominator flow. `portabilityrisk_figures.py` does the same for
figures — fifteen input digests verified before drawing, every output digest recorded.

## Failures are in the record

[`validation/PORTABILITYRISK_ERROR_AND_DEVIATION_REGISTER_V4.tsv`](validation/PORTABILITYRISK_ERROR_AND_DEVIATION_REGISTER_V4.tsv)
has 30 entries, including the ones that reached a reader before they were caught. The blinded
expert audit of the chromosomal MGE layer **failed at 62 of 120** on its first round; that round
is retained, because it diagnosed the instrument rather than the classifier, and the corrected
instrument then returned 0.9920. Deleting the failure would leave the pass looking unearned.

One sub-claim of the three-architecture analysis is registered **false**, and one host contrast is
registered **weak**. Both are labelled in the claim status matrix and in the manuscript.

## Redacted derivatives

Four files here are **redacted public derivatives**: infrastructure identifiers naming the compute
instance are replaced by `[REDACTED:…]` tokens. No scientific field, value, threshold or protocol
term was altered. Each carries an in-band marker and appears in
[`REDACTION_MANIFEST.json`](REDACTION_MANIFEST.json) with **both** its own SHA-256 and the
canonical private SHA-256.

The NMIS frozen-protocol digest quoted throughout the manuscript refers to the **private
canonical**. The public derivative hashes to something else, by design — silently editing the
canonical so the hash still appeared to match would be provenance theatre.

## What is deliberately not here

- **The upstream occurrence-level dataset is no longer withheld.** It is published in the
  Zenodo deposit under **CC BY 4.0** — https://doi.org/10.5281/zenodo.22116987 (concept DOI
  https://doi.org/10.5281/zenodo.22065541 resolves to the latest version) — together with its data
  dictionary, denominator dictionary and evidence-level dictionary. Earlier revisions of
  this README described it as `PRIVATE — pending deposit`; that was true then and is not
  true now.
- **Adjudication material.** Unblinding keys, casebooks and received adjudications would
  retrospectively unblind the expert audits. The **rubric is here**, because an audit whose scoring
  rules are unpublished cannot be checked; the keys are not.
- **Sequence caches and conda environments** — regenerable from accession receipts with per-file
  digests, and from the pinned lock in `protocols/environments/`.
- **Anything from the other papers in this programme.** Where another paper is named, it is to
  state a boundary.

**Executability, stated honestly:** scripts that consume the occurrence-level dataset can now be
run rather than only read, because that dataset is published. They resolve their input directory
from an environment variable and stop with a message if it is unset — see the `path_rewrites` block in [`REDACTION_MANIFEST.json`](REDACTION_MANIFEST.json).
[`SCOPE_AUDIT.txt`](SCOPE_AUDIT.txt) lists which is which.

## Licensing

Code is **Apache-2.0** ([`LICENSE`](LICENSE)). **The deposited dataset is CC BY 4.0** and may be
redistributed, built on and used commercially under attribution; see the Zenodo record. Tables
held in this repository that are derivatives of that deposit inherit CC BY 4.0. Anything in this
tree that is NOT part of the deposit remains under the terms in [`DATA_LICENSE.md`](DATA_LICENSE.md)
and [`DATA_NOTICE.md`](DATA_NOTICE.md). Underlying records are public NCBI data under their own terms.

## Status

The manuscript is a **draft**, not a submission, and has not been peer reviewed.

See [`AUDIT.md`](AUDIT.md) for the publication-safety audit, and
[`SCOPE_AUDIT.txt`](SCOPE_AUDIT.txt) for the per-file scope classification.
