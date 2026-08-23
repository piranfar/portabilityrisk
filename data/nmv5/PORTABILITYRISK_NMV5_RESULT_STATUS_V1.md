# NM-V5 — result status

**Version 1.0.0, 2026-08-23.** The publication-facing resource layer for Paper 1.

## Verdict

**NM-V5 RESOURCE COMPLETE.** The occurrence-level dataset is packaged for deposit with a full data
dictionary, denominator and evidence-level dictionaries, a provenance map, pinned software and
database versions, JSON Schemas, checksums, a manifest and a standalone verification script.

## What was produced

| | |
|---|---|
| deposited dataset | `portabilityrisk_occurrence_portability_v1.tsv` — 74,349 rows, 19 columns |
| canonical source | `determinant_portability_classes.tsv`, 74,349 x 20, retained **privately and unchanged** |
| difference | one column dropped: `plasmidcall_predicted_location`, empty in all 74,349 rows |
| proof | every retained cell identical, no row removed, reordered or renamed, every denominator unchanged |
| canonical SHA-256 | `ce0c99d61ba8435895135b2cfa6c8a48aa41a03c07e1540c87092d3c1332617b` |
| derivative SHA-256 | `957e529ae56c79593f85c45fe95cb8592fc111ffadb6dbfc1b22857faa88d564` |

The derivative digest differs from the canonical **by design**. It must never be quoted as the
canonical digest.

## Reconciliation, re-derived from the deposited file

74,349 occurrences; chromosome 35,140; plasmid 39,209; A 18,837; B 16,303; C 7,170; D 6,043;
E 25,996; A+B = 35,140; C+D+E = 39,209. Cohort: 6,288 genomes, 6,285 BioSamples, 2,283 BioProjects,
109 species, 12,811 replicons.

## Terminology corrected

**36.586% is `3,569 / 9,755` events after collapsing to at most one event per genome per
compartment.** It is **not** the arithmetic mean of per-genome percentages, which is
**35.932096%**. Neither number changed; the description of the first one did. The verification
script computes both and asserts they differ, so the confusion cannot recur silently.

## Not deposited

The occurrence-level dataset is **not** committed to GitHub. Restricted material — unblinding keys,
blinded casebooks, adjudication instruments, credentials, recovery artefacts — is deposited
nowhere. Raw sequence, caches and environments are excluded as regenerable.

## Licence

Dataset and metadata **CC-BY-NC-ND-4.0**. Code **Apache-2.0**. Third-party data and software retain
their own terms; no licence is applied to them by this deposit.
