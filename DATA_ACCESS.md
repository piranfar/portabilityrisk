# Getting the occurrence-level dataset

**It is not in this repository, and that is deliberate.** The 74,349-row occurrence-level dataset
is deposited in Zenodo, which gives it a persistent identifier and a durable archive. Ordinary
repository storage gives it neither.

What **is** here: its complete dictionary, its schemas, its provenance, its checksums, and these
instructions.

---

## Status

| | |
|---|---|
| deposit package | **built and verified** |
| Zenodo record | **not yet created** — `ZENODO_AUTHENTICATION_REQUIRED` |
| DOI | **none** — not reserved, not published, not fabricated |

**Do not cite a DOI for this dataset yet.** When the Zenodo record is published, its DOI and
version will be added here, to `CITATION.cff`, and to the manuscript's data-availability statement.
Until then the correct statement is that the dataset *has been prepared for deposition*.

## What the deposit contains

| file | what it is |
|---|---|
| `portabilityrisk_occurrence_portability_v1.tsv` | the dataset: **74,349 rows × 19 columns** |
| `PORTABILITYRISK_DATA_DICTIONARY_V1.tsv` | every field: type, nullability, units, permitted values, denominator, derivation, source, interpretation, caveats, observation status |
| `PORTABILITYRISK_DENOMINATOR_DICTIONARY_V1.tsv` | the six denominators |
| `PORTABILITYRISK_EVIDENCE_LEVEL_DICTIONARY_V1.tsv` | the five evidence levels and what each does **not** establish |
| `PORTABILITYRISK_SOFTWARE_AND_DATABASE_VERSIONS_V1.tsv` | every tool and database, pinned |
| `PORTABILITYRISK_PROVENANCE_MAP_V1.tsv` | artefact → script → frozen protocol |
| `PORTABILITYRISK_CLAIM_TO_EVIDENCE_CROSSWALK_V1.tsv` | headline claim → field |
| `PORTABILITYRISK_NMV5_TRANSFORMATION_RECEIPT_V1.json` | proof the deposited file equals the private canonical in every retained cell |
| `schemas/*.schema.json` | JSON Schema per table |
| `SHA256SUMS`, `MANIFEST.tsv` | checksums and per-file manifest |
| `verify_deposit.py` | re-derives every headline denominator from the archive alone |

The same dictionaries and schemas are mirrored in [`data/`](data/) here, so you can read the field
definitions before downloading anything.

## Verifying what you download

```bash
tar xzf portabilityrisk_nmv5_dataset_v1.tar.gz
cd portabilityrisk_nmv5_dataset_v1
sha256sum -c SHA256SUMS
python verify_deposit.py
```

Archive SHA-256 of the prepared package:
`990b13da5f0b8794a571eb06f969a17b07ac4e3e3c4a214f172ab838c58196eb`

`verify_deposit.py` needs only the Python standard library. It re-derives the denominators from the
dataset rather than reading them from a summary, validates each table against its schema, re-hashes
every file, and reports its own disagreement count. It should report **16 checks, 0 disagreements**.

## What it must reconcile to

74,349 occurrences · chromosome 35,140 · plasmid 39,209 · A 18,837 · B 16,303 · C 7,170 ·
D 6,043 · E 25,996 · A+B = 35,140 · C+D+E = 39,209.
Cohort: 6,288 genomes · 6,285 BioSamples · 2,283 BioProjects · 109 species · 12,811 replicons.

If your copy does not reconcile, it is not this dataset.

## Licence

The dataset and its metadata are **CC BY-NC-ND 4.0**. You may share them unmodified, for
non-commercial purposes, with attribution. **Verification is unaffected** — downloading,
re-running the checks and reporting a disagreement are reading, not adapting. Publishing a modified
or derived table, pooling it into an aggregated resource, or commercial use requires permission.

The **code in this repository is Apache-2.0**, which is a different licence with different terms.
Running it over the data does not place the output under Apache-2.0. Source NCBI records remain
governed by NCBI. See [`LICENSING_AND_REUSE_BOUNDARIES.md`](LICENSING_AND_REUSE_BOUNDARIES.md).

## Why the dataset is not committed here

Three reasons, in order of weight. A repository has no persistent identifier, so a citation to it
is a citation to a moving target. Git stores every version of a large table forever, and this one
is 17 MB. And a deposit carries a licence, a version and a landing page that a directory listing
does not.
