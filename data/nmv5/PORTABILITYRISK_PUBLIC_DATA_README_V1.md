# PortabilityRisk — NM-V5 occurrence-level dataset

**Version 1.0.0.** Everything needed to check the numbers in the manuscript without trusting the
manuscript.

---

## What this is

`portabilityrisk_occurrence_portability_v1.tsv` — **74,349 rows, 19 columns**. One row is one
**acquired antimicrobial-resistance gene occurrence**: one determinant call, at one coordinate
interval, on one replicon, in one closed complete bacterial genome.

For every occurrence the file records the genome, sample and project it came from, the species,
the replicon and **whether that replicon is a chromosome or a plasmid**, the determinant and its
family and drug class, its coordinates, and its **operational portability class A–E**.

The compartment assignment is **documentation, not prediction**. Every genome is closed, so NCBI
states the molecule type per replicon and the occurrence's coordinates fall inside it. That is why
the location layer has no error rate: **0 unmatched, 0 ambiguous, 0 missing coordinates**. It also
bounds the whole dataset — this holds for closed genomes and does not transfer to draft assemblies.

## Files

| file | what it is |
|---|---|
| `portabilityrisk_occurrence_portability_v1.tsv` | the dataset, 74,349 × 19 |
| `PORTABILITYRISK_DATA_DICTIONARY_V1.tsv` | every column: type, nullability, units, permitted values, denominator, derivation, source, interpretation, caveats, and whether it is observed, computed or operationally classified |
| `PORTABILITYRISK_DENOMINATOR_DICTIONARY_V1.tsv` | the six denominators and which question each answers |
| `PORTABILITYRISK_EVIDENCE_LEVEL_DICTIONARY_V1.tsv` | the five evidence levels and **what each does not establish** |
| `PORTABILITYRISK_SOFTWARE_AND_DATABASE_VERSIONS_V1.tsv` | every tool and database version, pinned |
| `PORTABILITYRISK_PROVENANCE_MAP_V1.tsv` | artefact → producing script → governing frozen protocol |
| `PORTABILITYRISK_CLAIM_TO_EVIDENCE_CROSSWALK_V1.tsv` | headline claim → field it is read from |
| `PORTABILITYRISK_NMV5_TRANSFORMATION_RECEIPT_V1.json` | proof this file equals the private canonical in every retained cell |
| `schemas/*.schema.json` | JSON Schema per table |
| `SHA256SUMS`, `MANIFEST.tsv` | checksums and a per-file manifest |
| `verify_deposit.py` | run it; it re-derives every headline denominator from this archive alone |

## Verify it before you use it

```bash
sha256sum -c SHA256SUMS
python verify_deposit.py
```

The second re-derives every denominator below from the extracted files and reports its own
disagreement count. It needs only Python's standard library.

## The reconciliation this dataset must satisfy

| quantity | value |
|---|---|
| acquired ARG occurrences | **74,349** |
| chromosome | **35,140** |
| plasmid | **39,209** |
| A — chromosomal, no MGE marker within ±10 kb | **18,837** |
| B — chromosomal, ≥1 MGE marker within ±10 kb | **16,303** |
| C — plasmid, no mobility marker | **7,170** |
| D — plasmid, relaxase | **6,043** |
| E — plasmid, relaxase and MPF | **25,996** |
| A + B | **35,140** |
| C + D + E | **39,209** |

Cohort: **6,288 genomes · 6,285 BioSamples · 2,283 BioProjects · 109 species · 12,811 replicons.**

## The one number readers get wrong

**36.586 % is `3,569 / 9,755` genome-collapsed events.**

A genome-collapsed event is **one event per (genome, compartment)**: each of the 6,288 genomes
contributes at most one plasmid event and at most one chromosomal event, however many resistance
genes it carries in that compartment. That gives 3,569 plasmid + 6,186 chromosomal = 9,755, and
3,569 / 9,755 = 36.586 %.

**It is a count of events after collapse. It is not the arithmetic mean of per-genome plasmid
percentages** — that statistic is about **35.93 %**, and it answers a different question. Both
numbers are correct; only the terminology was wrong in earlier drafts, and it is corrected here.

Two values circulate for that mean, and the difference is exact rather than sloppy:
**35.932096 %** is the mean of the *stored* `pct_plasmid` column, which is rounded to two decimals,
and **35.932101 %** is the mean recomputed from the raw counts. They differ by 5 parts in a
million because one averages a rounded column. `verify_deposit.py` computes both and asserts
neither equals the collapsed-event share.

Both differ again from the **occurrence-weighted** plasmid share, **52.736 %**, which counts genes
rather than genomes. Quote the denominator every time.

The same trap exists on the chromosomal side: **46.39 %** of chromosomal *occurrences* lie within
10 kb of a mobile-element marker, but only **30.14 %** of chromosomal *blocks* do, because
marker-positive blocks carry 2.46 resistance genes on average against 1.23 for marker-negative ones.

## What the evidence levels do and do not establish

| level | establishes | does **not** establish |
|---|---|---|
| documented replicon location | the gene **is** on a chromosome or a plasmid | anything about whether it can move |
| homology-based MGE proximity | a mobile-element **marker is nearby** | that an intact element is present |
| structurally resolved complete IS | a **complete element is there**, end to end | that it moved, moves, or can move |
| plasmid conjugative-potential | the **machinery is encoded** | that conjugation occurred or could |
| **experimental transfer** | — | **never measured. No mating, transposition or transfer assay exists in this study.** |

Absence is always detection-limited. A negative means *no qualifying marker or element was detected
under the frozen definition*, never *none is present*. **Class C is not "non-mobilizable" and class
A is not "immobile"** — both are statements about what a marker database returned.

## Two columns are empty by construction, not by omission

`evidence_layer_mobility` is empty for all 35,140 chromosomal rows (47.26 %); `evidence_layer_mge`
is empty for all 39,209 plasmid rows (52.74 %). The two counts are exactly complementary and sum to
74,349. Nothing is missing.

## Relationship to the canonical file

This is a **public deposit derivative**. The private canonical,
`determinant_portability_classes.tsv`, has 20 columns; this has 19. The single dropped column,
`plasmidcall_predicted_location`, was **empty in all 74,349 rows** and named a separate manuscript.

`PORTABILITYRISK_NMV5_TRANSFORMATION_RECEIPT_V1.json` records both digests and proves the rest:
every retained cell identical row-for-row, no row removed, reordered or renamed, every denominator
unchanged. **Its digest differs from the canonical by design and must not be quoted as the
canonical digest.**

## Licence

**CC BY-NC-ND 4.0** — Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International.
See `LICENSE`.

You may share this dataset unmodified, for non-commercial purposes, with attribution. **Verifying
it is unaffected**: downloading it, re-running the checks, re-deriving the denominators and
reporting a disagreement are all reading, not adapting. Publishing a **modified or derived** table,
pooling it into an aggregated public resource, or using it commercially requires permission from
the author.

The **analysis code** is separately licensed **Apache-2.0** and lives at
https://github.com/piranfar/portabilityrisk. Running Apache-2.0 code over this data does not place
the output under Apache-2.0.

**Source genome data** — NCBI assemblies, BioSample and BioProject records — remain governed by
NCBI's terms. This dataset contains accessions and computed values derived from them, not
redistributed sequence. No licence is applied to third-party data or software by this deposit.

## Citation

Cite this Zenodo record with its DOI and version, and the accompanying manuscript. See
`CITATION.cff`.

## Funding

This research received no specific grant from any funding agency in the public, commercial or
not-for-profit sectors. Cloud-computing resources were supported through promotional credits
provided by Oracle Cloud Infrastructure.

## What is not here

Raw sequence, FASTA caches, conda environments, the validation modules' intermediate outputs, and
all restricted material — unblinding keys, blinded casebooks, adjudication instruments, credentials.
The **scoring rubric** for the blinded expert audit is published in the code repository, because an
audit whose rules are secret cannot be checked; the keys that would unblind it are published
nowhere.
