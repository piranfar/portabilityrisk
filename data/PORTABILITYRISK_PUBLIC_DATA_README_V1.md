# PortabilityRisk public data resource

Everything needed to check the numbers in the manuscript without trusting the manuscript.

| file | what it is |
|---|---|
| `PORTABILITYRISK_DATA_DICTIONARY_V1.tsv` | every column of every published table: type, definition, denominator, nullability. **273 columns, 100 % curated** |
| `PORTABILITYRISK_TABLE_CATALOGUE_V1.tsv` | every published table: rows, columns, grain, SHA-256 |
| `PORTABILITYRISK_DENOMINATOR_DICTIONARY_V1.tsv` | the five denominators, and which question each answers |
| `PORTABILITYRISK_EVIDENCE_LEVEL_DICTIONARY_V1.tsv` | the five evidence levels, and what each does **not** establish |
| `PORTABILITYRISK_PROVENANCE_MAP_V1.tsv` | artefact → producing script → governing frozen protocol → inputs |
| `PORTABILITYRISK_CLAIM_TO_EVIDENCE_CROSSWALK_V1.tsv` | every headline claim → the table, field and digest it is read from |
| `schemas/*.schema.json` | JSON Schema per table, for machine validation |

---

## The one number readers get wrong

**36.586 % is `3,569 / 9,755` genome-collapsed events.**

A genome-collapsed event is **one event per (genome, compartment)**: each of the 6,288 genomes
contributes at most one plasmid event and at most one chromosomal event, however many resistance
genes it carries in that compartment. That gives 3,569 plasmid events + 6,186 chromosomal events =
9,755, and 3,569 / 9,755 = 36.586 %.

**It is not the mean of per-genome plasmid percentages.** That statistic is **35.932 %** — a
different number answering a different question. `genome_level_summary.tsv` has a `pct_plasmid`
column; averaging it does **not** reproduce 36.586 %, and it is not meant to.

Both differ again from the occurrence-weighted plasmid share, **52.736 %**, which counts genes
rather than genomes. All three are correct. Quote the denominator every time.

The same trap exists on the chromosomal side: **46.39 %** of chromosomal *occurrences* are within
10 kb of a mobile-element marker, but only **30.14 %** of chromosomal *blocks* are, because
marker-positive blocks carry 2.46 resistance genes on average against 1.23 for marker-negative ones.

---

## The three things that must never be conflated

| level | what it establishes | what it does **not** |
|---|---|---|
| **homology proximity** | a mobile-element marker is nearby | that an intact element is there |
| **structural IS detection** | a complete element, transposase ORF plus bilateral terminal inverted repeats, is fully inside the window | that it moved, moves, or can move |
| **experimental mobility** | — | **never measured.** No mating, transposition or transfer assay exists in this study |

Absence is always detection-limited. A censored row means *no qualifying element was detected
within 10 kb under the frozen definition*, never *no element is present*. Class C is not
"non-mobilizable"; class A is not "immobile". Both are statements about a marker database.

---

## The two structural positives that are homology-negative

12,034 occurrences carry a structurally complete, fully contained insertion sequence. 12,032 of
them are class B. **Two are class A** — the structural caller resolved a complete element where
profile homology found nothing within 10 kb:

| occurrence | species | structural distance |
|---|---|---:|
| `GCF_015135855.1 \| NZ_AP022446.1 \| blaACT-102` | *Enterobacter kobei* | 4,334 bp |
| `GCF_009738085.1 \| NZ_CP033102.1 \| oqxB` | *Enterobacter hormaechei* | 4,031 bp |

Both are homology-censored, so their transposases did not produce a profile hit passing the frozen
threshold — most plausibly divergent or unclassified transposases that the structural caller
resolves from element architecture instead of sequence similarity.

Two consequences, both worth stating. First, the structural set is **very nearly, but not exactly**,
a subset of the homology set: the methods disagree in 2 of 12,034 cases, and in both the structural
method is the more sensitive one. Second, **neither occurrence is in a headline species group** —
both fall in `other` — so the 73.80 % corroboration rate and every species contrast are unaffected.

---

## What is published here, and what is not

**Published.** All result and summary tables the manuscript quotes, including the two
occurrence-level endpoint tables (`nmis_occurrence_endpoints.tsv`, 35,140 rows;
`nmdist_occurrence_block_distances.tsv`, 35,140 rows), each with a full schema and dictionary
entry.

**Not published yet.** The upstream occurrence-level dataset — `determinant_occurrences.tsv`
(37 MB), `determinant_portability_classes.tsv` (17 MB), `arg_mge_neighbourhood.tsv` (8 MB) and
their companions, about 95 MB. These are listed in the provenance map with `PRIVATE — pending
deposit`. They are withheld until schema, dictionary, licence status and disclosure audit all pass;
publishing them ahead of that would substitute a file drop for a deposit.

Consequence, stated plainly: most scripts under `code/context/` **cannot run** against this
repository alone. They are published as the method of record. Of 61 scripts, 7 are executable from
the public tree; the rest need the deposit. `SCOPE_AUDIT.txt` lists which is which.

---

## Reuse

Underlying records are public NCBI data under their own terms. The derived tables here are
published for inspection, verification and reproducibility assessment. **No licence or permission
for redistribution, relicensing or commercial reuse is granted by this repository** — see
`DATA_LICENSE.md` and `DATA_NOTICE.md`. Provenance fields travel with the values wherever they go.
