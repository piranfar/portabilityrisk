# NMIS — structural IS reconstruction: manuscript interpretation

**Dated 2026-08-22.** Executed under `NMIS_FROZEN_PROTOCOL_V1.json` =
`5438045a3b73d123347fcd60b2456779f050c2a8fae9cd016782b8b6168b03a3` with
`NMIS_AMENDMENT_001_OCCURRENCE_WINDOW_CONTAINMENT.json` =
`72ce50c2cfd512c1bc6cb74d29830f70474ec5f1a12056a9e744883246e67f80`, both frozen before any block
was scored. Independent verification: **35 checks, 0 disagreements.**

**Boundary:** this is structural corroboration and sensitivity analysis of the homology-based
NM-DIST signal. It does not redefine, replace or reclassify class B, and no frozen quantity moved.

---

## 1. Census and endpoint accounting

**21,955 / 21,955 blocks**, 0 tool failures, 14,426 ISEScan elements (type `c` 10,769, `p` 3,657).

### Element categories — reconcile exactly to 14,426

| category | n |
|---|---:|
| complete and fully contained | **7,923** |
| complete in shared block but crossing the occurrence-window boundary | **2,384** |
| partial or structurally incomplete / boundary-limited | **4,119** |
| **total** | **14,426** |

### Occurrence endpoint states — reconcile exactly to 35,140

| state | n |
|---|---:|
| complete and fully contained | **12,034** |
| complete in shared block but crossing the occurrence-window boundary | **716** |
| partial or boundary-limited | **1,902** |
| no structurally resolved IS detected under the frozen ISEScan definition | **20,488** |
| tool failure | **0** |
| **total** | **35,140** |

Class A **18,837** and class B **16,303** unchanged; 21,955 blocks; every occurrence carries
exactly one state; nothing silently excluded.

---

## 2. Structural endpoint, block-balanced

| group | F(1 kb) | F(2 kb) | F(5 kb) | F(10 kb) | RMD | median |
|---|---:|---:|---:|---:|---:|---|
| ***A. baumannii*** | **0.3531** | **0.3956** | **0.4315** | **0.4555** | 5,873 bp | not reached |
| *P. aeruginosa* | 0.0276 | 0.0493 | 0.0817 | 0.0954 | 9,281 bp | not reached |
| *Klebsiella* group | 0.0633 | 0.0736 | 0.1021 | 0.1435 | 8,988 bp | not reached |

### Primary contrasts — BioProject bootstrap (B = 2,000, seed 20260822), Holm-corrected

| contrast | 1 kb | 2 kb | 5 kb | 10 kb |
|---|---:|---:|---:|---:|
| **N1** *A. baumannii* − *Klebsiella* | **+0.2898** | **+0.3221** | +0.3293 | +0.3119 |
| **N2** *A. baumannii* − *P. aeruginosa* | **+0.3255** | **+0.3463** | +0.3498 | +0.3601 |

All eight CIs exclude zero; all Holm-corrected **p = 0.001**.

> ### GATE VERDICT: **SUCCESS**
> Both primary contrasts are positive with CIs excluding zero at **both** 1 kb and 2 kb after Holm
> correction — the registered success condition.

---

## 3. Homology versus structural

| group | homology F(1 kb) | structural F(1 kb) | homology F(10 kb) | structural F(10 kb) |
|---|---:|---:|---:|---:|
| *A. baumannii* | 0.5224 | 0.3531 | 0.6323 | 0.4555 |
| *P. aeruginosa* | 0.0993 | 0.0276 | 0.1717 | 0.0954 |
| *Klebsiella* | 0.1039 | 0.0633 | 0.2390 | 0.1435 |

**The five registered questions:**

1. **Does the *A. baumannii* short-range excess remain positive at 1 and 2 kb?** **Yes** —
   +0.2898 and +0.3221 against *Klebsiella*, +0.3255 and +0.3463 against *P. aeruginosa*.
2. **Do both primary CIs exclude zero after Holm?** **Yes**, at every landmark.
3. **Is the NM-DIST ordering retained?** **Yes** — *A. baumannii* ≫ *Klebsiella* > *P. aeruginosa*
   at every landmark, identical to the homology ordering.
4. **How much of the homology signal is structurally corroborated?** **12,032 of 16,303 class-B
   occurrences — 73.80 %.** By species: *A. baumannii* **86.01 %**, *Klebsiella* 65.46 %,
   *P. aeruginosa* 64.24 %. Only **2** occurrences are structurally positive while
   homology-negative, so the structural set is very nearly a strict subset of the homology set —
   the two methods agree in direction and disagree only in stringency.
5. **Does IS6 dominance survive the structural gates?** **Yes, and it sharpens.** IS6 retains
   **90.4 %** of its elements through the full structural gate (4,385 of 4,848) — the highest of
   any abundant family. IS4 retains 70.4 %, IS5 85.2 %, IS91 53.7 %, IS3 41.3 %, and IS1380
   **0 %** (371 elements, none structurally complete). IS6 is therefore not merely the most
   abundant family but the most consistently completely resolved.

**Short-range concentration is preserved and is more distinctive than under homology.** Under the
structural endpoint **86.9 %** of *A. baumannii* detections fall within 2 kb, against **51.7 %**
for *P. aeruginosa* and **51.3 %** for *Klebsiella*.

---

## 4. Interpretation, stated separately

**A. What is directly established.** In 21,955 chromosomal context windows analysed under one
uniform pinned environment, structurally resolved complete insertion sequences — `type c` with a
complete transposase ORF and bilateral resolved terminal inverted repeats, fully contained within
the individual ARG occurrence's own ±10 kb window — are present for **12,034 of 35,140** chromosomal
acquired-ARG occurrences. Their distance distribution differs sharply by host: *A. baumannii*
reaches F(1 kb) = 0.3531 while both comparators remain below 0.07.

**B. What is strengthened.** The NM-DIST finding. The short-range chromosomal signature survives a
substantially stricter endpoint, retains the same species ordering, and is corroborated
structurally in **73.80 %** of class-B occurrences overall and **86.01 %** in *A. baumannii*. The
IS/transposase attribution from NM-DIST is confirmed at element level, and IS6 emerges as the
dominant structurally complete family.

**C. What remains homology-based.** Class B itself, the 46.39 % occurrence-level and 30.14 %
block-level estimates, and the NM-DIST curves. **NMIS does not replace them.** The 26.20 % of
class-B occurrences without a contained complete element are **not** reclassified: a homology
marker without a fully resolved element remains valid evidence of mobile-element context under the
frozen class definitions.

**D. What is not measured.** Transposition. Transfer. Horizontal gene transfer. Mobilisation.
Phenotype. Element activity. Copy-number dynamics. Any distance beyond 10 kb. Absence is reported
as *"no structurally resolved IS detected under the frozen ISEScan definition"*, never as
biological absence.

**E. Is class B unchanged?** **Yes — 16,303, verified.** Class A 18,837, total 35,140, blocks
21,955, all unchanged. NMIS added an orthogonal endpoint; it reclassified nothing.

**F. Does this license a new headline, or only a sensitivity statement?** **A new headline, with a
precisely bounded claim.** This is the first structurally resolved, census-scale demonstration that
the chromosomal mobile-element context of acquired resistance genes differs by host not only in
frequency but in *structural completeness and spatial proximity*. It is a headline because it
survives the strictest available endpoint at census scale with an independent verification of zero
disagreements — not merely a robustness footnote.

---

## 5. Permitted and prohibited wording

**Three distinct things must never be conflated:**

| level | what may be said |
|---|---|
| homology proximity | "within X bp of an IS/transposase homology marker" |
| **structural detection** | "a structurally resolved complete insertion sequence — complete transposase ORF and bilateral terminal inverted repeats — lies fully within the occurrence's ±10 kb window" |
| experimental mobility | **NOT MEASURED.** No transposition, transfer or activity assay exists anywhere in this study. |

**Permitted:** *structurally resolved*, *structurally complete*, *fully contained within the
occurrence window*, *conjugation-consistent* (plasmid side), *no structurally resolved IS detected
under the frozen ISEScan definition*.

**Prohibited:** *active transposon*, *demonstrated transposition*, *mobilised*, *transferred*,
*HGT event*, *proof of mobility*, and any statement treating non-detection as biological absence or
treating homology-only calls as false positives.

---

## 6. Limitations

1. The structural endpoint is **detection-limited**: ISEScan under the frozen definition, HMMER
   3.3.2 profiles, and a ±10 kb window. A shorter window truncates large elements — this is why
   2,384 complete elements were classed as crossing the occurrence-window boundary and excluded
   from the primary endpoint.
2. **716 occurrences** have a complete element in their shared block that crosses their own window
   and are therefore censored. They are reported separately and are available for a secondary
   sensitivity analysis.
3. No IS boundary was inferred and no flanking sequence retrieved; boundary-touching elements
   remain partial by construction.
4. Medians are **not reached** in any group under the structural endpoint, because F(10 kb) < 0.5
   everywhere. Restricted mean distance is reported instead.
5. The result concerns three frozen species groups; all other species remain descriptive.

## 7. Provenance

| artefact | SHA-256 |
|---|---|
| `NMIS_FROZEN_PROTOCOL_V1.json` | `5438045a3b73d123347fcd60b2456779f050c2a8fae9cd016782b8b6168b03a3` |
| `NMIS_AMENDMENT_001_…json` | `72ce50c2cfd512c1bc6cb74d29830f70474ec5f1a12056a9e744883246e67f80` |
| `nmis_occurrence_endpoints.tsv` | `1fc2a6748f2a15a156d36ddb69d33b6d4225dd26798c2e680341f0df83a9c827` |
| `nmis_primary_estimates.tsv` | `395ee5162f09bcdc10c8f078a8a55a4d73b115eb08f90df1c30d302096c22afd` |
| `nmis_species_contrasts.tsv` | `8e3602ba8e9d3a2a62f06b8b2822f60cb406144309943c4291e859ed8993e253` |
| `nmis_sensitivity_results.tsv` | `5ee44f9f7e746315111dcc8e03a60205fe1aec444ea3917a2080b01aad9d1d20` |
| `nmis_vs_nmdist_comparison.tsv` | `7eccd77f224ef056d29717d2ec553addef42041945daceb8c173553da9c51f9f` |
| `NMIS_RESULT_RECEIPT_V1.json` | `3725c18716ecd4228e5c5ca1cd261ffafb7d5e18d8ab2a93d72f765decc2b8cc` |

Environment reproduced NM-V1 exactly: ISEScan 1.7.3 (`isescan.py` digest `a6601aab…`, byte-identical),
HMMER **3.3.2** (pinned after the default solve returned 3.4), BLAST 2.17.0, FragGeneScan 1.32,
Biopython 1.88, Python 3.14.7. Command line unchanged:
`isescan.py --seqfile <block.fna> --output isescan_out --nthread 1`.
