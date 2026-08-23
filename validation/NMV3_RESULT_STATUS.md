# NM-V3 — tool and database robustness: result
> **REDACTED PUBLIC DERIVATIVE. Canonical private artefact: docs/nature_microbiology/NMV3_RESULT_STATUS.md, SHA-256 5850bfd54e2c35c1edc233f5e9f734c46a74ff8b16cccb08eb85a248283f9874. Infrastructure identifiers replaced by [REDACTED:...]; no scientific content altered.**

**Dated 2026-08-22.** Executed on the isolated instance `[REDACTED:INSTANCE_NAME]` under
`NMV3_FROZEN_DESIGN.json` (`8d96b72304a9580163132e07ecf16ede65cebe84760b0781a19d3e438f992a2d`)
and `NMV3_ARM_I_FROZEN_CONFIG.json` (`2351be401b5b80b727b48f1703018c6b487b8d102d70d6d1d0c96739590b7619`).
No frozen gate or classification threshold was changed.

---

## 1. Verdicts

| gate | verdict |
|---|---|
| **C10** — class C/D/E boundaries stable to mobility-marker database and version | **PASS** |
| **C09** — conjugative vs marker-negative convergence contrast | **PASS** |

Worst-case replicon category change across every executed arm: **0.0000 %**, against NM-0's
≤ 5 % bound. No headline direction reversed. The location layer is untouched — it is not a
MOB-suite output and was not recomputed.

## 2. The baseline reproduced the frozen run exactly

The original mobility layer was re-derived from scratch on a clean machine, from a fresh install
and a fresh database download:

| | |
|---|---|
| replicons | **6,621** |
| marker-triple mismatches vs the frozen 2026-08-21 run | **0** |
| mobility-category mismatches | **0** |

This is verification, not re-use. It is stronger than checking a surviving environment's digests,
because it demonstrates the layer is reconstructible rather than merely preserved.

Two provenance caveats. The recorded `mob_typer_sha256` `43408d89…` **cannot** match by
construction: `mob_typer` is a 248-byte shebang wrapper embedding the environment path. The
authoritative identity is the conda package, `mob_suite-3.1.9-pyhdfd78af_1`, which is the exact
build the frozen run used. Likewise the recorded `database_recursive_sha256` `aae85a3b…` uses an
unknown digest algorithm; database identity was instead established by per-file content hashes.

## 3. Arm T — tool version

**MOB-suite 3.1.8 vs 3.1.9, database held constant at v3.1.8.**

| | |
|---|---|
| marker-call mismatches | **0 / 6,621** |
| replicon category transitions | **0 (0.0000 %)** |
| occurrence transitions | **0 / 39,209** |
| transition matrix | diagonal |

Raw file digests differ; **sorted content is identical** across all three configurations
(`2c7e6e81…`). The difference is row order only — precisely the reproducibility fix documented in
the 3.1.9 changelog, which the frozen design named in advance as a known confound rather than
discovering after the fact.

## 4. Arm D — database version

**Database v2.0.0 vs v3.1.8, tool held constant at 3.1.9.** Four independent lines of evidence,
reported separately as required:

**(a) Analytic comparison of the marker databases.** Every file that determines C/D/E is
**byte-identical** across versions: `mob.proteins.faa` (relaxase), `mpf.proteins.faa` (MPF),
`orit.fas` (oriT), plus `rep.dna.fas`, `repetitive.dna.fas`, `ncbi_plasmid_full_seqs.fas` and
`host_range_literature_plasmidDB.txt`. **Only `clusters.txt` differs** — MOB-cluster assignments,
which drive cluster identifiers and host-range prediction, not marker detection. The database axis
is therefore *structurally incapable* of moving a C/D/E assignment.

**(b) Empirical smoke result.** 20 plasmids: output byte-identical to baseline.

**(c) Empirical full-census result.** 6,621 replicons: **0 marker mismatches, 0 transitions**,
output byte-identical to baseline.

**(d) Limitation.** This bounds the only two MOB-suite database versions that exist. It says
nothing about a *future* database, which could change the marker files. C10 must be worded to that
scope.

### Frozen-design departure, preserved explicitly

- the frozen design named Zenodo record **3785613**;
- that record is **incomplete** — `repetitive.dna.fas` is absent (13 files);
- the complete sister record **3786915** was used (15 files);
- the **13 shared files are byte-identical** between the two records;
- every C/D/E-determining marker database is byte-identical across 3785613, 3786915 and v3.1.8
  where present.

The frozen design was **not rewritten**. The departure is recorded here and in the receipt.

## 5. Arm I — independent implementation

**Prodigal 2.6.3 → MacSyFinder 2.1.6 / CONJScan 2.1.0.** 6,621 plasmids, **0 failures, 0
zero-CDS**, 1,006,053 proteins (5,300 partial, retained by design). Reported as **concordance
only** — neither tool is a truth standard, and a detection by one and not the other is a
disagreement between two homology-based predictors, not an error by either.

| axis | MOB-suite + | CONJScan + | both | MOB-only | CONJScan-only | raw agreement | Cohen's κ |
|---|---:|---:|---:|---:|---:|---:|---:|
| relaxase | 4,613 | 4,369 | 4,322 | 291 | 47 | 94.90 % | **0.883** |
| MPF | 4,378 | 3,657 | 3,602 | 776 | 55 | 87.45 % | **0.740** |
| **class E definition** | **3,937** | **3,657** | **3,594** | **343** | **63** | **93.87 %** | **0.875** |

**Occurrence-weighted:** 24,012 of 39,209 plasmid-borne occurrences — **61.24 %** — sit on class-E
replicons corroborated by *both* implementations, against 66.30 % under MOB-suite alone. 1,984
occurrences rest on class-E calls that only MOB-suite makes.

119 replicons carry only decayed (`dCONJ_*`) systems; these were recorded separately and never
counted as MPF evidence.

**Scope limit:** CONJScan has no oriT model, and no independent oriT scheme is installable offline.
Arm I therefore covers the **complete** definition of class E (relaxase AND MPF) but only
**partially** covers C and D, which depend on oriT. The oriT axis is reported as **NOT
INDEPENDENTLY TESTED**, not substituted.

**Parsing correction, disclosed.** A first extraction read only system-level model names and
counted relaxase evidence as a standalone `MOB` system. That undercounts badly, because CONJScan's
`T4SS_type*` systems carry the relaxase as a *mandatory* component (`T4SS_MOBF`, `T4SS_MOBQ`,
`T4SS_MOBP1`…). Evidence was re-extracted at **gene level**. No Prodigal or MacSyFinder setting was
changed — this was scoring, not configuration — and the error was caught before any result was
scored.

## 6. E1 / E2 evidence hierarchy — stable

| tier | evidence | occurrences | replicons |
|---|---|---:|---:|
| **E1** | relaxase + MPF, oriT not detected | **16,397** | 2,828 |
| **E2** | relaxase + MPF + oriT detected | **9,599** | 1,109 |
| **E** | frozen class E | **25,996** | 3,937 |

**66.301104 %** of plasmid-borne occurrences are class E; **24.481624 %** are the triply
corroborated E2 subset; E2 is **36.924912 %** of class E and E1 is **63.075088 %**. Identical under
both Arm T and Arm D. No detected oriT is not equivalent to biological absence of oriT.

## 7. C09

Conjugative minus marker-negative difference in ≥3-drug-class carriage, BioProject-clustered
bootstrap (B = 2,000, seed 20260822):

| configuration | difference | 95 % CI | direction | CI excludes 0 |
|---|---:|---|---|---|
| baseline | **+19.861 pp** | 14.160 – 25.414 | preserved | yes |
| Arm T | +19.861 pp | 14.007 – 25.438 | preserved | yes |
| Arm D | +19.861 pp | 14.188 – 25.655 | preserved | yes |

**PASS.** C09 remains Tier 2 and stays in the abstract.

## 8. Manuscript consequences

| statement | consequence |
|---|---|
| A/B/C/D/E = 74,349; plasmid share 52.736419 % | **unchanged**, verified |
| classes A (18,837) and B (16,303) | **untouched** |
| 66.3 % conjugation-consistent | **retained**, with the §6 operational wording; add that **61.24 %** is corroborated by an independent implementation |
| class E as Tier 1 descriptive, operationally defined | **supported** — zero version or database sensitivity, κ = 0.875 against an independent tool |
| C09 / Figure 6 | **retained**, Tier 2 |
| C10 | **PASS**, worded to the two existing database versions |
| version-stamping every C/D/E figure | **still required** — the result bounds existing versions, not future ones |
| *A. baumannii* discordance, NM-V1C, NM-V2, NM-V4 | **unaffected** — no exposure to the mobility layer |

## 9. Limitations

1. Only two MOB-suite database versions exist; a future release could change the marker files.
2. No newer tool than 3.1.9 exists, so the version axis is necessarily backward-looking.
3. The oriT axis has no independent implementation.
4. MOB-suite is more permissive on MPF than CONJScan (776 vs 55 one-sided calls); this is a
   documented difference between predictors, not evidence that either is wrong.
5. `FROZEN_PORTABILITY_CONTEXT_PROTOCOL_V1.json` still pins no mobility-tool digest and no stop
   condition gates it. NM-V3 bounds that gap but does not close it.

## 10. Provenance

| artefact | SHA-256 |
|---|---|
| `NMV3_FROZEN_DESIGN.json` | `8d96b72304a9580163132e07ecf16ede65cebe84760b0781a19d3e438f992a2d` |
| `NMV3_ARM_I_FROZEN_CONFIG.json` | `2351be401b5b80b727b48f1703018c6b487b8d102d70d6d1d0c96739590b7619` |
| `NMV3_RESULT_RECEIPT.json` | `e4de8950db7f8a16275b1395e361c43a569a2581f534c00afc640e7f53a6c247` |
| `NMV3_TRANSITION_MATRICES.tsv` | `dd48dcbff485b2ba5e64499d3b9607b92c5d15208953dc6ed278cce596a860f1` |
| `NMV3_VERIFICATION_REPORT.txt` | see repository |
| deliverable archive `nmv3_deliver.tar.gz` | `7a58e21c7b4b9d24f9c00ade189e85b2626e7b3f5200b6f01d6dba9bd1e3a691` |
| required-set aggregate (server = local) | `db9088e75a7bf94715d03ae1b68d2cc1d386b09f2392747612fb2e6fcd0b6803` |

Independent verification: **34 checks, 0 disagreements**, via a code path that imports nothing from
the scorer.

Environment: `mob_suite-3.1.9-pyhdfd78af_1` / `mob_suite-3.1.8-pyhdfd78af_1`,
`blast-2.17.0-h66d330f_0`, `mash-2.3-hf85e966_11`, `prodigal-2.6.3-h577a1d6_11`
(exe `b595f285…`), MacSyFinder 2.1.6, HMMER 3.4, CONJScan 2.1.0 (recursive `8c2c529d…`).
