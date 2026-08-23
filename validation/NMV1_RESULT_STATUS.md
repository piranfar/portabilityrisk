# NM-V1 — result status

**Dated 2026-08-22. Written after the corrected confirmatory audit was scored against the gate
frozen before any case was drawn.**

Module: NM-V1, *Stratified mobile-genetic-element validation*. Gates class **B**, claims **C05**
and **C06**, and the **46.39 %** chromosomal-MGE headline.

---

## Verdict

| | |
|---|---|
| **NM-V1C gate verdict** | **SUCCESS** |
| primary — design-weighted three-state agreement | **0.9920** |
| stratified bootstrap 95 % CI (B = 2,000, seed 20260822) | **0.9761 – 1.0000** |
| secondary — unweighted agreement | 0.9917 (**119 / 120**) |
| gate rule applied | agreement ≥ 0.90 **AND** bootstrap lower bound ≥ 0.80 |

Both limbs of the frozen SUCCESS rule are met with margin. Stop condition **S1 does not fire**:
class B stands, C05 and C06 keep headline status, and the 46.39 % figure is not demoted.

### Per machine state

| state | n | agreed | unweighted | 95 % CI | weighted |
|---|---:|---:|---:|---|---:|
| MOBILE | 60 | 59 | 0.9833 | 0.9500 – 1.0000 | 0.9825 |
| QUIESCENT | 60 | 60 | **1.0000** | 1.0000 – 1.0000 | 1.0000 |

### Per stratum

| stratum | n | design weight | agreed | unweighted |
|---|---:|---:|---:|---:|
| A | 25 | 4.2400 | 25 | 1.0000 |
| B | 34 | 4.6176 | 33 | 0.9706 |
| D | 1 | 1.0000 | 1 | 1.0000 |
| E | 60 | 5.2667 | 60 | 1.0000 |

### The single disagreement

`NMC-096`, stratum B, engine MOBILE → adjudicator QUIESCENT. The panel showed one Method Y
complete element at 17,922–19,125 with bilateral 20 bp inverted repeats (identity 19) and an
internal ORF, in a 20,858 bp circular block with no boundary problem and normal tool completion.
No note was left. One disagreement in 120 is within every bound the protocol set, and a failure
confined to one state would revise that state only; no state failed.

The adjudicator left **0 notes** and chose **NON_EVALUABLE 0 times**.

---

## Why this result is credible given the earlier FAIL

The registered R3 audit **FAILED at 62/120, three-state agreement 0.5167**. That result stands as
reported and is not erased, passed or rescored.

Amendment 005 diagnosed the cause as instrument contamination: Method X — the profile-based path
under test — was displayed to the adjudicator while the ground-truth engine was structurally
prohibited from using it. The evidence was a perfect dissociation in the quiescent stratum: of 47
engine-QUIESCENT cases, the **38 with a Method X marker drew 0 QUIESCENT calls** and the **9
without drew 9 of 9**.

NM-V1C removed Method X and changed nothing else — same engine, same F3 threshold, same gate
thresholds, 120 **fresh** blocks with zero overlap with any prior package, and **no adjudication
reused**. The quiescent stratum went from **0 / 38** to **60 / 60**.

That is as close to a controlled natural experiment as this design permits, and it confirms the
Amendment 005 diagnosis: the earlier failure measured the instrument, not the engine.

---

## Scope of what was validated — stated before the draw, restated here

Recorded in `NMV1C_FROZEN_PROTOCOL.json` under
`STRUCTURAL_LIMITATION_STATED_BEFORE_DRAWING`, before any case was selected:

> The unused population contains **zero NON_EVALUABLE blocks and zero integron-exclusive (rule C)
> blocks**, because the frozen rules routed all 114 indeterminate and all 8 rule-C blocks to
> mandatory adjudication, so all of them were consumed by the V1/V2 packages.

**Therefore NM-V1C validates the MOBILE versus QUIESCENT discrimination only.** It does **not**
validate the NON_EVALUABLE state, and it yields **no integron-specific estimate**. Any manuscript
sentence citing this result must carry that boundary.

The deferred partial-element question lives exactly in the untested region: 19 of 20 partial-only
cases in the failed audit were called MOBILE against the engine's F3 indeterminate. That remains a
**developmental finding for a future revised classifier**, which requires a new version number and
a separate untouched holdout. It was deliberately not acted on here, because revising F3 in light
of a failed audit would convert validation into model development.

---

## Provenance

| artefact | SHA-256 |
|---|---|
| `NMV1C_FROZEN_PROTOCOL.json` | `b2058877c0f7165f0ab7bf9be7adcb2cb3b5b4c8d2b81143d85ab03fe0f0c04c` |
| `NMV1C_ADJUDICATION_APP_R2.html` | `56581a89c3e7acf0ce39b831c122bacd8d165ebfbdd7aa58a9d8dafda1eee41b` |
| `NMV1C_UNBLINDING_KEY.tsv` | `40e142e27978631da9ec437b9cdd7aeb8c5fda52bdecbbcda5f2fa2db8501efe` |
| `NMV1C_ADJUDICATED_120_RECEIVED.json` | `a572945c9c02d316e00b076ac6f881472be554f5b3dc068cfee0bbdcf24c6b75` |
| `nmv1c_scored_cases.tsv` | `2dabb700bcb5d5722b8397c1afbbf7a4d8aea2b883a2b012609461ee80c864f0` |
| `NMV1C_SCORING_RECEIPT.json` | `a2ae6a8d8d7133574c1534172b8c2ccfcec8068ece4ac7b11ae6a1b48b1c42ae` |
| rule engine (unchanged) | `ed5db383bb0afe1a1a8433886d6666fe72c324975de99c6763a37824d51c2bee` |

Adjudicator: Vahhab Piranfar, blinded. Completed manually; the file was frozen and hashed before
the key was opened. The key was opened only after freezing, and only by the scorer, which verifies
all four input digests and aborts on mismatch.

Reproduce:

```
python audit/ingest/assay_aware_emergence/v2/nm_validation/nmv1c_score.py --dir docs/nature_microbiology
```

## Verification performed on the instrument before delivery

43 automated checks, zero failures (`NMV1C_LEAKAGE_VERIFICATION_REPORT.txt`), plus a functional
test in a browser with storage deliberately disabled, plus a coordinate audit of **608 mapped
elements across all 120 cases — 0 misplaced, 0 out of bounds**.

## Module status after this result

| module | status |
|---|---|
| NM-V1 | **PASS** (this document) |
| NM-V2 | pass |
| NM-V3 | not run |
| NM-V4 / V4C / V4D | pass |
| NM-V5 | not run |
| NM-V6 | design only; requires rewrite against published plasmid-classification tools |
