# NM-V1 Amendment 001 — QC arm restoration and adjudication evidence depth

**Dated 2026-08-21, America/New_York (UTC−04:00).**

Corrects the adjudication package only. **No tool was rerun.** The frozen sample, the frozen
biological rules and the reference tool calls are unchanged and were not recomputed.

---

## 1. What went wrong

### 1.1 The QC arm was short by 20 blocks

`NMV1_GROUNDTRUTH_RULES_FROZEN.json` fixes a quality-control arm of **100** blocks:
40 `AUTO_IS_MOBILE`, 40 `AUTO_QUIESCENT`, 20 `AUTO_INTEGRON`. Only **80** were produced.

**Cause.** The automatic-label rules are evaluated as an `elif` chain, and `AUTO_IS_MOBILE`
precedes `AUTO_INTEGRON`:

```
elif r["hmm_is_pos"] and r["isescan_complete"] > 0:   -> AUTO_IS_MOBILE
elif r["hmm_int_pos"] and r["if_complete"] > 0:       -> AUTO_INTEGRON
```

Any block carrying **both** a complete IS and a complete integron is captured by the first
branch. Because complete integrons in this cohort almost always co-occur with a complete IS,
the `AUTO_INTEGRON` cell was left with **zero** eligible blocks and the QC arm drew 40 + 40 + 0.

This is a precedence defect in the label assignment, not in the biology: the affected blocks
are correctly labelled mobile either way. What was lost is the ability to audit the
*integron* automatic rule specifically.

### 1.2 The package carried counts only

The V1 package exposed marker **counts** and nothing else — no coordinates, no element
boundaries, no inverted-repeat evidence, no attC positions, no distance to the resistance-gene
interval, no feature map. A biological adjudicator cannot judge whether a block sits in mobile
context from integer counts alone. The package was insufficient to support independent ground
truth, which is the entire purpose of NM-V1.

## 2. State at the time of correction

| | |
|---|---|
| aggregate ground-truth performance computed | **none** |
| adjudication begun | **no** |
| unblinding key opened | **no — V1 key remains sealed** |
| tool outputs altered | **none** |
| frozen sample altered | **none** |
| frozen biological rules altered | **none** |
| reference calls recomputed | **none** |

The correction affects **QC sampling and package presentation only**.

## 3. Superseded artefacts, preserved

| artefact | SHA-256 | status |
|---|---|---|
| `NMV1_ADJUDICATION_BLINDED_PACKAGE_V1_SUPERSEDED_INSUFFICIENT_EVIDENCE.tsv` | `436e295cbea45139e984c8f12e7df7e322eed2d24fefaad9e641412679bae975` | superseded, retained |
| `NMV1_ADJUDICATION_UNBLINDING_KEY_V1_SUPERSEDED_INSUFFICIENT_EVIDENCE.tsv` | `731fa74c820e2b3f15260d97be1d7ee9bc01e14ebd065a8cdd42133c4aa89c9e` | **sealed**, superseded, retained |
| `NMV1_ADJUDICATION_RUBRIC_V1_SUPERSEDED_INSUFFICIENT_EVIDENCE.md` | `a48494c56d86582f328bb16f78ad91de395ff730de5aa2886e3801173a36dd8f` | superseded, retained |

Nothing was overwritten or deleted.

## 4. The correction

### 4.1 QC arm restored to 100, three mutually exclusive cells

Selection is deterministic: candidates sorted by `SHA-256("<block_id>|20260821|<cell>")`, first
*n* taken. Cells are processed **INTEGRON → IS → QUIESCENT** so that a block claimed by an
earlier cell cannot reappear in a later one.

| cell | rule | eligible pool | selected |
|---|---|---:|---:|
| `QC_AUTO_INTEGRON` | auto-labelled **and** ≥1 complete IntegronFinder structure — complete-IS blocks admitted | 154 | **20** |
| `QC_AUTO_IS` | `AUTO_IS_MOBILE`, complete IS, **no** complete integron | 235 | **40** |
| `QC_AUTO_QUIESCENT` | `AUTO_QUIESCENT` | 399 | **40** |

Admitting a complete-IS block into the integron cell **does not change its biological label**.
It identifies which automatic evidence rule is under audit. Eligible pools, selection rule and
per-cell selected-token hashes were written to `NMV1_QC_MANIFEST_V2.json` **before** the
package was generated.

### 4.2 Reconciliation

| quantity | n |
|---|---:|
| unique sampled blocks | **1,283** |
| automatic-label blocks | 788 |
| mandatory-adjudication blocks | 495 |
| 788 + 495 | **1,283** ✔ |
| QC rows — a concealed review subset **of** the 788 | 100 |
| **package rows = 495 + 100** | **595** |

The 100 QC rows are **not additional unique blocks**. The sampled total remains 1,283.

### 4.3 Evidence depth

Every one of the 595 cases now carries: opaque token; topology and boundary warning; window
length; anonymised resistance-gene intervals; per-feature coordinates, strand and e-value;
complete vs partial IS status; terminal inverted-repeat coordinates, length and identity;
transposase ORF coordinates and length; integrase ORF coordinates; attC site count and
positions; integron structural class (complete / CALIN / In0); distance from each feature to
the nearest gene interval; tool execution status; and a standardised SVG feature map.

**Coordinate convention:** all coordinates shown are 1-based and relative to the block.
ISEScan and IntegronFinder report block-relative coordinates natively; HMM features and gene
intervals are chromosome-relative in the source tables and were converted by subtracting
`block_start`.

## 5. Blinding

Withheld: species, BioProject, accession, gene name, original portability class, sampling
stratum, tool identity, aggregate results. Methods appear only as **X**, **Y** and **Z**.

Automated leakage check across **1,034** identity-bearing strings: **0 leaks in the casebook,
0 in the spreadsheet.**

QC tokens are interleaved through the token ordering (positions 0–594 of 595), not clustered,
so concealed controls are indistinguishable from mandatory cases.

## 6. What is unchanged

The seven permitted adjudication outcomes, the decision rubric, the frozen sample of 1,283
blocks, the frozen strata and sampling fractions, the automatic-label rule definitions, and
every reference tool call. Amendment 001 changes **which** auto-labelled blocks are audited and
**how** the evidence is presented — nothing about what the evidence is.
