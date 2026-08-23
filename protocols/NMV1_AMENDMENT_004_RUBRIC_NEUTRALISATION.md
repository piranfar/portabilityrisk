# NM-V1 Amendment 004 — rubric neutralisation and the partial-element evidentiary boundary

**Dated 2026-08-22, America/New_York (UTC−04:00). Issued before any adjudication decision was
entered.**

Corrects an inconsistency between the rubric shown to the adjudicator and the frozen decision
engine. Changes **wording only**. No case, token, evidence value, feature map, rule label,
threshold or key is altered.

---

## 1. State at the time of this amendment

| | |
|---|---|
| adjudication decisions entered | **none** |
| agreement computed | **none** |
| unblinding keys opened | **none — all three sealed** |
| rule-engine labels altered | **none** |
| 120-case sample altered | **none** |
| success thresholds altered | **none** |
| tool rerun | **none** |

## 2. The defect

The rubric delivered to the adjudicator defined, verbatim:

> `IS_associated_supported` — "a transposase ORF with **at least one** resolved terminal
> inverted repeat, **or** a complete element reported by a boundary-based method"

The frozen engine requires, verbatim:

> rule **B** — `IS_strong AND NOT INT_strong`, where
> `IS_strong = IS_bilateral_TIR_n > 0 AND IS_complete_orf_n > 0`, and
> `IS_bilateral_TIR_n` = "count of **type = c** elements having **BOTH** terminal
> inverted-repeat coordinate pairs … with irLen > 0"

Everything else falls to rule **F3**, `biologically_indeterminate`.

Two bars differ. The engine counts bilateral TIRs **only on elements ISEScan typed `c`**, so a
`p`-typed element can never satisfy `IS_strong` however complete its inverted repeats are. The
rubric admitted exactly such elements.

### Magnitude

Computed from the blinded package alone, no key opened: **19 of 120 cases (15.8 %)** carry a
partial-only IS with both TIR arms reported, a transposase ORF, no integron evidence and no
boundary or tool problem.

`NMV1A-001` is the case the adjudicator raised. As presented it shows
`methodY_TIR_detail = L:20308-20339 R:21399-21430 len=32 id=21` and
`methodY_ORFs = ORF 20282-21381 len 1100`. Under the delivered rubric that is unambiguously
`IS_associated_supported`; under the engine it is F3.

`IS_associated_supported` maps to **MOBILE**; `biologically_indeterminate` maps to
**NON-EVALUABLE**. The disagreement therefore **crosses a state boundary**, so the three-state
collapse introduced by Amendment 003 does not absorb it. Had the adjudicator followed the
delivered rubric on all 19, primary agreement would have been capped at **101/120 = 0.842** —
inside REVISE — from an instrument inconsistency rather than any scientific disagreement.

## 3. Options considered

**(a) Align the rubric to the engine** — instruct that partial-only is indeterminate.
**Rejected**: it dictates the answer on 19 of 120 cases and makes the audit circular. Agreement
would be manufactured, not measured.

**(b) Align the engine to the rubric** — **rejected**: that would edit frozen ground truth
after inspecting the audit composition, which is tuning truth to the result.

**(c) Neutralise the rubric and return the judgement to the adjudicator** — **adopted, owner
approved.** The rubric no longer asserts any completeness threshold. If disagreement survives,
it is a real finding about where the engine draws its line, not an artefact of contradictory
instructions.

## 4. The correction

### 4.1 Neutral rubric wording

Every outcome description now names the evidence types without prescribing sufficiency, and the
instrument carries this instruction explicitly:

> **No completeness threshold is prescribed. Some elements are reported as complete and some as
> partial; some have one terminal inverted repeat resolved and some both. Where you draw the
> line between sufficient and insufficient evidence is precisely what this audit measures, so
> record your genuine judgement rather than trying to match any rule.**

`biologically_indeterminate` is now described as *"evidence is present but you judge it too
fragmentary, ambiguous or incomplete to establish or exclude mobile context"* — an available
judgement rather than a prescribed outcome for any specific evidence pattern.

The three-state scoring rule from Amendment 003 is disclosed in the instrument, so the
adjudicator knows that choosing a general mobile outcome over a specific one costs nothing.
Disclosing that removes a second source of artefact without revealing any case's label.

Verified after regeneration: **zero** engine-threshold terms leak into either file. (`F3`
matched only spreadsheet cell references such as `<c r="F3">` and a Japanese font codepoint in
the theme; **zero cells contain it as visible text**.)

### 4.2 `NMV1A-001` excluded from the primary denominator

Its machine label was discussed with the adjudicator before review, so its blinding is lost.
That was my error, not the adjudicator's.

It is **not replaced**: drawing a substitute would require reading the sealed key to learn its
stratum, and no key will be opened before adjudication is returned. The case is flagged in both
instruments, the adjudicator still answers it, and it is excluded from the primary agreement
denominator.

**Primary denominator: 119.** All other frozen quantities are unchanged.

### 4.3 Prespecified partial-element evidentiary-threshold analysis

Declared now, before any decision is seen:

- The 19 exposed cases are analysed as a named secondary stratum.
- If the adjudicator judges them **mobile**, the finding is that **the engine is conservative at
  the partial-element boundary**. Cohort-wide, 106 blocks carry rule F3, so real mobile context
  may sit inside the indeterminate class.
- **Scope of any resulting revision is confined to the partial-element boundary.** Under the
  Amendment 003 scope-of-failure rule, this revises the F3 evidentiary threshold only. It does
  **not** invalidate class B, whose 333 blocks rest on complete elements with bilateral TIRs and
  complete ORFs — evidence unaffected by where the partial boundary is drawn.
- If disagreement is instead spread evenly across states, the scope-limitation does not apply
  and the ordinary gate governs.

This is a hypothesis registered in advance, not a result. Nothing has been computed.

## 5. Files

### Use these

| file | SHA-256 |
|---|---|
| `NMV1_AUDIT_BLINDED_120_R2.xlsx` | `dd4d32edb023c792493c93598a4263dce8ded3a3d4f1b901cdceb7b02c5c4bee` |
| `NMV1_AUDIT_CASEBOOK_R2.html` | `6eef49f7fea5cd19f1d9847b58883d045a4e3db74a16af36300c2259862d6b47` |

Verified identical to the originals across all 120 case rows and columns 1–20; 120 dropdown
validations preserved; 120 casebook panels and 120 SVG feature maps preserved. The only
additions are the rewritten rubric sheet, the rewritten casebook rubric block, and a
`scoring_note` column carrying the exclusion note on one row.

### Preserved unchanged

| file | SHA-256 |
|---|---|
| `NMV1_AUDIT_BLINDED_120.xlsx` | `f66dfea818dc5002cdeab64321ba5db28bcba6fcad7cc453cead79228f2759a4` |
| `NMV1_AUDIT_CASEBOOK.html` | `4e5531f6b37434757459bae46701e5c5ac5d429a7694d6b051b9838410299e55` |
| `NMV1_RULE_ENGINE_FROZEN.json` | `ed5db383bb0afe1a1a8433886d6666fe72c324975de99c6763a37824d51c2bee` |
| `NMV1_RULE_BASED_GROUND_TRUTH.tsv` | `1beecaa39048f4df52a3235f2dbc538056af9adc6912ee057cbac8ca55b85897` |
| `NMV1_AUDIT_MANIFEST.json` | `7473b4a5a0ec2c964376ba07936c1e9a1e3fccb2fe869459251011528143171f` |
| `NMV1_AUDIT_UNBLINDING_KEY.tsv` | `b7042b95365a939132f9c190093a6a871326095ef0cb629b90e8a6ef77639929` — **sealed** |
| `NMV1_AUDIT_SCORING_SPECIFICATION.json` | `55ab34bd7d1d6c66b4f13c56eedff72e5eb07e53b0ba744a2d176c4b4bdda704` |

## 6. What this amendment does not do

It does not change any rule-engine label, any evidence value, the 120-case selection, the
seven permitted outcomes, the three-state mapping, the success/revise/failure thresholds, or
any sealed key. It changes the wording of the instructions and removes one compromised case
from the scoring denominator.

## 7. Still prohibited

Sensitivity, specificity, PPV, NPV, balanced accuracy, NM-V1 gate verdicts, and opening any
unblinding key — until the completed adjudication is returned, frozen and hashed.
