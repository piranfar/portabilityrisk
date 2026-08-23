# NM-V1 Amendment 002 — rule-based ground truth with blinded expert audit

**Dated 2026-08-22, America/New_York (UTC−04:00).**

Replaces exhaustive manual review of 595 cases with a deterministic structural decision engine
applied to all 1,283 blocks, audited by a blinded expert review of **120** cases.

**No tool was rerun.** The frozen sample, the frozen rubric, the reference tool calls and the
complete V2 package are unchanged and preserved.

---

## 1. State at the time of this amendment

| | |
|---|---|
| adjudication begun | **no** |
| aggregate ground-truth performance computed | **none** |
| sensitivity / specificity / PPV / NPV / balanced accuracy | **not computed** |
| NM-V1 gate verdict | **not computed** |
| V1 unblinding key opened | **no — sealed** |
| V2 unblinding key opened | **no — sealed** |
| tool outputs altered | **none** |
| frozen sample altered | **none** |

## 2. Why 595 manual cases was the wrong instrument

One adjudicator judging 595 evidence-rich cases is exposed to fatigue-related classification
drift: decisions made late in a long session are not exchangeable with those made early, and
that drift would be indistinguishable from genuine biological disagreement. It is also
unnecessary. Most cases are decidable by structural criteria already written into the frozen
rubric — complete element, bilateral terminal inverted repeats, complete transposase ORF,
integrase plus attC — and applying those criteria by hand 595 times measures endurance, not
biology.

The expert judgement is better spent where it is decisive: checking whether the rules agree
with an expert on a stratified sample that deliberately includes the classes most likely to
fail.

## 3. Preserved, superseded, sealed

| artefact | SHA-256 | status |
|---|---|---|
| `NMV1_ADJUDICATION_BLINDED_PACKAGE_V2.xlsx` | `1a12a246256a14c8a2bb85fb665b03f7b0b4037d6596200fccac124ce0bbb8f5` | preserved, superseded as the review instrument |
| `NMV1_ADJUDICATION_CASEBOOK_V2.html` | `3cee4b2a8c682547c05b16ffb66c06c4acaf393fdc1eacfb1e2ddddb4637ab02` | preserved |
| `NMV1_ADJUDICATION_UNBLINDING_KEY_V2.tsv` | `54409bc5f704b2de29bca15199b2a1fc718ac888aea4640dcf7472c5f7cb459e` | **sealed**, preserved |
| `NMV1_QC_MANIFEST_V2.json` | `787685a17323c2e1bb1fed3fc42c4e512148227647fe8572fbf886b3e1e79ff1` | preserved |
| V1 package / key / rubric | `436e295c…` / `731fa74c…` / `a48494c5…` | **sealed**, preserved |

Nothing deleted or overwritten. The V2 package remains the fallback instrument if the audit
fails and exhaustive review becomes necessary.

## 4. The frozen decision engine

`NMV1_RULE_ENGINE_FROZEN.json` — `ed5db383bb0afe1a1a8433886d6666fe72c324975de99c6763a37824d51c2bee`

Frozen **before** application. Rules evaluated in fixed order, first match wins, so every block
receives exactly one label.

| order | rule | label | condition |
|---|---|---|---|
| 1 | **F1** | indeterminate | tool failure or missing output |
| 2 | **F2** | indeterminate | truncated or origin-wrapped block |
| 3 | **A** | multiple MGE evidence | complete IS with bilateral TIRs **and** complete integron with intI + attC |
| 4 | **B** | IS-associated | complete IS, complete transposase ORF, bilateral TIRs, no boundary problem |
| 5 | **C** | integron-associated | complete integron with integrase + attC, no tool failure or boundary problem |
| 6 | **D** | chromosomal-mobile | complete-level mobile evidence not exclusively assignable |
| 7 | **F3** | indeterminate | partial-only element or incomplete integron |
| 8 | **E** | chromosomal-quiescent | no reference evidence, both tools completed, evaluable, no boundary warning |
| 9 | **F4** | indeterminate | conflicting evidence the rules cannot resolve |

**The original HMM annotation is not an input to any rule.** Using the method under test to
build the truth it is judged against would be circular. The HMM label is read once, to define
the `S_DISC` audit stratum by selection, and is never shown to the adjudicator. This is
verified two ways in §6.

## 5. Ground truth over all 1,283 blocks

| rule | label | n |
|---|---|---:|
| A | multiple MGE evidence supported | 154 |
| B | IS-associated supported | 333 |
| C | integron-associated supported | 8 |
| D | chromosomal-mobile supported | 1 |
| E | chromosomal-quiescent supported | 673 |
| F1 | indeterminate — tool failure | 1 |
| F2 | indeterminate — boundary | 7 |
| F3 | indeterminate — partial / incomplete | 106 |
| F4 | indeterminate — unresolved | 0 |
| | **total** | **1,283** |

**Evaluable 1,169 · indeterminate 114 (8.9 %).**

`NMV1_RULE_BASED_GROUND_TRUTH.tsv` — `1beecaa39048f4df52a3235f2dbc538056af9adc6912ee057cbac8ca55b85897`

Each row carries the label, rule id, the evidence fields that triggered it, and evaluable
status.

**Worth stating now, because it shapes what NM-V1 can conclude:** only **8** blocks are
exclusively integron-associated and only **1** is class D. Integron evidence in this cohort
almost always co-occurs with a complete IS and is therefore absorbed into rule A. The
integron-specific arm is consequently thin, and no strong class-specific claim about integrons
should be expected from this design.

## 6. Independent rule verification — **PASS, 27 checks, zero disagreements**

`NMV1_RULE_VERIFICATION_REPORT.txt`

A separate verifier, importing nothing from the application script, re-derives every predicate
from the raw hit tables and re-applies the rule order as an explicit ordered predicate table
rather than an if/elif chain, so a control-flow error in one implementation would not be
reproduced by the other.

- rule assignments disagreeing: **0 of 1,283**
- blocks with no matching rule: **0** (totality holds)
- per-rule label consistency: all nine rules **MATCH**
- indeterminate == F1+F2+F3+F4: **114 = 114**
- **HMM independence, tested two ways:** flipping every HMM field changed **0** labels, and
  the string `hmm` occurs **0** times inside the classifier function
- audit package: 120 cases, 120 unique tokens, 120 unique blocks, all six strata reconcile,
  key labels match ground truth, cap respected

The verifier also reports that 36 blocks satisfy more than one rule condition. That is expected
and not a defect: the engine declares first-match-wins, so the **assignment** is unique even
where raw conditions overlap. Stating it rather than hiding it is the point.

## 7. Blinded expert audit — 120 cases

Deterministic selection: candidates sorted by `SHA-256("<block_id>|20260821|<stratum>")`,
strata processed in fixed order so selections are mutually exclusive.

| stratum | definition | pool | target | allocated |
|---|---|---:|---:|---:|
| `S_IS` | rule B | 333 | 20 | **23** |
| `S_INT` | rule C | **8** | 20 | **8** |
| `S_MULTI` | rule A | 154 | 20 | **21** |
| `S_QUIET` | rule E | 673 | 20 | **25** |
| `S_DISC` | HMM disagrees with the rule reading, rule label not indeterminate | 274 | 20 | **22** |
| `S_INDET` | rules F1–F4 | 114 | 20 | **21** |
| | | | | **120** |

`S_INT` holds only 8 eligible cases, so all 8 were taken and the 12-case shortfall was
reallocated proportionally across the other five strata **before any human decision was seen**,
exactly as the frozen shortfall rule requires. Pools, rule and per-stratum token hashes are in
`NMV1_AUDIT_MANIFEST.json` — `7473b4a5a0ec2c964376ba07936c1e9a1e3fccb2fe869459251011528143171f`.

**Withheld from the adjudicator:** species, gene identity, BioProject, accession, original HMM
class, tool names, the rule-based label, the rule id, sampling stratum, audit stratum, and all
aggregate results. Automated leakage check over 255 identity strings: **0 leaks in the
casebook, 0 in the spreadsheet.**

## 8. Audit gate — frozen before review

Primary metric: overall agreement between the adjudicator's outcome and the rule-based label,
Wilson 95 % CI.

| verdict | condition |
|---|---|
| **SUCCESS** | agreement ≥ 0.90 **and** lower CI ≥ 0.80 |
| **REVISE** | agreement 0.80–0.899, **or** lower CI 0.65–0.799 |
| **FAIL** | agreement < 0.80, **or** systematic disagreement within a major class |

Class-specific agreement is reported alongside. **A failure confined to one evidence class
revises that class only and does not invalidate unaffected classes.** Indeterminate cases are
scored for agreement on *indeterminate status itself* and are never forced into a mobile or
quiescent label.

## 9. Still prohibited until the audit is returned and frozen

Sensitivity, specificity, PPV, NPV, balanced accuracy, NM-V1 gate verdicts, and opening either
unblinding key.
