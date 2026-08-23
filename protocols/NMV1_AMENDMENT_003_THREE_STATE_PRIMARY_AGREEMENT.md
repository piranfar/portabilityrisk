# NM-V1 Amendment 003 — three-state architecture as the primary agreement metric

**Dated 2026-08-22, America/New_York (UTC−04:00). Issued before any adjudication decision was
entered.**

Corrects **what the primary agreement metric measures**. It does not change the 120-case
sample, the rule-engine labels, the evidence, the sealed keys, or the success thresholds.

---

## 1. State at the time of this amendment

| | |
|---|---|
| adjudication decisions entered | **none** |
| agreement computed | **none** |
| sensitivity / specificity / PPV / NPV / balanced accuracy | **not computed** |
| NM-V1 gate verdict | **not computed** |
| unblinding keys opened | **none — V1, V2 and audit keys all sealed** |
| tool rerun | **none** |

## 2. The defect

`NMV1_RULE_ENGINE_FROZEN.json` → `audit_gate.primary_metric` reads, verbatim:

> "overall agreement between the adjudicator's outcome and the rule-based label across the
> audited cases, with a 95 per cent CI"

Both sides of that comparison are **seven-label** values. As frozen, the primary SUCCESS /
REVISE / FAIL gate is therefore driven by **exact seven-label agreement**. A search of
`NMV1_RULE_ENGINE_FROZEN.json`, `NMV1_GROUNDTRUTH_RULES_FROZEN.json` and
`NMV1_FROZEN_DESIGN.json` returns **zero** occurrences of any mapping from the seven outcomes
onto a collapsed architecture. The only partial collapse present is
`indeterminate_handling`, which covers indeterminate status alone.

### Why this biases the verdict downward

NM-V1 validates the **chromosomal-mobile versus chromosomal-quiescent** classification. That is
the quantity class B rests on and the quantity the manuscript claims.

The rule engine assigns subtype labels by strict structural criteria: rule B requires bilateral
terminal inverted repeats *and* a complete transposase ORF before it will say
`IS_associated_supported`. A biologist looking at the same panel may reasonably record the more
general `chromosomal_mobile_supported` — declining to commit to a subtype is good practice, not
error. Under exact seven-label scoring that is counted as a disagreement.

The failure mode is systematic rather than random: it would depress agreement uniformly across
the mobile strata, which are 52 of the 120 audited cases. A run in which the adjudicator agreed
with the architecture on every single case could still be scored REVISE, or conceivably FAIL,
purely on subtype granularity. That would be a measurement artefact reported as a biological
result.

## 3. The correction

### 3.1 Primary metric — three-state architecture

The seven permitted outcomes map onto three states. This mapping is **total and mutually
exclusive**: every outcome belongs to exactly one state.

| state | outcomes |
|---|---|
| **MOBILE** | `multiple_MGE_evidence_supported`, `IS_associated_supported`, `integron_associated_supported`, `chromosomal_mobile_supported` |
| **QUIESCENT** | `chromosomal_quiescent_supported` |
| **NON-EVALUABLE** | `neither_classification_supported`, `biologically_indeterminate` |

**Primary agreement** = the proportion of audited cases where the adjudicator's outcome and the
rule-based label fall in the **same state**, with a Wilson 95 % confidence interval.

### 3.2 Thresholds — unchanged

Carried over exactly as frozen; only the metric they are applied to changes.

| verdict | condition |
|---|---|
| **SUCCESS** | agreement ≥ 0.90 **and** Wilson lower CI ≥ 0.80 |
| **REVISE** | agreement 0.80–0.899, **or** lower CI 0.65–0.799 |
| **FAIL** | agreement < 0.80, **or** systematic disagreement within a major state |

### 3.3 Secondary — subtype resolution

Exact seven-label agreement is retained and reported as a **secondary subtype-resolution
analysis**. It carries no gate and cannot change the NM-V1 verdict.

A human call of `chromosomal_mobile_supported` against a rule-engine call of
`IS_associated_supported` is a subtype-resolution difference, **not** a primary biological
disagreement: both assert the same mobile architecture.

### 3.4 Integron cases — descriptive only, explicitly not a gate

Only **8** of the 120 audited cases are exclusively integron-associated (rule C), because
integron evidence in this cohort almost always co-occurs with a complete IS and is absorbed
into rule A.

Eight cases cannot support a class-specific success gate. Their agreement is reported
**descriptively** with an **exact Clopper–Pearson** 95 % interval — exact rather than Wilson,
because at n = 8 the normal-approximation family is unreliable — and is explicitly barred from
determining the overall NM-V1 verdict. With n = 8 the interval will be wide whatever the
result, and that width is the finding: this design cannot resolve integron-specific
performance.

### 3.5 State-specific agreement

Agreement within MOBILE, QUIESCENT and NON-EVALUABLE is reported separately. The frozen
scope-of-failure rule is unchanged and now applies at state level: **a failure confined to one
state revises that state only and does not invalidate unaffected states.**

Indeterminate handling is unchanged: NON-EVALUABLE cases are scored for agreement on
non-evaluable status itself and are never forced into MOBILE or QUIESCENT.

## 4. Audit composition under the corrected metric

Derived from the rule-id counts already reported for the 120-case audit set. **No key was
opened and no adjudication decision exists.**

| state | rules | n of 120 |
|---|---|---:|
| MOBILE | A 21 + B 23 + C 8 + D 0 | **52** |
| QUIESCENT | E 47 | **47** |
| NON-EVALUABLE | F2 2 + F3 19 | **21** |
| | | **120** |

All three states carry enough cases to be scored. MOBILE and QUIESCENT are well represented at
52 and 47; NON-EVALUABLE at 21 is adequate for a state-specific estimate but not for fine
subdivision, and will be reported with that limit stated.

## 5. What is explicitly unchanged

The 120-case sample and its selection; every rule-engine label and rule id; all evidence,
feature maps and casebook content; the seven permitted adjudication outcomes and the rubric
shown to the adjudicator; the success, revise and failure thresholds; all sealed unblinding
keys; every prior hash.

| artefact | SHA-256 | status |
|---|---|---|
| `NMV1_RULE_ENGINE_FROZEN.json` | `ed5db383bb0afe1a1a8433886d6666fe72c324975de99c6763a37824d51c2bee` | unchanged; `audit_gate.primary_metric` superseded by §3.1 of this file |
| `NMV1_RULE_BASED_GROUND_TRUTH.tsv` | `1beecaa39048f4df52a3235f2dbc538056af9adc6912ee057cbac8ca55b85897` | unchanged |
| `NMV1_AUDIT_BLINDED_120.xlsx` | `f66dfea818dc5002cdeab64321ba5db28bcba6fcad7cc453cead79228f2759a4` | unchanged |
| `NMV1_AUDIT_CASEBOOK.html` | `4e5531f6b37434757459bae46701e5c5ac5d429a7694d6b051b9838410299e55` | unchanged |
| `NMV1_AUDIT_MANIFEST.json` | `7473b4a5a0ec2c964376ba07936c1e9a1e3fccb2fe869459251011528143171f` | unchanged |
| `NMV1_AUDIT_UNBLINDING_KEY.tsv` | `b7042b95365a939132f9c190093a6a871326095ef0cb629b90e8a6ef77639929` | **sealed**, unchanged |

The adjudicator still records one of the seven detailed outcomes. The collapse to three states
happens only at scoring time, so no information is discarded and the secondary subtype analysis
remains possible.

## 6. Still prohibited

Sensitivity, specificity, PPV, NPV, balanced accuracy, NM-V1 gate verdicts, and opening any
unblinding key — until the completed adjudication is returned, frozen and hashed.
