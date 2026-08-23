# NM-V4C amendment 001 — co-primary gate clarification

**Dated 2026-08-21, America/New_York (UTC−04:00).**

This amendment corrects the **gate ledger** of NM-V4C. It withdraws no result, alters no
computation, and rewrites nothing.

| artefact | SHA-256 | status |
|---|---|---|
| `NMV4C_FROZEN_DESIGN.json` | `8a3c76b157cbf2cd5279342cb2752e1974a2baac976b62455ea0c66f3de4495d` | **unchanged** |
| `NMV4C_RESULT_RECEIPT.json` | `2ac346f7cc3be43bfa4b2bc92791e1428333693a094064ea64339fcba4cc1679` | **unchanged — superseded ledger only** |
| `nmv4c_family_host_vehicle.tsv` | `f13de8f47bc67c3170ce1d54b40b370fb1a943b8cd9bfd9c54e7a4b83af0fd98` | **unchanged** |
| `nmv4c_host_vehicle_summary.tsv` | `dab40b18db26b062076004586f9d5daef96aaee55f9b0bc9a1fca62db04c3903` | **unchanged** |

**Supersession pointer.** The `gates` and `gates_passed` fields inside
`NMV4C_RESULT_RECEIPT.json` are superseded by §5 of this file. Every numeric result in that
receipt stands. The receipt is not edited, because a receipt records what a run computed and
editing it would destroy that record.

---

## 1. The ambiguity in gate G5

G5 was frozen as:

> "the MH odds ratio for A. baumannii against Pseudomonas aeruginosa excludes 1 in the same
> direction, so the effect is not a generic low-plasmid artefact"

The frozen design separately declared **P1-full** and **P2-one-genome-per-BioProject** to be
**co-primary**, with the rule that disagreement means the result "is reported as
sampling-structure dependent and the universal claim is not made."

**G5 did not state whether co-primary agreement was required for every contrast or only for the
principal *A. baumannii*–*Klebsiella* contrast.** That ambiguity is the defect this amendment
records. It was in the design, not in the execution.

## 2. What the implementation did

`nmv4c_score.py` evaluated G5 against **P1-full only**:

```
"G5_pa_control": bool(abpa and abpa["excludes_1"] and abpa["or"] > 1)
```

where `abpa` is the P1-full result. The P2 estimator was computed, exported and reported, but
was not consulted by the gate expression.

## 3. What the two estimators showed

| contrast | P1-full | P2 BioProject-balanced | co-primary |
|---|---:|---:|---|
| *A. baumannii* vs *Klebsiella* | 50.29 [45.61, 55.45] | 21.95 [18.92, 25.47] | **agree** |
| *P. aeruginosa* vs *Klebsiella* | 28.94 [25.79, 32.46] | 16.06 [13.79, 18.70] | **agree** |
| ***A. baumannii* vs *P. aeruginosa*** | 1.26 [1.11, 1.43] | **1.18 [0.99, 1.41]** | **DISAGREE** |

P1 excluded 1. P2, which balances BioProjects, did not.

## 4. Why an unqualified "7 of 7 gates" is incorrect

Because the co-primary rule was violated for that one contrast, and a gate ledger that reports
7/7 without saying so implies a robustness the data do not carry for that comparison.
Reporting the count alone would have been misleading even though every individual number was
correct.

## 5. Corrected gate status — this table supersedes the receipt's `gates` block

| gate | subject | status |
|---|---|---|
| G1 | primary *A. baumannii*–*Klebsiella* host effect | **PASS** |
| G2 | family conditioning, ≥10 families, weight concentration ≤30 % | **PASS** (56 of 58 families; max weight share 15.5 %) |
| G3 | BioProject balancing on the primary contrast | **PASS** (OR 21.95, CI excludes 1) |
| — | cluster permutation of host labels within family | **PASS** (0 of 2000, p = 0.0005) |
| G4 | leave-one-family-out influence | **PASS** (max 2.68 %, threshold 20 %) |
| G7 | leave-one-BioProject-out influence | **PASS** (max 2.02 %, threshold 20 %) |
| G6 | five-class and focused-route coherence, primary contrast | **PASS** |
| **G5** | ***A. baumannii*–*P. aeruginosa* focused B-versus-plasmid separation** | **QUALIFIED / NOT ROBUST ACROSS CO-PRIMARY ESTIMATORS** |

## 6. What changes, and what does not

**Nothing is withdrawn.** All 58 family-level 2×2 tables, all odds ratios, the permutation
result and both influence analyses stand exactly as computed and independently verified.

**The interpretation changes.** The focused B-versus-plasmid contrast was built as a
*two-route* model — chromosomal-mobile against plasmid. That model cannot represent what the
data show, because it conditions away class A, and class A is precisely where *A. baumannii*
and *P. aeruginosa* diverge:

| host | A | B | C | D | E |
|---|---:|---:|---:|---:|---:|
| *A. baumannii* | 16.4 % | **69.9 %** | 3.7 % | 3.6 % | 6.4 % |
| *P. aeruginosa* | **50.5 %** | 36.2 % | 11.3 % | 1.1 % | 1.0 % |
| *Klebsiella* | 19.9 % | 11.3 % | 9.6 % | 9.4 % | **49.8 %** |

The failure of G5 is therefore **informative rather than adverse**: the two-route model is
underspecified, and a **three-architecture** model is required. That model is frozen and
tested in `NMV4D_FROZEN_THREE_ARCHITECTURE_DESIGN.json`, which exists because of this
amendment.

## 7. Terminology, binding on all downstream text

**Preferred central formulation:**

> Identical resistance determinants are routed into distinct plasmid-conjugative,
> chromosomal-mobile and chromosomal-quiescent architectures across bacterial hosts.

**Permitted formulation:**

> Portability architecture is a host-associated property of the determinant–host combination.

**"Host determines"** may be used descriptively **only** if immediately defined as statistical
association and not causal control. Preferred verbs are *routed into*, *segregated into*,
*associated with*. The word *routed* describes where determinants are **found**, not a process
that was observed.

**Still prohibited:** observed horizontal transfer, causal host control, demonstrated
conjugation, transfer rate, clinical risk.

## 8. Standing conditionality

Every class-B conclusion in NM-V4C and NM-V4D — including the 69.9 % chromosomal-mobile
fraction in *A. baumannii* — remains **conditional on NM-V1 validation of the MGE layer**,
which has not been run.
