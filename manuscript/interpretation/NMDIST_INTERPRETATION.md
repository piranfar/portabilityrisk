# NM-DIST — distance-to-MGE analysis: manuscript interpretation

**Dated 2026-08-22.** Executed under `NMDIST_FROZEN_PROTOCOL_V1.json` =
`226c0691cbd6ac9750fd8e816192420b3c0d40cadd4bb2ad05c7ca0089593a26`, frozen before any
species-specific distance outcome was computed. No frozen cohort, portability class, denominator
or NM result was altered. Independent verification: **32 checks, 0 disagreements.**

---

## 1. What was measured

All **35,140** PRIMARY chromosomal acquired-ARG occurrences across **21,955** frozen context
blocks. For each occurrence, the distance to the nearest qualifying MGE marker **within its own
±10 kb window**; overlap scored as 0; topology-aware on circular replicons. **18,837 occurrences
had no qualifying marker in-window and are right-censored at 10,000 bp** — never assigned a
distance beyond 10 kb. **±20 kb remains NOT EVALUABLE.**

**The 111 documented exceptions.** 111 occurrences sit in an MGE-positive *block* whose nearest
in-block feature lies **10,231–16,931 bp** away, i.e. outside their own window. Blocks are merged
shared-context regions and can span up to 62,349 bp, so block-positive is not window-positive.
These 111 are **censored at 10 kb**; using the out-of-window feature would manufacture a distance
the design declares unavailable. Confirming the arithmetic: 16,414 block-positive − 111 = **16,303
window-positive = frozen class B exactly**.

Primary estimand: **block-balanced** (each block contributes total weight 1; occurrence weight
1/m). Weights sum to **21,955.000000**. Every block nests within exactly one BioProject, so
recomputing weights inside a bootstrap resample returns the same value — verified.

---

## 2. Results

### Primary, block-balanced, nearest MGE of either type

| group | n | F(1 kb) | F(2 kb) | F(5 kb) | F(10 kb) | RMD (bp) | median |
|---|---:|---:|---:|---:|---:|---:|---|
| ***A. baumannii*** | 8,005 | **0.5224** | **0.5878** | **0.5964** | **0.6323** | **4,180** | **647 bp** |
| *P. aeruginosa* | 7,150 | 0.0993 | 0.1256 | 0.1410 | 0.1717 | 8,595 | not reached |
| *Klebsiella* group | 15,568 | 0.1039 | 0.1227 | 0.1533 | 0.2390 | 8,365 | not reached |

### Primary contrasts, BioProject-bootstrap (B = 2,000, seed 20260822), Holm-corrected

| contrast | 1 kb | 2 kb | 5 kb | 10 kb |
|---|---:|---:|---:|---:|
| **P1** *A. baumannii* − *Klebsiella* | **+0.4186** | **+0.4651** | **+0.4430** | +0.3932 |
| **P2** *A. baumannii* − *P. aeruginosa* | **+0.4231** | **+0.4622** | **+0.4554** | **+0.4606** |

All eight Holm-corrected **p = 0.001** (the bootstrap resolution floor). Restricted-mean distance:
P1 **−4,185 bp** [−4,588, −3,696]; P2 **−4,415 bp** [−4,805, −3,953].

---

## 3. The six required questions

### Does *A. baumannii* remain exceptional across the full spatial scale?
**Yes, and the excess is not a 10 kb artefact.** The contrast is present at every landmark and is
**largest at 2 kb (+0.4651), not at 10 kb (+0.3932)**. Against *P. aeruginosa* it is flat across
the whole window (+0.42 to +0.46). Panel **a** shows the *A. baumannii* curve rising almost
vertically and plateauing by ~2 kb, while both comparators climb slowly across the full window.

### At what distances does separation emerge?
**Immediately, and it is essentially complete by 2 kb.** *A. baumannii* reaches
0.5878 of its eventual 0.6323 by 2 kb — **93 % of all its detections occur within 2 kb**. For
*Klebsiella* the same ratio is 0.1227/0.2390 = **51 %**. Median nearest-marker distance is **647 bp**
for *A. baumannii* and **not reached** (>10 kb) for both comparators. Under the any-ARG block
distance (S3) the *A. baumannii* median falls to **279 bp**.

### Does the contrast survive block balancing?
**Yes.** Block balancing lowers every absolute level — *A. baumannii* F(10 kb) falls from 0.8091
occurrence-weighted to 0.6323 block-balanced — but the **contrast is barely affected**: P1 at 1 kb
is +0.4671 occurrence-weighted and **+0.4186 block-balanced**. It also survives the deterministic
one-occurrence-per-block sensitivity (**S2: +0.4209**), truncated-block exclusion (S4: +0.4187),
circular-wrapped exclusion (S5: +0.4191), and BioProject resampling. Leave-one-BioProject-out
across all 2,248 projects moves P1 at 1 kb by at most **0.0140** on a baseline of 0.4186.

### Is the result driven primarily by IS/transposase?
**Entirely.** The IS/transposase-only endpoint (S6) reproduces the primary result almost exactly
(P1 at 1 kb **+0.4186**; RMD 4,181 against 4,180). The integrase/integron-only endpoint shows **no
*A. baumannii* excess at all**: +0.0100 at 1 kb and **−0.0136 at 10 kb**, with all three groups
between 0.099 and 0.116 at 10 kb. Panel **c** makes this unmistakable. The integrase arm is
secondary by design (3,033 of 32,364 features) and NM-V1 did not independently validate every
NON-EVALUABLE boundary.

### Does this strengthen, qualify or leave unchanged the 46.39 % headline?
**It leaves the number unchanged and strengthens its interpretation, with one qualification.**

- **Unchanged.** 46.39 % is the occurrence-weighted ±10 kb figure over all species. No frozen
  artefact was touched; the analysis reproduces the frozen distance column with **0 mismatches**,
  and occurrence-weighted *A. baumannii* F(10 kb) = **0.8091**, matching the published 80.9 %.
- **Strengthened.** The ±10 kb threshold is **not load-bearing**. Signal is concentrated far inside
  the window, so the headline is not an artefact of where the boundary was drawn — the most
  obvious methodological objection to a fixed window is now answered with data.
- **Qualified.** Block balancing reduces absolute levels substantially (*A. baumannii* 80.9 % →
  63.2 %), consistent with the already-disclosed multiplicity effect behind 46.39 % versus
  30.14 %. The **contrast** is robust; the **absolute proportion** depends on the denominator, and
  both must continue to be reported together.

### Which manuscript claim and figure should change?
1. **C05/C06 gain a spatial sentence.** Not "within 10 kb" alone, but that the association is
   concentrated at short range — median 647 bp in *A. baumannii*, 93 % of detections within 2 kb.
2. **P4 (the *A. baumannii* claim) strengthens in kind, not just degree.** It is no longer only
   "more often near an MGE" but "**far nearer**": median 647 bp against not-reached in both
   comparators, a restricted-mean difference of ~4.2–4.4 kb.
3. **Figure 3** (chromosomal MGE, both denominators) should gain panel **a** of this figure, or
   this figure should enter as a new main figure. Panel **d** replaces any separate weighting
   supplementary.
4. **The IS-versus-integrase decomposition is new and should be stated explicitly** — the
   chromosomal signal is an IS/transposase phenomenon, with no integron contribution to the
   species contrast.
5. **No claim changes direction, and nothing is withdrawn.**

---

## 4. Wording constraints observed

No claim of intact transposons, demonstrated mobilization, horizontal transfer, conjugation or
phenotype. No distance beyond 10 kb is estimated or implied. Marker-negative observations are
reported as **right-censored at 10 kb**, never imputed. Permitted phrasing: *"the nearest detected
MGE marker lies within X bp"*, *"marker-negative within the ±10 kb window"*, *"concentrated at
short range"*.

## 5. Provenance

| artefact | SHA-256 |
|---|---|
| `NMDIST_FROZEN_PROTOCOL_V1.json` | `226c0691cbd6ac9750fd8e816192420b3c0d40cadd4bb2ad05c7ca0089593a26` |
| `nmdist_occurrence_block_distances.tsv` | `8d06ee836b5c7d0c05d3a58ba1f854005982d43a386bb9b9e8494c691a36916c` |
| `nmdist_primary_estimates.tsv` | `905ee498db419347e74fca453a5ac006baacaf818f54772dc4b0a5f3907af683` |
| `nmdist_species_contrasts.tsv` | `293fc60e3d5c4ed51933649d16617eefcc19b57c5e51c6cf3b2b490bce45e47c` |
| `nmdist_sensitivity_results.tsv` | `aaf2b41494e00124a720bc1d8956f8916af74bbeaf36b79f9a71593257e0478f` |
| `NMDIST_RESULT_RECEIPT_V1.json` | `a335adcafaf38027df3834948472864508599929e52f27c6b2b6518eb1ac3549` |

**Figure caption (draft).** Distance from chromosomal acquired resistance-gene occurrences to the
nearest mobile-element marker. **a**, Weighted cumulative detection *F*(*d*), block-balanced;
shaded bands are BioProject-bootstrap 95 % intervals (B = 2,000); marker-negative occurrences are
right-censored at 10 kb. **b**, Differences in *F*(*d*) at the four prespecified landmarks with
bootstrap 95 % intervals; Holm-corrected across the two primary contrasts at each landmark.
**c**, Decomposition by marker type at 10 kb. **d**, Occurrence-weighted versus block-balanced
weighting at 10 kb. *n* = 8,005 / 7,150 / 15,568 chromosomal occurrences.
