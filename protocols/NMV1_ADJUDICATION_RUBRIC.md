# NM-V1 blinded adjudication rubric

**Adjudicator:** Vahhab Piranfar

Frozen rules SHA-256 `e454873a89d1d56b20adb9ae157f20224076966daa0ed34bcbe44318063260f6`.

## What you are looking at

Each row is one genomic block. Three methods are shown, de-identified as
**Method X**, **Method Y** and **Method Z**. You are not told which is which,
nor the species, study, gene identity or original classification.

| column | meaning |
|---|---|
| `methodX_*` | marker counts from a protein-homology method |
| `methodY_*` | element counts from a boundary-based method, with complete vs partial and inverted-repeat evidence |
| `methodZ_*` | integron structures: complete, CALIN, In0 |

## Permitted outcomes

Enter exactly one in the `adjudication` column:

- `chromosomal_mobile_supported`
- `chromosomal_quiescent_supported`
- `integron_associated_supported`
- `IS_associated_supported`
- `multiple_MGE_evidence_supported`
- `neither_classification_supported`
- `biologically_indeterminate`

## Decision rules

**IS_associated_supported** — a transposase ORF with at least one resolved terminal inverted repeat, or a complete element reported by a boundary-based method

**integron_associated_supported** — an integrase with at least one attC site, or a complete integron structure

**multiple_MGE_evidence_supported** — both IS and integron evidence present in the block

**chromosomal_mobile_supported** — any of the three above; the block sits in mobile context

**chromosomal_quiescent_supported** — no credible mobile-element evidence within the block

**neither_classification_supported** — evidence is present but contradicts both readings

**biologically_indeterminate** — the evidence shown cannot decide; USE THIS RATHER THAN GUESSING. An honest indeterminate is more useful than a forced call and is reported as its own category.

**boundary_cases** — if an element appears truncated at a block edge, judge only what is visible and mark indeterminate if the visible part is insufficient

## Instruction

Judge only the structural evidence shown. Do not infer from gene identity, organism or study; none is provided.

Put a short free-text justification in `reason`. Leave no row blank.
