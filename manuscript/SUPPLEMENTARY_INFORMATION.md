# Supplementary Information

**Version 6.** Supersedes version 5; see `PORTABILITYRISK_CHANGELOG_V5_TO_V6.md`.

**Replicon-resolved portability of 74,349 acquired resistance-gene occurrences across
6,288 closed Gram-negative ESKAPE genomes**

Vahhab Piranfar

---

# Supplementary Methods

## Supplementary Method 1 | Blinded internal adjudication of the chromosomal mobile-element layer

**Design.** 120 chromosomal context blocks were drawn from a frozen 1,283-block sampling frame under
a stratified design registered before any block was drawn. The four strata pool 580 of those
1,283 blocks, not all of them: the frame is every block in the sampling universe, while the
strata contain only the blocks eligible for the state they represent. Strata and allocations were: M1 (106 available, 25 drawn), M2 (157, 34), M4 (1, 1) and Q (316, 60). The frozen protocol states the intended balance as MOBILE 60 (M1 25 + M2 34 + M4 1) against QUIESCENT 60 (Q 60), chosen for equal power in the two states that can be sampled; **these are rule-engine states, not portability classes** — earlier revisions wrote A, B, D and E here, which read as sampling the no-marker chromosomal class into MOBILE and the conjugative plasmid class into QUIESCENT; the delivered package contains exactly MOBILE 60, QUIESCENT 60,
NON_EVALUABLE 0.

**Weighting.** Strata are sampled at different fractions, so estimates are design-weighted by the
inverse probability of selection, weight = pool size / allocated: A 4.2400, B 4.617647, D 1.0000, E
5.266667.

**Blinding.** Methods were de-identified as X, Y and Z. Species, study, gene identity and any prior
classification were withheld. The scored case file was frozen and hashed before the key was opened.

**Adjudicator.** One adjudicator, Vahhab Piranfar, who is also the repository owner and this
manuscript's sole author. The registered eligibility requirement was "at least one person not
involved in building the pipeline"; the recorded eligibility basis is that the adjudicator "did not
develop the original MGE annotation pipeline or the NM-V1 scoring code, and has not inspected
block-level validation outcomes". No professional credential, affiliation or years of experience is
recorded, no second reader was used, and no inter-rater statistic exists for this module. These are
stated because they bound what the audit can support: it establishes that a blinded human applying
the frozen rubric reproduces the engine's state, not that two independent experts agree.

**Metric and gate.** The primary metric is design-weighted agreement between the adjudicator's
three-state call and the rule-engine state, under the mapping MOBILE = {M1, M2, M3, M4},
QUIESCENT = {Q}, NON_EVALUABLE = {NE1–NE4}. **These are rule-engine states, not portability
classes.** Earlier revisions labelled them A–F, which collided with the five portability classes
A–E and read as asserting that class A is mobile and class E quiescent — close to the opposite
of what this paper reports. The registered gate was SUCCESS at agreement ≥ 0.90 **and** bootstrap
lower bound ≥ 0.80; REVISE at 0.80–0.899 or lower bound 0.65–0.799; FAIL below 0.80 or on systematic
disagreement within a sampled state. These thresholds were carried over unchanged from the original
frozen gate, registered before any adjudication decision was made.

**Interval.** The confidence interval is a stratified bootstrap interval taken from the 2.5th and
97.5th order statistics — a percentile interval, although the NM-V1C record does not use that word.

**Result.** Confusion matrix MOBILE→MOBILE 59, MOBILE→QUIESCENT 1, QUIESCENT→QUIESCENT 60; raw
agreement 119 of 120; design-weighted agreement 0.9920 (95% CI 0.9761–1.0000); gate verdict SUCCESS.
The adjudicator left zero notes and selected NON_EVALUABLE zero times.

**The single disagreement.** Case NMC-096, stratum B: engine MOBILE, adjudicator QUIESCENT. The
panel showed one complete element at 17,922–19,125 with bilateral 20 bp inverted repeats and an
internal open reading frame, in a 20,858 bp circular block with no boundary problem and normal tool
completion.

**Scope limitation, binding.** This audit establishes rubric consistency for the MOBILE versus QUIESCENT discrimination only, and is not a validation: the adjudicator built neither the pipeline nor the scoring code but is the author. It cannot establish consistency for the NON_EVALUABLE state and cannot produce an integron-specific estimate, because
the unused population contained zero NON_EVALUABLE and zero integron-exclusive blocks. Any sentence
citing the 0.9920 figure carries that boundary.

**The superseded audit.** An earlier round of this audit scored 62/120 and returned a registered
FAIL. That result stands as a FAIL in the record and is not withdrawn. Its cause is documented: the
instrument displayed Method X — the profile-homology annotation path **under test** — to the
adjudicator, while the frozen ground-truth rule engine was defined exclusively on the two structural
methods Y and Z and was therefore prohibited from using Method X at all. The adjudicator was shown
evidence the engine could not use and asked to agree with it, creating an asymmetric comparison. The
contamination was deterministic rather than probabilistic: of 47 relevant cases, the adjudicator
called QUIESCENT in 0 of 38 where a Method X marker was displayed and in 9 of 9 where it was not.

The corrected confirmatory audit removed Method X from the instrument in every form, used a fresh
draw with **zero** reuse of any previous adjudication (machine-checked, zero overlap), and retained
the same rule engine and the same registered gate. The 62/120 result may not be quoted as the
engine's biological accuracy, its adjudications may not be reused, and it may not be rescored under
a changed rule; those prohibitions are recorded in the amendment that diagnosed it.


## Supplementary Method 2 | Terminology

These nine terms are used in one sense each throughout the manuscript and this Supplement. They are
not interchangeable, and several pairs that look similar are separated by an evidence boundary.

| term | means | does not mean |
|---|---|---|
| **documented replicon location** | the molecule type NCBI records for that replicon in the assembly's sequence report, read from `assigned_molecule_location_type` | a prediction from a plasmid-classification tool; nothing here is predicted |
| **mobile-element marker** | any insertion-sequence, transposase, integrase or integron feature in the frozen marker inventory | a complete or intact element |
| **transposase / integrase homology marker** | a profile-homology hit to a transposase or integrase family | that a functional protein is encoded, or that an element boundary exists |
| **structurally complete insertion sequence** | an element reported by ISEScan as type `c`, with a complete transposase reading frame and bilateral resolved terminal inverted repeats | that the element has moved, moves, or can move |
| **chromosomal mobile-element-associated context** | a chromosomal occurrence with ≥1 marker within its own ±10 kb window (class B) | that the determinant arrived by transposition |
| **mobility-marker-negative plasmid** | a plasmid on which MOB-suite detected no relaxase and no mating-pair formation system (class C) | non-mobilizable; it is a statement about a marker database |
| **mobilization-consistent plasmid** | a plasmid carrying a detected relaxase but no mating-pair formation system (class D) | that mobilization has occurred |
| **conjugation-consistent plasmid** | a plasmid carrying both a relaxase and a mating-pair formation system (class E) | a demonstrated conjugative plasmid |
| **observed transfer** | — | **NOT MEASURED.** No mating, transposition, transfer or element-activity assay exists anywhere in this study |

Two further distinctions carry the same weight. **Occurrence-weighted** and **block-weighted** are
different denominators over the same data and give different absolute levels (46.39% against
30.14%); the species ordering is preserved under both. **Homology endpoint** and **structural
endpoint** are two stringencies over the same 21,955 blocks, not two independent lines of evidence.

A third weighting is used and is easy to miss: **block-balanced** gives each occurrence weight
1/*m* inside a block of *m*, so weights sum to the occurrence count, while **block-weighted** takes
the block itself as the unit and they sum to the block count. Their values are close enough to be
mistaken for one another — *A. baumannii* is 0.6323 block-balanced and 0.6329 block-weighted — and
they answer different questions.

## Supplementary Method 3 | Cohort eligibility and retrieval

Assemblies were retrieved from the NCBI Assembly resource and were eligible if the assembly level
was reported as "Complete Genome"; every replicon in the assembly carried an
`assigned_molecule_location_type` of Chromosome or Plasmid; and the organism fell inside the frozen
Gram-negative ESKAPE scope. Assemblies whose replicon inventory contained an unplaced scaffold, an
unlocalised sequence, or a molecule of undeclared type were not eligible, because the location
guarantee that the study depends on would not hold for them.

Retrieval recorded the accession with its version suffix, the retrieval date and a SHA-256 digest of
each downloaded file. Version pinning matters here: RefSeq assemblies are re-annotated, and a
determinant call made against one annotation version cannot be reconciled with a replicon inventory
taken from another. Every downstream table therefore carries the accession.version, not the bare
accession.

The final cohort is 6,288 assemblies. No assembly was added or removed after any outcome column was
read.

## Supplementary Method 4 | Denominator derivation

AMRFinderPlus returned 184,538 records across the cohort. The primary denominator applies five
prespecified filters, in the order given, each of which is an arithmetic subset operation:

1. **Element type.** Only records of `Element type = AMR` are retained. Records typed STRESS,
BIOCIDE, METAL or VIRULENCE describe traits that are not antibiotic resistance and are held in a
separate layer, which is used only for the metal co-location analysis in Extended Data Fig. 2c.
2. **Point mutations.** Records of method type `POINT`, `POINTX`, `POINTN` or `POINTP`, and their
disrupted variants, are excluded. A resistance-conferring substitution in a housekeeping gene is not
an acquired gene, and it has no meaningful compartment: it is located wherever its host gene is.
Point mutations are also detected only where an `--organism` flag exists, so retaining them would
make the denominator depend on which organisms the tool happens to support.
3. **Efflux.** Records whose subclass is an efflux phenotype are removed from the primary
denominator and retained as sensitivity set **S2**. Efflux pumps are frequently core chromosomal
machinery, and including them inflates the chromosomal compartment with genes that were never
acquired.
4. **Scope.** Records with `Scope = plus` that are not efflux are removed and retained as
sensitivity set **S1**. The `plus` scope is a broader, less curated catalogue. One consequence is
worth naming explicitly because it is counter-intuitive: *mcr-9* is a `plus` record and therefore
appears only in S1 and never in a primary result.
5. **Deduplication.** An occurrence is keyed by (assembly accession.version, sequence accession,
start, end). The deposited table contains zero duplicate keys.

The result is **74,349 acquired resistance-gene occurrences**. Supplementary Table 1 gives the full
evidence-layer inventory: the layers are disjoint, they sum exactly to 184,538, and only the first
enters the primary denominator.

## Supplementary Method 5 | Context-block construction

Chromosomal occurrences were merged into shared context blocks so that neighbouring resistance genes
are not counted as independent neighbourhoods. Each block spans the union of the ±10 kb windows of
the occurrences it contains; two occurrences within 20 kb of one another therefore share a block.
This produced **21,955 blocks**, the largest spanning 62,349 bp.

Replicon topology was handled explicitly rather than by clipping. Where a ±10 kb window crossed the
origin of a circular replicon, the window was wrapped and the interval evaluated in both segments
(**57 blocks**). Where a window ran past the end of a linear replicon or an incompletely closed
molecule, the block was flagged truncated (**5 blocks**). Both categories are reported and both were
excluded in sensitivity sets S4 and S5 (Supplementary Table 9); neither changes any contrast beyond
the third decimal place.

## Supplementary Method 6 | The distance estimand

For each chromosomal occurrence the distance recorded is to the nearest qualifying marker lying
**within that occurrence's own ±10 kb window**. Direct overlaps score zero. Occurrences with no
qualifying marker in-window are right-censored at 10,000 bp; a censored distance is never imputed,
and no censored observation is converted to a point value.

The primary estimand is weighted cumulative detection *F*(*d*) — the block-balanced share of
occurrences whose nearest qualifying element lies within *d*. Block balancing assigns weight 1/*m*
to each occurrence in a block of *m* occurrences, so the weights sum to 21,955 and a gene-dense
neighbourhood contributes no more than a sparse one. Restricted mean distance is E[min(*D*,
10,000)]. Medians are reported as "not reached" wherever *F*(10 kb) < 0.5; they are never
extrapolated beyond the censoring horizon.

**The 111 out-of-window occurrences.** Blocks are shared, so a block can contain a marker that is
outside a given occurrence's own window. This happens for 111 occurrences, at marker distances of
10,231–16,931 bp. Using the out-of-block-mate marker would manufacture a distance the design
declares unavailable, so these occurrences are censored at 10 kb. The arithmetic reconciles exactly:
16,414 block-positive occurrences − 111 = 16,303 window-positive = class B.

## Supplementary Method 7 | Structural insertion-sequence census

All 21,955 blocks were processed with ISEScan 1.7.3 in a pinned environment: HMMER 3.3.2, BLAST+
2.17.0, FragGeneScan 1.32, Biopython 1.88. One thread per block; no block was rerun and no block
received a different parameter set. The census completed with **zero tool failures** and resolved
**14,426 elements** (ISEScan type `c` 10,769, type `p` 3,657).

HMMER was pinned to 3.3.2 deliberately. Transposase open-reading-frame detection is sensitive to the
profile search version, an earlier module of this study had used 3.3.2, and a default environment
solve returned 3.4. The environment was rebuilt to 3.3.2 **before** the census was run rather than
reconciled afterwards.

An element is **structurally complete** if ISEScan reports it as type `c` with a complete
transposase open reading frame and bilateral resolved terminal inverted repeats. The primary
structural endpoint additionally requires **full containment within the individual occurrence's own
±10 kb window**. This containment requirement was registered as a numbered amendment and hashed
before any block was scored, because it is the one design choice that could plausibly have been made
after seeing results. Its cost is explicit: 2,384 complete elements lie inside a shared block but
cross the occurrence-window boundary, and are excluded from the primary endpoint (Supplementary
Table 4). A window shorter than an element truncates that element by construction, so counting such
elements would make the endpoint a function of window width rather than of biology.

## Supplementary Method 8 | Statistical procedures

**Uncertainty.** Two interval methods are used and they are not interchangeable. Detection
fractions, contrasts, differences in proportions and the discordance residual carry **95%
percentile** confidence intervals from a cluster bootstrap, taken as the 2.5th and 97.5th
percentiles of the bootstrap distribution. The pooled Mantel–Haenszel odds ratios carry an
**analytic 95% Wald interval on ln(OR) with Robins–Breslow–Greenland variance**, not a bootstrap
interval. No bias-corrected, accelerated, basic or studentized interval is used anywhere in this
study. An earlier revision of this section said every interval was a cluster bootstrap; that was
never true of the pooled odds ratios.

The resampling unit is the **BioProject**, not the occurrence and not the block. Blocks nest within
exactly one BioProject, so the unit is unambiguous, and block weights are recomputed inside each
resample. B = 2,000 throughout. Seeds are recorded for reproducibility and carry no statistical
justification: NM-DIST and NM-IS use 20260822 with BioProject as the unit; NM-V4 uses 20260821 with
BioProject **within species** as the unit, because that analysis compares species.

Four BioProject counts appear in this study and are not interchangeable. The chromosomal analysis set
— the projects contributing context blocks — spans **2,248** BioProjects with an effective number of
**136.9** by the inverse Herfindahl–Hirschman index. The full 6,288-genome cohort spans **2,283** BioProjects at an effective **114.3**. The matched-family bootstrap clusters over **1,404** BioProjects, and the leave-one-BioProject-out arm behind the *A. baumannii*/*K. pneumoniae* ratio runs over the **1,195** projects that contribute to it. Each is a different set of projects selected by a different analysis, and no two are interchangeable. Both the nominal and the effective count are reported
wherever an interval is, because the nominal count materially overstates the independence available.

A Holm-corrected *P* of 0.002 is the **bootstrap resolution floor**, not a smaller value
rounded up. The smallest attainable two-sided bootstrap *P* is 2/B = 0.001 — twice the
one-sided floor, because a two-sided value doubles the smaller tail — and Holm correction
across the two primary contrasts doubles it again to 0.002. It is therefore written
*P* ≤ 0.002.

**Pooling.** Family-matched contrasts use Mantel–Haenszel pooling across matched gene families.
Cochran's *Q*, *I*² and the largest single-family weight share are reported alongside every pooled
estimate, so that an estimate dominated by one stratum is visible as such rather than presented as a
consensus. Where the largest weight share exceeds 50%, the analysis is reported as partial. **This 50% reporting convention is not the registered 30% gate.** NM-V4C's frozen design registers gate G2 — the largest single family carries no more than 30 per cent of the Mantel–Haenszel weight — for the matched-family module only, and a gate failure there is a stronger statement than a partial-reporting flag here. Both are applied to the same statistic and neither is retrofitted to the other: the three-architecture arm at 50.5% is reported partial under this convention, and the within-chromosome arm at 63.9% fails the registered gate.

**Multiplicity.** Within each registered contrast family, Holm's step-down procedure controls the
family-wise error rate. Per-gene-family compartment enrichment, which involves 158 simultaneous
tests, instead uses Fisher's exact test with Benjamini–Hochberg control of the false discovery rate;
family-wise control would be inappropriately conservative for a screen of that size.

**Discovery and confirmation.** The discordance analysis fits the relationship between plasmid share
and block-weighted chromosomal association on eight confirmation species only, excluding both
species in which the pattern was noticed, and evaluates the discovery species against that fit. The
split was registered before the fit was run. *A. baumannii*'s plasmid share falls inside the
confirmation-species range, so the evaluation is an interpolation.

## Supplementary Method 9 | Sensitivity set definitions

| set | definition |
|---|---|
| **S1** | `Scope = plus` non-efflux records restored to the denominator |
| **S2** | efflux records restored to the denominator |
| **S3** | distance measured to the nearest resistance-gene-bearing block rather than to a marker |
| **S4** | truncated blocks excluded (5 blocks) |
| **S5** | circular-wrapped blocks excluded (57 blocks) |
| **S6** | endpoint restricted to insertion-sequence and transposase markers only |
| **SEC_INT** | endpoint restricted to integrase and integron markers only |
| **one-per-block** | deterministic single occurrence per block, replacing 1/*m* weighting |

**Other labelled sets in this document, and where they are defined.** They are separate
namespaces and do not continue this one; S4 here is not BN4, and M1 here is not the M1 of the
adjudication rule engine.

| label | set | defined in |
|---|---|---|
| **BN1–BN7** | the seven background-null sensitivities | Supplementary Result 10 |
| **S8a–S8c** | the intrinsic-determinant exclusions | Supplementary Result 13 |
| **M1, M2** | the two matching units of the matched-family contrast | Supplementary Result 21 |
| **M1–M4, Q, NE1–NE4** | adjudication rule-engine states | Supplementary Method 1 |
| **G1–G4** | the registered background-null decision gates | Supplementary Result 10 |

---


# Supplementary Results

## Supplementary Result 1 | What is in the deposit, and what it answers

Generated from the deposit directory itself, so it cannot drift from what is published.
Row counts are read from the files. DOI 10.5281/zenodo.22116987, CC BY 4.0.

### The five data layers — four in the deposit, one in the repository

**`portabilityrisk_occurrence_portability_v1.tsv`** — 17,751,843 bytes, 74,349 rows.

- *scale*: 74,349 acquired resistance-gene occurrences
- *key columns*: assembly_version, replicon_accession, replicon_molecule_type, determinant_name, gene_family, drug_class, gene_start, gene_end, portability_class, mobility_category
- *the question it answers*: Where is each resistance gene, on a molecule NCBI documents, and how portable is the evidence for it? This is the layer everything else joins to.

**`genome_wide_is_elements.tsv`** — 48,361,465 bytes, 190,999 rows.

- *scale*: 190,999 insertion sequences on 6,190 chromosomes
- *key columns*: replicon_accession, family, cluster, isBegin, isEnd, type, orfLen, irLen, start1..end2
- *the question it answers*: Where is every structurally resolved insertion sequence, genome-wide rather than only in the windows around resistance genes? Supplies the background.

**`NMBG_CHROMOSOME_ACCESSION_MANIFEST_FROZEN_V1.tsv`** — 1,573,219 bytes, 6,190 rows.

- *scale*: 6,190 ARG-bearing chromosomes
- *key columns*: assembly_version, replicon_accession, replicon_length, topology, species, bioproject_accession, download_url
- *the question it answers*: Which chromosomes were annotated, and how does a reader refetch the sequence? Sequence itself is not redistributed; this makes that unnecessary.

**`chromosome_status.tsv`** — 353,149 bytes, 6,190 rows.

- *scale*: 6,190 chromosomes, one status each
- *key columns*: replicon_accession, status, elements
- *the question it answers*: Did the census complete on this chromosome? One row per chromosome so nothing is silently absent.

**`NM_V4C_MLST_CALLS_V1.tsv`** — 4,240 rows. **Not in the Zenodo deposit.** The sequence types
were called after version 2 of the deposit was published, so they are released in the code
repository and will enter the next deposit version rather than being backdated into this one.

- *scale*: 4,240 genomes typed
- *key columns*: assembly_version, organism, scheme, st, profile
- *the question it answers*: which lineage does each genome of the primary contrast belong to? It
  makes clonal structure visible rather than assumed away.

### The dictionaries and the provenance record

| file | what it is for |
|---|---|
| `PORTABILITYRISK_DATA_DICTIONARY_V1.tsv` | every column of every table above, defined |
| `PORTABILITYRISK_DENOMINATOR_DICTIONARY_V1.tsv` | which denominator each rate belongs to |
| `PORTABILITYRISK_EVIDENCE_LEVEL_DICTIONARY_V1.tsv` | what each evidence level does and does not establish |
| `PORTABILITYRISK_PROVENANCE_MAP_V1.tsv` | which artefact produced which published value |
| `PORTABILITYRISK_CLAIM_TO_EVIDENCE_CROSSWALK_V1.tsv` | which evidence supports which claim |
| `PORTABILITYRISK_SOFTWARE_AND_DATABASE_VERSIONS_V1.tsv` | the pinned toolchain |

The frozen design, the primary estimates, the seven background-null sensitivities, the
gate evaluation, the independent verification report and the reconciliation against the
window-limited endpoint are deposited as JSON receipts alongside these.

### What you can compute without re-running anything

The occurrence table is the join key for everything. `assembly_version` joins it to the
chromosome manifest and to the MLST calls; `replicon_accession` joins it to the
genome-wide insertion-sequence table. From those three joins alone a user can recompute
the five-class composition, any per-species or per-family rate, the distance from any
resistance gene to any annotated element, and any lineage-adjusted version of those.

What cannot be recomputed from the deposit: anything needing the **block** layer, because
the context-block table is not deposited — the per-block occupancy figures are the one
quantity in this paper a reader cannot independently derive. The permutation null also
requires the saved matrices rather than the tables here, and the raw chromosome sequence
must be refetched from NCBI using the manifest's `download_url`.

**Frozen versus derived.** The occurrence table, the chromosome manifest and the
insertion-sequence census are frozen records: they are the output of runs whose designs
were hashed beforehand and they are not regenerated. The MLST calls are frozen in the
same sense. Everything in the JSON receipts is derived from those and can be recomputed
from them.

---


## Supplementary Result 2 | Three levels of evidence, and the boundary of each

| level | what it establishes | what it does not establish |
|---|---|---|
| homology proximity | a mobile-element marker lies nearby | that an intact element is present |
| structural detection | a complete element — transposase ORF plus bilateral terminal inverted repeats — lies fully inside the window | that it moved, moves, or can move |
| experimental mobility | — | never measured; no mating, transposition or transfer assay was performed |

Absence is detection-limited throughout. A censored observation means that no qualifying element was
detected within 10 kb under the frozen definition — never that no element is present. Class C
plasmids are not "non-mobilizable" and class A occurrences are not "immobile"; both are statements
about what a marker database contains.


## Supplementary Result 3 | The nested evidence tier inside class E

Moved from the main text under the journal's word limit; not one word is altered.

Class E partitions into two evidence tiers, which are disjoint and not nested: **E2** (9,599 occurrences on 1,109 replicons) carries relaxase, mating-pair formation machinery and a detected origin of transfer; **E1** (16,397 occurrences on 2,828 replicons) is the remainder, carrying the first two without a detected origin of transfer (Fig. 3b). They sum to class E — 16,397 + 9,599 = 25,996 occurrences on 2,828 + 1,109 = 3,937 replicons. Three evidence layers — documented
location, predicted plasmid mobility, sequence-annotated chromosomal mobile-element context — are
kept separable and ranked, so that a revision of the mobility marker database moves C/D/E without
touching A/B or the location layer (Supplementary Result 2).

---


## Supplementary Result 4 | Full distance distribution under the homology endpoint

The species ordering is *A. baumannii* ≫ *Klebsiella* > *P. aeruginosa* at every landmark
(Supplementary Table 3). Two features of the distribution matter more than the headline threshold.

First, **the ±10 kb window is not load-bearing**. The *A. baumannii* excess over *Klebsiella* is
largest at 2 kb (+0.4651), not at 10 kb (+0.3932); against *P. aeruginosa* it is essentially flat
across the whole window (+0.42 to +0.46). The most obvious methodological objection to a fixed
window — that the result is an artefact of where the boundary was drawn — is therefore answered by
the data rather than argued away.

Second, **separation is essentially complete by 2 kb**. *A. baumannii* reaches 0.5878 of its
eventual 0.6323 by 2 kb, so 93% of its detections occur within 2 kb; for *Klebsiella* the same ratio
is 51%. Under the any-resistance-gene block distance (S3) the *A. baumannii* median falls further,
to 279 bp.

Restricted mean distance differences: *A. baumannii* − *Klebsiella* **−4,185 bp** (95% CI −4,588 to
−3,696); *A. baumannii* − *P. aeruginosa* **−4,415 bp** (95% CI −4,805 to −3,953).

## Supplementary Result 5 | Marker-type decomposition

The chromosomal signal is attributable to insertion sequences, not to integrons. Restricting the
endpoint to insertion-sequence and transposase markers (S6) reproduces the primary result almost
exactly: the *A. baumannii* − *Klebsiella* contrast at 1 kb is +0.4186 against +0.4186 in the
primary analysis, and restricted mean distance is 4,181 bp against 4,180 bp.

Restricting instead to integrase and integron markers (SEC_INT) shows **no *A. baumannii* excess at
all**: +0.0100 at 1 kb and −0.0136 at 10 kb, with all three species groups falling between 0.099 and
0.116 at 10 kb. The integrase arm is secondary by design — 3,033 of 32,364 features — and this
decomposition is what converts the finding from "resistance genes are near mobile elements" into a
specific statement about which class of element.

## Supplementary Result 6 | The structural endpoint

Under the stricter structural endpoint the absolute levels fall in every group but the ordering is
unchanged (Supplementary Table 5). Both registered contrasts are positive with 95% confidence
intervals excluding zero at all four landmarks, Holm-corrected *P* ≤ 0.002 throughout — 0.002 being the resolution
floor, not a point estimate.

Short-range concentration is not merely preserved under the structural endpoint but becomes more
distinctive: **86.9%** of *A. baumannii* structural detections fall within 2 kb, against 51.7% for
*P. aeruginosa* and 51.3% for *Klebsiella*. Medians are not reached in any group, because *F*(10 kb)
< 0.5 everywhere; restricted mean distance is reported instead.

## Supplementary Result 7 | Agreement between the homology and structural endpoints

**12,032 of 16,303 class-B occurrences (73.80%)** carry a structurally complete, fully contained
element. By species: *A. baumannii* **86.01%**, *Klebsiella* 65.46%, *P. aeruginosa* 64.24%.

The structural set is very nearly, but not exactly, a subset of the homology set. Of the 12,034
occurrences with a structurally complete contained element, **two** are class A — structurally
positive while homology-negative (Supplementary Table 7). Both are *Enterobacter* occurrences with a
complete element near 4 kb. Their transposases produced no profile hit passing the frozen threshold.
**Why they did not was not tested.** A divergent or unclassified transposase resolved from element
architecture rather than sequence similarity is one hypothesis consistent with the observation; it
is not evaluated in this study and no evidence here supports it over alternatives such as a
truncated or atypical profile match. Neither occurrence falls in a headline species group, so no
reported rate or contrast is affected.

The 26.20% of class-B occurrences without a contained complete element are **not** reclassified. The
structural endpoint is a stricter sequence-structural endpoint, not a correction to the homology
endpoint; a homology marker without a fully resolved element remains valid evidence of
mobile-element context.

The two endpoints are **not independent**. They are computed over the same 21,955 context blocks and
the same 35,140 chromosomal occurrences, under the same block-balancing weights and the same 10 kb
censoring horizon. ISEScan is itself a homology-based caller: its pinned dependencies are HMMER
3.3.2 for profile search, BLAST+ 2.17.0 and FragGeneScan 1.32 for open-reading-frame prediction.
What the structural endpoint adds is a requirement for element architecture — a complete transposase
reading frame with bilateral resolved terminal inverted repeats, fully inside the occurrence's own
window — and nothing else. It is sequence-structural corroboration using complete insertion
sequences, not independent biological validation, and it is reported as such throughout.

## Supplementary Result 8 | Insertion-sequence family retention

IS*6* retains **90.4%** of its elements through the full structural gate (4,385 of 4,848), the
highest of any abundant family. IS*5* retains 85.2%, IS*4* 70.4%, IS*91* 53.7%, IS*3* 41.3%, and
IS*1380* **0%** — 371 elements, none structurally complete under the frozen definition
(Supplementary Table 8). IS*6* is therefore not merely the most abundant family in this cohort but
the most consistently completely resolved, and a retention rate of zero for IS*1380* is a statement
about detectability under this definition, not about that family's biology.

## Supplementary Result 9 | Lineage-balanced weighting of the background null

Registered in `NM_C1_LINEAGE_WEIGHTED_AMENDMENT_008.json` before the arm was computed,
reusing the thresholds already registered for the composite-element and non-ST2 arms
rather than choosing new ones.

The paper argues that BioProject balancing is not clonal-lineage adjustment, and then types
every genome of the primary contrast. But the four registered weightings for the enrichment
— occurrence, one-per-block, genome-balanced, BioProject-balanced — contain no
lineage-balanced arm, and the non-ST2 result is a stratification rather than a weighting.
This applies the paper's own argument to the paper's own primary result.

Each occurrence is weighted by 1/*k*, where *k* is the number of occurrences sharing its
(host, sequence type), so a heavily sequenced clone contributes the same total weight as a
singleton. Untypeable genomes are singleton lineages, not pooled: pooling them would build
one artificial mega-lineage from the most divergent genomes.

| group | *n* | lineages | observed | expected | enrichment |
|---|---:|---:|---:|---:|---:|
| *A. baumannii* | 8,005 | 145 | 0.2719 | 0.0352 | **7.72** |
| *Klebsiella group* | 15,568 | 1128 | 0.0283 | 0.0156 | **1.81** |
| *P. aeruginosa* | — | — | — | — | **NOT EVALUABLE** |
| *Enterobacter group* | — | — | — | — | **NOT EVALUABLE** |

**The registered rule fired.** *A. baumannii*'s enrichment falls to 7.72, below the 8.46
floor, so by the criterion registered in advance the **absolute** enrichment is recorded as
substantially lineage-driven. The manuscript says so.

**The contrast is not.** Balancing lowers *Klebsiella* as far as it lowers *A. baumannii*
— 6.29 to 1.81 — so the ratio between them rises to 4.27, against a registered floor of 1.5.

| adjustment | *A. baumannii* / *Klebsiella* |
|---|---:|
| published | 2.69 |
| composite-flanked removed | 4.05 |
| ST2 removed | 2.03 |
| lineage-balanced | 4.27 |

The contrast is stable, not monotonically strengthened: removing ST2 lowers it to 2.03 and
lineage balancing raises it to 4.27. What matters is that across four adjustments acting on
four different confounders it never changes direction and never approaches 1.

**These lineage counts are per-occurrence, not per-genome.** 145 and 1,128 count the lineages of genomes carrying at least one *chromosomal* occurrence — 778 of the 780 *A. baumannii* genomes and 3,426 of the 3,460 *Klebsiella* genomes. The per-genome counts elsewhere in this document are larger, 120 + 26 = 146 and 715 + 447 = 1,162, because they include the genomes whose resistance genes are all plasmid-borne.

**Two limits.** Sequence types exist only for the 4,240 genomes of the primary contrast, so
*P. aeruginosa* and *Enterobacter* are NOT EVALUABLE here rather than omitted; and the
*Klebsiella* group carries 1,128 lineages against *A. baumannii*'s 145, so balancing
necessarily bites harder on the less diverse host. Neither is a reason to prefer the
unweighted figure; both are reasons to report the pair.

---


## Supplementary Result 10 | The genome-wide background null, in full

Design, permutation count, seed and decision thresholds were registered in
`NM_BACKGROUND_ENRICHMENT_AMENDMENT_001.json` before any species outcome was computed.
B = 2000, seed 20260824. Empirical *P* is (#{null ≥ observed} + 1)/(B + 1); the floor is 1/2001 and is
reported as ≤ 0.0005. *P* = 0 is never reported.

### Primary estimates, occurrence weighting

| group | *n* | landmark | observed | expected | null 95% | enrichment | *P* |
|---|---:|---:|---:|---:|---|---:|---:|
| A. baumannii | 8,005 | 1000 bp | 0.5786 | 0.0342 | 0.0304–0.0381 | **16.91** | ≤0.0005 |
|  |  | 2000 bp | 0.6786 | 0.0484 | 0.0436–0.0531 | **14.02** | ≤0.0005 |
|  |  | 5000 bp | 0.7493 | 0.0874 | 0.0816–0.0936 | **8.58** | ≤0.0005 |
|  |  | 10000 bp | 0.7725 | 0.1322 | 0.1252–0.1397 | **5.85** | ≤0.0005 |
| Klebsiella group | 15,568 | 1000 bp | 0.1206 | 0.0192 | 0.0171–0.0214 | **6.29** | ≤0.0005 |
|  |  | 2000 bp | 0.1497 | 0.0273 | 0.0247–0.0297 | **5.49** | ≤0.0005 |
|  |  | 5000 bp | 0.1987 | 0.0503 | 0.0469–0.0538 | **3.95** | ≤0.0005 |
|  |  | 10000 bp | 0.2466 | 0.0776 | 0.0735–0.0818 | **3.18** | ≤0.0005 |
| P. aeruginosa | 7,150 | 1000 bp | 0.0958 | 0.0106 | 0.0083–0.0130 | **9.07** | ≤0.0005 |
|  |  | 2000 bp | 0.1538 | 0.0150 | 0.0123–0.0179 | **10.24** | ≤0.0005 |
|  |  | 5000 bp | 0.2411 | 0.0274 | 0.0235–0.0312 | **8.81** | ≤0.0005 |
|  |  | 10000 bp | 0.2811 | 0.0409 | 0.0365–0.0457 | **6.87** | ≤0.0005 |
| Enterobacter group | 3,099 | 1000 bp | 0.0720 | 0.0144 | 0.0103–0.0187 | **4.99** | ≤0.0005 |
|  |  | 2000 bp | 0.1000 | 0.0199 | 0.0152–0.0248 | **5.02** | ≤0.0005 |
|  |  | 5000 bp | 0.1500 | 0.0350 | 0.0290–0.0416 | **4.29** | ≤0.0005 |
|  |  | 10000 bp | 0.1794 | 0.0523 | 0.0449–0.0600 | **3.43** | ≤0.0005 |

### Restricted mean distance

| group | observed (bp) | expected (bp) | difference (bp) |
|---|---:|---:|---:|
| A. baumannii | 2,930 | 9,151 | -6,221 |
| Klebsiella group | 8,111 | 9,508 | -1,397 |
| P. aeruginosa | 7,846 | 9,735 | -1,889 |
| Enterobacter group | 8,629 | 9,660 | -1,031 |

### All four weighting schemes, enrichment at 1 kb

| scheme | A. baumannii | Klebsiella group | P. aeruginosa | Enterobacter group |
|---|---:|---:|---:|---:|
| occurrence | 16.91 | 6.29 | 9.07 | 4.99 |
| one per block | 15.04 | 4.10 | 3.87 | 2.78 |
| genome balanced | 15.05 | 4.20 | 5.63 | 2.53 |
| bioproject balanced | 13.96 | 5.43 | 6.43 | 2.87 |

The ordering is preserved under every scheme, with *A. baumannii* first.

### The seven registered sensitivities

These are the **background-null** sensitivities, labelled BN1–BN7. They are a different set
from the distance-endpoint sensitivities S1–S6 and SEC_INT defined in the Supplementary Methods;
earlier revisions used S1–S7 for both, and four of the six labels collided with different
meanings.

**BN1 — chromosome-length stratification.** Enrichment at 1 kb within length tertiles.

| group | short | mid | long |
|---|---:|---:|---:|
| A. baumannii | 16.56 | 17.24 | 16.91 |
| Klebsiella group | 2.54 | 5.53 | 9.60 |
| P. aeruginosa | 4.21 | 10.13 | 9.98 |
| Enterobacter group | 6.15 | 3.79 | 5.29 |

The long-tertile value and the unstratified enrichment both round to 16.91 and are not the same
number: the receipt gives 16.9133 and 16.9073. The coincidence is arithmetic, not a transcription.

**BN2 — genome-wide structural IS density.** Reported beside every enrichment.
Enrichment is already density-normalised by construction: the null relocates within the
same chromosome, so that chromosome's density is preserved exactly.

| group | complete elements | chromosomal bp | per Mb | enrichment (1 kb) |
|---|---:|---:|---:|---:|
| A. baumannii | 26,120 | 3,075,391,409 | 8.493 | 16.91 |
| Klebsiella group | 81,984 | 18,436,925,108 | 4.447 | 6.29 |
| P. aeruginosa | 12,311 | 5,933,048,273 | 2.075 | 9.07 |
| Enterobacter group | 11,558 | 3,334,144,326 | 3.467 | 4.99 |

The four groups hold 131,973 of the 145,779 complete structural elements. The remaining
13,806 lie on ARG-bearing chromosomes of species outside the four groups, which the genome-wide
census covers and this table does not.

**BN3 — IS-family decomposition.** Families holding ≥5% of complete structural elements,
plus IS6 because the amendment names it. Selection fixed before results were read.

| family | elements | *A. baumannii* enrichment (1 kb) | excludes null |
|---|---:|---:|---|
| IS3 | 34,728 | 2.78 | yes |
| IS5 | 33,763 | 7.75 | yes |
| ISNCY | 21,320 | 0.00 | no |
| IS4 | 14,734 | 12.66 | yes |
| IS1 | 8,070 | 0.00 | no |
| IS6 | 4,627 | 35.09 | yes |

Four of 6 families show enrichment with an interval excluding the null. **Two, ISNCY and
IS1, show none within 1 kb.** The effect is carried by a subset of families and is reported
as such.

**BN4 — leave-one-BioProject-out** and **BN5 — leave-one-genome-out.** Registered ceiling 20%.

| group | BN4 refits | BN4 max excursion | BN5 refits | BN5 max excursion | all refits > 1 |
|---|---:|---:|---:|---:|---|
| A. baumannii | 268 | -1.95% | 200 | +0.82% | yes |
| Klebsiella group | 1153 | +9.63% | 200 | -0.77% | yes |
| P. aeruginosa | 425 | -3.93% | 200 | -2.07% | yes |
| Enterobacter group | 270 | +9.28% | 200 | -4.22% | yes |

**BN6 — exclusion of wrapped and truncated contexts.**

| group | excluded | enrichment, all | enrichment, retained |
|---|---:|---:|---:|
| A. baumannii | 39 | 16.9073 | 16.8722 |
| Klebsiella group | 28 | 6.2860 | 6.2541 |
| P. aeruginosa | 9 | 9.0712 | 9.0295 |
| Enterobacter group | 2 | 4.9893 | 4.9925 |

**BN7 — homology endpoint, secondary corroboration.**

**BN7 is occurrence-weighted and genome-wide**, and its numbers are therefore not comparable
row-for-row with Supplementary Tables 3 and 5, which are block-balanced and window-limited. Two
things differ at once — the weighting and whether the element must lie inside the occurrence's
own window — so *A. baumannii*'s *F*(1 kb) appears four times in this document at four values:

| estimate | endpoint | weighting | window | *F*(1 kb) |
|---|---|---|---|---:|
| Supplementary Table 3 | homology | block-balanced | window-limited | 0.5224 |
| Supplementary Table 5 | structural | block-balanced | window-limited | 0.3531 |
| primary background null | structural | occurrence | genome-wide | 0.5786 |
| BN7 | homology | occurrence | genome-wide | 0.6335 |

The same applies to the medians, 647 bp and 326 bp.

| group | *n* | F(1 kb) | F(2 kb) | median (bp) |
|---|---:|---:|---:|---:|
| A. baumannii | 8,005 | 0.6335 | 0.7761 | 326.0 |
| Klebsiella group | 15,568 | 0.1664 | 0.1995 | 1050.0 |
| P. aeruginosa | 7,150 | 0.2485 | 0.3280 | 810.0 |
| Enterobacter group | 3,099 | 0.1284 | 0.1562 | 869.0 |

Host ordering: A. baumannii > P. aeruginosa > Klebsiella group > Enterobacter group.

> **Limitation.** Observed only. No genome-wide homology-marker null exists: the census annotated structurally complete insertion sequences genome-wide, not homology markers, so the homology endpoint cannot be background-normalised here. It corroborates host ORDERING, not enrichment.

### The registered decision gates

| gate | condition | result |
|---|---|---|
| G1 | A. baumannii enrichment exceeds 1 with an interval excluding 1 at both 1 kb and 2 kb | **PASS** |
| G2 | the background-normalised A. baumannii contrast remains positive against both comparison groups | **PASS** |
| G3 | direction survives BioProject balancing | **PASS** |
| G4 | no single BioProject or genome dominates | **PASS** |

Verdict: **BACKGROUND_NORMALISED_HOST_CONDITIONING_SUCCESS**, one of the four registered in the amendment before any species
outcome was computed.

### Independent verification

A separate code path that does not import the primary analysis and re-implements every
mechanic brute force: **20 of 20 checks pass** (`C1_INDEPENDENTLY_VERIFIED`). The brute-force distance check
reproduced the primary on 400 of 400 sampled occurrences. The census itself was verified by
36 independent ISEScan re-runs, all of which reproduced the recorded element count exactly.


## Supplementary Result 11 | The genome-wide census against the window-limited endpoint

NMIS remains the frozen window-contained endpoint of the manuscript. C1 is the genome-wide background analysis. C1 does not replace NMIS and the two must not be presented as one measurement.

| | pooled |
|---|---:|
| occurrences | 35,140 |
| window-limited, complete and contained | 12,034 |
| genome-wide, detected | 13,238 |
| **net change** | **+1204** |
| gains — boundary | 1,154 |
| gains — context | 384 |
| losses — context | 334 |
| unmatched | 0 |

**BOUNDARY** — NMIS could not call the element complete-and-contained because the extracted block truncated it or it crossed the window boundary; whole-chromosome annotation removes the truncation. NMIS declared this limitation in its own protocol.

**CONTEXT** — the ISEScan prediction itself differs between block and whole-chromosome scope, because gene prediction and element clustering see different neighbours. Runs in both directions.

### Per species

The fourteen species with the most chromosomal occurrences, 33,100 of 35,140 (94.2%). The
cutoff falls between 149 and 137 occurrences; the 2,040 not shown are spread over species with
fewer than 149 each, three of the largest of which are unresolved `Pseudomonas sp.`,
`Enterobacter sp.` and `Acinetobacter sp.` labels.


| species | *n* | window-limited | genome-wide | net | boundary gains | context losses |
|---|---:|---:|---:|---:|---:|---:|
| *Klebsiella pneumoniae* | 13,082 | 3,217 | 3,527 | +310 | 273 | 33 |
| *Acinetobacter baumannii* | 8,005 | 5,571 | 6,184 | +613 | 496 | 96 |
| *Pseudomonas aeruginosa* | 7,150 | 1,838 | 2,010 | +172 | 278 | 176 |
| *Enterobacter hormaechei* | 1,749 | 320 | 326 | +6 | 13 | 9 |
| *Klebsiella variicola* | 646 | 16 | 30 | +14 | 13 | 0 |
| *Klebsiella quasipneumoniae* | 573 | 38 | 53 | +15 | 13 | 0 |
| *Klebsiella aerogenes* | 351 | 44 | 43 | -1 | 0 | 1 |
| *Enterobacter cloacae* | 312 | 81 | 89 | +8 | 8 | 1 |
| *Klebsiella michiganensis* | 298 | 95 | 101 | +6 | 3 | 2 |
| *Enterobacter asburiae* | 245 | 43 | 49 | +6 | 6 | 3 |
| *Klebsiella oxytoca* | 191 | 25 | 28 | +3 | 0 | 0 |
| *Acinetobacter pittii* | 179 | 31 | 37 | +6 | 5 | 0 |
| *Enterobacter roggenkampii* | 170 | 27 | 30 | +3 | 3 | 0 |
| *Enterobacter ludwigii* | 149 | 7 | 7 | +0 | 0 | 0 |

Three quarters of the gains are boundary artefacts the window-limited protocol declared in
its own limitations. The remainder differ because insertion-sequence prediction sees
different neighbours at whole-chromosome scope, and that difference runs in both
directions. The window-limited endpoint is not superseded: it remains the frozen structural
measure of this work, and the genome-wide census is the background against which it is
normalised.


## Supplementary Result 12 | Two explanations the enrichment is not

Registered in `NM_C1_COMPOSITE_ELEMENT_AMENDMENT_006.json` and
`NM_C1_NON_ST2_ENRICHMENT_AMENDMENT_007.json`, each before its own arm was computed. The
thresholds in the second were **reused from the first rather than chosen again**, because
picking a fresh threshold while knowing what the sibling arm returned is what a freeze
exists to prevent.

### Clonal replication

*A. baumannii* is dominated by one lineage. Sequence type 2 holds **48.2%** of its genomes and **63.0%** of its chromosomal
occurrences; the two are not the same number and an earlier revision of this section printed
the second in both places. If the enrichment were a property of that clone,
removing it would remove the effect.

| stratum | *n* | observed | expected | enrichment |
|---|---:|---:|---:|---:|
| A. baumannii, all | 8,005 | 0.5786 | 0.0342 | **16.91** |
| A. baumannii, ST2 only | 5,047 | 0.6788 | 0.0355 | **19.10** |
| A. baumannii, non-ST2 | 2,958 | 0.4077 | 0.0320 | **12.74** |
| Klebsiella group | 15,568 | 0.1206 | 0.0192 | **6.29** |
| P. aeruginosa | 7,150 | 0.0958 | 0.0106 | **9.07** |

It does not. Non-ST2 *A. baumannii* retains **12.74-fold** enrichment against a registered
floor of 8.46, and **2.02 times** the *Klebsiella* group against a floor of 1.5. The
ratio to *P. aeruginosa* falls to 1.41 from a published 1.86; no floor was registered for
that comparison and none is claimed for it here.

### Composite-transposon structure

Where a resistance gene is the cargo of a composite element, its distance to an insertion
sequence is near zero by construction: the flanking copies are what make it a transposon. A
uniform relocation null cannot address that, so the occurrences were stratified. An
occurrence counts as composite-flanked when a complete-structural element of family *F* lies
entirely upstream within *D* bp and another of the **same family** lies entirely downstream
within *D* bp, on its own replicon.

| distance | flanked | share of occurrences in the four analysis groups (*n* = 33,822) |
|---|---:|---:|
| *D* = 10 kb | 4,324 | 12.8% |
| *D* = 5 kb | 3,145 | 9.3% |

| group | flanked | non-flanked |
|---|---:|---:|
| *A. baumannii* | 23.64 (*n* = 2,544) | **13.58** (*n* = 5,461) |
| *Klebsiella group* | 29.02 (*n* = 1,316) | **3.35** (*n* = 14,252) |
| *P. aeruginosa* | 18.80 (*n* = 280) | **8.24** (*n* = 6,870) |
| *Enterobacter group* | 31.06 (*n* = 184) | **2.13** (*n* = 2,915) |

The flanked stratum is enriched 23–31-fold in every host, which is what a definitional
signal looks like and confirms the detector finds what it should. But it is a minority, and
among the occurrences that are **not** flanked *A. baumannii* retains **13.58-fold**
enrichment and **4.05 times** the *Klebsiella* group — higher than the unstratified ratio of
2.69, not lower. The objection is real for the flanked minority and does not explain the
effect.

A flanking pair of the same family is necessary but not sufficient for a composite
transposon. Some flanked occurrences will be coincidental, and true composites with one
degenerate copy will be missed. No transposon database was consulted and no element is named.

---

## Supplementary Result 13 | Intrinsic determinants do not carry the host contrast

Registered in `NM_V4C_INTRINSIC_SENSITIVITY_AMENDMENT_001.json` before any arm was computed.
The concern is that a determinant which is an intrinsic chromosomal gene of its host would
be chromosomal for reasons that have nothing to do with mobility.

The matched-family design already excludes most such genes, and not by intention. The
registered eligibility rule requires a family to appear in at least three species; a
species-core gene cannot satisfy that. *fosA*, *oqxA*, *oqxB*, *bla*PDC, *bla*ACT and
*bla*LEN are all present in the cohort and all absent from the 58 matched families for that
reason.

Two eligible families do carry intrinsic members: *bla*OXA, which in *A. baumannii* lumps the
intrinsic chromosomal OXA-51-like enzyme with acquired OXA-23, OXA-24 and OXA-58; and
*bla*SHV, which in *K. pneumoniae* is partly chromosomal. Allele-to-family assignment was
taken from the AMRFinderPlus reference protein database at version 2026-08-07.1, the version
that produced the calls. AMRFinderPlus does not subdivide *bla*SHV into intrinsic and
acquired nodes, so no allele rule for it was invented; the conservative arm removes the whole
family instead.

| arm | families | odds ratio | 95% CI | concordant |
|---|---:|---:|---|---:|
| baseline, as published | 58 | **50.29** | 45.61–55.45 | 56 |
| **S8a** OXA-51-like alleles removed | 58 | **48.45** | 43.94–53.41 | 56 |
| **S8b** *bla*OXA and *bla*SHV removed whole | 56 | **51.53** | 46.28–57.37 | 54 |
| **S8c** S8b plus a species-core screen | 55 | **51.51** | 46.26–57.35 | 53 |

Removing the intrinsic determinants does not shrink the contrast; the two conservative arms
sit above the unadjusted estimate. The objection is answered in the direction opposite to the
one it predicted.

This arm covers the matched-family contrast only. The distance analysis, the 46.39%
chromosomal association and the genome-wide background null are not family-matched and are
not addressed by it.

---


## Supplementary Result 14 | Accessory context: a local null and an island stratification

Registered in `NM_C1_ACCESSORY_CONTEXT_AMENDMENT_009.json` before either arm was computed,
reusing the contrast floor of 1.5 already registered for the composite-element, non-ST2
and lineage-balanced arms rather than choosing a new one.

The registered null relocates an occurrence anywhere on its own chromosome. That preserves the
chromosome's element density exactly, but not the occurrence's *position within* the chromosome:
resistance genes sit in accessory regions, insertion sequences accumulate in the same regions,
and a chromosome-wide null will mostly place an occurrence in core sequence where elements are
sparse. The enrichment could then be a statement about where accessory DNA is rather than about
where resistance genes are within it.

### ARM LOC — a local relocation null

Identical to the registered null in every respect except the relocation range: each
occurrence moves uniformly within ±*R* of its own start instead of anywhere on the
chromosome. Interval length, chromosome, genome, species and BioProject identity are
preserved; circular chromosomes wrap; linear chromosomes exclude positions whose ±10 kb
window would run past an end. *R* = 50 kb and *R* = 100 kb were registered together, so
neither could be chosen after seeing the other. 2,000 permutations, seed 20260827, seeded
per chromosome. ISEScan was not rerun.

One deliberate difference from the registered null: relocated intervals are not forced to be
mutually non-overlapping. Inside a ±50 kb window around a cluster of resistance genes that
constraint is often unsatisfiable, and it does not enter this estimator — each occurrence's
distance depends only on its own interval and the frozen element coordinates, never on where
the other occurrences landed.

| group | *n* | observed | expected | enrichment | excludes null |
|---|---:|---:|---:|---:|---|
| **±100 kb** | | | | | |
| *A. baumannii* | 8,005 | 0.5786 | 0.0804 | **7.19** | yes |
| *Klebsiella group* | 15,568 | 0.1206 | 0.0366 | **3.30** | yes |
| *P. aeruginosa* | 7,150 | 0.0958 | 0.0274 | **3.50** | yes |
| *Enterobacter group* | 3,099 | 0.0720 | 0.0282 | **2.55** | yes |
| **±50 kb** | | | | | |
| *A. baumannii* | 8,005 | 0.5786 | 0.1321 | **4.38** | yes |
| *Klebsiella group* | 15,568 | 0.1206 | 0.0465 | **2.59** | yes |
| *P. aeruginosa* | 7,150 | 0.0958 | 0.0442 | **2.17** | yes |
| *Enterobacter group* | 3,099 | 0.0720 | 0.0403 | **1.78** | yes |

**The accessory explanation is real and large, and it is not the whole effect.** The
expected detection fraction for *A. baumannii* rises from 0.034 under the chromosome-wide
null to 0.0804 at ±100 kb and 0.1321 at ±50 kb: the neighbourhoods around resistance genes are
far more element-dense than their chromosomes are. The enrichment falls with it, 16.91 to
7.19 to 4.38. It does not vanish. Every group still excludes its null interval, and the
*A. baumannii* / *Klebsiella* contrast holds at 2.18 and 1.69, above the registered floor of 1.5
but below every previous adjustment. **1.69 is the new minimum of the reported contrast range.**

### ARM ISL — a resistance-island stratification

Occurrences on a replicon were single-linkage clustered at a 10,000 bp gap; a cluster of ≥3
occurrences was treated as a resistance island. 9,498 of 35,140 occurrences (27.0%) fall in one.
The enrichment was then recomputed within each stratum from the existing null matrices, with
no new permutation.

| group | stratum | *n* | observed | expected | enrichment |
|---|---|---:|---:|---:|---:|
| *A. baumannii* | island | 4,384 | 0.6307 | 0.0344 | **18.32** |
| *A. baumannii* | isolated | 3,621 | 0.5156 | 0.0340 | **15.18** |
| *Klebsiella group* | island | 1,612 | 0.4721 | 0.0228 | **20.69** |
| *Klebsiella group* | isolated | 13,956 | 0.0800 | 0.0188 | **4.26** |
| *P. aeruginosa* | island | 2,457 | 0.2475 | 0.0138 | **17.92** |
| *P. aeruginosa* | isolated | 4,693 | 0.0164 | 0.0089 | **1.85** |
| *Enterobacter group* | island | 327 | 0.4618 | 0.0187 | **24.73** |
| *Enterobacter group* | isolated | 2,772 | 0.0260 | 0.0139 | **1.87** |

**The two hosts behave differently, and that is the result.** *Klebsiella*'s enrichment is
largely an island phenomenon — 20.69 inside multi-gene islands against 4.26 outside. *A. baumannii*'s
is not — 18.32 against 15.18. The contrast therefore *rises* outside islands, to 3.56, against
2.69 unstratified. *P. aeruginosa* and the *Enterobacter* group behave like *Klebsiella*
(1.85 and 1.87 outside islands).

### What this is not

ARM ISL is a **structural proxy**. It is not a core-versus-accessory partition, and it does
not identify AbaR or AbGRI islands by their backbone or their *comM* insertion site; it
identifies multi-gene resistance clusters, which is a broader and cruder category. ARM LOC is
an **approximation** to an accessory-restricted null, not that null: no core-genome
alignment exists in this project, so the true accessory partition cannot be drawn from any
tracked artefact. Both are reported as what they are.

Coding versus intergenic position, local gene density and recombination hotspots remain
unaddressed by either arm.

---


## Supplementary Result 15 | The matched-family analysis in full

Matched-family contrasts ask whether the *same* gene family occupies different routes in different
hosts, which removes gene-composition differences between species as an explanation. Two frozen
modules produce them: NM-V4C (two-way, chromosomal-and-mobile against plasmid-borne) and NM-V4D
(three-architecture, over the partition chromosomal_quiescent 18,837 / chromosomal_mobile 16,303 /
plasmid_borne 39,209).

### Estimation

Per-family effects are Woolf (log) odds ratios with a Haldane–Anscombe 0.5 continuity correction
applied where a cell is empty. The **pooled** estimate is Mantel–Haenszel and carries **no**
continuity correction. A stratum with an empty margin — no exposure or no outcome in that family —
is dropped rather than corrected. This policy exists only in the scoring source; it is stated here
for the first time.

### The primary contrast, *A. baumannii* against *Klebsiella*

| arm | families | OR | 95% CI | up / down | *I*² |
|---|---:|---:|---|---|---:|
| headline, occurrence-level | 58 | **50.29** | 45.61–55.45 | 56 / 2 | 83.6% |
| block-balanced | 58 | 24.43 | 21.12–28.25 | 55 / 3 | 72.9% |
| BioProject-balanced | 49 | **21.95** | 18.92–25.47 | 43 / 4 | 75.6% |

**The up/down counts exclude exact ties** and therefore need not sum to the family total. A
family whose Woolf log odds ratio is exactly zero counts in neither column. Two such families
sit in the BioProject-balanced arm here (43 + 4 + 2 = 49) and two in the *A. baumannii* against
*P. aeruginosa* contrast below (24 + 22 + 2 = 48). The convention was not stated in earlier
revisions, which made those two rows look like arithmetic errors.

Largest single-family Mantel–Haenszel weight share and family-weighted estimate, by arm: headline
15.47% (*bla*OXA) and 40.61; block-balanced 9.68% and 24.22; BioProject-balanced 15.46% and 16.72.
Cochran's *Q* is 348.17 on 57 degrees of freedom for the headline arm and 197.06 on 48 for the
BioProject-balanced arm. All three arms exclude 1 and point the same way; the effect is large under
every weighting.

The **block-balanced arm is the least weight-concentrated and the most internally concordant**: its
largest single family carries 9.68% rather than 15.47%, its top three families carry 26.46%
combined, and its Mantel–Haenszel estimate (24.43) and family-weighted estimate (24.22) agree
closely, which they do not in the headline arm (50.29 against 40.61). Fifty-five of its 58 families
point the same way. Reported together, the three arms bracket the effect at roughly 22–50-fold
rather than resting on a single weighting.

**Influence.** Leave-one-family-out moves the pooled log odds ratio by at most **2.68%** (worst
family *sul1*) and leave-one-BioProject-out by at most **2.02%** (worst project PRJNA789460),
against a **registered ceiling of 20%** relative. Zero families and zero BioProjects exceed it. In
the three-architecture module the same statistic is 3.83% (worst *sul1*) for the block-balanced
contrast. A permutation test on the primary contrast gave an empirical *P* of 0.0005 at 2,000
permutations (seed 20260821) — again a resolution floor.

**The two families pointing the other way** in the headline arm are *dfrA16* (OR 0.78, weight
0.0064%) and *tet(39)* (OR 0.067, weight 0.0000%); both carry essentially no pooled weight. The
heaviest family, *bla*OXA, has 938 chromosomal-and-mobile against 224 plasmid-borne occurrences in
*A. baumannii* and 123 against 1,301 in *Klebsiella*, giving OR 44.29 — close to the pooled estimate
rather than driving it away from the others.

### The qualified and partial contrasts

Against *P. aeruginosa* the same two-way analysis gives OR **1.26** (95% CI 1.11–1.43) over 48
families, with families splitting 24 up against 22 down; its BioProject-balanced co-primary gives
1.18 (95% CI 0.98–1.41), which does **not** exclude 1 and splits 17 up against 20 down. This
contrast is reported as **qualified** for that reason: the interval in one arm excludes 1, the
direction is inconsistent across families, and the co-primary arm does not replicate it.

The three-architecture partition PC1 (*A. baumannii* against *P. aeruginosa*, chromosomal_mobile
against chromosomal_quiescent) gives OR **16.54** (95% CI 11.88–23.01) over 45 families, reproducing
the direction. One family carries **50.5%** of the pooled weight and the family-weighted estimate is
1.41 against a Mantel–Haenszel 16.54, so this sub-analysis is reported as **partial** and its
pre-registered weight-ceiling sub-claim is **false**, not merely unsupported.

### What is not available

Per-family confidence intervals are not stored as such: the exported table carries the odds ratio, its log and its variance but no interval column. They are **derivable** from what is stored, and are now derived and drawn — see the next Supplementary Result. An earlier revision said no file in the study held them, which confused *not stored* with *not obtainable*. Per-family
leave-one-family-out values are likewise not persisted — only the maximum and the count over
threshold. Cochran's *Q* is reported by NM-V4C but not by the three-architecture module, which
stores *I*² only. These are limitations of the frozen runs and are stated rather than filled in.


## Supplementary Result 16 | Per-family intervals for the matched-family contrast

An earlier revision of Figure 1 drew no per-family interval and said they had never been computed
or stored. That was wrong. The exported table `nmv4c_family_host_vehicle.tsv` carries `ln_or` and
`var_ln_or`, so the Wald interval is exp(ln(OR) ± 1.96·√var) — arithmetic on an existing column,
not a new analysis, and no part of the frozen record is altered by computing it.

The check that licenses drawing them: the stored variance reproduces the Woolf variance computed
from the table's own 2×2 counts for **58 of 58** families.

**44 of 58 families exclude 1 individually.** The 14 that do not are
*aac(6')-Ib4*, *aph(4)*, *blaCARB*, *blaGES*, *blaIMP*, *blaPER*, *blaSHV*, *catB2*, *dfrA16*, *dfrB1*, *erm(B)*, *rmtC*, *tet(39)*, *tet(M)*.

They are drawn in grey in Figure 1a rather than omitted. Individual families being uncertain is
expected at these counts and is not an argument against the pooled estimate; hiding which ones
they are would have been.

---


## Supplementary Result 17 | Why the within-chromosome contrast is not reported

Registered in `NM_V4C_REFEREE_ANALYSES_AMENDMENT_004.json` and
`NM_V4C_MARGINAL_ADJACENCY_AMENDMENT_005.json`. This section reports an analysis that
**failed**, because the question it asks is the right one and a reader is entitled to know
why it has no answer here.

The two-way endpoint asks whether a determinant is chromosomal-and-mobile *rather than*
plasmid-borne. The question the resource's own framing raises is narrower: given that a
determinant is chromosomal, is it marker-adjacent? That is class B against class A.

| arm | odds ratio | 95% CI | largest family weight | informative families |
|---|---:|---|---:|---|
| W1 as the referee specified it | **0.48** | 0.33–0.70 | blaOXA 63.9% | 11 of 49 |
| W1 minus intrinsic Acinetobacter OXA (NM-V4C-001 exclusion) | **3.68** | 2.13–6.37 | aph(3') 39.4% | 11 of 49 |
| W1 minus the whole blaOXA family | **6.11** | 2.93–12.75 | aph(3') 46.4% | 10 of 48 |
| W1 minus intrinsic OXA and minus ST2 | **2.51** | 1.41–4.48 | aph(3') 43.8% | 9 of 46 |

**Every arm fails gate G2**, registered in `NMV4C_FROZEN_DESIGN.json` long before this
analysis was requested: *the single largest family carries no more than 30 per cent of the
MH weight*. And the estimate rests on far fewer families than the count suggests — a family
in which both hosts have zero class-A occurrences contributes zero to both Mantel–Haenszel
sums and disappears silently.

The sign is not stable either. It depends entirely on how one family is handled: with
*bla*OXA included the direction favours *Klebsiella*; with the intrinsic OXA-51-like alleles
removed it favours *A. baumannii*. In this cohort that comparison is between the intrinsic
chromosomal enzyme of one host and the integron-borne acquired enzymes of the other, which
is not an architectural comparison.

### The marginal form is stable but compositional

The rate is estimable where the pooled odds ratio is not: **80.9%** of *A. baumannii*
chromosomal occurrences are marker-adjacent against **34.2%** in the *Klebsiella* group, a
difference of +46.7 percentage points whose interval excludes zero under both
BioProject-clustered and lineage-clustered resampling, and which survives collapsing to one
genome per sequence type (65.2% against 24.9%).

Direct standardisation then removes it. Applying each host's family-specific rates to the
other's family distribution accounts for **136%** of the crude gap. The reason is visible
in the data: most of the *Klebsiella* class-A pool sits in *fosA*, *bla*SHV, *oqxA* and
*oqxB* — species-core genes that do not occur in *A. baumannii* at all. The two hosts'
chromosomal resistomes are built from different families, and that, rather than a difference
in how shared families sit, is what the marginal gap measures.

**Neither form of the within-chromosome contrast is reported as a result.** The pooled form
fails a registered gate; the marginal form measures composition. Both are here so that the
next person does not have to rediscover it.

---


## Supplementary Result 18 | The discordance analysis in full

### The registered design

Ten species met a floor of ≥40 genomes **and** ≥200 chromosomal resistance-gene occurrences. Two —
*Acinetobacter baumannii* and *Klebsiella pneumoniae* — were designated **discovery** species,
because the discordance pattern was noticed in them, and were excluded from the fit. The remaining
**eight** are the confirmation panel. The design, including this split, was frozen at
2026-08-22T00:48:01Z with `frozen_before_any_outcome_was_computed: true`, and the run began 93
seconds later; the frozen design file deliberately excludes the outcome columns *P* and *M*.

### The model

Both proportions are transformed to the logit scale and an unweighted ordinary least-squares line is
fitted across the eight confirmation species:

> **logit(*M*) = −1.421185 + 0.080631 × logit(*P*)**,  in-sample *R*² = **0.027524**

where **P** is the occurrence-weighted plasmid share of that species' resistance-gene occurrences
(*n* plasmid occurrences / *n* occurrences) and **M** is the **block-weighted** chromosomal
mobile-element association (blocks carrying ≥1 marker / all context blocks on that species'
chromosomes). *M* is block-weighted rather than occurrence-weighted by prior registration, because
marker-positive blocks carry 2.46 resistance genes on average against 1.23 for marker-negative
blocks, and occurrence weighting would manufacture discordance.

These two are **means over blocks**, not ratios of totals: the per-block occurrence count averaged
over blocks of each kind. The ratio of totals is a different and slightly larger quantity —
16,414/6,617 = 2.48 against (35,140 − 16,414)/(21,955 − 6,617) = 1.22 — and the two should not be
read as the same number. The block-level table this is computed from is held in the working tree
and is not part of the deposit, so a reader cannot recompute it from the published data; that is a
gap in the deposit rather than in the analysis, and it is stated here rather than left to be found.

### The eight confirmation species, complete

| species | occurrences | blocks | plasmid share *P* | chromosomal association *M* | logit *P* | logit *M* | LOSO residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| *Enterobacter asburiae* | 466 | 180 | 0.4742 | 0.2722 | −0.1031 | −0.9834 | +0.5099 |
| *Enterobacter cloacae* | 783 | 221 | 0.6015 | 0.2805 | +0.4119 | −0.9418 | +0.5426 |
| *Enterobacter hormaechei* | 5,260 | 1,173 | 0.6675 | 0.1586 | +0.6969 | −1.6689 | −0.4071 |
| *Klebsiella aerogenes* | 672 | 245 | 0.4777 | 0.1347 | −0.0893 | −1.8601 | −0.4934 |
| *Klebsiella michiganensis* | 718 | 177 | 0.5850 | 0.2825 | +0.3432 | −0.9322 | +0.5521 |
| *Klebsiella quasipneumoniae* | 1,391 | 430 | 0.5881 | 0.1419 | +0.3560 | −1.7999 | −0.4889 |
| *Klebsiella variicola* | 1,075 | 488 | 0.3991 | 0.1537 | −0.4093 | −1.7060 | −0.2946 |
| *Pseudomonas aeruginosa* | 8,151 | 4,102 | 0.1228 | 0.1767 | −1.9661 | −1.5386 | +0.2613 |

The two excluded discovery species, for reference: *A. baumannii* 9,304 occurrences, 4,138 blocks,
*P* = 0.1396, *M* = 0.6329; *K. pneumoniae* 39,918 occurrences, 9,075 blocks, *P* = 0.6723, *M* =
0.2575.

**Per-species uncertainty is not available.** The bootstrap resamples every species but retains only
five derived aggregates, so no per-species confidence interval for *P* or *M* was computed or stored
and none is reported here. The only per-species uncertainty quantity computed at all is the binomial
standard error of logit(*M*), and only its mean across the eight (0.1320) survives into the receipt.
This is a limitation of the frozen run, not an omission from this table.

### The three registered tests

**Occurrence-weighted associations, for the four species the manuscript names.** These are a
different denominator from the block-weighted *M* above and the two must not be read down one
column: occurrence-weighted counts every chromosomal occurrence once, block-weighted gives each
context block weight 1. Each value is class B over classes A and B, and every one recomputes
from the deposited occurrence table.

| species | class A | class B | association, occurrence-weighted | plasmid share |
|---|---:|---:|---:|---:|
| *A. baumannii* | 1,528 | 6,477 | **80.91%** | 13.96% |
| *K. pneumoniae* | 8,276 | 4,806 | **36.74%** | 67.23% |
| *P. aeruginosa* | 4,289 | 2,861 | **40.01%** | 12.28% |
| *E. hormaechei* | 1,329 | 420 | **24.01%** | 66.75% |

*A. baumannii*'s block-weighted association is 0.6329 against this 80.91%; the gap is the
denominator, not a disagreement.

**T2 — the discovery-species residual.** *A. baumannii*'s plasmid share (0.1396) lies inside the
confirmation range [0.1228, 0.6675], so evaluating it against the fit is interpolation. Predicted
logit(*M*) = −1.5678 (*M* = 0.1725); observed logit(*M*) = +0.5447 (*M* = 0.6329); residual
**+2.1126** logits, bootstrap median 2.1194, 95% CI **1.8065 to 2.4376**, excluding zero. The
corresponding *K. pneumoniae* residual is +0.3044, an order of magnitude smaller and reported
without an interval because none was stored.

**The slope.** Bootstrap median 0.0849, 95% CI **−0.0588 to 0.2287** — including zero. Together with
*R*² = 0.0275 this is the substantive finding: across the confirmation panel, plasmid share carries
almost no information about chromosomal mobile-element association.

**T3 — the low-plasmid contrast.** *A. baumannii* (*P* = 0.1396) against *P. aeruginosa* (*P* =
0.1228): gap in *M* of **0.4548**, 95% CI 0.4021 to 0.5039, excluding zero. That figure is the
**bootstrap median** of the gap (`M_gap_median` in the result receipt), not the difference of
the two point estimates, which is 0.6329 − 0.1767 = 0.4562. The two differ by 0.0014 and the
median is the registered quantity. Two species with almost
identical plasmid share differ by 45 percentage points in chromosomal association.

**T1 — leave-one-species-out, reported here for the first time.** This test was computed in the
frozen NM-V4 run and is **load-bearing for that run's SUCCESS verdict**, but was not reported in any
previous version of the manuscript or Supplement. It is exposed here in full. The procedure
regresses logit(*M*) on logit(*P*) over the confirmation set minus one species, predicts the
held-out species, and asks whether the residuals are systematically larger than sampling noise; if
plasmid share alone predicted chromosomal association, they would not be. The eight residuals are in
the table above — four positive, four negative, none small. Their standard deviation is **0.4867**
against a sampling-noise floor of **0.1320** (the mean within-species binomial standard error of
logit *M*), a ratio of **3.69**. The 95% bootstrap CI for that standard deviation is **0.3840 to
0.8665**, whose lower bound exceeds the noise floor; that comparison is the registered T1 gate
condition, and it passes.

### What the panel does and does not establish

Plasmid fraction and chromosomal mobile-element association are **non-redundant across the
registered species panel**: the slope cannot be distinguished from zero, *R*² is 2.75%, the
leave-one-species-out residuals are 3.69 times sampling noise, and two species with nearly identical
plasmid share differ by 45 percentage points in chromosomal association. The panel is eight species
and the fit is unweighted ordinary least squares on eight points; it establishes non-redundancy
across these species and does not extrapolate to Gram-negative bacteria at large.


## Supplementary Result 19 | Leverage in the eight-species fit

The discordance argument rests on an ordinary least-squares fit of logit(*M*) on logit(*P*) over
eight confirmation species. *P. aeruginosa* sits at logit(*P*) = −1.97 while the other seven lie
between −0.41 and +0.70, so the fit was checked for a single influential point. It has one.

| species | hat | in-sample residual | Cook's *D* | slope without it |
|---|---:|---:|---:|---:|
| *Enterobacter asburiae* | 0.1250 | +0.4461 | 0.0870 | +0.0815 |
| *Enterobacter cloacae* | 0.1777 | +0.4462 | 0.1400 | +0.0243 |
| *Enterobacter hormaechei* | 0.2535 | -0.3039 | 0.1125 | +0.1467 |
| *Klebsiella aerogenes* | 0.1250 | -0.4317 | 0.0815 | +0.0812 |
| *Klebsiella michiganensis* | 0.1643 | +0.4613 | 0.1341 | +0.0311 |
| *Klebsiella quasipneumoniae* | 0.1667 | -0.4074 | 0.1066 | +0.1258 |
| *Klebsiella variicola* | 0.1452 | -0.2518 | 0.0337 | +0.0617 |
| *Pseudomonas aeruginosa* | **0.8425** | +0.0412 | 0.1540 | +0.1809 |

**One point carries 84% of the leverage.** With eight observations the average hat value is
0.25; *P. aeruginosa* has 0.8425. The slope is +0.0806 with all eight and +0.1809 without it, a
factor of 2.2, and across the eight leave-one-out fits the slope ranges from +0.024 to +0.181.

**The conclusion is unchanged and the fit is not stable, and both are true.** *R*² is 0.0275 with
all eight and 0.0252 without *P. aeruginosa*; the slope interval covers zero either way. What the
diagnostic removes is any reading of the slope itself as an estimate. The published T1 test
reported leave-one-species-out residuals only, which is why this was not visible before.

The fit reproduces from the published eight-species table to four decimal places, so this
diagnostic describes the fit the paper reports and not a re-fit of it.

---

## Supplementary Result 20 | Cargo convergence is confounded with replicon length

Conjugation-consistent replicons are larger, and larger replicons carry more of everything.
The unstratified difference in the share carrying three or more drug classes is **+19.86
percentage points** (95% CI 14.16–25.41). Stratifying by replicon length and pooling across
strata leaves **+8.10 percentage points**.

| stratum | length range (bp) | *n* marker-negative | difference (pp) |
|---:|---|---:|---:|
| 0 | 0–53097 | 351 | -13.84 |
| 1 | 53098–87948 | 194 | -14.3 |
| 2 | 87949–128536 | 551 | 34.06 |
| 3 | 128537–205642 | 186 | -4.69 |
| 4 | 205643–max | 191 | 8.71 |

**The strata disagree in sign.** Three of five are negative and one, stratum 2, carries
+34.06 points on its own. A pooled difference over strata that point in opposite
directions does not describe a common effect, and the unstratified +19.86 should not be
read as a size-independent one. The manuscript reports both.

This is a length stratification, not a length adjustment: the strata are the boundaries
the phase-3 analysis had already fixed, no interval was computed on the pooled value, and
no model of the length effect is fitted. It establishes that the crude figure is
confounded, not by how much.

---


## Supplementary Result 21 | Lineage adjustment of the matched-family contrast

Registered in `NM_V4C_LINEAGE_ADJUSTMENT_AMENDMENT_002.json` before any arm was computed,
including the three outcomes the result could take and what each would require of the
manuscript. The outcome that occurred was the middle one, and the title changed because the
rule said it would.

### Cohort composition by sequence type

| host | genomes | sequence types | untypeable | effective STs (1/HHI) | largest ST | its share |
|---|---:|---:|---:|---:|---|---:|
| *A. baumannii* | 780 | 120 | 26 | **3.88** | ST2 | **48.2%** |
| *Klebsiella group* | 3,460 | 715 | 447 | **20.76** | ST11 | **17.3%** |

The asymmetry is the point. *A. baumannii* contributes 780 genomes but an effective **3.88**
lineages, because sequence type 2 — global clone GC2 — accounts for 48.2% of them. Any
estimate that treats those genomes as independent overstates its own precision.

### The registered arms

| arm | families | odds ratio | 95% CI | concordant |
|---|---:|---:|---|---:|
| baseline, as published | 58 | **50.29** | 45.61–55.45 | 56 |
| **L1** one genome per sequence type | 49 | **21.76** | 17.49–27.07 | 47 |
| **L1u** L1, untypeable genomes dropped | 47 | **21.33** | 16.37–27.80 | 45 |
| **L3** lineage-balanced weighting | 58 | **22.05** | 17.74–27.41 | 54 |

| influence arm | lineages | max relative change in ln(OR) | ceiling | within |
|---|---:|---:|---:|---|
| **L2** leave-one-ST-out | 1,308 | **0.2544** | 0.15 | **no** |
| **L2u** untypeable pooled | 837 | 0.2544 | 0.15 | **no** |

The lineage that breaches the ceiling is `AB|2`. Removing that one sequence type moves
ln(OR) by 25.4%, against a ceiling of 15% borrowed from the BioProject arm of the
original freeze. The ceiling was borrowed rather than chosen here, because choosing a
threshold after seeing the data is precisely what a freeze exists to prevent.

### What this does and does not establish

It establishes that the contrast is not an artefact of clonal replication: collapsing to one
genome per sequence type leaves **21.76-fold (95% CI 17.49–27.07)**, an interval that excludes
1 and a direction unchanged in every arm.

It also establishes that the crude estimate overstates the effect. The adjusted interval does
not overlap the unadjusted one, and roughly 43% of the crude odds ratio survives
adjustment. Both statements are true at once and the manuscript reports both.

### The dominant clone removed, and the interval the design actually supports

Two further arms were run under `NM_V4C_REFEREE_ANALYSES_AMENDMENT_004.json`, frozen before
either was computed.

**Sequence type 2 removed entirely.** Dropping every *A. baumannii* genome of the dominant
lineage, with *Klebsiella* unchanged, leaves **18.56 (95% CI 16.66–20.69)** across 55 families,
52 of them concordant. The contrast is smaller than the unadjusted 50.29 and is not a property
of that clone.

**The interval on the unadjusted estimate.** Every other interval in this work is a cluster
bootstrap; the published 50.29 carried an analytic Wald interval that treats 74,349
occurrences as independent. Resampling clusters instead:

| resampling unit | clusters | 95% CI |
|---|---:|---|
| none (analytic, Robins–Breslow–Greenland) | — | 45.61–55.45 |
| BioProject | 1,404 | 34.26–74.08 |
| **sequence-type lineage** | 1,308 | **11.93–130.52** |

The lineage-clustered interval is an order of magnitude wider than the analytic one. That is
not a defect of the bootstrap; it is what an effective 3.88 lineages in the extreme host buys.
The analytic interval should not have been the one reported, and is not.

Sequence type is a seven-locus surrogate. Two genomes sharing a type are not clones and two
differing at one locus may be near-identical, so this removes clonal replication rather than
controlling descent. A core-genome or SNP-distance analysis would be a stronger control and
was not performed.

---

## Supplementary Result 22 | Robustness and internal audits

**Blinded adjudication.** See Supplementary Method 1 for the full design. Design-weighted agreement
between the rule-engine state and the adjudication was **0.9920** (95% CI 0.9761–1.0000) on 120
stratified blocks, against a registered gate of ≥0.90 with a bootstrap lower bound ≥0.80. Raw
agreement was 119 of 120. This audit establishes rubric consistency for the MOBILE versus QUIESCENT discrimination only; it is not an independent validation, and it cannot establish consistency for the NON_EVALUABLE state and cannot produce an integron-specific estimate.

**Tool and database version invariance.** Re-running the plasmid mobility layer under a different
MOB-suite version and a different marker database produced **no class transitions**.

**Independent implementation.** CONJScan within MacSyFinder agrees with MOB-suite on the class E
definition at **κ = 0.875**; 61.24% of plasmid-borne occurrences sit on replicons that both tools
call conjugation-consistent. Neither tool is treated as a truth standard, so this is reported as
concordance rather than as accuracy.

**Cluster robustness.** The *A. baumannii*/*K. pneumoniae* block-weighted ratio is 2.4577 at
baseline and 2.1504 when collapsed to one genome per BioProject. Leave-one-BioProject-out across
1,195 projects moves it between 2.3224 and 2.4832; no single project moves it by more than 15%. For
the distance analysis, leave-one-BioProject-out across all 2,248 projects moves the primary contrast
at 1 kb by at most 0.0140 on a baseline of 0.4186.

**Independent re-derivation.** 35 checks re-derived the published structural quantities from raw
tool output; 15 further checks re-derived the denominator flow. Both returned zero disagreements.

## Supplementary Result 23 | Denominators, and why three correct numbers differ

Three plasmid-share figures appear in this study. All three are correct and they answer different
questions; the denominator must be quoted with the number.

- **52.736%** — occurrence-weighted. Counts genes. Describes the resistome.
- **36.586%** — genome-collapsed events, 3,569/9,755. Counts (genome, compartment) events, at most
one per genome per compartment. Describes how often a compartment is used at all.
- **35.932%** — the mean of per-genome plasmid percentages. A per-genome average, which weights a
genome carrying two resistance genes equally with one carrying forty.

Averaging the `pct_plasmid` column of the genome-level summary does not reproduce 36.586%, and is
not intended to.

The same distinction recurs on the chromosomal side. **46.39%** of chromosomal *occurrences* lie
within 10 kb of a marker, but only **30.14%** of chromosomal *blocks* do, because marker-positive
blocks carry 2.46 resistance genes on average against 1.23 for marker-negative blocks. Occurrence
weighting counts the same neighbourhood once per gene in it; block weighting counts it once. Both
are reported throughout, and the contrast between species is robust to the choice while the absolute
level is not.

## Supplementary Result 24 | The chromosomal mobile compartment, under four denominators

This section stood in the main text of V3. It moves here unchanged under the journal's
word limit; not one word is altered, and no claim it makes is withdrawn.

**46.39% of chromosomal acquired resistance-gene occurrences (16,303/35,140) lie within 10 kb of an
insertion sequence, transposase, integrase or integron marker** (Extended Data Fig. 3a). At ±5 kb the figure is
40.51%; 0.26% of occurrences overlap a marker directly. The marker inventory holds 32,364 features:
29,331 insertion sequence or transposase [@siguier2014; @siguier2006] and 3,033 integrase or
integron [@gillings2014; @escudero2015].

This figure is occurrence-weighted, and the weighting is not incidental. Computed over the 21,955
unique context blocks, 6,617 blocks (30.14%) contain at least one marker. The gap arises because
marker-positive blocks carry 2.46 resistance genes on average against 1.23 for marker-negative
blocks: mobile-element context and multi-gene context co-occur, so occurrence weighting counts the
same neighbourhood more than once. Both weightings preserve the host ordering (Extended Data Fig. 3b).

A one-block discrepancy inside that figure is worth stating. Blocks containing at least one marker
number 6,617; blocks containing at least one class-B occurrence number 6,616 (Extended Data Fig. 3c). Exactly one
block carries a marker lying more than 10 kb from every resistance-gene occurrence in it. The
occurrence-level counterpart is larger: 111 occurrences sit in a block containing a marker outside
their own window, and are censored rather than counted.

The ±10 kb figure is an operational classification window, not a biological boundary, and it is not
load-bearing. The *A. baumannii* excess over *Klebsiella* is present at 1, 2, 5 and 10 kb and is
**largest at 2 kb** (+0.4651) rather than at 10 kb (+0.3932); 93% of *A. baumannii* homology
detections fall within 2 kb; and the host ordering survives the stricter structural endpoint. Nearly
half of what a plasmid-fraction summary scores as stable is in mobile-associated context.
Chromosomal location is not fixity.

## Supplementary Result 25 | Re-analysis of Jia et al. 2026

Jia et al. conclude that AMR mobility is dictated by gene function rather than by the host
bacterium. Their deposited data (OSF `10.17605/OSF.IO/WE3TX`, MIT licence) were re-analysed
under **their own mobility definition** — plasmid, or genomic island predicted by
IslandViewer 4 — so that any difference is attributable to the comparison rather than to the
definition. Families required ≥5 occurrences in each host. The design was frozen and hashed
before any outcome column was inspected.

All six pairwise contrasts among the four ESKAPE Gram-negative genera are reported; none is
omitted. The exhaustive set is the safeguard against selection.

| contrast | families | MH odds ratio | 95% CI | *I*² | direction |
|---|---:|---:|---|---:|---|
| Klebsiella vs Pseudomonas | 27 | **16.2145** | 14.3369–18.3381 | 95.7% | 22 up / 5 down |
| Acinetobacter vs Pseudomonas | 25 | **2.7445** | 2.4566–3.0661 | 96.7% | 14 up / 11 down |
| Klebsiella vs Enterobacter | 38 | **0.5129** | 0.4464–0.5892 | 93.0% | 19 up / 19 down |
| Acinetobacter vs Enterobacter | 23 | **0.2244** | 0.1839–0.2739 | 89.2% | 9 up / 14 down |
| Acinetobacter vs Klebsiella | 25 | **0.1657** | 0.1505–0.1825 | 97.0% | 10 up / 15 down |
| Pseudomonas vs Enterobacter | 26 | **0.0775** | 0.0633–0.0950 | 83.8% | 2 up / 24 down |

**All six intervals exclude 1**, and the point estimates span a 209-fold range. That span is
computed across **non-identical strata**: each pairwise contrast is restricted to the families
present in both of its own hosts, so the six use between 23 and 38 families and are not nested.
The consequence is that they are not transitive — (*A* vs *E*)/(*K* vs *E*) is 0.44 against a
directly estimated *A* vs *K* of 0.17, and (*K* vs *E*)/(*P* vs *E*) is 6.6 against a directly
estimated 16.2. The argument these six support is that the contrast is host-dependent at all,
not that the 209-fold span is a single well-defined quantity. Under a
gene-intrinsic model they should sit near unity.

**No pooled estimate is reported.** *I*² lies between 83.8% and 97.0% in every contrast; a
pooled odds ratio would not describe its own strata. The primary species pair
(*A. baumannii* against *K. pneumoniae*) gives 0.1203 over 21 families with *I*² = 96.3%
and a 9/12 direction split, and dropping OXA β-lactamase moves it by 51.65% — three
independent reasons not to pool it.

The clearest single case is OXA β-lactamase, the largest family in that contrast at 1,957
occurrences: chromosomal and non-mobile in 75.8% of *A. baumannii* occurrences (811 of
1,070) and mobile in 94.8% of *K. pneumoniae* occurrences (841 of 887).

> **Limitations.** This re-analysis cannot corroborate the magnitude reported in the main
> text: the determinant caller (RGI with CARD, against AMRFinderPlus), the mobility
> definition and the outcome contrast all differ. It also carries no study-level control —
> the deposited metadata records BioProject PRJNA224116 for 66,865 of 66,889 replicons, the
> NCBI annotation umbrella rather than the submitting study — so clonal oversampling cannot
> be collapsed at study level in that dataset by us or by anyone.


## Supplementary Result 26 | The matching unit, and pooling under heterogeneity

Registered in `NM_V4C_MATCHING_UNIT_AMENDMENT_003.json` before either question was computed.

### What a gene family is here

It is not the AMRFinderPlus family field. A family is the `Element symbol` with a trailing
allele designator removed at the last hyphen, when that tail is a number, a short roman
numeral, or one to three letters:

```
head, _, tail = symbol.rpartition("-")
family = head if head and re.match(r"^(?:\d+[A-Za-z'\"]*|[IVX]+[a-z]?|[A-Za-z]{1,3})$", tail)
         else symbol
```

so `blaTEM-1` becomes `blaTEM`, `aac(6')-Ib3` becomes `aac(6')-Ib`, and `sul1` is unchanged.
The rule reproduces the `gene_family` column of the deposited occurrence table exactly: **0
mismatches across all 74,349 rows**. It is stated here because it was not stated anywhere
before, and because it is a modelling choice rather than a database fact.

The choice matters in one direction that works against us. Collapsing puts the intrinsic
chromosomal *bla*OXA-51-like enzyme and the transposon-borne *bla*OXA-23 into a single
`blaOXA` family. So the contrast was repeated with no collapsing at all.

| arm | matching unit | units | odds ratio | 95% CI | concordant |
|---|---|---:|---:|---|---:|
| baseline | gene family | 58 | **50.29** | 45.61–55.45 | 56 |
| **M1** | exact `Element symbol` | 68 | **56.40** | 50.37–63.15 | 65 |
| **M2** | exact symbol, one genome per ST | 58 | **23.18** | 18.04–29.78 | 55 |

The effect does not depend on the collapsing rule. At allele resolution it is larger, not
smaller, and under the lineage adjustment it lands where the family-matched analysis lands.

### Pooling under heterogeneity

The family-matched contrast carries *I*² = 83.6%. A fixed-effect Mantel–Haenszel summary
under that much heterogeneity describes the pooled odds, not a common effect, so a
random-effects estimate is reported beside it wherever the pooled value appears.

| pooling | cohort | strata | odds ratio | 95% CI | τ² | *I*² |
|---|---|---:|---:|---|---:|---:|
| fixed effect (MH) | full | 58 | **50.29** | 45.61–55.45 | — | 83.6% |
| random effects (DL) | full | 58 | **42.78** | 30.80–59.41 | 0.795 | 83.6% |
| random effects (DL) | one genome per ST | 49 | **19.50** | 14.39–26.43 | 0.263 | 30.2% |

**Most of the heterogeneity was clonal replication.** Adjusting for lineage takes *I*² from
83.6% to 30.2% and τ² from 0.795 to 0.263, and the random-effects and fixed-effect estimates
then agree closely (19.50 against the 21.76 of the same gene-family matching unit; the
exact-symbol arm is not comparable to it and an earlier revision paired the two). The families were never as different from one another
as the unadjusted *I*² suggested; they were unevenly sampled from the same few lineages.

DerSimonian–Laird is a moment estimator and underestimates between-stratum variance when
strata are few. It is reported because omitting it while reporting *I*² = 83.6% would be
the more misleading choice, not because it is the better estimator.

---


## Supplementary Result 27 | Provenance of the headline quantities

Each value below is traced to the file that holds it, the arithmetic that produces it, and the
script that computed it. Two entries carry an explicit provenance caveat rather than a clean chain;
both are stated rather than smoothed over.

| quantity | value | weighting or collapse rule | source artefact |
| --- | --- | --- | --- |
| genome-collapsed plasmid share | 36.586% (3,569 / 9,755) | one event per (genome, compartment); each genome contributes ≤1 plasmid and ≤1 chromosomal event | `genome_level_summary.tsv` |
| mean per-genome plasmid share | 35.932% (mean of the per-genome `pct_plasmid` column) | per-genome average; weights a 2-gene genome equally with a 40-gene genome | `genome_level_summary.tsv`, column `pct_plasmid` |
| occurrence-weighted plasmid share | 52.736% (39,209 / 74,349) | occurrence-weighted | `determinant_portability_classes.tsv` |
| chromosomal association, occurrence-weighted | 46.39% (16,303 / 35,140) | occurrence-weighted, own ±10 kb window | `arg_mge_neighbourhood.tsv` |
| chromosomal association, block-weighted | 30.14% (6,617 / 21,955) | one block, one count | `shared_context_blocks.tsv` |
| *A. baumannii* logit residual | 2.1126 (95% CI 1.8065–2.4376) (observed logit *M* 0.5447 − predicted −1.5678) | fit on 8 confirmation species; BioProject-within-species cluster bootstrap | `NMV4_RESULT_RECEIPT.json` → `T2` |
| *A. baumannii*/*K. pneumoniae* ratio, baseline | 2.4577 (0.6329 / 0.2575 block-weighted associations) | block-weighted, full cohort | `NMV2_RESULT_RECEIPT.json` line 7 |
| the same ratio, BioProject-collapsed | 2.1504 (0.5778 / 0.2687 after collapse) | one genome per BioProject (269 *A. baumannii*, 945 *K. pneumoniae*) | `NMV2_RESULT_RECEIPT.json` lines 35–43 |
| combined conjugation-consistent and multi-class cargo architecture | 1,455 of 6,621 (21.98%) (1,455 / 6,621) | conjugation-consistent AND ≥3 drug classes AND (metal OR virulence) on the same replicon | `plasmid_convergence.tsv` (6,621 rows) · `pr_context_convergence.py` lines 252–255 |
| conjugative − marker-negative, ≥3-class share | +19.86 pp (95% CI 14.16–25.41) (61.54% − 41.68%) | BioProject-clustered percentile bootstrap | `NMV3_RESULT_RECEIPT.json` → `C09.baseline` |
| MOB-suite / CONJScan agreement | κ = 0.875 (two-tool concordance on the class E definition) | replicon-level, both tools run on the same replicons | `NMV3_RESULT_RECEIPT.json` → arm I |

Producing scripts, in table order: `pr_context_build_join.py` for the first three rows and the two
chromosomal-association rows via `pr_context_neighbourhood.py`; `nmv4_run.py`, `nmv2_run.py`,
`pr_context_convergence.py` and `nmv3_score.py` for the remainder.

**Two provenance caveats, stated plainly.**

*First*, the two per-genome plasmid figures differ in the fifth decimal place depending on which
column is averaged: 35.932096% is the mean of the stored two-decimal `pct_plasmid` column, and
35.932101% the mean computed from raw counts. Neither is wrong and neither was changed; the
manuscript quotes 35.932%, which is common to both.

*Second*, the count 1,455 was never written to a hashed deliverable. The named output table holds
only the top 40 replicon-type/drug-class combinations, and the total was printed to standard output
and transcribed by hand into the report text. For this revision the stored Boolean predicate was
re-applied to the stored 6,621-row table and returns **1,455 (21.9755%)**, matching the published
value exactly. That is a re-derivation of an existing number, not a new analysis, and it closes the
gap; the underlying count remains absent from any hashed artefact and should be added to the deposit
before publication.

---

# Supplementary Tables

**Supplementary Table 1 | Evidence-layer inventory of the 184,538 determinant records.** Layer
counts are disjoint and sum exactly to 184,538. Only the first layer enters the primary denominator.

| evidence layer | records | in primary denominator | reason |
|---|---:|---|---|
| acquired AMR, core scope, non-efflux | **74,349** | **yes** | the frozen primary layer |
| metal / stress tolerance | 65,355 | no | metal tolerance is not antibiotic resistance |
| point mutation | 11,424 | no | not an acquired gene; detected only where an organism flag exists |
| virulence | 11,408 | no | virulence is not antibiotic resistance |
| acquired AMR, efflux class | 9,124 | no | sensitivity set S2 only |
| biocide tolerance | 5,544 | no | biocide tolerance is not antibiotic resistance |
| other stress (heat, acid) | 5,300 | no | not antibiotic resistance |
| acquired AMR, `plus` scope, non-efflux | 2,034 | no | sensitivity set S1 only |
| **total** | **184,538** | | |

The three acquired-AMR layers together number 85,507 records; removing the efflux (9,124) and
`plus`-scope (2,034) layers leaves the primary denominator of 74,349.

| primary denominator | occurrences |
|---|---:|
| plasmid-borne | 39,209 |
| chromosomal | 35,140 |
| unmatched, ambiguous, or missing coordinates | **0** |
| **total** | **74,349** |

**Supplementary Table 2 | Portability classes.**

| class | definition | occurrences |
|---|---|---:|
| A | chromosomal, no mobile-element marker within ±10 kb | 18,837 |
| B | chromosomal, at least one marker within ±10 kb | 16,303 |
| C | plasmid, no detected mobility marker | 7,170 |
| D | plasmid, relaxase detected | 6,043 |
| E | plasmid, relaxase and mating-pair formation detected | 25,996 |
| — E1 | relaxase and mating-pair formation, **no** detected origin of transfer (2,828 replicons) | 16,397 |
| — E2 | relaxase, mating-pair formation **and** a detected origin of transfer (1,109 replicons) | 9,599 |
| **total** | | **74,349** |

**Supplementary Table 3 | Weighted cumulative detection under the homology endpoint,
block-balanced.** Three groups appear here and in Supplementary Table 5 because three were registered:
`groups_and_contrasts` in `NMIS_FROZEN_PROTOCOL_V1.json` names *A. baumannii*, the *Klebsiella*
group and *P. aeruginosa*, and the two primary contrasts among them. The *Enterobacter* group's
occurrences are inside the 35,140 of Supplementary Table 4 and inside the background-null
sensitivities, but no distance-endpoint estimate was registered for it, so none is reported.

| group | n | *F*(1 kb) | *F*(2 kb) | *F*(5 kb) | *F*(10 kb) | RMD (bp) | median |
|---|---:|---:|---:|---:|---:|---:|---|
| *A. baumannii* | 8,005 | 0.5224 | 0.5878 | 0.5964 | 0.6323 | 4,180 | 647 bp |
| *Klebsiella* group | 15,568 | 0.1039 | 0.1227 | 0.1533 | 0.2390 | 8,365 | not reached |
| *P. aeruginosa* | 7,150 | 0.0993 | 0.1256 | 0.1410 | 0.1717 | 8,595 | not reached |

**Supplementary Table 4 | Structural census: element and occurrence accounting.**

| element category | n |
|---|---:|
| complete and fully contained | 7,923 |
| complete in shared block, crossing the occurrence-window boundary | 2,384 |
| partial, structurally incomplete or boundary-limited | 4,119 |
| **total elements resolved** | **14,426** |

| occurrence endpoint state | n |
|---|---:|
| complete and fully contained | 12,034 |
| complete in shared block, crossing the occurrence-window boundary | 716 |
| partial or boundary-limited | 1,902 |
| no structurally resolved element under the frozen definition | 20,488 |
| tool failure | 0 |
| **total chromosomal occurrences** | **35,140** |

**Supplementary Table 5 | Weighted cumulative detection under the structural endpoint,
block-balanced.**

| group | n | *F*(1 kb) | *F*(2 kb) | *F*(5 kb) | *F*(10 kb) | RMD (bp) | median |
|---|---:|---:|---:|---:|---:|---:|---|
| *A. baumannii* | 8,005 | 0.3531 | 0.3956 | 0.4315 | 0.4555 | 5,873 | not reached |
| *Klebsiella* group | 15,568 | 0.0633 | 0.0736 | 0.1021 | 0.1435 | 8,988 | not reached |
| *P. aeruginosa* | 7,150 | 0.0276 | 0.0493 | 0.0817 | 0.0954 | 9,281 | not reached |

**Supplementary Table 6 | Registered contrasts at four landmarks.** Cluster bootstrap over
BioProjects, B = 2,000, seed 20260822; 95% percentile intervals; Holm-corrected within each endpoint
family. All eight intervals in each panel exclude zero and all Holm-corrected *P* ≤ 0.002 (the
resolution floor, 2/B doubled by Holm across two contrasts).

| endpoint | contrast | 1 kb | 2 kb | 5 kb | 10 kb |
|---|---|---:|---:|---:|---:|
| homology | *A. baumannii* − *Klebsiella* | +0.4186 | +0.4651 | +0.4430 | +0.3932 |
| homology | *A. baumannii* − *P. aeruginosa* | +0.4231 | +0.4622 | +0.4554 | +0.4606 |
| structural | *A. baumannii* − *Klebsiella* | +0.2898 | +0.3221 | +0.3293 | +0.3119 |
| structural | *A. baumannii* − *P. aeruginosa* | +0.3255 | +0.3463 | +0.3498 | +0.3601 |

**Supplementary Table 7 | The two structurally positive, homology-negative occurrences.**

| assembly | replicon | determinant | species | structural distance |
|---|---|---|---|---:|
| GCF_015135855.1 | NZ_AP022446.1 | *bla*ACT-102 | *Enterobacter kobei* | 4,334 bp |
| GCF_009738085.1 | NZ_CP033102.1 | *oqxB* | *Enterobacter hormaechei* | 4,031 bp |

**Supplementary Table 8 | Retention of insertion-sequence elements through the structural gate, by
family.** The eight most abundant families in the window census, which are the eight plotted in
Extended Data Fig. 1b; families below IS*1380*'s 371 elements are not shown. `unassigned` is
ISEScan's `new` label, carried through unchanged: an element resolved structurally but not
matched to a named family profile. It is a label, not a family.

| family | elements called | complete after the structural gate | retention |
|---|---:|---:|---:|
| IS*6* | 4,848 | 4,385 | 90.4% |
| IS*30* | 384 | 342 | 89.1% |
| IS*5* | 1,278 | 1,089 | 85.2% |
| IS*4* | 2,655 | 1,869 | 70.4% |
| IS*91* | 1,176 | 631 | 53.7% |
| unassigned | 829 | 359 | 43.3% |
| IS*3* | 1,202 | 497 | 41.3% |
| IS*1380* | 371 | 0 | 0.0% |

**Why IS*6* carries two different counts, and why no ratio is taken between them.** The 4,848 above counts every IS*6* call inside the ±10 kb windows, complete or not, of which 4,385 pass the structural gate. Sensitivity BN3 reports 4,627 *complete* IS*6* elements genome-wide. Neither figure exceeds its genome-wide total — the genome-wide census holds 5,249 IS*6* calls of which 4,627 are complete — so the two are not arithmetically inconsistent, but they are also not two views of one count.

Restricting the genome-wide census to the same ±10 kb windows and counting directly gives **3,167** complete IS*6*, against the window census's 4,385. **The two ISEScan runs disagree on the same sequence.** They were run on different inputs — 22 kb extracted blocks in one, whole chromosomes in the other — and ISEScan's copy-number estimate, ORF context and family assignment all depend on what else is in the input. No quantity in this paper is computed across the two censuses, and the near-equality of 4,385 and 4,627 must not be read as a detection rate: the correct within-census statement is that **68.4% of complete genome-wide IS*6* (3,167 of 4,627) lies inside the 1.57% of chromosomal sequence covered by the windows**, against 8.2% for complete elements of all families (12,004 of 145,779).

**Supplementary Table 9 | Sensitivity analyses for the primary contrast at 1 kb** (*A. baumannii* −
*Klebsiella*, homology endpoint; baseline +0.4186).

| analysis | value |
|---|---:|
| primary, block-balanced | +0.4186 |
| occurrence-weighted | +0.4671 |
| deterministic one occurrence per block | +0.4209 |
| truncated blocks excluded (S4) | +0.4187 |
| circular-wrapped blocks excluded (S5) | +0.4191 |
| insertion sequence / transposase only (S6) | +0.4186 |
| integrase / integron only (SEC_INT) | +0.0100 |
| leave-one-BioProject-out, maximum absolute shift | 0.0140 |

**Supplementary Table 10 | Software and reference data versions.**

| component | version |
|---|---|
| AMRFinderPlus | 4.2.7 |
| AMRFinderPlus reference gene catalogue | 2026-08-07.1 |
| MOB-suite | 3.1.9 |
| MOB-suite database | 3.1.8 |
| ISEScan | 1.7.3 |
| HMMER | 3.3.2 |
| BLAST+ | 2.17.0 |
| FragGeneScan | 1.32 |
| Biopython | 1.88 |

---

# Supplementary Discussion

## Limitations of the structural endpoint

The structural endpoint is detection-limited in four specific ways, each of which sets a boundary on
what Supplementary Results 5 and 6 can support.

1. **Definition-limited.** Completeness is ISEScan's type `c` under one frozen parameterisation and
one HMMER version. A different caller, or a different profile set, would move the absolute levels.
2. **Window-limited.** A ±10 kb window truncates large elements by construction. This is why 2,384
complete elements that cross the occurrence-window boundary are excluded from the primary endpoint,
and why 716 occurrences are censored despite having a complete element in their shared block. Both
categories are reported rather than absorbed.
3. **Boundary-limited.** No insertion-sequence boundary was inferred and no flanking sequence was
retrieved beyond the window, so elements touching a boundary remain partial by construction.
4. **Horizon-limited.** No median is reached in any group under the structural endpoint, because
*F*(10 kb) < 0.5 everywhere. Restricted mean distance is reported instead, and no median is
extrapolated past the censoring horizon.

## Scope of the species claims

Three species groups were frozen for contrast: *A. baumannii*, the *Klebsiella* group and *P.
aeruginosa*. Every other species in the cohort is described but not contrasted. The *P. aeruginosa*
matched-family odds ratio of 1.26 (95% CI 1.11–1.43) is reported as qualified because, although the
interval excludes 1, the effect is small and the constituent families split 24 up against 22 down —
a pattern that is not consistent with a single dominant mechanism.

The three-architecture partition reproduces the direction of the host effect (odds ratio 16.54, 95%
CI 11.88–23.01) but one gene family carries 50.5% of the pooled weight. That sub-analysis is
therefore reported as partial, and its pre-registered weight-ceiling sub-claim is false rather than
merely unsupported.

## What co-location does and does not show

Extended Data Fig. 2 reports that conjugation-consistent replicons carry more drug classes and co-locate metal
tolerance more often. Co-location on a documented replicon is exactly that. It is not co-transfer:
no transfer was observed. It is not co-selection: no selective history was reconstructed. And it is
not evidence that conjugative architecture causes accumulation, since independent acquisition of
several determinants onto one replicon produces an identical genomic record. The association is
reported because it is a property of the deposited sequences; the mechanism behind it is not
addressed by this design.

## Generalisability

The cohort is a convenience sample of publicly deposited closed genomes. Closed assemblies are not a
random sample of clinical isolates: outbreak investigations, reference collections and
methodologically motivated sequencing are all over-represented, and the 2,248 BioProjects have an
effective count of 136.9. No fraction reported here is a prevalence estimate, and none should be
read as one. The species contrasts are internally valid comparisons within this cohort; their
external validity to any defined clinical population is not established by this design.
