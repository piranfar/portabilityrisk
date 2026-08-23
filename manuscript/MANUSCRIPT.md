# Portability, not presence: replicon-resolved evidence that plasmid fraction underestimates the mobile resistome

**PortabilityRisk — Paper 1. Manuscript draft V1, 2026-08-22.**

Assembled from frozen results only. Every quantitative statement below traces to a hashed
artefact registered in `PORTABILITYRISK_CANONICAL_ARTEFACT_MANIFEST_V2.tsv`; the mapping is
`PORTABILITYRISK_MANUSCRIPT_ASSET_INDEX_V2.tsv`. No number in this draft was computed while
writing it. Claim-by-claim permission is `PORTABILITYRISK_CLAIM_STATUS_MATRIX_V5.tsv`.

**Author.** Vahhab Piranfar.

---

## Abstract

Antimicrobial-resistance surveillance is built on gene-presence tables, which record whether a
resistance gene is there but not whether it can move. We assembled a replicon-resolved
portability dataset from 6,288 closed complete genomes of Gram-negative ESKAPE pathogens, in
which each of 74,349 acquired resistance-gene occurrences is assigned to a documented chromosome
or plasmid using NCBI molecule designations rather than prediction: 100 % resolution, 0
unmatched, 0 ambiguous, 0 missing coordinates. Every occurrence receives one of five
evidence-ranked portability classes combining documented location, plasmid mobilization
machinery, and chromosomal mobile-element context. To our knowledge, this is the first closed,
occurrence-level framework to resolve the portability of acquired resistance determinants across
both chromosomal and plasmid architectures using directly documented complete replicons.

We show for the first time that plasmid carriage and chromosomal mobile-element association are
distinct, largely orthogonal dimensions of portability; that bacterial host identity governs
which route a matched determinant takes; and that the low plasmid resistance-gene fraction of
*Acinetobacter baumannii* conceals an exceptionally mobile-element-associated chromosomal
resistome. One hundred and sixty-four gene families occur in both compartments, so portability
is a property of the occurrence rather than of the gene. Of chromosomal occurrences, 46.39 % lie
within 10 kb of a mobile-element marker, and that association is spatially tight rather than
diffuse: in *A. baumannii* the median distance to the nearest marker is 647 bp. Restricting the
endpoint to structurally resolved complete insertion sequences — a complete transposase open
reading frame with bilateral terminal inverted repeats, fully contained in the occurrence's own
window — reproduces the same host ordering across a 21,955-block census, and corroborates 73.80 %
of the homology-based calls.

We report no transfer event. The claim is about genomic architecture, not about demonstrated
mobilization.

---

## Introduction

Genomic surveillance of antimicrobial resistance answers one question well and a second question
badly. Whether a resistance determinant is present in an isolate is now routine to establish.
Whether that determinant is positioned to leave the isolate is not, and it is the second question
that governs consequence. A carbapenemase fixed in a chromosomal core spreads only as fast as its
host lineage. The same carbapenemase on a conjugative plasmid, or flanked by an intact insertion
sequence, is a different epidemiological object.

The gap has persisted for a structural reason rather than a conceptual one. Assigning a gene to a
chromosome or a plasmid requires knowing which molecule the contig carrying it belongs to, and in
draft assemblies that molecule is unknown. It must be predicted, by tools whose outputs carry an
error rate that propagates into every downstream statement about mobility. Large mobilome–resistome
association studies are consequently built on a location layer that is an estimate. Plasmid-centric
studies avoid the problem by making plasmids the unit of analysis, which is rigorous about one
compartment and silent about the other. Gene-presence resistome tables sidestep it entirely.

Closed complete genomes remove the estimate. For a closed assembly, NCBI states
`assigned_molecule_location_type` per molecule; a determinant whose coordinates fall inside that
molecule is on that molecule. This is documentation, not inference. It also imposes the boundary
of everything reported here: the results hold for closed genomes and cannot be extrapolated to
draft-assembly collections.

Because plasmid location is easy to summarise, plasmid fraction became the default proxy for
mobilization potential. Its failure mode is the subject of this paper. Plasmid fraction measures
one of two routes by which a resistance gene can be portable, and it is blind to the other:
chromosomal insertion in mobile-element context. If the two routes were correlated across hosts,
the proxy would be adequate. They are not.

We report a five-class, occurrence-level portability framework applied to 74,349 acquired
resistance-gene occurrences in 6,288 closed genomes; evidence that portability is a property of
the occurrence and not of the gene; a quantification of the chromosomal mobile compartment under
four denominators, a spatial distance distribution, and a structural insertion-sequence endpoint;
and an independently confirmed demonstration that the two portability routes are non-redundant
and can invert.

---

## Results

### R1 — A location layer with zero missingness

The cohort comprises 6,288 closed complete genomes of Gram-negative ESKAPE pathogens. Determinant
calling and denominator construction reduce 184,538 raw determinant records to a PRIMARY
denominator of **74,349 acquired resistance-gene occurrences** (Fig. 1). Each occurrence is joined
to the replicon whose coordinate interval contains it.

Every occurrence resolved: **39,209 plasmid-borne and 35,140 chromosomal, 0 unmatched, 0
ambiguous, 0 missing coordinates, 100.000 %**. Fifteen independent verification checks on the
denominator flow returned 0 disagreements.

Plasmid share is **52.736 %** occurrence-weighted. Collapsed to genome-level events it is
**36.586 %**, and that figure needs its denominator stated precisely, because two plausible
readings give different numbers. A *genome-collapsed event* is one event per (genome,
compartment): each of the 6,288 genomes contributes at most one plasmid event and at most one
chromosomal event, however many resistance genes it carries in that compartment. That yields
**3,569 plasmid events and 6,186 chromosomal events, 9,755 in total, and 3,569 / 9,755 =
36.586 %**. It is **not** the mean of per-genome plasmid percentages, which is 35.932 % — a
different statistic answering a different question.

Both figures are correct, and the divergence is the first of four denominator effects that recur
throughout this work. An occurrence-weighted figure describes the resistome; a genome-collapsed
figure describes how often a compartment is used at all.

### R2 — Five evidence-ranked portability classes

Each occurrence receives exactly one class (Fig. 2):

| class | evidence | n | % of 74,349 |
|---|---|---:|---:|
| A | chromosomal, no MGE marker within ±10 kb | 18,837 | 25.34 |
| B | chromosomal, ≥1 MGE marker within ±10 kb | 16,303 | 21.93 |
| C | plasmid, no mobility marker | 7,170 | 9.64 |
| D | plasmid, mobilization-consistent (relaxase) | 6,043 | 8.13 |
| E | plasmid, conjugation-consistent (relaxase and MPF) | 25,996 | 34.97 |

A + B + C + D + E = 74,349 = 100.000 %, reconciled independently.

Class E subdivides into a nested evidence tier: **E1 (16,397 occurrences on 2,828 replicons)**
carries relaxase and MPF; **E2 (9,599 occurrences on 1,109 replicons)** additionally carries a
detected *oriT*. E2 is 36.92 % of class E and 24.48 % of all plasmid-borne occurrences. The tiers
are nested, not competing: E2 is a subset of E1's evidential requirement plus one more marker.

Three separately labelled evidence layers — documented location, predicted plasmid mobility,
sequence-annotated chromosomal MGE context — are kept separable and ranked, so that a revision of
the mobility marker database moves C/D/E without touching A/B or the location layer.

### R3 — Portability is a property of the occurrence, not of the gene

**One hundred and sixty-four gene families occur in both compartments** (Fig. 3). Of the 158
families with ≥20 occurrences, 144 show a significant compartment preference at BH *q* < 0.05,
and **95 (66.0 %)** retain a Mantel–Haenszel interval excluding 1 after species adjustment.

The direction of these preferences is biologically correct and was recovered without any prior:
*qnr*, *bla*CMY, *sul3* and *tmexCD* fall out plasmid-side; *bla*PDC, *bla*ACT, *bla*LEN, *fosA*
and *oqxAB* fall out chromosome-side. That the pipeline reproduces known determinant biology blind
is the internal-validity argument for trusting it where the answer is not already known.

The 49 families that do not survive species adjustment are composition-driven, not refuted; 43 are
flagged lineage-dominated and are reported as flagged.

The consequence is structural. A gene-presence table cannot express this result, because the same
family name appears in both compartments and the table has no field in which the difference could
be written down.

### R4 — The chromosomal mobile compartment, under four denominators

**46.39 % of chromosomal acquired resistance-gene occurrences (16,303 / 35,140) lie within 10 kb
of an IS, transposase, integrase or integron marker** (Fig. 4). At ±5 kb the figure is 40.51 %;
0.26 % of occurrences overlap a marker directly. The marker inventory holds 32,364 features:
29,331 IS/transposase and 3,033 integrase/integron.

This figure is occurrence-weighted, and the weighting is not incidental. Computed over the 21,955
unique context blocks, **6,617 blocks (30.14 %) contain at least one MGE feature**. The gap arises
because marker-positive blocks carry 2.46 ARGs on average against 1.23 for marker-negative blocks,
a 2.01× ratio: mobile-element context and multi-gene context co-occur, so occurrence weighting
counts the same neighbourhood more than once.

A one-block discrepancy inside that figure is worth stating, because it is the block-level
signature of a real phenomenon rather than an error. Blocks containing ≥1 MGE feature number
6,617; blocks containing ≥1 class-B occurrence number **6,616**. Exactly one block carries a
marker that lies more than 10 kb from every resistance-gene occurrence in it. The
occurrence-level counterpart is larger: 111 occurrences sit in a block that contains a marker
which falls outside that occurrence's own ±10 kb window, and those are censored rather than
counted (Methods).

Both weightings are reported throughout. A single-number version of this result would be the
manuscript's largest referee risk, and the two numbers answer different questions: 46.39 % is the
share of the chromosomal resistome in mobile context, 30.14 % is the share of chromosomal
neighbourhoods that are mobile.

Nearly half of what a plasmid-fraction summary scores as "stable" is in demonstrably
mobile-associated context. Chromosomal location is not fixity.

### R5 — The association is short-range, and it is insertion sequences

Converting the ±10 kb threshold into a distance distribution changes the interpretation
qualitatively (Fig. 5). Using weighted cumulative detection *F*(*d*) with right-censoring at
10 kb and block-balanced weights:

| group | *F*(1 kb) | *F*(2 kb) | *F*(10 kb) | restricted mean distance | median |
|---|---:|---:|---:|---:|---|
| ***A. baumannii*** | **0.5224** | 0.5878 | 0.6323 | 4,180 bp | **647 bp** |
| *Klebsiella* group | 0.1039 | 0.1227 | 0.2390 | 8,365 bp | not reached |
| *P. aeruginosa* | 0.0993 | 0.1256 | 0.1717 | 8,595 bp | not reached |

A median of 647 bp means that in *A. baumannii* the typical acquired resistance gene is not merely
"in a neighbourhood containing" a mobile-element marker; it is adjacent to one.

The decomposition is what makes the result mechanistic rather than descriptive. An
**IS/transposase-only** endpoint reproduces the primary result. An **integrase/integron-only**
endpoint shows no *A. baumannii* excess at all (−0.0136 at 10 kb). The chromosomal mobile context
of this cohort is not a general "some MGE is nearby" effect. It is insertion sequences, close in,
in one host.

### R6 — The signature survives a structural endpoint at census scale

Homology markers establish that a transposase-like protein is nearby. They do not establish that
an intact element is there. We therefore ran a full structural reconstruction over all **21,955
context blocks** — 0 tool failures, 14,426 elements resolved — and restricted the endpoint to a
**structurally complete insertion sequence: a complete transposase open reading frame with
bilateral resolved terminal inverted repeats, fully contained within that occurrence's own ±10 kb
window** (Fig. 6).

Under this strictest available endpoint the host ordering is unchanged:

| group | *F*(1 kb) | *F*(2 kb) | *F*(10 kb) |
|---|---:|---:|---:|
| ***A. baumannii*** | **0.3531** | 0.3956 | 0.4555 |
| *Klebsiella* group | 0.0633 | 0.0736 | 0.1435 |
| *P. aeruginosa* | 0.0276 | 0.0493 | 0.0954 |

Both registered contrasts are positive with confidence intervals excluding zero at 1 kb and 2 kb
and at every landmark: *A. baumannii* − *Klebsiella* **+0.2898** at 1 kb and +0.3221 at 2 kb;
*A. baumannii* − *P. aeruginosa* **+0.3255** and **+0.3463**. All eight Holm-adjusted *p* = 0.001.

Across the whole chromosomal compartment, **12,034 of 35,140 occurrences** carry a structurally
complete, fully contained insertion sequence. **12,032 of those are class B**, so **73.80 % of
class-B occurrences (12,032 / 16,303) are structurally corroborated**, rising to **86.01 % in
*A. baumannii*** against 65.46 % in *Klebsiella* and 64.24 % in *P. aeruginosa*.

The difference between 12,034 and 12,032 is **exactly two occurrences that are structurally
positive while homology-negative**, and they are worth naming rather than rounding away:
`blaACT-102` in *Enterobacter kobei* (assembly GCF_015135855.1, replicon NZ_AP022446.1), with a
complete element at 4,334 bp, and `oqxB` in *Enterobacter hormaechei* (GCF_009738085.1,
NZ_CP033102.1), at 4,031 bp. Both are class A — homology-censored — so their transposases produced
no profile hit passing the frozen threshold, most plausibly because they are divergent or
unclassified transposases that the structural caller resolves from element architecture rather
than sequence similarity. Two consequences follow. The structural call set is *very nearly, but
not exactly*, a subset of the homology set, and in both discordant cases the structural method is
the more sensitive one. And neither occurrence falls in a headline species group — both are in
`other` — so the 73.80 % corroboration rate and every species contrast are unaffected.

**IS6 retains 90.4 % of its elements through the structural gate** (4,385 / 4,848), the highest of
any abundant family, while IS1380 retains none.

Short-range concentration sharpens rather than dilutes: **86.9 % of *A. baumannii* structural
detections fall within 2 kb**, against 51.7 % for *P. aeruginosa* and 51.3 % for *Klebsiella*.

The 26.20 % of class-B occurrences without a contained complete element are **not** reclassified.
A homology marker without a fully resolved element remains valid evidence of mobile-element
context under the frozen class definitions; the structural endpoint is an orthogonal, stricter
measurement, not a correction.

### R7 — Plasmid fraction and chromosomal MGE association are non-redundant

The two portability routes diverge across hosts (Fig. 7):

| species | plasmid share (occ.) | chromosomal MGE (occ.) | chromosomal MGE (block) |
|---|---:|---:|---:|
| *A. baumannii* | **13.96 %** | **80.91 %** | 63.29 % |
| *K. pneumoniae* | 67.23 % | 36.74 % | 25.75 % |
| *P. aeruginosa* | 12.28 % | 40.01 % | 17.65 % |
| *E. hormaechei* | 66.75 % | 24.01 % | 15.86 % |

*P. aeruginosa* is the decisive comparator. It has a plasmid share as low as *A. baumannii*
(12.28 % against 13.96 %) but a mid-range chromosomal MGE association. The *A. baumannii*
architecture is therefore not a mechanical consequence of low plasmid fraction.

Because the principle was discovered in *A. baumannii*, testing it on *A. baumannii* would be
circular. We therefore fitted the relationship between plasmid share and chromosomal MGE
association on **eight confirmation species only**, none of which was used to articulate the
principle, and then asked where *A. baumannii* falls. Its plasmid share lies inside the
confirmation range, so it is not an extrapolation. Its observed block-weighted MGE association is
**0.6329**; the value predicted from its plasmid share is **0.1725**. The residual is **2.1126 on
the logit scale, bootstrap CI 1.8065–2.4376, excluding zero**.

The fit itself is the second half of the argument. Across the eight confirmation species, plasmid
share explains almost none of the variance in chromosomal MGE association: in-sample *R*² =
**0.0275**, and the slope's bootstrap interval **includes zero** (−0.0588 to 0.2287). A low *R*²
here is not a weak result — it is the result. If the two axes measured one latent quantity, the
slope would be steep and the residual small. Both are the opposite.

Direct contrast of the two low-plasmid species gives a chromosomal MGE gap of **0.4548
(CI 0.4021–0.5039)**, excluding zero.

### R8 — Host identity governs which route a matched determinant takes

The species-level result could still be a composition artefact if *A. baumannii* simply carries
different gene families. It does not. Restricting to families present in both hosts and pooling by
Mantel–Haenszel across matched families, the odds that a determinant of the same family is
chromosomal-and-mobile rather than plasmid-borne is **50.29-fold higher in *A. baumannii* than in
*Klebsiella* (CI 45.61–55.45)** across 58 eligible families, with 56 of 58 families pointing the
same way. Collapsing to one genome per BioProject the estimate is **21.95 (CI 18.92–25.47)** over
49 families.

Against *P. aeruginosa* the same analysis gives **OR 1.26 (CI 1.11–1.43)**: the interval excludes
1, but the effect is small and families split 24 up against 22 down. This contrast is reported as
**qualified**. The host-conditioning result is strong against *Klebsiella* and weak against
*P. aeruginosa*, and the manuscript states it that way.

A three-architecture partition (chromosomal-quiescent 18,837, chromosomal-mobile 16,303,
plasmid-borne 39,209) reproduces the direction — *A. baumannii* versus *P. aeruginosa*, OR 16.54
(CI 11.88–23.01) — but one family accounts for 50.5 % of the pooled weight. That sub-analysis is
therefore reported as **partial**: the direction is consistent, the pooled magnitude is not
family-robust, and the pre-registered sub-claim that no family would exceed a weight ceiling is
**false**.

### R9 — Conjugation-consistent replicons carry the convergent cargo

Of 6,621 ARG-bearing plasmids typed, **3,937 are conjugation-consistent, 1,211
mobilization-consistent, and 1,473 marker-negative**. Cargo convergence tracks that architecture:
plasmids carrying ≥3 drug classes are **61.54 %** of conjugative replicons against 42.94 % of
mobilizable and 41.68 % of marker-negative ones (conjugative − marker-negative **+19.86
percentage points, BioProject-clustered CI 14.16–25.41**). Metal-resistance co-location follows
the same order (44.42 % / 27.99 % / 25.93 %), as does median ARG count (5 / 4 / 2). **1,455
replicons (21.98 %)** meet the high-concern architecture definition.

"Convergence" here means present on the same documented replicon. It is not co-transfer, not
co-selection, and not evidence that conjugative plasmids cause accumulation.

### R10 — Robustness

**The chromosomal MGE layer was validated by blinded expert adjudication.** On 120 stratified
blocks presented with methods de-identified, design-weighted agreement between the automated
MOBILE/QUIESCENT state and the adjudicator was **0.9920 (CI 0.9761–1.0000)**; 119 of 120 agreed,
with perfect agreement on QUIESCENT (60/60) and one disagreement on MOBILE (59/60). The scope
validated is MOBILE versus QUIESCENT and nothing wider. An earlier round of the same audit
**failed at 62/120** and is retained in the record: it diagnosed the instrument, not the
classifier — the casebook had presented evidence too thin to adjudicate — and the corrected
instrument produced the result above.

**Class assignments are invariant to tool and database version.** Re-running the plasmid mobility
layer under a different MOB-suite version and a different marker database produced **0.0000 %
class transitions**: C→C 1,473, D→D 1,211, E→E 3,937, and E1/E2 identical to the digit.

**An independent implementation corroborates class E.** CONJScan, which shares no code with
MOB-suite, agrees on the class E definition at **κ = 0.875** (relaxase κ = 0.883, MPF κ = 0.740);
**61.24 % of plasmid-borne occurrences sit on replicons both tools call conjugation-consistent.**
Neither tool is treated as a truth standard, and the comparison is reported as concordance.
CONJScan has no *oriT* model, so this arm covers class E fully and classes C and D only in part.

**The species contrast survives clustering and lineage collapse.** The *A. baumannii* /
*K. pneumoniae* block-weighted ratio is **2.4577** at baseline, 2.4590 collapsing to one genome
per BioSample and **2.1504** collapsing to one genome per BioProject. Leave-one-BioProject-out
across 1,195 BioProjects moves it between 2.3224 and 2.4832 — a maximum relative change of
**5.50 %**, with **no BioProject moving it by more than 15 %**. BioProject cluster bootstrap gives
2.4378 (CI 2.1565–2.8021) at full weighting and 2.1497 (CI 2.0286–2.2821) collapsed.

**The structural census was independently verified.** Thirty-five checks re-deriving the published
NMIS quantities from raw tool output returned **0 disagreements**; a further 43 checks on the
consolidated result tables returned 0.

---

## Discussion

**Portability is measurable, and the unit is the occurrence.** The central result of this work is
not a number but a change of unit. Once 164 gene families are observed in both compartments, a
per-gene mobility score is no longer well defined: it averages over occurrences that differ in
exactly the property being scored. Gene-level frameworks are the closest prior art, and this is
the limitation they carry; occurrence-level resolution is what removes it.

**The chromosomal mobile compartment is half the resistome and is under-instrumented.** Plasmid
biology has decades of dedicated method development. The chromosomal half — insertion sequences,
transposons, integrons in resistance context — is measured mostly as a by-product. Yet 46.39 % of
chromosomal occurrences here sit within 10 kb of a marker, in *A. baumannii* at a median of 647 bp,
and three quarters of those calls survive a structural endpoint requiring a complete element with
resolved terminal inverted repeats. Instrumenting only the plasmid compartment discards a large,
spatially tight, structurally real fraction of the mobile resistome.

**Two proxies, not one latent quantity.** The strongest evidence in this paper is negative:
plasmid share explains 2.75 % of the variance in chromosomal MGE association across the eight
confirmation species, and the slope relating them cannot be distinguished from zero. Against that
background the *A. baumannii* residual of 2.11 logits is not an outlier to be explained away; it
is what a genuinely two-dimensional space looks like when it is projected onto one axis. The
practical corollary is direct: ranking species, lineages or isolates by plasmid fraction alone
systematically misranks exactly the organisms whose resistance genes travel by the other route.

**Host, not gene, allocates the vehicle.** The matched-family analysis shows that the same
determinant family is routed differently depending on the host that carries it — 50-fold between
*A. baumannii* and *Klebsiella*. This positions the result against gene-level mobility scores: what
those scores average away is not noise but a host effect. The contrast against *P. aeruginosa* is
weak, and we do not claim it is general across all Gram-negative hosts.

**What this design cannot claim.** No transfer, conjugation or horizontal-gene-transfer event was
observed. Co-location on a replicon is not co-transfer, and independent acquisition produces an
identical genomic record. No transposition or element activity was assayed; a structurally complete
insertion sequence is a structure, not an event. Absence of a marker is a statement about a marker
database, not about biology: class C plasmids are not "non-mobilizable" and class A occurrences are
not "immobile". There is no phenotype, no MIC and no outcome data here, so no clinical risk claim
is available. The cohort is a convenience sample of public closed genomes and supports no
prevalence estimate.

**The closed-genome boundary.** Everything reported rests on documented replicon assignment, which
exists only for closed assemblies. That restriction is what gives the location layer zero
missingness, and it is equally what forbids extrapolating any of these fractions to draft-assembly
collections. Applying this framework at draft scale requires a prediction step, with an error rate
that must be measured — which is a different paper.

---

## Limitations

1. **Closed genomes only.** No result generalises to draft assemblies.
2. **Detection-limited endpoints.** The MGE layer is bounded by profile homology; the structural
   layer by ISEScan under a frozen definition, HMMER 3.3.2 profiles and a ±10 kb window. A shorter
   window truncates large elements, which is why 2,384 structurally complete elements were classed
   as crossing the occurrence window and excluded from the primary structural endpoint.
3. **Censored by construction.** 716 occurrences have a complete element in their block that
   crosses their own window; 111 occurrences have an in-block marker beyond their own window. Both
   are censored, reported separately, and not counted as positive.
4. **Medians unreachable under the structural endpoint** in every group, because *F*(10 kb) < 0.5
   everywhere; restricted mean distance is reported instead.
5. **The *P. aeruginosa* matched-family contrast is weak** (OR 1.26) and the three-architecture
   sub-analysis is family-dominated (one family, 50.5 % of weight). Both are labelled accordingly.
6. **Convenience sampling.** BioProject structure is uneven; every species-level estimate is
   accompanied by a clustered interval and a collapse analysis for that reason.
7. **No lineage typing.** cgMLST or clonal-complex control for the *A. baumannii* result is not
   included; BioProject collapse and leave-one-out are the available substitutes.

---

## Methods

**Cohort and frozen protocol.** 6,288 closed complete Gram-negative ESKAPE genomes. The context
protocol was frozen and hashed before any outcome column was read
(`FROZEN_PORTABILITY_CONTEXT_PROTOCOL_V1.json`,
`ddb155707dea1331b82d6ee5bf59ca2bf188fc8ad7bcb35ead2d7b08af8142ef`), as were the class definitions
(`frozen_portability_class_definitions.json`).

**Denominator construction.** 184,538 raw determinant records reduce to 85,507 and then to the
PRIMARY denominator of 74,349 acquired occurrences. Point mutations, intrinsic determinants and
`Scope=plus` records are excluded from PRIMARY; *mcr-9* is `Scope=plus` and is analysed only in
sensitivity set S1.

**Direct replicon join.** Each occurrence is assigned to the replicon whose coordinate interval
contains it, using NCBI `assigned_molecule_location_type`. No prediction step is involved. 0
unmatched, 0 ambiguous, 0 missing coordinates.

**Plasmid mobility typing.** MOB-suite relaxase, MPF and *oriT* models define classes C (no
marker), D (relaxase) and E (relaxase and MPF), with E2 the nested subset additionally carrying
*oriT*.

**Neighbourhood construction.** Chromosomal occurrences are merged into shared context blocks; each
block spans the union of the ±10 kb windows of the occurrences it contains. Circular replicons are
handled topologically, with wrapped coordinates resolved rather than clipped (57 wrapped blocks, 5
truncated). 21,955 blocks in total.

**MGE annotation.** Profile-homology annotation of IS/transposase and integrase/integron features
within blocks: 32,364 features (29,331 IS/transposase, 3,033 integrase/integron). An occurrence is
class B if ≥1 feature lies within its own ±10 kb window. Block-level marker positivity is a
different quantity and is reported separately: 6,617 blocks contain a feature; 6,616 contain a
class-B occurrence.

**Distance analysis.** Weighted cumulative detection *F*(*d*) with right-censoring at 10 kb;
block-balanced weights (1/*m* per occurrence, summing to 21,955); restricted mean distance computed
as E[min(*D*, *L*)]; BioProject cluster bootstrap, *B* = 2,000, seed 20260822, 2,248 BioProjects,
effective count by 1/HHI = 136.9. Holm correction across the registered contrast family.

**Structural insertion-sequence census.** ISEScan 1.7.3 under a pinned environment (HMMER 3.3.2,
BLAST 2.17.0, FragGeneScan 1.32) over all 21,955 blocks, one thread per block, no reruns. An
element is structurally complete if it is `type c` with a complete transposase ORF and bilateral
resolved terminal inverted repeats. The primary endpoint additionally requires **full containment
within the individual occurrence's own ±10 kb window**; elements complete in the shared block but
crossing that boundary form a separate, explicitly reported category. This containment rule was
registered as a numbered amendment, hashed before any benchmark block was scored.

**Statistics.** Family-matched contrasts use Mantel–Haenszel pooling with Cochran's *Q* and *I*²
reported; species-level contrasts use BioProject cluster bootstrap. Confirmation-species fits
exclude both discovery species. No outcome column was inspected to choose a feature, threshold,
species or hyperparameter.

**Validation programme.** Blinded expert adjudication of the MGE layer (n = 120, stratified);
lineage and BioProject robustness; tool-version, database-version and independent-implementation
arms; a discovery/confirmation split for the discordance principle; a distance analysis; and a
structural census. Each ran under its own protocol, frozen and hashed before outcomes were read,
with an independent verifier that re-derives the published numbers and reports its own
disagreement count.

**Reproducibility.** All frozen protocols, receipts, result tables and analysis code carry full
64-character SHA-256 digests and are versioned rather than overwritten. Superseded versions are
retained.

---

## Data and code availability

All inputs are public NCBI records and remain governed by NCBI's terms.

**Data.** The occurrence-level dataset — 74,349 acquired resistance-gene occurrences with their
replicon assignment and portability class — is deposited in Zenodo at
**https://doi.org/10.5281/zenodo.22065542** (version 1.0.0; concept DOI
https://doi.org/10.5281/zenodo.22065541 resolves to the latest version). The deposit includes a
full data dictionary defining every field, denominator and evidence-level dictionaries, JSON
Schemas for every table, a provenance map, pinned software and database versions, per-file SHA-256
checksums, and a standalone verification script that re-derives every headline denominator from the
archive alone and reports its own disagreement count.

Result and summary tables, frozen protocols, receipts and validation records are available now at
https://github.com/piranfar/portabilityrisk, which also carries the deposit's dictionaries,
schemas, checksums and access instructions. The occurrence-level dataset itself is not held in
that repository.

**Code.** Analysis code is at the same repository, versioned and hashed, with the exact commands.
Each module follows a freeze / score / verify structure: the protocol is hashed before any outcome
is read, the scorer aborts on a digest mismatch, and an independent verifier re-derives the
published numbers. Of 61 published scripts, 7 run against that repository alone; the remainder
require the Zenodo deposit and are published as the method of record.

Sequence caches and conda environments are not deposited; both are regenerable, from receipts
carrying accession.version with per-file digests and from the pinned environment lock.

**Licence.** The dataset and its metadata are licensed CC BY-NC-ND 4.0. Analysis code is licensed
Apache-2.0. Third-party data and software retain their own terms; no licence is applied to them by
this work.

**Outstanding before submission:** confirmation of journal-specific data and code availability
wording, which could not be retrieved from an official page.

## Licence

This preprint is made available under the Creative Commons Attribution-NonCommercial-NoDerivatives
4.0 International licence.

## Funding

This research received no specific grant from any funding agency in the public, commercial or
not-for-profit sectors. Cloud-computing resources were supported through promotional credits
provided by Oracle Cloud Infrastructure.

## Competing interests

The author declares no competing interests.

---

## Figures

All eight are drawn, in SVG, PNG and PDF, under
`docs/nature_microbiology/figures/manuscript/`. Figures 1–4 and 6–8 are produced by
`portabilityrisk_figures.py`, which verifies all fifteen input digests before drawing anything and
writes every output digest to `PORTABILITYRISK_FIGURE_RECEIPT_V1.json`. Figure 5 is the NM-DIST
four-panel figure, produced under the NM-DIST frozen protocol and reused unchanged.

| Figure | title | file stem |
|---|---|---|
| 1 | Cohort, denominator flow and the resolution guarantee | `figure1_cohort_and_resolution` |
| 2 | The five-class portability architecture, with E1/E2 nesting | `figure2_five_class_architecture` |
| 3 | Portability is a property of the occurrence, not the gene | `figure3_occurrence_not_gene` |
| 4 | The chromosomal mobile compartment under four denominators | `figure4_chromosomal_mobile_compartment` |
| 5 | Distance to the nearest marker, by host, with the IS/integron decomposition | `figure5_distance_to_mge` |
| 6 | The structural insertion-sequence endpoint at census scale | `figure6_structural_is_endpoint` |
| 7 | The discordance principle and its independent confirmation | `figure7_discordance_principle` |
| 8 | Conjugation-consistent replicons carry the convergent cargo | `figure8_convergence_by_mobility` |

Three quantities in the text were recomputed by the figure generator from the frozen tables, and
agree exactly: 144 families significant after BH, 95 surviving species adjustment, and the
6,617 / 6,616 block pair of Result R4.

Demoted to supplementary: plasmid replicon-type distribution; geography, source and clinical
status (irreducibly confounded with sequencing programme); evidence-hierarchy integration figures
belonging to a different manuscript's scope.

---

## Prohibited claims

Reproduced here so that no revision reintroduces one: demonstrated conjugation or mobilization ·
observed HGT or transfer · demonstrated transposition or element activity · causal co-selection ·
clinical risk prediction · population prevalence · "class C plasmids are non-mobilizable" · "class
A occurrences are immobile" · homology-only detections described as false positives · any
extrapolation of these fractions to draft assemblies · any *mcr-9* statement on the PRIMARY
denominator · any predictive-model result from a different paper in this programme.

Cohort size is never used as a novelty claim: several published closed-genome collections are
larger. The priority claims concern the framework and the propositions.
