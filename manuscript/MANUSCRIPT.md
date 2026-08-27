# Replicon-resolved portability of 74,349 acquired resistance-gene occurrences
across 6,288 closed Gram-negative ESKAPE genomes

**Vahhab Piranfar**^1,2^*

^1^ Independent Researcher, Jersey City, NJ, USA

^2^ Farname Inc, Ontario, Canada

\* Corresponding author: Vahhab Piranfar, Independent Researcher, Jersey City, NJ, USA.
Email: vahab.p@gmail.com. ORCID: 0000-0003-3653-5739

**Version 6.** Supersedes version 5. New: a lineage-balanced weighting of the background null, accessory-context arms, per-family intervals and a leverage diagnostic; the claim now rests on the host contrast, 1.69 to 4.27 across five adjustments, not on the absolute enrichment. Twelve corrections of fact and the removal of one unsupported count are listed in `PORTABILITYRISK_CHANGELOG_V5_TO_V6.md`.

---

## Abstract

Whether a resistance gene can leave its host depends on where it sits, and in draft assemblies that
location is predicted. We assembled an occurrence-level resource over **6,288 closed Gram-negative
ESKAPE genomes**: **74,349 acquired resistance-gene occurrences**, each assigned to a chromosome or
plasmid documented by NCBI, ranked into five evidence classes, with a genome-wide census of
**190,999 insertion sequences, 145,779 of them structurally complete**, across **6,190 ARG-bearing
chromosomes** and sequence types for the 4,240 genomes behind the primary contrast. Against a null
relocating each occurrence within its own chromosome, resistance genes lie closer to a complete
insertion sequence than chance allows in every group. *Acinetobacter baumannii* is extreme:
**1.7 to 4.3 times the *Klebsiella* group's enrichment** at 1 kb, a contrast that holds under adjustment for element density, local accessory context, composite-transposon structure, intrinsic determinants and clonal replication,
though the absolute enrichment does not.

Antimicrobial resistance was associated with an estimated 4.95 million deaths in 2019 [@murray2022],
and the greatest clinical burden falls on Gram-negative species [@who2024bppl] of the group
named for their capacity to escape antibiotic action [@rice2008]. Genomic surveillance answers one question about them well and a second badly. Whether a
determinant is present is routine to establish; whether it is positioned to leave the isolate is
not, and the second governs consequence. A carbapenemase fixed in a chromosomal core spreads only as fast as its host lineage;
the same carbapenemase on a conjugative plasmid, or flanked by an intact insertion sequence, is a
different epidemiological object [@partridge2018].

The gap has persisted for a structural reason rather than a conceptual one. Assigning a gene to a chromosome or a plasmid requires knowing which molecule its contig belongs
to, and in a draft assembly that molecule must be predicted, by tools whose error rate propagates
into every downstream statement about mobility [@arredondoalonso2017; @teixeira2025]; a controlled
benchmark places plasmid identification precision at 0.57 and insertion-sequence sensitivity at
0.58 [@kerkvliet2024]. Large mobilome–resistome studies are consequently built on a location layer that is an estimate
[@khedkar2022]. Plasmid-centric studies make plasmids the unit, rigorous about one compartment
and silent about the other [@redondo2020; @smillie2010; @conjsurvival2026]. Gene-presence tables sidestep the problem, and the limits of reading risk from inventories alone are under discussion [@larsson2026; @klumper2025]; the largest such analysis works from predicted genomic islands [@jia2026], not documented replicon location.

Closed complete genomes remove the estimate. Long-read and hybrid assembly resolves chromosomes and plasmids as separate circular molecules [@wick2017; @wick2023], and NCBI states the molecule type per replicon of a closed assembly [@kitts2016; @oleary2016]. A determinant whose coordinates fall inside that molecule is on it: documentation, not inference. That also bounds everything reported here — the results hold for closed genomes and cannot be extrapolated to draft-assembly collections.

Because plasmid location is easy to summarise, plasmid fraction became the default proxy for mobilization potential; its failure mode is this paper's subject. It measures one of two routes to portability and is blind to the other — chromosomal insertion in mobile-element context, across which exchange between compartments is structured rather than free [@wangdagan2021]. Frameworks that do quantify mobility aggregate contexts into a per-gene score [@ellabaan2022; @jia2026]; that would be harmless if the two routes were correlated across hosts. They are not.

In this study, **portability** denotes an evidence-ranked genomic property inferred from documented
replicon location, mobile-element-associated chromosomal context, structurally complete insertion
sequences and plasmid mobility-marker architecture. It does not denote observed transfer,
conjugation, transposition or element activity, none of which was measured here. We report a five-class, occurrence-level portability framework over 74,349 occurrences in 6,288 closed genomes; evidence that portability is a property of the occurrence, not the gene; the chromosomal mobile compartment under four denominators, a distance
distribution and a structural insertion-sequence endpoint; and a demonstration, fitted where the
pattern was not noticed, that the two routes are non-redundant.

---

## Results

### What the resource contains

In a closed genome that assignment is documented rather than predicted, and the resource is built on the distinction. Across **6,288 closed complete genomes** of Gram-negative ESKAPE
pathogens, **74,349 acquired resistance-gene occurrences** were each assigned to a documented
replicon, with 0 unmatched (Methods; Fig. 1; Supplementary Table 1). Every occurrence carries a five-class portability rank, a plasmid mobility call and, for the
35,140 chromosomal occurrences, the distance to the nearest marker inside its own ±10 kb window.

Three layers were built here and existed in no prior public record. A **genome-wide structural 
insertion-sequence census** annotated all **6,190 ARG-bearing chromosomes** under one pinned toolchain,
recovering **190,999 elements**, of which 145,779 meet a complete-structural definition: a complete transposase
open reading frame with bilateral resolved terminal inverted repeats. A **within-chromosome
permutation null** supplies the expectation against which those distances are read. **Sequence types**
were called for the 4,240 genomes contributing to the primary species contrast, which makes clonal structure visible.

The layer is deposited under CC BY 4.0, one row per file with its scale, its key columns and
the question it answers (Data availability; Supplementary Result 1). Of 7,216 assemblies
fetched, 6,288 met the eligibility rules; the count excluded by each individual rule was not
retained in any receipt and is not reconstructed here.

### Portability is a property of the occurrence, not the gene

**One hundred and sixty-four gene families occur in both compartments** (Fig. 2a). Of the 158
families with ≥20 occurrences, 144 show a significant compartment preference at Benjamini–Hochberg
*q* < 0.05 [@benjamini1995], and **95 (66.0%)** retain a Mantel–Haenszel interval excluding 1 after
species adjustment [@mantel1959] (Fig. 2b,c).

The direction of these preferences agrees with what is independently known of these
determinants, and was recovered without any prior:
*qnr*, *bla*CMY, *sul3* and *tmexCD* fall out plasmid-side; *bla*PDC, *bla*ACT, *bla*LEN, *fosA* and
*oqxAB* fall out chromosome-side. Reproducing known determinant biology blind is the internal-validity argument for trusting the pipeline where the answer is not known — weaker than it looks, since the same literature informs database curation, so the agreement is not fully independent of the annotation. The 49 families that do not survive species adjustment are composition-driven, not refuted.

The consequence is structural. A gene-presence table cannot express this result: the same family name appears in both compartments and the table has no field in which the difference could be written.

### Five evidence-ranked portability classes

Each occurrence receives exactly one class (Fig. 3a): **A**, chromosomal with no mobile-element
marker within ±10 kb (18,837); **B**, chromosomal with at least one marker (16,303); **C**, plasmid
with no mobility marker (7,170); **D**, plasmid carrying a relaxase (6,043); **E**, plasmid carrying
relaxase and mating-pair formation machinery (25,996). A + B + C + D + E = 74,349, reconciled
independently (Supplementary Table 2). Plasmid mobility was typed with MOB-suite [@robertson2018]
under the relaxase classification scheme [@garcillanbarcia2009].

Class E subdivides by mobility evidence, and the three evidence layers — documented
location, predicted plasmid mobility, sequence-annotated chromosomal context — are kept separable
and ranked, so a revision of the mobility marker database moves C/D/E without touching A/B or
the location layer (Fig. 3b; Supplementary Results 2 and 3).

### The association is short-range, structural, and specific to insertion sequences

Converting the ±10 kb threshold into a distance distribution changes the interpretation qualitatively. Under weighted cumulative detection *F*(*d*), right-censored at 10 kb with block-balanced weights, *A. baumannii* reaches *F*(1 kb) = 0.5224 against 0.1039 for the
*Klebsiella* group and 0.0993 for *P. aeruginosa*, with a **median distance of 647 bp** — not
reached within 10 kb in either comparator (Fig. 4d; Supplementary Table 3, Supplementary Result 4). The typical acquired resistance gene in *A. baumannii* is not merely near a marker; it is adjacent to one.

The decomposition narrows the result to one marker class; it stays observational, and proximity is not evidence of mobilisation. An insertion sequence–only
endpoint reproduces the primary result; an integrase/integron-only endpoint shows no *A. baumannii*
excess (Extended Data Fig. 1a; Supplementary Result 5). The chromosomal mobile context of this
cohort is insertion sequences, close in, in one host.

Homology markers establish that a transposase-like protein is nearby; they do not establish that an
intact element is present. We therefore applied a stricter sequence-structural endpoint,
reconstructing elements with ISEScan [@xie2017] over all **21,955 context blocks** — no failures, 14,426 elements resolved (Supplementary Table 4, Supplementary Result 6) — and restricted
the endpoint to a structurally complete insertion sequence: a complete transposase open reading
frame with bilateral resolved terminal inverted repeats, fully contained within the occurrence's own
±10 kb window.

Under this stricter endpoint the host ordering is unchanged (Fig. 4a; Supplementary Table 5): *A. baumannii* *F*(1 kb) =
0.3531 against 0.0633 and 0.0276. Both registered contrasts are positive with confidence intervals
excluding zero at every landmark — *A. baumannii* − *Klebsiella* +0.2898 at 1 kb, *A. baumannii* −
*P. aeruginosa* +0.3255 — all Holm-adjusted *P* ≤ 0.002 [@holm1979] (Fig. 4b;
Supplementary Table 6). Intervals here are 95% percentile intervals from a BioProject cluster bootstrap.

Across the chromosomal compartment, 12,034 of 35,140 occurrences carry a structurally complete,
fully contained insertion sequence; 12,032 of those are class B, so **73.80% of class-B occurrences
are structurally corroborated**, rising to **86.01% in *A. baumannii*** (Fig. 4c; Supplementary
Result 7). Two *Enterobacter* occurrences are structurally positive while homology-negative; neither falls in
a headline species group (Supplementary Table 7). IS*6*
retains 90.4% of its elements through the structural gate, the highest of any abundant family
(Extended Data Fig. 1b; Supplementary Table 8, Supplementary Result 8).

The 26.20% of class-B occurrences without a contained complete element are not reclassified: a homology marker is still evidence of context. The structural endpoint
is stricter, not a correction, and it is not independent — both run over the same blocks and
occurrences, and ISEScan itself uses profile homology. What it adds is a requirement for element architecture (Supplementary Result 6, Supplementary
Table 5).

### Resistance genes sit closer to intact insertion sequences than chance allows

Testing whether a host simply carries more insertion sequences requires a background, and none 
existed: elements had only ever been searched inside the ±10 kb windows, 1.57% of chromosomal sequence. Each occurrence was compared with a null relocating it uniformly
**within its own chromosome**, which preserves that chromosome's element density exactly, so
enrichment is density-normalised by construction. Design, permutation count, seed and thresholds
were registered before any species outcome existed.

**Resistance genes lie closer to a complete insertion sequence than chance allows in every group
tested** (Fig. 5). Within 1 kb *A. baumannii* detects at 0.579 against 0.034 expected — **16.91-fold**, outside the null 95% interval at *P* ≤ 0.0005, the 1/2001 floor; restricted mean distance 2,930 bp against 9,151 bp. **That absolute figure is
inflated by clonal sampling**: weighting each sequence type equally leaves 7.72-fold, and *Klebsiella* falls the same way, 6.29 to 1.81. The contrast between them does not change direction: 4.27 lineage-balanced against 2.69 unweighted (Supplementary
Result 9). All seven registered sensitivities support the enrichment, and the census reconciles against
the window-limited endpoint species by species (Supplementary Results 10 and 11).

**Four further explanations were tested and none accounts for the effect.** *Clonal replication*: 48.2%
of the *A. baumannii* genomes belong to sequence type 2, and 63% of its chromosomal occurrences
sit on those genomes; removing them leaves **12.74-fold**, still 2.03 times the *Klebsiella* group
(Supplementary Result 12). *Composite-transposon structure*: cargo of a composite element sits at near-zero distance by construction, but only 4,324 of the 33,822 chromosomal occurrences in the four groups (12.8%) are flanked by a same-family pair of complete elements; among *A. baumannii* occurrences that are not, enrichment is **13.58-fold**, 4.05 times the *Klebsiella* group — higher than unstratified, not lower. *Intrinsic determinants*: species-core genes are excluded from the
matched analyses by the registered eligibility rule, and removing the intrinsic *bla*OXA alleles at allele level does not reduce any contrast
(Supplementary Result 13). *Accessory-genome location*: the null may place an occurrence anywhere on its chromosome, mostly in element-sparse core. Confining relocation to its own ±100 kb neighbourhood leaves **7.19-fold**, and ±50 kb **4.38-fold**; both exclude the null, with the contrast at 2.18 and 1.69. Outside multi-gene resistance islands the contrast rises to 3.56 — *Klebsiella*'s enrichment is largely an island phenomenon and *A. baumannii*'s is not (Supplementary Result 14).

### A species contrast that the same data do not support

The resource also permits a matched-family contrast: among gene families present in both hosts,
the odds that a determinant is chromosomal-and-mobile **rather than plasmid-borne**. Pooled by
Mantel–Haenszel it is **50.29-fold higher in *A. baumannii* than in *Klebsiella*** across 58 families, 56 of which point the same way and 44 of which
exclude 1 individually (Supplementary Results 15 and 16). We report it, and that it does not carry the weight the number suggests.

Three things qualify it. Between-family heterogeneity is high (*I*² = 83.6%), so the fixed-effect
summary is reported with its random-effects companion, **42.78 (95% CI 30.80–59.41)**. The analytic
interval assumes 74,349 independent occurrences; resampling lineages instead gives **11.93–130.52**, an
order of magnitude wider, because *A. baumannii* carries an effective **3.88** lineages. And the
endpoint multiplies two effects — how often a host puts determinants on plasmids at all, and how
often a chromosomal determinant sits near a marker — so a host extreme on both produces a large odds
ratio partly by construction.

Separating them does not work here. The within-chromosome contrast, class B against class A,
fails in every variant the weight-concentration gate registered in the original frozen design —
one family carries 64% of the pooled weight against a 30% ceiling, only 11 of 49 families are
informative, and the sign moves between 0.48 and 6.11 with the handling of *bla*OXA alone. The
marginal rates are stable, 80.9% against 34.2%, but standardisation attributes the whole gap to
which families each host carries, not how shared families sit. **The within-chromosome question is not estimable at this project's registered
standard** (Supplementary Result 17).

What does survive is narrower and worth stating plainly: the contrast keeps its direction under
lineage adjustment (**21.76**, 17.49–27.07), with the dominant clone removed (**18.56**, 16.66–20.69), and with
*bla*OXA and *bla*SHV removed outright (**51.53**). *A. baumannii* keeps acquired resistance on the
chromosome where *Klebsiella* keeps it on plasmids, and plasmid fraction alone does not see that.

### Plasmid fraction and chromosomal mobile-element association are non-redundant

The two portability routes diverge across hosts. *A. baumannii* has one of the cohort's two lowest plasmid shares (13.96%) and its highest chromosomal mobile-element association (80.91% occurrence-weighted);
*K. pneumoniae* is its mirror (67.23%, 36.74%); *E. hormaechei* resembles *K. pneumoniae* (66.75%, 24.01%; Supplementary Result 18). *P. aeruginosa* is the decisive comparator: its plasmid share is lower still (12.28%) but its chromosomal association is mid-range (40.01%), so the *A. baumannii* architecture
is not a mechanical consequence of low plasmid fraction.

Because the pattern was noticed in *A. baumannii*, testing it there would be circular. We fitted the
relationship between plasmid share and chromosomal mobile-element association on **eight
confirmation species only**, excluding both species in which it was noticed, and then asked where
*A. baumannii* falls (Fig. 6a). Its plasmid share lies inside the confirmation range, so this is
interpolation rather than extrapolation. Its observed block-weighted association is 0.6329; the
value predicted from its plasmid share is 0.1725. The residual is **2.1126 on the logit scale (95%
bootstrap CI 1.8065–2.4376)** (Fig. 6c).

The fit itself is the second half of the argument. Across the eight confirmation species, plasmid
share explains almost none of the variance in chromosomal association: in-sample *R*² = **0.0275**, and the slope's bootstrap interval includes zero (−0.0588 to 0.2287);
one species carries 84% of the leverage and the slope doubles without it, so the slope is not an
estimate here (Supplementary Result 19) (Fig. 6b). A low *R*² is what a two-axis structure would produce, but with eight species it shows non-redundancy within this panel, not independence in general; a slope indistinguishable from zero is not evidence that it is zero. One latent quantity would give a steep slope and a small residual; both are the opposite. Direct contrast of the two low-plasmid
species gives a chromosomal association gap of 0.4548 (95% CI 0.4021–0.5039). The fit, the residual
and every diagnostic behind it are reported in full in Supplementary Result 18.

### Conjugation-consistent replicons carry the convergent cargo

Of 6,621 resistance-gene-bearing plasmids typed, **3,937 are conjugation-consistent, 1,211
mobilization-consistent and 1,473 marker-negative** (Extended Data Fig. 2a). Cargo convergence tracks that architecture: plasmids carrying ≥3 drug classes are 61.54% of conjugative replicons against 42.94% of mobilizable and 41.68% of marker-negative ones (+19.86 percentage points, 95% BioProject-clustered CI 14.16–25.41) (Extended Data Fig. 2b,d). **That difference is confounded with plasmid length**: conjugative replicons are larger, and larger replicons carry more of everything. Length stratification leaves +8.10 points and the strata disagree in sign, so the crude figure is not a size-independent effect (Supplementary Result 20). Metal-resistance co-location follows the same order
(44.42%/27.99%/25.93%) (Extended Data Fig. 2c) [@cooccurrence2026], as does median resistance-gene
count (5/4/2), and 1,455 of 6,621 replicons (21.98%) combine conjugation-consistency, three or more
drug classes and a metal or virulence determinant. Convergence here means present on the same documented replicon: not co-transfer, co-selection, or cause.

### Verification, and what the resource cannot be used for

The chromosomal mobile-element layer was audited for rubric consistency by blinded manual
adjudication. **This is an internal audit and not a validation.**
**The adjudicator was the author**, working from panels with the methods de-identified and blind to
species, study, gene identity and prior classification. The frozen design had required an
adjudicator not involved in building the pipeline; none was recruited, so that criterion is
recorded as **not met** and this audit is internal rather than independent. On 120 stratified blocks, design-weighted agreement between the automated state and the adjudicator was **0.9920 (95% CI 0.9761–1.0000)** against a registered gate
of ≥ 0.90, with one disagreement, and it establishes rubric consistency for the MOBILE versus QUIESCENT
discrimination only. It would become a validation if any uninvolved reader re-adjudicated the
same package, which is frozen and hashed for that purpose. An earlier round failed at 62/120 and stands as a FAIL: the annotation path under test was
visible to the adjudicator but excluded from the rule engine. The corrected audit removed that
asymmetry and drew fresh cases (Supplementary Method 1).

Class assignments are invariant to tool and database version: re-running the plasmid mobility
layer under a different MOB-suite version and marker database produced **no class transitions**,
and CONJScan [@cury2020; @cury2017] within MacSyFinder [@neron2023; @abby2014] agrees on class E
at **κ = 0.875**, neither tool a truth standard. Leave-one-BioProject-out across the 1,195 projects behind the *A. baumannii*/*K. pneumoniae*
ratio moves it by at most 5.50% against a prespecified 15% ceiling. **BioProject balancing addresses project-level sampling structure but is not clonal-lineage adjustment**: a BioProject is a submission unit, a clone a descent unit. We therefore typed every genome behind the headline contrast. *A. baumannii* is dominated by one lineage — 780 genomes in 120 sequence types and 26 untypeable, but sequence type 2 alone holds **48.2%** and the effective number of lineages is **3.88** by inverse Herfindahl–Hirschman index; the *Klebsiella* group is broader at 3,460 genomes, 715 sequence types and an effective **20.76**.

**The contrast survives lineage adjustment but is roughly halved by it.** Collapsing to one genome
per sequence type gives an odds ratio of **21.76 (95% CI 17.49–27.07)** across 49 families, 47 of them
concordant; dropping untypeable genomes gives 21.33 (16.37–27.80) and lineage-balanced weighting 22.05
(17.74–27.41). Every arm keeps its direction with an interval excluding 1, and the adjusted interval
does not overlap the crude one. Leave-one-sequence-type-out across 1,308 lineages moves ln(OR) by up
to **25.4%**, above the 15% ceiling borrowed from the BioProject arm, and the lineage responsible is
sequence type 2. The effect is therefore real and lineage-sensitive, and is reported as an
association with host rather than as a host effect (Supplementary Result 21).

The complementary objection, that these determinants might be intrinsic chromosomal genes rather
than acquired ones, does not hold. Six of the determinants most often named as intrinsic —
*fosA*, *oqxA*, *oqxB*, *bla*PDC, *bla*ACT and *bla*LEN — are absent from the matched set because a species-core gene
cannot meet the registered three-species eligibility rule. Two eligible families do carry intrinsic
members, and removing them changes little: excluding the *bla*OXA-51-like alleles at allele level
gives 48.45 (43.94–53.41), and excluding the *bla*OXA and *bla*SHV families outright gives 51.53
(46.28–57.37), above the unadjusted estimate rather than below it (Supplementary Result 13). The verification programme is reported in Supplementary Result 22.

---

The chromosomal mobile-compartment estimate depends on its denominator, and the host ordering is
preserved under each of them (Extended Data Fig. 3); all four denominators and their reconciliation
are reported in Supplementary Results 23 and 24.

## Discussion

**Portability-relevant genomic architecture is measurable, and the unit is the occurrence.** Once 164 gene families are observed
in both compartments, a per-gene mobility score is not well defined: it averages over occurrences
differing in exactly the property being scored. Gene-level frameworks carry that limitation
[@ellabaan2022], as do studies that resolve one family but predict its location [@liu2026] and
plasmid-only surveys that never allocate a chromosomal occurrence [@hou2026; @coluzzi2025;
@ecoli9700]. Occurrence-level resolution removes it.

**The chromosomal compartment is 47.26% of the resistome — 35,140 of 74,349 occurrences — and its
mobile fraction is under-instrumented.** Plasmid biology has decades of dedicated method development [@smillie2010; @redondo2020]; the
chromosomal half is mostly a by-product [@efaecium2026]. Yet 46.39% of chromosomal occurrences sit within 10 kb of a marker, in
*A. baumannii* at a median of 647 bp, and roughly three quarters survive a stricter structural
endpoint. That *A. baumannii* carries chromosomal resistance islands is long established
[@fournier2006; @hamidian2018; @nigro2012; @nigro2013]; new here is the multi-species quantification
showing this architecture is invisible to the plasmid-based proxy.

**Two proxies, not one latent quantity, and allocation is host-associated.** The strongest
evidence here is negative: across eight confirmation species plasmid share explains 2.75% of the
variance in chromosomal association and the slope cannot be distinguished from zero, on a fit
where one species carries 84% of the leverage. The *A. baumannii* residual of 2.11 logits is
large, which is what a two-dimensional space projected onto one axis produces; a mechanistic
duality between routes has been reported for carbapenem resistance [@carbapenemMobilome], and
draft-assembly work reaches related conclusions by prediction rather than documentation
[@fullassembly2026]. Ranking by plasmid fraction alone omits a distinct chromosomal-context axis,
so an organism whose resistance genes sit in that context is ranked as though it had none. No
external true ranking of portability was defined, and none is claimed here.

**What this design cannot claim, and where it stops.** No transfer event was observed; co-location
is not co-transfer, and independent acquisition produces an identical record. No transposition was
assayed: a complete insertion sequence is a structure, not an event. Absence of a marker is a
statement about a database, so class C plasmids are not non-mobilizable and class A occurrences are
not immobile. There is no phenotype and no outcome data, so no clinical risk claim is available, and
the cohort supports no prevalence estimate. Everything rests on documented replicon assignment,
which exists only for closed assemblies: that gives zero missingness and equally forbids
extrapolation to draft-assembly collections [@closed16622; @fullassembly2026].

The chromosome-wide density objection is answered on its own terms: a null relocating each
occurrence within its own chromosome preserves that chromosome's insertion-sequence density
exactly, leaves the *A. baumannii* contrast intact, and is supported by every registered
sensitivity. A uniform relocation null preserves neither coding position, accessory status, gene density, recombination hotspots nor resistance-island structure. Two were tested directly: relocating within the occurrence's own ±50–100 kb neighbourhood absorbs most of the enrichment, 16.91 to 4.38–7.19, without abolishing it, and the contrast survives outside resistance islands. The rest remain unaddressed, and no core-genome alignment exists here from which to build a true accessory-restricted null. It cannot establish causation — host species is confounded with lineage,
clinical context and sampling source, and proximity to a complete element is not evidence that it
has moved. The genome-wide census also differs from the frozen window-limited endpoint by a net 1,204 occurrences — 1,538 gains against 334 losses, and three quarters of the gains (1,154) are elements the extracted windows previously truncated, a limitation the earlier protocol declared. That window-limited endpoint remains this work's frozen
structural measure; the genome-wide census is the background against which it is normalised, not a
replacement.

A recent large-scale analysis reached the opposite conclusion, holding that AMR mobility is "largely
dictated by gene function rather than the host bacterium" across 39,089 RefSeq accessions
[@jia2026]. That holds pooled across thousands of species but not within their own data: re-analysed under
their own mobility definition, all six pairwise contrasts among the four ESKAPE Gram-negative
genera exclude 1 and span 0.078 to 16.2 (Supplementary Result 25), where a gene-intrinsic model
predicts unity. OXA β-lactamase is non-mobile in 75.8% of *A. baumannii* occurrences and mobile in 94.8% of
*K. pneumoniae* ones. The divergence is instrumental: genomic-island prediction "will not predict
very small mobile genetic entities" [@jia2026], and an ISAba1-flanked gene is far below that
scale. The two analyses measure different objects and neither refutes the other.

## Methods

### Definitions used throughout

Four terms carry the analysis and are used in one sense only.

An **occurrence** is one qualifying determinant record at one coordinate interval on one replicon of
one assembly. It is the unit of every primary count.

A **compartment** is chromosome or plasmid, taken from the molecule type NCBI documents for the
replicon. It is never predicted.

A **context block** is the union of the ±10 kb windows of a set of chromosomal occurrences that are
within 20 kb of one another. It is the unit of block weighting and the unit that nests within a
BioProject.

**Portability** is an evidence-ranked genomic property inferred from documented replicon location,
mobile-element-associated chromosomal context, structurally complete insertion sequences and plasmid
mobility-marker architecture, expressed as one of five ranked classes. It does not denote observed
transfer, conjugation, transposition or element activity, none of which was measured. It is not
transfer capability, and it is not a rate.

Three plasmid-share figures appear and are not interchangeable: 52.736% is occurrence-weighted and
counts genes; 36.586% is genome-collapsed events and counts (genome, compartment) pairs; 35.932% is
the mean of per-genome percentages. Each is reported with its denominator throughout, and the same
distinction separates the 46.39% occurrence-weighted and 30.14% block-weighted chromosomal figures
(Supplementary Method 2).

### Cohort and frozen protocol

The cohort comprises 6,288 closed complete genome assemblies of Gram-negative ESKAPE pathogens
retrieved from NCBI [@kitts2016; @oleary2016] with accession.version recorded and a SHA-256 digest
stored per downloaded file. An assembly was eligible if its assembly level was "Complete Genome", if
every replicon carried an `assigned_molecule_location_type` of Chromosome or Plasmid, and if the
organism fell within the frozen ESKAPE Gram-negative scope. No assembly was excluded after any
outcome column was inspected (Supplementary Method 3).

The analysis protocol and the class definitions were written and hashed **before any outcome column
was read**, and are published with the deposit. Every subsequent module follows the same three-stage
structure: a freeze stage writes a protocol and hashes it; a scoring stage verifies those digests
and aborts on mismatch; an independent verification stage re-derives the published quantities and
reports its own disagreement count. Independence here means a separate code path that does not
import the primary analysis; it was written and run by the same investigator and is not
independent scientific replication. Changes after a freeze are numbered amendments, themselves
hashed before use. No feature, threshold, species grouping or hyperparameter was chosen by
inspecting an outcome.

### Determinant calling and denominator construction

Acquired resistance determinants were called with AMRFinderPlus 4.2.7, reference gene catalogue
version 2026-08-07.1 [@feldgarden2021; @feldgarden2022]. The tool returned 184,538 records across
the cohort. The primary denominator retains only records of type AMR with scope "core" and excludes,
in this order: point mutations and disrupted point mutations, which are not acquired genes and were
searched only where an organism flag existed; efflux-class records, retained as sensitivity set S2;
`Scope=plus` non-efflux records, retained as sensitivity set S1; and the metal, biocide, stress and
virulence layers, which are not antibiotic resistance. This leaves **74,349 acquired resistance-gene
occurrences**. Each exclusion is arithmetic and prespecified; the layer inventory with row counts is
published with the deposit (Supplementary Method 4, Supplementary Table 1). *mcr-9* is `Scope=plus`
and therefore appears only in S1, never in a primary result.

**The matching unit is the gene family, and it is not the AMRFinderPlus family field.** A family is
derived from the `Element symbol` by removing a trailing allele designator at the last hyphen
when that tail is a number, a short roman numeral or one to three letters, so `blaTEM-1` becomes
`blaTEM` and `sul1` is unchanged. The rule exists because the question is whether the same
determinant takes different routes in different hosts, and at allele resolution most alleles occur
in one host only. It is also a choice that could merge determinants of genuinely different
mobility, so it was tested rather than assumed: repeating the contrast on the raw symbol, with no
collapsing at all, gives **56.40 (95% CI 50.37–63.15)** over 68 symbols, and the same analysis with one
genome per sequence type gives 23.18 (18.04–29.78) (Supplementary Result 26).

One occurrence is one qualifying record, identified by assembly accession, sequence accession and
coordinate interval. This composite is unique across the cohort: the deposited dataset has zero
duplicate keys.

### Direct replicon assignment

Each occurrence was assigned to the replicon whose coordinate interval contains it. Replicon
metadata was retrieved from the NCBI Datasets v2alpha `sequence_reports` endpoint (run receipt
2026-08-21T00:57:40Z; 7,216 assemblies, 23,129 replicon records, zero fetch failures), and the
molecule assignment is the `assigned_molecule_location_type` field, accepted only where it reads
exactly `Chromosome` or `Plasmid`. A direct location call requires all four of: the genome is in the
frozen complete-genome cohort; the AMRFinderPlus contig identifier maps exactly to one documented
replicon record; that replicon is explicitly designated Chromosome or Plasmid; and the gene
coordinates lie within that sequence's stated length. Matching tries the full accession.version
first and falls back to the accession with the version suffix stripped, recording which rule fired
in a per-row `join_method` field. Containment is 1-based and inclusive at both bounds, evaluated on
the minimum and maximum of Start and Stop so that reversed coordinates cannot cause a spurious
failure. The replicon unit is a unique sequence accession **within a versioned assembly**, so a
replicon accession is never matched across assemblies; most assemblies carry several plasmid
replicons. A replicon record found but carrying an empty molecule type becomes unclassified and
unresolved — it is never defaulted to chromosome — and an unmatched identifier is never inferred to
be chromosomal. No prediction step is involved and no plasmid-classification tool contributes to the
location layer. Within the closed-genome cohort the join produced zero unmatched, zero ambiguous and
zero missing-coordinate assignments: 39,209 plasmid and 35,140 chromosome. Location evidence is
recorded per occurrence in a dedicated field so that the layer can be audited separately from the
mobility and context layers.

**Denominator flow and completeness of assignment.**

The cohort comprises 6,288 closed complete genomes of Gram-negative ESKAPE pathogens. Determinant
calling with AMRFinderPlus [@feldgarden2019; @feldgarden2021; @feldgarden2022] and denominator
construction reduce 184,538 raw determinant records to a primary denominator of **74,349 acquired
resistance-gene occurrences** (Extended Data Fig. 4a). Each occurrence is joined to the replicon whose coordinate
interval contains it.

Every occurrence received a documented assignment: **39,209 plasmid-borne and 35,140 chromosomal, 0
unmatched, 0 ambiguous, 0 missing coordinates within the closed-genome cohort** (Extended Data Fig. 4b). This is
complete computational assignment to a documented replicon, not biological certainty about any
individual locus. Fifteen independent verification checks on the denominator flow returned no
disagreements.

Plasmid share is **52.736%** occurrence-weighted. Collapsed to genome-level events it is
**36.586%**, and that figure needs its denominator stated precisely. A genome-collapsed event is one
event per (genome, compartment): each genome contributes at most one plasmid event and at most one
chromosomal event, however many resistance genes it carries in that compartment. That yields 3,569
plasmid and 6,186 chromosomal events, 9,755 in total, and 3,569/9,755 = 36.586% (Extended Data Fig. 4c). It is a
count of events after collapse, not the mean of per-genome percentages, which is 35.93%. Both
figures are correct and answer different questions: an occurrence-weighted figure describes the
resistome, a genome-collapsed figure describes how often a compartment is used at all.

### Plasmid mobility typing

Resistance-gene-bearing plasmid replicons (n = 6,621) were typed with MOB-suite 3.1.9 using database
v3.1.8 [@robertson2018], which detects relaxase families under the MOB classification scheme
[@garcillanbarcia2009], mating-pair formation systems, and origins of transfer. Class C is a plasmid
with no detected mobility marker, class D carries a relaxase, and class E carries both relaxase and
mating-pair formation machinery. E1 and E2 partition class E and are not nested: E2 carries a detected origin of transfer in addition to relaxase and mating-pair formation machinery, E1 is the remainder. A negative result is a statement about the marker
database and is reported as such throughout.

### Chromosomal context construction and mobile-element annotation

Chromosomal occurrences were merged into shared context blocks. Each block spans the union of the
±10 kb windows of the occurrences it contains, so an occurrence with a neighbour within 20 kb shares
a block with it. Circular replicons were handled topologically: windows crossing the origin were
resolved by wrapping rather than clipping (57 wrapped blocks), and windows truncated by a replicon
end were recorded as truncated (5 blocks). This yielded **21,955 blocks** spanning up to 62,349 bp
(Supplementary Method 5).

Blocks were annotated by profile-homology search for insertion sequence and transposase markers
[@siguier2014; @siguier2006] and for integrase and integron markers [@gillings2014; @escudero2015],
giving **32,364 features** (29,331 and 3,033 respectively). An occurrence is class B if at least one
feature lies within **its own** ±10 kb window; block-level marker positivity is a different quantity
and is reported separately throughout, because the two differ (6,617 versus 6,616 blocks) and
conflating them was the single most likely misreading of this analysis.

### Distance analysis

For each chromosomal occurrence the distance to the nearest qualifying marker within its own window
was computed, with direct overlaps scored as zero. Occurrences with no qualifying marker in-window
were right-censored at 10,000 bp and never imputed. The primary estimand is weighted cumulative
detection *F*(*d*), the block-balanced share of occurrences whose nearest element lies within *d*.
Block balancing assigns weight 1/*m* to each occurrence in a block of *m*, so weights sum to 21,955
and a gene-dense neighbourhood does not count more than a sparse one. Restricted mean distance is
E[min(*D*, 10,000)]; medians are reported as "not reached" wherever *F*(10 kb) < 0.5 rather than
extrapolated (Supplementary Method 6).

A documented exception: 111 occurrences lie in a marker-positive block whose nearest in-block
feature falls outside that occurrence's own window (10,231–16,931 bp away). Using the out-of-window
feature would manufacture a distance the design declares unavailable, so these are censored at 10
kb. The arithmetic reconciles exactly: 16,414 block-positive − 111 = 16,303 window-positive = class
B.

### Structural insertion-sequence census

All 21,955 blocks were processed with ISEScan 1.7.3 [@xie2017] under a pinned environment — HMMER
3.3.2 [@eddy2011], BLAST+ 2.17.0 [@camacho2009], FragGeneScan 1.32 [@rho2010], Biopython 1.88
[@cock2009] — one thread per block, with no reruns and no per-block parameter variation. HMMER was
pinned to 3.3.2 because transposase open-reading-frame detection is sensitive to the profile search
version, and an earlier module had used 3.3.2; a default environment solve returned 3.4 and the
environment was rebuilt before the census rather than after. The census completed with no tool
failures and resolved 14,426 elements.

Three cohort counts differ and are not inconsistent. The cohort holds **6,288 genomes**. Of these,
**6,186 carry at least one chromosomal resistance gene** and contribute a genome-collapsed
chromosomal event; the remaining 102 carry resistance genes only on plasmids. Those 6,186
assemblies hold **6,190 ARG-bearing chromosomes**, because four of them carry two such
chromosomes rather than one.

An element is structurally complete if ISEScan reports it as type `c` with a complete transposase
open reading frame and bilateral resolved terminal inverted repeats. The primary structural endpoint
additionally requires **full containment within the individual occurrence's own ±10 kb window**.
Elements complete in the shared block but crossing that boundary form a separate, explicitly
reported category (2,384 elements, 716 occurrences) and are censored, because a window shorter than
an element truncates it by construction. This containment rule was registered as a numbered
amendment and hashed before any block was scored (Supplementary Method 7).

### Genome-wide background null

Insertion sequences were annotated across all 6,190 ARG-bearing chromosomes with ISEScan 1.7.3 under
the pinned environment described above, one chromosome at a time, yielding 190,999 elements. Element
completeness follows the definition used for the window-limited census without change: type `c` with
a complete transposase open reading frame and bilateral resolved terminal inverted repeats.

The null relocates each chromosomal occurrence uniformly at random within its own chromosome,
preserving chromosome, genome, species and BioProject identity, the number of occurrences on that
chromosome, and each occurrence's interval length. Relocated intervals may not overlap. On circular
chromosomes the ±10 kb window wraps; on linear chromosomes a start whose window would run past
either end is not permissible, matching how truncated windows are treated in the observed data.
Permutations never rerun ISEScan: they score against the frozen coordinate table. B = 2,000, seed
20260824, seeded per chromosome as seed XOR hash(accession) so that a chromosome replays identically
regardless of scheduling order.

Empirical *P* is (#{null ≥ observed} + 1)/(B + 1). With B = 2,000 the floor is 1/2001 ≈ 0.0005 and
is reported as ≤ 0.0005; *P* = 0 is never reported. Enrichment is the observed detection fraction
divided by its null expectation at the same landmark, with a 95% interval taken as the 2.5th and
97.5th percentiles of the null distribution.

The design, the permutation count, the seed, the four weighting schemes, the seven sensitivities and
the decision thresholds were registered in a frozen amendment before any species outcome was
computed, and the order in which the cohort was fixed and the outcome read is set out in
Extended Data Fig. 5. The analysis was verified through a separate code path that re-implements every mechanic
independently, including topology-aware distance on circular chromosomes, censoring, permutation
logic, all four weightings and the empirical *P* calculation; 20 of 20 checks passed. The census
itself was verified by 36 independent ISEScan re-runs, all of which reproduced the recorded element
count exactly.

**Accessory-context arms.** A chromosome-wide relocation null preserves a chromosome's element
density but not an occurrence's position within it, so it cannot separate the enrichment from the
fact that resistance genes and insertion sequences both concentrate in accessory regions. Two arms
were registered in `NM_C1_ACCESSORY_CONTEXT_AMENDMENT_009.json` before either was computed, reusing
the contrast floor of 1.5 already registered for the composite-element, non-ST2 and lineage-balanced
arms. In the **local relocation null**, each occurrence is relocated uniformly within ±*R* of its own
start rather than anywhere on its chromosome, with *R* = 50 kb and *R* = 100 kb registered together
so that neither could be chosen after the other was seen; everything else — interval length,
identity, wrapping, the complete-structural element definition, the 10 kb horizon, B = 2,000 — is
unchanged, and ISEScan is not rerun. One departure is deliberate: relocated intervals are not forced
to be mutually non-overlapping, because inside a local window around a resistance-gene cluster that
constraint is often unsatisfiable and it does not enter this estimator, each occurrence's distance
depending only on its own interval and the frozen coordinates. In the **resistance-island
stratification**, occurrences on a replicon are single-linkage clustered at a 10 kb gap and a cluster
of three or more is treated as an island; the enrichment is recomputed within each stratum from the
existing null matrices, with no new permutation. Neither arm is an accessory-restricted null: no
core-genome alignment exists in this project, and both are reported as approximations.
### Cargo convergence on documented replicons

The unit of this analysis is the replicon, not the occurrence. Every plasmid carrying at least one
acquired resistance gene was typed for mobility as above, giving 6,621 plasmids in three mutually
exclusive categories: conjugation-consistent (3,937), mobilization-consistent (1,211) and
marker-negative (1,473).

Cargo breadth is the number of distinct AMRFinderPlus drug classes represented among the acquired
resistance-gene occurrences lying on that replicon, and a replicon is multi-class if it carries at
least three. A determinant contributes its own class only, so a replicon carrying five genes of one
class is not multi-class. Metal-tolerance and virulence co-location use the STRESS, BIOCIDE, METAL
and VIRULENCE records that the primary denominator excludes; those layers are retained for this
analysis alone, and a replicon counts as co-locating when at least one such record lies on it.
Median resistance-gene count is the per-replicon count of primary-denominator occurrences.

Differences between mobility categories are reported with 95% percentile intervals from a cluster
bootstrap over BioProjects, B = 2,000, seed 20260822; replicon-clustered intervals were computed
alongside and are recorded in the result receipt. Four mobility-rule variants were run, differing in
whether mating-pair formation alone, or an origin of transfer, is sufficient evidence; the ordering
of the three categories is preserved under all four. These are associations on a shared replicon.
No transfer was observed, and co-location is not co-transfer.

### Statistics

Family-matched contrasts use Mantel–Haenszel pooling across matched gene families [@mantel1959],
with Cochran's *Q* and *I*² reported and the largest single-family weight share stated, so that a
pooled estimate dominated by one family is visible as such. Per-family compartment enrichment uses
Fisher's exact test with Benjamini–Hochberg control of the false discovery rate [@benjamini1995].

Uncertainty is estimated by a cluster bootstrap [@efron1979; @field2007; @davison1997] resampling
**BioProjects**, not occurrences or blocks, with B = 2,000. Blocks nest within exactly one
BioProject, so the resampling unit is well defined. Two interval methods are used and they are not interchangeable. Detection fractions, contrasts,
differences in proportions and the discordance residual carry **95% percentile** intervals from
the cluster bootstrap, taken as the 2.5th and 97.5th percentiles; no bias-corrected, accelerated
or studentized interval is used anywhere. The pooled Mantel–Haenszel odds ratios carry an
**analytic 95% Wald interval on ln(OR) with Robins–Breslow–Greenland variance**, not a bootstrap
interval. Each reported interval names its own method. Seeds are recorded for
reproducibility and carry no statistical justification. The chromosomal analysis set spans **2,248**
BioProjects whose effective number by the inverse Herfindahl–Hirschman index is **136.9**; the full
6,288-genome cohort spans 2,283 BioProjects at an effective 114.3. Both the nominal and the
effective count are reported wherever an interval is, because the nominal count overstates the
independence available. Multiplicity within each registered contrast family is controlled by Holm's
procedure [@holm1979]. The smallest attainable two-sided bootstrap *P* is 2/B = 0.001 — twice the one-sided floor,
because a two-sided value doubles the smaller tail — so a Holm-adjusted value of 0.002 across
the two primary contrasts is the resolution floor and is reported as *P* ≤ 0.002. An earlier
version of this manuscript quoted 1/B and a floor of 0.001; that overstated the attainable
resolution by a factor of two (Supplementary Method 8).

The discordance analysis fits the relationship between plasmid share and block-weighted chromosomal
association on the eight confirmation species only, excluding both species in which the pattern was
noticed, and evaluates the discovery species against that fit. This discovery/confirmation split was
registered before the fit was run.

The family-matched contrast is two-way. The numerator state is a chromosomal occurrence with a
mobile-element marker inside its own ±10 kb window (class B); the comparator is any plasmid-borne
occurrence (classes C, D and E). Chromosomal occurrences with no marker (class A) are excluded
rather than pooled with either state, because the endpoint asks how a determinant that is already
chromosomal is positioned, not whether it is chromosomal. The exclusion makes the contrast
asymmetric: every plasmid occurrence qualifies as a comparator, whereas a chromosomal occurrence
qualifies as a case only when a marker is detected. The three-state decomposition reported
alongside it does not impose that asymmetry, and the two are not interchangeable.

### Lineage typing and lineage adjustment

Sequence types were called with mlst 2.35.0 against the PubMLST schemes bundled with that
release, on every assembly contributing to the headline contrast (4,240 genomes; 780 in
*A. baumannii*, 3460 in the *Klebsiella* group). Reconciliation was exact: 4240 accessions in,
4240 typed, no failures. The scheme was **forced per host** — Pasteur `abaumannii_2` for
*A. baumannii* and `klebsiella` for the *Klebsiella* species complex — because *A. baumannii*
has two PubMLST schemes and per-genome auto-detection would mix Oxford with Pasteur types and
return sequence types that are not comparable with one another. A genome the scheme cannot
resolve is untypeable; untypeable genomes are not silently dropped, and both handlings are
reported.

Three adjustment arms and one influence arm were registered before any of them was computed
(`NM_V4C_LINEAGE_ADJUSTMENT_AMENDMENT_002.json`), together with the three outcomes the result
could take and what each would require of the manuscript. One outcome specified that the
claim be relabelled from host-conditioned to host-associated; that is the outcome that
occurred, and the relabelling follows the rule rather than a judgement made afterwards. The
estimator is unchanged from the primary analysis: only which occurrences enter, and with what
weight, differs. Sequence type is a seven-locus surrogate for lineage, not a phylogeny — two
genomes sharing a type are not clones, and the adjustment should be read as removing clonal
replication, not as controlling descent.

The intrinsic-determinant sensitivity was registered the same way
(`NM_V4C_INTRINSIC_SENSITIVITY_AMENDMENT_001.json`). Family membership for allele-level
exclusion was taken from the AMRFinderPlus reference protein database at version 2026-08-07.1,
the version that produced the calls, rather than assigned by hand.

### Verification programme

Five modules were executed under separate frozen protocols: blinded manual adjudication of the
chromosomal mobile-element layer (n = 120 stratified blocks, methods de-identified as X/Y/Z, the
adjudicator blind to species, study, gene identity and prior classification, and the adjudicator
being the author — the registered criterion of an uninvolved adjudicator was not met); lineage and
BioProject robustness including leave-one-BioProject-out; tool-version, database-version and
independent-implementation arms for the plasmid mobility layer, the last using CONJScan [@cury2020;
@cury2017] within MacSyFinder [@neron2023; @abby2014] and reported as concordance because neither
tool is a truth standard; a discovery/confirmation split for the discordance principle; and the
structural census described above. Each emitted an independent verification report stating its own
disagreement count. The sensitivity sets applied to the window-limited endpoint are defined in Supplementary
Method 9 and their results tabulated in Supplementary Table 9. The primary contrast at 1 kb stays positive
under all of them, within 0.05 of its +0.4186 baseline in every set but one: restricting the
marker to integrases and integrons alone leaves +0.0100, which is the point of that set — the
contrast is carried by insertion sequences and transposases, not by integrons.

### Reproducibility

All frozen protocols, receipts, result tables and analysis code carry full 64-character SHA-256
digests and are versioned rather than overwritten; superseded versions are retained. Figures are
regenerated from digest-verified inputs by a generator that verifies every input hash before drawing
and records every output hash. Three quantities quoted in this manuscript were recomputed
independently during figure generation and agreed exactly. Each headline quantity is traced to its source artefact in Supplementary Result 27, and tool and reference database versions are
listed in Supplementary Table 10.

### Use of AI tools

Large language models were used to assist with code development, quality-control documentation and
manuscript drafting under the author's direction. All analyses, outputs, citations and scientific
claims were checked by the author, who accepts full responsibility for the work.

No figure in this paper was drawn by a generative model. An early draft of the schematic in Extended Data Fig. 5 was, from a written specification, and it was replaced: the submitted figure is redrawn
by a script that reads every value from the result receipts, so it is reproducible rather than
trusted, and the model-drawn draft is retained separately as a record. No estimate, statistic or
scientific conclusion in this paper was produced by a generative model.

---

## Data availability

The occurrence-level dataset — 74,349 acquired resistance-gene occurrences with their replicon
assignment and portability class — and the genome-wide insertion-sequence annotation of 190,999 insertion sequences, 145,779 of them
structurally complete, across 6,190 chromosomes are deposited in Zenodo at
https://doi.org/10.5281/zenodo.22116987 (version 2.0.0; the concept DOI
https://doi.org/10.5281/zenodo.22065541 resolves to the latest version, and version 1.0.0 at
https://doi.org/10.5281/zenodo.22065542 remains available under its original terms). The deposit
includes the frozen chromosome accession manifest, the frozen background-null design, the primary
estimates, the seven registered sensitivities, the gate evaluation, the independent verification
report and the reconciliation against the window-limited endpoint, together with a data dictionary
defining every field, denominator and evidence-level dictionaries, JSON Schemas for every table, a
provenance map, pinned software and database versions, per-file SHA-256 checksums, and a standalone
verification script that re-derives every headline denominator from the archive alone. Raw
chromosome sequence is not redistributed: the accession manifest names every chromosome by
accession.version with its expected length and NCBI source, which regenerates the exact input.
Result and summary tables, frozen protocols, receipts and verification records are available at https://github.com/piranfar/portabilityrisk, which is also where the 4,240 sequence-type calls (`NM_V4C_MLST_CALLS_V1.tsv`) are released: they were produced after version 2 of the Zenodo deposit and will enter the next deposit version rather than being backdated into it. The deposited data are licensed CC BY 4.0. The primary
inputs are public NCBI records and remain governed by NCBI's own terms, which nothing here alters.

## Code availability

Analysis code is available at https://github.com/piranfar/portabilityrisk under the Apache
License 2.0, versioned and hashed, with the exact commands. The occurrence-level dataset the
scripts consume is the Zenodo deposit above; because that deposit is now published under
CC BY 4.0, the scripts that require it can be run and not merely read. The release carries the
frozen protocol amendments and result receipts for the genome-wide background null, the
lineage adjustment, the intrinsic-determinant sensitivity and the matching-unit sensitivity,
each with the script that produced it.

## Acknowledgements

Cloud-computing resources were supported through promotional credits provided by Oracle Cloud
Infrastructure. This work was performed without collaborators.

## Author contributions

V.P. conceived the study, designed and implemented the analysis, performed the internal audits, and wrote
the manuscript.

## Competing interests

The author declares no competing interests.

## Funding

This research received no specific grant from any funding agency in the public, commercial or
not-for-profit sectors.

---

## Figure legends

**Fig. 1 | The same gene family takes different routes in different hosts.** **a**, per-family log
odds ratio for chromosomal-and-mobile against plasmid-borne across the 58 matched families, point
area proportional to Mantel–Haenszel weight share, with 95% Wald intervals; the dotted line is the
pooled estimate. Intervals covering zero are drawn in grey — the axis is ln(odds ratio) — and 44 of 58 families exclude zero individually.
**b**, the same contrast under all three weighting arms with analytic 95% Wald intervals
(Robins–Breslow–Greenland variance), not bootstrap intervals; the fraction
beneath each point is the number of families favouring *A. baumannii*.

**Fig. 2 | Portability is a property of the occurrence, not the gene.** **a**, 164 gene families
occur in both compartments, distributed by the share of occurrences in their minority compartment.
**b**, compartment preference by family against family size; colour denotes significance after
Benjamini–Hochberg correction and survival of Mantel–Haenszel species adjustment. **c**, attrition
from 164 both-context families to the 95 that survive adjustment.

**Fig. 3 | The five-class portability architecture.** **a**, each occurrence takes exactly one of
five evidence-ranked classes; the five sum to 74,349. **b**, class E partitions into two evidence tiers, E2 carrying a detected origin of transfer and E1 not.

**Fig. 4 | Spatial and structural evidence for the chromosomal route.** **a**, weighted cumulative
detection against distance; solid lines the structural endpoint, dashed the homology endpoint, over
identical rows. **b**, both registered structural contrasts at four landmarks, with BioProject
cluster-bootstrap intervals; all exclude zero. **c**, share of class-B occurrences corroborated by a
structurally complete, fully contained element. **d**, median distance to the nearest homology
marker; not reached within 10 kb in either comparator.

**Fig. 5 | Resistance genes sit closer to a complete insertion sequence than the genome background
allows.** **a**, weighted cumulative detection *F*(*d*) at 1, 2, 5 and 10 kb, observed (bars) against
the permutation-null expectation (horizontal marks, with the null 95% interval), per group. **b**,
enrichment (observed / expected) at the two registered gate landmarks; the dashed line marks no
enrichment. **c**, the null distribution of *F*(1 kb) for *A. baumannii* over 2,000
within-chromosome permutations with the observed value marked; no permutation reaches it, so the
empirical *P* sits at the 1/2001 resolution floor and is reported as ≤ 0.0005. **d**, genome-wide
density of complete structural insertion sequences against enrichment; with four groups this cannot
establish independence, but *P. aeruginosa* carries the fewest elements per Mb and ranks second by
enrichment. The null relocates each occurrence within its own chromosome, so that chromosome's own
insertion-sequence density is preserved exactly in every permutation and enrichment is
density-normalised by construction.

**Fig. 6 | Plasmid fraction and chromosomal mobile-element association are non-redundant.** **a**,
the two axes, with the fit on eight confirmation species (grey) and the two discovery species
coloured; the arrow marks the *A. baumannii* residual. **b**, variance in chromosomal association
explained by plasmid share. **c**, both registered effects with bootstrap intervals.

**Extended Data Fig. 1 | Marker-type decomposition and insertion-sequence family retention.** **a**,
an insertion sequence–only endpoint reproduces the primary result while an integrase/integron-only
endpoint shows no *A. baumannii* excess. **b**, retention of elements through the structural gate by
insertion-sequence family.

**Extended Data Fig. 2 | Mobilization architecture tracks cargo convergence on documented
replicons.** **a**, 6,621 resistance-gene-bearing plasmids typed by mobility class:
conjugation-consistent 3,937, mobilization-consistent 1,211, marker-negative 1,473. **b**, share
carrying ≥3 drug classes and **c**, share co-locating a metal-tolerance determinant — in both panels
the denominator is the plasmid count of that mobility class, printed inside the bar, not 6,621.
**d**, BioProject-clustered difference between conjugation-consistent and marker-negative replicons,
95% percentile interval. Association, not co-transfer.

**Extended Data Fig. 3 | The chromosomal mobile compartment under four denominators.** **a**, the
estimate depends on the denominator: 46.39% occurrence-weighted at ±10 kb against 30.14%
block-weighted. **b**, the host ordering is preserved under both weightings. **c**, two block-level
quantities that differ by exactly one block.

**Extended Data Fig. 4 | Cohort and replicon resolution.** **a**, 184,538 determinant records reduce
to the primary denominator of 74,349 acquired resistance-gene occurrences. **b**, every occurrence
resolves to a documented replicon: 39,209 plasmid, 35,140 chromosome, none unresolved. **c**, plasmid
share on two denominators — occurrence-weighted and genome-collapsed.

**Extended Data Fig. 5 | Registration and verification.** The order in which the analysis was
governed, from frozen design to verdict. **Stage 1**, the protocol — null model, permutation count,
seed, weightings, the seven sensitivities and the four possible verdicts — was registered and
hashed before any outcome was computed; the analysis aborts if the hash does not match. **Stage 2**,
the genome-wide census was run under a pinned toolchain and checked by 36 independent re-runs on a
seeded selection written to disk before the first re-run, all of which reproduced the recorded
element count exactly. **Stage 3**, the cohort was built and its receipt written to disk before any
estimate existed; the dashed line marks the point beyond which the outcome could not be read, and
the amber card marks where reading it became permitted. **Stage 4**, the four gates were applied
mechanically against thresholds unchanged from registration, and every mechanic was re-implemented
through a separate code path that does not import the primary analysis. The figure reports governance, not results: no value in it is an estimate, and the verdict
shown is one of four registered in advance. That verdict string is the frozen name of the
registered outcome and is reproduced unchanged; the claim it belongs to is reported as
host-*associated* in this paper, for the reasons given in the Results.
