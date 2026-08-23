# NM-V4C — interpretation qualifications, and the D-NM4 clonality plan

Companion to `NMV4C_FROZEN_DESIGN.json` (`8a3c76b157cbf2cd…`) and `NMV4C_RESULT_RECEIPT.json`.
Written after scoring. Alters no frozen artefact and no receipt.

---

## 1. The qualification the gate count does not show

All seven frozen gates passed. One result nevertheless constrains what may be claimed, and it
must be stated beside the verdict rather than below it.

The frozen design declares **P1-full** and **P2-one-genome-per-BioProject** to be
**co-primary**, and rules that if they disagree in direction or significance, "the result is
reported as sampling-structure dependent and the universal claim is not made."

| contrast | P1 OR | P2 OR | co-primary |
|---|---:|---:|---|
| *A. baumannii* vs *Klebsiella* | 50.29 [45.61, 55.45] | 21.95 [18.92, 25.47] | **agree** |
| *P. aeruginosa* vs *Klebsiella* | 28.94 [25.79, 32.46] | 16.06 [13.79, 18.70] | **agree** |
| ***A. baumannii* vs *P. aeruginosa*** | 1.26 [1.11, 1.43] | **1.18 [0.99, 1.41]** | **DISAGREE** |

Gate G5 was written to test the *A. baumannii* vs *P. aeruginosa* contrast without naming
which co-primary it applied to; my implementation used P1, so G5 registered as passed. Under
the co-primary rule the honest reading is different: **the two low-plasmid hosts are not
reliably separated by the B-versus-plasmid route contrast** once BioProjects are balanced.

### Why that is not a failure, and what it actually shows

The five-class representation — which the design required precisely so class A could not be
silently discarded — explains it:

| host | A | B | C | D | E |
|---|---:|---:|---:|---:|---:|
| *A. baumannii* | 16.4 % | **69.9 %** | 3.7 % | 3.6 % | 6.4 % |
| *P. aeruginosa* | **50.5 %** | 36.2 % | 11.3 % | 1.1 % | 1.0 % |
| *Klebsiella* | 19.9 % | 11.3 % | 9.6 % | 9.4 % | **49.8 %** |

*A. baumannii* and *P. aeruginosa* are both plasmid-poor, so **B / (C+D+E)** is similar for
both and the focused contrast has little to detect. They differ instead on the axis the
focused contrast conditions away: *P. aeruginosa* puts **half** its acquired ARGs in class A —
chromosomal with no mobile-element marker within 10 kb — against 16.4 % for *A. baumannii*.

So the three hosts occupy three distinct architectures, not two:

- *Klebsiella* — plasmid route, conjugation-consistent (49.8 % class E)
- *A. baumannii* — chromosomal **mobile** route (69.9 % class B)
- *P. aeruginosa* — chromosomal **quiescent** route (50.5 % class A)

This is a **stronger** result than the one the focused contrast was designed to find, and it
is exactly the Mulkern-style structure: the determinant does not specify the route; the
determinant–host combination does. Low plasmid fraction is not one phenotype.

### What may and may not be claimed

**Permitted:** *"Portability architecture is a property of the determinant–host combination:
identical ARG families occupy systematically different chromosomal-MGE and plasmid mobility
contexts across bacterial hosts."*

**Permitted, and worth stating:** the two low-plasmid hosts diverge on the A/B axis rather
than on the plasmid axis.

**Not permitted without an amendment:** a formal A-versus-B inferential test. The frozen design
fixed the 2×2 as B vs C+D+E and reports class A descriptively. Testing A vs B is a new
contrast; running it now and reporting it as confirmatory would be choosing a test after
seeing the data. It requires a numbered amendment written before it is fitted.

**Not permitted at all:** observed horizontal transfer, causal host control, demonstrated
conjugation, transfer rate, clinical risk.

---

## 2. D-NM4 — plan for controlling chromosomal clonality

**Not executed.** No sequence was retrieved, no clustering was run, no threshold was applied.

### Why BioProject adjustment is not clonal correction

They are different confounders and neither substitutes for the other.

- A **BioProject** is a submission unit. Adjusting for it removes correlation induced by one
  laboratory sequencing many isolates together.
- A **clone** is a descent unit. Two isolates from the same global clone submitted by two
  laboratories on two continents are BioProject-independent and clonally identical.

NM-V2 showed the *A. baumannii* result survives BioProject adjustment decisively — 269
BioProjects, largest contributing 6.3 %, no single removal moving the ratio more than 5.50 %.
That is genuine evidence against a *submission* artefact. It is **not** evidence against a
*clonal* artefact. Global *A. baumannii* is dominated by a small number of international clonal
lineages (GC1/GC2), and the AbaR-type islands that produce class B are a lineage property. If
most cohort genomes belong to two lineages, the finding could be a statement about two clones
rather than about the species.

**This is currently the single largest unquantified risk to claim C07**, and it must not be
described as resolved by NM-V2.

### Candidate methods

| method | resolution | inputs | cost | independent threshold available |
|---|---|---|---|---|
| **7-locus MLST** | sequence type | assemblies | low | yes — PubMLST scheme definitions |
| **cgMLST** | clonal complex | assemblies | medium | yes — published scheme thresholds |
| **Mash / skani ANI** | fine clusters | assemblies | low | **no** — must be derived |
| **core-genome SNP** | transmission-level | assemblies + alignment | high | partly |

**Recommendation: cgMLST as primary, skani ANI as a cross-check.** cgMLST carries published,
externally defined clonal-complex thresholds, which is the only way to satisfy the protocol's
rule that a threshold be justified independently of any outcome. ANI clustering alone would
force us to invent a cutoff, which is exactly what the protocol prohibits.

The project already has precedent and tooling: **FastANI 1.34 and skani 0.2.2** were installed
and used for PR-ANI-1, which adjudicated 419 genomes. The method is not new to this programme.

### Required sequence retrieval and storage

| item | quantity |
|---|---|
| genomes needing retrieval | **6,288** (all ARG-bearing cohort genomes) |
| minimum viable subset | 780 *A. baumannii* + 2,822 *Klebsiella* + 884 *P. aeruginosa* = **4,486** |
| mean assembly size | ≈ 5.5 Mb |
| uncompressed | ≈ 35 GB full cohort, ≈ 25 GB minimum subset |
| gzipped | ≈ 10 GB full cohort |
| cgMLST scheme databases | ≈ 1–2 GB |
| working space | ≈ 60 GB recommended |

Retrieval is from NCBI by exact `assembly_version`, the same route already used for plasmids
and windows, with per-file digests recorded. **This is a genuine new download and requires
explicit authorisation** — it is far larger than anything NM-0 permitted.

### Expected server shape and duration

| stage | shape | duration |
|---|---|---|
| assembly retrieval, 4,486–6,288 genomes | network-bound, 8–16 parallel streams | 3–6 h |
| cgMLST allele calling | CPU-bound, 16–32 cores | 8–20 h |
| skani ANI all-vs-all cross-check | CPU-bound, embarrassingly parallel | 1–3 h |
| clustering and leave-one-lineage-out re-analysis | minutes on a laptop | < 1 h |

**A server is required.** This does not fit comfortably on the laptop, and the former analysis
server no longer exists — provisioning is **D-NM5**.

### Threshold justification, fixed in advance

1. **Primary:** adopt the **published** clonal-complex threshold of the chosen cgMLST scheme,
   cited to its source. Not derived from our data at all.
2. **Cross-check:** if an ANI threshold is needed, derive it from the **within-BioProject**
   distance distribution — pairs known to be sampling-related — examined **without reference to
   any class, MGE or plasmid quantity**, and freeze it before recomputing any outcome.
3. **Prohibited:** selecting the threshold that maximises or minimises the *A. baumannii*
   effect; trying several and reporting the best; adjusting after seeing the result.

### What the analysis then delivers

Leave-one-major-lineage-out on the *A. baumannii* chromosomal-MGE fraction, plus lineage-aware
resampling as an additional clustering level above BioProject. **Prespecified reading:** if the
*A. baumannii* vs *Klebsiella* ratio survives removal of the largest lineage with its CI still
excluding 1, the finding is a species-level architecture. If it collapses, it is a lineage
property and must be reported as such — which would still be publishable, but as a different
and narrower claim.

---

## 3. Remaining blockers

| blocker | status | blocks |
|---|---|---|
| **NM-V1** stratified MGE validation | **NOT RUN** — needs ISEScan, IntegronFinder, sequences, 4–9.5 h | claims C05, C06, and every class-B number including the 69.9 % above |
| **NM-V3** tool/database robustness | **NOT RUN** — needs a second MOB-suite database | claim C10, and the C/D/E boundaries underpinning the plasmid arm |
| **D-NM4** clonal structure | **NOT RUN** — plan above; needs D-NM5 and a large download | claim C07 at species level |

**The manuscript is not Nature-ready.** All three blockers are open, and NM-V4C depends on the
same unvalidated MGE layer as C05: every class-B count in this module inherits NM-V1's risk.
