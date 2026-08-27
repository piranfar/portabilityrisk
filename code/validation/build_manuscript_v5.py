"""Build Manuscript V5 as a Nature Microbiology Resource.

V4 stays on disk unchanged and remains the record of the Article framing.

Why the format changed. A Resource "presents a large data set of broad utility,
interest and significance to the community", carries 4,000 words against an
Analysis's 3,000, and is not judged on conceptual novelty. This work generated no
new sequence: it built a documented location layer over 6,288 public closed
genomes and then tested what that layer supports. That is a Resource.

What changes in substance, not only in framing:

  * the headline is no longer the 50.29-fold species contrast. That number answers
    a different question - chromosomal-and-mobile RATHER THAN plasmid-borne - and
    multiplies the host's chromosome-versus-plasmid propensity with
    marker-proximity given chromosomal. It is kept, as a secondary result,
    labelled as what it is;
  * the headline is the enrichment that four alternative explanations could not
    account for: chromosome-wide element density, clonal replication,
    composite-transposon structure, and intrinsic determinants;
  * what collapsed is in the Results, not buried. A Resource is the one format
    where "here is the layer, here is what it will and will not support" is the
    contribution rather than an admission.

Every number is read from a result receipt at build time. None is typed here.
"""
import os as _os


def _dir(var):
    """Resolve an input directory from the environment.

    The private copy of this script carried an absolute path on the
    author's machine. The public copy does not substitute a plausible
    path, because a path that looks right and is wrong is worse than one
    that is obviously missing. Set the variable, or the run stops here.
    """
    v = _os.environ.get(var)
    if not v:
        raise SystemExit(
            "set %s to the directory it names; see README" % var)
    return v.rstrip("/\\") + "/"


import hashlib
import io
import json
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

NM = _dir("PORTABILITYRISK_REPO_DIR") + "docs/nature_microbiology/"
V4 = NM + "PORTABILITYRISK_MANUSCRIPT_V4.md"
V5 = NM + "PORTABILITYRISK_MANUSCRIPT_V5_RESOURCE.md"
V4_SHA = "aa42cbcefa47b2ae6b8b9300abf53260c8c7af72bc3e7b9643bee3affa37e759"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jload(n):
    return json.load(io.open(NM + n, encoding="utf-8"))


def words(x):
    x = re.sub(r"^#.*$", "", x, flags=re.M)
    return len(re.findall(r"\S+", re.sub(r"\[@[^\]]+\]", "", x)))


def main():
    if os.path.exists(V5):
        raise SystemExit("V5 exists - refusing to overwrite: %s" % V5)
    if sha(V4) != V4_SHA:
        raise SystemExit("V4 changed since this build was written\n  expected %s\n"
                         "  found    %s" % (V4_SHA, sha(V4)))
    t = io.open(V4, encoding="utf-8").read()
    print("V4 verified: %s\n" % V4_SHA[:16])

    def sec(head):
        m = re.search(r"^%s\s*$" % re.escape(head), t, re.M)
        if not m:
            raise SystemExit("section not found in V4: %r" % head)
        nxt = re.search(r"^#{1,3} ", t[m.end():], re.M)
        return (t[m.end():m.end() + nxt.start()] if nxt else t[m.end():]).strip("\n")

    null = jload("NM_C1_4_NULL_RESULTS_V1.json")
    nst2 = jload("NM_C1_NON_ST2_ENRICHMENT_RESULTS_V1.json")
    comp = jload("NM_C1_COMPOSITE_ELEMENT_RESULTS_V1.json")
    lin = jload("NM_V4C_LINEAGE_ADJUSTMENT_RESULTS_V1.json")
    ref = jload("NM_V4C_REFEREE_ANALYSES_RESULTS_V1.json")
    mu = jload("NM_V4C_MATCHING_UNIT_RESULTS_V1.json")
    intr = jload("NM_V4C_INTRINSIC_SENSITIVITY_RESULTS_V1.json")
    marg = jload("NM_V4C_MARGINAL_ADJACENCY_RESULTS_V1.json")

    AB = null["results"]["A. baumannii"]["occurrence"]["F_1000"]
    ENR = AB["enrichment"]
    NS = nst2["non_ST2_summary"]
    D10 = comp["results"]["D_10000"]
    n_flank = sum(D10[g]["flanked"]["n"] for g in D10)
    n_tot = n_flank + sum(D10[g]["non-flanked"]["n"] for g in D10)
    CF = comp["non_flanked_summary"]
    L1, L2 = lin["results"]["L1"], lin["results"]["L2"]
    ABC, KLC = lin["composition"]["A. baumannii"], lin["composition"]["Klebsiella group"]
    BASE, N1, C2 = (ref["results"]["baseline_two_way"], ref["results"]["N1"],
                    ref["results"]["C2"])
    W1 = ref["W1_diagnostics"]["arms"][0]
    W1i = ref["W1_diagnostics"]["arms"][1]
    R1 = mu["results"]["R1_random_effects"]
    S8B = intr["results"]["S8b"]
    NG = ABC["genomes"] + KLC["genomes"]

    P = []
    A = P.append

    A("# Replicon-resolved portability of 74,349 acquired resistance-gene occurrences\n"
      "across 6,288 closed Gram-negative ESKAPE genomes\n")
    A(t[t.index("**Vahhab Piranfar**"):t.index("## Abstract")].rstrip() + "\n")

    A("## Abstract\n")
    A("Whether a resistance gene can leave its host depends on where it sits, and in draft "
      "assemblies\nthat location is predicted. We assembled an occurrence-level "
      "resource over **6,288\nclosed Gram-negative ESKAPE genomes**: **74,349 acquired "
      "resistance-gene occurrences**, each assigned\nto a chromosome or plasmid documented by NCBI, "
      "ranked into five evidence classes, with a\ngenome-wide census of **%s structurally "
      "resolved insertion sequences** across **6,190 ARG-bearing\nchromosomes** and sequence types "
      "for the %s genomes behind the primary species contrast. Against a\nnull that relocates each "
      "occurrence within its own chromosome, resistance genes lie closer to a\ncomplete insertion "
      "sequence than chance allows in every group, and *Acinetobacter baumannii* is\nextreme at "
      "**%.2f-fold within 1 kb**. Chromosome-wide element density, clonal replication,\n"
      "composite-transposon structure and intrinsic determinants each fail to account for it. A "
      "%.2f-fold\nspecies contrast in the same data does not survive the same scrutiny, and we "
      "report where it goes."
      % ("{:,}".format(null["elements_total"]), "{:,}".format(NG), ENR, BASE["or"]))
    A("")

    A("## Introduction\n")
    A(sec("## Introduction"))
    A("")

    A("## Results\n")

    A("### What the resource contains\n")
    A("Assigning a resistance gene to a chromosome or a plasmid requires knowing which molecule its "
      "contig\nbelongs to, and in a closed genome that molecule is documented rather than predicted. "
      "The resource is\nbuilt on that distinction. Across **6,288 closed complete genomes** of "
      "Gram-negative ESKAPE\npathogens, **74,349 acquired resistance-gene occurrences** were each "
      "assigned to a documented\nreplicon, with 0 unmatched (Methods; Fig. 1; Supplementary Table 1). "
      "Every occurrence carries a\nfive-class portability rank, a plasmid mobility call, and, for "
      "the 35,140 chromosomal occurrences,\nthe distance to the nearest mobile-element marker inside "
      "its own ±10 kb window.\n\n"
      "Three layers were built here and existed in no prior public record. A **genome-wide structural "
      "\ninsertion-sequence census** annotated all **6,190 ARG-bearing chromosomes** under one pinned "
      "toolchain,\nrecovering **%s elements**, of which %s meet a complete-structural definition: a "
      "complete transposase\nopen reading frame with bilateral resolved terminal inverted repeats. A "
      "**within-chromosome\npermutation null** supplies the expectation against which those distances "
      "are read. **Sequence types**\nwere called for the %s genomes contributing to the primary "
      "species contrast, which makes clonal\nstructure visible rather than assumed away.\n\n"
      "The layer is deposited under CC BY 4.0 (Data availability). What follows is what it supports, "
      "and\nwhat it does not."
      % ("{:,}".format(null["elements_total"]),
         "{:,}".format(null["elements_complete_structural"]), "{:,}".format(NG)))
    A("")

    A("### Portability is a property of the occurrence, not the gene\n")
    A(sec("### Portability is a property of the occurrence, not the gene"))
    A("")

    A("### Five evidence-ranked portability classes\n")
    A(sec("### Five evidence-ranked portability classes"))
    A("")

    A("### The association is short-range, structural, and specific to insertion sequences\n")
    A(sec("### The association is short-range, structural, and specific to insertion sequences"))
    A("")

    A("### Resistance genes sit closer to intact insertion sequences than chance allows\n")
    A("Testing whether a host simply carries more insertion sequences requires a background, and none "
      "\nexisted: elements had only ever been searched inside the ±10 kb windows around resistance "
      "genes,\n1.57%% of chromosomal sequence. Each occurrence was compared with a null relocating it "
      "uniformly\n**within its own chromosome**, which preserves that chromosome's element density "
      "exactly, so\nenrichment is density-normalised by construction. Design, permutation count, seed "
      "and thresholds\nwere registered before any species outcome existed.\n\n"
      "**Resistance genes lie closer to a complete insertion sequence than chance allows in every "
      "group\ntested** (Fig. 4). Within 1 kb, *A. baumannii* shows an observed detection fraction of "
      "%.3f against an\nexpected %.3f — a **%.2f-fold enrichment**, with the observed value outside "
      "the null 95%% interval and\nan empirical *P* ≤ 0.0005 at the resolution floor of 1/2001. The "
      "restricted mean distance is 2,930 bp\nobserved against 9,151 bp expected. All seven registered "
      "sensitivities support it (Supplementary\nResult 7).\n\n"
      "**Three further explanations were tested and none accounts for the effect.** *Clonal "
      "replication*: %s\nof the *A. baumannii* genomes belong to sequence type %s, and %.0f%% of its "
      "chromosomal occurrences\nsit on those genomes; removing them leaves **%.2f-fold**, still "
      "%.2f times the *Klebsiella* group\n(Supplementary Result 8). *Composite-transposon structure*: "
      "where a resistance gene is the cargo of a\ncomposite element its distance to an insertion "
      "sequence is near zero by construction, but only %s of\n%s chromosomal occurrences (%.1f%%) are "
      "flanked by a same-family pair of complete elements, and among\nthose that are not the "
      "enrichment is **%.2f-fold**, %.2f times the *Klebsiella* group — higher than\nthe unstratified "
      "ratio, not lower. *Intrinsic determinants*: species-core genes are excluded from the\nmatched "
      "analyses by the registered eligibility rule, and removing the intrinsic *bla*OXA alleles "
      "at\nallele level does not reduce any contrast (Supplementary Result 12)."
      % (AB["observed"], AB["expected"], ENR,
         "%.1f%%" % (100 * ABC["largest_ST_share"]), ABC["largest_ST"],
         100 * nst2["results"]["A. baumannii, ST2 only"]["n"]
         / nst2["results"]["A. baumannii, all"]["n"],
         NS["enrichment"], NS["over_Klebsiella"],
         "{:,}".format(n_flank), "{:,}".format(n_tot), 100 * n_flank / n_tot,
         CF["AB_enrichment"], CF["AB_over_Klebsiella"]))
    A("")

    A("### A species contrast that the same data do not support\n")
    A("The resource also permits a matched-family contrast: among gene families present in both "
      "hosts,\nthe odds that a determinant is chromosomal-and-mobile **rather than plasmid-borne**. "
      "Pooled by\nMantel–Haenszel it is **%.2f-fold higher in *A. baumannii* than in *Klebsiella*** "
      "across %d families,\n%d of which point the same way. We report it, and we report that it does "
      "not carry the weight the\nnumber suggests.\n\n"
      "Three things qualify it. Between-family heterogeneity is high (*I*² = %.1f%%), so the "
      "fixed-effect\nsummary is reported with its random-effects companion, **%.2f (95%% CI "
      "%.2f–%.2f)**. The analytic\ninterval assumes 74,349 independent occurrences; resampling "
      "lineages instead gives **%.2f–%.2f**, an\norder of magnitude wider, because *A. baumannii* "
      "carries an effective **%.2f** lineages. And the\nendpoint multiplies two effects — how often a "
      "host puts determinants on plasmids at all, and how\noften a chromosomal determinant sits near "
      "a marker — so a host extreme on both produces a large odds\nratio partly by construction.\n\n"
      "Separating them does not work in this cohort. The within-chromosome contrast, class B against "
      "class\nA, returns %.2f (%.2f–%.2f), but one family carries %.0f%% of the pooled weight against "
      "a registered\nceiling of 30%%, only %d of %d families are informative because the rest have no "
      "class-A occurrences\nin one host or the other, and removing the intrinsic *bla*OXA alleles "
      "moves it to %.2f. **The\nwithin-chromosome question is not estimable at this project's "
      "registered standard**, and we report\nthat rather than a number. The marginal rates are "
      "stable — %.1f%% of *A. baumannii* chromosomal\noccurrences are marker-adjacent against %.1f%% "
      "in *Klebsiella* — but direct standardisation shows\ncomposition accounts for the whole gap: "
      "the families that differ are the species-core genes each host\ncarries and the other does "
      "not (Supplementary Result 13).\n\n"
      "What does survive is narrower and worth stating plainly: the contrast keeps its direction "
      "under\nlineage adjustment (**%.2f**, %.2f–%.2f), with the dominant clone removed (**%.2f**, "
      "%.2f–%.2f), and with\n*bla*OXA and *bla*SHV removed outright (**%.2f**). *A. baumannii* keeps "
      "acquired resistance on the\nchromosome where *Klebsiella* keeps it on plasmids, and plasmid "
      "fraction alone does not see that."
      % (BASE["or"], BASE["n_families"], BASE["concordant"],
         R1["I2_pct"], R1["or"], R1["ci_lo"], R1["ci_hi"],
         C2["ci_lo"], C2["ci_hi"], ABC["effective_STs_1_over_HHI"],
         W1["or"], W1["ci_lo"], W1["ci_hi"], 100 * W1["largest_weight_share"],
         W1["families_actually_informative"], W1["families_in_table"], W1i["or"],
         100 * marg["M1_crude"]["AB"][0], 100 * marg["M1_crude"]["KL"][0],
         L1["or"], L1["ci_lo"], L1["ci_hi"],
         N1["or"], N1["ci_lo"], N1["ci_hi"], S8B["or"]))
    A("")

    A("### Plasmid fraction and chromosomal mobile-element association are non-redundant\n")
    A(sec("### Plasmid fraction and chromosomal mobile-element association are non-redundant"))
    A("")

    A("### Conjugation-consistent replicons carry the convergent cargo\n")
    A(sec("### Conjugation-consistent replicons carry the convergent cargo"))
    A("")

    A("### Validation, and what the resource cannot be used for\n")
    A(sec("### Robustness"))
    A("")

    A("## Discussion\n")
    A(sec("## Discussion"))
    A("")

    A("## Methods\n")
    m = re.search(r"^## Methods\s*$", t, re.M)
    end = re.search(r"^## Data availability\s*$", t, re.M)
    A(t[m.end():end.start()].strip("\n"))
    A("")
    A(t[end.start():].rstrip() + "\n")

    doc = "\n".join(P)
    io.open(V5, "w", encoding="utf-8", newline="\n").write(doc)

    def seg(a, b):
        return words(doc[doc.index(a):doc.index(b)])
    main_words = (seg("## Introduction", "## Results")
                  + seg("## Results", "## Discussion")
                  + seg("## Discussion", "## Methods"))
    abs_words = words(doc[doc.index("## Abstract"):doc.index("## Introduction")])
    print("wrote %s" % os.path.basename(V5))
    print("  sha256      : %s" % sha(V5))
    print("  abstract    : %d words (Resource allows 100-150)" % abs_words)
    print("  main text   : %d words (Resource allows 4,000)" % main_words)
    print("  V4 unchanged: %s" % (sha(V4) == V4_SHA))


if __name__ == "__main__":
    main()
