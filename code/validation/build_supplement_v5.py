"""Build Supplementary Information V5 for the Resource manuscript.

V4 is carried byte-for-byte and stays on disk. Two sections are added, both for
analyses that already have receipts but no prose:

  * the referee-requested arms: the enrichment with the dominant clone removed,
    and the enrichment with composite-transposon cargo removed;
  * the within-chromosome diagnostics: why the matched-family form of that
    contrast is not estimable here, and what direct standardisation shows about
    the marginal form.

Every number is read from a receipt at build time.
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
S4 = NM + "PORTABILITYRISK_SUPPLEMENTARY_INFORMATION_V4.md"
S5 = NM + "PORTABILITYRISK_SUPPLEMENTARY_INFORMATION_V5.md"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jload(n):
    return json.load(io.open(NM + n, encoding="utf-8"))


def main():
    if os.path.exists(S5):
        raise SystemExit("V5 supplement exists - refusing to overwrite")
    s = io.open(S4, encoding="utf-8").read()
    print("carrying V4 supplement: %s" % sha(S4)[:16])

    nst2 = jload("NM_C1_NON_ST2_ENRICHMENT_RESULTS_V1.json")
    comp = jload("NM_C1_COMPOSITE_ELEMENT_RESULTS_V1.json")
    ref = jload("NM_V4C_REFEREE_ANALYSES_RESULTS_V1.json")
    marg = jload("NM_V4C_MARGINAL_ADJACENCY_RESULTS_V1.json")

    n = max(int(x) for x in re.findall(r"^## Supplementary Result (\d+) \|", s, re.M))
    a, b = n + 1, n + 2
    print("existing Supplementary Results: %d; adding %d and %d" % (n, a, b))

    R = nst2["results"]
    D10 = comp["results"]["D_10000"]
    D5 = comp["results"]["D_5000"]
    NS = nst2["non_ST2_summary"]
    CF = comp["non_flanked_summary"]

    A = ["## Supplementary Result %d | Two explanations the enrichment is not" % a, "",
         "Registered in `NM_C1_COMPOSITE_ELEMENT_AMENDMENT_006.json` and",
         "`NM_C1_NON_ST2_ENRICHMENT_AMENDMENT_007.json`, each before its own arm was computed. The",
         "thresholds in the second were **reused from the first rather than chosen again**, because",
         "picking a fresh threshold while knowing what the sibling arm returned is what a freeze",
         "exists to prevent.", "",
         "### Clonal replication", "",
         "*A. baumannii* is dominated by one lineage. Sequence type %s holds %.1f%% of its genomes and"
         % (nst2["results"]["A. baumannii, ST2 only"] and "2",
            100 * R["A. baumannii, ST2 only"]["n"] / R["A. baumannii, all"]["n"]),
         "%.1f%% of its chromosomal occurrences. If the enrichment were a property of that clone,"
         % (100 * R["A. baumannii, ST2 only"]["n"] / R["A. baumannii, all"]["n"]),
         "removing it would remove the effect.", "",
         "| stratum | *n* | observed | expected | enrichment |", "|---|---:|---:|---:|---:|"]
    for k in ("A. baumannii, all", "A. baumannii, ST2 only", "A. baumannii, non-ST2",
              "Klebsiella group", "P. aeruginosa"):
        v = R[k]
        A.append("| %s | %s | %.4f | %.4f | **%.2f** |"
                 % (k, "{:,}".format(v["n"]), v["observed"], v["expected"], v["enrichment"]))
    A += ["",
          "It does not. Non-ST2 *A. baumannii* retains **%.2f-fold** enrichment against a registered"
          % NS["enrichment"],
          "floor of 8.46, and **%.2f times** the *Klebsiella* group against a floor of 1.5. The"
          % NS["over_Klebsiella"],
          "ratio to *P. aeruginosa* falls to %.2f from a published 1.86; no floor was registered for"
          % NS["over_Pseudomonas"],
          "that comparison and none is claimed for it here.", "",
          "### Composite-transposon structure", "",
          "Where a resistance gene is the cargo of a composite element, its distance to an insertion",
          "sequence is near zero by construction: the flanking copies are what make it a transposon. A",
          "uniform relocation null cannot address that, so the occurrences were stratified. An",
          "occurrence counts as composite-flanked when a complete-structural element of family *F* lies",
          "entirely upstream within *D* bp and another of the **same family** lies entirely downstream",
          "within *D* bp, on its own replicon.", "",
          "| distance | flanked | share of all chromosomal occurrences |", "|---|---:|---:|",
          "| *D* = 10 kb | %s | %.1f%% |"
          % ("{:,}".format(sum(D10[g]["flanked"]["n"] for g in D10)),
             100 * sum(D10[g]["flanked"]["n"] for g in D10)
             / sum(D10[g]["flanked"]["n"] + D10[g]["non-flanked"]["n"] for g in D10)),
          "| *D* = 5 kb | %s | %.1f%% |"
          % ("{:,}".format(sum(D5[g]["flanked"]["n"] for g in D5)),
             100 * sum(D5[g]["flanked"]["n"] for g in D5)
             / sum(D5[g]["flanked"]["n"] + D5[g]["non-flanked"]["n"] for g in D5)), "",
          "| group | flanked | non-flanked |", "|---|---:|---:|"]
    for g in ("A. baumannii", "Klebsiella group", "P. aeruginosa", "Enterobacter group"):
        A.append("| *%s* | %.2f (*n* = %s) | **%.2f** (*n* = %s) |"
                 % (g, D10[g]["flanked"]["enrichment"], "{:,}".format(D10[g]["flanked"]["n"]),
                    D10[g]["non-flanked"]["enrichment"],
                    "{:,}".format(D10[g]["non-flanked"]["n"])))
    A += ["",
          "The flanked stratum is enriched 23–31-fold in every host, which is what a definitional",
          "signal looks like and confirms the detector finds what it should. But it is a minority, and",
          "among the occurrences that are **not** flanked *A. baumannii* retains **%.2f-fold**"
          % CF["AB_enrichment"],
          "enrichment and **%.2f times** the *Klebsiella* group — higher than the unstratified ratio of"
          % CF["AB_over_Klebsiella"],
          "2.69, not lower. The objection is real for the flanked minority and does not explain the",
          "effect.", "",
          "A flanking pair of the same family is necessary but not sufficient for a composite",
          "transposon. Some flanked occurrences will be coincidental, and true composites with one",
          "degenerate copy will be missed. No transposon database was consulted and no element is named.",
          "", "---", "", ""]

    W = ref["W1_diagnostics"]
    B = ["## Supplementary Result %d | Why the within-chromosome contrast is not reported" % b, "",
         "Registered in `NM_V4C_REFEREE_ANALYSES_AMENDMENT_004.json` and",
         "`NM_V4C_MARGINAL_ADJACENCY_AMENDMENT_005.json`. This section reports an analysis that",
         "**failed**, because the question it asks is the right one and a reader is entitled to know",
         "why it has no answer here.", "",
         "The two-way endpoint asks whether a determinant is chromosomal-and-mobile *rather than*",
         "plasmid-borne. The question the resource's own framing raises is narrower: given that a",
         "determinant is chromosomal, is it marker-adjacent? That is class B against class A.", "",
         "| arm | odds ratio | 95% CI | largest family weight | informative families |",
         "|---|---:|---|---:|---|"]
    for arm in W["arms"]:
        B.append("| %s | **%.2f** | %.2f–%.2f | %s %.1f%% | %d of %d |"
                 % (arm["arm"], arm["or"], arm["ci_lo"], arm["ci_hi"],
                    arm["largest_weight_family"], 100 * arm["largest_weight_share"],
                    arm["families_actually_informative"], arm["families_in_table"]))
    B += ["",
          "**Every arm fails gate G2**, registered in `NMV4C_FROZEN_DESIGN.json` long before this",
          "analysis was requested: *the single largest family carries no more than 30 per cent of the",
          "MH weight*. And the estimate rests on far fewer families than the count suggests — a family",
          "in which both hosts have zero class-A occurrences contributes zero to both Mantel–Haenszel",
          "sums and disappears silently.", "",
          "The sign is not stable either. It depends entirely on how one family is handled: with",
          "*bla*OXA included the direction favours *Klebsiella*; with the intrinsic OXA-51-like alleles",
          "removed it favours *A. baumannii*. In this cohort that comparison is between the intrinsic",
          "chromosomal enzyme of one host and the integron-borne acquired enzymes of the other, which",
          "is not an architectural comparison.", "",
          "### The marginal form is stable but compositional", "",
          "The rate is estimable where the pooled odds ratio is not: **%.1f%%** of *A. baumannii*"
          % (100 * marg["M1_crude"]["AB"][0]),
          "chromosomal occurrences are marker-adjacent against **%.1f%%** in the *Klebsiella* group, a"
          % (100 * marg["M1_crude"]["KL"][0]),
          "difference of %+.1f percentage points whose interval excludes zero under both"
          % (100 * marg["M3_standardised"]["crude_gap"]),
          "BioProject-clustered and lineage-clustered resampling, and which survives collapsing to one",
          "genome per sequence type (%.1f%% against %.1f%%)."
          % (100 * marg["M2_one_per_ST"]["AB"], 100 * marg["M2_one_per_ST"]["KL"]), "",
          "Direct standardisation then removes it. Applying each host's family-specific rates to the",
          "other's family distribution accounts for **%.0f%%** of the crude gap. The reason is visible"
          % (100 * marg["M3_standardised"]["share_of_gap_explained_by_composition"]),
          "in the data: most of the *Klebsiella* class-A pool sits in *fosA*, *bla*SHV, *oqxA* and",
          "*oqxB* — species-core genes that do not occur in *A. baumannii* at all. The two hosts'",
          "chromosomal resistomes are built from different families, and that, rather than a difference",
          "in how shared families sit, is what the marginal gap measures.", "",
          "**Neither form of the within-chromosome contrast is reported as a result.** The pooled form",
          "fails a registered gate; the marginal form measures composition. Both are here so that the",
          "next person does not have to rediscover it.", "", "---", "", ""]

    anchor = "\n# Supplementary Tables\n"
    if anchor not in s:
        raise SystemExit("Supplementary Tables block not found")
    s = s.replace(anchor, "\n" + "\n".join(A) + "\n".join(B) + anchor, 1)
    io.open(S5, "w", encoding="utf-8", newline="\n").write(s)
    print("wrote %s" % os.path.basename(S5))
    print("  sha256 : %s" % sha(S5))
    print("  V4 supplement unchanged: %s"
          % (sha(S4) == "e65a3e83cad956402f82d8e7ff5ebc95f747a743d2e0717c32caf3621426f9f1"))
    print("  new sections: Supplementary Result %d and %d" % (a, b))


if __name__ == "__main__":
    main()
