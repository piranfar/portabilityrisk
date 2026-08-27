"""Freeze the within-chromosome marker-adjacency estimand before computing it.

The matched-family Mantel-Haenszel contrast cannot answer the within-chromosome
question at this project's registered standard: one family carries 64% of the
weight, only 9-11 of 46-49 families carry any weight at all, and the sign moves
between 0.48 and 6.11 with the handling of blaOXA alone.

The quantity that does not move is the marginal one: among chromosomal acquired
resistance-gene occurrences, the share with a mobile-element marker inside their
own +-10 kb window. That is a rate, not a pooled odds ratio, and no family can
flip it.

It buys robustness by giving something up, and the give-up is registered here
rather than discovered later: a marginal rate is NOT adjusted for gene-family
composition. If A. baumannii's chromosomal resistome is concentrated in families
that are intrinsically IS-associated, the marginal gap counts that as part of
the effect. Direct standardisation is therefore run alongside, precisely to
measure how much of the gap is composition, and BOTH are reported.

Known when this was written: the crude rates are 80.9% for A. baumannii and 34.2%
for the Klebsiella group; the matched-family arms above; and the lineage
structure (A. baumannii effective 3.88 lineages, ST2 at 48.2%). NOT known: any
standardised rate, any interval on any of these, and any lineage-adjusted rate.
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
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

NM = _dir("PORTABILITYRISK_REPO_DIR") + "docs/nature_microbiology/"
OUT = NM + "NM_V4C_MARGINAL_ADJACENCY_AMENDMENT_005.json"

BODY = {
  "amendment": "NM-V4C-005",
  "title": "Within-chromosome marker adjacency as a marginal rate, with direct "
           "standardisation for gene-family composition",
  "why_this_replaces_the_matched_contrast": (
      "the matched-family odds ratio for the within-chromosome endpoint fails the "
      "weight-concentration gate G2 registered in NMV4C_FROZEN_DESIGN.json in "
      "every variant, and its sign depends on one family. A marginal rate is "
      "estimable where the pooled odds ratio is not."),
  "frozen_before_any_outcome_was_scored": True,
  "post_hoc_and_labelled_as_such": (
      "this estimand was chosen AFTER seeing that the matched-family arm was "
      "unstable. That is a post-hoc change of estimand and is labelled post hoc "
      "wherever it appears. What is not post hoc is the direction of the "
      "conclusion: the rule below fixes what each outcome requires before any "
      "interval is computed."),

  "estimand": {
    "definition": "among chromosomal acquired resistance-gene occurrences of a "
                  "host, the proportion whose own +-10 kb window contains at "
                  "least one mobile-element marker; that is class B / (class A + "
                  "class B)",
    "unit": "occurrence",
    "hosts": ["Acinetobacter baumannii", "Klebsiella group"],
    "not_adjusted_for": "gene-family composition, by construction"
  },

  "arms": {
    "M1": {"name": "crude rate per host",
           "intervals": "percentile bootstrap, B = 2,000, seed 20260827, "
                        "resampling BioProjects; and separately resampling "
                        "(host, sequence type) lineages"},
    "M2": {"name": "rate under one genome per sequence type",
           "rule": "the deterministic selection registered in NM-V4C-002"},
    "M3": {"name": "directly standardised rates",
           "rule": "apply each host's family-specific rates to the OTHER host's "
                   "family distribution, over families present in both, and "
                   "report both directions",
           "purpose": "to measure how much of the crude gap is family "
                      "composition rather than within-family difference. This is "
                      "the caveat quantified, not a second estimate of the effect"},
    "M4": {"name": "rate difference and rate ratio",
           "rule": "A. baumannii minus the Klebsiella group, and their ratio, "
                   "with the same two bootstrap clusterings"}
  },

  "decision_rule": {
    "claim_supported": "the within-chromosome claim is supported if the rate "
                       "difference keeps its direction with a 95% interval "
                       "excluding zero under BOTH clusterings AND under the "
                       "one-genome-per-sequence-type arm",
    "claim_is_compositional_if": "direct standardisation removes more than half "
                                 "the crude gap, in which case the manuscript "
                                 "must say the difference is substantially a "
                                 "difference in which families each host carries, "
                                 "not in how those families sit",
    "claim_withdrawn_if": "the difference loses direction or its interval covers "
                          "zero under any registered clustering",
    "reported_regardless": "all four arms, both standardisation directions, and "
                           "the failed matched-family arms are reported. Removing "
                           "an analysis because it failed is the omission this "
                           "manuscript has already been criticised for",
    "no_threshold_may_move": "these criteria may not be changed after any arm runs"
  },

  "known_limits": [
    "a marginal rate cannot separate architecture from composition; that is what "
    "M3 measures and why M3 is not optional",
    "the two hosts differ in how much of their resistome is chromosomal at all "
    "(86% against 34%), so these rates are computed on different-sized and "
    "differently-composed denominators",
    "class B means a marker within 10 kb. It is not evidence that the element "
    "mobilises the gene"
  ]
}


def main():
    if os.path.exists(OUT):
        raise SystemExit("amendment exists - refusing to overwrite")
    body = json.dumps(BODY, indent=1, ensure_ascii=False, sort_keys=True)
    h = hashlib.sha256(body.encode("utf-8")).hexdigest()
    doc = dict(BODY)
    doc["sha256_of_body"] = h
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(doc, indent=1, ensure_ascii=False, sort_keys=True) + "\n")
    print("frozen: %s" % os.path.basename(OUT))
    print("  sha256_of_body : %s" % h)
    print("  arms: M1 crude, M2 lineage-adjusted, M3 standardised, M4 difference")
    print("  the rule names in advance the outcome that calls this compositional")
    print("  no arm has been computed")


if __name__ == "__main__":
    main()
