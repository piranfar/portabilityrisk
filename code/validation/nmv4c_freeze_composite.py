"""Freeze the composite-element test before computing it (referee point 4).

The objection: blaOXA-23 travels in Tn2006/Tn2008, delivered by ISAba1;
aph(3'')/aph(6) travel in Tn5393. Where an ARG is the cargo of a composite
transposon, its distance to the nearest insertion sequence is near zero BY
CONSTRUCTION, because the flanking copies are what make it a transposon. A null
that relocates the gene uniformly within the chromosome cannot address that: it
asks whether the gene sits near an IS more often than chance, when the answer is
definitional for that subset. If most of the A. baumannii signal comes from
cargo occurrences, the 16.91-fold enrichment is a restatement of composite
transposon structure and not a statement about genome architecture.

The test is a stratification, not a new estimator. The published null matrices
are reused unchanged, so observed and expected are computed exactly as in the
primary analysis and only the subset differs.

Known when this was written: the published enrichment is 16.91-fold at 1 kb for
A. baumannii, 2.69 times the Klebsiella group and 1.86 times P. aeruginosa; the
census holds 190,999 elements of which 145,779 are complete-structural. NOT
known: how many occurrences are composite-flanked, or the enrichment in either
stratum.
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
OUT = NM + "NM_C1_COMPOSITE_ELEMENT_AMENDMENT_006.json"

BODY = {
  "amendment": "NM-C1-006",
  "title": "Composite-element stratification of the genome-wide background null",
  "raised_by": "referee point 4: the IS-proximity enrichment is partly tautological "
               "where the resistance gene is the cargo of a composite transposon",
  "frozen_before_any_outcome_was_scored": True,

  "definition": {
    "composite_flanked": (
        "a chromosomal resistance-gene occurrence is composite-flanked when, on its "
        "own replicon, there is at least one complete-structural insertion sequence "
        "of family F lying entirely upstream of the gene within D bp of its start, "
        "AND at least one complete-structural insertion sequence of the SAME family "
        "F lying entirely downstream of the gene within D bp of its end"),
    "why_same_family": (
        "a composite transposon is bounded by two copies of the same element. "
        "Requiring the same family is what separates 'flanked by a pair' from "
        "'sits in an IS-dense neighbourhood', and the second is the thing the "
        "paper claims to measure"),
    "complete_structural": (
        "unchanged from the primary analysis: ISEScan type c, complete transposase "
        "open reading frame, bilateral resolved terminal inverted repeats"),
    "D_primary": 10000,
    "D_secondary": 5000,
    "not_claimed": (
        "this identifies a flanking ARRANGEMENT consistent with a composite "
        "element. It does not demonstrate that the element is intact, mobile, or "
        "that it moved. Tn2006 is not called by name here and no transposon "
        "database is consulted")
  },

  "arms": {
    "P1": {"name": "enrichment among composite-flanked occurrences",
           "estimator": "unchanged; observed and expected read from the saved "
                        "null matrices c1_4_matrices.npz, subset by column"},
    "P2": {"name": "enrichment among occurrences that are NOT composite-flanked"},
    "P3": {"name": "the host contrast within the non-flanked stratum",
           "quantity": "A. baumannii enrichment divided by the Klebsiella group's, "
                       "and by P. aeruginosa's, at 1 kb"},
    "P4": {"name": "the same stratification at D = 5 kb"}
  },

  "decision_rule": {
    "enrichment_is_substantially_composite_structure_if": (
        "A. baumannii's enrichment among non-flanked occurrences falls below half "
        "the published 16.91, i.e. below 8.46"),
    "host_contrast_is_composite_structure_if": (
        "the A. baumannii / Klebsiella enrichment ratio in the non-flanked stratum "
        "falls below 1.5, against the published 2.69, or loses direction"),
    "objection_answered_if": (
        "both the non-flanked enrichment and the non-flanked host ratio stay above "
        "those values, in which case the signal is not carried by cargo occurrences "
        "and the referee's objection does not hold for this cohort"),
    "reported_regardless": "both strata, both distances, and the number of "
                           "occurrences in each, whatever they show",
    "no_threshold_may_move": "8.46 and 1.5 are fixed here and may not be changed "
                             "after any arm runs"
  },

  "known_limits": [
    "a flanking pair of the same family is necessary but not sufficient for a "
    "composite transposon; some flanked occurrences will be coincidental and some "
    "true composites will be missed where one copy is degenerate and fails the "
    "complete-structural filter",
    "the two strata are not independent samples of one population: cargo status is "
    "a property of the determinant as much as of the host",
    "this tests the genome-wide null only. It says nothing about the matched-family "
    "contrasts"
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
    print("  thresholds fixed in advance: 8.46 enrichment, 1.5 host ratio")
    print("  no arm has been computed")


if __name__ == "__main__":
    main()
