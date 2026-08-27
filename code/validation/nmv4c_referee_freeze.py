"""Freeze the referee-requested analyses before any of them is computed.

The referee makes three claims that, if right, change what this paper is:

  * the headline two-way endpoint multiplies the host's chromosome-versus-plasmid
    propensity with marker-proximity-given-chromosomal, so a host extreme on both
    produces a large odds ratio partly by construction. The contrast that tests
    the TITLE - chromosomal mobility architecture - is class B against class A
    WITHIN the chromosome, and it has never been reported for A. baumannii
    against Klebsiella;
  * leave-one-sequence-type-out already breached the registered ceiling on ST2,
    which is GC2, the clone defined by AbaR resistance islands. The decisive
    analysis is the contrast computed in non-ST2 A. baumannii alone. If it
    collapses, the finding is "GC2 carries resistance islands", which is known;
  * the headline interval is analytic while every other interval in the paper is
    a BioProject cluster bootstrap, and the effective lineage count is 3.88. A
    +-10% interval on ~4 effective lineages is not defensible.

All three are testable with data already in hand. This registers what will be
run and what each outcome will require, before any of it is computed.

Known at the time of writing: baseline 50.29 (45.61-55.45); lineage-adjusted
21.76 (17.49-27.07); leave-one-ST-out 25.4% on ST2 against a 15% ceiling;
intrinsic sensitivity survives; exact-symbol 56.40; random effects 42.78 with
I2 83.6%. NOT known: any within-chromosome odds ratio for the headline pair, any
non-ST2 estimate, and any cluster-bootstrap interval for the headline.
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
OUT = NM + "NM_V4C_REFEREE_ANALYSES_AMENDMENT_004.json"

BODY = {
  "amendment": "NM-V4C-004",
  "title": "Within-chromosome endpoint, non-ST2 restriction, and cluster-bootstrap "
           "intervals for the headline contrast",
  "raised_by": "referee report, major points 1, 2 and 3",
  "frozen_before_any_outcome_was_scored": True,

  "arms": {
    "W1": {"name": "within-chromosome matched-family contrast, headline pair",
           "endpoint": "class B against class A, chromosomal occurrences only; "
                       "plasmid classes C, D and E excluded entirely",
           "contrast": "A. baumannii against the Klebsiella group",
           "why": "this is the contrast the title makes. The published two-way "
                  "endpoint asks a different question: whether a determinant is "
                  "chromosomal-and-mobile rather than plasmid-borne, which is "
                  "sensitive to how often the host puts determinants on plasmids "
                  "at all"},
    "W2": {"name": "W1, lineage-adjusted",
           "rule": "W1 restricted to one genome per sequence type, the same "
                   "deterministic selection as NM-V4C-002"},
    "N1": {"name": "headline two-way contrast, non-ST2 A. baumannii only",
           "rule": "drop every A. baumannii genome whose sequence type is 2; "
                   "Klebsiella unchanged",
           "why": "ST2 is GC2, the lineage defined by AbaR islands, and it "
                  "breached the registered leave-one-out ceiling"},
    "N2": {"name": "within-chromosome contrast, non-ST2 A. baumannii only",
           "rule": "W1 with the same ST2 exclusion"},
    "C1": {"name": "cluster-bootstrap interval for the headline two-way contrast",
           "rule": "resample BioProjects with replacement, B = 2,000, seed "
                   "20260827, recomputing the matched-family Mantel-Haenszel odds "
                   "ratio in each replicate; report the 2.5th and 97.5th "
                   "percentiles",
           "why": "every other interval in the paper is a cluster bootstrap; the "
                  "headline is not, and it is the one estimate where the "
                  "independence assumption is least defensible"},
    "C2": {"name": "cluster-bootstrap interval, resampling sequence types",
           "rule": "as C1 but resampling (host, ST) lineages, which is the "
                   "clustering the leave-one-out arm showed to matter"}
  },

  "decision_rule": {
    "primary_becomes": "if W1 is estimable, the within-chromosome odds ratio "
                       "becomes the primary result for the title's claim and the "
                       "two-way estimate is reported beside it as a different "
                       "question, not as the headline",
    "abstract_must_change_if": "W1 is materially smaller than the two-way "
                               "estimate, defined in advance as a point estimate "
                               "below half of 50.29",
    "finding_is_reframed_if": "N1 or N2 loses direction or its interval covers 1, "
                              "in which case the result is a property of ST2/GC2 "
                              "and the paper must say so in the title and abstract",
    "interval_replaced_if": "the cluster-bootstrap interval from C1 or C2 is wider "
                            "than the analytic one, in which case the widest of "
                            "the three becomes the reported interval",
    "reported_regardless": "every arm is reported whatever it shows, including any "
                           "arm that is not estimable and why",
    "no_threshold_may_move": "these criteria may not be changed after any arm runs"
  },

  "known_limits": [
    "W1 conditions on being chromosomal, so it answers the architecture question "
    "and NOT the question of how much of a host's resistome is chromosomal; both "
    "are needed and neither substitutes for the other",
    "removing ST2 removes roughly half the A. baumannii genomes, so N1 and N2 will "
    "have wider intervals; width is expected and is not itself evidence against "
    "the effect",
    "a cluster bootstrap over ~4 effective lineages cannot manufacture precision "
    "that the design does not contain. A very wide interval here is an honest "
    "result, not a failed analysis"
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
    print("  arms: W1 W2 (within-chromosome), N1 N2 (non-ST2), C1 C2 (cluster bootstrap)")
    print("  the rule names, in advance, the outcome that renames the paper")
    print("  no arm has been computed")


if __name__ == "__main__":
    main()
