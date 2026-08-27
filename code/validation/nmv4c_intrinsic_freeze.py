"""Freeze the intrinsic-determinant sensitivity BEFORE any of its results exist.

Reviewer 1 argues that the 50.29-fold matched-family contrast could partly reflect
determinants that are intrinsic chromosomal genes of the host rather than acquired
ones. The objection is legitimate and this registers the test of it.

What was already known when this design was written, stated in full because the
project's rule is that an amendment must say what was known at the time:

  * the baseline is published: 58 matched families, MH OR 50.29 (95% CI
    45.61-55.45), 56 families concordant, and it was reproduced exactly from the
    deposited occurrence table before this design was written;
  * the 58 family names were read. Five of the determinants the reviewer names -
    fosA, oqxA, oqxB, blaPDC, blaACT, blaLEN - are NOT among them. They exist in
    the cohort but fail the frozen eligibility rule of >=3 species, because a
    species-intrinsic gene is by definition not spread across species. The
    reviewer could not have known this;
  * two eligible families do carry intrinsic members: blaOXA, which in
    A. baumannii lumps the intrinsic chromosomal OXA-51-like enzyme with acquired
    OXA-23/24/58; and blaSHV, which in K. pneumoniae is partly chromosomal;
  * from AMRFinderPlus database 2026-08-07.1, 736 of the 1,750 blaOXA
    occurrences in A. baumannii carry the family node blaOXA-51_fam;
  * AMRFinderPlus does NOT subdivide blaSHV into intrinsic and acquired nodes, so
    no allele-level rule for SHV is derivable from the database and none is
    invented here.

What was NOT known: the odds ratio under any exclusion. No arm below has been
computed. That is the point of freezing this file first.
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
OUT = NM + "NM_V4C_INTRINSIC_SENSITIVITY_AMENDMENT_001.json"

BODY = {
  "amendment": "NM-V4C-001",
  "title": "Intrinsic-determinant sensitivity for the matched-family host contrast",
  "raised_by": "Reviewer 1, major concern 1: intrinsic versus acquired determinants "
               "are not separated, so host-specific architecture could reflect "
               "species-core chromosomal genes",
  "frozen_before_any_outcome_was_scored": True,
  "written_after": "the baseline was reproduced (58 families, MH OR 50.2912) and the "
                   "58 family names and their AMRFinderPlus family nodes were read",
  "written_before": "any sensitivity odds ratio was computed",

  "estimator": {
    "method": "Mantel-Haenszel across matched gene families, Robins-Breslow-Greenland "
              "variance, 95% Wald interval on ln(OR)",
    "source": "imported unchanged from nmv4c_score.mh, the module that produced the "
              "published baseline, so only the row filter differs between baseline "
              "and sensitivity",
    "concordance_count": "families with Woolf log odds ratio > 0 under the "
                         "Haldane-Anscombe 0.5 correction, as in the baseline"
  },

  "cohort": {
    "source_table": "portabilityrisk_occurrence_portability_v1.tsv (deposit v2)",
    "contrast": "Acinetobacter baumannii against the Klebsiella group",
    "endpoint": "two-way: class B against classes C, D and E; class A excluded",
    "family_eligibility": ">=3 species AND >=10 BioProjects AND >=20 occurrences AND "
                          "class B observable, unchanged from the original freeze"
  },

  "arms": {
    "S8a": {
      "name": "allele-level intrinsic exclusion",
      "rule": "drop every occurrence whose AMRFinderPlus family node is a recognised "
              "species-intrinsic chromosomal node of an Acinetobacter species",
      "nodes_excluded": ["blaOXA-51_fam", "blaOXA-134_fam", "blaOXA-213_fam"],
      "why_these": "blaOXA-51_fam is the intrinsic chromosomal class D beta-lactamase "
                   "of A. baumannii; blaOXA-134_fam and blaOXA-213_fam are the "
                   "corresponding intrinsic families of other Acinetobacter species. "
                   "Node membership is taken from AMRProt.fa in AMRFinderPlus "
                   "database 2026-08-07.1, the same version that produced the calls.",
      "unmapped_alleles": "bare 'blaOXA' calls carry no node and are RETAINED in this "
                          "arm; they are removed in S8b, which drops the family whole"
    },
    "S8b": {
      "name": "whole-family exclusion, conservative",
      "rule": "drop the blaOXA and blaSHV families entirely from the matched set",
      "why": "it needs no allele taxonomy at all and removes the acquired members "
             "along with the intrinsic ones, so it can only understate the effect. "
             "It is the arm that answers the objection without relying on any "
             "classification of ours."
    },
    "S8c": {
      "name": "species-core screen",
      "rule": "S8b, and additionally drop any matched family whose name matches a "
              "known species-core determinant of the four study hosts",
      "screen": ["fosA", "fosA3", "fosA5", "fosA6", "oqxA", "oqxB", "blaADC",
                 "blaPDC", "blaACT", "blaMIR", "blaLEN", "blaOKP", "blaOXY",
                 "catB7", "aph(3')-IIb"],
      "expected": "no change relative to S8b, because none of these except fosA3 is "
                  "among the 58. The arm is run to demonstrate that rather than to "
                  "assert it. fosA3 is the acquired plasmid-borne variant and is "
                  "screened here only because its name collides with intrinsic fosA."
    }
  },

  "decision_rule": {
    "survives": "the claim that a matched determinant family occupies different "
                "portability routes across hosts survives if, in EVERY arm, the MH "
                "odds ratio keeps its direction and its 95% interval excludes 1",
    "fails": "if any arm loses the direction or its interval covers 1, the headline "
             "claim is revised, not defended",
    "reported_regardless": "all three arms, with family counts, concordance and "
                           "interval, are reported whatever they show",
    "no_threshold_may_move": "these criteria may not be changed after the arms are "
                             "computed"
  },

  "known_limits": [
    "blaSHV cannot be split into intrinsic and acquired alleles from the database, "
    "so S8a does not attempt it and S8b removes the family whole",
    "this sensitivity addresses the matched-family contrast only. The distance "
    "analysis, the 46.39% chromosomal association and the genome-wide null are not "
    "family-matched and are not covered by these arms",
    "excluding a family shrinks the matched set, so a wider interval in an arm is "
    "expected and is not by itself evidence against the effect"
  ]
}


def main():
    if os.path.exists(OUT):
        raise SystemExit("amendment already exists - refusing to overwrite: %s" % OUT)
    body = json.dumps(BODY, indent=1, ensure_ascii=False, sort_keys=True)
    body_sha = hashlib.sha256(body.encode("utf-8")).hexdigest()
    doc = dict(BODY)
    doc["sha256_of_body"] = body_sha
    text = json.dumps(doc, indent=1, ensure_ascii=False, sort_keys=True) + "\n"
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(text)
    file_sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    print("frozen: %s" % os.path.basename(OUT))
    print("  sha256_of_body : %s" % body_sha)
    print("  sha256_of_file : %s" % file_sha)
    print("\narms registered: S8a allele-level, S8b whole-family, S8c species-core screen")
    print("no arm has been computed")


if __name__ == "__main__":
    main()
