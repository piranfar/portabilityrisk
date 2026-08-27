"""Freeze the lineage-adjustment analysis before any lineage-adjusted estimate exists.

Reviewer 1, major concern 2: the manuscript says "host-conditioned" while itself
recording that the cohort carries no MLST, no SNP distance and no lineage
assignment, so leave-one-major-lineage-out was reported NOT EVALUABLE. Species is
confounded with clone, and BioProject balancing is not a substitute: a BioProject
is a submission unit, a clone is a descent unit.

That gap is now closable. Sequence types are being called from the assemblies with
mlst 2.35.0 against PubMLST, which makes the NOT EVALUABLE arm evaluable for the
first time. This registers what will be done with them.

What was known when this was written:

  * the baseline: 58 matched families, MH OR 50.29 (95% CI 45.61-55.45), reproduced
    exactly from the deposited table;
  * the intrinsic-determinant sensitivity NM-V4C-001 has been run and survived
    (S8a 48.45, S8b 51.53, S8c 51.51, every interval excluding 1);
  * MLST is running. Sequence types are genotype, not outcome, so calling them
    reveals nothing about the estimand;
  * A. baumannii contributes 780 genomes and the Klebsiella group 3,460;
  * from a four-genome pilot, A. baumannii ST1 and ST2 are present - the two
    global clones - so lineage concentration in that host is expected and is
    exactly what the reviewer suspects.

What was NOT known: any odds ratio computed within, across or adjusted for
sequence type. No arm below has been run.
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
OUT = NM + "NM_V4C_LINEAGE_ADJUSTMENT_AMENDMENT_002.json"

BODY = {
  "amendment": "NM-V4C-002",
  "title": "Lineage adjustment of the matched-family host contrast using MLST",
  "raised_by": "Reviewer 1, major concern 2: host effects are not separated from "
               "lineage and database-sampling structure, and 'host-conditioned' "
               "asserts more than the design supports",
  "supersedes_status": "the manuscript reports leave-one-major-lineage-out as NOT "
                       "EVALUABLE because no lineage assignment existed. This "
                       "amendment makes it evaluable and commits to reporting the "
                       "result whichever way it falls",
  "frozen_before_any_outcome_was_scored": True,

  "typing": {
    "tool": "mlst 2.35.0, PubMLST schemes bundled with that release",
    "scheme_forced_per_host": {"Acinetobacter baumannii": "abaumannii_2 (Pasteur)",
                               "Klebsiella group": "klebsiella (KpSC, Pasteur)"},
    "why_forced": "A. baumannii has two PubMLST schemes. Auto-detection would pick "
                  "per genome and mix Oxford with Pasteur types, producing sequence "
                  "types that are not comparable with one another.",
    "untypeable": "a genome whose scheme returns ST '-' is untypeable. It is NOT "
                  "dropped by default; both handlings are registered as L1u and L2u "
                  "below, because dropping them silently would remove exactly the "
                  "divergent genomes least likely to belong to a dominant clone."
  },

  "estimator": "nmv4c_score.mh, Mantel-Haenszel with Robins-Breslow-Greenland "
               "variance, imported unchanged from the module that produced the "
               "baseline. Only the row filter or weighting differs between arms.",

  "arms": {
    "L1": {
      "name": "one genome per sequence type per host",
      "rule": "within each (host, ST), retain exactly one genome, chosen "
              "deterministically as the smallest SHA-256 of its assembly accession - "
              "the same tie-break rule the original freeze used for collapsing to "
              "one genome per BioProject",
      "purpose": "removes clonal replication entirely; the strongest single test"
    },
    "L1u": {"name": "L1, untypeable genomes excluded",
            "rule": "L1 after dropping every genome with ST '-'"},
    "L2": {
      "name": "leave-one-sequence-type-out",
      "rule": "drop each ST in turn from both hosts and recompute; report the "
              "maximum relative change in ln(OR)",
      "ceiling": 0.15,
      "why_this_ceiling": "it is the ceiling the original freeze set for "
                          "leave-one-BioProject-out. ST is the analogous clustering "
                          "unit, so the threshold is borrowed rather than chosen "
                          "now. Borrowing is stated explicitly because choosing a "
                          "threshold after seeing data is what freezing forbids."
    },
    "L2u": {"name": "L2, untypeable genomes treated as one pooled lineage",
            "rule": "every ST '-' genome is assigned to a single synthetic lineage "
                    "and that lineage is subject to leave-one-out like any other"},
    "L3": {
      "name": "lineage-balanced weighting",
      "rule": "each occurrence is weighted by 1/g, where g is the number of genomes "
              "of its own (host, ST) in the cohort, so a heavily sequenced clone "
              "contributes no more than a singleton",
      "note": "weights are applied to the 2x2 cell counts before pooling; the "
              "Mantel-Haenszel estimator is unchanged"
    }
  },

  "descriptive_reporting_required": [
    "the number of genomes, distinct STs, untypeable genomes and effective number "
    "of STs by inverse Herfindahl-Hirschman index, per host",
    "the share of each host's genomes held by its single largest ST",
    "these are reported whatever the arms show, because the reviewer also asked "
    "for cohort composition to be visible in the main text"
  ],

  "decision_rule": {
    "survives_fully": "every arm keeps direction with a 95% interval excluding 1, "
                      "AND leave-one-ST-out maximum relative change stays at or "
                      "below 0.15",
    "survives_partly": "arms keep direction and exclude 1 but leave-one-ST-out "
                       "exceeds 0.15. Then the effect is real but lineage-sensitive, "
                       "and the title must move from 'host-conditioned' to "
                       "'host-associated'",
    "fails": "any arm loses direction or covers 1. Then the host claim is withdrawn "
             "and rewritten as an association within this closed-genome cohort",
    "no_threshold_may_move": "these three outcomes and the 0.15 ceiling may not be "
                             "changed after any arm is computed"
  },

  "known_limits": [
    "MLST is a seven-locus scheme. Two genomes sharing an ST are not clones, and "
    "two genomes differing at one locus may be near-identical. ST is a lineage "
    "surrogate, not a phylogeny, and the result must be reported as such",
    "the Klebsiella group spans K. pneumoniae, K. quasipneumoniae and K. variicola; "
    "the KpSC scheme covers all three but ST numbering is shared, so an ST is "
    "interpreted within its species",
    "this arm covers the matched-family contrast only. The distance analysis and "
    "the genome-wide null are not lineage-adjusted by it",
    "no MLST scheme existed in this project before this amendment, so nothing here "
    "revises a previously reported lineage result - it replaces a NOT EVALUABLE"
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
    print("frozen: %s" % os.path.basename(OUT))
    print("  sha256_of_body : %s" % body_sha)
    print("  sha256_of_file : %s"
          % hashlib.sha256(text.encode("utf-8")).hexdigest())
    print("\narms registered: L1, L1u, L2, L2u, L3")
    print("decision rule has three named outcomes, including one that renames the title")
    print("no arm has been computed")


if __name__ == "__main__":
    main()
