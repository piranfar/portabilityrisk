"""Does the genome-wide IS-proximity enrichment survive removing GC2?

The decisive test under the owner's decision framework, and the one gap left in
this programme. The matched-family contrast has been run without ST2 (18.56).
The enrichment has been run without composite elements (13.58). It has never been
run without ST2, and ST2 is 48.2% of the A. baumannii genomes.

Design frozen here, run below, thresholds fixed before either.

Thresholds are the ones already registered for the composite-element test, reused
rather than chosen again: enrichment must stay above 8.46, which is half the
published 16.91, and the A. baumannii / Klebsiella ratio above 1.5 against the
published 2.69. Reusing a threshold from a sibling arm is deliberate - picking a
fresh one now, knowing what the composite arm returned, is the thing a freeze
exists to prevent.

Known when written: 16.91 published; 13.58 among non-composite-flanked; ST2 is
48.2% of A. baumannii genomes and breached the matched-family leave-one-out
ceiling. NOT known: any enrichment computed on non-ST2 genomes.
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


import collections
import csv
import hashlib
import io
import json
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

NM = _dir("PORTABILITYRISK_REPO_DIR") + "docs/nature_microbiology/"
AMEND = NM + "NM_C1_NON_ST2_ENRICHMENT_AMENDMENT_007.json"
MATS = _dir("PORTABILITYRISK_C1_DIR") + "c1_4_matrices.npz"
MLST = NM + "NM_V4C_MLST_CALLS_V1.tsv"
OUT = NM + "NM_C1_NON_ST2_ENRICHMENT_RESULTS_V1.json"

PUBLISHED_AB, FLOOR, RATIO_FLOOR = 16.91, 8.46, 1.5

BODY = {
  "amendment": "NM-C1-007",
  "title": "Genome-wide enrichment with the GC2 lineage removed",
  "raised_by": "the owner's decision framework: whether the effect holds in "
               "non-ST2 A. baumannii decides where this paper goes",
  "frozen_before_any_outcome_was_scored": True,
  "arms": {"E1": "enrichment at 1 kb for A. baumannii genomes whose sequence type "
                 "is not 2, with Klebsiella and P. aeruginosa unchanged",
           "E2": "E1 restricted further to occurrences that are not "
                 "composite-flanked, so both alternative explanations are removed "
                 "at once"},
  "estimator": "unchanged; observed and expected read from the saved null "
               "matrices c1_4_matrices.npz, subset by column",
  "decision_rule": {
      "effect_is_GC2_specific_if": "non-ST2 A. baumannii enrichment falls below "
                                   "%s, or the A. baumannii / Klebsiella ratio "
                                   "falls below %s" % (FLOOR, RATIO_FLOOR),
      "thresholds_reused_from": "NM-C1-006, deliberately, rather than chosen now",
      "reported_regardless": "both arms and the occurrence counts, whatever they show",
      "no_threshold_may_move": True}
}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def freeze():
    if os.path.exists(AMEND):
        return
    body = json.dumps(BODY, indent=1, ensure_ascii=False, sort_keys=True)
    d = dict(BODY)
    d["sha256_of_body"] = hashlib.sha256(body.encode("utf-8")).hexdigest()
    io.open(AMEND, "w", encoding="utf-8", newline="\n").write(
        json.dumps(d, indent=1, ensure_ascii=False, sort_keys=True) + "\n")
    print("frozen: %s  body %s" % (os.path.basename(AMEND), d["sha256_of_body"][:16]))


def main():
    freeze()
    doc = json.load(io.open(AMEND, encoding="utf-8"))
    claimed = doc.pop("sha256_of_body")
    if hashlib.sha256(json.dumps(doc, indent=1, ensure_ascii=False,
                                 sort_keys=True).encode("utf-8")).hexdigest() != claimed:
        raise SystemExit("REFUSING: amendment body changed since freeze")

    z = np.load(MATS, allow_pickle=False)
    species, assembly = z["species"], z["assembly"]
    obs, null = z["obs"].astype(np.int32), z["null"].astype(np.int32)
    st = {r["assembly_version"]: r["st"].strip()
          for r in csv.DictReader(io.open(MLST, encoding="utf-8", newline=""),
                                  delimiter="\t")}

    def grp(s):
        if s == "Acinetobacter baumannii":
            return "A. baumannii"
        if s.startswith("Klebsiella"):
            return "Klebsiella group"
        if s == "Pseudomonas aeruginosa":
            return "P. aeruginosa"
        return None

    g = np.array([grp(s) for s in species])
    is_ab = g == "A. baumannii"
    ab_st = np.array([st.get(a, "?") for a in assembly])
    st2 = is_ab & (ab_st == "2")
    typed = is_ab & np.isin(ab_st, list({v for v in st.values()}))
    print("A. baumannii chromosomal occurrences: %s | on ST2 genomes: %s (%.1f%%)"
          % ("{:,}".format(int(is_ab.sum())), "{:,}".format(int(st2.sum())),
             100 * st2.sum() / is_ab.sum()))
    print("  untyped A. baumannii occurrences (no MLST call): %s"
          % "{:,}".format(int((is_ab & (ab_st == "?")).sum())))

    def enrich(mask, d=1000):
        if mask.sum() == 0:
            return None
        o = float(((obs[mask] >= 0) & (obs[mask] <= d)).mean())
        e = float(((null[:, mask] >= 0) & (null[:, mask] <= d)).mean())
        return {"n": int(mask.sum()), "observed": round(o, 4),
                "expected": round(e, 4), "enrichment": round(o / e, 2) if e else None}

    rows = [("A. baumannii, all", is_ab),
            ("A. baumannii, ST2 only", st2),
            ("A. baumannii, non-ST2", is_ab & ~st2),
            ("Klebsiella group", g == "Klebsiella group"),
            ("P. aeruginosa", g == "P. aeruginosa")]
    res = {}
    print("\n  %-26s %9s %10s %10s %12s"
          % ("stratum", "n", "observed", "expected", "enrichment"))
    for lab, m in rows:
        r = enrich(m)
        res[lab] = r
        print("  %-26s %9s %10.4f %10.4f %12.2f"
              % (lab, "{:,}".format(r["n"]), r["observed"], r["expected"],
                 r["enrichment"]))

    ab_nst2 = res["A. baumannii, non-ST2"]["enrichment"]
    ratio_kl = ab_nst2 / res["Klebsiella group"]["enrichment"]
    ratio_pa = ab_nst2 / res["P. aeruginosa"]["enrichment"]
    print("\n  non-ST2 A. baumannii enrichment      %.2f   (published all: %.2f, floor %.2f)"
          % (ab_nst2, PUBLISHED_AB, FLOOR))
    print("  non-ST2 / Klebsiella                 %.2f   (published 2.69, floor %.2f)"
          % (ratio_kl, RATIO_FLOOR))
    print("  non-ST2 / P. aeruginosa              %.2f   (published 1.86)" % ratio_pa)

    gc2 = ab_nst2 < FLOOR or ratio_kl < RATIO_FLOOR
    print("\nfrozen decision rule -> effect is GC2-specific: %s" % gc2)

    rec = {"receipt": "NM-C1-007 enrichment with GC2 removed",
           "amendment_sha256_of_body": claimed, "amendment_sha256_of_file": sha(AMEND),
           "inputs": {"null_matrices": sha(MATS), "mlst_calls": sha(MLST)},
           "results": res,
           "non_ST2_summary": {"enrichment": ab_nst2,
                               "over_Klebsiella": round(ratio_kl, 3),
                               "over_Pseudomonas": round(ratio_pa, 3)},
           "decision_rule": doc["decision_rule"],
           "effect_is_GC2_specific": bool(gc2)}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("receipt: %s" % os.path.basename(OUT))


if __name__ == "__main__":
    main()
