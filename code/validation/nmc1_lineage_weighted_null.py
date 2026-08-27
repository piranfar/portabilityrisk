"""A lineage-balanced weighting for the background null (item 4.2).

The paper argues that BioProject balancing is not clonal-lineage adjustment - a
BioProject is a submission unit and a clone is a descent unit - and then types
every genome in the headline contrast. But the four registered weightings for the
enrichment are occurrence, one-per-block, genome-balanced and BioProject-balanced.
None of them is lineage-balanced. The non-ST2 stratification already reported is a
stratification, not a weighting.

So the paper's own argument had not been applied to the paper's own primary
result. This applies it.

Scope, stated before the arm runs: sequence types exist only for the 4,240 genomes
of the primary contrast. A. baumannii and the Klebsiella group can be
lineage-balanced; P. aeruginosa and Enterobacter cannot, because they were never
typed. Those two are reported as NOT EVALUABLE rather than quietly omitted.

Threshold reused rather than chosen: the enrichment must stay above 8.46 and the
A. baumannii / Klebsiella ratio above 1.5, exactly as registered for the
composite-element and non-ST2 arms.
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
AMEND = NM + "NM_C1_LINEAGE_WEIGHTED_AMENDMENT_008.json"
MATS = _dir("PORTABILITYRISK_C1_DIR") + "c1_4_matrices.npz"
MLST = NM + "NM_V4C_MLST_CALLS_V1.tsv"
OUT = NM + "NM_C1_LINEAGE_WEIGHTED_RESULTS_V1.json"

FLOOR, RATIO_FLOOR = 8.46, 1.5

BODY = {
  "amendment": "NM-C1-008",
  "title": "Lineage-balanced weighting of the genome-wide background null",
  "raised_by": "the paper argues BioProject balancing is not lineage adjustment, then does not "
               "offer a lineage-balanced weighting for its own primary enrichment",
  "frozen_before_any_outcome_was_scored": True,
  "arm": {"name": "lineage-balanced",
          "rule": "each chromosomal occurrence is weighted by 1/k, where k is the number of "
                  "occurrences sharing its (host, sequence type). A heavily sequenced clone "
                  "then contributes the same total weight as a singleton lineage.",
          "estimator": "unchanged; observed and expected read from the saved null matrices "
                       "c1_4_matrices.npz, only the weights differ",
          "untypeable": "a genome whose scheme returns ST '-' is its own singleton lineage, "
                        "not pooled, because pooling them would create one artificial "
                        "mega-lineage out of the most divergent genomes"},
  "scope_limit": "sequence types exist only for the 4,240 genomes of the primary contrast. "
                 "P. aeruginosa and Enterobacter were never typed and are reported NOT "
                 "EVALUABLE for this weighting.",
  "decision_rule": {
      "thresholds_reused_from": "NM-C1-006 and NM-C1-007, deliberately, rather than chosen now",
      "effect_is_lineage_artefact_if": "A. baumannii lineage-balanced enrichment falls below "
                                       "%s, or the A. baumannii / Klebsiella ratio below %s"
                                       % (FLOOR, RATIO_FLOOR),
      "reported_regardless": "the arm and the two NOT EVALUABLE groups, whatever it shows",
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
        if s.startswith("Enterobacter"):
            return "Enterobacter group"
        return None

    g = np.array([grp(s) for s in species])
    lin = np.array(["%s|%s" % (g[i], st[assembly[i]] if st.get(assembly[i], "-") != "-"
                               else "UNTYPED:" + assembly[i])
                    if assembly[i] in st else "" for i in range(len(g))])
    typed = lin != ""
    print("occurrences: %s | typed: %s (%.1f%%)"
          % ("{:,}".format(len(g)), "{:,}".format(int(typed.sum())),
             100 * typed.mean()))

    size = collections.Counter(lin[typed])
    wt = np.array([1.0 / size[lin[i]] if typed[i] else 0.0 for i in range(len(g))])

    def enrich(mask, weights, d=1000):
        if mask.sum() == 0 or weights[mask].sum() == 0:
            return None
        w = weights[mask]
        o = float((((obs[mask] >= 0) & (obs[mask] <= d)) * w).sum() / w.sum())
        nn = ((null[:, mask] >= 0) & (null[:, mask] <= d)).astype(np.float64)
        e = float((nn * w).sum() / (nn.shape[0] * w.sum()))
        return {"n": int(mask.sum()), "effective_lineages": len(set(lin[mask])),
                "observed": round(o, 4), "expected": round(e, 4),
                "enrichment": round(o / e, 2) if e else None}

    res = {}
    print("\n  %-22s %9s %8s %10s %10s %11s"
          % ("group", "n", "lineages", "observed", "expected", "enrichment"))
    for name in ("A. baumannii", "Klebsiella group", "P. aeruginosa", "Enterobacter group"):
        m = (g == name)
        if not (m & typed).any():
            res[name] = {"status": "NOT EVALUABLE",
                         "why": "no sequence types were called for this group"}
            print("  %-22s %9s   NOT EVALUABLE - never typed" % (name, "{:,}".format(int(m.sum()))))
            continue
        r = enrich(m & typed, wt)
        res[name] = r
        print("  %-22s %9s %8d %10.4f %10.4f %11.2f"
              % (name, "{:,}".format(r["n"]), r["effective_lineages"],
                 r["observed"], r["expected"], r["enrichment"]))

    ab = res["A. baumannii"]["enrichment"]
    kl = res["Klebsiella group"]["enrichment"]
    ratio = ab / kl
    print("\n  lineage-balanced A. baumannii enrichment  %.2f  (floor %.2f)" % (ab, FLOOR))
    print("  lineage-balanced A. baumannii / Klebsiella %.2f  (floor %.2f)" % (ratio, RATIO_FLOOR))
    artefact = ab < FLOOR or ratio < RATIO_FLOOR
    print("\nfrozen decision rule -> effect is a lineage artefact: %s" % artefact)

    rec = {"receipt": "NM-C1-008 lineage-balanced weighting of the background null",
           "amendment_sha256_of_body": claimed, "amendment_sha256_of_file": sha(AMEND),
           "inputs": {"null_matrices": sha(MATS), "mlst_calls": sha(MLST)},
           "results": res,
           "summary": {"AB_enrichment": ab, "AB_over_Klebsiella": round(ratio, 3)},
           "decision_rule": doc["decision_rule"],
           "effect_is_lineage_artefact": bool(artefact)}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("receipt: %s" % os.path.basename(OUT))


if __name__ == "__main__":
    main()
