"""Run the frozen intrinsic-determinant sensitivity (NM-V4C-001).

Refuses to run unless the amendment's own body hash still matches, so the design
cannot be edited after seeing an arm and re-run as if it had not been.

The estimator is imported from nmv4c_score, the module that produced the published
baseline. Only the row filter differs between the baseline and any arm, so a change
in the odds ratio cannot come from a change in the estimator.

The baseline is recomputed first. If it does not reproduce 58 families and
50.2912, the pipeline is wrong and no arm is reported.
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


import csv
import hashlib
import io
import json
import math
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = _dir("PORTABILITYRISK_REPO_DIR")
NM = REPO + "docs/nature_microbiology/"
AMEND = NM + "NM_V4C_INTRINSIC_SENSITIVITY_AMENDMENT_001.json"
TABLE = (_dir("PORTABILITYRISK_DEPOSIT_DIR")
         "portabilityrisk_occurrence_portability_v1.tsv")
AMRPROT = REPO + ".cache/model2_g_a/amrfinder_db/2026-08-07.1/AMRProt.fa"
OUT = NM + "NM_V4C_INTRINSIC_SENSITIVITY_RESULTS_V1.json"

sys.path.insert(0, REPO + "audit/ingest/assay_aware_emergence/v2/nm_validation")
from nmv4c_score import mh, woolf                                  # noqa: E402

BASE_FAMILIES, BASE_OR = 58, 50.2912


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_amendment():
    doc = json.load(io.open(AMEND, encoding="utf-8"))
    claimed = doc.pop("sha256_of_body")
    body = json.dumps(doc, indent=1, ensure_ascii=False, sort_keys=True)
    actual = hashlib.sha256(body.encode("utf-8")).hexdigest()
    if actual != claimed:
        raise SystemExit("REFUSING: the amendment body has changed since it was frozen\n"
                         "  frozen %s\n  now    %s" % (claimed, actual))
    doc["sha256_of_body"] = claimed
    print("frozen design verified: body %s" % claimed[:16])
    return doc


def node_map():
    """allele symbol -> AMRFinderPlus family node, from the pinned database."""
    m = {}
    for line in io.open(AMRPROT, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            p = line[1:].rstrip("\n").split("|")
            if len(p) > 5 and p[3]:
                m.setdefault(p[3], p[4])
    return m


def group(o):
    if o == "Acinetobacter baumannii":
        return "AB"
    if o.startswith("Klebsiella"):
        return "KL"
    return None


def eligible(rows):
    fam = {}
    for r in rows:
        d = fam.setdefault(r["gene_family"], {"sp": set(), "bp": set(), "n": 0, "B": 0})
        d["sp"].add(r["organism_harmonized"])
        d["bp"].add(r["bioproject_accession"])
        d["n"] += 1
        if r["portability_class"] == "B":
            d["B"] += 1
    return {f for f, d in fam.items()
            if len(d["sp"]) >= 3 and len(d["bp"]) >= 10 and d["n"] >= 20 and d["B"] > 0}


def contrast(rows, elig, keep_row=lambda r: True, drop_families=frozenset()):
    tab = {}
    for r in rows:
        g = group(r["organism_harmonized"])
        if g is None or r["portability_class"] == "A":
            continue
        f = r["gene_family"]
        if f not in elig or f in drop_families or not keep_row(r):
            continue
        cell = tab.setdefault(f, {"AB": [0, 0], "KL": [0, 0]})
        cell[g][0 if r["portability_class"] == "B" else 1] += 1
    tables, fams, conc = [], [], 0
    for f, c in sorted(tab.items()):
        a, b = c["AB"]
        cc, d = c["KL"]
        if (a + b) == 0 or (cc + d) == 0:
            continue
        tables.append((a, b, cc, d))
        fams.append(f)
        if woolf(a, b, cc, d)[0] > 0:
            conc += 1
    m = mh(tables)
    if m is None:
        return {"n_families": len(fams), "estimable": False}
    return {"n_families": len(fams), "estimable": True,
            "or": round(m["or"], 4),
            "ci_lo": round(m["ci_lo"], 4), "ci_hi": round(m["ci_hi"], 4),
            "excludes_1": bool(m["ci_lo"] > 1 or m["ci_hi"] < 1),
            "direction": "A. baumannii" if m["or"] > 1 else "Klebsiella",
            "families_concordant": conc,
            "families": fams}


def main():
    doc = load_amendment()
    rows = list(csv.DictReader(io.open(TABLE, encoding="utf-8", newline=""),
                              delimiter="\t"))
    nodes = node_map()
    elig = eligible(rows)
    print("occurrences %d | eligible families %d | alleles with a node %d"
          % (len(rows), len(elig), len(nodes)))

    base = contrast(rows, elig)
    print("\nbaseline reproduction")
    print("  families %d (expected %d) | OR %.4f (expected %.4f) | CI %.2f-%.2f | concordant %d"
          % (base["n_families"], BASE_FAMILIES, base["or"], BASE_OR,
             base["ci_lo"], base["ci_hi"], base["families_concordant"]))
    if base["n_families"] != BASE_FAMILIES or abs(base["or"] - BASE_OR) > 0.001:
        raise SystemExit("baseline did not reproduce - no arm will be reported")

    A = doc["arms"]
    drop_nodes = set(A["S8a"]["nodes_excluded"])
    s8a = contrast(rows, elig,
                   keep_row=lambda r: nodes.get(r["determinant_name"]) not in drop_nodes)
    s8b = contrast(rows, elig, drop_families={"blaOXA", "blaSHV"})
    screen = set(A["S8c"]["screen"])
    s8c = contrast(rows, elig, drop_families={"blaOXA", "blaSHV"} | screen)

    print("\n%-38s %9s %6s %-18s %s" % ("arm", "families", "OR", "95% CI", "concordant"))
    res = {"baseline": base, "S8a": s8a, "S8b": s8b, "S8c": s8c}
    for k, label in (("baseline", "baseline, as published"),
                     ("S8a", "S8a allele-level intrinsic removed"),
                     ("S8b", "S8b blaOXA and blaSHV removed whole"),
                     ("S8c", "S8c + species-core screen")):
        v = res[k]
        print("%-38s %9d %6.2f %-18s %d"
              % (label, v["n_families"], v["or"],
                 "%.2f-%.2f" % (v["ci_lo"], v["ci_hi"]), v["families_concordant"]))

    survives = all(res[k]["excludes_1"] and res[k]["direction"] == "A. baumannii"
                   for k in ("S8a", "S8b", "S8c"))
    print("\nfrozen decision rule: every arm must keep direction and exclude 1")
    print("verdict: %s" % ("SURVIVES" if survives else "FAILS - the headline claim "
                           "must be revised, not defended"))

    rec = {"receipt": "NM-V4C-001 intrinsic-determinant sensitivity",
           "amendment_sha256_of_body": doc["sha256_of_body"],
           "amendment_sha256_of_file": sha_file(AMEND),
           "inputs": {"occurrence_table": sha_file(TABLE),
                      "amrfinder_AMRProt_fa": sha_file(AMRPROT),
                      "amrfinder_db_version": "2026-08-07.1"},
           "estimator": "nmv4c_score.mh, Mantel-Haenszel with RBG variance, imported "
                        "unchanged from the module that produced the baseline",
           "baseline_reproduced": True,
           "results": {k: {kk: vv for kk, vv in v.items() if kk != "families"}
                       for k, v in res.items()},
           "families_by_arm": {k: v["families"] for k, v in res.items()},
           "decision_rule": doc["decision_rule"],
           "verdict": "SURVIVES" if survives else "FAILS"}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("\nreceipt: %s" % os.path.basename(OUT))


if __name__ == "__main__":
    main()
