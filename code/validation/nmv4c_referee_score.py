"""Run the frozen referee analyses (NM-V4C-004)."""
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
import random
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = "" + _dir("PORTABILITYRISK_REPO_DIR") + ""
NM = REPO + "docs/nature_microbiology/"
AMEND = NM + "NM_V4C_REFEREE_ANALYSES_AMENDMENT_004.json"
MLST = NM + "NM_V4C_MLST_CALLS_V1.tsv"
TABLE = ("" + _dir("PORTABILITYRISK_DEPOSIT_DIR") + ""
         "portabilityrisk_occurrence_portability_v1.tsv")
OUT = NM + "NM_V4C_REFEREE_ANALYSES_RESULTS_V1.json"

sys.path.insert(0, REPO + "audit/ingest/assay_aware_emergence/v2/nm_validation")
from nmv4c_score import mh, woolf                                  # noqa: E402

BASE_N, BASE_OR = 58, 50.2912
B_BOOT, SEED = 2000, 20260827


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def host(o):
    if o == "Acinetobacter baumannii":
        return "AB"
    if o.startswith("Klebsiella"):
        return "KL"
    return None


def eligible(rows):
    agg = {}
    for r in rows:
        d = agg.setdefault(r["gene_family"], {"sp": set(), "bp": set(), "n": 0, "B": 0})
        d["sp"].add(r["organism_harmonized"])
        d["bp"].add(r["bioproject_accession"])
        d["n"] += 1
        if r["portability_class"] == "B":
            d["B"] += 1
    return {k for k, d in agg.items()
            if len(d["sp"]) >= 3 and len(d["bp"]) >= 10 and d["n"] >= 20 and d["B"] > 0}


def build(rows, elig, mode, keep=lambda r: True):
    """2x2 tables per family. mode 'two_way' = B vs C+D+E; 'within' = B vs A."""
    tab = {}
    for r in rows:
        h = host(r["organism_harmonized"])
        if h is None or r["gene_family"] not in elig or not keep(r):
            continue
        c = r["portability_class"]
        if mode == "two_way":
            if c == "A":
                continue
            idx = 0 if c == "B" else 1
        else:
            if c not in ("A", "B"):
                continue
            idx = 0 if c == "B" else 1
        tab.setdefault(r["gene_family"], {"AB": [0, 0], "KL": [0, 0]})[h][idx] += 1
    out, names = [], []
    for k, c in sorted(tab.items()):
        a, b = c["AB"]
        cc, d = c["KL"]
        if (a + b) == 0 or (cc + d) == 0:
            continue
        out.append((a, b, cc, d))
        names.append(k)
    return out, names


def summarise(tables, names, label):
    m = mh(tables)
    if m is None:
        return {"arm": label, "estimable": False,
                "why": "no family has both a case and a comparator in both hosts"}
    return {"arm": label, "estimable": True, "n_families": len(names),
            "or": round(m["or"], 4), "ci_lo": round(m["ci_lo"], 4),
            "ci_hi": round(m["ci_hi"], 4),
            "excludes_1": bool(m["ci_lo"] > 1 or m["ci_hi"] < 1),
            "direction": "A. baumannii" if m["or"] > 1 else "Klebsiella",
            "concordant": sum(1 for t in tables if woolf(*t)[0] > 0)}


def cluster_boot(rows, elig, mode, unit_of, label):
    """Percentile interval from resampling clusters with replacement."""
    by = {}
    for r in rows:
        if host(r["organism_harmonized"]) is None:
            continue
        by.setdefault(unit_of(r), []).append(r)
    units = sorted(by)
    rng = random.Random(SEED)
    est = []
    for _ in range(B_BOOT):
        draw = [by[units[rng.randrange(len(units))]] for _ in range(len(units))]
        flat = [r for g in draw for r in g]
        t, n = build(flat, elig, mode)
        m = mh(t)
        if m:
            est.append(m["or"])
    est.sort()
    if len(est) < 100:
        return {"arm": label, "replicates": len(est), "estimable": False}
    lo = est[int(0.025 * len(est))]
    hi = est[int(0.975 * len(est))]
    return {"arm": label, "clusters": len(units), "replicates": len(est),
            "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
            "median": round(est[len(est) // 2], 4),
            "excludes_1": bool(lo > 1 or hi < 1)}


def main():
    doc = json.load(io.open(AMEND, encoding="utf-8"))
    claimed = doc.pop("sha256_of_body")
    actual = hashlib.sha256(json.dumps(doc, indent=1, ensure_ascii=False,
                                       sort_keys=True).encode("utf-8")).hexdigest()
    if actual != claimed:
        raise SystemExit("REFUSING: amendment body changed since freeze")
    print("frozen design verified: body %s\n" % claimed[:16])

    rows = list(csv.DictReader(io.open(TABLE, encoding="utf-8", newline=""),
                              delimiter="\t"))
    st = {r["assembly_version"]: r["st"].strip()
          for r in csv.DictReader(io.open(MLST, encoding="utf-8", newline=""),
                                  delimiter="\t")}
    for r in rows:
        if host(r["organism_harmonized"]):
            r["_st"] = st.get(r["assembly_version"], "?")
            r["_lin"] = "%s|%s" % (host(r["organism_harmonized"]),
                                   r["_st"] if r["_st"] != "-"
                                   else "UNTYPED:" + r["assembly_version"])
    elig = eligible(rows)

    base = summarise(*build(rows, elig, "two_way"), "baseline two-way")
    print("baseline reproduction: %d families, OR %.4f" % (base["n_families"], base["or"]))
    if base["n_families"] != BASE_N or abs(base["or"] - BASE_OR) > 0.001:
        raise SystemExit("baseline did not reproduce - no arm reported")

    chosen = {}
    for r in rows:
        if host(r["organism_harmonized"]) is None:
            continue
        d = hashlib.sha256(r["assembly_version"].encode()).hexdigest()
        if r["_lin"] not in chosen or d < chosen[r["_lin"]][0]:
            chosen[r["_lin"]] = (d, r["assembly_version"])
    one_per_st = {a for _d, a in chosen.values()}

    def not_st2(r):
        return not (host(r["organism_harmonized"]) == "AB" and r["_st"] == "2")

    W1 = summarise(*build(rows, elig, "within"), "W1 within-chromosome, B vs A")
    W2 = summarise(*build(rows, elig, "within",
                          keep=lambda r: r["assembly_version"] in one_per_st),
                   "W2 within-chromosome, one genome per ST")
    N1 = summarise(*build(rows, elig, "two_way", keep=not_st2),
                   "N1 two-way, non-ST2 A. baumannii")
    N2 = summarise(*build(rows, elig, "within", keep=not_st2),
                   "N2 within-chromosome, non-ST2 A. baumannii")

    print("\n%-44s %8s %8s %-18s %s" % ("arm", "families", "OR", "95% CI", "concordant"))
    for m in (base, W1, W2, N1, N2):
        if not m["estimable"]:
            print("%-44s  NOT ESTIMABLE: %s" % (m["arm"], m["why"]))
            continue
        print("%-44s %8d %8.2f %-18s %d"
              % (m["arm"], m["n_families"], m["or"],
                 "%.2f-%.2f" % (m["ci_lo"], m["ci_hi"]), m["concordant"]))

    print("\ncluster bootstrap on the headline two-way contrast, B = %d, seed %d"
          % (B_BOOT, SEED))
    C1 = cluster_boot(rows, elig, "two_way", lambda r: r["bioproject_accession"],
                      "C1 BioProject-clustered")
    C2 = cluster_boot(rows, elig, "two_way", lambda r: r["_lin"],
                      "C2 lineage-clustered")
    for c in (C1, C2):
        print("  %-26s %5d clusters | 95%% CI %.2f-%.2f | median %.2f"
              % (c["arm"], c["clusters"], c["ci_lo"], c["ci_hi"], c["median"]))
    print("  analytic RBG interval for comparison: %.2f-%.2f"
          % (base["ci_lo"], base["ci_hi"]))

    dr = doc["decision_rule"]
    half = W1["estimable"] and W1["or"] < BASE_OR / 2
    reframe = any(m["estimable"] and (not m["excludes_1"]
                                      or m["direction"] != "A. baumannii")
                  for m in (N1, N2))
    widest = max((base["ci_hi"] - base["ci_lo"],
                  C1["ci_hi"] - C1["ci_lo"], C2["ci_hi"] - C2["ci_lo"]))
    print("\nfrozen decision rule:")
    print("  abstract must change (W1 below half of 50.29): %s" % half)
    print("  finding reframed as ST2/GC2-specific            : %s" % reframe)
    print("  widest interval is the cluster bootstrap        : %s"
          % (widest > base["ci_hi"] - base["ci_lo"]))

    rec = {"receipt": "NM-V4C-004 referee-requested analyses",
           "amendment_sha256_of_body": claimed,
           "amendment_sha256_of_file": sha(AMEND),
           "inputs": {"occurrence_table": sha(TABLE), "mlst_calls": sha(MLST)},
           "baseline_reproduced": True,
           "results": {"baseline_two_way": base, "W1": W1, "W2": W2,
                       "N1": N1, "N2": N2, "C1": C1, "C2": C2},
           "decision_rule": dr,
           "abstract_must_change": bool(half),
           "reframe_as_ST2_specific": bool(reframe)}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("\nreceipt: %s" % os.path.basename(OUT))


if __name__ == "__main__":
    main()
