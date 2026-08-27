"""Run the frozen lineage-adjustment arms (NM-V4C-002).

Refuses to run unless the amendment's body hash still matches, and refuses to
report any arm unless the published baseline reproduces first.

The estimator is imported from nmv4c_score, the module that produced the baseline,
so a change between baseline and arm cannot come from a change in the estimator.
Only which occurrences enter, and with what weight, differs.
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
AMEND = NM + "NM_V4C_LINEAGE_ADJUSTMENT_AMENDMENT_002.json"
MLST = NM + "NM_V4C_MLST_CALLS_V1.tsv"
TABLE = (_dir("PORTABILITYRISK_DEPOSIT_DIR")
         "portabilityrisk_occurrence_portability_v1.tsv")
OUT = NM + "NM_V4C_LINEAGE_ADJUSTMENT_RESULTS_V1.json"

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
    actual = hashlib.sha256(json.dumps(doc, indent=1, ensure_ascii=False,
                                       sort_keys=True).encode("utf-8")).hexdigest()
    if actual != claimed:
        raise SystemExit("REFUSING: the amendment body changed since it was frozen")
    doc["sha256_of_body"] = claimed
    print("frozen design verified: body %s" % claimed[:16])
    return doc


def host(o):
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


def pool(rows, elig, keep=lambda r: True, weight=lambda r: 1.0):
    """Matched-family Mantel-Haenszel over the AB/KL contrast."""
    tab = {}
    for r in rows:
        h = host(r["organism_harmonized"])
        if h is None or r["portability_class"] == "A":
            continue
        if r["gene_family"] not in elig or not keep(r):
            continue
        c = tab.setdefault(r["gene_family"], {"AB": [0.0, 0.0], "KL": [0.0, 0.0]})
        c[h][0 if r["portability_class"] == "B" else 1] += weight(r)
    tables, conc, fams = [], 0, []
    for f, c in sorted(tab.items()):
        a, b = c["AB"]
        cc, d = c["KL"]
        if (a + b) <= 0 or (cc + d) <= 0:
            continue
        tables.append((a, b, cc, d))
        fams.append(f)
        if woolf(a, b, cc, d)[0] > 0:
            conc += 1
    m = mh(tables)
    if m is None:
        return None
    return {"n_families": len(fams), "or": round(m["or"], 4), "ln_or": m["ln_or"],
            "ci_lo": round(m["ci_lo"], 4), "ci_hi": round(m["ci_hi"], 4),
            "excludes_1": bool(m["ci_lo"] > 1 or m["ci_hi"] < 1),
            "direction": "A. baumannii" if m["or"] > 1 else "Klebsiella",
            "families_concordant": conc}


def main():
    doc = load_amendment()
    rows = list(csv.DictReader(io.open(TABLE, encoding="utf-8", newline=""),
                              delimiter="\t"))
    st = {}
    for r in csv.DictReader(io.open(MLST, encoding="utf-8", newline=""),
                            delimiter="\t"):
        st[r["assembly_version"]] = r["st"].strip()

    # attach a lineage id to every occurrence of the two hosts
    miss = set()
    for r in rows:
        h = host(r["organism_harmonized"])
        if h is None:
            continue
        s = st.get(r["assembly_version"])
        if s is None:
            miss.add(r["assembly_version"])
            s = "?"
        r["_st"] = s
        r["_lin"] = "%s|%s" % (h, s if s != "-" else "UNTYPED:" + r["assembly_version"])
        r["_linpool"] = "%s|%s" % (h, s)          # untypeable pooled into one lineage
    if miss:
        raise SystemExit("REFUSING: %d genomes have no MLST call: %s"
                         % (len(miss), sorted(miss)[:5]))

    elig = eligible(rows)
    base = pool(rows, elig)
    print("\nbaseline reproduction: %d families | OR %.4f | CI %.2f-%.2f"
          % (base["n_families"], base["or"], base["ci_lo"], base["ci_hi"]))
    if base["n_families"] != BASE_FAMILIES or abs(base["or"] - BASE_OR) > 0.001:
        raise SystemExit("baseline did not reproduce - no arm reported")

    # ---- descriptive composition, required by the amendment --------------
    desc = {}
    for h, name in (("AB", "A. baumannii"), ("KL", "Klebsiella group")):
        gen = {r["assembly_version"]: r["_st"] for r in rows
               if host(r["organism_harmonized"]) == h}
        cnt = {}
        for s in gen.values():
            cnt[s] = cnt.get(s, 0) + 1
        n = len(gen)
        typed = {k: v for k, v in cnt.items() if k != "-"}
        tot_typed = sum(typed.values())
        hhi = sum((v / tot_typed) ** 2 for v in typed.values()) if tot_typed else 0
        top = max(typed.items(), key=lambda kv: kv[1]) if typed else ("-", 0)
        desc[name] = {"genomes": n, "distinct_STs": len(typed),
                      "untypeable_genomes": cnt.get("-", 0),
                      "effective_STs_1_over_HHI": round(1 / hhi, 2) if hhi else None,
                      "largest_ST": top[0], "largest_ST_genomes": top[1],
                      "largest_ST_share": round(top[1] / n, 4)}
        d = desc[name]
        print("  %-18s %4d genomes | %3d STs | %2d untypeable | effective %.2f | "
              "largest ST %s = %.1f%%"
              % (name, d["genomes"], d["distinct_STs"], d["untypeable_genomes"],
                 d["effective_STs_1_over_HHI"], d["largest_ST"],
                 100 * d["largest_ST_share"]))

    # ---- L1 / L1u : one genome per sequence type -------------------------
    def one_per(field, drop_untyped):
        chosen = {}
        for r in rows:
            if host(r["organism_harmonized"]) is None:
                continue
            if drop_untyped and r["_st"] == "-":
                continue
            k = r[field]
            a = r["assembly_version"]
            dig = hashlib.sha256(a.encode()).hexdigest()
            if k not in chosen or dig < chosen[k][0]:
                chosen[k] = (dig, a)
        return {a for _d, a in chosen.values()}

    keepL1 = one_per("_lin", False)
    keepL1u = one_per("_lin", True)
    L1 = pool(rows, elig, keep=lambda r: r["assembly_version"] in keepL1)
    L1u = pool(rows, elig, keep=lambda r: r["assembly_version"] in keepL1u)

    # ---- L3 : lineage-balanced weighting ---------------------------------
    gsize = {}
    for r in rows:
        if host(r["organism_harmonized"]) is None:
            continue
        gsize.setdefault(r["_lin"], set()).add(r["assembly_version"])
    L3 = pool(rows, elig, weight=lambda r: 1.0 / len(gsize[r["_lin"]]))

    # ---- L2 / L2u : leave-one-lineage-out --------------------------------
    def loo(field):
        worst = (None, 0.0)
        lins = sorted({r[field] for r in rows
                       if host(r["organism_harmonized"]) is not None})
        for lin in lins:
            m = pool(rows, elig, keep=lambda r, L=lin: r[field] != L)
            if m is None:
                continue
            rel = abs(m["ln_or"] - base["ln_or"]) / abs(base["ln_or"])
            if rel > worst[1]:
                worst = (lin, rel)
        return len(lins), worst

    n_l2, w_l2 = loo("_lin")
    n_l2u, w_l2u = loo("_linpool")

    print("\n%-42s %8s %7s %-16s %s" % ("arm", "families", "OR", "95% CI", "concordant"))
    for label, m in (("baseline, as published", base),
                     ("L1  one genome per ST", L1),
                     ("L1u L1, untypeable dropped", L1u),
                     ("L3  lineage-balanced weighting", L3)):
        print("%-42s %8d %7.2f %-16s %d"
              % (label, m["n_families"], m["or"],
                 "%.2f-%.2f" % (m["ci_lo"], m["ci_hi"]), m["families_concordant"]))
    print("\nL2  leave-one-ST-out over %d lineages : max relative change in ln(OR) "
          "%.4f  (ceiling %.2f) on %s" % (n_l2, w_l2[1], 0.15, w_l2[0]))
    print("L2u untypeable pooled, %d lineages    : max relative change %.4f on %s"
          % (n_l2u, w_l2u[1], w_l2u[0]))

    ceiling = doc["arms"]["L2"]["ceiling"]
    arms_ok = all(m["excludes_1"] and m["direction"] == "A. baumannii"
                  for m in (L1, L1u, L3))
    loo_ok = w_l2[1] <= ceiling and w_l2u[1] <= ceiling
    verdict = ("SURVIVES_FULLY" if arms_ok and loo_ok else
               "SURVIVES_PARTLY" if arms_ok else "FAILS")
    print("\nfrozen decision rule -> %s" % verdict)
    print({"SURVIVES_FULLY": "  every arm keeps direction and excludes 1, and "
                             "leave-one-ST-out stays within the ceiling. The title "
                             "stands.",
           "SURVIVES_PARTLY": "  arms hold but leave-one-ST-out exceeds the ceiling: "
                              "the effect is real and lineage-sensitive. The title "
                              "must become 'host-associated'.",
           "FAILS": "  an arm lost direction or covers 1. The host claim is "
                    "withdrawn and rewritten."}[verdict])

    rec = {"receipt": "NM-V4C-002 lineage adjustment by MLST",
           "amendment_sha256_of_body": doc["sha256_of_body"],
           "amendment_sha256_of_file": sha_file(AMEND),
           "inputs": {"occurrence_table": sha_file(TABLE),
                      "mlst_calls": sha_file(MLST),
                      "mlst_version": "2.35.0",
                      "schemes": doc["typing"]["scheme_forced_per_host"]},
           "baseline_reproduced": True,
           "composition": desc,
           "results": {"baseline": base, "L1": L1, "L1u": L1u, "L3": L3,
                       "L2": {"n_lineages": n_l2, "max_relative_ln_or_change": w_l2[1],
                              "worst_lineage": w_l2[0], "ceiling": ceiling,
                              "within_ceiling": bool(w_l2[1] <= ceiling)},
                       "L2u": {"n_lineages": n_l2u, "max_relative_ln_or_change": w_l2u[1],
                               "worst_lineage": w_l2u[0], "ceiling": ceiling,
                               "within_ceiling": bool(w_l2u[1] <= ceiling)}},
           "decision_rule": doc["decision_rule"],
           "verdict": verdict}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("\nreceipt: %s" % os.path.basename(OUT))


if __name__ == "__main__":
    main()
