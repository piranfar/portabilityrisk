"""Run the frozen marginal within-chromosome adjacency arms (NM-V4C-005)."""
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
import random
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = _dir("PORTABILITYRISK_REPO_DIR")
NM = REPO + "docs/nature_microbiology/"
AMEND = NM + "NM_V4C_MARGINAL_ADJACENCY_AMENDMENT_005.json"
MLST = NM + "NM_V4C_MLST_CALLS_V1.tsv"
TABLE = (_dir("PORTABILITYRISK_DEPOSIT_DIR")
         "portabilityrisk_occurrence_portability_v1.tsv")
OUT = NM + "NM_V4C_MARGINAL_ADJACENCY_RESULTS_V1.json"

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


def rate(rows):
    """class B / (class A + class B), per host."""
    c = collections.Counter()
    for r in rows:
        h = r["_h"]
        if r["portability_class"] in ("A", "B"):
            c[(h, r["portability_class"])] += 1
    out = {}
    for h in ("AB", "KL"):
        a, b = c[(h, "A")], c[(h, "B")]
        out[h] = (b / (a + b)) if (a + b) else None
        out[h + "_n"] = a + b
        out[h + "_B"] = b
    return out


def boot(rows, unit_of, keep=lambda r: True):
    by = collections.defaultdict(list)
    for r in rows:
        if keep(r):
            by[unit_of(r)].append(r)
    units = sorted(by)
    rng = random.Random(SEED)
    ab, kl, diff, ratio = [], [], [], []
    for _ in range(B_BOOT):
        flat = []
        for _i in range(len(units)):
            flat.extend(by[units[rng.randrange(len(units))]])
        s = rate(flat)
        if s["AB"] is None or s["KL"] is None:
            continue
        ab.append(s["AB"])
        kl.append(s["KL"])
        diff.append(s["AB"] - s["KL"])
        ratio.append(s["AB"] / s["KL"] if s["KL"] else float("nan"))

    def ci(v):
        v = sorted(x for x in v if x == x)
        return [round(v[int(0.025 * len(v))], 4), round(v[int(0.975 * len(v))], 4)]
    return {"clusters": len(units), "replicates": len(diff),
            "AB_ci": ci(ab), "KL_ci": ci(kl), "difference_ci": ci(diff),
            "ratio_ci": ci(ratio)}


def standardise(rows, target, source):
    """Apply `source`'s family-specific rates to `target`'s family distribution."""
    num = collections.defaultdict(lambda: [0, 0])   # family -> [B, A+B] for source
    dist = collections.Counter()                    # family -> A+B for target
    for r in rows:
        if r["portability_class"] not in ("A", "B"):
            continue
        f = r["gene_family"]
        if r["_h"] == source:
            num[f][1] += 1
            if r["portability_class"] == "B":
                num[f][0] += 1
        if r["_h"] == target:
            dist[f] += 1
    shared = [f for f in dist if num.get(f, [0, 0])[1] > 0]
    w = sum(dist[f] for f in shared)
    if not w:
        return None
    std = sum(dist[f] * (num[f][0] / num[f][1]) for f in shared) / w
    return {"standardised_rate": round(std, 4), "families_shared": len(shared),
            "target_occurrences_covered": w,
            "target_coverage": round(w / sum(dist.values()), 4)}


def main():
    doc = json.load(io.open(AMEND, encoding="utf-8"))
    claimed = doc.pop("sha256_of_body")
    actual = hashlib.sha256(json.dumps(doc, indent=1, ensure_ascii=False,
                                       sort_keys=True).encode("utf-8")).hexdigest()
    if actual != claimed:
        raise SystemExit("REFUSING: amendment body changed since freeze")
    print("frozen design verified: body %s\n" % claimed[:16])

    allrows = list(csv.DictReader(io.open(TABLE, encoding="utf-8", newline=""),
                                 delimiter="\t"))
    st = {r["assembly_version"]: r["st"].strip()
          for r in csv.DictReader(io.open(MLST, encoding="utf-8", newline=""),
                                  delimiter="\t")}
    rows = []
    for r in allrows:
        h = host(r["organism_harmonized"])
        if h is None:
            continue
        r["_h"] = h
        s = st.get(r["assembly_version"], "?")
        r["_lin"] = "%s|%s" % (h, s if s != "-" else "UNTYPED:" + r["assembly_version"])
        rows.append(r)

    # ---- M1 crude -------------------------------------------------------
    m1 = rate(rows)
    print("M1  crude within-chromosome marker adjacency")
    for h, name in (("AB", "A. baumannii"), ("KL", "Klebsiella group")):
        print("      %-18s %6.1f%%   (%s of %s chromosomal occurrences)"
              % (name, 100 * m1[h], "{:,}".format(m1[h + "_B"]),
                 "{:,}".format(m1[h + "_n"])))
    print("      difference %+.1f pp | ratio %.2f"
          % (100 * (m1["AB"] - m1["KL"]), m1["AB"] / m1["KL"]))

    bp = boot(rows, lambda r: r["bioproject_accession"])
    ln = boot(rows, lambda r: r["_lin"])
    print("\n      BioProject-clustered  difference 95%% CI %+.4f to %+.4f  (%d clusters)"
          % (bp["difference_ci"][0], bp["difference_ci"][1], bp["clusters"]))
    print("      lineage-clustered     difference 95%% CI %+.4f to %+.4f  (%d clusters)"
          % (ln["difference_ci"][0], ln["difference_ci"][1], ln["clusters"]))

    # ---- M2 one genome per sequence type --------------------------------
    chosen = {}
    for r in rows:
        d = hashlib.sha256(r["assembly_version"].encode()).hexdigest()
        if r["_lin"] not in chosen or d < chosen[r["_lin"]][0]:
            chosen[r["_lin"]] = (d, r["assembly_version"])
    keep1 = {a for _d, a in chosen.values()}
    sub = [r for r in rows if r["assembly_version"] in keep1]
    m2 = rate(sub)
    m2boot = boot(sub, lambda r: r["_lin"])
    print("\nM2  one genome per sequence type")
    print("      A. baumannii %6.1f%%  |  Klebsiella %6.1f%%  |  difference %+.1f pp"
          % (100 * m2["AB"], 100 * m2["KL"], 100 * (m2["AB"] - m2["KL"])))
    print("      lineage-clustered difference 95%% CI %+.4f to %+.4f"
          % (m2boot["difference_ci"][0], m2boot["difference_ci"][1]))

    # ---- M3 direct standardisation --------------------------------------
    ab_on_kl = standardise(rows, target="KL", source="AB")
    kl_on_ab = standardise(rows, target="AB", source="KL")
    print("\nM3  direct standardisation, to measure how much of the gap is composition")
    print("      A. baumannii rates applied to Klebsiella's family mix : %6.1f%%  "
          "(covers %.0f%% of Klebsiella occurrences)"
          % (100 * ab_on_kl["standardised_rate"], 100 * ab_on_kl["target_coverage"]))
    print("      Klebsiella rates applied to A. baumannii's family mix : %6.1f%%  "
          "(covers %.0f%% of A. baumannii occurrences)"
          % (100 * kl_on_ab["standardised_rate"], 100 * kl_on_ab["target_coverage"]))
    crude_gap = m1["AB"] - m1["KL"]
    resid_gap = m1["AB"] - kl_on_ab["standardised_rate"]
    removed = 1 - resid_gap / crude_gap
    print("      crude gap %+.1f pp | gap remaining after standardising to "
          "A. baumannii's mix %+.1f pp" % (100 * crude_gap, 100 * resid_gap))
    print("      composition explains %.0f%% of the crude gap" % (100 * removed))

    compositional = removed > 0.5
    holds = (bp["difference_ci"][0] > 0 and ln["difference_ci"][0] > 0
             and m2boot["difference_ci"][0] > 0)
    print("\nfrozen decision rule:")
    print("  claim supported (difference excludes zero under every clustering): %s" % holds)
    print("  must be described as substantially compositional (>50%% explained) : %s"
          % compositional)

    rec = {"receipt": "NM-V4C-005 within-chromosome marker adjacency as a rate",
           "amendment_sha256_of_body": claimed,
           "amendment_sha256_of_file": sha(AMEND),
           "inputs": {"occurrence_table": sha(TABLE), "mlst_calls": sha(MLST)},
           "post_hoc": True,
           "M1_crude": {h: (round(m1[h], 4), m1[h + "_B"], m1[h + "_n"])
                        for h in ("AB", "KL")},
           "M1_bootstrap_bioproject": bp, "M1_bootstrap_lineage": ln,
           "M2_one_per_ST": {"AB": round(m2["AB"], 4), "KL": round(m2["KL"], 4),
                             "difference": round(m2["AB"] - m2["KL"], 4),
                             "bootstrap": m2boot},
           "M3_standardised": {"AB_rates_on_KL_mix": ab_on_kl,
                               "KL_rates_on_AB_mix": kl_on_ab,
                               "crude_gap": round(crude_gap, 4),
                               "gap_after_standardisation": round(resid_gap, 4),
                               "share_of_gap_explained_by_composition": round(removed, 4)},
           "decision_rule": doc["decision_rule"],
           "claim_supported": bool(holds),
           "substantially_compositional": bool(compositional)}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("\nreceipt: %s" % os.path.basename(OUT))


if __name__ == "__main__":
    main()
