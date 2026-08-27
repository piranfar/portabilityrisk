"""Diagnostics for the within-chromosome arm W1, and the pre-registered gate it fails.

W1 came back at 0.485 - the direction reversed. Reporting that as "Klebsiella
chromosomal determinants are more IS-adjacent" would be wrong, and checking why
before reporting is not rescuing the result: the two checks applied here were
both registered before this arm existed.

  * G2, from the ORIGINAL frozen design (NMV4C_FROZEN_DESIGN.json): "the single
    largest family carries no more than 30 per cent of the MH weight". W1 puts
    63.9% of its weight on blaOXA, so W1 fails a gate this project registered
    long before the referee asked for the arm.
  * NM-V4C-001 already registered the exclusion of the intrinsic Acinetobacter
    OXA nodes. Applying an existing registered exclusion to a new endpoint is
    not a new analysis choice.

The other thing W1 hides is sparsity: a family where BOTH hosts have zero class-A
occurrences contributes exactly zero to both MH sums and vanishes without
comment. Counting those is the difference between "49 families" and what the
estimate actually rests on.
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
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = "" + _dir("PORTABILITYRISK_REPO_DIR") + ""
NM = REPO + "docs/nature_microbiology/"
TABLE = ("" + _dir("PORTABILITYRISK_DEPOSIT_DIR") + ""
         "portabilityrisk_occurrence_portability_v1.tsv")
MLST = NM + "NM_V4C_MLST_CALLS_V1.tsv"
AMRPROT = REPO + ".cache/model2_g_a/amrfinder_db/2026-08-07.1/AMRProt.fa"
REC = NM + "NM_V4C_REFEREE_ANALYSES_RESULTS_V1.json"

sys.path.insert(0, REPO + "audit/ingest/assay_aware_emergence/v2/nm_validation")
from nmv4c_score import mh, woolf                                  # noqa: E402

INTRINSIC = {"blaOXA-51_fam", "blaOXA-134_fam", "blaOXA-213_fam"}
G2_CEILING = 0.30


def host(o):
    if o == "Acinetobacter baumannii":
        return "AB"
    if o.startswith("Klebsiella"):
        return "KL"
    return None


def main():
    rows = list(csv.DictReader(io.open(TABLE, encoding="utf-8", newline=""),
                               delimiter="\t"))
    node = {}
    for line in io.open(AMRPROT, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            p = line[1:].split("|")
            if len(p) > 5 and p[3]:
                node.setdefault(p[3], p[4])
    st = {r["assembly_version"]: r["st"].strip()
          for r in csv.DictReader(io.open(MLST, encoding="utf-8", newline=""),
                                  delimiter="\t")}

    fam = collections.defaultdict(lambda: {"sp": set(), "bp": set(), "n": 0, "B": 0})
    for r in rows:
        d = fam[r["gene_family"]]
        d["sp"].add(r["organism_harmonized"])
        d["bp"].add(r["bioproject_accession"])
        d["n"] += 1
        if r["portability_class"] == "B":
            d["B"] += 1
    elig = {k for k, d in fam.items()
            if len(d["sp"]) >= 3 and len(d["bp"]) >= 10 and d["n"] >= 20 and d["B"] > 0}

    def arm(keep, label):
        tab = collections.defaultdict(lambda: {"AB": [0, 0], "KL": [0, 0]})
        for r in rows:
            h = host(r["organism_harmonized"])
            c = r["portability_class"]
            if h and c in ("A", "B") and r["gene_family"] in elig and keep(r):
                tab[r["gene_family"]][h][0 if c == "B" else 1] += 1
        tables, names, wts = [], [], []
        zero = 0
        for f, v in sorted(tab.items()):
            a, b = v["AB"]
            c2, d = v["KL"]
            n = a + b + c2 + d
            if n == 0 or (a + b) == 0 or (c2 + d) == 0:
                continue
            w = a * d / n + b * c2 / n
            if w == 0:
                zero += 1
            tables.append((a, b, c2, d))
            names.append(f)
            wts.append(w)
        m = mh(tables)
        tw = sum(wts)
        order = sorted(zip(names, wts), key=lambda kv: -kv[1])
        top_f, top_w = order[0] if order else ("-", 0.0)
        share = top_w / tw if tw else 0.0
        out = {"arm": label, "families_in_table": len(names),
               "families_carrying_zero_weight": zero,
               "families_actually_informative": len(names) - zero,
               "largest_weight_family": top_f,
               "largest_weight_share": round(share, 4),
               "passes_registered_G2_ceiling": bool(share <= G2_CEILING),
               "or": round(m["or"], 4) if m else None,
               "ci_lo": round(m["ci_lo"], 4) if m else None,
               "ci_hi": round(m["ci_hi"], 4) if m else None,
               "direction": ("A. baumannii" if m and m["or"] > 1 else "Klebsiella")
                            if m else None,
               "concordant": sum(1 for t in tables if woolf(*t)[0] > 0)}
        return out

    def not_st2(r):
        return not (host(r["organism_harmonized"]) == "AB"
                    and st.get(r["assembly_version"], "?") == "2")

    arms = [
        arm(lambda r: True, "W1 as the referee specified it"),
        arm(lambda r: node.get(r["determinant_name"]) not in INTRINSIC,
            "W1 minus intrinsic Acinetobacter OXA (NM-V4C-001 exclusion)"),
        arm(lambda r: r["gene_family"] != "blaOXA",
            "W1 minus the whole blaOXA family"),
        arm(lambda r: not_st2(r) and node.get(r["determinant_name"]) not in INTRINSIC,
            "W1 minus intrinsic OXA and minus ST2"),
    ]

    print("%-58s %7s %-16s %6s %s" % ("arm", "OR", "95% CI", "G2", "informative families"))
    for a in arms:
        print("%-58s %7.2f %-16s %-6s %d of %d"
              % (a["arm"], a["or"], "%.2f-%.2f" % (a["ci_lo"], a["ci_hi"]),
                 "pass" if a["passes_registered_G2_ceiling"] else "FAIL",
                 a["families_actually_informative"], a["families_in_table"]))
        print("%-58s        largest family %s at %.1f%% of the weight"
              % ("", a["largest_weight_family"], 100 * a["largest_weight_share"]))

    r = json.load(io.open(REC, encoding="utf-8"))
    r["W1_diagnostics"] = {
        "why": ("W1 reversed direction. Before reporting that, two checks "
                "registered BEFORE this arm existed were applied: the G2 weight "
                "ceiling from the original frozen design, and the intrinsic-OXA "
                "exclusion from NM-V4C-001."),
        "G2_ceiling": G2_CEILING,
        "G2_source": "NMV4C_FROZEN_DESIGN.json, gates.G2_multiple_families",
        "sparsity_note": ("a family in which BOTH hosts have zero class-A "
                          "occurrences contributes zero to both Mantel-Haenszel "
                          "sums and disappears without comment; those are counted "
                          "here rather than left inside the family total"),
        "arms": arms}
    io.open(REC, "w", encoding="utf-8", newline="\n").write(
        json.dumps(r, indent=1, ensure_ascii=False) + "\n")
    print("\ndiagnostics appended to %s" % REC.split("/")[-1])


if __name__ == "__main__":
    main()
