"""Run the frozen composite-element stratification (NM-C1-006)."""
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
AMEND = NM + "NM_C1_COMPOSITE_ELEMENT_AMENDMENT_006.json"
MATS = _dir("PORTABILITYRISK_C1_DIR") + "c1_4_matrices.npz"
ISTAB = (_dir("PORTABILITYRISK_DEPOSIT_DIR")
         "genome_wide_is_elements.tsv")
OUT = NM + "NM_C1_COMPOSITE_ELEMENT_RESULTS_V1.json"

GROUPS = {"Acinetobacter baumannii": "A. baumannii",
          "Pseudomonas aeruginosa": "P. aeruginosa"}
PUBLISHED_AB = 16.91


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def group_of(sp):
    if sp in GROUPS:
        return GROUPS[sp]
    if sp.startswith("Klebsiella"):
        return "Klebsiella group"
    if sp.startswith("Enterobacter"):
        return "Enterobacter group"
    return None


def complete_structural(r):
    """The primary analysis definition, applied to a census row."""
    try:
        if (r.get("type") or "").strip() != "c":
            return False
        if int(float(r.get("orfLen") or 0)) <= 0:
            return False
        if int(float(r.get("irLen") or 0)) <= 0:
            return False
        for k in ("start1", "end1", "start2", "end2"):
            if not (r.get(k) or "").strip():
                return False
        return True
    except (TypeError, ValueError):
        return False


def main():
    doc = json.load(io.open(AMEND, encoding="utf-8"))
    claimed = doc.pop("sha256_of_body")
    actual = hashlib.sha256(json.dumps(doc, indent=1, ensure_ascii=False,
                                       sort_keys=True).encode("utf-8")).hexdigest()
    if actual != claimed:
        raise SystemExit("REFUSING: amendment body changed since freeze")
    print("frozen design verified: body %s\n" % claimed[:16])

    z = np.load(MATS, allow_pickle=False)
    ids, acc, species = z["ids"], z["acc"], z["species"]
    obs, null = z["obs"].astype(np.int32), z["null"].astype(np.int32)
    n = len(ids)
    print("occurrences in the null matrices: %s | permutations: %d" % ("{:,}".format(n), null.shape[0]))

    starts = np.zeros(n, dtype=np.int64)
    ends = np.zeros(n, dtype=np.int64)
    for i, s in enumerate(ids):
        p = s.split("|")
        starts[i] = int(p[-2])
        ends[i] = int(p[-1])

    # ---- complete-structural elements, indexed by replicon and family ----
    by = collections.defaultdict(list)
    kept = 0
    with io.open(ISTAB, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not complete_structural(r):
                continue
            try:
                b, e = int(float(r["isBegin"])), int(float(r["isEnd"]))
            except (TypeError, ValueError):
                continue
            by[(r["replicon_accession"], r["family"])].append((b, e))
            kept += 1
    print("complete-structural elements indexed: %s" % "{:,}".format(kept))

    fams = collections.defaultdict(set)
    for (rep, fam) in by:
        fams[rep].add(fam)

    def flanked_mask(D):
        m = np.zeros(n, dtype=bool)
        for i in range(n):
            rep = acc[i]
            gs, ge = starts[i], ends[i]
            for fam in fams.get(rep, ()):
                up = dn = False
                for b, e in by[(rep, fam)]:
                    if e < gs and gs - e <= D:
                        up = True
                    elif b > ge and b - ge <= D:
                        dn = True
                    if up and dn:
                        break
                if up and dn:
                    m[i] = True
                    break
        return m

    def enrich(mask, d=1000):
        if mask.sum() == 0:
            return None
        o = float(((obs[mask] >= 0) & (obs[mask] <= d)).mean())
        e = float(((null[:, mask] >= 0) & (null[:, mask] <= d)).mean())
        return {"n": int(mask.sum()), "observed": round(o, 4),
                "expected": round(e, 4),
                "enrichment": round(o / e, 2) if e > 0 else None}

    results = {}
    for D in (doc["definition"]["D_primary"], doc["definition"]["D_secondary"]):
        fl = flanked_mask(D)
        print("\nD = %s bp : %s of %s chromosomal occurrences are composite-flanked (%.1f%%)"
              % ("{:,}".format(D), "{:,}".format(int(fl.sum())), "{:,}".format(n),
                 100 * fl.mean()))
        print("  %-20s %9s %10s %10s %12s"
              % ("group", "n", "observed", "expected", "enrichment"))
        per = {}
        for g in ("A. baumannii", "Klebsiella group", "P. aeruginosa",
                  "Enterobacter group"):
            gm = np.array([group_of(s) == g for s in species])
            for lab, m in (("flanked", gm & fl), ("non-flanked", gm & ~fl)):
                r = enrich(m)
                per.setdefault(g, {})[lab] = r
                if r:
                    print("  %-20s %9s %10.4f %10.4f %12.2f"
                          % ("%s, %s" % (g[:11], lab), "{:,}".format(r["n"]),
                             r["observed"], r["expected"], r["enrichment"]))
        results["D_%d" % D] = per

    prim = results["D_%d" % doc["definition"]["D_primary"]]
    ab_nf = prim["A. baumannii"]["non-flanked"]["enrichment"]
    kl_nf = prim["Klebsiella group"]["non-flanked"]["enrichment"]
    pa_nf = prim["P. aeruginosa"]["non-flanked"]["enrichment"]
    ratio_kl = ab_nf / kl_nf if kl_nf else None
    ratio_pa = ab_nf / pa_nf if pa_nf else None

    print("\nnon-flanked stratum, at D = 10 kb:")
    print("  A. baumannii enrichment        %.2f   (published, all occurrences: %.2f)"
          % (ab_nf, PUBLISHED_AB))
    print("  A. baumannii / Klebsiella      %.2f   (published: 2.69)" % ratio_kl)
    print("  A. baumannii / P. aeruginosa   %.2f   (published: 1.86)" % ratio_pa)

    substantially = ab_nf < 8.46
    host_structural = (ratio_kl is None) or ratio_kl < 1.5
    print("\nfrozen decision rule:")
    print("  enrichment is substantially composite structure (< 8.46) : %s" % substantially)
    print("  host contrast is composite structure (ratio < 1.5)       : %s" % host_structural)
    print("  referee objection answered                                : %s"
          % (not substantially and not host_structural))

    rec = {"receipt": "NM-C1-006 composite-element stratification",
           "amendment_sha256_of_body": claimed,
           "amendment_sha256_of_file": sha(AMEND),
           "inputs": {"null_matrices": sha(MATS), "is_elements": sha(ISTAB)},
           "complete_structural_elements_indexed": kept,
           "results": results,
           "non_flanked_summary": {"AB_enrichment": ab_nf,
                                   "AB_over_Klebsiella": round(ratio_kl, 3),
                                   "AB_over_Pseudomonas": round(ratio_pa, 3)},
           "decision_rule": doc["decision_rule"],
           "substantially_composite_structure": bool(substantially),
           "host_contrast_is_composite_structure": bool(host_structural),
           "referee_objection_answered": bool(not substantially and not host_structural)}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("\nreceipt: %s" % os.path.basename(OUT))


if __name__ == "__main__":
    main()
