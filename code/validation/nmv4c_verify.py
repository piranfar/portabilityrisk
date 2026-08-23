"""NM-V4C independent verifier.

Imports nothing from the scoring script and touches no in-memory model object. It re-reads the
exported per-family 2x2 counts from nmv4c_family_host_vehicle.tsv, re-implements the
Mantel-Haenszel estimator from the published formula, and checks the result against the
summary TSV and the receipt. It also re-derives the per-family counts from the frozen class
table by an independent path, so a shared aggregation bug would not survive both.
"""
import argparse, collections, csv, hashlib, json, math, os, sys

VERSION = "nmv4c_verify_v1.0.0"
AB, PA, KL = "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Klebsiella"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def mh_or(tabs):
    R = S = 0.0
    pr = ps = qr = qs = 0.0
    for a, b, c, d in tabs:
        n = a + b + c + d
        if not n:
            continue
        Ri = a * d / n; Si = b * c / n
        P = (a + d) / n; Q = (b + c) / n
        R += Ri; S += Si
        pr += P * Ri; ps += P * Si; qr += Q * Ri; qs += Q * Si
    o = R / S
    v = pr / (2 * R * R) + (ps + qr) / (2 * R * S) + qs / (2 * S * S)
    se = math.sqrt(v)
    return o, math.exp(math.log(o) - 1.96 * se), math.exp(math.log(o) + 1.96 * se)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    D = a.dir
    fam = list(csv.DictReader(open(os.path.join(D, "nmv4c_family_host_vehicle.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    summ = list(csv.DictReader(open(os.path.join(D, "nmv4c_host_vehicle_summary.tsv"),
                                    encoding="utf-8"), delimiter="\t"))
    rec = json.load(open(os.path.join(D, "NMV4C_RESULT_RECEIPT.json"), encoding="utf-8"))
    print("%s\n  re-reading %d family rows from the exported TSV" % (VERSION, len(fam)))

    tabs = [(int(r["h1_classB"]), int(r["h1_plasmid_CDE"]),
             int(r["h2_classB"]), int(r["h2_plasmid_CDE"])) for r in fam]
    o, lo, hi = mh_or(tabs)
    srow = [s for s in summ if s["contrast"] == "%s vs %s" % (AB, KL)
            and s["analysis"] == "P1_full"][0]
    key = "Acinetobacter_baumannii__vs__Klebsiella__P1_full"
    print("\n=== headline recomputed from the TSV, not from memory ===")
    print("  %-34s %12s %12s %12s" % ("source", "OR", "CI lo", "CI hi"))
    print("  %-34s %12.4f %12.4f %12.4f" % ("verifier, from family TSV", o, lo, hi))
    print("  %-34s %12.4f %12.4f %12.4f" % ("summary TSV", float(srow["or"]),
                                            float(srow["ci_lo"]), float(srow["ci_hi"])))
    print("  %-34s %12.4f %12.4f %12.4f" % ("receipt", rec["results"][key]["or"],
                                            rec["results"][key]["ci_lo"],
                                            rec["results"][key]["ci_hi"]))
    ok = (abs(o - float(srow["or"])) < 1e-3
          and abs(o - rec["results"][key]["or"]) < 1e-3
          and abs(lo - float(srow["ci_lo"])) < 1e-3
          and abs(hi - float(srow["ci_hi"])) < 1e-3)
    print("  AGREEMENT: %s" % ("MATCH" if ok else "*** MISMATCH ***"))

    nfam = len(fam)
    up = sum(1 for r in fam if float(r["ln_or"]) > 0)
    dn = sum(1 for r in fam if float(r["ln_or"]) < 0)
    mxw = max(float(r["mh_weight_share"]) for r in fam)
    print("\n  %-40s verifier %6s | receipt %6s | %s"
          % ("n families", nfam, srow["n_families"],
             "MATCH" if nfam == int(srow["n_families"]) else "MISMATCH"))
    print("  %-40s verifier %6d | receipt %6s | %s"
          % ("families with OR > 1", up, srow["families_or_gt_1"],
             "MATCH" if up == int(srow["families_or_gt_1"]) else "MISMATCH"))
    print("  %-40s verifier %6d | receipt %6s | %s"
          % ("families with OR < 1", dn, srow["families_or_lt_1"],
             "MATCH" if dn == int(srow["families_or_lt_1"]) else "MISMATCH"))
    print("  %-40s verifier %6.4f | receipt %6.4f | %s"
          % ("max family weight share", mxw, float(srow["max_family_weight_share"]),
             "MATCH" if abs(mxw - float(srow["max_family_weight_share"])) < 1e-4 else "MISMATCH"))
    print("  %-40s %s" % ("weight shares sum to 1",
                          "MATCH" if abs(sum(float(r["mh_weight_share"]) for r in fam) - 1) < 1e-3
                          else "*** MISMATCH ***"))

    # ---- independent re-derivation of the per-family counts from the frozen class table ----
    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    cls = list(csv.DictReader(open(os.path.join(O, "determinant_portability_classes.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    print("\n=== independent re-derivation from the frozen class table ===")
    print("  class table rows: %d %s" % (len(cls), "MATCH" if len(cls) == 74349 else "MISMATCH"))
    want = {r["gene_family"] for r in fam}
    agg = collections.defaultdict(collections.Counter)
    for r in cls:
        f = r["gene_family"]
        if f not in want:
            continue
        sp, gn = r["organism_harmonized"], r["genus"]
        h = AB if sp == AB else (KL if (gn == KL or sp.startswith(KL)) else None)
        if h:
            agg[(f, h)][r["portability_class"]] += 1
    bad = 0
    for r in fam:
        f = r["gene_family"]
        c1 = agg[(f, AB)]; c2 = agg[(f, KL)]
        if (c1["B"] != int(r["h1_classB"])
                or c1["C"] + c1["D"] + c1["E"] != int(r["h1_plasmid_CDE"])
                or c2["B"] != int(r["h2_classB"])
                or c2["C"] + c2["D"] + c2["E"] != int(r["h2_plasmid_CDE"])):
            bad += 1
            if bad <= 3:
                print("     MISMATCH on %s" % f)
    print("  per-family 2x2 cells re-derived independently: %d of %d disagree" % (bad, len(fam)))
    tabs2 = []
    for r in fam:
        f = r["gene_family"]
        c1 = agg[(f, AB)]; c2 = agg[(f, KL)]
        tabs2.append((c1["B"], c1["C"] + c1["D"] + c1["E"],
                      c2["B"], c2["C"] + c2["D"] + c2["E"]))
    o2, lo2, hi2 = mh_or(tabs2)
    print("  MH OR from the independent re-derivation: %.4f [%.4f, %.4f]  %s"
          % (o2, lo2, hi2, "MATCH" if abs(o2 - o) < 1e-6 else "*** MISMATCH ***"))

    print("\n=== co-primary agreement audit ===")
    for con in ("%s vs %s" % (AB, KL), "%s vs %s" % (AB, PA), "%s vs %s" % (PA, KL)):
        rows = {s["analysis"]: s for s in summ if s["contrast"] == con}
        p1, p2 = rows.get("P1_full"), rows.get("P2_one_genome_per_bioproject")
        if not (p1 and p2):
            continue
        agree = (p1["excludes_1"] == p2["excludes_1"]
                 and (float(p1["or"]) > 1) == (float(p2["or"]) > 1))
        print("  %-58s P1 OR %8.3f excl1=%-5s | P2 OR %8.3f excl1=%-5s | co-primary %s"
              % (con, float(p1["or"]), p1["excludes_1"], float(p2["or"]),
                 p2["excludes_1"], "AGREE" if agree else "*** DISAGREE ***"))

    allok = ok and bad == 0 and abs(o2 - o) < 1e-6 and len(cls) == 74349
    print("\n  VERIFIER VERDICT: %s" % ("PASS" if allok else "*** FAIL ***"))
    sys.exit(0 if allok else 5)


if __name__ == "__main__":
    main()
