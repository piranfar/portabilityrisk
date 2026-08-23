"""NM-V4D independent verifier (success criterion SC7).

Imports nothing from the scoring script. Recomputes every headline effect from the exported
TSVs, and independently re-derives the BioProject-balanced per-family cells from the frozen
class table by a separate code path.
"""
import argparse, collections, csv, hashlib, json, math, os, sys

VERSION = "nmv4d_verify_v1.0.0"
AB, PA, KL = "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Klebsiella"
ARCH = {"chromosomal_quiescent": ["A"], "chromosomal_mobile": ["B"],
        "plasmid_borne": ["C", "D", "E"]}


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def mh(tabs):
    R = S = 0.0; pr = ps = qr = qs = 0.0
    for a, b, c, d in tabs:
        n = a + b + c + d
        if n <= 0:
            continue
        Ri = a * d / n; Si = b * c / n
        P = (a + d) / n; Q = (b + c) / n
        R += Ri; S += Si
        pr += P * Ri; ps += P * Si; qr += Q * Ri; qs += Q * Si
    o = R / S
    se = math.sqrt(pr / (2 * R * R) + (ps + qr) / (2 * R * S) + qs / (2 * S * S))
    return o, math.exp(math.log(o) - 1.96 * se), math.exp(math.log(o) + 1.96 * se)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    D = a.dir
    fam = list(csv.DictReader(open(os.path.join(D, "nmv4d_family_three_architecture.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    summ = list(csv.DictReader(open(os.path.join(D, "nmv4d_host_three_architecture_summary.tsv"),
                                    encoding="utf-8"), delimiter="\t"))
    rec = json.load(open(os.path.join(D, "NMV4D_RESULT_RECEIPT.json"), encoding="utf-8"))
    print("%s\n  family rows re-read from TSV: %d" % (VERSION, len(fam)))

    ok = True
    print("\n=== headline effects recomputed from the exported TSV ===")
    print("  %-5s %12s %12s %12s | %12s %12s %12s | %s"
          % ("id", "OR(verif)", "lo", "hi", "OR(summary)", "lo", "hi", "match"))
    for pc in ("PC1", "PC2", "PC3"):
        rows = [r for r in fam if r["comparison"] == pc]
        tabs = [(float(r["h1_arch1"]), float(r["h1_arch2"]),
                 float(r["h2_arch1"]), float(r["h2_arch2"])) for r in rows]
        o, lo, hi = mh(tabs)
        s = [x for x in summ if x["comparison"] == pc
             and x["analysis"] == "PRIMARY_bp_balanced"][0]
        m = (abs(o - float(s["or"])) < 1e-3 and abs(lo - float(s["ci_lo"])) < 1e-3
             and abs(hi - float(s["ci_hi"])) < 1e-3
             and abs(o - rec["results"]["%s__PRIMARY_bp_balanced" % pc]["or"]) < 1e-3)
        ok &= m
        print("  %-5s %12.4f %12.4f %12.4f | %12.4f %12.4f %12.4f | %s"
              % (pc, o, lo, hi, float(s["or"]), float(s["ci_lo"]), float(s["ci_hi"]),
                 "MATCH" if m else "*** MISMATCH ***"))
        nf = len(rows)
        m2 = nf == int(s["n_families"])
        up = sum(1 for r in rows if float(r["ln_or"]) > 0)
        m3 = up == int(s["families_gt_1"])
        mw = max(float(r["mh_weight_share"]) for r in rows)
        m4 = abs(mw - float(s["max_family_weight_share"])) < 1e-4
        ws = sum(float(r["mh_weight_share"]) for r in rows)
        m5 = abs(ws - 1.0) < 1e-3
        ok &= (m2 and m3 and m4 and m5)
        print("        n_fam %d/%s %s | OR>1 %d/%s %s | maxw %.4f/%.4f %s | weights sum %.4f %s"
              % (nf, s["n_families"], "OK" if m2 else "BAD", up, s["families_gt_1"],
                 "OK" if m3 else "BAD", mw, float(s["max_family_weight_share"]),
                 "OK" if m4 else "BAD", ws, "OK" if m5 else "BAD"))

    # ---- independent re-derivation of the BioProject-balanced cells ----
    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    cls = list(csv.DictReader(open(os.path.join(O, "determinant_portability_classes.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    print("\n=== independent re-derivation from the frozen class table ===")
    print("  class table rows %d %s" % (len(cls), "MATCH" if len(cls) == 74349 else "MISMATCH"))
    ok &= len(cls) == 74349
    tot3 = collections.Counter()
    for r in cls:
        for k, v in ARCH.items():
            if r["portability_class"] in v:
                tot3[k] += 1
    print("  partition: quiescent %d | mobile %d | plasmid %d | sum %d %s"
          % (tot3["chromosomal_quiescent"], tot3["chromosomal_mobile"], tot3["plasmid_borne"],
             sum(tot3.values()), "MATCH" if sum(tot3.values()) == 74349 else "MISMATCH"))
    ok &= sum(tot3.values()) == 74349

    def host(sp, gn):
        if sp == AB: return AB
        if sp == PA: return PA
        if gn == KL or sp.startswith(KL): return KL
        return None

    raw = collections.defaultdict(collections.Counter)
    for r in cls:
        h = host(r["organism_harmonized"], r["genus"])
        if h:
            raw[(r["gene_family"], h, r["bioproject_accession"])][r["portability_class"]] += 1

    bad = 0; checked = 0
    for pc in ("PC1", "PC2", "PC3"):
        rows = [r for r in fam if r["comparison"] == pc]
        if not rows:
            continue
        h1, h2 = rows[0]["host1"], rows[0]["host2"]
        a1, a2 = rows[0]["arch1"], rows[0]["arch2"]
        acc = collections.defaultdict(lambda: [0.0, 0.0])
        for (f, h, bp), c in raw.items():
            if h not in (h1, h2):
                continue
            x = sum(c[z] for z in ARCH[a1]); y = sum(c[z] for z in ARCH[a2])
            if x + y == 0:
                continue
            v = acc[(f, h)]
            v[0] += x / (x + y); v[1] += y / (x + y)
        tabs = []
        for r in rows:
            f = r["gene_family"]
            p = acc.get((f, h1), [0, 0]); q = acc.get((f, h2), [0, 0])
            checked += 1
            for got, want in ((p[0], "h1_arch1"), (p[1], "h1_arch2"),
                              (q[0], "h2_arch1"), (q[1], "h2_arch2")):
                if abs(got - float(r[want])) > 1e-6:
                    bad += 1
                    if bad <= 3:
                        print("     MISMATCH %s %s %s: %.6f vs %s" % (pc, f, want, got, r[want]))
                    break
            tabs.append((p[0], p[1], q[0], q[1]))
        o2, lo2, hi2 = mh(tabs)
        s = [x for x in summ if x["comparison"] == pc
             and x["analysis"] == "PRIMARY_bp_balanced"][0]
        m = abs(o2 - float(s["or"])) < 1e-3
        ok &= m
        print("  %-5s re-derived OR %12.4f [%10.4f, %12.4f]  %s"
              % (pc, o2, lo2, hi2, "MATCH" if m else "*** MISMATCH ***"))
    print("  family cells re-derived independently: %d checked, %d disagree" % (checked, bad))
    ok &= bad == 0

    print("\n=== success-criteria audit ===")
    sc = rec["success_criteria"]
    for k in ("SC1", "SC2", "SC3", "SC4", "SC5", "SC6"):
        print("  %-5s %s" % (k, sc[k]))
    print("  SC7   %s" % ("PASS" if ok else "FAIL"))
    print("\n  per-comparison robustness detail:")
    for pc in ("PC1", "PC2", "PC3"):
        i = rec["influence"][pc]
        s = [x for x in summ if x["comparison"] == pc
             and x["analysis"] == "PRIMARY_bp_balanced"][0]
        print("  %-5s maxFamWeight %.3f %-4s | LOFO %.4f %-4s | LOBO %.4f %-4s | fam OR>1 %s/%s"
              % (pc, float(s["max_family_weight_share"]),
                 "OK" if float(s["max_family_weight_share"]) <= 0.30 else "OVER",
                 i["lofo_max"], "OK" if i["lofo_max"] <= 0.20 else "OVER",
                 i["lobo_max"], "OK" if i["lobo_max"] <= 0.20 else "OVER",
                 s["families_gt_1"], s["n_families"]))
    print("\n  VERIFIER VERDICT: %s" % ("PASS" if ok else "*** FAIL ***"))
    sys.exit(0 if ok else 6)


if __name__ == "__main__":
    main()
