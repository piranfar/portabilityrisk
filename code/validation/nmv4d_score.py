"""NM-V4D scoring -- three-architecture analysis. Runs only after the freeze.

Primary is BioProject-BALANCED: each BioProject contributes unit mass to its (family, host)
stratum, so effective n is the number of BioProjects, not the number of occurrences. That is
deliberately conservative and was fixed before scoring.
"""
import argparse, collections, csv, datetime, hashlib, json, math, os, sys
import numpy as np

VERSION = "nmv4d_score_v1.0.0"
FROZEN_SHA = "bfa88f6c55ec348c7f158af40326aeb17d18623f7f8d8a154bbbad4c9eeaafc5"
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
    R = S = 0.0; pr = ps = qr = qs = 0.0; w = []
    for a, b, c, d in tabs:
        n = a + b + c + d
        if n <= 0:
            w.append(0.0); continue
        Ri = a * d / n; Si = b * c / n
        P = (a + d) / n; Q = (b + c) / n
        R += Ri; S += Si
        pr += P * Ri; ps += P * Si; qr += Q * Ri; qs += Q * Si
        w.append(Ri + Si)
    if R <= 0 or S <= 0:
        return None
    o = R / S
    v = pr / (2 * R * R) + (ps + qr) / (2 * R * S) + qs / (2 * S * S)
    se = math.sqrt(v)
    return {"or": o, "ln_or": math.log(o), "se": se,
            "lo": math.exp(math.log(o) - 1.96 * se), "hi": math.exp(math.log(o) + 1.96 * se),
            "w": w, "tw": R + S}


def woolf(a, b, c, d):
    if min(a, b, c, d) <= 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    return math.log((a * d) / (b * c)), 1 / a + 1 / b + 1 / c + 1 / d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--frozen", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--perm", type=int, default=2000)
    a = ap.parse_args()
    if sha256_file(a.frozen) != FROZEN_SHA:
        print("REFUSING: frozen design digest mismatch"); sys.exit(1)
    F = json.load(open(a.frozen, encoding="utf-8"))
    if not F.get("frozen_before_any_outcome_was_scored"):
        print("REFUSING: design does not assert pre-outcome freeze"); sys.exit(1)
    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    for k, v in F["input_digests"].items():
        if sha256_file(os.path.join(O, k)) != v:
            print("REFUSING: input %s changed since the freeze" % k); sys.exit(2)
    print("%s | frozen design %s | %d inputs verified"
          % (VERSION, FROZEN_SHA[:16], len(F["input_digests"])))

    ELIG = set(F["eligibility"]["eligible_families"])
    THR = F["success_criteria"]["influence_threshold_relative_log_or"]
    MAXW = F["success_criteria"]["max_single_family_weight_share"]

    cls = list(csv.DictReader(open(os.path.join(O, "determinant_portability_classes.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    if len(cls) != 74349:
        print("REFUSING: PRIMARY is %d, expected 74349" % len(cls)); sys.exit(3)
    rec5 = collections.Counter(r["portability_class"] for r in cls)
    if sum(rec5.values()) != 74349 or set(rec5) != set("ABCDE"):
        print("REFUSING: A-E reconciliation failed"); sys.exit(4)
    tot3 = {k: sum(rec5[c] for c in v) for k, v in ARCH.items()}
    if sum(tot3.values()) != 74349:
        print("REFUSING: three-architecture partition does not sum to 74349"); sys.exit(5)
    print("  PRIMARY 74,349 | A-E %s" % " ".join("%s=%d" % (k, rec5[k]) for k in "ABCDE"))
    print("  three-architecture partition: quiescent %d | mobile %d | plasmid %d | sum %d"
          % (tot3["chromosomal_quiescent"], tot3["chromosomal_mobile"],
             tot3["plasmid_borne"], sum(tot3.values())))

    gbp = {}
    for r in csv.DictReader(open(os.path.join(O, "genome_level_summary.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        gbp[r["assembly_version"]] = r["bioproject_accession"]
    best = {}
    for asm, bp in gbp.items():
        h = hashlib.sha256(asm.encode()).hexdigest()
        if bp not in best or h < best[bp][0]:
            best[bp] = (h, asm)
    keep1 = {v[1] for v in best.values()}

    def hg(sp, gn):
        if sp == AB: return AB
        if sp == PA: return PA
        if gn == KL or sp.startswith(KL): return KL
        return None

    cell = collections.defaultdict(collections.Counter)
    cell2 = collections.defaultdict(collections.Counter)
    excl = collections.Counter()
    for r in cls:
        f = r["gene_family"]
        if f not in ELIG:
            excl["family_not_eligible"] += 1; continue
        h = hg(r["organism_harmonized"], r["genus"])
        if h is None:
            excl["host_not_a_comparison_group"] += 1; continue
        k = (f, h, r["bioproject_accession"])
        cell[k][r["portability_class"]] += 1
        if r["assembly_version"] in keep1:
            cell2[k][r["portability_class"]] += 1
    print("\n=== DENOMINATOR FLOW ===")
    print("  PRIMARY 74349")
    for k, v in excl.most_common():
        print("  excluded, %-38s %d" % (k, v))
    print("  retained %d | one-genome-per-BioProject %d"
          % (sum(sum(c.values()) for c in cell.values()),
             sum(sum(c.values()) for c in cell2.values())))

    def arch_count(c, name):
        return sum(c[x] for x in ARCH[name])

    def build(cells, h1, h2, a1, a2, balanced, drop_fam=None, drop_bp=None):
        """balanced: each BioProject contributes unit mass to its (family,host) stratum."""
        acc = collections.defaultdict(lambda: [0.0, 0.0])
        for (f, h, bp), c in cells.items():
            if f not in ELIG or h not in (h1, h2):
                continue
            if drop_fam is not None and f == drop_fam:
                continue
            if drop_bp is not None and bp == drop_bp:
                continue
            x = arch_count(c, a1); y = arch_count(c, a2)
            if x + y == 0:
                continue
            if balanced:
                x, y = x / (x + y), y / (x + y)
            v = acc[(f, h)]
            v[0] += x; v[1] += y
        T = []; fams = []
        for f in sorted({k[0] for k in acc}):
            p = acc.get((f, h1)); q = acc.get((f, h2))
            if not p or not q:
                continue
            if (p[0] + p[1]) <= 0 or (q[0] + q[1]) <= 0:
                continue
            T.append((p[0], p[1], q[0], q[1])); fams.append(f)
        return T, fams

    PCS = F["primary_comparisons"]
    results = {}
    famrows = []
    print("\n=== PRIMARY: BioProject-balanced, family-stratified MH ===")
    print("  %-5s %-58s %6s %10s %20s" % ("id", "comparison", "n_fam", "OR", "95% CI"))
    for pc in PCS:
        h1, h2 = pc["hosts"]; a1, a2 = pc["contrast"]
        for lab, cells, bal in (("PRIMARY_bp_balanced", cell, True),
                                ("S1_occurrence_level", cell, False),
                                ("S2_one_genome_per_bp", cell2, False)):
            T, fams = build(cells, h1, h2, a1, a2, bal)
            m = mh(T)
            key = "%s__%s" % (pc["id"], lab)
            if m is None:
                results[key] = None; continue
            lo = [woolf(*t) for t in T]
            lns = [l for l, _ in lo]; ws = [1 / v for _, v in lo]
            fx = sum(w * l for w, l in zip(ws, lns)) / sum(ws)
            Q = sum(w * (l - fx) ** 2 for w, l in zip(ws, lns)); df = len(T) - 1
            I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
            fam_w = float(np.mean(lns))
            results[key] = {"pc": pc["id"], "hosts": [h1, h2], "contrast": [a1, a2],
                            "analysis": lab, "n_families": len(T), "or": m["or"],
                            "ci_lo": m["lo"], "ci_hi": m["hi"], "ln_or": m["ln_or"],
                            "excludes_1": bool(m["lo"] > 1 or m["hi"] < 1),
                            "direction_positive": bool(m["or"] > 1),
                            "families_gt_1": sum(1 for l in lns if l > 0),
                            "families_lt_1": sum(1 for l in lns if l < 0),
                            "I2_pct": I2,
                            "max_family_weight_share": max(m["w"]) / m["tw"] if m["tw"] else 1.0,
                            "family_weighted_mean_ln_or": fam_w,
                            "family_weighted_or": math.exp(fam_w)}
            if lab == "PRIMARY_bp_balanced":
                print("  %-5s %-58s %6d %10.3f  [%7.3f, %8.3f] %s"
                      % (pc["id"], "%s vs %s : %s vs %s" % (h1, h2, a1, a2), len(T),
                         m["or"], m["lo"], m["hi"],
                         "*" if (m["lo"] > 1 or m["hi"] < 1) else ""))
                for i, f in enumerate(fams):
                    x1, y1, x2, y2 = T[i]
                    famrows.append({"comparison": pc["id"], "gene_family": f,
                                    "host1": h1, "host2": h2, "arch1": a1, "arch2": a2,
                                    "h1_arch1": round(x1, 6), "h1_arch2": round(y1, 6),
                                    "h2_arch1": round(x2, 6), "h2_arch2": round(y2, 6),
                                    "ln_or": round(lo[i][0], 6),
                                    "or": round(math.exp(lo[i][0]), 6),
                                    "mh_weight_share": round(m["w"][i] / m["tw"], 6)})
    print("\n  sensitivities:")
    for pc in PCS:
        for lab in ("S1_occurrence_level", "S2_one_genome_per_bp"):
            r = results.get("%s__%s" % (pc["id"], lab))
            if r:
                print("  %-5s %-24s n_fam %4d  OR %9.3f [%7.3f, %9.3f] %s"
                      % (pc["id"], lab, r["n_families"], r["or"], r["ci_lo"], r["ci_hi"],
                         "*" if r["excludes_1"] else ""))
    print("\n  family-weighted summary (each family counts once, S3):")
    for pc in PCS:
        r = results["%s__PRIMARY_bp_balanced" % pc["id"]]
        print("  %-5s MH OR %9.3f | family-weighted OR %9.3f | I2 %5.1f%% | max fam weight %.3f"
              % (pc["id"], r["or"], r["family_weighted_or"], r["I2_pct"],
                 r["max_family_weight_share"]))

    # ---------------- influence ----------------
    print("\n=== INFLUENCE (primary BioProject-balanced) ===")
    infl = {}
    for pc in PCS:
        h1, h2 = pc["hosts"]; a1, a2 = pc["contrast"]
        base = results["%s__PRIMARY_bp_balanced" % pc["id"]]["ln_or"]
        T0, f0 = build(cell, h1, h2, a1, a2, True)
        lof = []
        for f in f0:
            T, _ = build(cell, h1, h2, a1, a2, True, drop_fam=f)
            m = mh(T)
            if m:
                lof.append((f, abs(m["ln_or"] - base) / abs(base)))
        bps = sorted({bp for (f, h, bp) in cell if h in (h1, h2) and f in ELIG})
        lob = []
        for bp in bps:
            T, _ = build(cell, h1, h2, a1, a2, True, drop_bp=bp)
            m = mh(T)
            if m:
                lob.append((bp, abs(m["ln_or"] - base) / abs(base)))
        lof.sort(key=lambda x: -x[1]); lob.sort(key=lambda x: -x[1])
        infl[pc["id"]] = {"lofo_max": lof[0][1], "lofo_worst": lof[0][0],
                          "lofo_over": sum(1 for _, v in lof if v > THR), "lofo_n": len(lof),
                          "lobo_max": lob[0][1], "lobo_worst": lob[0][0],
                          "lobo_over": sum(1 for _, v in lob if v > THR), "lobo_n": len(lob)}
        print("  %-5s LOFO n=%3d max %.4f (%s) over=%d | LOBO n=%4d max %.4f (%s) over=%d"
              % (pc["id"], len(lof), lof[0][1], lof[0][0], infl[pc["id"]]["lofo_over"],
                 len(lob), lob[0][1], lob[0][0], infl[pc["id"]]["lobo_over"]))

    # ---------------- permutation ----------------
    print("\n=== PERMUTATION, host labels over BioProjects within family (S6) ===")
    rng = np.random.default_rng(20260821)
    perms = {}
    for pc in PCS:
        h1, h2 = pc["hosts"]; a1, a2 = pc["contrast"]
        obs = abs(results["%s__PRIMARY_bp_balanced" % pc["id"]]["ln_or"])
        units = collections.defaultdict(list)
        for (f, h, bp), c in cell.items():
            if f in ELIG and h in (h1, h2):
                x = arch_count(c, a1); y = arch_count(c, a2)
                if x + y > 0:
                    units[f].append((h, x / (x + y), y / (x + y)))
        ge = done = 0
        for _ in range(a.perm):
            T = []
            for f, lst in units.items():
                labs = [h for h, _, _ in lst]
                rng.shuffle(labs)
                p = [0.0, 0.0]; q = [0.0, 0.0]
                for (_, x, y), nh in zip(lst, labs):
                    t = p if nh == h1 else q
                    t[0] += x; t[1] += y
                if (p[0] + p[1]) > 0 and (q[0] + q[1]) > 0:
                    T.append((p[0], p[1], q[0], q[1]))
            m = mh(T)
            if m:
                done += 1
                if abs(m["ln_or"]) >= obs:
                    ge += 1
        p_ = (ge + 1) / (done + 1)
        perms[pc["id"]] = {"n": done, "ge": ge, "p": p_}
        print("  %-5s permutations %4d | >= observed %4d | empirical p = %.5f"
              % (pc["id"], done, ge, p_))

    # ---------------- five-class coherence ----------------
    print("\n=== SECONDARY FULL FIVE-CLASS (class A and C both retained) ===")
    five = {}
    for h in (AB, PA, KL):
        c = collections.Counter()
        for (f, hh, bp), cc in cell.items():
            if hh == h and f in ELIG:
                c.update(cc)
        n = sum(c.values()); five[h] = {k: c[k] for k in "ABCDE"}
        print("  %-26s n=%6d  A %5.1f%%  B %5.1f%%  C %5.1f%%  D %5.1f%%  E %5.1f%%"
              % (h, n, *[100 * c[k] / n for k in "ABCDE"]))
    def frac(h, name):
        t = sum(five[h].values())
        return sum(five[h][x] for x in ARCH[name]) / t
    coh = {
      "PC1": frac(AB, "chromosomal_mobile") / max(frac(AB, "chromosomal_quiescent"), 1e-9)
             > frac(PA, "chromosomal_mobile") / max(frac(PA, "chromosomal_quiescent"), 1e-9),
      "PC2": frac(AB, "chromosomal_mobile") / max(frac(AB, "plasmid_borne"), 1e-9)
             > frac(KL, "chromosomal_mobile") / max(frac(KL, "plasmid_borne"), 1e-9),
      "PC3": frac(PA, "chromosomal_quiescent") / max(frac(PA, "plasmid_borne"), 1e-9)
             > frac(KL, "chromosomal_quiescent") / max(frac(KL, "plasmid_borne"), 1e-9)}
    print("  five-class coherence with three-class direction: %s" % coh)

    # ---------------- criteria ----------------
    def R(pc, lab="PRIMARY_bp_balanced"):
        return results["%s__%s" % (pc, lab)]
    SC = {
      "SC1": bool(R("PC1")["direction_positive"] and R("PC1")["excludes_1"]),
      "SC2": bool(R("PC2")["direction_positive"] and R("PC2")["excludes_1"]),
      "SC3": bool(R("PC3")["direction_positive"] and R("PC3")["excludes_1"]),
      "SC4": bool(all(R(p, l)["direction_positive"] and R(p, l)["excludes_1"]
                      for p in ("PC1", "PC2", "PC3")
                      for l in ("PRIMARY_bp_balanced", "S1_occurrence_level",
                                "S2_one_genome_per_bp"))),
      "SC5": bool(all(infl[p]["lofo_max"] <= THR and infl[p]["lobo_max"] <= THR
                      and R(p)["max_family_weight_share"] <= MAXW
                      for p in ("PC1", "PC2", "PC3"))),
      "SC6": bool(all(coh.values())),
      "SC7": "deferred to nmv4d_verify.py"}
    print("\n=== SUCCESS CRITERIA ===")
    for k in ("SC1", "SC2", "SC3", "SC4", "SC5", "SC6", "SC7"):
        print("  %-5s %s" % (k, SC[k]))
    core = all(SC[k] for k in ("SC1", "SC2", "SC3", "SC4", "SC5", "SC6"))
    verdict = ("NMV4D_SUCCESS_THREE_ARCHITECTURES_SUPPORTED" if core else
               "NMV4D_PARTIAL_ARCHITECTURE_MODEL"
               if any(SC[k] for k in ("SC1", "SC2", "SC3")) else "NMV4D_FAIL")
    print("  VERDICT (pending SC7): %s" % verdict)

    os.makedirs(a.outdir, exist_ok=True)
    p1 = os.path.join(a.outdir, "nmv4d_family_three_architecture.tsv")
    c1 = ["comparison", "gene_family", "host1", "host2", "arch1", "arch2",
          "h1_arch1", "h1_arch2", "h2_arch1", "h2_arch2", "or", "ln_or", "mh_weight_share"]
    with open(p1, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(c1) + "\n")
        for r in famrows:
            fh.write("\t".join(str(r[c]) for c in c1) + "\n")
    p2 = os.path.join(a.outdir, "nmv4d_host_three_architecture_summary.tsv")
    c2 = ["comparison", "hosts", "contrast", "analysis", "n_families", "or", "ci_lo", "ci_hi",
          "excludes_1", "families_gt_1", "families_lt_1", "I2_pct",
          "max_family_weight_share", "family_weighted_or"]
    with open(p2, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(c2) + "\n")
        for k, v in results.items():
            if not v:
                continue
            fh.write("\t".join(str(x) for x in [
                v["pc"], " vs ".join(v["hosts"]), " vs ".join(v["contrast"]), v["analysis"],
                v["n_families"], round(v["or"], 6), round(v["ci_lo"], 6), round(v["ci_hi"], 6),
                v["excludes_1"], v["families_gt_1"], v["families_lt_1"], round(v["I2_pct"], 3),
                round(v["max_family_weight_share"], 6),
                round(v["family_weighted_or"], 6)]) + "\n")
    rec = {"builder": VERSION,
           "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "frozen_design_sha256": FROZEN_SHA, "primary_denominator_verified": 74349,
           "class_reconciliation_five": dict(rec5),
           "three_architecture_partition": tot3,
           "denominator_flow": {"primary": 74349, "exclusions": dict(excl),
                                "retained": sum(sum(c.values()) for c in cell.values())},
           "results": results, "influence": infl, "permutation": perms,
           "five_class_by_host": five, "five_class_coherence": coh,
           "success_criteria": SC, "verdict_pending_verification": verdict,
           "outputs": {"family_tsv": sha256_file(p1), "summary_tsv": sha256_file(p2)},
           "statements": ["No transfer, conjugation or HGT event was observed or is claimed.",
                          "Host association is statistical, not causal control.",
                          "Class A and class C were both retained; A and B were never merged.",
                          "Every chromosomal_mobile conclusion is conditional on NM-V1.",
                          "No PortabilityEvent or PlasmidCall artefact was read."]}
    rp = os.path.join(a.outdir, "NMV4D_RESULT_RECEIPT.json")
    json.dump(rec, open(rp, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  %s  %s" % (os.path.basename(p1), sha256_file(p1)))
    print("  %s  %s" % (os.path.basename(p2), sha256_file(p2)))
    print("  %s  %s" % (os.path.basename(rp), sha256_file(rp)))


if __name__ == "__main__":
    main()
