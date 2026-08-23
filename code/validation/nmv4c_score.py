"""NM-V4C scoring -- within-determinant host-vehicle analysis. Runs only after the freeze.

Emits family-level effects, pooled effects, heterogeneity and influence diagnostics, plus two
TSVs from which an independent verifier recomputes the headline values. Writes no prose.
"""
import argparse, collections, csv, datetime, hashlib, json, math, os, sys
import numpy as np

VERSION = "nmv4c_score_v1.0.0"
FROZEN_SHA = "8a3c76b157cbf2cd5279342cb2752e1974a2baac976b62455ea0c66f3de4495d"
AB, PA, KL = "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Klebsiella"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def mh(tables):
    """Mantel-Haenszel OR with Robins-Breslow-Greenland variance. tables: list of (a,b,c,d)."""
    R = S = 0.0; pr = ps = qr = qs = 0.0; wts = []
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            wts.append(0.0); continue
        Ri = a * d / n; Si = b * c / n
        P = (a + d) / n; Q = (b + c) / n
        R += Ri; S += Si
        pr += P * Ri; ps += P * Si; qr += Q * Ri; qs += Q * Si
        wts.append(Ri + Si)
    if R <= 0 or S <= 0:
        return None
    orr = R / S
    var = pr / (2 * R * R) + (ps + qr) / (2 * R * S) + qs / (2 * S * S)
    se = math.sqrt(var)
    return {"or": orr, "ln_or": math.log(orr), "se_ln": se,
            "ci_lo": math.exp(math.log(orr) - 1.96 * se),
            "ci_hi": math.exp(math.log(orr) + 1.96 * se),
            "weights": wts, "total_weight": R + S}


def woolf(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    l = math.log((a * d) / (b * c))
    v = 1 / a + 1 / b + 1 / c + 1 / d
    return l, v


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
    print("%s | frozen design verified %s" % (VERSION, FROZEN_SHA[:16]))

    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    for k, v in F["input_digests"].items():
        if sha256_file(os.path.join(O, k)) != v:
            print("REFUSING: input %s changed since the freeze" % k); sys.exit(2)
    print("  all %d frozen input digests verified" % len(F["input_digests"]))

    ELIG = set(F["eligibility"]["eligible_families"])
    THR = F["gates"]["influence_threshold_relative_log_or"]
    MINFAM = F["ab_specific_test"]["minimum_supporting_families"]
    SEED = F["primary_test"]["global_test"]["seed"]

    cls = list(csv.DictReader(open(os.path.join(O, "determinant_portability_classes.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    if len(cls) != 74349:
        print("REFUSING: PRIMARY is %d, expected 74349" % len(cls)); sys.exit(3)
    recon = collections.Counter(r["portability_class"] for r in cls)
    if sum(recon.values()) != 74349 or set(recon) != set("ABCDE"):
        print("REFUSING: A-E reconciliation failed"); sys.exit(4)
    print("  PRIMARY 74,349 verified | A-E reconciliation %s"
          % " ".join("%s=%d" % (k, recon[k]) for k in "ABCDE"))

    # deterministic one-genome-per-BioProject selection
    gen_bp = {}
    for r in csv.DictReader(open(os.path.join(O, "genome_level_summary.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        gen_bp[r["assembly_version"]] = r["bioproject_accession"]
    best = {}
    for asm, bp in gen_bp.items():
        k = (bp,)
        h = hashlib.sha256(asm.encode()).hexdigest()
        if k not in best or h < best[k][0]:
            best[k] = (h, asm)
    keep_one = {v[1] for v in best.values()}

    def hostgroup(sp, genus):
        if sp == AB: return AB
        if sp == PA: return PA
        if genus == KL or sp.startswith(KL): return KL
        return None

    # cell[(family, host, bioproject)] -> Counter over classes ; and the P2 variant
    cell = collections.defaultdict(collections.Counter)
    cell2 = collections.defaultdict(collections.Counter)
    excl = collections.Counter()
    for r in cls:
        fam = r["gene_family"]
        if fam not in ELIG:
            excl["family_not_eligible"] += 1; continue
        hg = hostgroup(r["organism_harmonized"], r["genus"])
        if hg is None:
            excl["host_not_a_comparison_group"] += 1; continue
        key = (fam, hg, r["bioproject_accession"])
        cell[key][r["portability_class"]] += 1
        if r["assembly_version"] in keep_one:
            cell2[key][r["portability_class"]] += 1
    kept = sum(sum(c.values()) for c in cell.values())
    print("\n=== DENOMINATOR FLOW ===")
    print("  PRIMARY occurrences                         : 74349")
    for k, v in excl.most_common():
        print("  excluded, %-40s: %d" % (k, v))
    print("  retained in NM-V4C                          : %d" % kept)
    print("  retained after one-genome-per-BioProject    : %d"
          % sum(sum(c.values()) for c in cell2.values()))

    def agg(cells, fams, host):
        out = collections.Counter()
        for (f, h, bp), c in cells.items():
            if h == host and f in fams:
                out.update(c)
        return out

    def tables_for(cells, h1, h2, fams, drop_bp=None, drop_fam=None):
        per = collections.defaultdict(lambda: collections.Counter())
        for (f, h, bp), c in cells.items():
            if f not in fams or h not in (h1, h2):
                continue
            if drop_bp is not None and bp == drop_bp:
                continue
            if drop_fam is not None and f == drop_fam:
                continue
            per[(f, h)].update(c)
        T = []; fl = []
        for f in sorted(fams):
            if drop_fam is not None and f == drop_fam:
                continue
            c1 = per.get((f, h1), collections.Counter())
            c2 = per.get((f, h2), collections.Counter())
            a_ = c1["B"]; b_ = c1["C"] + c1["D"] + c1["E"]
            c_ = c2["B"]; d_ = c2["C"] + c2["D"] + c2["E"]
            if (a_ + b_) == 0 or (c_ + d_) == 0:
                continue
            T.append((a_, b_, c_, d_)); fl.append(f)
        return T, fl

    CONTRASTS = [(AB, KL), (AB, PA), (PA, KL)]
    results = {}
    famrows = []
    print("\n=== PRIMARY: family-stratified Mantel-Haenszel, B vs C+D+E ===")
    for h1, h2 in CONTRASTS:
        for lab, cells in (("P1_full", cell), ("P2_one_genome_per_bioproject", cell2)):
            T, fl = tables_for(cells, h1, h2, ELIG)
            m = mh(T)
            key = "%s__vs__%s__%s" % (h1.replace(" ", "_"), h2.replace(" ", "_"), lab)
            if m is None:
                results[key] = None; continue
            lo = [woolf(*t) for t in T]
            wf = [1 / v for _, v in lo]
            lnf = [l for l, _ in lo]
            fixed = sum(w * l for w, l in zip(wf, lnf)) / sum(wf)
            Q = sum(w * (l - fixed) ** 2 for w, l in zip(wf, lnf))
            df = len(T) - 1
            I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
            up = sum(1 for l in lnf if l > 0); dn = sum(1 for l in lnf if l < 0)
            mxw = max(m["weights"]) / m["total_weight"] if m["total_weight"] else 1.0
            results[key] = {"host1": h1, "host2": h2, "analysis": lab,
                            "n_families": len(T), "or": m["or"],
                            "ci_lo": m["ci_lo"], "ci_hi": m["ci_hi"],
                            "ln_or": m["ln_or"], "se_ln": m["se_ln"],
                            "excludes_1": bool(m["ci_lo"] > 1 or m["ci_hi"] < 1),
                            "families_or_gt_1": up, "families_or_lt_1": dn,
                            "cochran_Q": Q, "df": df, "I2_pct": I2,
                            "max_family_weight_share": mxw,
                            "largest_weight_family": fl[int(np.argmax(m["weights"]))] if fl else None}
            print("  %-28s vs %-24s %-30s n_fam %3d  OR %7.3f [%6.3f, %7.3f] %s"
                  % (h1, h2, lab, len(T), m["or"], m["ci_lo"], m["ci_hi"],
                     "*" if (m["ci_lo"] > 1 or m["ci_hi"] < 1) else ""))
            if lab == "P1_full" and (h1, h2) == (AB, KL):
                for i, f in enumerate(fl):
                    aa, bb, cc, dd = T[i]
                    l, v = lo[i]
                    famrows.append({"gene_family": f, "host1": h1, "host2": h2,
                                    "h1_classB": aa, "h1_plasmid_CDE": bb,
                                    "h2_classB": cc, "h2_plasmid_CDE": dd,
                                    "ln_or": round(l, 6), "var_ln_or": round(v, 6),
                                    "or": round(math.exp(l), 6),
                                    "mh_weight_share": round(m["weights"][i] / m["total_weight"], 6)})

    # ---------- five-class composition and class A, per family x host ----------
    comp = {}
    for host in (AB, PA, KL):
        for f in sorted(ELIG):
            c = agg(cell, {f}, host)
            n = sum(c.values())
            if n:
                comp[(f, host)] = (c, n)
    print("\n=== FIVE-CLASS REPRESENTATION (class A retained, reported separately) ===")
    for host in (AB, PA, KL):
        tot = collections.Counter()
        for f in ELIG:
            if (f, host) in comp:
                tot.update(comp[(f, host)][0])
        n = sum(tot.values())
        print("  %-26s n=%6d  A %5.1f%%  B %5.1f%%  C %5.1f%%  D %5.1f%%  E %5.1f%%"
              % (host, n, *[100 * tot[k] / n for k in "ABCDE"]))

    # total variation distance between hosts, per family, on the full 5-class vector
    tvd = {}
    for f in sorted(ELIG):
        for (x, y) in CONTRASTS:
            if (f, x) in comp and (f, y) in comp:
                cx, nx = comp[(f, x)]; cy, ny = comp[(f, y)]
                tvd[(f, x, y)] = 0.5 * sum(abs(cx[k] / nx - cy[k] / ny) for k in "ABCDE")

    # ---------- influence: leave-one-family-out and leave-one-BioProject-out ----------
    print("\n=== INFLUENCE DIAGNOSTICS (A. baumannii vs Klebsiella, P1_full) ===")
    base = results["%s__vs__%s__P1_full" % (AB.replace(" ", "_"), KL.replace(" ", "_"))]
    T0, fl0 = tables_for(cell, AB, KL, ELIG)
    lof = []
    for f in fl0:
        T, _ = tables_for(cell, AB, KL, ELIG, drop_fam=f)
        m = mh(T)
        if m:
            lof.append((f, abs(m["ln_or"] - base["ln_or"]) / abs(base["ln_or"])))
    lof.sort(key=lambda x: -x[1])
    print("  leave-one-family-out: %d families | max relative ln(OR) change %.4f (%.2f%%) on %s"
          % (len(lof), lof[0][1], 100 * lof[0][1], lof[0][0]))
    print("  families moving ln(OR) by more than %.0f%%: %d"
          % (100 * THR, sum(1 for _, v in lof if v > THR)))
    bps = sorted({bp for (f, h, bp) in cell if h in (AB, KL) and f in ELIG})
    lob = []
    for bp in bps:
        T, _ = tables_for(cell, AB, KL, ELIG, drop_bp=bp)
        m = mh(T)
        if m:
            lob.append((bp, abs(m["ln_or"] - base["ln_or"]) / abs(base["ln_or"])))
    lob.sort(key=lambda x: -x[1])
    print("  leave-one-BioProject-out: %d | max relative ln(OR) change %.4f (%.2f%%) on %s"
          % (len(lob), lob[0][1], 100 * lob[0][1], lob[0][0]))
    print("  BioProjects moving ln(OR) by more than %.0f%%: %d"
          % (100 * THR, sum(1 for _, v in lob if v > THR)))

    # ---------- cluster-level permutation of host labels within family ----------
    print("\n=== GLOBAL PERMUTATION (host labels permuted over BioProjects within family) ===")
    rng = np.random.default_rng(SEED)
    units = collections.defaultdict(list)
    for (f, h, bp), c in cell.items():
        if h in (AB, KL) and f in ELIG:
            units[f].append((h, c))
    obs = base["ln_or"]; ge = 0; done = 0
    for _ in range(a.perm):
        T = []
        for f, lst in units.items():
            labs = [h for h, _ in lst]
            rng.shuffle(labs)
            c1 = collections.Counter(); c2 = collections.Counter()
            for (_, c), nh in zip(lst, labs):
                (c1 if nh == AB else c2).update(c)
            a_ = c1["B"]; b_ = c1["C"] + c1["D"] + c1["E"]
            c_ = c2["B"]; d_ = c2["C"] + c2["D"] + c2["E"]
            if (a_ + b_) and (c_ + d_):
                T.append((a_, b_, c_, d_))
        m = mh(T)
        if m:
            done += 1
            if abs(m["ln_or"]) >= abs(obs):
                ge += 1
    pperm = (ge + 1) / (done + 1)
    print("  permutations completed: %d | |ln OR| >= observed: %d | empirical p = %.5f"
          % (done, ge, pperm))

    # ---------- secondary classifications ----------
    S = F["secondary_classifications"]
    bp_share = collections.defaultdict(collections.Counter)
    fam_tot = collections.Counter()
    for (f, h, bp), c in cell.items():
        bp_share[f][bp] += sum(c.values()); fam_tot[f] += sum(c.values())
    sec = {}
    for f in sorted(ELIG):
        top = max(bp_share[f].values()) / fam_tot[f] if fam_tot[f] else 0
        if top >= 0.70:
            sec[f] = "lineage_or_project_private"; continue
        bf = {}
        for h in (AB, PA, KL):
            if (f, h) in comp:
                c, n = comp[(f, h)]
                den = c["B"] + c["C"] + c["D"] + c["E"]
                if den:
                    bf[h] = c["B"] / den
        if len(bf) >= 2:
            mx = max(bf.values()) - min(bf.values())
            if mx >= 0.30:
                sec[f] = "host_dependent"; continue
            if mx < 0.15:
                sec[f] = "vehicle_stable"; continue
        nobs = sum(1 for k in "ABCDE" if any(comp.get((f, h), (collections.Counter(), 0))[0][k]
                                             for h in (AB, PA, KL)))
        sec[f] = "broad_route_portable" if nobs >= 4 else "unclassified"
    sc = collections.Counter(sec.values())
    print("\n=== SECONDARY CLASSIFICATIONS (frozen thresholds) ===")
    for k, v in sc.most_common():
        print("  %-30s %3d" % (k, v))

    # ---------- write TSVs ----------
    os.makedirs(a.outdir, exist_ok=True)
    p1 = os.path.join(a.outdir, "nmv4c_family_host_vehicle.tsv")
    cols = ["gene_family", "host1", "host2", "h1_classB", "h1_plasmid_CDE",
            "h2_classB", "h2_plasmid_CDE", "or", "ln_or", "var_ln_or", "mh_weight_share",
            "secondary_class", "tvd_ab_kl"]
    with open(p1, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in famrows:
            r["secondary_class"] = sec.get(r["gene_family"], "unclassified")
            r["tvd_ab_kl"] = round(tvd.get((r["gene_family"], AB, KL), float("nan")), 6)
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    p2 = os.path.join(a.outdir, "nmv4c_host_vehicle_summary.tsv")
    scols = ["contrast", "analysis", "n_families", "or", "ci_lo", "ci_hi", "excludes_1",
             "families_or_gt_1", "families_or_lt_1", "I2_pct", "max_family_weight_share"]
    with open(p2, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(scols) + "\n")
        for k, v in results.items():
            if not v:
                continue
            fh.write("\t".join(str(x) for x in [
                "%s vs %s" % (v["host1"], v["host2"]), v["analysis"], v["n_families"],
                round(v["or"], 6), round(v["ci_lo"], 6), round(v["ci_hi"], 6),
                v["excludes_1"], v["families_or_gt_1"], v["families_or_lt_1"],
                round(v["I2_pct"], 3), round(v["max_family_weight_share"], 6)]) + "\n")

    # ---------- gates ----------
    abkl = results["%s__vs__%s__P1_full" % (AB.replace(" ", "_"), KL.replace(" ", "_"))]
    abkl2 = results["%s__vs__%s__P2_one_genome_per_bioproject" % (AB.replace(" ", "_"), KL.replace(" ", "_"))]
    abpa = results["%s__vs__%s__P1_full" % (AB.replace(" ", "_"), PA.replace(" ", "_"))]
    cAB = collections.Counter(); cKL = collections.Counter()
    for f in ELIG:
        if (f, AB) in comp: cAB.update(comp[(f, AB)][0])
        if (f, KL) in comp: cKL.update(comp[(f, KL)][0])
    fiveclass_dir = (cAB["B"] / sum(cAB.values())) > (cKL["B"] / sum(cKL.values()))
    G = {
      "G1_pooled_effect": bool(abkl["excludes_1"] and abkl["or"] > 1),
      "G2_multiple_families": bool(abkl["families_or_gt_1"] >= MINFAM
                                   and abkl["max_family_weight_share"] <= 0.30),
      "G3_bioproject_balancing": bool(abkl2 and abkl2["excludes_1"] and abkl2["or"] > 1),
      "G4_leave_one_family_out": bool(lof[0][1] <= THR),
      "G5_pa_control": bool(abpa and abpa["excludes_1"] and abpa["or"] > 1),
      "G6_representation_coherence": bool(fiveclass_dir and abkl["or"] > 1),
      "G7_bioproject_influence": bool(lob[0][1] <= THR)}
    print("\n=== GATES ===")
    for k, v in G.items():
        print("  %-34s %s" % (k, v))
    npass = sum(G.values())
    verdict = ("NMV4C_SUCCESS_HOST_VEHICLE_PRINCIPLE_SUPPORTED" if npass == 7 else
               "NMV4C_PARTIAL_FAMILY_SPECIFIC_HETEROGENEITY" if G["G1_pooled_effect"] else
               "NMV4C_FAIL_NO_WITHIN_FAMILY_HOST_EFFECT")
    print("  gates passed: %d of 7" % npass)
    print("  VERDICT: %s" % verdict)

    rec = {"builder": VERSION,
           "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "frozen_design_sha256": FROZEN_SHA,
           "primary_denominator_verified": 74349,
           "class_reconciliation": dict(recon),
           "denominator_flow": {"primary": 74349, "exclusions": dict(excl),
                                "retained": kept,
                                "retained_one_genome_per_bioproject":
                                    sum(sum(c.values()) for c in cell2.values())},
           "n_eligible_families": len(ELIG),
           "results": results,
           "five_class_by_host": {h: {k: (agg(cell, ELIG, h)[k]) for k in "ABCDE"}
                                  for h in (AB, PA, KL)},
           "influence": {"leave_one_family_out_max_rel": lof[0][1],
                         "leave_one_family_out_worst": lof[0][0],
                         "n_families_over_threshold": sum(1 for _, v in lof if v > THR),
                         "leave_one_bioproject_out_max_rel": lob[0][1],
                         "leave_one_bioproject_out_worst": lob[0][0],
                         "n_bioprojects_over_threshold": sum(1 for _, v in lob if v > THR),
                         "threshold_relative_log_or": THR},
           "permutation": {"n_requested": a.perm, "n_completed": done,
                           "seed": SEED, "empirical_p": pperm},
           "secondary_classification_counts": dict(sc),
           "gates": G, "gates_passed": npass, "verdict": verdict,
           "outputs": {"family_tsv": sha256_file(p1), "summary_tsv": sha256_file(p2)},
           "statements": ["No transfer, conjugation or HGT event was observed or is claimed.",
                          "No causal host control is claimed.",
                          "Class A was retained and reported separately, never discarded.",
                          "No PortabilityEvent or PlasmidCall artefact was read."]}
    rp = os.path.join(a.outdir, "NMV4C_RESULT_RECEIPT.json")
    json.dump(rec, open(rp, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  %s  sha256 %s" % (p1, sha256_file(p1)))
    print("  %s  sha256 %s" % (p2, sha256_file(p2)))
    print("  %s  sha256 %s" % (rp, sha256_file(rp)))


if __name__ == "__main__":
    main()
