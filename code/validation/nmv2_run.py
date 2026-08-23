"""NM-V2 -- lineage and sampling-structure robustness.

Tests whether the headline findings, and the Acinetobacter baumannii chromosomal-MGE result in
particular, survive collapsing the cohort onto independent sampling units.

Gates come from NM0_VALIDATION_PROTOCOL_V1.json and are not restated as choices here:
  the A. baumannii / K. pneumoniae chromosomal-MGE ratio stays at or above 1.5x, with a lower
  bootstrap CI bound above 1.0 under BioProject-clustered resampling AND after collapsing to
  one genome per BioProject, AND no single BioProject removal moves the ratio by more than
  15 per cent relative.

One requirement of the module cannot be met with the current artefacts and is reported as such
rather than approximated: leave-one-major-lineage-out for the CHROMOSOMAL outcome. This cohort
has no MLST, no SNP distance and no chromosomal lineage assignment. Plasmid MOB clusters exist
but do not partition chromosomes, so using them here would invent a key.
"""
import argparse, collections, csv, datetime, hashlib, json, os, sys
import numpy as np

VERSION = "nmv2_run_v1.0.0"
PROTOCOL_SHA = "c968fb6d16a528a64d064d6f8bbac745390804df4ad0897c09b12de84ca3fbff"
AB, KP = "Acinetobacter baumannii", "Klebsiella pneumoniae"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def digest_rank(s):
    return hashlib.sha256(s.encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    if sha256_file(a.protocol) != PROTOCOL_SHA:
        print("REFUSING: protocol digest mismatch"); sys.exit(1)
    P = json.load(open(a.protocol, encoding="utf-8"))
    G = P["modules"]["NM-V2"]["prespecified_interpretations"]
    SEED = 20260821
    print("%s | protocol verified %s | seed %d" % (VERSION, PROTOCOL_SHA[:16], SEED))

    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    gen = {}
    for r in csv.DictReader(open(os.path.join(O, "genome_level_summary.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        gen[r["assembly_version"]] = {
            "sp": r["organism"], "bp": r["bioproject_accession"],
            "bs": r["biosample_accession"], "occ": int(r["n_arg_occurrences"]),
            "pl": int(r["n_plasmid_args"]), "ch": int(r["n_chromosomal_args"]),
            "nrep": int(r["n_arg_bearing_replicons"]), "npl": int(r["n_arg_bearing_plasmids"])}
    rep2asm = {}
    for r in csv.DictReader(open(os.path.join(O, "replicon_level_summary.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        if r["replicon_molecule_type"].lower().startswith("chrom"):
            rep2asm[r["replicon_accession"]] = r["assembly_version"]
    pos = set(r["block_id"] for r in csv.DictReader(
        open(os.path.join(O, "mge_feature_inventory.tsv"), encoding="utf-8"), delimiter="\t"))
    blk = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(O, "shared_context_blocks.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        asm = rep2asm.get(r["replicon_accession"])
        if asm:
            blk[asm].append(r["block_id"] in pos)
    # occurrence-level chromosomal MGE, for the second weighting
    occ_mge = collections.defaultdict(lambda: [0, 0])
    for r in csv.DictReader(open(os.path.join(O, "arg_mge_neighbourhood.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        v = occ_mge[r["assembly_version"]]
        v[0] += 1
        if int(r["n_mge_features"]) > 0:
            v[1] += 1

    by_sp = collections.defaultdict(list)
    for asm, g in gen.items():
        by_sp[g["sp"]].append(asm)

    def mge_block(asms):
        t = p_ = 0
        for x in asms:
            for b in blk.get(x, ()):
                t += 1; p_ += 1 if b else 0
        return (p_ / t if t else float("nan")), t

    def mge_occ(asms):
        t = p_ = 0
        for x in asms:
            v = occ_mge.get(x)
            if v:
                t += v[0]; p_ += v[1]
        return (p_ / t if t else float("nan")), t

    def ratio(asms_ab, asms_kp, fn):
        a1, n1 = fn(asms_ab); a2, n2 = fn(asms_kp)
        return (a1 / a2 if a2 else float("nan")), a1, a2, n1, n2

    print("\n=== BASELINE, full cohort ===")
    for name, fn in (("block-weighted", mge_block), ("occurrence-weighted", mge_occ)):
        r, m1, m2, n1, n2 = ratio(by_sp[AB], by_sp[KP], fn)
        print("  %-20s A.b %.4f (n=%d) | K.p %.4f (n=%d) | ratio %.3f"
              % (name, m1, n1, m2, n2, r))
    base_ratio = ratio(by_sp[AB], by_sp[KP], mge_block)[0]

    # ---------------- collapse rules ----------------
    def one_per(field, asms):
        best = {}
        for x in asms:
            k = gen[x][field]
            if k not in best or digest_rank(x) < digest_rank(best[k]):
                best[k] = x
        return sorted(best.values())

    print("\n=== COLLAPSE TO INDEPENDENT UNITS (deterministic: smallest SHA-256 of accession) ===")
    collapses = {}
    for label, field in (("full cohort", None), ("one genome per BioSample", "bs"),
                         ("one genome per BioProject", "bp")):
        A = by_sp[AB] if field is None else one_per(field, by_sp[AB])
        K = by_sp[KP] if field is None else one_per(field, by_sp[KP])
        r, m1, m2, n1, n2 = ratio(A, K, mge_block)
        collapses[label] = {"n_genomes_ab": len(A), "n_genomes_kp": len(K),
                            "ab_M": m1, "kp_M": m2, "ratio": r,
                            "ab_blocks": n1, "kp_blocks": n2}
        print("  %-28s genomes A.b %4d K.p %4d | A.b %.4f K.p %.4f | ratio %.3f"
              % (label, len(A), len(K), m1, m2, r))

    # ---------------- leave-one-BioProject-out ----------------
    print("\n=== LEAVE-ONE-BIOPROJECT-OUT ===")
    bps = sorted({gen[x]["bp"] for x in by_sp[AB] + by_sp[KP]})
    worst = None; rows = []
    for b in bps:
        A = [x for x in by_sp[AB] if gen[x]["bp"] != b]
        K = [x for x in by_sp[KP] if gen[x]["bp"] != b]
        if not A or not K:
            continue
        r = ratio(A, K, mge_block)[0]
        rel = abs(r - base_ratio) / base_ratio
        rows.append((b, r, rel))
        if worst is None or rel > worst[2]:
            worst = (b, r, rel)
    rels = np.array([x[2] for x in rows])
    print("  BioProjects tested: %d" % len(rows))
    print("  ratio range        : %.3f to %.3f (baseline %.3f)"
          % (min(x[1] for x in rows), max(x[1] for x in rows), base_ratio))
    print("  max relative change: %.4f (%.2f%%) on %s" % (worst[2], 100 * worst[2], worst[0]))
    print("  BioProjects moving the ratio by more than 15%%: %d"
          % int((rels > 0.15).sum()))

    # ---------------- BioProject-clustered bootstrap ----------------
    print("\n=== BIOPROJECT-CLUSTERED BOOTSTRAP, %d resamples ===" % a.boot)
    bp_by_sp = collections.defaultdict(lambda: collections.defaultdict(list))
    for x, g in gen.items():
        bp_by_sp[g["sp"]][g["bp"]].append(x)
    rng = np.random.default_rng(SEED)
    RB, RC = [], []
    for it in range(a.boot):
        out = {}
        for sp in (AB, KP):
            keys = list(bp_by_sp[sp])
            pick = rng.choice(len(keys), size=len(keys), replace=True)
            out[sp] = [x for j in pick for x in bp_by_sp[sp][keys[j]]]
        r = ratio(out[AB], out[KP], mge_block)[0]
        if r == r:
            RB.append(r)
        # same resample, collapsed to one genome per BioProject
        cA = one_per("bp", out[AB]); cK = one_per("bp", out[KP])
        r2 = ratio(cA, cK, mge_block)[0]
        if r2 == r2:
            RC.append(r2)
        if (it + 1) % 500 == 0:
            print("  %d/%d" % (it + 1, a.boot), flush=True)
    ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    RB = np.array(RB); RC = np.array(RC)
    b_lo, b_hi = ci(RB); c_lo, c_hi = ci(RC)
    print("\n  full-cohort ratio        median %.3f  95%% CI [%.3f, %.3f]"
          % (float(np.median(RB)), b_lo, b_hi))
    print("  one-genome-per-BioProject median %.3f  95%% CI [%.3f, %.3f]"
          % (float(np.median(RC)), c_lo, c_hi))

    # ---------------- four denominators, cohort-wide plasmid share ----------------
    print("\n=== FOUR DENOMINATORS, cohort-wide plasmid share ===")
    occ_t = sum(g["occ"] for g in gen.values()); occ_p = sum(g["pl"] for g in gen.values())
    rep_t = sum(g["nrep"] for g in gen.values()); rep_p = sum(g["npl"] for g in gen.values())
    gen_ev = sum((1 if g["pl"] else 0) + (1 if g["ch"] else 0) for g in gen.values())
    gen_pl = sum(1 for g in gen.values() if g["pl"])
    bpagg = collections.defaultdict(lambda: [0, 0])
    for g in gen.values():
        v = bpagg[g["bp"]]
        v[0] += 1 if g["pl"] else 0
        v[1] += 1 if g["ch"] else 0
    bp_ev = sum((1 if v[0] else 0) + (1 if v[1] else 0) for v in bpagg.values())
    bp_pl = sum(1 for v in bpagg.values() if v[0])
    dens = [("occurrence-weighted", occ_p, occ_t), ("replicon-weighted", rep_p, rep_t),
            ("genome-weighted", gen_pl, gen_ev), ("BioProject-weighted", bp_pl, bp_ev)]
    for lab, p_, t_ in dens:
        print("  %-22s %7d / %7d = %6.3f%%" % (lab, p_, t_, 100 * p_ / t_))

    # ---------------- balanced species weighting ----------------
    print("\n=== BALANCED SPECIES WEIGHTING (each species counts once) ===")
    sps = [s for s in by_sp if len(by_sp[s]) >= 40]
    pv = []
    for s in sps:
        t = sum(gen[x]["occ"] for x in by_sp[s]); p_ = sum(gen[x]["pl"] for x in by_sp[s])
        if t:
            pv.append(p_ / t)
    print("  species with >=40 genomes: %d" % len(sps))
    print("  occurrence-weighted cohort plasmid share : %.3f%%" % (100 * occ_p / occ_t))
    print("  unweighted mean of species shares        : %.3f%%" % (100 * float(np.mean(pv))))
    print("  median of species shares                 : %.3f%%" % (100 * float(np.median(pv))))

    # ---------------- gate ----------------
    ok_ratio = base_ratio >= 1.5
    ok_boot = b_lo > 1.0
    ok_collapse = collapses["one genome per BioProject"]["ratio"] >= 1.5 and c_lo > 1.0
    ok_loo = float(rels.max()) <= 0.15
    verdict = "SUCCESS" if (ok_ratio and ok_boot and ok_collapse and ok_loo) else "FAILURE"
    print("\n=== GATE ===")
    print("  baseline ratio >= 1.5x                            : %s (%.3f)" % (ok_ratio, base_ratio))
    print("  bootstrap lower CI > 1.0                          : %s (%.3f)" % (ok_boot, b_lo))
    print("  survives one-genome-per-BioProject with CI > 1.0  : %s (%.3f, CI lo %.3f)"
          % (ok_collapse, collapses["one genome per BioProject"]["ratio"], c_lo))
    print("  no single BioProject moves the ratio > 15%%        : %s (max %.2f%%)"
          % (ok_loo, 100 * float(rels.max())))
    print("  VERDICT: %s" % verdict)
    print("\n  leave-one-major-lineage-out (chromosomal): NOT EVALUABLE - this cohort has no")
    print("  MLST, no SNP distance and no chromosomal lineage assignment. Reported, not faked.")

    rec = {"builder": VERSION,
           "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "protocol_sha256": PROTOCOL_SHA, "seed": SEED,
           "baseline": {"block_weighted_ratio": base_ratio,
                        "ab_M_block": collapses["full cohort"]["ab_M"],
                        "kp_M_block": collapses["full cohort"]["kp_M"],
                        "occurrence_weighted": {k: v for k, v in zip(
                            ("ratio", "ab_M", "kp_M"), ratio(by_sp[AB], by_sp[KP], mge_occ)[:3])}},
           "collapses": collapses,
           "leave_one_bioproject_out": {
               "n_tested": len(rows), "baseline_ratio": base_ratio,
               "ratio_min": float(min(x[1] for x in rows)),
               "ratio_max": float(max(x[1] for x in rows)),
               "max_relative_change": float(rels.max()),
               "worst_bioproject": worst[0],
               "n_moving_more_than_15pct": int((rels > 0.15).sum())},
           "bootstrap": {"unit": "BioProject within species", "n_requested": a.boot,
                         "n_completed_full": int(len(RB)), "n_completed_collapsed": int(len(RC)),
                         "full_median": float(np.median(RB)), "full_ci": [b_lo, b_hi],
                         "collapsed_median": float(np.median(RC)), "collapsed_ci": [c_lo, c_hi]},
           "four_denominators": {lab: {"plasmid": p_, "total": t_, "pct": 100 * p_ / t_}
                                 for lab, p_, t_ in dens},
           "balanced_species_weighting": {
               "n_species_ge_40_genomes": len(sps),
               "occurrence_weighted_pct": 100 * occ_p / occ_t,
               "unweighted_species_mean_pct": 100 * float(np.mean(pv)),
               "species_median_pct": 100 * float(np.median(pv))},
           "not_evaluable": {
               "leave_one_major_lineage_out_chromosomal":
                 "NOT EVALUABLE - no MLST, no SNP distance, no chromosomal lineage assignment "
                 "exists in this cohort. Plasmid MOB clusters do not partition chromosomes and "
                 "using them here would invent a key. Requires owner decision D-NM4."},
           "gate": {"baseline_ratio_ge_1_5": bool(ok_ratio),
                    "bootstrap_lower_ci_gt_1": bool(ok_boot),
                    "survives_bioproject_collapse": bool(ok_collapse),
                    "no_single_bioproject_moves_gt_15pct": bool(ok_loo),
                    "verdict": verdict},
           "statements": ["No transfer, conjugation or HGT event was observed or is claimed.",
                          "No PortabilityEvent or PlasmidCall artefact was read.",
                          "The frozen PRIMARY denominator of 74,349 was not altered."]}
    json.dump(rec, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  receipt: %s\n  sha256 %s" % (a.out, sha256_file(a.out)))


if __name__ == "__main__":
    main()
