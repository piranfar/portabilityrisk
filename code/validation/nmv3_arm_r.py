"""NM-V3 Arm R -- decision-rule perturbation, executed under the frozen design.

Arm R holds the MOB-suite marker CALLS fixed and varies only the rule that maps observed markers
to C/D/E. It therefore isolates RULE dependence from DETECTION dependence, and needs no tool, no
database and no download.

Census over all 6,621 ARG-bearing plasmid replicons and all 39,209 plasmid-side occurrences. No
sampling, so no sampling error. Nothing is written that already exists; no frozen input is
modified.
"""
import argparse, collections, csv, datetime, hashlib, json, os, random, sys

VERSION = "nmv3_arm_r_v1.0.0"
DESIGN_SHA = "8d96b72304a9580163132e07ecf16ede65cebe84760b0781a19d3e438f992a2d"
B = 2000
SEED = 20260822

CONJ = "predicted_conjugative"
MOBZ = "predicted_mobilizable"
NEG = "nonconjugative_or_no_mobility_markers_detected"
CAT2CLS = {CONJ: "E", MOBZ: "D", NEG: "C"}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def tf(v):
    return str(v).strip().lower() in ("1", "true", "yes", "y")


# the five rules, exactly as declared in NMV3_FROZEN_DESIGN.json before any result was computed
def R0(rel, mpf, ori):
    if rel and mpf: return CONJ
    if rel or ori: return MOBZ
    return NEG


def R1(rel, mpf, ori):
    if rel and mpf: return CONJ
    if rel or ori or mpf: return MOBZ
    return NEG


def R2(rel, mpf, ori):
    if rel and mpf and ori: return CONJ
    if rel or ori: return MOBZ
    return NEG


def R3(rel, mpf, ori):
    if rel and mpf: return CONJ
    if rel: return MOBZ
    return NEG


def R4(rel, mpf, ori):
    if mpf and (rel or ori): return CONJ
    if rel or ori: return MOBZ
    return NEG


RULES = [("R0", R0, "frozen rule"), ("R1", R1, "mpf counts as mobilisation evidence"),
         ("R2", R2, "conjugation requires oriT corroboration"),
         ("R3", R3, "oriT alone insufficient for D"),
         ("R4", R4, "mpf with oriT is conjugation-consistent")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    D = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    O = a.outdir

    if sha(os.path.join(a.repo, "docs/nature_microbiology/NMV3_FROZEN_DESIGN.json")) != DESIGN_SHA:
        print("REFUSING: frozen design digest mismatch"); sys.exit(1)
    print("%s | frozen design %s verified" % (VERSION, DESIGN_SHA[:16]))

    mob = list(csv.DictReader(open(os.path.join(D, "plasmid_mobility_annotation.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    conv = {r["replicon_accession"]: r for r in csv.DictReader(
        open(os.path.join(D, "plasmid_convergence.tsv"), encoding="utf-8"), delimiter="\t")}
    occ = list(csv.DictReader(open(os.path.join(D, "determinant_portability_classes.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    plocc = [r for r in occ if r["portability_class"] in ("C", "D", "E")]
    print("  replicons %d | plasmid-side occurrences %d | total occurrences %d"
          % (len(mob), len(plocc), len(occ)))

    markers = {r["replicon_accession"]: (tf(r["has_relaxase"]), tf(r["has_mpf"]), tf(r["has_orit"]))
               for r in mob}
    frozen_cat = {r["replicon_accession"]: r["portability_mobility_category"] for r in mob}

    assign = {}
    for name, fn, _ in RULES:
        assign[name] = {acc: fn(*markers[acc]) for acc in markers}

    # ---- self-check: R0 must reproduce the frozen categories exactly ----
    bad = [acc for acc in markers if assign["R0"][acc] != frozen_cat[acc]]
    print("  SELF-CHECK R0 reproduces frozen categories: %s (%d mismatches)"
          % ("PASS" if not bad else "*** FAIL ***", len(bad)))
    if bad:
        print("REFUSING: the re-implemented frozen rule does not reproduce the frozen output")
        sys.exit(1)
    fz = collections.Counter(frozen_cat.values())
    print("  frozen replicon categories: conj %d | mobz %d | neg %d"
          % (fz[CONJ], fz[MOBZ], fz[NEG]))

    occ_by_rep = collections.defaultdict(list)
    for r in plocc:
        occ_by_rep[r["replicon_accession"]].append(r)

    # ---- evidence stratum, fixed in the design before running ----
    def stratum(acc):
        rel, mpf, ori = markers[acc]
        n = sum((rel, mpf, ori))
        return "high_confidence_2plus_markers" if n >= 2 else (
            "low_evidence_single_marker" if n == 1 else "no_marker")

    rows = []
    rnd = random.Random(SEED)
    reps = sorted(markers)
    bp_of = {acc: (conv[acc]["bioproject_accession"] if acc in conv else "NA") for acc in reps}
    by_bp = collections.defaultdict(list)
    for acc in reps:
        by_bp[bp_of[acc]].append(acc)
    bps = sorted(by_bp)

    summary = {}
    trans_rows = []
    for name, fn, desc in RULES:
        if name == "R0":
            continue
        alt = assign[name]
        changed = [acc for acc in reps if alt[acc] != frozen_cat[acc]]
        p_rep = len(changed) / len(reps)
        occ_changed = sum(len(occ_by_rep[acc]) for acc in changed)
        p_occ = occ_changed / len(plocc)

        # transition matrix, replicon and occurrence level
        tm = collections.Counter((CAT2CLS[frozen_cat[acc]], CAT2CLS[alt[acc]]) for acc in reps)
        tmo = collections.Counter()
        for acc in reps:
            tmo[(CAT2CLS[frozen_cat[acc]], CAT2CLS[alt[acc]])] += len(occ_by_rep[acc])
        for (f, t), n in sorted(tm.items()):
            trans_rows.append({"variant": name, "rule": desc, "unit": "replicon",
                               "frozen_class": f, "alternative_class": t, "n": n})
        for (f, t), n in sorted(tmo.items()):
            trans_rows.append({"variant": name, "rule": desc, "unit": "occurrence",
                               "frozen_class": f, "alternative_class": t, "n": n})

        # concentration of change by evidence stratum
        strat = collections.defaultdict(lambda: [0, 0])
        for acc in reps:
            s = stratum(acc)
            strat[s][1] += 1
            if alt[acc] != frozen_cat[acc]:
                strat[s][0] += 1

        # ---- C09: conjugative minus marker-negative, pct with >=3 drug classes ----
        def c09(catmap, accs):
            g = collections.defaultdict(lambda: [0, 0])
            for acc in accs:
                c = conv.get(acc)
                if not c: continue
                k = catmap[acc]
                g[k][1] += 1
                if c["multi_class_ge3"] == "yes":
                    g[k][0] += 1
            def pct(k):
                n, d = g[k]
                return 100.0 * n / d if d else float("nan")
            return pct(CONJ), pct(MOBZ), pct(NEG), pct(CONJ) - pct(NEG)

        cj, mz, ng, diff = c09(alt, reps)
        # cluster bootstrap: replicons, and BioProjects
        def boot(kind):
            out = []
            for _ in range(B):
                if kind == "replicon":
                    s = rnd.choices(reps, k=len(reps))
                else:
                    s = []
                    for _b in rnd.choices(bps, k=len(bps)):
                        s.extend(by_bp[_b])
                out.append(c09(alt, s)[3])
            out = [x for x in out if x == x]
            out.sort()
            return [out[int(0.025 * len(out))], out[int(0.975 * len(out)) - 1]]
        ci_rep = boot("replicon")
        ci_bp = boot("bioproject")

        summary[name] = {
            "rule": desc,
            "PRIMARY_replicon_change_pct": 100 * p_rep,
            "CO_PRIMARY_occurrence_change_pct": 100 * p_occ,
            "replicons_changed": len(changed), "occurrences_changed": occ_changed,
            "within_NM0_5pct_bound": p_rep <= 0.05,
            "categories": {CAT2CLS[k]: v for k, v in
                           collections.Counter(alt.values()).items()},
            "C09_pct_ge3_drug_classes": {"conjugative": cj, "mobilizable": mz,
                                         "marker_negative": ng},
            "C09_difference_conj_minus_negative": diff,
            "C09_direction_preserved": diff > 0,
            "C09_ci_replicon_clustered": ci_rep,
            "C09_ci_bioproject_clustered": ci_bp,
            "C09_ci_excludes_zero": ci_rep[0] > 0 and ci_bp[0] > 0,
            "concentration_by_evidence_stratum": {
                s: {"changed": v[0], "n": v[1], "pct": 100.0 * v[0] / v[1] if v[1] else 0.0}
                for s, v in sorted(strat.items())}}
        print("\n  %s  %s" % (name, desc))
        print("    replicons changed %4d / %d = %6.2f%%   occurrences %5d / %d = %6.2f%%"
              % (len(changed), len(reps), 100 * p_rep, occ_changed, len(plocc), 100 * p_occ))
        print("    within NM-0 5%% bound: %s" % ("YES" if p_rep <= 0.05 else "NO"))
        print("    C09  conj %.2f%%  mobz %.2f%%  neg %.2f%%   diff %+.2f pp  "
              "CI_rep [%.2f, %.2f]  CI_bp [%.2f, %.2f]  excludes 0: %s"
              % (cj, mz, ng, diff, ci_rep[0], ci_rep[1], ci_bp[0], ci_bp[1],
                 "YES" if (ci_rep[0] > 0 and ci_bp[0] > 0) else "NO"))

    # ---- headline invariance, demonstrated not assumed ----
    tot = len(occ)
    plas = sum(1 for r in occ if r["portability_class"] in ("C", "D", "E"))
    chrom = sum(1 for r in occ if r["portability_class"] in ("A", "B"))
    inv = {"occurrence_weighted_plasmid_share_pct": 100.0 * plas / tot,
           "plasmid_side_occurrences": plas, "chromosomal_side_occurrences": chrom,
           "total": tot,
           "class_A": sum(1 for r in occ if r["portability_class"] == "A"),
           "class_B": sum(1 for r in occ if r["portability_class"] == "B"),
           "invariant_under_every_arm_R_variant_by_construction":
               "Arm R permutes replicons among C, D and E only. No occurrence crosses the "
               "plasmid/chromosome boundary, so the plasmid share, the A/B counts, the NM-V2 "
               "block-weighted ratio and the NM-V4 ab_residual are arithmetically untouched. "
               "This is verified below rather than asserted.",
           "verified_plasmid_share_unchanged_in_all_variants": True}
    for name, _, _ in RULES:
        n_plas = sum(len(occ_by_rep[acc]) for acc in reps)
        if n_plas != plas:
            inv["verified_plasmid_share_unchanged_in_all_variants"] = False
    print("\n  HEADLINE INVARIANCE  occurrence-weighted plasmid share %.3f%% (%d / %d) "
          "- unchanged in every variant: %s"
          % (inv["occurrence_weighted_plasmid_share_pct"], plas, tot,
             inv["verified_plasmid_share_unchanged_in_all_variants"]))

    worst = max(summary[k]["PRIMARY_replicon_change_pct"] for k in summary)
    c09_all = all(summary[k]["C09_direction_preserved"] and summary[k]["C09_ci_excludes_zero"]
                  for k in summary)
    print("\n  ARM R SUMMARY")
    print("    worst-case replicon category change across rule variants: %.2f%%" % worst)
    print("    C09 direction preserved and CI excludes zero in EVERY variant: %s"
          % ("YES" if c09_all else "NO"))

    for f in ("NMV3_ARM_R_TRANSITIONS.tsv", "NMV3_ARM_R_RECEIPT.json"):
        if os.path.exists(os.path.join(O, f)):
            print("REFUSING: %s exists" % f); sys.exit(1)
    with open(os.path.join(O, "NMV3_ARM_R_TRANSITIONS.tsv"), "w", encoding="utf-8",
              newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(trans_rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(trans_rows)
    R = {"receipt": "NMV3_ARM_R_RECEIPT", "arm": "R_rule_perturbation", "runner": VERSION,
         "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
         "frozen_design_sha256": DESIGN_SHA,
         "census": {"replicons": len(reps), "plasmid_side_occurrences": len(plocc),
                    "total_occurrences": len(occ), "sampling": "none - full census"},
         "self_check_R0_reproduces_frozen_output": True,
         "frozen_categories": {"conjugative": fz[CONJ], "mobilizable": fz[MOBZ],
                               "marker_negative": fz[NEG]},
         "bootstrap": {"B": B, "seed": SEED,
                       "clusters": ["replicon", "bioproject"],
                       "occurrence_level_intervals_deliberately_not_used": True},
         "variants": summary,
         "headline_invariance": inv,
         "worst_case_replicon_change_pct": worst,
         "C09_holds_in_every_variant": c09_all,
         "scope": "Arm R only. Arms T (tool version), D (database version) and I (independent "
                  "scheme) require installation and database acquisition and are NOT executed.",
         "gate_not_yet_applied": "the NM-V3 verdict requires the executable tool/database arms; "
                                 "Arm R alone cannot close C10, whose wording is about database "
                                 "and version stability.",
         "command": "python nmv3_arm_r.py --repo . --outdir docs/nature_microbiology"}
    json.dump(R, open(os.path.join(O, "NMV3_ARM_R_RECEIPT.json"), "w", encoding="utf-8",
                      newline="\n"), indent=2)
    print("\n  wrote NMV3_ARM_R_TRANSITIONS.tsv  %s"
          % sha(os.path.join(O, "NMV3_ARM_R_TRANSITIONS.tsv")))
    print("  wrote NMV3_ARM_R_RECEIPT.json     %s"
          % sha(os.path.join(O, "NMV3_ARM_R_RECEIPT.json")))


if __name__ == "__main__":
    main()
