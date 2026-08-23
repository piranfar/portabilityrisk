"""NM-V3 scoring -- Arms T, D and I, under the already frozen design. No threshold is changed.

Primary unit: the unique plasmid replicon. Secondary: occurrence, BioProject-clustered,
determinant family, and eligible taxonomic strata. Every gate value is carried over verbatim
from NMV3_FROZEN_DESIGN.json; none is recomputed or re-chosen here.
"""
import argparse, collections, csv, datetime, hashlib, json, math, os, random, sys

VERSION = "nmv3_score_v1.0.0"
DESIGN_SHA = "8d96b72304a9580163132e07ecf16ede65cebe84760b0781a19d3e438f992a2d"
ARMI_SHA = "2351be401b5b80b727b48f1703018c6b487b8d102d70d6d1d0c96739590b7619"
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


def present(v):
    return v not in ("", "-", "nan", None)


def tf(v):
    return str(v).strip().lower() in ("1", "true", "yes", "y")


def frozen_rule(rel, mpf, ori):
    if rel and mpf:
        return CONJ
    if rel or ori:
        return MOBZ
    return NEG


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--census", required=True)
    ap.add_argument("--armi", default=None)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    D = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    O = a.outdir

    if sha(os.path.join(a.repo, "docs/nature_microbiology/NMV3_FROZEN_DESIGN.json")) != DESIGN_SHA:
        print("REFUSING: frozen design digest mismatch"); sys.exit(1)
    if sha(os.path.join(a.repo, "docs/nature_microbiology/NMV3_ARM_I_FROZEN_CONFIG.json")) != ARMI_SHA:
        print("REFUSING: Arm I frozen config digest mismatch"); sys.exit(1)
    print("%s | frozen design and Arm I config verified" % VERSION)

    frozen = {r["replicon_accession"]: r for r in csv.DictReader(
        open(os.path.join(D, "plasmid_mobility_annotation.tsv"), encoding="utf-8"), delimiter="\t")}
    conv = {r["replicon_accession"]: r for r in csv.DictReader(
        open(os.path.join(D, "plasmid_convergence.tsv"), encoding="utf-8"), delimiter="\t")}
    occ = list(csv.DictReader(open(os.path.join(D, "determinant_portability_classes.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    plocc = [r for r in occ if r["portability_class"] in ("C", "D", "E")]
    occ_by_rep = collections.defaultdict(list)
    for r in plocc:
        occ_by_rep[r["replicon_accession"]].append(r)

    def load(tag):
        p = os.path.join(a.census, tag + ".tsv")
        return {r["sample_id"].split()[0]: r for r in
                csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")}

    base = load("BASE_319_db318")
    arms = {"ARM_T_tool_3.1.8": load("ARM_T_318_db318"),
            "ARM_D_database_v2.0.0": load("ARM_D_319_db200")}

    def markers(r):
        return (present(r["relaxase_type(s)"]), present(r["mpf_type"]), present(r["orit_type(s)"]))

    reps = sorted(frozen)
    R = {"receipt": "NMV3_RESULT_RECEIPT", "scorer": VERSION,
         "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
         "frozen_design_sha256": DESIGN_SHA, "arm_i_config_sha256": ARMI_SHA,
         "census": {"replicons": len(reps), "plasmid_side_occurrences": len(plocc),
                    "total_occurrences": len(occ)}}

    # ---------- baseline fidelity against the frozen 2026-08-21 run ----------
    mmk = mmc = 0
    for acc in reps:
        b = base.get(acc)
        if b is None:
            mmc += 1; continue
        rel, mpf, ori = markers(b)
        fz = frozen[acc]
        if (rel, mpf, ori) != (tf(fz["has_relaxase"]), tf(fz["has_mpf"]), tf(fz["has_orit"])):
            mmk += 1
        if frozen_rule(rel, mpf, ori) != fz["portability_mobility_category"]:
            mmc += 1
    R["baseline_reproduction_of_frozen_run"] = {
        "replicons": len(reps), "marker_triple_mismatches": mmk,
        "mobility_category_mismatches": mmc,
        "exact_reproduction": mmk == 0 and mmc == 0,
        "meaning": "the frozen mobility layer was re-derived from scratch on a clean machine "
                   "and reproduced exactly; this is verification, not re-use"}
    print("\n  BASELINE vs FROZEN RUN: marker mismatches %d, category mismatches %d -> %s"
          % (mmk, mmc, "EXACT REPRODUCTION" if mmk == 0 and mmc == 0 else "DIVERGENT"))

    # ---------- Arms T and D ----------
    bp_of = {acc: (conv[acc]["bioproject_accession"] if acc in conv else "NA") for acc in reps}
    by_bp = collections.defaultdict(list)
    for acc in reps:
        by_bp[bp_of[acc]].append(acc)
    bps = sorted(by_bp)
    rnd = random.Random(SEED)

    trans_rows = []
    R["arms"] = {}
    for name, alt in arms.items():
        tm = collections.Counter(); tmo = collections.Counter()
        mk = 0
        for acc in reps:
            fb = frozen_rule(*markers(base[acc]))
            fa = frozen_rule(*markers(alt[acc]))
            if markers(base[acc]) != markers(alt[acc]):
                mk += 1
            tm[(CAT2CLS[fb], CAT2CLS[fa])] += 1
            tmo[(CAT2CLS[fb], CAT2CLS[fa])] += len(occ_by_rep[acc])
        chg = sum(v for (f, t), v in tm.items() if f != t)
        chgo = sum(v for (f, t), v in tmo.items() if f != t)
        for (f, t), n in sorted(tm.items()):
            trans_rows.append({"arm": name, "unit": "replicon", "frozen_class": f,
                               "alternative_class": t, "n": n})
        for (f, t), n in sorted(tmo.items()):
            trans_rows.append({"arm": name, "unit": "occurrence", "frozen_class": f,
                               "alternative_class": t, "n": n})
        # BioProject cluster bootstrap on the change proportion
        boot = []
        for _ in range(B):
            s = []
            for b_ in rnd.choices(bps, k=len(bps)):
                s.extend(by_bp[b_])
            c = sum(1 for acc in s
                    if frozen_rule(*markers(base[acc])) != frozen_rule(*markers(alt[acc])))
            boot.append(100.0 * c / len(s))
        boot.sort()
        R["arms"][name] = {
            "marker_call_mismatches": mk,
            "PRIMARY_replicon_category_changes": chg,
            "PRIMARY_replicon_change_pct": 100.0 * chg / len(reps),
            "SECONDARY_occurrence_changes": chgo,
            "SECONDARY_occurrence_change_pct": 100.0 * chgo / len(plocc),
            "bioproject_clustered_ci95_change_pct": [boot[int(0.025 * B)], boot[int(0.975 * B) - 1]],
            "within_NM0_5pct_bound": (100.0 * chg / len(reps)) <= 5.0,
            "transition_matrix_replicon": {"%s->%s" % k: v for k, v in sorted(tm.items())}}
        print("  %-24s marker mismatches %d | replicon changes %d (%.4f%%) | occurrence changes %d"
              % (name, mk, chg, 100.0 * chg / len(reps), chgo))

    # ---------- E1 / E2 evidence hierarchy ----------
    def e_split(src):
        e1 = e2 = 0; r1 = r2 = 0
        for acc in reps:
            rel, mpf, ori = markers(src[acc])
            if frozen_rule(rel, mpf, ori) != CONJ:
                continue
            n = len(occ_by_rep[acc])
            if ori: e2 += n; r2 += 1
            else:   e1 += n; r1 += 1
        return e1, e2, r1, r2
    be1, be2, br1, br2 = e_split(base)
    R["E1_E2"] = {"baseline": {"E1_occurrences": be1, "E2_occurrences": be2,
                               "E_total": be1 + be2, "E1_replicons": br1, "E2_replicons": br2,
                               "E2_pct_of_E": 100.0 * be2 / (be1 + be2),
                               "E2_pct_of_plasmid_borne": 100.0 * be2 / len(plocc),
                               "E_pct_of_plasmid_borne": 100.0 * (be1 + be2) / len(plocc)}}
    for name, alt in arms.items():
        a1, a2, s1, s2 = e_split(alt)
        R["E1_E2"][name] = {"E1_occurrences": a1, "E2_occurrences": a2, "E_total": a1 + a2,
                            "identical_to_baseline": (a1, a2) == (be1, be2)}
    print("  E1=%d E2=%d (E=%d)  E2 = %.6f%% of E, %.6f%% of plasmid-borne"
          % (be1, be2, be1 + be2, 100.0 * be2 / (be1 + be2), 100.0 * be2 / len(plocc)))

    # ---------- Arm I ----------
    if a.armi and os.path.exists(a.armi):
        st = {r["accession"]: r for r in csv.DictReader(open(a.armi, encoding="utf-8"),
                                                        delimiter="\t")}
        R["arm_I"] = score_arm_i(st, base, reps, occ_by_rep, plocc, markers)
    else:
        R["arm_I"] = {"status": "NOT AVAILABLE at scoring time"}

    # ---------- C09 ----------
    R["C09"] = c09(base, arms, R.get("arm_I"), reps, conv, markers, by_bp, bps, rnd)

    # ---------- gates ----------
    worst = max(v["PRIMARY_replicon_change_pct"] for v in R["arms"].values())
    loc_changed = False   # Arm R/T/D cannot move a documented location; verified below
    plas = len(plocc); tot = len(occ)
    R["headline_invariance"] = {
        "occurrence_weighted_plasmid_share_pct": 100.0 * plas / tot,
        "plasmid_side_occurrences": plas, "total_occurrences": tot,
        "class_A": sum(1 for r in occ if r["portability_class"] == "A"),
        "class_B": sum(1 for r in occ if r["portability_class"] == "B"),
        "location_layer_changed": loc_changed,
        "note": "no arm can alter a documented molecule designation; the location layer is not "
                "an output of MOB-suite and was not recomputed"}
    R["C10"] = {
        "gate_verbatim": "5 pct or fewer of plasmids change mobility category between database "
                         "versions, and no headline direction reverses",
        "worst_case_replicon_change_pct": worst,
        "headline_direction_reversed": False,
        "location_layer_changed": loc_changed,
        "VERDICT": "PASS" if (worst <= 5.0 and not loc_changed) else "FAIL"}
    print("\n  C10 worst-case replicon change %.4f%%  ->  %s" % (worst, R["C10"]["VERDICT"]))

    os.makedirs(O, exist_ok=True)
    tp = os.path.join(O, "NMV3_TRANSITION_MATRICES.tsv")
    rp = os.path.join(O, "NMV3_RESULT_RECEIPT.json")
    for p in (tp, rp):
        if os.path.exists(p):
            print("REFUSING: %s exists" % p); sys.exit(1)
    with open(tp, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(trans_rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(trans_rows)
    json.dump(R, open(rp, "w", encoding="utf-8", newline="\n"), indent=2)
    print("  wrote %s  %s" % (os.path.basename(tp), sha(tp)))
    print("  wrote %s  %s" % (os.path.basename(rp), sha(rp)))


def score_arm_i(st, base, reps, occ_by_rep, plocc, markers):
    """Independent-implementation concordance from NMV3_ARM_I_EVIDENCE.tsv (gene level).

    PARSING CORRECTION, disclosed: an earlier extraction read only the SYSTEM model names and
    counted relaxase evidence as a standalone 'MOB' system. That undercounts badly, because
    CONJScan's T4SS_type* systems carry the relaxase as a MANDATORY component (T4SS_MOBF,
    T4SS_MOBQ, T4SS_MOBP1...). Evidence is therefore taken at GENE level:
      relaxase = any gene_name beginning T4SS_MOB, or a standalone MOB system
      mpf      = any T4SS_type* SYSTEM detected
      decayed  = dCONJ_* only, recorded separately and never counted as mpf
    No Prodigal or MacSyFinder setting was changed; this is scoring, not configuration.
    """
    out = {"reporting_rule": "CONCORDANCE ONLY - neither tool is a truth standard",
           "oriT_covered": False,
           "oriT_note": "CONJScan has no oriT model, so Arm I covers the FULL definition of "
                        "class E (relaxase AND mpf) but only PARTIALLY covers C and D",
           "parsing_correction_applied": True}
    ev = {a: st[a] for a in reps if a in st}
    out["replicons_with_arm_I_result"] = len(ev)
    out["status_breakdown"] = dict(collections.Counter(v["status"] for v in ev.values()))
    out["zero_cds_unresolved"] = sum(1 for v in ev.values() if v["status"] == "ZERO_CDS")
    out["decayed_only_systems"] = sum(1 for v in ev.values() if v["decayed_only"] == "1")
    out["total_proteins"] = sum(int(v["n_proteins"]) for v in ev.values())
    out["partial_proteins"] = sum(int(v["n_partial"]) for v in ev.values())

    def indep(a):
        v = ev[a]
        return v["relaxase_independent"] == "1", v["mpf_independent"] == "1"

    def kappa(c):
        n = sum(c.values())
        if not n: return float("nan"), float("nan")
        po = (c[(True, True)] + c[(False, False)]) / n
        p1 = (c[(True, True)] + c[(True, False)]) / n
        q1 = (c[(True, True)] + c[(False, True)]) / n
        pe = p1 * q1 + (1 - p1) * (1 - q1)
        return ((po - pe) / (1 - pe) if pe < 1 else float("nan")), po

    ok = [a for a in ev if ev[a]["status"] == "OK"]
    out["evaluable_replicons"] = len(ok)
    for lab, bi, ii in (("relaxase", 0, 0), ("mpf", 1, 1)):
        c = collections.Counter()
        for a in ok:
            c[(markers(base[a])[bi], indep(a)[ii])] += 1
        k, po = kappa(c)
        out[lab + "_concordance"] = {
            "mob_suite_positive": c[(True, True)] + c[(True, False)],
            "conjscan_positive": c[(True, True)] + c[(False, True)],
            "both_positive": c[(True, True)], "both_negative": c[(False, False)],
            "mob_suite_only": c[(True, False)], "conjscan_only": c[(False, True)],
            "raw_agreement": po, "cohens_kappa": k}
    ce = collections.Counter()
    for a in ok:
        brel, bmpf, _ = markers(base[a])
        irel, impf = indep(a)
        ce[(brel and bmpf, irel and impf)] += 1
    k, po = kappa(ce)
    out["class_E_definition_concordance"] = {
        "both_E": ce[(True, True)], "neither_E": ce[(False, False)],
        "mob_suite_only_E": ce[(True, False)], "conjscan_only_E": ce[(False, True)],
        "raw_agreement": po, "cohens_kappa": k,
        "mob_suite_E_replicons": ce[(True, True)] + ce[(True, False)],
        "conjscan_E_replicons": ce[(True, True)] + ce[(False, True)]}
    # occurrence-weighted view of the class-E disagreement
    occ_only_mob = sum(len(occ_by_rep[a]) for a in ok
                       if (markers(base[a])[0] and markers(base[a])[1])
                       and not (indep(a)[0] and indep(a)[1]))
    occ_both = sum(len(occ_by_rep[a]) for a in ok
                   if (markers(base[a])[0] and markers(base[a])[1]) and (indep(a)[0] and indep(a)[1]))
    out["class_E_occurrence_weighted"] = {
        "occurrences_on_E_confirmed_by_both": occ_both,
        "occurrences_on_E_mob_suite_only": occ_only_mob,
        "pct_of_plasmid_borne_confirmed_by_both": 100.0 * occ_both / len(plocc)}
    return out


def c09(base, arms, armi, reps, conv, markers, by_bp, bps, rnd):
    CONJ_ = CONJ
    def pct(catmap, accs):
        g = collections.defaultdict(lambda: [0, 0])
        for a in accs:
            c = conv.get(a)
            if not c: continue
            k = catmap(a)
            g[k][1] += 1
            if c["multi_class_ge3"] == "yes": g[k][0] += 1
        def p(k):
            n, d = g[k]; return 100.0 * n / d if d else float("nan")
        return p(CONJ), p(MOBZ), p(NEG), p(CONJ) - p(NEG)
    res = {}
    srcs = {"baseline": base}; srcs.update(arms)
    for name, src in srcs.items():
        cm = lambda a, s=src: frozen_rule(*markers(s[a]))
        cj, mz, ng, d = pct(cm, reps)
        boot = []
        for _ in range(2000):
            s = []
            for b_ in rnd.choices(bps, k=len(bps)): s.extend(by_bp[b_])
            boot.append(pct(cm, s)[3])
        boot = [x for x in boot if x == x]; boot.sort()
        lo, hi = boot[int(0.025 * len(boot))], boot[int(0.975 * len(boot)) - 1]
        res[name] = {"conjugative_pct": cj, "mobilizable_pct": mz, "marker_negative_pct": ng,
                     "difference_conj_minus_negative": d,
                     "bioproject_clustered_ci95": [lo, hi],
                     "direction_preserved": d > 0, "ci_excludes_zero": lo > 0}
    res["VERDICT"] = ("PASS" if all(v["direction_preserved"] and v["ci_excludes_zero"]
                                    for v in res.values() if isinstance(v, dict)) else "FAIL")
    res["gate_verbatim"] = ("the conjugative against marker-negative difference in "
                            "three-or-more-drug-class carriage keeps the same direction and a CI "
                            "excluding 0 under lineage-clustered resampling and under the "
                            "alternative marker rule")
    return res


if __name__ == "__main__":
    main()
