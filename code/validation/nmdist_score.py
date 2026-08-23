"""NM-DIST Phase 2 -- execution under the frozen protocol. No threshold or group is redefined."""
import argparse, collections, csv, datetime, hashlib, json, os, sys
import numpy as np

VERSION = "nmdist_score_v1.0.0"
PROTO_SHA = "226c0691cbd6ac9750fd8e816192420b3c0d40cadd4bb2ad05c7ca0089593a26"
L = 10000                      # window horizon / administrative censoring point
LAND = [1000, 2000, 5000, 10000]
B = 2000
SEED = 20260822


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def rd(p):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))


def sep(a1, a2, b1, b2):
    if not (b2 < a1 or b1 > a2):
        return 0
    return (a1 - b2) if b2 < a1 else (b1 - a2)


def circ_sep(a1, a2, b1, b2, Ln):
    best = sep(a1, a2, b1, b2)
    for s in (-Ln, Ln):
        best = min(best, sep(a1, a2, b1 + s, b2 + s))
    return max(0, best)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    D = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    O = a.outdir
    proto = os.path.join(a.repo, "docs/nature_microbiology/NMDIST_FROZEN_PROTOCOL_V1.json")
    if sha(proto) != PROTO_SHA:
        print("REFUSING: frozen protocol digest mismatch"); sys.exit(1)
    P = json.load(open(proto, encoding="utf-8"))
    for rel, meta in P["input_digests"].items():
        if sha(os.path.join(a.repo, rel)) != meta["sha256"]:
            print("REFUSING: input digest changed: %s" % rel); sys.exit(1)
    print("%s | frozen protocol and all %d input digests verified"
          % (VERSION, len(P["input_digests"])))

    win = rd(os.path.join(D, "arg_neighbourhood_windows.tsv"))
    mge = rd(os.path.join(D, "mge_feature_inventory.tsv"))
    blk = {r["block_id"]: r for r in rd(os.path.join(D, "shared_context_blocks.tsv"))}
    occ = {(r["assembly_version"], r["replicon_accession"], r["determinant_name"],
            r["gene_start"], r["gene_end"]): r
           for r in rd(os.path.join(D, "determinant_portability_classes.tsv"))
           if r["portability_class"] in ("A", "B")}
    nb = {(r["assembly_version"], r["replicon_accession"], r["determinant_name"],
           r["gene_start"], r["gene_end"]): r
          for r in rd(os.path.join(D, "arg_mge_neighbourhood.tsv"))}

    MB = collections.defaultdict(list)
    for r in mge:
        MB[r["block_id"]].append(r)

    # ---------- build the occurrence-level distance table ----------
    rows = []
    inblock_beyond_window = 0
    for r in win:
        k = (r["assembly_version"], r["replicon_accession"], r["determinant_name"],
             r["gene_start"], r["gene_end"])
        o = occ[k]
        b = r["block_id"]
        Ln = int(r["replicon_length"])
        circ = r["topology"].strip().lower().startswith("circ")
        gs, ge = int(r["gene_start"]), int(r["gene_end"])
        d_any = d_is = d_int = None
        any_in_block_beyond = False
        for f in MB.get(b, []):
            fs, fe = int(f["chrom_start"]), int(f["chrom_end"])
            d = circ_sep(gs, ge, fs, fe, Ln) if circ else sep(gs, ge, fs, fe)
            if d > L:
                any_in_block_beyond = True
                continue                       # strictly window-scoped, per the frozen protocol
            d_any = d if d_any is None else min(d_any, d)
            if f["feature_class"] == "IS_or_transposase":
                d_is = d if d_is is None else min(d_is, d)
            elif f["feature_class"] == "integrase_or_integron":
                d_int = d if d_int is None else min(d_int, d)
        if d_any is None and any_in_block_beyond:
            inblock_beyond_window += 1
        oid = "|".join(k)
        rows.append({
            "occurrence_id": oid, "block_id": b,
            "assembly_version": k[0], "replicon_accession": k[1],
            "bioproject_accession": o["bioproject_accession"],
            "organism_harmonized": o["organism_harmonized"], "genus": o["genus"],
            "determinant_name": k[2], "gene_family": o["gene_family"],
            "portability_class": o["portability_class"],
            "topology": r["topology"], "wrapped_circular": r["wrapped_circular"],
            "truncated": r["truncated"],
            "dist_any_bp": "" if d_any is None else d_any,
            "dist_is_bp": "" if d_is is None else d_is,
            "dist_integrase_bp": "" if d_int is None else d_int,
            "censored_any": int(d_any is None), "censored_is": int(d_is is None),
            "censored_integrase": int(d_int is None),
            "overlap_zero": int(d_any == 0),
            "in_block_feature_beyond_window": int(any_in_block_beyond and d_any is None)})

    # ---------- stop conditions ----------
    nblk = collections.Counter(r["block_id"] for r in rows)
    stop = []
    if len(rows) != 35140: stop.append("occurrences != 35,140")
    if len(nblk) != 21955: stop.append("blocks != 21,955")
    pos = sum(1 for r in rows if r["censored_any"] == 0)
    if pos != 16303: stop.append("window-positive != 16,303 (got %d)" % pos)
    clsB = sum(1 for r in rows if r["portability_class"] == "B")
    if pos != clsB: stop.append("window-positive != class B")
    if sum(1 for b in blk.values() if b["truncated"].strip().lower() in ("yes", "true", "1")) != 5:
        stop.append("truncated != 5")
    if inblock_beyond_window != 111:
        stop.append("in-block-beyond-window != 111 (got %d)" % inblock_beyond_window)
    # frozen distance column agreement
    mism = 0
    for r in rows:
        k = tuple(r["occurrence_id"].split("|"))
        fz = nb[k]["nearest_mge_distance_bp"].strip()
        got = r["dist_any_bp"]
        if fz in ("", "NA", "nan"):
            if got != "": mism += 1
        elif got == "" or int(float(fz)) != got:
            mism += 1
    if mism: stop.append("frozen distance disagreement: %d" % mism)
    if stop:
        print("*** STOP CONDITIONS TRIGGERED ***"); [print("   ", s) for s in stop]; sys.exit(9)
    print("  stop conditions: all clear | window-positive %d = class B | censored %d | "
          "in-block-beyond-window %d | frozen-column mismatches %d"
          % (pos, len(rows) - pos, inblock_beyond_window, mism))

    # ---------- weights ----------
    for r in rows:
        r["block_n_args"] = nblk[r["block_id"]]
        r["weight_block_balanced"] = 1.0 / nblk[r["block_id"]]
    wsum = sum(r["weight_block_balanced"] for r in rows)
    if abs(wsum - len(nblk)) > 1e-6:
        print("*** STOP: block-balanced weights sum to %.6f, expected %d" % (wsum, len(nblk)))
        sys.exit(9)
    print("  block-balanced weights sum to %.6f = %d blocks  PASS" % (wsum, len(nblk)))

    # blocks nest within BioProjects -> weights are invariant under BioProject resampling
    bpb = collections.defaultdict(set)
    for r in rows:
        bpb[r["block_id"]].add(r["bioproject_accession"])
    nest = all(len(v) == 1 for v in bpb.values())
    print("  every block nests within exactly one BioProject: %s "
          "(so recomputing weights inside a resample returns the same value)" % nest)
    if not nest:
        print("*** STOP: blocks do not nest within BioProjects"); sys.exit(9)

    # ---------- S2 deterministic one-occurrence-per-block ----------
    pick = {}
    for r in rows:
        h = hashlib.sha256(("%s|%s|20260822" % (r["block_id"], r["occurrence_id"]))
                           .encode()).hexdigest()
        if r["block_id"] not in pick or h < pick[r["block_id"]][0]:
            pick[r["block_id"]] = (h, r["occurrence_id"])
    s2 = {v[1] for v in pick.values()}
    for r in rows:
        r["s2_selected"] = int(r["occurrence_id"] in s2)
    print("  S2 deterministic selection: %d occurrences, one per block" % len(s2))

    # ---------- arrays ----------
    def grp(r):
        if r["organism_harmonized"] == "Acinetobacter baumannii": return "A. baumannii"
        if r["organism_harmonized"] == "Pseudomonas aeruginosa": return "P. aeruginosa"
        if r["genus"] == "Klebsiella": return "Klebsiella group"
        return "other"
    for r in rows:
        r["frozen_group"] = grp(r)
    G = ["A. baumannii", "P. aeruginosa", "Klebsiella group"]

    idx = {g: np.array([i for i, r in enumerate(rows) if r["frozen_group"] == g]) for g in G}
    def arr(field, cen):
        return np.array([L if r[cen] else int(r[field]) for r in rows], dtype=float)
    d_any = arr("dist_any_bp", "censored_any")
    d_is = arr("dist_is_bp", "censored_is")
    d_int = arr("dist_integrase_bp", "censored_integrase")
    obs_any = np.array([r["censored_any"] == 0 for r in rows])
    obs_is = np.array([r["censored_is"] == 0 for r in rows])
    obs_int = np.array([r["censored_integrase"] == 0 for r in rows])
    w_bb = np.array([r["weight_block_balanced"] for r in rows])
    w_occ = np.ones(len(rows))
    w_s2 = np.array([float(r["s2_selected"]) for r in rows])
    bps = np.array([r["bioproject_accession"] for r in rows])
    trunc = np.array([r["truncated"].strip().lower() in ("yes", "true", "1") for r in rows])
    wrapped = np.array([r["wrapped_circular"].strip().lower() in ("yes", "true", "1")
                        for r in rows])

    def summarise(d, obs, w, sel):
        """Weighted F(d) at landmarks, RMD through 10 kb, AUC, median."""
        if sel.sum() == 0 or w[sel].sum() == 0:
            return None
        dd, oo, ww = d[sel], obs[sel], w[sel]
        W = ww.sum()
        out = {"n_occurrences": int(sel.sum()), "weight_total": float(W),
               "n_observed": int(oo.sum()), "n_censored_at_10kb": int((~oo).sum())}
        for t in LAND:
            out["F_%d" % t] = float(ww[(dd <= t) & oo].sum() / W)
        rmd = float((ww * np.minimum(dd, L)).sum() / W)
        out["restricted_mean_distance_bp"] = rmd
        out["auc_detection_bp"] = L - rmd
        order = np.argsort(dd)
        cw = np.cumsum(ww[order]) / W
        do = dd[order]; oo2 = oo[order]
        med = None
        for i in range(len(do)):
            if oo2[i] and cw[i] >= 0.5:
                med = float(do[i]); break
        out["median_distance_bp"] = med if med is not None else "not reached"
        return out

    # ---------- primary estimates ----------
    prim = {}
    for g in G:
        sel = np.zeros(len(rows), bool); sel[idx[g]] = True
        prim[g] = summarise(d_any, obs_any, w_bb, sel)
    print("\n  PRIMARY (block-balanced, nearest MGE of either type)")
    for g in G:
        p = prim[g]
        print("    %-16s n=%5d  F1k=%.4f F2k=%.4f F5k=%.4f F10k=%.4f  RMD=%7.1f  median=%s"
              % (g, p["n_occurrences"], p["F_1000"], p["F_2000"], p["F_5000"], p["F_10000"],
                 p["restricted_mean_distance_bp"], p["median_distance_bp"]))

    # ---------- BioProject cluster bootstrap ----------
    ubp, inv = np.unique(bps, return_inverse=True)
    bp_idx = [np.where(inv == i)[0] for i in range(len(ubp))]
    rng = np.random.default_rng(SEED)
    print("\n  bootstrap: %d BioProjects, %d resamples, seed %d" % (len(ubp), B, SEED))

    def boot(d, obs, w, groups):
        keys = ["F_%d" % t for t in LAND] + ["restricted_mean_distance_bp", "auc_detection_bp"]
        acc = {g: {k: np.empty(B) for k in keys} for g in groups}
        gid = np.array([{"A. baumannii": 0, "P. aeruginosa": 1,
                         "Klebsiella group": 2}.get(r["frozen_group"], 3) for r in rows])
        for bi in range(B):
            take = rng.integers(0, len(ubp), len(ubp))
            ii = np.concatenate([bp_idx[t] for t in take])
            dd, oo, ww, gg = d[ii], obs[ii], w[ii], gid[ii]
            for gi, g in enumerate(groups):
                m = gg == gi
                W = ww[m].sum()
                if W == 0:
                    for k in keys: acc[g][k][bi] = np.nan
                    continue
                for t in LAND:
                    acc[g]["F_%d" % t][bi] = ww[m & (dd <= t) & oo].sum() / W
                rmd = (ww[m] * np.minimum(dd[m], L)).sum() / W
                acc[g]["restricted_mean_distance_bp"][bi] = rmd
                acc[g]["auc_detection_bp"][bi] = L - rmd
        return acc

    acc = boot(d_any, obs_any, w_bb, G)
    def ci(v):
        v = v[~np.isnan(v)]
        return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    for g in G:
        for k in acc[g]:
            prim[g][k + "_ci95"] = ci(acc[g][k])

    # ---------- contrasts ----------
    CON = [("P1", "A. baumannii", "Klebsiella group", "primary"),
           ("P2", "A. baumannii", "P. aeruginosa", "primary"),
           ("P3", "P. aeruginosa", "Klebsiella group", "secondary")]
    contrasts = []
    for tag, ga, gb, kind in CON:
        for t in LAND:
            k = "F_%d" % t
            diff = prim[ga][k] - prim[gb][k]
            bd = acc[ga][k] - acc[gb][k]
            bd = bd[~np.isnan(bd)]
            lo, hi = np.percentile(bd, 2.5), np.percentile(bd, 97.5)
            p = 2 * min((bd <= 0).mean(), (bd >= 0).mean())
            p = min(1.0, max(p, 1.0 / len(bd)))
            contrasts.append({"contrast": tag, "kind": kind, "group_a": ga, "group_b": gb,
                              "landmark_bp": t, "estimand": "F(d) block-balanced",
                              "value_a": prim[ga][k], "value_b": prim[gb][k],
                              "difference": diff, "ci_lo": float(lo), "ci_hi": float(hi),
                              "p_bootstrap": float(p), "p_holm": ""})
        dr = prim[ga]["restricted_mean_distance_bp"] - prim[gb]["restricted_mean_distance_bp"]
        bd = acc[ga]["restricted_mean_distance_bp"] - acc[gb]["restricted_mean_distance_bp"]
        bd = bd[~np.isnan(bd)]
        contrasts.append({"contrast": tag, "kind": kind, "group_a": ga, "group_b": gb,
                          "landmark_bp": "RMD_0_10kb", "estimand": "restricted mean distance",
                          "value_a": prim[ga]["restricted_mean_distance_bp"],
                          "value_b": prim[gb]["restricted_mean_distance_bp"],
                          "difference": dr,
                          "ci_lo": float(np.percentile(bd, 2.5)),
                          "ci_hi": float(np.percentile(bd, 97.5)),
                          "p_bootstrap": float(min(1.0, max(
                              2 * min((bd <= 0).mean(), (bd >= 0).mean()), 1.0 / len(bd)))),
                          "p_holm": ""})
    # Holm across the two PRIMARY contrasts, separately at each landmark
    for t in LAND:
        cell = [c for c in contrasts if c["kind"] == "primary" and c["landmark_bp"] == t]
        order = sorted(cell, key=lambda c: c["p_bootstrap"])
        m = len(order)
        prev = 0.0
        for i, c in enumerate(order):
            adj = min(1.0, max(prev, (m - i) * c["p_bootstrap"]))
            c["p_holm"] = adj
            prev = adj
    print("\n  PRIMARY CONTRASTS (Holm-corrected across the two primary contrasts per landmark)")
    for c in contrasts:
        if c["kind"] != "primary": continue
        print("    %s %-14s vs %-16s %-11s diff %+.4f  CI [%+.4f, %+.4f]  p=%.4g holm=%s"
              % (c["contrast"], c["group_a"], c["group_b"], str(c["landmark_bp"]),
                 c["difference"], c["ci_lo"], c["ci_hi"], c["p_bootstrap"],
                 ("%.4g" % c["p_holm"]) if c["p_holm"] != "" else "-"))

    # ---------- sensitivities ----------
    sens = []
    def add(tag, desc, d, obs, w, mask=None):
        for g in G:
            sel = np.zeros(len(rows), bool); sel[idx[g]] = True
            if mask is not None: sel &= mask
            s = summarise(d, obs, w, sel)
            if s is None: continue
            sens.append({"sensitivity": tag, "description": desc, "group": g, **{
                k: v for k, v in s.items()}})
    add("PRIMARY", "block-balanced, nearest MGE either type", d_any, obs_any, w_bb)
    add("S1", "occurrence-weighted, every occurrence weight 1", d_any, obs_any, w_occ)
    add("S2", "one deterministic occurrence per block", d_any, obs_any, w_s2)
    keep = np.zeros(len(rows), bool)
    bmin = {}
    for i, r in enumerate(rows):
        b = r["block_id"]
        v = d_any[i] if obs_any[i] else np.inf
        if b not in bmin or v < bmin[b][1]: bmin[b] = (i, v)
    for i, _ in bmin.values(): keep[i] = True
    add("S3", "minimum ARG-to-MGE distance per block (any-ARG block distance; FAVOURS multi-ARG "
        "blocks)", d_any, obs_any, w_occ, keep)
    add("S4", "excluding the five truncated blocks", d_any, obs_any, w_bb, ~trunc)
    add("S5", "excluding circular-wrapped blocks", d_any, obs_any, w_bb, ~wrapped)
    add("S6", "IS/transposase-only endpoint", d_is, obs_is, w_bb)
    add("SEC_INT", "integrase/integron-only endpoint (SECONDARY: sparse evidence arm)",
        d_int, obs_int, w_bb)
    print("\n  SENSITIVITIES computed: %s" % sorted({s["sensitivity"] for s in sens}))

    # ---------- leave-one-BioProject-out influence on P1 at 1 kb ----------
    selA = np.zeros(len(rows), bool); selA[idx["A. baumannii"]] = True
    selK = np.zeros(len(rows), bool); selK[idx["Klebsiella group"]] = True
    def f1k(sel, mask):
        m = sel & mask
        W = w_bb[m].sum()
        return float(w_bb[m & (d_any <= 1000) & obs_any].sum() / W) if W > 0 else np.nan
    base = f1k(selA, np.ones(len(rows), bool)) - f1k(selK, np.ones(len(rows), bool))
    infl = []
    for i, bp in enumerate(ubp):
        m = inv != i
        v = f1k(selA, m) - f1k(selK, m)
        if not np.isnan(v): infl.append((abs(v - base), bp, v))
    infl.sort(reverse=True)
    print("  LOBO on P1 at 1 kb: baseline %+.4f | max |change| %+.4f (%s) | n tested %d"
          % (base, infl[0][0], infl[0][1], len(infl)))

    # ---------- outputs ----------
    os.makedirs(O, exist_ok=True)
    def wtsv(name, data, fields):
        p = os.path.join(O, name)
        if os.path.exists(p): print("REFUSING: %s exists" % name); sys.exit(1)
        with open(p, "w", encoding="utf-8", newline="\n") as fh:
            w_ = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
            w_.writeheader(); w_.writerows(data)
        print("  wrote %-46s %s" % (name, sha(p)))
        return sha(p)

    occ_fields = ["occurrence_id", "block_id", "assembly_version", "replicon_accession",
                  "bioproject_accession", "organism_harmonized", "genus", "frozen_group",
                  "determinant_name", "gene_family", "portability_class", "topology",
                  "wrapped_circular", "truncated", "block_n_args", "weight_block_balanced",
                  "s2_selected", "dist_any_bp", "dist_is_bp", "dist_integrase_bp",
                  "censored_any", "censored_is", "censored_integrase", "overlap_zero",
                  "in_block_feature_beyond_window"]
    h1 = wtsv("nmdist_occurrence_block_distances.tsv", rows, occ_fields)
    pe = [{"group": g, **{k: v for k, v in prim[g].items()}} for g in G]
    pf = ["group", "n_occurrences", "weight_total", "n_observed", "n_censored_at_10kb"] + \
         ["F_%d" % t for t in LAND] + ["F_%d_ci95" % t for t in LAND] + \
         ["restricted_mean_distance_bp", "restricted_mean_distance_bp_ci95",
          "auc_detection_bp", "auc_detection_bp_ci95", "median_distance_bp"]
    h2 = wtsv("nmdist_primary_estimates.tsv", pe, pf)
    h3 = wtsv("nmdist_species_contrasts.tsv", contrasts,
              ["contrast", "kind", "group_a", "group_b", "landmark_bp", "estimand",
               "value_a", "value_b", "difference", "ci_lo", "ci_hi", "p_bootstrap", "p_holm"])
    sf = ["sensitivity", "description", "group", "n_occurrences", "weight_total", "n_observed",
          "n_censored_at_10kb"] + ["F_%d" % t for t in LAND] + \
         ["restricted_mean_distance_bp", "auc_detection_bp", "median_distance_bp"]
    h4 = wtsv("nmdist_sensitivity_results.tsv", sens, sf)

    # full curves for the figure
    curves = {}
    for g in G:
        sel = np.zeros(len(rows), bool); sel[idx[g]] = True
        W = w_bb[sel].sum()
        gridd = np.arange(0, L + 1, 25)
        curves[g] = {"grid": gridd.tolist(),
                     "F": [float(w_bb[sel & (d_any <= t) & obs_any].sum() / W) for t in gridd]}
    cb = {g: {"lo": [], "hi": []} for g in G}
    gid = np.array([{"A. baumannii": 0, "P. aeruginosa": 1,
                     "Klebsiella group": 2}.get(r["frozen_group"], 3) for r in rows])
    gridd = np.arange(0, L + 1, 250)
    rng2 = np.random.default_rng(SEED)
    bootc = {g: np.empty((B, len(gridd))) for g in G}
    for bi in range(B):
        take = rng2.integers(0, len(ubp), len(ubp))
        ii = np.concatenate([bp_idx[t] for t in take])
        dd, oo, ww, gg = d_any[ii], obs_any[ii], w_bb[ii], gid[ii]
        for gi, g in enumerate(G):
            m = gg == gi
            W = ww[m].sum()
            bootc[g][bi] = [ww[m & (dd <= t) & oo].sum() / W for t in gridd] if W > 0 else np.nan
    for g in G:
        cb[g]["grid"] = gridd.tolist()
        cb[g]["lo"] = np.nanpercentile(bootc[g], 2.5, axis=0).tolist()
        cb[g]["hi"] = np.nanpercentile(bootc[g], 97.5, axis=0).tolist()
    json.dump({"curves": curves, "bands": cb},
              open(os.path.join(O, "nmdist_curves.json"), "w", encoding="utf-8"), indent=1)

    def hhi(c):
        n = sum(c.values()); return 1.0 / sum((v / n) ** 2 for v in c.values())
    R = {"receipt": "NMDIST_RESULT_RECEIPT", "version": "1.0.0", "scorer": VERSION,
         "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
         "frozen_protocol_sha256": PROTO_SHA,
         "population": {"occurrences": len(rows), "blocks": len(nblk),
                        "window_positive": pos, "censored_at_10kb": len(rows) - pos,
                        "equals_class_B": pos == clsB,
                        "in_block_feature_beyond_window_censored": inblock_beyond_window,
                        "overlap_zero": int(sum(r["overlap_zero"] for r in rows)),
                        "truncated_blocks": 5, "wrapped_blocks": 57},
         "weights": {"block_balanced_sum": wsum, "n_blocks": len(nblk),
                     "blocks_nest_within_bioprojects": bool(nest)},
         "bootstrap": {"n_bioprojects": int(len(ubp)), "effective_1_over_HHI":
                       hhi(collections.Counter(bps.tolist())), "B": B, "seed": SEED},
         "primary_estimates": prim, "contrasts": contrasts,
         "influence_LOBO_P1_at_1kb": {"baseline_difference": base,
                                      "max_abs_change": infl[0][0],
                                      "most_influential_bioproject": infl[0][1],
                                      "n_tested": len(infl)},
         "frozen_column_mismatches": mism,
         "output_digests": {"nmdist_occurrence_block_distances.tsv": h1,
                            "nmdist_primary_estimates.tsv": h2,
                            "nmdist_species_contrasts.tsv": h3,
                            "nmdist_sensitivity_results.tsv": h4},
         "prohibited_claims_observed": ["no intact transposon claim", "no demonstrated "
                                        "mobilization", "no HGT", "no transfer", "no phenotype",
                                        "no distance beyond 10 kb estimated"]}
    rp = os.path.join(O, "NMDIST_RESULT_RECEIPT_V1.json")
    if os.path.exists(rp): print("REFUSING: receipt exists"); sys.exit(1)
    json.dump(R, open(rp, "w", encoding="utf-8", newline="\n"), indent=2, default=str)
    print("  wrote %-46s %s" % ("NMDIST_RESULT_RECEIPT_V1.json", sha(rp)))


if __name__ == "__main__":
    main()
