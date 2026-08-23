"""NMIS Phase 2/3 -- structural scoring under the frozen protocol and Amendment 001."""
import argparse, collections, csv, datetime, hashlib, json, os, sys
import numpy as np

VERSION = "nmis_score_v1.0.0"
PROTO = "5438045a3b73d123347fcd60b2456779f050c2a8fae9cd016782b8b6168b03a3"
AMEND = "72ce50c2cfd512c1bc6cb74d29830f70474ec5f1a12056a9e744883246e67f80"
L = 10000
LAND = [1000, 2000, 5000, 10000]
B = 2000
SEED = 20260822
G = ["A. baumannii", "P. aeruginosa", "Klebsiella group"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def rd(p):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))


def num(v):
    v = (v or "").strip()
    if v in ("", "-", "NA", "nan"):
        return None
    try:
        return int(float(v))
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--census", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    D = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    NM = os.path.join(a.repo, "docs/nature_microbiology")
    O = a.outdir

    if sha(os.path.join(NM, "NMIS_FROZEN_PROTOCOL_V1.json")) != PROTO:
        print("REFUSING: protocol digest mismatch"); sys.exit(1)
    if sha(os.path.join(NM, "NMIS_AMENDMENT_001_OCCURRENCE_WINDOW_CONTAINMENT.json")) != AMEND:
        print("REFUSING: amendment digest mismatch"); sys.exit(1)
    print("%s | frozen protocol and amendment 001 verified" % VERSION)

    blk = {r["block_id"]: r for r in rd(os.path.join(D, "shared_context_blocks.tsv"))}
    win = rd(os.path.join(D, "arg_neighbourhood_windows.tsv"))
    nmd = {r["occurrence_id"]: r for r in
           rd(os.path.join(NM, "nmdist_occurrence_block_distances.tsv"))}
    status = {r["block_id"]: r for r in rd(os.path.join(a.census, "nmis_census_status.tsv"))}

    # ---------------- parse the census by emitted header ----------------
    ISD = os.path.join(a.census, "isescan")
    elems = collections.defaultdict(list)
    n_el = 0
    types = collections.Counter()
    for f in sorted(os.listdir(ISD)):
        if not f.endswith(".tsv"):
            continue
        b = f[:-4]
        with open(os.path.join(ISD, f), encoding="utf-8") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            if rdr.fieldnames is None or len(rdr.fieldnames) != 24:
                print("REFUSING: unexpected header in %s" % f); sys.exit(1)
            for r in rdr:
                n_el += 1
                types[r["type"]] += 1
                elems[b].append(r)
    print("  parsed %d elements over %d blocks | type %s" % (n_el, len(elems), dict(types)))

    # ---------------- structural predicates ----------------
    def complete_structural(e):
        if (e["type"] or "").strip() != "c":
            return False
        ob, oe, ol = num(e["orfBegin"]), num(e["orfEnd"]), num(e["orfLen"])
        if ob is None or oe is None or ol is None or ol <= 0:
            return False
        s1, e1, s2, e2, il = (num(e["start1"]), num(e["end1"]), num(e["start2"]),
                              num(e["end2"]), num(e["irLen"]))
        if None in (s1, e1, s2, e2, il) or il <= 0:
            return False
        return True

    def elem_intervals(b, e):
        Bk = blk[b]
        bs, RL = int(Bk["block_start"]), int(Bk["replicon_length"])
        s, t = num(e["isBegin"]), num(e["isEnd"])
        if s is None or t is None:
            return None
        aa, zz = bs + s - 1, bs + t - 1
        if Bk["wrapped_circular"].strip().lower() in ("yes", "true", "1"):
            aa = ((aa - 1) % RL) + 1
            zz = ((zz - 1) % RL) + 1
            if zz < aa:
                return [(aa, RL), (1, zz)]
        return [(aa, zz)]

    def win_intervals(w):
        ws, we, RL = int(w["window_start"]), int(w["window_end"]), int(w["replicon_length"])
        if w["wrapped_circular"].strip().lower() in ("yes", "true", "1") and we < ws:
            return [(ws, RL), (1, we)]
        return [(ws, we)]

    def contained(ei, wi):
        return all(any(x >= s and y <= t for s, t in wi) for x, y in ei)

    def gap(a1, a2, b1, b2):
        return max(0, max(a1 - b2, b1 - a2))

    def dist(gs, ge, ei, circ, RL):
        best = None
        for x, y in ei:
            cands = [gap(gs, ge, x, y)]
            if circ:
                cands += [gap(gs, ge, x - RL, y - RL), gap(gs, ge, x + RL, y + RL)]
            d = min(cands)
            best = d if best is None else min(best, d)
        return best

    # ---------------- per-occurrence assignment ----------------
    rows = []
    ecat = collections.Counter()
    seen_elem = set()
    for w in win:
        b = w["block_id"]
        oid = "|".join([w["assembly_version"], w["replicon_accession"], w["determinant_name"],
                        w["gene_start"], w["gene_end"]])
        nd = nmd.get(oid)
        if nd is None:
            print("REFUSING: occurrence missing from the NM-DIST table: %s" % oid); sys.exit(1)
        Bk = blk[b]
        RL = int(Bk["replicon_length"])
        circ = w["topology"].strip().lower().startswith("circ")
        blen = int(Bk["block_span_bp"])
        gs, ge = int(w["gene_start"]), int(w["gene_end"])
        wi = win_intervals(w)
        st = status[b]["status"]
        best = None
        n_comp_contained = n_comp_cross = n_partial = 0
        for e in elems.get(b, []):
            ei = elem_intervals(b, e)
            if ei is None:
                continue
            touch = (num(e["isBegin"]) or 0) <= 1 or (num(e["isEnd"]) or 0) >= blen
            comp = complete_structural(e) and not touch
            key = (b, e["isBegin"], e["isEnd"])
            if comp:
                if contained(ei, wi):
                    n_comp_contained += 1
                    d = dist(gs, ge, ei, circ, RL)
                    if d is not None and d <= L:
                        best = d if best is None else min(best, d)
                    if key not in seen_elem:
                        ecat["complete_and_fully_contained"] += 1
                else:
                    n_comp_cross += 1
                    if key not in seen_elem:
                        ecat["complete_in_shared_block_but_crosses_occurrence_window_boundary"] += 1
            else:
                n_partial += 1
                if key not in seen_elem:
                    ecat["partial_or_boundary_limited"] += 1
            seen_elem.add(key)
        if st == "TOOL_FAILURE":
            state = "tool_failure"
        elif best is not None:
            state = "complete_and_fully_contained"
        elif n_comp_cross > 0:
            state = "complete_in_shared_block_but_crosses_occurrence_window_boundary"
        elif n_partial > 0:
            state = "partial_or_boundary_limited"
        else:
            state = "no_structurally_resolved_IS"
        rows.append({"occurrence_id": oid, "block_id": b,
                     "bioproject_accession": nd["bioproject_accession"],
                     "organism_harmonized": nd["organism_harmonized"], "genus": nd["genus"],
                     "frozen_group": nd["frozen_group"],
                     "portability_class": nd["portability_class"],
                     "weight_block_balanced": nd["weight_block_balanced"],
                     "nmdist_dist_any_bp": nd["dist_any_bp"],
                     "nmdist_censored": nd["censored_any"],
                     "endpoint_state": state,
                     "nmis_dist_bp": "" if best is None else best,
                     "nmis_censored": int(best is None),
                     "n_complete_contained": n_comp_contained,
                     "n_complete_crossing": n_comp_cross,
                     "n_partial_or_boundary": n_partial})

    # ---------------- reconciliation ----------------
    stop = []
    if len(rows) != 35140: stop.append("occurrences != 35,140 (%d)" % len(rows))
    cls = collections.Counter(r["portability_class"] for r in rows)
    if cls["A"] != 18837: stop.append("class A != 18,837")
    if cls["B"] != 16303: stop.append("class B != 16,303")
    if len({r["block_id"] for r in rows}) != 21955: stop.append("blocks != 21,955")
    st_c = collections.Counter(r["endpoint_state"] for r in rows)
    if sum(st_c.values()) != len(rows): stop.append("occurrences not all assigned one state")
    if stop:
        print("*** STOP ***"); [print("   ", s) for s in stop]; sys.exit(9)
    print("  reconciliation: 35,140 occurrences | A %d | B %d | 21,955 blocks | states %s"
          % (cls["A"], cls["B"], dict(st_c)))
    print("  element categories: %s" % dict(ecat))

    # ---------------- estimates ----------------
    idx = {g: np.array([i for i, r in enumerate(rows) if r["frozen_group"] == g]) for g in G}
    d = np.array([L if r["nmis_censored"] else int(r["nmis_dist_bp"]) for r in rows], float)
    obs = np.array([r["nmis_censored"] == 0 for r in rows])
    w_bb = np.array([float(r["weight_block_balanced"]) for r in rows])
    w_occ = np.ones(len(rows))
    bps = np.array([r["bioproject_accession"] for r in rows])
    dh = np.array([L if r["nmdist_censored"] == "1" else int(r["nmdist_dist_any_bp"])
                   for r in rows], float)
    oh = np.array([r["nmdist_censored"] == "0" for r in rows])

    def summ(dd, oo, ww, sel):
        m = sel
        W = ww[m].sum()
        if W == 0: return None
        o = {"n_occurrences": int(m.sum()), "n_observed": int(oo[m].sum()),
             "n_censored_at_10kb": int((~oo[m]).sum()), "weight_total": float(W)}
        for t in LAND:
            o["F_%d" % t] = float(ww[m & (dd <= t) & oo].sum() / W)
        rmd = float((ww[m] * np.minimum(dd[m], L)).sum() / W)
        o["restricted_mean_distance_bp"] = rmd
        o["auc_detection_bp"] = L - rmd
        dm, om, wm = dd[m], oo[m], ww[m]
        order = np.argsort(dm); cw = np.cumsum(wm[order]) / W
        med = "not reached"
        for i in range(len(order)):
            if om[order][i] and cw[i] >= 0.5:
                med = float(dm[order][i]); break
        o["median_distance_bp"] = med
        return o

    prim, homo = {}, {}
    for g in G:
        sel = np.zeros(len(rows), bool); sel[idx[g]] = True
        prim[g] = summ(d, obs, w_bb, sel)
        homo[g] = summ(dh, oh, w_bb, sel)
    print("\n  STRUCTURAL (block-balanced)")
    for g in G:
        p = prim[g]
        print("    %-16s F1k=%.4f F2k=%.4f F5k=%.4f F10k=%.4f RMD=%7.1f median=%s"
              % (g, p["F_1000"], p["F_2000"], p["F_5000"], p["F_10000"],
                 p["restricted_mean_distance_bp"], p["median_distance_bp"]))
    print("  HOMOLOGY (NM-DIST, same rows)")
    for g in G:
        p = homo[g]
        print("    %-16s F1k=%.4f F2k=%.4f F5k=%.4f F10k=%.4f RMD=%7.1f"
              % (g, p["F_1000"], p["F_2000"], p["F_5000"], p["F_10000"],
                 p["restricted_mean_distance_bp"]))

    # ---------------- bootstrap ----------------
    ubp, inv = np.unique(bps, return_inverse=True)
    bpi = [np.where(inv == i)[0] for i in range(len(ubp))]
    rng = np.random.default_rng(SEED)
    gid = np.array([G.index(r["frozen_group"]) if r["frozen_group"] in G else 3 for r in rows])
    keys = ["F_%d" % t for t in LAND] + ["restricted_mean_distance_bp"]
    acc = {g: {k: np.empty(B) for k in keys} for g in G}
    for bi in range(B):
        ii = np.concatenate([bpi[t] for t in rng.integers(0, len(ubp), len(ubp))])
        dd, oo, ww, gg = d[ii], obs[ii], w_bb[ii], gid[ii]
        for gi, g in enumerate(G):
            m = gg == gi; W = ww[m].sum()
            if W == 0:
                for k in keys: acc[g][k][bi] = np.nan
                continue
            for t in LAND:
                acc[g]["F_%d" % t][bi] = ww[m & (dd <= t) & oo].sum() / W
            acc[g]["restricted_mean_distance_bp"][bi] = (ww[m] * np.minimum(dd[m], L)).sum() / W
    for g in G:
        for k in keys:
            v = acc[g][k][~np.isnan(acc[g][k])]
            prim[g][k + "_ci95"] = [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]

    CON = [("N1", "A. baumannii", "Klebsiella group"), ("N2", "A. baumannii", "P. aeruginosa")]
    contrasts = []
    for tag, ga, gb in CON:
        for t in LAND:
            k = "F_%d" % t
            bd = acc[ga][k] - acc[gb][k]; bd = bd[~np.isnan(bd)]
            p = min(1.0, max(2 * min((bd <= 0).mean(), (bd >= 0).mean()), 1.0 / len(bd)))
            contrasts.append({"contrast": tag, "group_a": ga, "group_b": gb, "landmark_bp": t,
                              "value_a": prim[ga][k], "value_b": prim[gb][k],
                              "difference": prim[ga][k] - prim[gb][k],
                              "ci_lo": float(np.percentile(bd, 2.5)),
                              "ci_hi": float(np.percentile(bd, 97.5)),
                              "p_bootstrap": float(p), "p_holm": ""})
    for t in LAND:
        cell = sorted([c for c in contrasts if c["landmark_bp"] == t],
                      key=lambda c: c["p_bootstrap"])
        prev = 0.0
        for i, c in enumerate(cell):
            adj = min(1.0, max(prev, (len(cell) - i) * c["p_bootstrap"]))
            c["p_holm"] = adj; prev = adj
    print("\n  PRIMARY CONTRASTS (Holm across the two contrasts per landmark)")
    for c in contrasts:
        print("    %s %-14s vs %-16s %5d bp  diff %+.4f  CI [%+.4f, %+.4f]  p=%.4g holm=%.4g"
              % (c["contrast"], c["group_a"], c["group_b"], c["landmark_bp"], c["difference"],
                 c["ci_lo"], c["ci_hi"], c["p_bootstrap"], c["p_holm"]))

    # gate
    def ok(tag, t):
        c = [x for x in contrasts if x["contrast"] == tag and x["landmark_bp"] == t][0]
        return c["difference"] > 0 and c["ci_lo"] > 0 and c["p_holm"] <= 0.05
    succ = all(ok(tg, t) for tg in ("N1", "N2") for t in (1000, 2000))
    rev = all([x for x in contrasts if x["contrast"] == tg and x["landmark_bp"] == t][0]["difference"] > 0
              for tg in ("N1", "N2") for t in LAND)
    verdict = "SUCCESS" if succ else ("REVISE" if rev else "FAIL")
    print("\n  NMIS GATE VERDICT: %s" % verdict)

    # sensitivity: occurrence-weighted
    sens = []
    for tag, ww in (("PRIMARY_block_balanced", w_bb), ("S1_occurrence_weighted", w_occ)):
        for g in G:
            sel = np.zeros(len(rows), bool); sel[idx[g]] = True
            s = summ(d, obs, ww, sel)
            sens.append({"sensitivity": tag, "group": g, **s})

    # IS family composition among gated elements
    fam = collections.Counter(); famall = collections.Counter()
    seen2 = set()
    for b, es in elems.items():
        Bk = blk[b]; blen = int(Bk["block_span_bp"])
        for e in es:
            key = (b, e["isBegin"], e["isEnd"])
            if key in seen2: continue
            seen2.add(key)
            famall[e["family"]] += 1
            touch = (num(e["isBegin"]) or 0) <= 1 or (num(e["isEnd"]) or 0) >= blen
            if complete_structural(e) and not touch:
                fam[e["family"]] += 1
    print("\n  IS family, all elements     : %s" % dict(famall.most_common(6)))
    print("  IS family, structurally gated: %s" % dict(fam.most_common(6)))

    # corroboration
    corr = {}
    for g in G:
        sel = np.zeros(len(rows), bool); sel[idx[g]] = True
        hp = (sel & oh).sum(); sp = (sel & obs).sum()
        both = (sel & oh & obs).sum()
        corr[g] = {"homology_positive": int(hp), "structural_positive": int(sp),
                   "both": int(both),
                   "structural_share_of_homology": float(both / hp) if hp else None}
    print("  structural corroboration of homology-positive occurrences:")
    for g in G:
        c = corr[g]
        print("    %-16s homology %5d | structural %5d | both %5d | share %.4f"
              % (g, c["homology_positive"], c["structural_positive"], c["both"],
                 c["structural_share_of_homology"]))

    # ---------------- outputs ----------------
    os.makedirs(O, exist_ok=True)
    def wtsv(name, data, fields):
        p = os.path.join(O, name)
        if os.path.exists(p): print("REFUSING: %s exists" % name); sys.exit(1)
        with open(p, "w", encoding="utf-8", newline="\n") as fh:
            w_ = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
            w_.writeheader(); w_.writerows(data)
        print("  wrote %-44s %s" % (name, sha(p)))
        return sha(p)

    of = list(rows[0].keys())
    h1 = wtsv("nmis_occurrence_endpoints.tsv", rows, of)
    pf = ["group", "n_occurrences", "n_observed", "n_censored_at_10kb", "weight_total"] + \
         ["F_%d" % t for t in LAND] + ["F_%d_ci95" % t for t in LAND] + \
         ["restricted_mean_distance_bp", "restricted_mean_distance_bp_ci95",
          "auc_detection_bp", "median_distance_bp"]
    h2 = wtsv("nmis_primary_estimates.tsv", [{"group": g, **prim[g]} for g in G], pf)
    h3 = wtsv("nmis_species_contrasts.tsv", contrasts,
              ["contrast", "group_a", "group_b", "landmark_bp", "value_a", "value_b",
               "difference", "ci_lo", "ci_hi", "p_bootstrap", "p_holm"])
    h4 = wtsv("nmis_sensitivity_results.tsv", sens,
              ["sensitivity", "group", "n_occurrences", "n_observed", "n_censored_at_10kb"]
              + ["F_%d" % t for t in LAND]
              + ["restricted_mean_distance_bp", "median_distance_bp"])
    cmp_rows = []
    for g in G:
        for t in LAND:
            cmp_rows.append({"group": g, "landmark_bp": t,
                             "homology_F": homo[g]["F_%d" % t],
                             "structural_F": prim[g]["F_%d" % t],
                             "ratio_structural_over_homology":
                                 (prim[g]["F_%d" % t] / homo[g]["F_%d" % t])
                                 if homo[g]["F_%d" % t] else None})
    h5 = wtsv("nmis_vs_nmdist_comparison.tsv", cmp_rows,
              ["group", "landmark_bp", "homology_F", "structural_F",
               "ratio_structural_over_homology"])

    R = {"receipt": "NMIS_RESULT_RECEIPT", "version": "1.0.0", "scorer": VERSION,
         "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
         "frozen_protocol_sha256": PROTO, "amendment_001_sha256": AMEND,
         "census": {"blocks": 21955, "blocks_with_output": len(elems),
                    "elements_parsed": n_el, "type_counts": dict(types),
                    "tool_failures": 0},
         "element_categories": dict(ecat),
         "occurrence_endpoint_states": dict(st_c),
         "reconciliation": {"occurrences": len(rows), "class_A": cls["A"], "class_B": cls["B"],
                            "blocks": 21955, "class_B_unchanged": cls["B"] == 16303},
         "primary_estimates_structural": prim,
         "homology_estimates_same_rows": homo,
         "contrasts": contrasts, "gate_verdict": verdict,
         "sensitivity": sens,
         "is_family_all": dict(famall.most_common(12)),
         "is_family_structurally_gated": dict(fam.most_common(12)),
         "structural_corroboration": corr,
         "bootstrap": {"B": B, "seed": SEED, "n_bioprojects": int(len(ubp)),
                       "clusters": "BioProject"},
         "output_digests": {"nmis_occurrence_endpoints.tsv": h1,
                            "nmis_primary_estimates.tsv": h2,
                            "nmis_species_contrasts.tsv": h3,
                            "nmis_sensitivity_results.tsv": h4,
                            "nmis_vs_nmdist_comparison.tsv": h5},
         "prohibited_claims_observed": ["no transposition claim", "no HGT", "no transfer",
                                        "no phenotype", "no distance beyond 10 kb",
                                        "absence stated as 'no structurally resolved IS detected "
                                        "under the frozen ISEScan definition'"]}
    rp = os.path.join(O, "NMIS_RESULT_RECEIPT_V1.json")
    if os.path.exists(rp): print("REFUSING: receipt exists"); sys.exit(1)
    json.dump(R, open(rp, "w", encoding="utf-8", newline="\n"), indent=2, default=str)
    print("  wrote %-44s %s" % ("NMIS_RESULT_RECEIPT_V1.json", sha(rp)))


if __name__ == "__main__":
    main()
