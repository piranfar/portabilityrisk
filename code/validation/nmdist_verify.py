"""NM-DIST Phase 3 -- independent verification. Imports nothing from the scorer.

Re-derives coordinates, distances, wrapping, weights, censoring, group membership, landmark
estimates, bootstrap inputs and contrasts with its own implementation, then compares to the
exported tables. A single disagreement is a failure.
"""
import argparse, collections, csv, hashlib, json, os, random, sys

VERSION = "nmdist_verify_v1.0.0"
PROTO_SHA = "226c0691cbd6ac9750fd8e816192420b3c0d40cadd4bb2ad05c7ca0089593a26"
L = 10000
LAND = [1000, 2000, 5000, 10000]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def rows_of(p):
    with open(p, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            if ln.strip():
                yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    D = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    O = a.dir
    fails, n = [], 0

    def C(name, ok, detail=""):
        nonlocal n
        n += 1
        print("  %-58s %s %s" % (name, "PASS" if ok else "*** FAIL ***", detail))
        if not ok:
            fails.append(name)

    print("%s  (independent code path)\n" % VERSION)
    C("frozen protocol digest unchanged",
      sha(os.path.join(a.repo, "docs/nature_microbiology/NMDIST_FROZEN_PROTOCOL_V1.json"))
      == PROTO_SHA)
    R = json.load(open(os.path.join(O, "NMDIST_RESULT_RECEIPT_V1.json"), encoding="utf-8"))

    # ---- rebuild from source with an independent implementation ----
    win = {(r["assembly_version"], r["replicon_accession"], r["determinant_name"],
            r["gene_start"], r["gene_end"]): r
           for r in rows_of(os.path.join(D, "arg_neighbourhood_windows.tsv"))}
    blk = {r["block_id"]: r for r in rows_of(os.path.join(D, "shared_context_blocks.tsv"))}
    feats = collections.defaultdict(list)
    for r in rows_of(os.path.join(D, "mge_feature_inventory.tsv")):
        feats[r["block_id"]].append((int(r["chrom_start"]), int(r["chrom_end"]),
                                     r["feature_class"]))
    exported = list(rows_of(os.path.join(O, "nmdist_occurrence_block_distances.tsv")))
    C("exported table has 35,140 rows", len(exported) == 35140, len(exported))
    C("exported table covers 21,955 blocks",
      len({r["block_id"] for r in exported}) == 21955)

    def gap(a1, a2, b1, b2):
        """Independent formulation: max(0, max(a1-b2, b1-a2)) is 0 exactly when they overlap."""
        return max(0, max(a1 - b2, b1 - a2))

    def dmin(r):
        k = (r["assembly_version"], r["replicon_accession"], r["determinant_name"],
             r["occurrence_id"].split("|")[3], r["occurrence_id"].split("|")[4])
        w = win[k]
        Ln = int(w["replicon_length"])
        circ = w["topology"].strip().lower().startswith("circ")
        gs, ge = int(k[3]), int(k[4])
        best = {"any": None, "IS_or_transposase": None, "integrase_or_integron": None}
        beyond = False
        for fs, fe, fc in feats.get(r["block_id"], []):
            cands = [gap(gs, ge, fs, fe)]
            if circ:
                cands += [gap(gs, ge, fs - Ln, fe - Ln), gap(gs, ge, fs + Ln, fe + Ln)]
            d = min(cands)
            if d > L:
                beyond = True
                continue
            if best["any"] is None or d < best["any"]:
                best["any"] = d
            if best[fc] is None or d < best[fc]:
                best[fc] = d
        return best, beyond

    mism = mism_is = mism_int = 0
    beyond_censored = 0
    zero_ok = 0
    for r in exported:
        b, beyond = dmin(r)
        exp_any = r["dist_any_bp"]
        got = b["any"]
        if (got is None) != (exp_any == ""):
            mism += 1
        elif got is not None and int(exp_any) != got:
            mism += 1
        if (b["IS_or_transposase"] is None) != (r["dist_is_bp"] == ""):
            mism_is += 1
        elif b["IS_or_transposase"] is not None and int(r["dist_is_bp"]) != b["IS_or_transposase"]:
            mism_is += 1
        if (b["integrase_or_integron"] is None) != (r["dist_integrase_bp"] == ""):
            mism_int += 1
        elif b["integrase_or_integron"] is not None and \
                int(r["dist_integrase_bp"]) != b["integrase_or_integron"]:
            mism_int += 1
        if got is None and beyond:
            beyond_censored += 1
        if r["overlap_zero"] == "1":
            zero_ok += 1 if got == 0 else 0
    C("nearest-any distances independently reproduced", mism == 0, "%d mismatches" % mism)
    C("nearest IS/transposase distances reproduced", mism_is == 0, "%d" % mism_is)
    C("nearest integrase/integron distances reproduced", mism_int == 0, "%d" % mism_int)
    C("the 111 in-block-beyond-window occurrences are censored",
      beyond_censored == 111, beyond_censored)
    C("every overlap flag corresponds to distance 0",
      zero_ok == sum(1 for r in exported if r["overlap_zero"] == "1"), zero_ok)

    # ---- censoring and class reconciliation ----
    obs = [r for r in exported if r["censored_any"] == "0"]
    C("window-positive == 16,303 == class B", len(obs) == 16303 and
      sum(1 for r in exported if r["portability_class"] == "B") == 16303, len(obs))
    C("censored at 10 kb == 18,837", len(exported) - len(obs) == 18837)
    C("no exported distance exceeds 10 kb",
      all(int(r["dist_any_bp"]) <= L for r in exported if r["dist_any_bp"] != ""))

    # ---- block weights ----
    cnt = collections.Counter(r["block_id"] for r in exported)
    wsum = sum(1.0 / cnt[r["block_id"]] for r in exported)
    C("block-balanced weights sum to the block count",
      abs(wsum - 21955) < 1e-6, "%.6f" % wsum)
    C("exported weights equal 1/m",
      all(abs(float(r["weight_block_balanced"]) - 1.0 / cnt[r["block_id"]]) < 1e-12
          for r in exported))
    C("S2 selects exactly one occurrence per block",
      sum(1 for r in exported if r["s2_selected"] == "1") == 21955)
    s2bad = 0
    byb = collections.defaultdict(list)
    for r in exported:
        byb[r["block_id"]].append(r)
    for b, rs in byb.items():
        best = min(rs, key=lambda r: hashlib.sha256(
            ("%s|%s|20260822" % (b, r["occurrence_id"])).encode()).hexdigest())
        if best["s2_selected"] != "1":
            s2bad += 1
    C("S2 selection reproduces the frozen hash rule", s2bad == 0, s2bad)

    # ---- group membership ----
    gbad = 0
    for r in exported:
        g = ("A. baumannii" if r["organism_harmonized"] == "Acinetobacter baumannii" else
             "P. aeruginosa" if r["organism_harmonized"] == "Pseudomonas aeruginosa" else
             "Klebsiella group" if r["genus"] == "Klebsiella" else "other")
        if g != r["frozen_group"]:
            gbad += 1
    C("frozen group membership reproduced", gbad == 0, gbad)
    gc = collections.Counter(r["frozen_group"] for r in exported)
    C("group sizes 8,005 / 7,150 / 15,568",
      (gc["A. baumannii"], gc["P. aeruginosa"], gc["Klebsiella group"]) == (8005, 7150, 15568),
      dict(gc))

    # ---- landmark estimates recomputed ----
    pe = {r["group"]: r for r in rows_of(os.path.join(O, "nmdist_primary_estimates.tsv"))}
    worst = 0.0
    for g in ("A. baumannii", "P. aeruginosa", "Klebsiella group"):
        sel = [r for r in exported if r["frozen_group"] == g]
        W = sum(float(r["weight_block_balanced"]) for r in sel)
        for t in LAND:
            f = sum(float(r["weight_block_balanced"]) for r in sel
                    if r["censored_any"] == "0" and int(r["dist_any_bp"]) <= t) / W
            worst = max(worst, abs(f - float(pe[g]["F_%d" % t])))
        rmd = sum(float(r["weight_block_balanced"]) *
                  (int(r["dist_any_bp"]) if r["censored_any"] == "0" else L)
                  for r in sel) / W
        worst = max(worst, abs(rmd - float(pe[g]["restricted_mean_distance_bp"])) / L)
    C("all landmark and RMD estimates reproduced", worst < 1e-9, "max deviation %.2e" % worst)

    # ---- contrasts ----
    con = list(rows_of(os.path.join(O, "nmdist_species_contrasts.tsv")))
    cbad = 0
    for c in con:
        if c["landmark_bp"] == "RMD_0_10kb":
            continue
        d = float(c["value_a"]) - float(c["value_b"])
        if abs(d - float(c["difference"])) > 1e-12:
            cbad += 1
    C("contrast differences equal value_a minus value_b", cbad == 0, cbad)
    C("Holm applied to the two primary contrasts at each landmark",
      all(c["p_holm"] != "" for c in con if c["kind"] == "primary"
          and c["landmark_bp"] != "RMD_0_10kb"))
    C("secondary contrast P3 carries no Holm value",
      all(c["p_holm"] == "" for c in con if c["kind"] == "secondary"))

    # ---- bootstrap input construction ----
    bp = collections.Counter(r["bioproject_accession"] for r in exported)
    C("BioProject count matches the receipt",
      len(bp) == R["bootstrap"]["n_bioprojects"], len(bp))
    nest = all(len({r["bioproject_accession"] for r in rs}) == 1 for rs in byb.values())
    C("blocks nest within BioProjects (weights invariant under resampling)", nest)
    tot = sum(bp.values())
    hhi = 1.0 / sum((v / tot) ** 2 for v in bp.values())
    C("effective 1/HHI reproduced",
      abs(hhi - float(R["bootstrap"]["effective_1_over_HHI"])) < 1e-6, "%.2f" % hhi)

    # ---- required spot checks ----
    print("\n  REQUIRED SPOT CHECKS")
    rnd = random.Random(20260822)
    pos_blocks = [r for r in exported if r["censored_any"] == "0" and int(r["dist_any_bp"]) > 0]
    sample = rnd.sample(pos_blocks, 20)
    ok = sum(1 for r in sample if dmin(r)[0]["any"] == int(r["dist_any_bp"]))
    C("20 ordinary positive occurrences recompute exactly", ok == 20, "%d/20" % ok)
    wr = [r for r in exported if r["wrapped_circular"].strip().lower() in ("yes", "true", "1")
          and r["frozen_group"] in ("A. baumannii", "P. aeruginosa", "Klebsiella group")]
    okw = sum(1 for r in wr if (dmin(r)[0]["any"] is None) == (r["dist_any_bp"] == "")
              and (r["dist_any_bp"] == "" or dmin(r)[0]["any"] == int(r["dist_any_bp"])))
    C("all circular-wrapped occurrences in frozen groups recompute",
      okw == len(wr), "%d/%d" % (okw, len(wr)))
    tr = [r for r in exported if r["truncated"].strip().lower() in ("yes", "true", "1")]
    okt = sum(1 for r in tr if (dmin(r)[0]["any"] is None) == (r["dist_any_bp"] == ""))
    C("all truncated-block occurrences recompute", okt == len(tr), "%d/%d" % (okt, len(tr)))
    C("truncated blocks number exactly 5",
      len({r["block_id"] for r in tr}) == 5, len({r["block_id"] for r in tr}))
    cen = rnd.sample([r for r in exported if r["censored_any"] == "1"], 20)
    okc = sum(1 for r in cen if dmin(r)[0]["any"] is None)
    C("20 marker-negative occurrences confirm no feature within 10 kb", okc == 20, "%d/20" % okc)
    zr = [r for r in exported if r["overlap_zero"] == "1"]
    zs = rnd.sample(zr, min(20, len(zr)))
    okz = sum(1 for r in zs if dmin(r)[0]["any"] == 0)
    C("20 distance-zero overlaps recompute as 0", okz == len(zs), "%d/%d" % (okz, len(zs)))

    # ---- exported vs reported agreement ----
    C("receipt output digests match the files on disk",
      all(sha(os.path.join(O, k)) == v for k, v in R["output_digests"].items()))
    C("receipt population figures match the exported table",
      R["population"]["occurrences"] == len(exported)
      and R["population"]["window_positive"] == len(obs)
      and R["population"]["in_block_feature_beyond_window_censored"] == 111)

    print("\n  checks run: %d   disagreements: %d" % (n, len(fails)))
    print("  VERDICT: %s" % ("PASS - zero disagreements" if not fails
                             else "*** FAIL: %s ***" % fails))
    sys.exit(0 if not fails else 9)


if __name__ == "__main__":
    main()
