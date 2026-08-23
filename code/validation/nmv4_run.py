"""NM-V4 -- generalisation of the discordance principle. Step 4 of the frozen order.

Refuses to run unless NMV4_FROZEN_DESIGN.json matches the digest fixed before any outcome was
computed. Every threshold, species list and gate below comes from that file, not from here.

The confirmation is non-circular by construction: the relationship between plasmid fraction and
chromosomal MGE fraction is fitted on species that were never used to articulate the principle,
and only then applied to Acinetobacter baumannii.
"""
import argparse, collections, csv, datetime, hashlib, json, os, sys
import numpy as np

VERSION = "nmv4_run_v1.0.0"
FROZEN_SHA = "d06e1e68e9dd5a9303ed12c86748f687e71de7afbc594753da947a4a8c7357d9"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def logit(p, n):
    lo, hi = 0.5 / n, 1 - 0.5 / n
    p = min(max(p, lo), hi)
    return float(np.log(p / (1 - p)))


def ols(x, y):
    X = np.column_stack([np.ones(len(x)), x])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--frozen", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    if sha256_file(a.frozen) != FROZEN_SHA:
        print("REFUSING: frozen design digest mismatch"); sys.exit(1)
    F = json.load(open(a.frozen, encoding="utf-8"))
    if not F.get("frozen_before_any_outcome_was_computed"):
        print("REFUSING: design does not assert pre-outcome freeze"); sys.exit(1)
    CONF = F["confirmation_species"]
    DISC = F["discovery_species"]
    SEED = F["uncertainty"]["seed"]
    print("%s | frozen design verified %s" % (VERSION, FROZEN_SHA[:16]))
    print("  confirmation species: %d | discovery species: %d | seed %d"
          % (len(CONF), len(DISC), SEED))

    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    for k, v in F["input_digests"].items():
        if sha256_file(os.path.join(O, k)) != v:
            print("REFUSING: input %s changed since the freeze" % k); sys.exit(2)
    print("  all four frozen input digests verified")

    # ---------- genome level: species, bioproject, plasmid/chromosome counts ----------
    gen = {}
    for r in csv.DictReader(open(os.path.join(O, "genome_level_summary.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        gen[r["assembly_version"]] = {
            "sp": r["organism"], "bp": r["bioproject_accession"],
            "occ": int(r["n_arg_occurrences"]), "pl": int(r["n_plasmid_args"])}

    # ---------- replicon -> assembly, restricted to chromosomes carrying ARGs ----------
    rep2asm = {}
    for r in csv.DictReader(open(os.path.join(O, "replicon_level_summary.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        if r["replicon_molecule_type"].lower().startswith("chrom"):
            rep2asm[r["replicon_accession"]] = r["assembly_version"]

    # ---------- blocks and their marker status ----------
    posblocks = set()
    for r in csv.DictReader(open(os.path.join(O, "mge_feature_inventory.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        posblocks.add(r["block_id"])
    blocks_by_asm = collections.defaultdict(list)   # assembly -> list of bool(positive)
    nblk = 0
    for r in csv.DictReader(open(os.path.join(O, "shared_context_blocks.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        asm = rep2asm.get(r["replicon_accession"])
        if asm is None:
            continue
        blocks_by_asm[asm].append(r["block_id"] in posblocks)
        nblk += 1
    print("  blocks mapped to a genome: %d of 21955" % nblk)

    by_sp = collections.defaultdict(list)
    for asm, g in gen.items():
        by_sp[g["sp"]].append(asm)

    def measure(asms):
        occ = pl = 0; bt = bp_ = 0
        for asm in asms:
            g = gen[asm]; occ += g["occ"]; pl += g["pl"]
            for isp in blocks_by_asm.get(asm, ()):
                bt += 1; bp_ += 1 if isp else 0
        P = pl / occ if occ else float("nan")
        M = bp_ / bt if bt else float("nan")
        return P, M, occ, bt

    rows = {}
    for sp in CONF + DISC:
        P, M, occ, bt = measure(by_sp[sp])
        rows[sp] = {"P": P, "M": M, "n_occ": occ, "n_blocks": bt,
                    "lP": logit(P, occ), "lM": logit(M, bt)}

    print("\n=== per-species measurements ===")
    print("  %-32s %8s %8s %9s %8s %8s" % ("species", "P", "M", "n_occ", "blocks", "D"))
    for sp in CONF:
        r = rows[sp]
        print("  %-32s %8.4f %8.4f %9d %8d %+8.3f"
              % (sp, r["P"], r["M"], r["n_occ"], r["n_blocks"], r["lM"] - r["lP"]))
    for sp in DISC:
        r = rows[sp]
        print("  %-32s %8.4f %8.4f %9d %8d %+8.3f  [DISCOVERY]"
              % (sp, r["P"], r["M"], r["n_occ"], r["n_blocks"], r["lM"] - r["lP"]))

    Pconf = np.array([rows[s]["P"] for s in CONF])
    print("\n  confirmation-set plasmid-fraction range: %.4f to %.4f" % (Pconf.min(), Pconf.max()))
    ab = rows["Acinetobacter baumannii"]["P"]
    print("  Acinetobacter baumannii P = %.4f -> %s the confirmation range"
          % (ab, "INSIDE" if Pconf.min() <= ab <= Pconf.max() else "OUTSIDE (extrapolation)"))

    # ---------------- T2: fit on confirmation only, predict A. baumannii ----------------
    x = np.array([rows[s]["lP"] for s in CONF])
    y = np.array([rows[s]["lM"] for s in CONF])
    b = ols(x, y)
    pred_ab = b[0] + b[1] * rows["Acinetobacter baumannii"]["lP"]
    resid_ab = rows["Acinetobacter baumannii"]["lM"] - pred_ab
    pred_kp = b[0] + b[1] * rows["Klebsiella pneumoniae"]["lP"]
    resid_kp = rows["Klebsiella pneumoniae"]["lM"] - pred_kp
    yhat = b[0] + b[1] * x
    ss_res = float(np.sum((y - yhat) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
    print("\n=== T2: fit on the 8 confirmation species, predict the discovery species ===")
    print("  fitted line: logit(M) = %+.4f %+.4f * logit(P)   in-sample R2 = %.4f"
          % (b[0], b[1], 1 - ss_res / ss_tot if ss_tot else float("nan")))
    print("  A. baumannii: observed logit(M) %+.3f | predicted %+.3f | RESIDUAL %+.3f"
          % (rows["Acinetobacter baumannii"]["lM"], pred_ab, resid_ab))
    print("     on the proportion scale: observed M %.4f | predicted %.4f"
          % (rows["Acinetobacter baumannii"]["M"], 1 / (1 + np.exp(-pred_ab))))
    print("  K. pneumoniae: observed %+.3f | predicted %+.3f | residual %+.3f"
          % (rows["Klebsiella pneumoniae"]["lM"], pred_kp, resid_kp))

    # ---------------- T1: leave-one-species-out within the confirmation set -------------
    print("\n=== T1: leave-one-species-out within the confirmation set ===")
    loso = []
    for i, s in enumerate(CONF):
        m = np.ones(len(CONF), bool); m[i] = False
        bb = ols(x[m], y[m])
        p = bb[0] + bb[1] * x[i]
        loso.append(y[i] - p)
        print("  %-32s observed %+7.3f  predicted %+7.3f  residual %+7.3f"
              % (s, y[i], p, y[i] - p))
    loso = np.array(loso)
    # sampling noise floor: binomial SE of logit(M) is 1/sqrt(n*M*(1-M))
    se = np.array([1.0 / np.sqrt(rows[s]["n_blocks"] * rows[s]["M"] * (1 - rows[s]["M"]))
                   for s in CONF])
    print("  LOSO residual SD          : %.4f" % float(np.std(loso, ddof=1)))
    print("  mean within-species SE    : %.4f" % float(se.mean()))
    print("  ratio SD/SE               : %.2f  (>>1 means P alone does not determine M)"
          % (float(np.std(loso, ddof=1)) / float(se.mean())))

    # ---------------- T3: the low-plasmid control --------------------------------------
    pa = rows["Pseudomonas aeruginosa"]
    abr = rows["Acinetobacter baumannii"]
    print("\n=== T3: low-plasmid control ===")
    print("  P. aeruginosa   P %.4f  M %.4f" % (pa["P"], pa["M"]))
    print("  A. baumannii    P %.4f  M %.4f" % (abr["P"], abr["M"]))
    print("  plasmid fractions differ by %.4f; chromosomal MGE fractions differ by %.4f"
          % (abs(abr["P"] - pa["P"]), abs(abr["M"] - pa["M"])))

    # ---------------- BioProject-clustered bootstrap -----------------------------------
    print("\n=== BioProject-clustered bootstrap within species, %d resamples ===" % a.boot)
    bp_by_sp = collections.defaultdict(lambda: collections.defaultdict(list))
    for asm, g in gen.items():
        bp_by_sp[g["sp"]][g["bp"]].append(asm)
    rng = np.random.default_rng(SEED)
    R_ab, R_kp, SLOPE, GAP, LOSOSD = [], [], [], [], []
    for it in range(a.boot):
        rs = {}
        ok = True
        for sp in CONF + DISC:
            bps = list(bp_by_sp[sp])
            pick = rng.choice(len(bps), size=len(bps), replace=True)
            asms = [asm for j in pick for asm in bp_by_sp[sp][bps[j]]]
            P, M, occ, bt = measure(asms)
            if not occ or not bt or M <= 0 or M >= 1:
                ok = False; break
            rs[sp] = (logit(P, occ), logit(M, bt), P, M)
        if not ok:
            continue
        xb = np.array([rs[s][0] for s in CONF]); yb = np.array([rs[s][1] for s in CONF])
        bb = ols(xb, yb)
        R_ab.append(rs["Acinetobacter baumannii"][1] - (bb[0] + bb[1] * rs["Acinetobacter baumannii"][0]))
        R_kp.append(rs["Klebsiella pneumoniae"][1] - (bb[0] + bb[1] * rs["Klebsiella pneumoniae"][0]))
        SLOPE.append(bb[1])
        GAP.append(rs["Acinetobacter baumannii"][3] - rs["Pseudomonas aeruginosa"][3])
        lo = []
        for i in range(len(CONF)):
            m = np.ones(len(CONF), bool); m[i] = False
            c = ols(xb[m], yb[m]); lo.append(yb[i] - (c[0] + c[1] * xb[i]))
        LOSOSD.append(float(np.std(lo, ddof=1)))
        if (it + 1) % 500 == 0:
            print("  %d/%d" % (it + 1, a.boot), flush=True)

    ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    R_ab = np.array(R_ab); R_kp = np.array(R_kp); GAP = np.array(GAP); SLOPE = np.array(SLOPE)
    LOSOSD = np.array(LOSOSD)
    a_lo, a_hi = ci(R_ab); k_lo, k_hi = ci(R_kp); g_lo, g_hi = ci(GAP)
    s_lo, s_hi = ci(SLOPE); l_lo, l_hi = ci(LOSOSD)
    print("\n  resamples completed: %d" % len(R_ab))
    print("  A. baumannii residual   median %+.3f  95%% CI [%+.3f, %+.3f]  excludes 0: %s"
          % (float(np.median(R_ab)), a_lo, a_hi, "YES" if a_lo > 0 or a_hi < 0 else "NO"))
    print("  K. pneumoniae residual  median %+.3f  95%% CI [%+.3f, %+.3f]  excludes 0: %s"
          % (float(np.median(R_kp)), k_lo, k_hi, "YES" if k_lo > 0 or k_hi < 0 else "NO"))
    print("  M gap A.b minus P.a     median %+.3f  95%% CI [%+.3f, %+.3f]  excludes 0: %s"
          % (float(np.median(GAP)), g_lo, g_hi, "YES" if g_lo > 0 or g_hi < 0 else "NO"))
    print("  fitted slope on logit(P) median %+.3f 95%% CI [%+.3f, %+.3f]"
          % (float(np.median(SLOPE)), s_lo, s_hi))
    print("  LOSO residual SD        median %.3f   95%% CI [%.3f, %.3f]"
          % (float(np.median(LOSOSD)), l_lo, l_hi))

    # ---------------- gate evaluation --------------------------------------------------
    g_success = (a_lo > 0) and (g_lo > 0) and (l_lo > float(se.mean()))
    g_direction = float(np.median(R_ab)) > 0
    verdict = ("SUCCESS" if g_success else
               ("REVISE" if g_direction else "FAILURE"))
    print("\n=== GATE ===")
    print("  T2 A. baumannii residual positive and CI excludes 0 : %s" % (a_lo > 0))
    print("  T3 A.b vs P.a chromosomal MGE gap CI excludes 0     : %s" % (g_lo > 0))
    print("  T1 LOSO residual SD exceeds sampling noise          : %s"
          % (l_lo > float(se.mean())))
    print("  VERDICT: %s" % verdict)

    rec = {"builder": VERSION,
           "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "frozen_design_sha256": FROZEN_SHA,
           "protocol_sha256": F["protocol_sha256"],
           "confirmation_species": CONF, "discovery_species": DISC,
           "per_species": {s: {k: (None if isinstance(v, float) and np.isnan(v) else v)
                               for k, v in rows[s].items()} for s in rows},
           "confirmation_P_range": [float(Pconf.min()), float(Pconf.max())],
           "ab_P_inside_confirmation_range": bool(Pconf.min() <= ab <= Pconf.max()),
           "T2": {"intercept": float(b[0]), "slope": float(b[1]),
                  "in_sample_r2": float(1 - ss_res / ss_tot) if ss_tot else None,
                  "ab_observed_logitM": rows["Acinetobacter baumannii"]["lM"],
                  "ab_predicted_logitM": float(pred_ab),
                  "ab_residual": float(resid_ab),
                  "ab_observed_M": rows["Acinetobacter baumannii"]["M"],
                  "ab_predicted_M": float(1 / (1 + np.exp(-pred_ab))),
                  "kp_residual": float(resid_kp),
                  "bootstrap_median": float(np.median(R_ab)), "bootstrap_ci": [a_lo, a_hi],
                  "excludes_zero": bool(a_lo > 0 or a_hi < 0)},
           "T1": {"loso_residuals": {CONF[i]: float(loso[i]) for i in range(len(CONF))},
                  "loso_residual_sd": float(np.std(loso, ddof=1)),
                  "mean_within_species_se": float(se.mean()),
                  "sd_over_se": float(np.std(loso, ddof=1) / se.mean()),
                  "bootstrap_sd_ci": [l_lo, l_hi]},
           "T3": {"pa_P": pa["P"], "pa_M": pa["M"], "ab_P": abr["P"], "ab_M": abr["M"],
                  "M_gap_median": float(np.median(GAP)), "M_gap_ci": [g_lo, g_hi],
                  "excludes_zero": bool(g_lo > 0 or g_hi < 0)},
           "bootstrap": {"unit": "BioProject within species", "seed": SEED,
                         "n_requested": a.boot, "n_completed": int(len(R_ab)),
                         "slope_median": float(np.median(SLOPE)), "slope_ci": [s_lo, s_hi]},
           "gate_verdict": verdict,
           "statements": [
             "The fit uses only the 8 confirmation species, none of which was used to "
             "articulate the discordance principle.",
             "M is block-weighted, not occurrence-weighted.",
             "No transfer, conjugation or HGT event was observed or is claimed.",
             "No PortabilityEvent or PlasmidCall artefact was read."]}
    json.dump(rec, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  receipt: %s\n  sha256 %s" % (a.out, sha256_file(a.out)))


if __name__ == "__main__":
    main()
