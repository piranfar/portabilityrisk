"""NM-V1 scoring -- stratified MGE validation.

Recall is DESIGN-WEIGHTED. The four strata were sampled at very different rates (S2 at 100 %,
S4 at 2.6 %), so an unweighted recall would be badly biased towards the over-sampled
HMM-positive strata. Every block carries weight 1/(sampling fraction of its stratum), and the
confidence interval comes from a stratified bootstrap rather than a binomial formula that
would assume simple random sampling.

Homology-only detections are never called false positives: no independent truth exists for
them outside the adjudicated subset, which is prepared here but not adjudicated.
"""
import argparse, collections, csv, datetime, hashlib, json, math, os, sys
import numpy as np

VERSION = "nmv1_score_v1.0.0"
FROZEN_SHA = "c2aea6cb583c24b997ab376861acc600295e94d59bfa0ef2d55cf2bcc424bb20"
IS_POS_STRATA = {"S1_hmm_pos_IS_only", "S3_hmm_pos_both"}
INT_POS_STRATA = {"S2_hmm_pos_integron_only", "S3_hmm_pos_both"}
NEG_STRATUM = "S4_hmm_negative"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z / d * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return p, max(0.0, c - h), min(1.0, c + h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--frozen", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    if sha256_file(a.frozen) != FROZEN_SHA:
        print("REFUSING: frozen design digest mismatch"); sys.exit(1)
    F = json.load(open(a.frozen, encoding="utf-8"))
    if not F.get("frozen_before_any_reference_tool_was_run"):
        print("REFUSING: design does not assert pre-reference freeze"); sys.exit(1)
    print("%s | frozen design verified %s" % (VERSION, FROZEN_SHA[:16]))

    man = list(csv.DictReader(open(os.path.join(
        a.repo, "docs/nature_microbiology/nmv1_sample_manifest.tsv"),
        encoding="utf-8"), delimiter="\t"))
    if sha256_file(os.path.join(a.repo, "docs/nature_microbiology/nmv1_sample_manifest.tsv")) \
            != F["sample"]["manifest_sha256"]:
        print("REFUSING: sample manifest changed since the freeze"); sys.exit(2)
    print("  sample manifest verified: %d blocks" % len(man))

    iss = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(a.results, "isescan_per_block.tsv"), encoding="utf-8"),
        delimiter="\t")}
    inf = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(a.results, "integronfinder_per_block.tsv"), encoding="utf-8"),
        delimiter="\t")}
    miss = [r["block_id"] for r in man if r["block_id"] not in iss or r["block_id"] not in inf]
    if miss:
        print("  WARNING: %d sampled blocks have no reference result" % len(miss))
        print("     first few: %s" % miss[:5])
    print("  reference results: ISEScan %d | IntegronFinder %d" % (len(iss), len(inf)))

    strata = F["strata"]
    frac = {k: strata[k]["allocation"] / strata[k]["population"] for k in strata}
    print("\n=== SAMPLING FRACTIONS (why recall must be design-weighted) ===")
    for k in sorted(strata):
        print("  %-28s %5d of %6d = %6.2f%%   weight %8.2f"
              % (k, strata[k]["allocation"], strata[k]["population"],
                 100 * frac[k], 1 / frac[k]))

    rows = []
    for r in man:
        b = r["block_id"]
        if b not in iss or b not in inf:
            continue
        st = r["stratum"]
        rows.append({
            "block_id": b, "stratum": st, "weight": 1.0 / frac[st],
            "species": r["species"], "bioproject": r["bioproject"],
            "topology": r["topology"], "span_bp": int(r["span_bp"]),
            "hmm_is_pos": int(r["hmm_is"]) > 0,
            "hmm_int_pos": int(r["hmm_integron"]) > 0,
            "hmm_any_pos": r["hmm_positive"] == "yes",
            "isescan_pos": iss[b]["isescan_positive"] == "yes",
            "isescan_n": int(iss[b]["isescan_n_is"]),
            "if_pos": inf[b]["if_positive"] == "yes",
            "if_complete": int(inf[b]["if_complete"]),
            "if_calin": int(inf[b]["if_calin"]), "if_in0": int(inf[b]["if_in0"])})
    print("  blocks scored: %d" % len(rows))

    def arm(name, hmm_key, ref_key):
        R = rows
        a11 = [r for r in R if r[hmm_key] and r[ref_key]]
        a10 = [r for r in R if r[hmm_key] and not r[ref_key]]
        a01 = [r for r in R if not r[hmm_key] and r[ref_key]]
        a00 = [r for r in R if not r[hmm_key] and not r[ref_key]]
        W = lambda L: sum(x["weight"] for x in L)
        w11, w10, w01, w00 = W(a11), W(a10), W(a01), W(a00)
        rec_w = w11 / (w11 + w01) if (w11 + w01) else float("nan")
        rec_u, lo_u, hi_u = wilson(len(a11), len(a11) + len(a01))
        agree_w = (w11 + w00) / (w11 + w10 + w01 + w00)
        rng = np.random.default_rng(20260821)
        bystr = collections.defaultdict(list)
        for r in R:
            bystr[r["stratum"]].append(r)
        boot = []
        for _ in range(a.boot):
            n11 = n01 = 0.0
            for st, L in bystr.items():
                idx = rng.integers(0, len(L), len(L))
                for i in idx:
                    x = L[i]
                    if x[ref_key]:
                        if x[hmm_key]:
                            n11 += x["weight"]
                        else:
                            n01 += x["weight"]
            if n11 + n01 > 0:
                boot.append(n11 / (n11 + n01))
        boot = np.array(boot)
        lo, hi = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))) \
            if len(boot) else (float("nan"), float("nan"))
        return {"arm": name, "n_blocks": len(R),
                "cells_unweighted": {"hmm+ref+": len(a11), "hmm+ref-": len(a10),
                                     "hmm-ref+": len(a01), "hmm-ref-": len(a00)},
                "cells_weighted": {"hmm+ref+": w11, "hmm+ref-": w10,
                                   "hmm-ref+": w01, "hmm-ref-": w00},
                "recall_design_weighted": rec_w, "recall_ci_lo": lo, "recall_ci_hi": hi,
                "recall_unweighted": rec_u, "recall_unweighted_wilson": [lo_u, hi_u],
                "agreement_weighted": agree_w,
                "discordant_hmm_only": len(a10), "discordant_ref_only": len(a01),
                "discordant_total": len(a10) + len(a01)}

    print("\n=== ARM A: IS / transposase, reference = ISEScan 1.7.3 ===")
    A = arm("IS", "hmm_is_pos", "isescan_pos")
    print("\n=== ARM C: integron, reference = IntegronFinder 2.0.6 ===")
    C = arm("integron", "hmm_int_pos", "if_pos")
    for r in (A, C):
        c = r["cells_unweighted"]
        print("  %s arm 2x2 (unweighted block counts):" % r["arm"])
        print("     HMM+ ref+ %5d | HMM+ ref- %5d" % (c["hmm+ref+"], c["hmm+ref-"]))
        print("     HMM- ref+ %5d | HMM- ref- %5d" % (c["hmm-ref+"], c["hmm-ref-"]))
        print("     recall, DESIGN-WEIGHTED : %.4f  95%% CI [%.4f, %.4f]  (stratified bootstrap)"
              % (r["recall_design_weighted"], r["recall_ci_lo"], r["recall_ci_hi"]))
        print("     recall, unweighted      : %.4f  Wilson [%.4f, %.4f]  (biased, shown for contrast)"
              % (r["recall_unweighted"], *r["recall_unweighted_wilson"]))
        print("     weighted agreement      : %.4f" % r["agreement_weighted"])
        print("     discordant blocks       : %d (HMM-only %d, reference-only %d)"
              % (r["discordant_total"], r["discordant_hmm_only"], r["discordant_ref_only"]))
        print()

    # ---- gates ----
    def gate(r):
        if r["recall_design_weighted"] >= 0.95 and r["recall_ci_lo"] >= 0.90:
            return "SUCCESS"
        if r["recall_design_weighted"] >= 0.85 or r["recall_ci_lo"] >= 0.80:
            return "REVISE"
        return "FAILURE"
    gA, gC = gate(A), gate(C)
    overall = "SUCCESS" if gA == "SUCCESS" and gC == "SUCCESS" else \
              ("FAILURE" if "FAILURE" in (gA, gC) else "REVISE")
    print("=== GATES (evaluated per arm, never pooled) ===")
    print("  IS arm       : %s" % gA)
    print("  integron arm : %s" % gC)
    print("  OVERALL      : %s" % overall)

    # ---- blinded adjudication package ----
    disc = [r for r in rows
            if (r["hmm_is_pos"] != r["isescan_pos"]) or (r["hmm_int_pos"] != r["if_pos"])]
    rng = np.random.default_rng(20260821)
    if len(disc) > 120:
        pick = sorted(rng.choice(len(disc), 120, replace=False))
        adj = [disc[i] for i in pick]
    else:
        adj = disc
    os.makedirs(a.outdir, exist_ok=True)
    tok = {}
    for i, r in enumerate(sorted(adj, key=lambda x: x["block_id"])):
        tok[r["block_id"]] = "ADJ%04d" % (i + 1)
    pk = os.path.join(a.outdir, "nmv1_adjudication_blinded.tsv")
    with open(pk, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("token\tspan_bp\ttopology\tn_candidate_features\tadjudication\treason\n")
        for b, t in sorted(tok.items(), key=lambda x: x[1]):
            r = [x for x in adj if x["block_id"] == b][0]
            fh.write("%s\t%d\t%s\t%d\t\t\n"
                     % (t, r["span_bp"], r["topology"],
                        r["isescan_n"] + r["if_complete"] + r["if_calin"] + r["if_in0"]))
    key = os.path.join(a.outdir, "nmv1_adjudication_unblinding_key.tsv")
    with open(key, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("token\tblock_id\tstratum\thmm_is\thmm_integron\tisescan\tintegronfinder\n")
        for b, t in sorted(tok.items(), key=lambda x: x[1]):
            r = [x for x in adj if x["block_id"] == b][0]
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                     % (t, b, r["stratum"], r["hmm_is_pos"], r["hmm_int_pos"],
                        r["isescan_pos"], r["if_pos"]))
    print("\n=== BLINDED ADJUDICATION PACKAGE ===")
    print("  discordant blocks total : %d" % len(disc))
    print("  in package              : %d" % len(adj))
    print("  blinded file  : %s  sha256 %s" % (os.path.basename(pk), sha256_file(pk)))
    print("  unblinding key: %s  sha256 %s" % (os.path.basename(key), sha256_file(key)))
    print("  ADJUDICATION NOT PERFORMED - requires a person not involved in building the")
    print("  pipeline. Claude built the analysis and cannot adjudicate. D-NM3 needs a name.")

    p1 = os.path.join(a.outdir, "nmv1_block_level_results.tsv")
    cols = ["block_id", "stratum", "weight", "species", "bioproject", "topology", "span_bp",
            "hmm_is_pos", "hmm_int_pos", "isescan_pos", "isescan_n", "if_pos",
            "if_complete", "if_calin", "if_in0"]
    with open(p1, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: x["block_id"]):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    rec = {"builder": VERSION,
           "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "frozen_design_sha256": FROZEN_SHA,
           "environment": F["environment"],
           "n_sampled": len(man), "n_scored": len(rows),
           "blocks_without_reference_result": len(miss),
           "sampling_fractions": frac,
           "arm_IS": A, "arm_integron": C,
           "gates": {"IS": gA, "integron": gC, "overall": overall},
           "adjudication": {"discordant_total": len(disc), "package_size": len(adj),
                            "blinded_sha256": sha256_file(pk),
                            "unblinding_key_sha256": sha256_file(key),
                            "performed": False,
                            "blocker": "D-NM3 adjudicator not named; must not be the analysis author"},
           "outputs": {"block_level_tsv": sha256_file(p1)},
           "statements": [
             "Recall is design-weighted by inverse sampling fraction; unweighted recall is "
             "reported only for contrast and is biased.",
             "Homology-only detections are NOT called false positives; no truth exists for "
             "them outside the adjudicated subset.",
             "The two arms are evaluated separately and never pooled.",
             "No transfer, conjugation or HGT event was observed or is claimed."]}
    rp = os.path.join(a.outdir, "NMV1_RESULT_RECEIPT.json")
    json.dump(rec, open(rp, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  %s  %s" % (os.path.basename(p1), sha256_file(p1)))
    print("  %s  %s" % (os.path.basename(rp), sha256_file(rp)))


if __name__ == "__main__":
    main()
