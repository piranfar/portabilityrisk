"""NM-V1 ground-truth construction and blinded adjudication package.

Verifies the frozen design and the frozen rules, audits the sample for substitution, builds
the per-block evidence table, applies the automatic label rules exactly as frozen, routes every
other case to adjudication, draws the QC arm, and writes the blinded package plus a hashed
unblinding key.

Computes tool-vs-tool concordance. Does NOT compute sensitivity, specificity, PPV or NPV
against ground truth: ground truth does not exist until adjudication is returned and frozen.
"""
import argparse, collections, csv, datetime, hashlib, json, os, sys
import numpy as np

VERSION = "nmv1_build_groundtruth_v1.0.0"
DESIGN_SHA = "c2aea6cb583c24b997ab376861acc600295e94d59bfa0ef2d55cf2bcc424bb20"
RULES_SHA = "e454873a89d1d56b20adb9ae157f20224076966daa0ed34bcbe44318063260f6"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--rules", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    if sha256_file(a.design) != DESIGN_SHA:
        print("REFUSING: design digest mismatch"); sys.exit(1)
    if sha256_file(a.rules) != RULES_SHA:
        print("REFUSING: rules digest mismatch"); sys.exit(1)
    D = json.load(open(a.design, encoding="utf-8"))
    R = json.load(open(a.rules, encoding="utf-8"))
    print("%s | design %s | rules %s verified" % (VERSION, DESIGN_SHA[:12], RULES_SHA[:12]))

    MAN = os.path.join(a.repo, "docs/nature_microbiology/nmv1_sample_manifest.tsv")
    if sha256_file(MAN) != D["sample"]["manifest_sha256"]:
        print("REFUSING: sample manifest changed since the freeze"); sys.exit(2)
    man = list(csv.DictReader(open(MAN, encoding="utf-8"), delimiter="\t"))
    print("  manifest verified: %d blocks" % len(man))

    # ---------- sample integrity audit ----------
    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    blk = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(O, "shared_context_blocks.tsv"), encoding="utf-8"), delimiter="\t")}
    print("\n=== SAMPLE INTEGRITY AUDIT ===")
    sub = 0; hdrbad = 0; nofasta = 0
    for r in man:
        b = r["block_id"]
        src = blk.get(b)
        if src is None or src["replicon_accession"] != r["replicon"]:
            sub += 1; continue
        fp = os.path.join(a.fasta, b + ".fna")
        if not os.path.exists(fp):
            nofasta += 1; continue
        with open(fp, encoding="utf-8") as fh:
            h = fh.readline().strip()
        want = ">%s|%s|%s-%s" % (b, src["replicon_accession"], src["block_start"], src["block_end"])
        if h != want:
            hdrbad += 1
    print("  blocks substituted or replicon-mismatched : %d" % sub)
    print("  blocks with FASTA header mismatch         : %d" % hdrbad)
    print("  blocks with no FASTA                      : %d" % nofasta)
    strata = collections.Counter(r["stratum"] for r in man)
    print("  strata represented : %d of 4 -> %s" % (len(strata), dict(strata)))
    print("  BioProjects        : %d (frozen expectation %d)"
          % (len({r["bioproject"] for r in man}), D["sample"]["n_bioprojects"]))
    print("  species            : %d (frozen expectation %d)"
          % (len({r["species"] for r in man}), D["sample"]["n_species"]))
    integrity_ok = (sub == 0 and hdrbad == 0 and nofasta == 0
                    and len(strata) == 4
                    and len({r["bioproject"] for r in man}) == D["sample"]["n_bioprojects"])
    print("  INTEGRITY: %s" % ("PASS" if integrity_ok else "*** FAIL ***"))

    # ---------- reference tool outputs ----------
    hits = collections.defaultdict(list)
    hp = os.path.join(a.results, "isescan_hits.tsv")
    if os.path.exists(hp):
        for r in csv.DictReader(open(hp, encoding="utf-8"), delimiter="\t"):
            hits[r["block_id"]].append(r)
    iss = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(a.results, "isescan_per_block.tsv"), encoding="utf-8"),
        delimiter="\t")} if os.path.exists(
        os.path.join(a.results, "isescan_per_block.tsv")) else {}
    inf = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(a.results, "integronfinder_per_block.tsv"), encoding="utf-8"),
        delimiter="\t")} if os.path.exists(
        os.path.join(a.results, "integronfinder_per_block.tsv")) else {}

    # ---------- evidence table ----------
    strat_pop = D["strata"]
    frac = {k: strat_pop[k]["allocation"] / strat_pop[k]["population"] for k in strat_pop}
    rows = []
    for r in man:
        b = r["block_id"]
        hh = hits.get(b, [])
        i = iss.get(b); f = inf.get(b)
        # prefer the per-block columns; fall back to recounting the hit rows
        ncomp = int(i["isescan_complete"]) if (i and "isescan_complete" in i) \
            else sum(1 for x in hh if x.get("type") == "c")
        npart = int(i["isescan_partial"]) if (i and "isescan_partial" in i) \
            else sum(1 for x in hh if x.get("type") == "p")
        irs = [x for x in hh if x.get("irLen") not in (None, "", "0")]
        status = "ok"
        if i is None or f is None:
            status = "tool_failure_or_missing_output"
        elif f.get("if_status", "ok") != "ok":
            status = f["if_status"]
        rows.append({
            "block_id": b, "stratum": r["stratum"], "weight": 1.0 / frac[r["stratum"]],
            "species": r["species"], "bioproject": r["bioproject"],
            "topology": r["topology"], "span_bp": int(r["span_bp"]),
            "wrapped": r["wrapped"], "truncated": r["truncated"],
            "hmm_is": int(r["hmm_is"]), "hmm_integron": int(r["hmm_integron"]),
            "hmm_is_pos": int(r["hmm_is"]) > 0,
            "hmm_int_pos": int(r["hmm_integron"]) > 0,
            "isescan_n": int(i["isescan_n_is"]) if i else 0,
            "isescan_pos": (i["isescan_positive"] == "yes") if i else False,
            "isescan_complete": ncomp, "isescan_partial": npart,
            "isescan_with_ir": len(irs),
            "if_complete": int(f["if_complete"]) if f else 0,
            "if_calin": int(f["if_calin"]) if f else 0,
            "if_in0": int(f["if_in0"]) if f else 0,
            "if_pos": (f["if_positive"] == "yes") if f else False,
            "tool_status": status})

    # ---------- frozen automatic label rules ----------
    for r in rows:
        lab = None; route = None
        if r["tool_status"] != "ok":
            route = "TOOL_FAILURE"
        elif r["truncated"] == "yes" or r["wrapped"] == "yes":
            route = "BOUNDARY_AMBIGUOUS"
        elif r["hmm_is_pos"] != r["isescan_pos"]:
            route = "DISCORDANT_IS"
        elif r["hmm_int_pos"] != r["if_pos"]:
            route = "DISCORDANT_INTEGRON"
        elif r["isescan_pos"] and r["isescan_complete"] == 0:
            route = "PARTIAL_ONLY"
        elif (r["if_calin"] + r["if_in0"]) > 0 and r["if_complete"] == 0:
            route = "INCOMPLETE_INTEGRON"
        elif r["hmm_is_pos"] and r["isescan_complete"] > 0:
            lab = "AUTO_IS_MOBILE"
        elif r["hmm_int_pos"] and r["if_complete"] > 0:
            lab = "AUTO_INTEGRON"
        elif (not r["hmm_is_pos"] and not r["hmm_int_pos"] and r["isescan_n"] == 0
              and (r["if_complete"] + r["if_calin"] + r["if_in0"]) == 0):
            lab = "AUTO_QUIESCENT"
        else:
            route = "REFERENCE_DISAGREEMENT"
        r["auto_label"] = lab or ""
        r["adjudication_route"] = route or ""
    auto = collections.Counter(r["auto_label"] for r in rows if r["auto_label"])
    adj = collections.Counter(r["adjudication_route"] for r in rows if r["adjudication_route"])
    print("\n=== FROZEN RULE APPLICATION ===")
    print("  auto-labelled: %d" % sum(auto.values()))
    for k, v in auto.most_common():
        print("     %-20s %5d" % (k, v))
    print("  routed to adjudication: %d" % sum(adj.values()))
    for k, v in adj.most_common():
        print("     %-24s %5d" % (k, v))

    # ---------- QC arm ----------
    rng = np.random.default_rng(R["quality_control_arm"]["seed"])
    qc = []
    for lab, n in R["quality_control_arm"]["allocation"].items():
        pool = sorted([r for r in rows if r["auto_label"] == lab],
                      key=lambda x: x["block_id"])
        take = min(n, len(pool))
        if take:
            idx = rng.choice(len(pool), take, replace=False)
            for j in sorted(idx):
                qc.append(pool[j])
        print("  QC arm %-18s requested %3d  available %5d  drawn %3d"
              % (lab, n, len(pool), take))
    qcset = {r["block_id"] for r in qc}
    pkg = [r for r in rows if r["adjudication_route"]] + qc
    pkg.sort(key=lambda x: x["block_id"])
    print("  adjudication package: %d discordant/ambiguous + %d QC = %d blocks"
          % (sum(adj.values()), len(qc), len(pkg)))

    # ---------- blinded package, rubric, unblinding key ----------
    os.makedirs(a.outdir, exist_ok=True)
    order = sorted(pkg, key=lambda x: hashlib.sha256(
        (x["block_id"] + "|nmv1blind").encode()).hexdigest())
    tok = {r["block_id"]: "ADJ%04d" % (i + 1) for i, r in enumerate(order)}

    pk = os.path.join(a.outdir, "NMV1_ADJUDICATION_BLINDED_PACKAGE.tsv")
    cols = ["token", "topology", "block_span_bp", "truncated_or_wrapped",
            "methodX_is_markers", "methodX_integron_markers",
            "methodY_elements_total", "methodY_complete_elements", "methodY_partial_elements",
            "methodY_elements_with_inverted_repeat",
            "methodZ_complete_integrons", "methodZ_calin", "methodZ_in0",
            "tool_status", "adjudication", "reason"]
    with open(pk, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in order:
            fh.write("\t".join(str(x) for x in [
                tok[r["block_id"]], r["topology"], r["span_bp"],
                "yes" if (r["truncated"] == "yes" or r["wrapped"] == "yes") else "no",
                r["hmm_is"], r["hmm_integron"],
                r["isescan_n"], r["isescan_complete"], r["isescan_partial"],
                r["isescan_with_ir"],
                r["if_complete"], r["if_calin"], r["if_in0"],
                r["tool_status"], "", ""]) + "\n")

    rub = os.path.join(a.outdir, "NMV1_ADJUDICATION_RUBRIC.md")
    RB = R["adjudication_rubric"]
    with open(rub, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# NM-V1 blinded adjudication rubric\n\n")
        fh.write("**Adjudicator:** %s\n\n" % R["adjudicator"]["name"])
        fh.write("Frozen rules SHA-256 `%s`.\n\n" % RULES_SHA)
        fh.write("## What you are looking at\n\n")
        fh.write("Each row is one genomic block. Three methods are shown, de-identified as\n")
        fh.write("**Method X**, **Method Y** and **Method Z**. You are not told which is which,\n")
        fh.write("nor the species, study, gene identity or original classification.\n\n")
        fh.write("| column | meaning |\n|---|---|\n")
        fh.write("| `methodX_*` | marker counts from a protein-homology method |\n")
        fh.write("| `methodY_*` | element counts from a boundary-based method, with complete "
                 "vs partial and inverted-repeat evidence |\n")
        fh.write("| `methodZ_*` | integron structures: complete, CALIN, In0 |\n\n")
        fh.write("## Permitted outcomes\n\nEnter exactly one in the `adjudication` column:\n\n")
        for o in R["permitted_adjudication_outcomes"]:
            fh.write("- `%s`\n" % o)
        fh.write("\n## Decision rules\n\n")
        for k in ("IS_associated_supported", "integron_associated_supported",
                  "multiple_MGE_evidence_supported", "chromosomal_mobile_supported",
                  "chromosomal_quiescent_supported", "neither_classification_supported",
                  "biologically_indeterminate", "boundary_cases"):
            fh.write("**%s** — %s\n\n" % (k, RB[k]))
        fh.write("## Instruction\n\n%s\n\n" % RB["instruction"])
        fh.write("Put a short free-text justification in `reason`. Leave no row blank.\n")

    key = os.path.join(a.outdir, "NMV1_ADJUDICATION_UNBLINDING_KEY.tsv")
    with open(key, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("token\tblock_id\tstratum\tspecies\tbioproject\troute\tauto_label\tis_qc_control\n")
        for r in order:
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                     % (tok[r["block_id"]], r["block_id"], r["stratum"], r["species"],
                        r["bioproject"], r["adjudication_route"] or "QC",
                        r["auto_label"], "yes" if r["block_id"] in qcset else "no"))

    # ---------- design-weighted tool concordance (NOT ground-truth performance) ----------
    def conc(hk, rk, label):
        W = lambda L: sum(x["weight"] for x in L)
        a11 = [r for r in rows if r[hk] and r[rk]]
        a10 = [r for r in rows if r[hk] and not r[rk]]
        a01 = [r for r in rows if not r[hk] and r[rk]]
        a00 = [r for r in rows if not r[hk] and not r[rk]]
        w11, w10, w01, w00 = W(a11), W(a10), W(a01), W(a00)
        rec = w11 / (w11 + w01) if (w11 + w01) else float("nan")
        rng2 = np.random.default_rng(20260821)
        bys = collections.defaultdict(list)
        for r in rows:
            bys[r["stratum"]].append(r)
        bt = []
        for _ in range(2000):
            n11 = n01 = 0.0
            for st, L in bys.items():
                for j in rng2.integers(0, len(L), len(L)):
                    x = L[j]
                    if x[rk]:
                        if x[hk]:
                            n11 += x["weight"]
                        else:
                            n01 += x["weight"]
            if n11 + n01 > 0:
                bt.append(n11 / (n11 + n01))
        bt = np.array(bt)
        return {"arm": label,
                "unweighted": {"hmm+ref+": len(a11), "hmm+ref-": len(a10),
                               "hmm-ref+": len(a01), "hmm-ref-": len(a00)},
                "weighted": {"hmm+ref+": w11, "hmm+ref-": w10,
                             "hmm-ref+": w01, "hmm-ref-": w00},
                "recall_design_weighted": rec,
                "recall_ci": [float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))],
                "recall_unweighted": len(a11) / (len(a11) + len(a01)) if (len(a11) + len(a01)) else float("nan"),
                "agreement_weighted": (w11 + w00) / (w11 + w10 + w01 + w00)}
    print("\n=== TOOL-VS-TOOL CONCORDANCE (design-weighted; NOT ground-truth performance) ===")
    A = conc("hmm_is_pos", "isescan_pos", "IS")
    C = conc("hmm_int_pos", "if_pos", "integron")
    for x in (A, C):
        u = x["unweighted"]
        print("  %s arm: HMM+ref+ %d | HMM+ref- %d | HMM-ref+ %d | HMM-ref- %d"
              % (x["arm"], u["hmm+ref+"], u["hmm+ref-"], u["hmm-ref+"], u["hmm-ref-"]))
        print("     recall design-weighted %.4f  95%% CI [%.4f, %.4f] | unweighted %.4f | agreement %.4f"
              % (x["recall_design_weighted"], x["recall_ci"][0], x["recall_ci"][1],
                 x["recall_unweighted"], x["agreement_weighted"]))

    ev = os.path.join(a.outdir, "nmv1_block_evidence_table.tsv")
    ecols = ["block_id", "stratum", "weight", "species", "bioproject", "topology", "span_bp",
             "wrapped", "truncated", "hmm_is", "hmm_integron", "isescan_n",
             "isescan_complete", "isescan_partial", "isescan_with_ir",
             "if_complete", "if_calin", "if_in0", "tool_status",
             "auto_label", "adjudication_route"]
    with open(ev, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(ecols) + "\n")
        for r in sorted(rows, key=lambda x: x["block_id"]):
            fh.write("\t".join(str(r[c]) for c in ecols) + "\n")

    rec = {"builder": VERSION,
           "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "design_sha256": DESIGN_SHA, "rules_sha256": RULES_SHA,
           "environment": D["environment"],
           "integrity": {"blocks_in_manifest": len(man), "substituted": sub,
                         "header_mismatch": hdrbad, "missing_fasta": nofasta,
                         "strata_represented": dict(strata),
                         "bioprojects": len({r["bioproject"] for r in man}),
                         "species": len({r["species"] for r in man}),
                         "pass": bool(integrity_ok)},
           "auto_labels": dict(auto), "adjudication_routes": dict(adj),
           "qc_arm": {"allocation": R["quality_control_arm"]["allocation"],
                      "drawn": len(qc), "seed": R["quality_control_arm"]["seed"]},
           "adjudication_package": {"n_blocks": len(pkg),
                                    "n_discordant_or_ambiguous": sum(adj.values()),
                                    "n_qc_controls": len(qc),
                                    "blinded_sha256": sha256_file(pk),
                                    "rubric_sha256": sha256_file(rub),
                                    "unblinding_key_sha256": sha256_file(key),
                                    "key_hashed_before_delivery": True,
                                    "adjudication_performed": False},
           "tool_concordance_not_ground_truth": {"IS": A, "integron": C},
           "ground_truth_status": "NOT YET ESTABLISHED - requires the returned, frozen and "
                                  "hashed adjudication file",
           "outputs": {"evidence_table_sha256": sha256_file(ev)},
           "statements": [
             "Ground truth is the frozen adjudicated labels, not any single tool output.",
             "Sensitivity, specificity, PPV and NPV are NOT computed here; they require the "
             "returned adjudication.",
             "Recall figures above are HMM-versus-reference concordance only.",
             "Design weights are mandatory because strata were sampled at 2.6 to 100 per cent.",
             "No transfer, conjugation or HGT event was observed or is claimed."]}
    rp = os.path.join(a.outdir, "NMV1_GROUNDTRUTH_BUILD_RECEIPT.json")
    json.dump(rec, open(rp, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  evidence table : %s" % sha256_file(ev))
    print("  BLINDED PACKAGE: %s" % sha256_file(pk))
    print("  RUBRIC         : %s" % sha256_file(rub))
    print("  UNBLINDING KEY : %s" % sha256_file(key))
    print("  build receipt  : %s" % sha256_file(rp))


if __name__ == "__main__":
    main()
