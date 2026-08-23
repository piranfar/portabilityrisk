"""NM-V1 freeze and stratified draw.

Fixes the sampling design, then draws the sample deterministically from it. The draw is a
function of the frozen seed and of strata defined by the CURRENT HMM classification, which is
the method under test -- never by the reference method, which has not been run.

Sampling unit is the unique context BLOCK. 21,955 blocks carry 35,140 chromosomal ARG
occurrences; marker-positive blocks carry 2.46 ARGs against 1.23 for marker-negative, so
sampling occurrences would oversample marker-positive sequence about two-fold.
"""
import argparse, collections, csv, datetime, hashlib, json, os, random, sys

VERSION = "nmv1_freeze_and_sample_v1.0.0"
PROTOCOL_SHA = "c968fb6d16a528a64d064d6f8bbac745390804df4ad0897c09b12de84ca3fbff"
SEED = 20260821
ALLOC = {"S1_hmm_pos_IS_only": 350, "S2_hmm_pos_integron_only": 183,
         "S3_hmm_pos_both": 350, "S4_hmm_negative": 400}
LINEAR_FLOOR = 25
MAX_PER_BIOPROJECT_PER_STRATUM = 2


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--env-manifest", required=True)
    ap.add_argument("--out-design", required=True)
    ap.add_argument("--out-sample", required=True)
    a = ap.parse_args()
    if sha256_file(a.protocol) != PROTOCOL_SHA:
        print("REFUSING: protocol digest mismatch"); sys.exit(1)
    for p in (a.out_design, a.out_sample):
        if os.path.exists(p):
            print("REFUSING: %s exists; never overwritten" % p); sys.exit(1)
    print("%s | protocol verified %s | seed %d" % (VERSION, PROTOCOL_SHA[:16], SEED))

    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    rep2asm = {}
    for r in csv.DictReader(open(os.path.join(O, "replicon_level_summary.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        if r["replicon_molecule_type"].lower().startswith("chrom"):
            rep2asm[r["replicon_accession"]] = (r["assembly_version"], r["organism"])
    asm2bp = {r["assembly_version"]: r["bioproject_accession"]
              for r in csv.DictReader(open(os.path.join(O, "genome_level_summary.tsv"),
                                           encoding="utf-8"), delimiter="\t")}
    feat = collections.defaultdict(lambda: [0, 0])
    for r in csv.DictReader(open(os.path.join(O, "mge_feature_inventory.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        cl = r["feature_class"].lower()
        feat[r["block_id"]][0 if ("is" in cl or "transpos" in cl) else 1] += 1

    blocks = []
    for r in csv.DictReader(open(os.path.join(O, "shared_context_blocks.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        asm, sp = rep2asm.get(r["replicon_accession"], (None, None))
        f = feat.get(r["block_id"])
        if f is None:
            st = "S4_hmm_negative"
        elif f[0] and f[1]:
            st = "S3_hmm_pos_both"
        elif f[0]:
            st = "S1_hmm_pos_IS_only"
        else:
            st = "S2_hmm_pos_integron_only"
        blocks.append({"block_id": r["block_id"], "replicon": r["replicon_accession"],
                       "topology": r["topology"], "span_bp": int(r["block_span_bp"]),
                       "n_args": int(r["n_args_in_block"]),
                       "wrapped": r["wrapped_circular"], "truncated": r["truncated"],
                       "assembly": asm or "", "species": sp or "unknown",
                       "bioproject": asm2bp.get(asm, "unknown"), "stratum": st,
                       "hmm_is": f[0] if f else 0, "hmm_integron": f[1] if f else 0,
                       "hmm_positive": "yes" if f else "no"})
    pop = collections.Counter(b["stratum"] for b in blocks)
    print("\n=== STRATA (disjoint, defined by the HMM path under test) ===")
    for k in ALLOC:
        print("  %-28s population %6d  allocation %4d" % (k, pop[k], ALLOC[k]))
    print("  %-28s population %6d" % ("TOTAL", len(blocks)))

    rng = random.Random(SEED)
    chosen = []
    for st, want in ALLOC.items():
        pool = sorted([b for b in blocks if b["stratum"] == st],
                      key=lambda b: b["block_id"])
        if want >= len(pool):
            chosen.extend(pool)
            print("  %-28s took all %d" % (st, len(pool)))
            continue
        by_sp = collections.defaultdict(list)
        for b in pool:
            by_sp[b["species"]].append(b)
        for v in by_sp.values():
            rng.shuffle(v)
        order = sorted(by_sp, key=lambda s: (-len(by_sp[s]), s))
        picked = []
        bpcount = collections.Counter()
        # proportional-by-species round robin, capped per BioProject
        quota = {s: max(1, round(want * len(by_sp[s]) / len(pool))) for s in order}
        for s in order:
            for b in by_sp[s]:
                if len(picked) >= want or quota[s] <= 0:
                    break
                if bpcount[b["bioproject"]] >= MAX_PER_BIOPROJECT_PER_STRATUM:
                    continue
                picked.append(b); bpcount[b["bioproject"]] += 1; quota[s] -= 1
        # top up if quotas under-filled
        if len(picked) < want:
            got = {b["block_id"] for b in picked}
            for b in pool:
                if len(picked) >= want:
                    break
                if b["block_id"] in got:
                    continue
                if bpcount[b["bioproject"]] >= MAX_PER_BIOPROJECT_PER_STRATUM:
                    continue
                picked.append(b); bpcount[b["bioproject"]] += 1
        chosen.extend(picked)
        print("  %-28s drew %4d from %6d  (%d BioProjects)"
              % (st, len(picked), len(pool), len(bpcount)))

    got = {b["block_id"] for b in chosen}
    lin = [b for b in chosen if b["topology"] == "linear"]
    if len(lin) < LINEAR_FLOOR:
        extra = sorted([b for b in blocks
                        if b["topology"] == "linear" and b["block_id"] not in got],
                       key=lambda b: b["block_id"])
        rng.shuffle(extra)
        add = extra[:LINEAR_FLOOR - len(lin)]
        chosen.extend(add)
        print("  linear-topology top-up: added %d (floor %d)" % (len(add), LINEAR_FLOOR))

    chosen.sort(key=lambda b: b["block_id"])
    print("\n  TOTAL SAMPLED BLOCKS: %d" % len(chosen))
    print("  distinct BioProjects : %d" % len({b["bioproject"] for b in chosen}))
    print("  distinct species     : %d" % len({b["species"] for b in chosen}))
    print("  circular %d | linear %d" % (sum(1 for b in chosen if b["topology"] == "circular"),
                                         sum(1 for b in chosen if b["topology"] == "linear")))
    print("  total sequence span  : %.1f Mb" % (sum(b["span_bp"] for b in chosen) / 1e6))
    sp = collections.Counter(b["species"] for b in chosen)
    print("  top species:")
    for k, v in sp.most_common(5):
        print("     %-34s %4d" % (k, v))

    cols = ["block_id", "stratum", "hmm_positive", "hmm_is", "hmm_integron", "replicon",
            "assembly", "species", "bioproject", "topology", "span_bp", "n_args",
            "wrapped", "truncated"]
    with open(a.out_sample, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for b in chosen:
            fh.write("\t".join(str(b[c]) for c in cols) + "\n")

    D = {
     "design": "NMV1_FROZEN_DESIGN", "version": "1.0.0", "builder": VERSION,
     "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
     "authorised_by": "owner, D-NM2 and D-NM3 approved 2026-08-21",
     "protocol_sha256": PROTOCOL_SHA,
     "frozen_before_any_reference_tool_was_run": True,
     "sampling_unit": {"unit": "unique shared context block",
       "not": "ARG occurrence",
       "why": "21,955 blocks carry 35,140 occurrences; marker-positive blocks carry 2.46 ARGs "
              "against 1.23 for marker-negative, so occurrence sampling would oversample "
              "marker-positive sequence about two-fold",
       "each_block_enters_at_most_once": True},
     "seed": SEED, "seed_declared_before_sampling": True,
     "strata": {k: {"population": pop[k], "allocation": ALLOC[k]} for k in ALLOC},
     "strata_are_disjoint": True,
     "strata_defined_by": "the HMM path under test, never by the reference method",
     "additional_stratification": {
       "species": "proportional within stratum",
       "bioproject": "no more than %d blocks per BioProject per stratum" % MAX_PER_BIOPROJECT_PER_STRATUM,
       "topology": "linear oversampled to a floor of %d blocks" % LINEAR_FLOOR},
     "sample_size_justification": {
       "formula": "n = z^2 p(1-p)/d^2, z=1.96, d=0.03",
       "planning_p": 0.90,
       "why_0_90_not_0_95": "0.90 is the revise/failure boundary. Planning at the hoped-for "
                            "0.95 would under-power the very comparison that matters.",
       "required_reference_positive_blocks_per_arm": 384,
       "interval_method": "Wilson, not normal approximation",
       "explicit_rejection": "the round number 1,000 was not adopted; it corresponds to no "
                             "stated objective"},
     "arms": {
       "IS_arm": {"reference": "ISEScan 1.7.3",
                  "hmm_positive_strata": ["S1_hmm_pos_IS_only", "S3_hmm_pos_both"],
                  "hmm_negative_stratum": "S4_hmm_negative"},
       "integron_arm": {"reference": "IntegronFinder 2.0.6",
                        "hmm_positive_strata": ["S2_hmm_pos_integron_only", "S3_hmm_pos_both"],
                        "hmm_negative_stratum": "S4_hmm_negative"}},
     "arms_evaluated_separately": "pooling would let the abundant IS signal mask an integron "
                                  "detection failure",
     "method_non_equivalence": {
       "HMM_path": "transposase and integrase PROTEIN homology; cannot resolve element "
                   "boundaries or terminal inverted repeats",
       "ISEScan": "insertion sequences with resolvable terminal inverted repeats; cannot "
                  "detect degraded or truncated transposase proteins",
       "IntegronFinder": "integron integrases plus attC sites; the HMM path uses only 3 "
                         "integrase models and no attC detection, so the integron arm has "
                         "never been evaluated before",
       "no_method_here_is_truth": True},
     "metrics": {
       "always": ["percent agreement", "discordance count and direction per arm",
                  "recall of reference-positive blocks by the HMM path with Wilson 95% CI"],
       "only_on_adjudicated_subset": ["positive predictive value", "adjudicated sensitivity"],
       "prohibited": ["labelling HMM-only detections false positives outside the adjudicated "
                      "subset", "specificity against a non-truth reference",
                      "a single pooled figure across the IS and integron arms"]},
     "gates": {
       "success": "recall >= 0.95 with Wilson lower bound >= 0.90 in BOTH arms separately, "
                  "and no systematic HMM failure mode in the adjudicated subset",
       "revise": "recall 0.85-0.95, or lower bound 0.80-0.90, or the arms disagree",
       "failure": "recall < 0.85 in either arm, or a systematic failure mode"},
     "blinded_adjudication": {
       "trigger": "every block where the HMM path and the reference disagree",
       "size": "all discordant blocks up to 120, else a seeded subsample of 120",
       "blinding": "method labels removed, block ids replaced by opaque tokens; the "
                   "unblinding key is written and hashed before adjudication begins",
       "adjudicator": "OWNER DECISION REQUIRED - must be a person not involved in building "
                      "the pipeline. Claude built the analysis and therefore cannot adjudicate.",
       "status": "package prepared, adjudication NOT performed"},
     "environment": {"platform": "WSL2 Ubuntu 26.04 on the local laptop",
                     "isescan": "1.7.3", "integron_finder": "2.0.6",
                     "hmmer": "3.3.2", "prodigal": "2.6.3", "infernal": "1.1.4",
                     "blast": "2.17.0",
                     "manifest_sha256": sha256_file(a.env_manifest)},
     "sample": {"n_blocks": len(chosen),
                "n_bioprojects": len({b["bioproject"] for b in chosen}),
                "n_species": len({b["species"] for b in chosen}),
                "total_span_mb": round(sum(b["span_bp"] for b in chosen) / 1e6, 3),
                "manifest_sha256": sha256_file(a.out_sample)},
     "input_digests": {f: sha256_file(os.path.join(O, f)) for f in
                       ("shared_context_blocks.tsv", "mge_feature_inventory.tsv",
                        "arg_mge_neighbourhood.tsv", "replicon_level_summary.tsv",
                        "genome_level_summary.tsv")},
    }
    json.dump(D, open(a.out_design, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  sample manifest : %s\n    sha256 %s" % (a.out_sample, sha256_file(a.out_sample)))
    print("  FROZEN DESIGN   : %s\n    sha256 %s" % (a.out_design, sha256_file(a.out_design)))


if __name__ == "__main__":
    main()
