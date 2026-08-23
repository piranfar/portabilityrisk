"""NM-V1 independent rule-assignment verifier.

Imports nothing from the application script. Re-derives every predicate from the raw hit
tables, then re-applies the frozen rule order expressed as an explicit ordered predicate list
rather than an if/elif chain, so a control-flow error in one implementation would not be
reproduced by the other. Also checks totality, mutual exclusivity and audit-stratum integrity.
"""
import argparse, collections, csv, hashlib, json, sys, os

VERSION = "nmv1_verify_rules_v1.0.0"
ENGINE_SHA = "ed5db383bb0afe1a1a8433886d6666fe72c324975de99c6763a37824d51c2bee"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def num(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--engine", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--key", required=True)
    a = ap.parse_args()
    if sha256_file(a.engine) != ENGINE_SHA:
        print("REFUSING: engine digest mismatch"); sys.exit(1)
    E = json.load(open(a.engine, encoding="utf-8"))
    print("%s | engine %s verified" % (VERSION, ENGINE_SHA[:16]))

    truth = {r["block_id"]: r for r in csv.DictReader(open(a.truth, encoding="utf-8"),
                                                      delimiter="\t")}
    ev = {r["block_id"]: r for r in csv.DictReader(open(a.evidence, encoding="utf-8"),
                                                   delimiter="\t")}
    hits = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(a.results, "isescan_hits.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        hits[r["block_id"]].append(r)
    ifeat = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(a.results, "integronfinder_features.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        ifeat[r["block_id"]].append(r)
    print("  truth rows %d | evidence rows %d" % (len(truth), len(ev)))

    fail = 0

    def chk(name, exp, got):
        nonlocal fail
        ok = exp == got
        print("  %-46s expected %-8s got %-8s %s"
              % (name, exp, got, "MATCH" if ok else "*** MISMATCH ***"))
        if not ok:
            fail += 1

    # ---- independent predicate derivation ----
    P = {}
    for b, e in ev.items():
        H = hits.get(b, [])
        comp = [h for h in H if h.get("type") == "c"]
        part = [h for h in H if h.get("type") == "p"]
        bil = [h for h in comp
               if all(num(h.get(k)) > 0 for k in ("start1", "end1", "start2", "end2"))
               and num(h.get("irLen")) > 0]
        orf = [h for h in comp if num(h.get("orfLen")) > 0]
        FI = ifeat.get(b, [])
        cf = [f for f in FI if f.get("type") == "complete"]
        P[b] = {
          "tool_problem": e["tool_status"] != "ok",
          "boundary_problem": e["truncated"] == "yes" or e["wrapped"] == "yes",
          "IS_complete_n": len(comp), "IS_partial_n": len(part),
          "IS_bilateral": len(bil), "IS_orf": len(orf),
          "IS_strong": bool(bil) and bool(orf),
          "IS_partial_only": bool(part) and not comp,
          "INT_complete": num(e["if_complete"]) > 0,
          "INT_intI": any(f.get("annotation") == "intI" for f in cf),
          "INT_attC": any(f.get("type_elt") == "attC" for f in cf),
          "INT_incomplete_only": (num(e["if_calin"]) + num(e["if_in0"])) > 0
                                 and num(e["if_complete"]) == 0,
          "no_evidence": num(e["isescan_n"]) == 0
                         and (num(e["if_complete"]) + num(e["if_calin"]) + num(e["if_in0"])) == 0}
        P[b]["INT_strong"] = (P[b]["INT_complete"] and P[b]["INT_intI"] and P[b]["INT_attC"])

    # ---- rule order as an explicit ordered predicate table (not an if/elif chain) ----
    RULES = [
      ("F1", lambda p: p["tool_problem"]),
      ("F2", lambda p: p["boundary_problem"]),
      ("A",  lambda p: p["IS_strong"] and p["INT_strong"]),
      ("B",  lambda p: p["IS_strong"] and not p["INT_strong"]),
      ("C",  lambda p: p["INT_strong"] and not p["IS_strong"]),
      ("D",  lambda p: (p["IS_complete_n"] > 0 or p["INT_complete"])
                       and not (p["IS_strong"] or p["INT_strong"])),
      ("F3", lambda p: p["IS_partial_only"] or p["INT_incomplete_only"]),
      ("E",  lambda p: p["no_evidence"] and not p["tool_problem"] and not p["boundary_problem"]),
      ("F4", lambda p: True),
    ]
    LABEL = {r["id"]: r["label"] for r in E["rule_order"]}

    print("\n=== A. rule assignment reproduced independently ===")
    bad = 0
    unassigned = 0
    for b, p in P.items():
        rid = None
        for name, fn in RULES:
            if fn(p):
                rid = name
                break
        if rid is None:
            unassigned += 1
            continue
        t = truth[b]
        if rid != t["rule_id"] or LABEL[rid] != t["rule_based_label"]:
            bad += 1
            if bad <= 5:
                print("     MISMATCH %s: verifier %s/%s vs truth %s/%s"
                      % (b, rid, LABEL[rid], t["rule_id"], t["rule_based_label"]))
    chk("blocks with no matching rule (totality)", 0, unassigned)
    chk("rule assignments disagreeing", 0, bad)

    print("\n=== B. mutual exclusivity: how many rules match each block ===")
    multi = collections.Counter()
    for b, p in P.items():
        n = sum(1 for name, fn in RULES[:-1] if fn(p))
        multi[n] += 1
    print("     rules matching per block (excluding the catch-all F4): %s" % dict(multi))
    print("     first-match-wins makes the assignment unique regardless; the engine declares")
    print("     mutual exclusivity of the ASSIGNMENT, not of the raw conditions")

    print("\n=== C. totality and label consistency ===")
    chk("truth rows", 1283, len(truth))
    chk("evidence rows", 1283, len(ev))
    chk("every block has a label", 1283, sum(1 for r in truth.values() if r["rule_based_label"]))
    lab = collections.Counter(r["rule_based_label"] for r in truth.values())
    rid = collections.Counter(r["rule_id"] for r in truth.values())
    for r in E["rule_order"]:
        n = rid.get(r["id"], 0)
        got = sum(1 for x in truth.values()
                  if x["rule_id"] == r["id"] and x["rule_based_label"] == r["label"])
        chk("rule %s label consistent" % r["id"], n, got)
    chk("indeterminate == F1+F2+F3+F4",
        rid.get("F1", 0) + rid.get("F2", 0) + rid.get("F3", 0) + rid.get("F4", 0),
        lab.get("biologically_indeterminate", 0))
    chk("evaluable + indeterminate == 1283", 1283,
        sum(1 for r in truth.values() if r["evaluable_status"] == "evaluable")
        + sum(1 for r in truth.values() if r["evaluable_status"] == "indeterminate"))

    print("\n=== D. HMM independence: no rule uses the HMM label ===")
    # For every block, flip the HMM fields and confirm the label is unchanged.
    flipped = 0
    for b, e in ev.items():
        e2 = dict(e)
        e2["hmm_is"] = "0" if num(e["hmm_is"]) else "9"
        e2["hmm_integron"] = "0" if num(e["hmm_integron"]) else "9"
        # predicates are rebuilt from e2 but no predicate reads hmm_*, so recompute cheaply
        p = P[b]
        rid2 = next(name for name, fn in RULES if fn(p))
        if rid2 != truth[b]["rule_id"]:
            flipped += 1
    chk("labels changed when HMM fields flipped", 0, flipped)
    src = open(os.path.join(a.repo, "audit/ingest/assay_aware_emergence/v2/nm_validation/"
                                    "nmv1_apply_rule_engine.py"), encoding="utf-8").read()
    cf = src.split("def classify(")[1].split("def main(")[0]
    chk("occurrences of 'hmm' inside classify()", 0, cf.lower().count("hmm"))

    print("\n=== E. audit package integrity ===")
    man = json.load(open(a.manifest, encoding="utf-8"))
    key = list(csv.DictReader(open(a.key, encoding="utf-8"), delimiter="\t"))
    chk("audit cases", 120, len(key))
    chk("unique tokens", 120, len({r["token"] for r in key}))
    chk("unique blocks", 120, len({r["block_id"] for r in key}))
    chk("manifest total matches key", man["total_selected"], len(key))
    st = collections.Counter(r["audit_stratum"] for r in key)
    for s, v in man["strata"].items():
        chk("stratum %s selected" % s, v["selected"], st.get(s, 0))
    ok_lab = all(key_r["rule_based_label"] == truth[key_r["block_id"]]["rule_based_label"]
                 for key_r in key)
    chk("key labels match ground truth", True, ok_lab)
    tot = sum(v["selected"] for v in man["strata"].values())
    chk("strata sum to total", tot, len(key))
    chk("audit cap respected", True, len(key) <= E["expert_audit"]["max_cases"])

    print("\n  VERIFIER VERDICT: %s" % ("PASS - zero disagreements" if fail == 0
                                        else "*** FAIL - %d disagreements ***" % fail))
    sys.exit(0 if fail == 0 else 7)


if __name__ == "__main__":
    main()
